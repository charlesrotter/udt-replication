#!/usr/bin/env python3
"""
fig01 -- Type Ia supernova Hubble diagram with Pantheon+ data.

Upper panel: distance modulus mu vs redshift z.
Lower panel: residuals.
Uses REAL Pantheon+SH0ES data, not synthetic.

Saves: manuscript/figures/fig01_sne_hubble.{pdf,png}
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from lib.constants import COSMO_K, COSMO_BETA, COSMO_GAMMA
from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import brentq

# === Publication style ===
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['CMU Serif', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'font.size': 9,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.3,
    'ytick.minor.width': 0.3,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'savefig.dpi': 300,
})

CB_BLUE = '#0072B2'
CB_RED  = '#D55E00'
CB_GRAY = '#999999'

# === Load Pantheon+ data ===
data_paths = [
    os.path.join(os.path.dirname(__file__), '..', 'data', 'external',
                 'Pantheon+SH0ES.dat'),
    os.path.join(os.path.dirname(__file__), '..', 'data', 'external',
                 'pantheon_plus', 'Pantheon+SH0ES.dat'),
]

data_file = None
for p in data_paths:
    if os.path.exists(p):
        data_file = p
        break

if data_file is None:
    raise FileNotFoundError("Pantheon+SH0ES.dat not found. Run download_data.sh.")

print(f"  Loading: {data_file}")

# Parse the whitespace-separated file
z_hd_all, mu_all, mu_err_all = [], [], []
with open(data_file) as f:
    header = f.readline()  # skip header
    for line in f:
        cols = line.split()
        if len(cols) < 12:
            continue
        z_hd = float(cols[2])       # zHD (Hubble-flow corrected)
        m_b = float(cols[8])        # m_b_corr
        m_b_err = float(cols[9])    # m_b_corr_err_DIAG
        z_hd_all.append(z_hd)
        mu_all.append(m_b)
        mu_err_all.append(m_b_err)

z_hd_all = np.array(z_hd_all)
mu_all = np.array(mu_all)
mu_err_all = np.array(mu_err_all)

# Filter: z > 0.01 and reasonable errors
mask = (z_hd_all > 0.01) & (mu_err_all < 1.0) & (mu_err_all > 0)
z_data = z_hd_all[mask]
mu_data = mu_all[mask]
mu_err = mu_err_all[mask]
print(f"  SNe after cuts: {len(z_data)}")

# === UDT model ===
def phi_cosmo(r):
    return COSMO_K * r + COSMO_BETA * r**2 + COSMO_GAMMA * r**3

def mu_udt(z_arr, M_offset=0.0):
    """Distance modulus from UDT polynomial. M_offset is the free parameter."""
    result = np.full_like(z_arr, np.nan, dtype=float)
    for i, z in enumerate(z_arr):
        target = np.log(1.0 + z)
        try:
            r_sol = brentq(lambda r: phi_cosmo(r) - target, 1e-6, 50.0)
        except ValueError:
            continue
        D_L = (1 + z) * r_sol  # Gpc
        D_L_pc = D_L * 1e9
        if D_L_pc > 0:
            result[i] = 5.0 * np.log10(D_L_pc / 10.0) + M_offset
    return result

# Fit M_offset (absolute magnitude) to minimize residuals
from scipy.optimize import minimize_scalar
def chi2_M(M):
    model = mu_udt(z_data, M)
    valid = ~np.isnan(model)
    return np.sum(((mu_data[valid] - model[valid]) / mu_err[valid])**2)

res = minimize_scalar(chi2_M, bounds=(-25, -15), method='bounded')
M_best = res.x
print(f"  Best-fit M offset: {M_best:.4f}")

# Model curve
z_model = np.linspace(0.01, 2.3, 500)
mu_model = mu_udt(z_model, M_best)

# Residuals
mu_pred_data = mu_udt(z_data, M_best)
residuals = mu_data - mu_pred_data

# Bin data for cleaner visualization
z_bins = np.linspace(0.01, 2.3, 40)
z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])
mu_binned = np.full(len(z_centers), np.nan)
mu_binned_err = np.full(len(z_centers), np.nan)
res_binned = np.full(len(z_centers), np.nan)
res_binned_err = np.full(len(z_centers), np.nan)

for i in range(len(z_centers)):
    in_bin = (z_data >= z_bins[i]) & (z_data < z_bins[i+1])
    if np.sum(in_bin) < 3:
        continue
    w = 1.0 / mu_err[in_bin]**2
    mu_binned[i] = np.average(mu_data[in_bin], weights=w)
    mu_binned_err[i] = 1.0 / np.sqrt(np.sum(w))
    res_binned[i] = np.average(residuals[in_bin], weights=w)
    res_binned_err[i] = mu_binned_err[i]

valid = ~np.isnan(mu_binned)

# RMS
rms = np.sqrt(np.nanmean(residuals**2))
print(f"  RMS residual: {rms:.3f} mag")

# === Create figure ===
fig, (ax_main, ax_res) = plt.subplots(
    2, 1, figsize=(3.4, 3.8),
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05},
    sharex=True
)

# --- Main panel ---
# Individual SNe as small dots (the data cloud)
ax_main.scatter(z_data, mu_data, s=1.5, c=CB_GRAY, zorder=1, rasterized=True,
                edgecolors='none', alpha=0.4,
                label=f'Pantheon+ ({len(z_data)} SNe)')

# UDT curve — on top of scatter cloud
ax_main.plot(z_model, mu_model, '-', color=CB_RED, linewidth=1.3,
             label='UDT geometric polynomial', zorder=5)

ax_main.set_ylabel(r'Distance modulus $\mu$ (mag)')
ax_main.legend(loc='lower right', frameon=True, fancybox=False,
               edgecolor='0.6', framealpha=0.95, fontsize=6.5,
               handlelength=1.5)
# Auto-scale y to data range
mu_lo = np.nanmin(mu_model) - 0.5
mu_hi = np.nanmax(mu_model) + 0.5
ax_main.set_xlim(0, 2.35)
ax_main.set_ylim(mu_lo, mu_hi)
ax_main.yaxis.set_minor_locator(AutoMinorLocator())
ax_main.tick_params(labelbottom=False)

# Polynomial label
ax_main.text(0.04, 0.95,
             r'$\phi(r) = \frac{3}{2}\mu_g r - \cos\frac{\pi}{5}\,\mu_g^2 r^2'
             r' + \frac{2}{3}\mu_g^3 r^3$',
             transform=ax_main.transAxes, fontsize=6, va='top',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       edgecolor='0.6', alpha=0.9, linewidth=0.4))

# --- Residual panel ---
ax_res.axhline(0, color='0.5', linewidth=0.5, zorder=1)
ax_res.axhspan(-0.1, 0.1, color='#E8F4E8', alpha=0.4, zorder=0)
ax_res.errorbar(z_centers[valid], res_binned[valid], yerr=res_binned_err[valid],
                fmt='o', color=CB_BLUE, markersize=3.5, elinewidth=0.6,
                capsize=1.5, capthick=0.4, zorder=3)

ax_res.set_xlabel(r'Redshift $z$')
ax_res.set_ylabel(r'$\Delta\mu$ (mag)')
ax_res.set_ylim(-0.4, 0.4)
ax_res.yaxis.set_minor_locator(AutoMinorLocator())
ax_res.xaxis.set_minor_locator(AutoMinorLocator())

# RMS label
ax_res.text(0.96, 0.90, f'RMS = {rms:.3f} mag',
            transform=ax_res.transAxes, fontsize=6.5,
            ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='0.6', alpha=0.9, linewidth=0.4))

fig.align_ylabels([ax_main, ax_res])

# === Save ===
savefig(fig, 'fig01_sne_hubble')
plt.close(fig)
print("fig01_sne_hubble: done")
