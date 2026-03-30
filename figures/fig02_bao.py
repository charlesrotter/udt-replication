#!/usr/bin/env python3
"""
fig02 -- BAO distance-redshift diagram.

D_V/r_d vs redshift with BOSS + DESI data.
r_d fitted as single free parameter (tests shape of D_V(z)).
UDT does NOT use r_d = 147 Mpc (LCDM sound horizon).

Saves: manuscript/figures/fig02_bao.{pdf,png}
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
from scipy.optimize import brentq, minimize_scalar

# === Publication style ===
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['CMU Serif', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'font.size': 9,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
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

CB_BLUE   = '#0072B2'
CB_GREEN  = '#009E73'
CB_ORANGE = '#E69F00'
CB_RED    = '#D55E00'

# === UDT cosmological model ===
def phi_cosmo(r):
    return COSMO_K * r + COSMO_BETA * r**2 + COSMO_GAMMA * r**3

def phi_prime(r):
    return COSMO_K + 2 * COSMO_BETA * r + 3 * COSMO_GAMMA * r**2

def r_of_z(z):
    target = np.log(1.0 + z)
    try:
        return brentq(lambda r: phi_cosmo(r) - target, 1e-6, 60.0)
    except ValueError:
        return np.nan

def D_V_Gpc(z):
    """Volume-averaged distance D_V(z) in Gpc."""
    r = r_of_z(z)
    pp = phi_prime(r)
    if pp <= 0 or r <= 0 or np.isnan(r):
        return np.nan
    return (z * r**2 / pp) ** (1.0 / 3.0)

# === BAO data ===
# BOSS DR12 (Alam et al. 2017)
boss = np.array([
    [0.38, 10.27, 0.15],
    [0.51, 13.38, 0.18],
    [0.61, 15.45, 0.25],
])

# DESI Y1 (2024)
desi = np.array([
    [0.30,  7.93, 0.15],
    [0.51, 13.62, 0.25],
    [0.71, 16.85, 0.32],
    [0.93, 21.71, 0.28],
    [1.32, 27.79, 0.69],
    [1.49, 30.69, 0.80],
])

# Ly-alpha (DESI, least LCDM-contaminated)
lya = np.array([
    [2.33, 39.71, 0.95],
])

# Combine all
all_z = np.concatenate([boss[:, 0], desi[:, 0], lya[:, 0]])
all_dv_rd = np.concatenate([boss[:, 1], desi[:, 1], lya[:, 1]])
all_err = np.concatenate([boss[:, 2], desi[:, 2], lya[:, 2]])

# Compute UDT D_V at data redshifts
DV_model = np.array([D_V_Gpc(z) for z in all_z])

# Fit r_d (single parameter): minimize chi2 of D_V_model/r_d vs data D_V/r_d
def chi2_rd(rd_Gpc):
    pred = DV_model / rd_Gpc
    return np.sum(((pred - all_dv_rd) / all_err)**2)

res = minimize_scalar(chi2_rd, bounds=(0.05, 0.3), method='bounded')
rd_best = res.x
rd_Mpc = rd_best * 1000
rms_pct = np.sqrt(np.mean(((DV_model / rd_best - all_dv_rd) / all_dv_rd)**2)) * 100
print(f"  Best-fit r_d = {rd_Mpc:.1f} Mpc ({rd_best:.4f} Gpc)")
print(f"  RMS: {rms_pct:.1f}%")

# Model curve
z_curve = np.linspace(0.15, 2.5, 300)
dv_curve = np.array([D_V_Gpc(z) / rd_best for z in z_curve])

# === Create figure ===
fig, ax = plt.subplots(figsize=(3.4, 3.0))

# UDT curve
ax.plot(z_curve, dv_curve, '-', color=CB_RED, linewidth=1.2,
        label='UDT prediction', zorder=3)

# BOSS
ax.errorbar(boss[:, 0], boss[:, 1], yerr=boss[:, 2], fmt='s',
            color=CB_BLUE, markersize=4.5, elinewidth=0.6,
            capsize=2, capthick=0.4, label='BOSS DR12', zorder=4)

# DESI
ax.errorbar(desi[:, 0], desi[:, 1], yerr=desi[:, 2], fmt='D',
            color=CB_GREEN, markersize=3.5, elinewidth=0.6,
            capsize=2, capthick=0.4, label='DESI Y1', zorder=4)

# Ly-alpha
ax.errorbar(lya[:, 0], lya[:, 1], yerr=lya[:, 2], fmt='^',
            color=CB_ORANGE, markersize=5, elinewidth=0.6,
            capsize=2, capthick=0.4, label=r'Ly$\alpha$ forest', zorder=4)

ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$D_V / r_d$')
ax.set_xlim(0, 2.5)
ax.set_ylim(0, 45)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.legend(loc='upper left', frameon=True, fancybox=False,
          edgecolor='0.6', framealpha=0.95, fontsize=7,
          handlelength=1.5)

# Annotation
ax.text(0.96, 0.08,
        f'$r_d = {rd_Mpc:.0f}$ Mpc (fitted)\nRMS = {rms_pct:.1f}%',
        transform=ax.transAxes, fontsize=6.5,
        ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFFDE7',
                  edgecolor='0.6', alpha=0.9, linewidth=0.4))

fig.tight_layout()

# === Save ===
savefig(fig, 'fig02_bao')
plt.close(fig)
print("fig02_bao: done")
