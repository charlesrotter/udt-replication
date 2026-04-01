#!/usr/bin/env python3
"""UDT Joint Likelihood Analysis: SNe + BAO + CMB simultaneously.

One free cosmological parameter: mu_g (screening mass).
Three geometric parameters locked by Diophantine identity.

SNe: d_L = r*(1+z), M offset analytically marginalized
BAO: D_M/r_d (transverse), native r_d = pi*r_CMB*1000/ell_A_micro
CMB: WKB2 peaks with Diophantine q1=2/9, q2=-1/20, peak 1 anchored

Addresses: "No global likelihood across SNe+BAO+CMB simultaneously."
"""
import os, sys
import numpy as np
from scipy.optimize import brentq, minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from lib.constants import PHI_CMB

pi = np.pi

# === Locked microphysical constants ===
R_STAR_MICRO = 6.9875
E2_E1 = 5.9017
I2 = 0.82296
ell_A_micro = 2 * pi * R_STAR_MICRO * E2_E1 / I2  # 314.8

# WKB2 corrections from Diophantine quantum numbers
q1 = 2.0 / 9.0     # (2j+1)/(2l+1)^2
q2 = -1.0 / 20.0   # -1/((2j+1)^2*(2kmax-1))

mu_g_micro = pi * np.sqrt(pi / 3) / 13  # 0.24730

print("=" * 60)
print("UDT Joint Likelihood Analysis")
print("=" * 60)
print(f"  mu_g_micro = {mu_g_micro:.5f} Gpc^-1")
print(f"  ell_A_micro = {ell_A_micro:.1f}")
print(f"  q1 = 2/9 = {q1:.6f},  q2 = -1/20 = {q2:.6f}")

# === Cosmological profile (parameterized by mu_g) ===
def phi_cosmo(r, mg):
    k = 1.5 * mg
    b = -np.cos(pi / 5) * mg**2
    g = (2.0 / 3.0) * mg**3
    return k * r + b * r**2 + g * r**3

def dphi_cosmo(r, mg):
    k = 1.5 * mg
    b = -np.cos(pi / 5) * mg**2
    g = (2.0 / 3.0) * mg**3
    return k + 2 * b * r + 3 * g * r**2

def r_of_z(z, mg):
    try:
        return brentq(lambda r: np.exp(phi_cosmo(r, mg)) - 1 - z, 1e-6, 60, xtol=1e-10)
    except:
        return np.nan

def distance_modulus(dL_Gpc):
    if np.isnan(dL_Gpc) or dL_Gpc <= 0:
        return np.nan
    return 5.0 * np.log10(dL_Gpc) + 40.0

# === Load Pantheon+ ===
data_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'external')
pantheon_file = os.path.join(data_dir, 'Pantheon+SH0ES.dat')

sne_z, sne_mu, sne_err = [], [], []
if os.path.isfile(pantheon_file):
    with open(pantheon_file) as f:
        f.readline()  # header
        for line in f:
            cols = line.split()
            if len(cols) < 12: continue
            z = float(cols[2])
            mb = float(cols[8])
            mbe = float(cols[9])
            if z > 0.01 and 0 < mbe < 5:
                sne_z.append(z); sne_mu.append(mb); sne_err.append(mbe)
    print(f"  Loaded {len(sne_z)} Pantheon+ SNe")
else:
    print("  WARNING: Pantheon+ not found, using representative data")
    sne_z = [0.01,0.05,0.1,0.2,0.4,0.7,1.0,1.5,2.0]
    sne_mu = [33.1,36.6,38.2,39.8,41.4,42.5,43.2,43.9,44.4]
    sne_err = [0.15]*len(sne_z)
sne_z = np.array(sne_z); sne_mu = np.array(sne_mu); sne_err = np.array(sne_err)

# === BAO data (8 surveys, VR §3) ===
bao_data = np.array([
    [0.11, 3.047, 0.137], [0.38, 10.27, 0.15], [0.51, 13.38, 0.18],
    [0.61, 15.33, 0.22], [0.70, 17.86, 0.25], [0.51, 13.62, 0.18],
    [1.48, 30.69, 0.50], [2.33, 39.71, 0.70],
])

# === CMB peaks (Planck 2018) ===
planck_peaks = np.array([220.0, 537.5, 810.8, 1120.9, 1444.2, 1776.0, 2081.0])
cmb_sigma = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0])

# === Chi-squared: SNe ===
def chi2_sne(mg):
    mu_pred = np.array([distance_modulus(r_of_z(z, mg) * (1 + z)) for z in sne_z])
    valid = np.isfinite(mu_pred)
    if valid.sum() < 10: return 1e10, np.inf
    mp, md, me = mu_pred[valid], sne_mu[valid], sne_err[valid]
    w = 1.0 / me**2
    M_eff = np.sum(w * (md - mp)) / np.sum(w)
    resid = md - mp - M_eff
    return float(np.sum(resid**2 / me**2)), float(np.sqrt(np.mean(resid**2)))

# === Chi-squared: BAO (D_M/r_d) ===
def chi2_bao(mg):
    r_cmb = r_of_z(1090.0, mg)
    if np.isnan(r_cmb): return 1e10, np.inf
    rd = pi * r_cmb * 1000 / ell_A_micro  # Mpc
    errs = []
    chi2 = 0.0
    for row in bao_data:
        z, obs, sig = row
        r = r_of_z(z, mg)
        if np.isnan(r): return 1e10, np.inf
        dm_rd = r * 1000 / rd
        chi2 += ((dm_rd - obs) / sig)**2
        errs.append((dm_rd / obs - 1))
    rms = np.sqrt(np.mean(np.array(errs)**2)) * 100
    return float(chi2), float(rms)

# === Chi-squared: CMB peaks (WKB2 + Diophantine q, peak 1 anchored) ===
def chi2_cmb(mg):
    r_cmb = r_of_z(1090.0, mg)
    if np.isnan(r_cmb): return 1e10, np.inf
    dp = dphi_cosmo(r_cmb, mg)
    eta = 2.0 / r_cmb - 2.0 * dp
    c_eff = np.arctan(eta) / pi
    # Anchor peak 1
    nu1 = 1.0 + c_eff
    F1 = nu1 + q1 / nu1 + q2 / nu1**3
    ell_A = planck_peaks[0] / F1
    chi2 = 0.0
    errs = []
    for n in range(1, 8):
        nu = n + c_eff
        ell_pred = ell_A * (nu + q1 / nu + q2 / nu**3)
        chi2 += ((ell_pred - planck_peaks[n-1]) / cmb_sigma[n-1])**2
        errs.append((ell_pred / planck_peaks[n-1] - 1))
    rms = np.sqrt(np.mean(np.array(errs)**2)) * 100
    return float(chi2), float(rms)

def chi2_joint(mg):
    return chi2_sne(mg)[0] + chi2_bao(mg)[0] + chi2_cmb(mg)[0]

# === Sanity check at mu_g_micro ===
print("\n--- Sanity checks at mu_g_micro ---")
_, rms_s = chi2_sne(mu_g_micro)
_, rms_b = chi2_bao(mu_g_micro)
_, rms_c = chi2_cmb(mu_g_micro)
r_cmb = r_of_z(1090.0, mu_g_micro)
rd = pi * r_cmb * 1000 / ell_A_micro
print(f"  r_CMB = {r_cmb:.4f} Gpc,  r_d = {rd:.1f} Mpc")
print(f"  RMS SNe = {rms_s:.3f} mag,  BAO = {rms_b:.1f}%,  CMB = {rms_c:.2f}%")

# === Scan mu_g ===
print("\n--- Scanning mu_g ---")
mg_scan = np.linspace(0.20, 0.30, 200)
c2_s = np.array([chi2_sne(mg)[0] for mg in mg_scan])
c2_b = np.array([chi2_bao(mg)[0] for mg in mg_scan])
c2_c = np.array([chi2_cmb(mg)[0] for mg in mg_scan])
c2_j = c2_s + c2_b + c2_c

# === Individual & joint best fits ===
res_s = minimize_scalar(lambda mg: chi2_sne(mg)[0], bounds=(0.20, 0.30), method='bounded')
res_b = minimize_scalar(lambda mg: chi2_bao(mg)[0], bounds=(0.20, 0.30), method='bounded')
res_c = minimize_scalar(lambda mg: chi2_cmb(mg)[0], bounds=(0.20, 0.30), method='bounded')
res_j = minimize_scalar(chi2_joint, bounds=(0.20, 0.30), method='bounded')

mg_s, mg_b, mg_c, mg_j = res_s.x, res_b.x, res_c.x, res_j.x
c2j_min = res_j.fun

cs_s, rs_s = chi2_sne(mg_s)
cs_b, rs_b = chi2_bao(mg_b)
cs_c, rs_c = chi2_cmb(mg_c)

cj_s, rj_s = chi2_sne(mg_j)
cj_b, rj_b = chi2_bao(mg_j)
cj_c, rj_c = chi2_cmb(mg_j)

cm_s, rm_s = chi2_sne(mu_g_micro)
cm_b, rm_b = chi2_bao(mu_g_micro)
cm_c, rm_c = chi2_cmb(mu_g_micro)
c2_micro = cm_s + cm_b + cm_c

N_s, N_b, N_c = len(sne_z), len(bao_data), len(planck_peaks)
N_tot = N_s + N_b + N_c

# === 1-sigma interval ===
try:
    mg_lo = brentq(lambda mg: chi2_joint(mg) - c2j_min - 1, 0.20, mg_j, xtol=1e-6)
except: mg_lo = np.nan
try:
    mg_hi = brentq(lambda mg: chi2_joint(mg) - c2j_min - 1, mg_j, 0.30, xtol=1e-6)
except: mg_hi = np.nan

micro_in_1s = (not np.isnan(mg_lo)) and (not np.isnan(mg_hi)) and (mg_lo <= mu_g_micro <= mg_hi)

# === Results ===
print("\n" + "=" * 60)
print("RESULTS")
print("=" * 60)
print(f"{'Dataset':>14}  {'mu_g':>8}  {'chi2/N':>8}  {'RMS':>12}")
print("-" * 50)
print(f"{'SNe ('+str(N_s)+')':>14}  {mg_s:.5f}  {cs_s/N_s:.3f}  {rs_s:.3f} mag")
print(f"{'BAO ('+str(N_b)+')':>14}  {mg_b:.5f}  {cs_b/N_b:.3f}  {rs_b:.1f}%")
print(f"{'CMB ('+str(N_c)+')':>14}  {mg_c:.5f}  {cs_c/N_c:.3f}  {rs_c:.2f}%")
print(f"{'micro':>14}  {mu_g_micro:.5f}  {c2_micro/N_tot:.3f}  ---")

print(f"\nJoint fit (1 free: mu_g):")
sig_lo = f"{mg_j-mg_lo:.5f}" if not np.isnan(mg_lo) else "---"
sig_hi = f"{mg_hi-mg_j:.5f}" if not np.isnan(mg_hi) else "---"
print(f"  mu_g = {mg_j:.5f} +{sig_hi} -{sig_lo}")
print(f"  chi2/N = {c2j_min/N_tot:.3f}  ({N_tot} points)")
print(f"  SNe: {rj_s:.3f} mag  |  BAO: {rj_b:.1f}%  |  CMB: {rj_c:.2f}%")

print(f"\nMicrophysics mu_g (0 fitted):")
print(f"  chi2/N = {c2_micro/N_tot:.3f}")
print(f"  SNe: {rm_s:.3f} mag  |  BAO: {rm_b:.1f}%  |  CMB: {rm_c:.2f}%")

spread = max(mg_s, mg_b, mg_c) - min(mg_s, mg_b, mg_c)
print(f"\nDataset consistency:")
print(f"  Optima: SNe={mg_s:.5f}  BAO={mg_b:.5f}  CMB={mg_c:.5f}")
print(f"  Spread = {spread:.5f} ({spread/mu_g_micro*100:.1f}% of mu_g_micro)")
print(f"\nmu_g_micro within 1-sigma of joint: {'YES' if micro_in_1s else 'NO'}")
if not np.isnan(mg_lo):
    print(f"  1-sigma: [{mg_lo:.5f}, {mg_hi:.5f}]")

print(f"\nParameter economy:")
print(f"  {'Model':>18}  {'Params':>6}  {'SNe':>10}  {'BAO':>8}  {'CMB':>8}")
print(f"  {'LCDM':>18}  {'6':>6}  {'0.155 mag':>10}  {'~1%':>8}  {'<1%':>8}")
print(f"  {'UDT (joint)':>18}  {'1':>6}  {rj_s:.3f} mag  {rj_b:.1f}%  {rj_c:.2f}%")
print(f"  {'UDT (micro)':>18}  {'0':>6}  {rm_s:.3f} mag  {rm_b:.1f}%  {rm_c:.2f}%")

# === Plot ===
out_dir = os.path.join(os.path.dirname(__file__), 'outputs')
os.makedirs(out_dir, exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
labels = [f'SNe ({N_s})', f'BAO ({N_b})', f'CMB ({N_c})']
colors = ['#E69F00', '#56B4E9', '#009E73']
arrays = [c2_s, c2_b, c2_c]

for ax, arr, lab, col in zip(axes.flat[:3], arrays, labels, colors):
    ax.plot(mg_scan, arr - arr.min(), color=col, lw=2)
    ax.axvline(mg_j, color='k', ls='--', alpha=0.5, label=f'joint={mg_j:.4f}')
    ax.axvline(mu_g_micro, color='r', ls=':', alpha=0.5, label=f'micro={mu_g_micro:.4f}')
    ax.set_xlabel(r'$\mu_g$ (Gpc$^{-1}$)'); ax.set_ylabel(r'$\Delta\chi^2$')
    ax.set_title(lab); ax.legend(fontsize=8)
    ymax = min(float(arr.max() - arr.min()), 100)
    ax.set_ylim(0, max(ymax, 5))

ax = axes[1, 1]
ax.plot(mg_scan, c2_j - c2j_min, 'k', lw=2)
ax.axhline(1, color='gray', ls='--', alpha=0.5)
ax.axvline(mg_j, color='k', ls='--', alpha=0.5)
ax.axvline(mu_g_micro, color='r', ls=':', alpha=0.5, label=f'micro={mu_g_micro:.4f}')
if not np.isnan(mg_lo):
    ax.axvspan(mg_lo, mg_hi, alpha=0.15, color='blue', label=r'$1\sigma$')
ax.set_xlabel(r'$\mu_g$ (Gpc$^{-1}$)'); ax.set_ylabel(r'$\Delta\chi^2$')
ax.set_title('Joint'); ax.legend(fontsize=8); ax.set_ylim(0, 20)

plt.tight_layout()
outpath = os.path.join(out_dir, 'udt_joint_chi2.png')
plt.savefig(outpath, dpi=150)
print(f"\n  Saved: {outpath}")
plt.close()
print("Done.")
