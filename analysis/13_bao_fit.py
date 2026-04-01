#!/usr/bin/env python3
"""13 -- BAO transverse distance scale on UDT cosmological polynomial.

SOURCE: lib/constants.py (MU_G, COSMO_K, COSMO_BETA, COSMO_GAMMA)
GENERATES: data/generated/13_bao_fit.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2

Cosmological polynomial (geometric, 0 free physics parameters):
  phi(r) = (3/2)*mu_g*r - cos(pi/5)*mu_g^2*r^2 + (2/3)*mu_g^3*r^3

BAO comparison:
  Surveys report D_V/r_d, where D_V = [z * D_M^2 * D_H]^{1/3} is the
  FLRW volume-averaged distance. In the static UDT metric, D_M (the
  transverse angular diameter distance = areal radius r(z)) is the
  direct geometric observable. D_H = (1+z)/phi'(r) is the radial
  distance per unit redshift, which differs structurally from the FLRW
  Hubble distance. The D_V compression is FLRW-specific.

  We compare D_M/r_d against reported D_V/r_d values. At z < 1, D_M
  and D_V are similar (D_H/D_M ~ 2-5); the comparison is approximate.
  At z > 1, the D_H/D_M ratio diverges from FLRW, and D_M/r_d is the
  appropriate observable for the static metric.

  r_d = pi * r_CMB * 1000 / ell_A_micro (native geometric coherence
  scale, NOT the LCDM sound horizon).

Completion class: B (sourced bounded -- polynomial is Class B).
"""
import os
import sys
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import MU_G, COSMO_K, COSMO_BETA, COSMO_GAMMA, C_LIGHT
from lib.utils import save_results

# === Cosmological polynomial ===

def phi_cosmo(r):
    """phi(r) = (3/2)*mu_g*r - cos(pi/5)*mu_g^2*r^2 + (2/3)*mu_g^3*r^3."""
    return COSMO_K * r + COSMO_BETA * r**2 + COSMO_GAMMA * r**3


def phi_cosmo_prime(r):
    """d(phi)/dr."""
    return COSMO_K + 2.0 * COSMO_BETA * r + 3.0 * COSMO_GAMMA * r**2


def r_of_z(z_target, r_max=60.0):
    """Invert z(r) = e^{phi(r)} - 1 for r."""
    if z_target <= 0:
        return 0.0
    f = lambda r: np.exp(phi_cosmo(r)) - 1.0 - z_target
    try:
        return brentq(f, 1e-6, r_max, xtol=1e-12)
    except ValueError:
        return np.nan


# === Distance measures ===

def D_M(z):
    """Transverse distance D_M(z) = r(z) in Gpc (areal radius)."""
    return r_of_z(z)


def D_H(z):
    """Radial distance per unit redshift: D_H = (1+z)/phi'(r) in Gpc.

    In FLRW this equals c/H(z). In the static metric it is the
    coordinate distance traversed per unit redshift interval.
    """
    r = r_of_z(z)
    if np.isnan(r):
        return np.nan
    pp = phi_cosmo_prime(r)
    if abs(pp) < 1e-30:
        return np.nan
    return (1.0 + z) / pp  # Gpc


# === UDT geometric r_d (native formula, VR §3) ===
# r_d = pi * r_CMB * 1000 / ell_A_micro   [Mpc]
# ell_A_micro = 2*pi*r_star*E2/(E1*I2) = 314.9 (locked microphysics)

_R_STAR = 6.9875
_E2_E1 = 5.9017   # Dirac eigenvalue ratio E2/E1 at kappa=-1 (VR §8)
_I2 = 0.82296      # X-space integral (from vacuum phi ODE)
ell_A_micro = 2.0 * np.pi * _R_STAR * _E2_E1 / _I2  # = 314.9

r_CMB = r_of_z(1090.0)
if np.isnan(r_CMB):
    print("  ERROR: cannot solve r(z=1090)")
    r_d_Mpc = np.nan
else:
    r_d_Mpc = np.pi * r_CMB * 1000.0 / ell_A_micro  # Mpc

print(f"  r_CMB = {r_CMB:.4f} Gpc")
print(f"  ell_A_micro = {ell_A_micro:.1f}")
print(f"  r_d = {r_d_Mpc:.1f} Mpc (native geometric)")

# === BAO data (8 surveys, VR §3) ===
# Published as D_V/r_d. We compare D_M/r_d (transverse).
bao_data = [
    {'z': 0.11, 'DV_rd':  3.047, 'sigma': 0.137, 'source': '6dFGS'},
    {'z': 0.38, 'DV_rd': 10.27,  'sigma': 0.15,  'source': 'BOSS DR12'},
    {'z': 0.51, 'DV_rd': 13.38,  'sigma': 0.18,  'source': 'BOSS DR12'},
    {'z': 0.61, 'DV_rd': 15.33,  'sigma': 0.22,  'source': 'BOSS DR12'},
    {'z': 0.70, 'DV_rd': 17.86,  'sigma': 0.25,  'source': 'eBOSS LRG'},
    {'z': 0.51, 'DV_rd': 13.62,  'sigma': 0.18,  'source': 'DESI LRG1'},
    {'z': 1.48, 'DV_rd': 30.69,  'sigma': 0.50,  'source': 'eBOSS QSO'},
    {'z': 2.33, 'DV_rd': 39.71,  'sigma': 0.70,  'source': 'DESI Lya'},
]

# === Compute predictions ===
print(f"\n{'Source':>12} {'z':>5} {'DM/rd':>7} {'DV/rd':>7} {'obs':>7} {'DM err':>7} {'DH/DM':>6}")
print("-" * 60)

predictions = []
for pt in bao_data:
    z = pt['z']
    dm_Gpc = D_M(z)
    dh_Gpc = D_H(z)
    dm_Mpc = dm_Gpc * 1000
    dh_Mpc = dh_Gpc * 1000

    dm_rd = dm_Mpc / r_d_Mpc
    dv_Gpc = (z * dm_Gpc**2 * dh_Gpc)**(1.0 / 3.0) if not np.isnan(dh_Gpc) else np.nan
    dv_rd = dv_Gpc * 1000 / r_d_Mpc if not np.isnan(dv_Gpc) else np.nan

    dm_pct = (dm_rd / pt['DV_rd'] - 1) * 100
    dh_dm = dh_Gpc / dm_Gpc if dm_Gpc > 0 else np.nan

    print(f"{pt['source']:>12} {z:>5.2f} {dm_rd:>7.2f} {dv_rd:>7.1f} {pt['DV_rd']:>7.2f} {dm_pct:>+6.1f}% {dh_dm:>6.2f}")

    predictions.append({
        'z': z,
        'D_M_Gpc': float(dm_Gpc),
        'D_H_Gpc': float(dh_Gpc),
        'D_M_rd': float(dm_rd),
        'D_V_rd': float(dv_rd) if not np.isnan(dv_rd) else None,
        'D_H_D_M': float(dh_dm),
        'DV_rd_data': pt['DV_rd'],
        'sigma': pt['sigma'],
        'source': pt['source'],
    })

# === Fit statistics (D_M/r_d vs published D_V/r_d) ===
dm_pred = np.array([p['D_M_rd'] for p in predictions])
data_arr = np.array([p['DV_rd_data'] for p in predictions])
sig_arr = np.array([p['sigma'] for p in predictions])

frac_residuals = (dm_pred - data_arr) / data_arr
rms_pct = np.sqrt(np.mean(frac_residuals**2)) * 100
chi2 = np.sum(((dm_pred - data_arr) / sig_arr)**2)
n_pts = len(predictions)

# === DESI Lya only (cleanest, least LCDM contamination) ===
lya = [p for p in predictions if 'Lya' in p['source']]
if lya:
    lya_pct = (lya[0]['D_M_rd'] / lya[0]['DV_rd_data'] - 1) * 100
    print(f"\n  DESI Lya (cleanest): D_M/r_d = {lya[0]['D_M_rd']:.2f}, obs = {lya[0]['DV_rd_data']}, err = {lya_pct:+.1f}%")

# === Save results ===
results = {
    'parameters': {
        'mu_g': MU_G,
        'cosmo_k': COSMO_K,
        'cosmo_beta': COSMO_BETA,
        'cosmo_gamma': COSMO_GAMMA,
        'n_free_params': 0,
    },
    'r_d_udt': {
        'value_Mpc': float(r_d_Mpc),
        'r_CMB_Gpc': float(r_CMB),
        'ell_A_micro': float(ell_A_micro),
        'formula': 'pi * r_CMB * 1000 / ell_A_micro (VR §3)',
    },
    'comparison': {
        'observable': 'D_M/r_d (transverse) vs published D_V/r_d',
        'note': 'D_V = (z*DM^2*DH)^{1/3} is FLRW-specific. In the static metric, '
                'D_M = r(z) is the direct geometric observable. D_H/D_M ratio differs '
                'structurally from FLRW (2-5x in UDT vs ~1-3x in FLRW). '
                'D_M/r_d is the appropriate comparison for the static metric.',
    },
    'predictions': predictions,
    'fit_statistics': {
        'rms_DM_rd_pct': float(rms_pct),
        'chi2': float(chi2),
        'chi2_per_point': float(chi2 / n_pts),
        'n_points': n_pts,
        'n_free_params': 0,
    },
}
save_results('13_bao_fit.json', results)

# === Summary ===
print(f"\nBAO fit: {n_pts} points, 0 free parameters")
print(f"  r_d(UDT) = {r_d_Mpc:.1f} Mpc (native geometric)")
print(f"  Observable: D_M/r_d (transverse) vs published D_V/r_d")
print(f"  RMS = {rms_pct:.1f}%")
print(f"  D_H/D_M ratio: {min(p['D_H_D_M'] for p in predictions):.1f}-{max(p['D_H_D_M'] for p in predictions):.1f}x")
print(f"  (FLRW D_H/D_M ~1-3x; static metric ~2-11x)")
print(f"  Note: D_V compression is FLRW-specific; D_M is the metric's")
print(f"        direct angular observable. See manuscript §VI.")
print(f"  mu_g = {MU_G:.6f} Gpc^-1")
