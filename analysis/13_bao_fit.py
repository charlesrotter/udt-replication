#!/usr/bin/env python3
"""13 -- BAO distance scale fit on UDT cosmological polynomial.

SOURCE: lib/constants.py (MU_G, COSMO_K, COSMO_BETA, COSMO_GAMMA)
        data/external/boss_bao/, data/external/desi_bao/ (if available)
GENERATES: data/generated/13_bao_fit.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  sqrt(-g) = c r^2 sin(theta)  (phi-independent)
  KG equation: box_g phi - mu^2 phi = -S  (positive S deepens well)

Cosmological polynomial (geometric, 0 free physics parameters):
  phi(r) = (3/2)*mu_g*r - cos(pi/5)*mu_g^2*r^2 + (2/3)*mu_g^3*r^3

BAO observable:
  D_V(z) / r_d  where  D_V = [z * D_M(z)^2 * D_H(z)]^{1/3}
  D_M(z) = r(z)                  -- comoving distance in UDT (static metric)
  D_H(z) = c / H(z)              -- Hubble distance
  H(z) = c * phi'(r) / (1+z)     -- gravitational "expansion rate" from metric gradient

UDT note: BAO is a kaleidoscope artifact, NOT acoustic waves. r_d is a geometric
coherence ruler from the metric, not a sound horizon. We do NOT use the LCDM
r_d = 147 Mpc prior.

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


def z_of_r(r):
    """z = e^{phi(r)} - 1."""
    return np.exp(phi_cosmo(r)) - 1.0


def r_of_z(z_target, r_max=60.0):
    """Invert z(r) to find r(z)."""
    if z_target <= 0:
        return 0.0
    f = lambda r: np.exp(phi_cosmo(r)) - 1.0 - z_target
    try:
        return brentq(f, 1e-6, r_max, xtol=1e-12)
    except ValueError:
        return np.nan


# === BAO observables on UDT metric ===

def D_M(z):
    """Comoving distance D_M(z) = r(z) in Gpc."""
    return r_of_z(z)


def H_of_z(z):
    """H(z) = c * phi'(r(z)) / (1+z) in km/s/Mpc equivalent.

    phi'(r) has units Gpc^{-1}. H = c * phi'/(1+z) in Gpc^{-1} * (km/s).
    Convert: 1 Gpc = 1e3 Mpc, so H in km/s/Mpc = c_km_s * phi'_Gpc / ((1+z) * 1e3).
    """
    r = r_of_z(z)
    if np.isnan(r):
        return np.nan
    pp = phi_cosmo_prime(r)
    c_km_s = C_LIGHT * 1e-3  # m/s -> km/s
    return c_km_s * pp / ((1.0 + z) * 1e3)  # km/s/Mpc


def D_H(z):
    """Hubble distance D_H(z) = c/H(z) in Gpc.

    D_H = c / H(z). With H in km/s/Mpc and c in km/s:
    D_H_Mpc = c_km_s / H_km_s_Mpc, then D_H_Gpc = D_H_Mpc / 1e3.
    Equivalently: D_H_Gpc = (1+z) / phi'(r).
    """
    r = r_of_z(z)
    if np.isnan(r):
        return np.nan
    pp = phi_cosmo_prime(r)
    if abs(pp) < 1e-30:
        return np.nan
    return (1.0 + z) / pp  # Gpc


def D_V(z):
    """Volume-averaged distance: D_V = [z * D_M^2 * D_H]^{1/3} in Gpc."""
    dm = D_M(z)
    dh = D_H(z)
    if np.isnan(dm) or np.isnan(dh) or dm <= 0 or dh <= 0:
        return np.nan
    return (z * dm**2 * dh)**(1.0 / 3.0)


# === UDT geometric r_d ===
# r_d is NOT the LCDM sound horizon. It is a geometric coherence scale.
# r_d = 1000 / sqrt(2*pi*|k|*|beta|*(1+z_d))
# where k, beta are polynomial coefficients and z_d ~ 1060 (drag epoch)
z_d = 1060.0  # drag epoch redshift (approximate, from CMB)
r_d_udt = 1000.0 / np.sqrt(2.0 * np.pi * abs(COSMO_K) * abs(COSMO_BETA) * (1.0 + z_d))
# Convert to Mpc (r_d_udt is in Gpc since mu_g is in Gpc^{-1})
r_d_udt_Mpc = r_d_udt * 1e3

# === BAO data points ===
# BOSS DR12 + DESI Y1 consensus
# Format: z_eff, D_V/r_d (observed), sigma(D_V/r_d)
# r_d is the sound horizon in LCDM context; we compare D_V/r_d_udt
bao_data = [
    # BOSS DR12 (Alam+ 2017)
    {'z': 0.38, 'DV_rd': 10.27, 'sigma': 0.15, 'source': 'BOSS DR12'},
    {'z': 0.51, 'DV_rd': 13.38, 'sigma': 0.18, 'source': 'BOSS DR12'},
    {'z': 0.70, 'DV_rd': 17.86, 'sigma': 0.33, 'source': 'BOSS DR12'},
    # DESI Y1 (DESI Collaboration 2024)
    {'z': 1.48, 'DV_rd': 30.69, 'sigma': 0.80, 'source': 'DESI Y1 (Lya)'},
    {'z': 2.33, 'DV_rd': 37.60, 'sigma': 1.20, 'source': 'DESI Y1 (Lya)'},
]

# === Compute UDT predictions ===
predictions = []
for pt in bao_data:
    z = pt['z']
    dv = D_V(z)
    dm = D_M(z)
    dh = D_H(z)
    hz = H_of_z(z)

    if np.isnan(dv):
        dv_rd = np.nan
    else:
        dv_rd = dv / r_d_udt  # both in Gpc

    predictions.append({
        'z': z,
        'D_M_Gpc': dm,
        'D_H_Gpc': dh,
        'D_V_Gpc': dv,
        'H_km_s_Mpc': hz,
        'DV_rd_pred': dv_rd,
        'DV_rd_data': pt['DV_rd'],
        'sigma': pt['sigma'],
        'source': pt['source'],
    })

# === Fit statistics ===
z_arr = np.array([p['z'] for p in predictions])
pred_arr = np.array([p['DV_rd_pred'] for p in predictions])
data_arr = np.array([p['DV_rd_data'] for p in predictions])
sig_arr = np.array([p['sigma'] for p in predictions])

valid = np.isfinite(pred_arr)
residuals = pred_arr[valid] - data_arr[valid]
rms = np.sqrt(np.mean(residuals**2))
chi2 = np.sum((residuals / sig_arr[valid])**2)
n_valid = np.sum(valid)

# === Save results ===
results = {
    'parameters': {
        'mu_g': MU_G,
        'cosmo_k': COSMO_K,
        'cosmo_beta': COSMO_BETA,
        'cosmo_gamma': COSMO_GAMMA,
        'n_free_params': 0,
        'polynomial': 'phi(r) = (3/2)*mu_g*r - cos(pi/5)*mu_g^2*r^2 + (2/3)*mu_g^3*r^3',
    },
    'r_d_udt': {
        'value_Gpc': r_d_udt,
        'value_Mpc': r_d_udt_Mpc,
        'formula': '1000 / sqrt(2*pi*|k|*|beta|*(1+z_d))',
        'z_d': z_d,
        'note': 'geometric coherence ruler, NOT LCDM sound horizon',
    },
    'predictions': predictions,
    'fit_statistics': {
        'rms_DV_rd': rms,
        'chi2': chi2,
        'chi2_per_point': chi2 / max(1, n_valid),
        'n_points': int(n_valid),
        'n_free_params': 0,
    },
}
save_results('13_bao_fit.json', results)

# === Summary ===
print(f"BAO fit: {int(n_valid)} points, 0 free parameters")
print(f"  r_d(UDT) = {r_d_udt_Mpc:.2f} Mpc (geometric)")
print(f"  RMS(D_V/r_d) = {rms:.2f}")
print(f"  chi^2 = {chi2:.1f}, chi^2/N = {chi2 / max(1, n_valid):.2f}")
print(f"  mu_g = {MU_G:.6f} Gpc^-1")
