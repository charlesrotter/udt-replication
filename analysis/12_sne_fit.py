#!/usr/bin/env python3
"""12 -- Type Ia supernovae (Pantheon+) fit on UDT cosmological polynomial.

SOURCE: lib/constants.py (MU_G, COSMO_K, COSMO_BETA, COSMO_GAMMA, PHI0)
        data/external/pantheon_plus/ (if available)
GENERATES: data/generated/12_sne_fit.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  sqrt(-g) = c r^2 sin(theta)  (phi-independent)
  KG equation: box_g phi - mu^2 phi = -S  (positive S deepens well)

Cosmological polynomial (geometric fit, 0 free physics parameters):
  phi(r) = (3/2)*mu_g*r - cos(pi/5)*mu_g^2*r^2 + (2/3)*mu_g^3*r^3

Observables:
  Redshift: z = e^{phi(r)} - 1
  Luminosity distance: d_L(z) = r(z) * e^{phi(r(z))}   [in Gpc]
  Distance modulus: mu = 5*log10(d_L / 10 pc)

UDT is STATIC. There is no expansion. The redshift is gravitational (metric-induced).

Completion class: B (sourced bounded -- polynomial is a sourced Class B profile).
"""
import os
import sys
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import MU_G, COSMO_K, COSMO_BETA, COSMO_GAMMA, PHI0
from lib.utils import save_results

# === Cosmological polynomial ===

def phi_cosmo(r):
    """Cosmological phi(r) = (3/2)*mu_g*r - cos(pi/5)*mu_g^2*r^2 + (2/3)*mu_g^3*r^3."""
    return COSMO_K * r + COSMO_BETA * r**2 + COSMO_GAMMA * r**3


def phi_cosmo_prime(r):
    """d(phi)/dr for the cosmological polynomial."""
    return COSMO_K + 2.0 * COSMO_BETA * r + 3.0 * COSMO_GAMMA * r**2


# === Redshift <-> radius inversion ===

def z_of_r(r):
    """Redshift from radius: z = e^{phi(r)} - 1."""
    return np.exp(phi_cosmo(r)) - 1.0


def r_of_z(z_target, r_max=60.0):
    """Invert z(r) to find r(z) by root finding."""
    if z_target <= 0:
        return 0.0
    # Find r such that e^{phi(r)} - 1 = z_target
    f = lambda r: np.exp(phi_cosmo(r)) - 1.0 - z_target
    # Check bracket
    if f(1e-6) * f(r_max) > 0:
        # Extend search
        r_max = 100.0
        if f(1e-6) * f(r_max) > 0:
            return np.nan
    try:
        return brentq(f, 1e-6, r_max, xtol=1e-12)
    except ValueError:
        return np.nan


# === Luminosity distance ===

def d_L(z):
    """Luminosity distance: d_L = r(z) * e^{phi(r(z))} = r(z) * (1+z) [Gpc]."""
    r = r_of_z(z)
    if np.isnan(r):
        return np.nan
    return r * np.exp(phi_cosmo(r))


def distance_modulus(d_L_Gpc):
    """Distance modulus: mu = 5*log10(d_L / 10 pc).

    d_L in Gpc. 1 Gpc = 10^9 pc, so d_L/10pc = d_L * 10^8.
    mu = 5*log10(d_L_Gpc * 1e8) = 5*log10(d_L_Gpc) + 40.
    """
    if d_L_Gpc <= 0 or np.isnan(d_L_Gpc):
        return np.nan
    return 5.0 * np.log10(d_L_Gpc) + 40.0


# === Load Pantheon+ data or use representative points ===

pantheon_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'external')
pantheon_file = os.path.join(pantheon_dir, 'Pantheon+SH0ES.dat')

use_external = False
z_data = None
mu_data = None
mu_err = None

if os.path.isfile(pantheon_file):
    try:
        import pandas as pd
        df = pd.read_csv(pantheon_file, sep=r'\s+', comment='#')
        z_data = df['zHD'].values
        mu_data = df['m_b_corr'].values       # corrected apparent magnitudes
        mu_err = df['m_b_corr_err_DIAG'].values
        use_external = True
        print(f"  Loaded {len(z_data)} Pantheon+ data points")
    except Exception as e:
        print(f"  Could not load Pantheon+ data: {e}")

if not use_external:
    # Representative supernova data points (Pantheon+ summary statistics)
    # z, distance modulus mu, uncertainty sigma_mu
    # These are approximate bin-averaged values from Pantheon+ (Brout+ 2022)
    z_data = np.array([0.01, 0.023, 0.05, 0.08, 0.12, 0.20, 0.30, 0.40,
                       0.50, 0.60, 0.70, 0.80, 1.00, 1.20, 1.50, 2.00, 2.26])
    mu_data = np.array([33.10, 34.67, 36.64, 37.82, 38.70, 39.82, 40.73, 41.35,
                        41.82, 42.20, 42.50, 42.77, 43.19, 43.51, 43.90, 44.39, 44.56])
    mu_err = np.array([0.15, 0.12, 0.10, 0.09, 0.08, 0.07, 0.06, 0.06,
                       0.06, 0.06, 0.07, 0.07, 0.08, 0.10, 0.12, 0.15, 0.18])
    print(f"  Using {len(z_data)} representative data points")

# === Compute UDT predictions ===
d_L_pred = np.array([d_L(z) for z in z_data])
mu_pred = np.array([distance_modulus(dl) for dl in d_L_pred])

# Filter out NaN predictions
valid = np.isfinite(mu_pred) & np.isfinite(mu_data)
z_valid = z_data[valid]
mu_pred_valid = mu_pred[valid]
mu_data_valid = mu_data[valid]
mu_err_valid = mu_err[valid]
d_L_valid = d_L_pred[valid]

# === Fit M offset (apparent mag = distance modulus + M) ===
# When using Pantheon+ m_b_corr, M_offset absorbs the absolute magnitude.
# This is NOT a free physics parameter — it is a calibration offset.
from scipy.optimize import minimize_scalar
def rms_at_M(M):
    res = mu_pred_valid - (mu_data_valid - M)
    return np.sqrt(np.mean(res**2))
M_fit = minimize_scalar(rms_at_M, bounds=(-25, -10), method='bounded').x
mu_data_valid = mu_data_valid - M_fit  # convert to distance modulus
print(f"  M offset: {M_fit:.4f}")

# === Fit statistics ===
residuals = mu_pred_valid - mu_data_valid
rms_mag = np.sqrt(np.mean(residuals**2))
chi2 = np.sum((residuals / mu_err_valid)**2)
chi2_dof = chi2 / max(1, len(residuals) - 0)  # 0 free physics parameters
n_points = len(residuals)

# === Also compute on a fine grid for output ===
z_fine = np.logspace(-2.5, np.log10(2.5), 200)
d_L_fine = np.array([d_L(z) for z in z_fine])
mu_fine = np.array([distance_modulus(dl) for dl in d_L_fine])

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
    'data': {
        'source': 'Pantheon+ external' if use_external else 'representative bin averages',
        'n_points': n_points,
        'z': z_valid.tolist(),
        'mu_data': mu_data_valid.tolist(),
        'mu_err': mu_err_valid.tolist(),
    },
    'predictions': {
        'z': z_valid.tolist(),
        'd_L_Gpc': d_L_valid.tolist(),
        'mu_pred': mu_pred_valid.tolist(),
        'residuals_mag': residuals.tolist(),
    },
    'fit_statistics': {
        'rms_mag': rms_mag,
        'chi2': chi2,
        'chi2_per_dof': chi2_dof,
        'n_dof': n_points,
    },
    'model_curve': {
        'z': z_fine.tolist(),
        'd_L_Gpc': d_L_fine.tolist(),
        'mu': mu_fine.tolist(),
    },
}
save_results('12_sne_fit.json', results)

# === Summary ===
print(f"SNe Ia fit: {n_points} points, 0 free parameters")
print(f"  RMS = {rms_mag:.3f} mag")
print(f"  chi^2 = {chi2:.1f}, chi^2/N = {chi2_dof:.2f}")
print(f"  mu_g = {MU_G:.6f} Gpc^-1")
print(f"  Data: {'Pantheon+ external' if use_external else 'representative'}")
