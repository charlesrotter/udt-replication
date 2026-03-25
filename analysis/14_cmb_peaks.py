#!/usr/bin/env python3
"""14 -- CMB acoustic peak positions from WKB2 law.

SOURCE: lib/constants.py (MU_G, COSMO_K, COSMO_BETA, COSMO_GAMMA, PHI0, MU, R_STAR)
        lib/vacuum_phi.py (solve, compute_I2)
        lib/dirac_formT.py (find_eigenvalues)
GENERATES: data/generated/14_cmb_peaks.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  sqrt(-g) = c r^2 sin(theta)  (phi-independent)
  KG equation: box_g phi - mu^2 phi = -S

WKB2 peak law (second-order WKB):
  ell_n = ell_A * (nu_n + q1/nu_n + q2/nu_n^3)
  nu_n  = n + c_eff
  c_eff = arctan(eta) / pi
  eta   = 2/r* - 2*phi'(r*)   [cosmological boundary]

  q1, q2 are second-order WKB corrections from the potential shape.
  At lowest order: q1 = q2 = 0 (pure harmonic).

Zero-parameter ell_A from microphysics:
  ell_A = 2*pi*r_star*E_2 / (E_1*I_2)

Planck 2018 observed peaks:
  [220, 537.5, 810.8, 1120.9, 1444.2, 1776.0, 2081.0]

Completion class: B (sourced bounded -- polynomial is Class B; micro quantities Class A).
"""
import os
import sys
import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import (MU_G, COSMO_K, COSMO_BETA, COSMO_GAMMA,
                            PHI0, MU, R_STAR, N_GRID, R_MIN, PHI_CMB)
from lib.vacuum_phi import solve, compute_I2, extract_phi2
from lib.dirac_formT import find_eigenvalues
from lib.utils import save_results, gate_check

# === Planck 2018 observed peak positions ===
planck_peaks = np.array([220.0, 537.5, 810.8, 1120.9, 1444.2, 1776.0, 2081.0])
n_peaks = len(planck_peaks)

# === Cosmological polynomial ===

def phi_cosmo(r):
    """phi(r) = (3/2)*mu_g*r - cos(pi/5)*mu_g^2*r^2 + (2/3)*mu_g^3*r^3."""
    return COSMO_K * r + COSMO_BETA * r**2 + COSMO_GAMMA * r**3


def phi_cosmo_prime(r):
    """d(phi)/dr."""
    return COSMO_K + 2.0 * COSMO_BETA * r + 3.0 * COSMO_GAMMA * r**2


# === Find cosmological r* where phi(r*) = phi_CMB ===
from scipy.optimize import brentq

def find_r_star_cosmo():
    """Find r such that phi(r) = ln(1 + z_CMB)."""
    f = lambda r: phi_cosmo(r) - PHI_CMB
    return brentq(f, 0.1, 100.0, xtol=1e-12)


r_star_cosmo = find_r_star_cosmo()
phi_prime_rstar = phi_cosmo_prime(r_star_cosmo)

# === WKB2 peak law ===

def wkb2_peaks(ell_A, eta, q1=0.0, q2=0.0, n_max=7):
    """Compute peak positions from WKB2 law.

    ell_n = ell_A * (nu_n + q1/nu_n + q2/nu_n^3)
    nu_n = n + c_eff
    c_eff = arctan(eta) / pi
    """
    c_eff = np.arctan(eta) / np.pi
    peaks = []
    for n in range(1, n_max + 1):
        nu_n = n + c_eff
        ell_n = ell_A * (nu_n + q1 / nu_n + q2 / nu_n**3)
        peaks.append(ell_n)
    return np.array(peaks), c_eff


# === Compute eta from cosmological boundary ===
eta_cosmo = 2.0 / r_star_cosmo - 2.0 * phi_prime_rstar

# === Method 1: Anchor ell_A to peak 1 = 220 ===
# Solve: 220 = ell_A * (1 + c_eff) => ell_A = 220 / (1 + c_eff)
c_eff_0 = np.arctan(eta_cosmo) / np.pi
ell_A_anchored = planck_peaks[0] / (1.0 + c_eff_0)

peaks_anchored, _ = wkb2_peaks(ell_A_anchored, eta_cosmo, q1=0, q2=0, n_max=n_peaks)

# === Method 2: Scan mu_g near 0.248 for best fit ===
def rms_peaks(mu_g_trial):
    """Compute peak RMS for a given mu_g."""
    k = 1.5 * mu_g_trial
    beta = -np.cos(np.pi / 5) * mu_g_trial**2
    gamma = (2.0 / 3.0) * mu_g_trial**3

    def phi_trial(r):
        return k * r + beta * r**2 + gamma * r**3

    def phi_trial_prime(r):
        return k + 2.0 * beta * r + 3.0 * gamma * r**2

    # Find r* where phi = phi_CMB
    try:
        f = lambda r: phi_trial(r) - PHI_CMB
        if f(0.1) * f(100.0) > 0:
            return 1e10
        r_s = brentq(f, 0.1, 100.0, xtol=1e-12)
    except (ValueError, RuntimeError):
        return 1e10

    pp = phi_trial_prime(r_s)
    eta = 2.0 / r_s - 2.0 * pp
    c_eff = np.arctan(eta) / np.pi
    ell_A = planck_peaks[0] / (1.0 + c_eff)

    peaks, _ = wkb2_peaks(ell_A, eta, n_max=n_peaks)
    residuals = (peaks - planck_peaks) / planck_peaks
    return np.sqrt(np.mean(residuals**2))


# Scan
mu_g_scan = np.linspace(0.220, 0.280, 1000)
rms_scan = np.array([rms_peaks(mg) for mg in mu_g_scan])
best_idx = np.argmin(rms_scan)
mu_g_best = mu_g_scan[best_idx]

# Refine with scalar minimizer
result = minimize_scalar(rms_peaks, bounds=(0.230, 0.270), method='bounded')
mu_g_best = result.x
rms_best = result.fun

# Compute best-fit peaks
k_best = 1.5 * mu_g_best
beta_best = -np.cos(np.pi / 5) * mu_g_best**2
gamma_best = (2.0 / 3.0) * mu_g_best**3

def phi_best(r):
    return k_best * r + beta_best * r**2 + gamma_best * r**3

def phi_best_prime(r):
    return k_best + 2.0 * beta_best * r + 3.0 * gamma_best * r**2

r_star_best = brentq(lambda r: phi_best(r) - PHI_CMB, 0.1, 100.0, xtol=1e-12)
eta_best = 2.0 / r_star_best - 2.0 * phi_best_prime(r_star_best)
c_eff_best = np.arctan(eta_best) / np.pi
ell_A_best = planck_peaks[0] / (1.0 + c_eff_best)
peaks_best, _ = wkb2_peaks(ell_A_best, eta_best, n_max=n_peaks)

# === Method 3: Zero-parameter ell_A from microphysics ===
# ell_A = 2*pi*r_star*E_2 / (E_1*I_2)
# Need E_1, E_2 from Dirac solver and I_2 from vacuum
r, phi, J, phip, overflow = solve(PHI0, R_STAR, N_GRID, MU, r_min=R_MIN)
assert overflow is None, f"Overflow at r={overflow}"
e2phi = np.exp(np.clip(2 * phi, -400, 400))
phi2 = extract_phi2(r, phi, PHI0)
I_2 = compute_I2(r, phi)

evals_km1 = find_eigenvalues(-1, r, phip, e2phi, phi0=PHI0, phi2=phi2,
                              E_min=0.1, E_max=50.0, n_scan=20000, n_modes=3)
if len(evals_km1) >= 2:
    E1 = evals_km1[0]
    E2 = evals_km1[1]
    ell_A_micro = 2.0 * np.pi * R_STAR * E2 / (E1 * I_2)
    peaks_micro, _ = wkb2_peaks(ell_A_micro, eta_cosmo, n_max=n_peaks)
    rms_micro = np.sqrt(np.mean(((peaks_micro - planck_peaks) / planck_peaks)**2))
else:
    E1, E2 = np.nan, np.nan
    ell_A_micro = np.nan
    peaks_micro = np.full(n_peaks, np.nan)
    rms_micro = np.nan

# === Per-peak comparison ===
peak_comparison = []
for i in range(n_peaks):
    peak_comparison.append({
        'n': i + 1,
        'planck': planck_peaks[i],
        'wkb2_anchored': float(peaks_anchored[i]),
        'wkb2_best_mug': float(peaks_best[i]),
        'wkb2_micro': float(peaks_micro[i]) if np.isfinite(peaks_micro[i]) else None,
        'error_anchored_pct': float((peaks_anchored[i] - planck_peaks[i]) / planck_peaks[i] * 100),
        'error_best_pct': float((peaks_best[i] - planck_peaks[i]) / planck_peaks[i] * 100),
    })

# === RMS for anchored method ===
rms_anchored = np.sqrt(np.mean(((peaks_anchored - planck_peaks) / planck_peaks)**2))

# === Save results ===
results = {
    'planck_peaks': planck_peaks.tolist(),
    'cosmological_boundary': {
        'r_star_cosmo_Gpc': r_star_cosmo,
        'phi_prime_rstar': phi_prime_rstar,
        'eta': eta_cosmo,
        'c_eff': c_eff_0,
        'mu_g_locked': MU_G,
    },
    'method_1_anchored': {
        'ell_A': ell_A_anchored,
        'peaks': peaks_anchored.tolist(),
        'rms_fractional': rms_anchored,
        'rms_pct': rms_anchored * 100,
        'note': 'ell_A anchored to peak 1 = 220, q1=q2=0',
    },
    'method_2_best_mug': {
        'mu_g_best': mu_g_best,
        'ell_A': ell_A_best,
        'r_star_cosmo_Gpc': r_star_best,
        'eta': eta_best,
        'c_eff': c_eff_best,
        'peaks': peaks_best.tolist(),
        'rms_fractional': rms_best,
        'rms_pct': rms_best * 100,
        'note': f'mu_g scanned; best = {mu_g_best:.6f}',
    },
    'method_3_microphysics': {
        'E1': E1,
        'E2': E2,
        'E2_E1': E2 / E1 if E1 > 0 else None,
        'I_2': I_2,
        'ell_A_micro': ell_A_micro,
        'peaks': peaks_micro.tolist() if np.all(np.isfinite(peaks_micro)) else None,
        'rms_fractional': rms_micro,
        'rms_pct': rms_micro * 100 if np.isfinite(rms_micro) else None,
        'note': 'ell_A = 2*pi*r_star*E_2/(E_1*I_2), zero free parameters',
    },
    'peak_comparison': peak_comparison,
}
save_results('14_cmb_peaks.json', results)

# === Summary ===
print(f"CMB peaks: 7 peaks, Planck 2018 comparison")
print(f"  Anchored (ell_A={ell_A_anchored:.1f}): RMS = {rms_anchored*100:.2f}%")
print(f"  Best mu_g={mu_g_best:.5f}: RMS = {rms_best*100:.2f}%")
print(f"  Micro ell_A={ell_A_micro:.1f}: RMS = {rms_micro*100:.2f}%" if np.isfinite(rms_micro) else "  Micro: FAIL (eigenvalues not found)")
print(f"  eta = {eta_cosmo:.4f}, c_eff = {c_eff_0:.4f}")
