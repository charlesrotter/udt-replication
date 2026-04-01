#!/usr/bin/env python3
"""14 -- CMB acoustic peak positions from WKB2 law.

SOURCE: lib/constants.py (MU_G, COSMO_K, COSMO_BETA, COSMO_GAMMA, PHI0, MU, R_STAR)
        lib/vacuum_phi.py (solve, compute_I2)
        lib/dirac_formT.py (find_eigenvalues)
GENERATES: data/generated/14_cmb_peaks.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2

WKB2 peak law (second-order WKB on the sourced phi profile):
  ell_n = ell_A * (nu_n + q1/nu_n + q2/nu_n^3)
  nu_n  = n + c_eff
  c_eff = arctan(eta) / pi
  eta   = 2/r*_cosmo - 2*phi'(r*_cosmo)   [cosmological boundary]

  WKB2 corrections from the angular sector:
    q1 = (2j+1) / (2l+1)^2 = 2/9
    q2 = -1 / ((2j+1)^2 * (2|kmax|-1)) = -1/20

  ell_A determined by anchoring peak 1 to observed value (220):
    ell_A = 220 / F1,  F1 = nu_1 + q1/nu_1 + q2/nu_1^3

  This is NOT a fit: eta comes from the polynomial at derived mu_g,
  q1 and q2 are Diophantine quantum numbers, and peak 1 sets the
  angular scale (1 anchored observable, 0 fitted physics parameters).

Planck 2018 observed peaks:
  [220, 537.5, 810.8, 1120.9, 1444.2, 1776.0, 2081.0]

Completion class: B (sourced bounded -- polynomial is Class B; micro quantities Class A).
"""
import os
import sys
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import (MU_G, COSMO_K, COSMO_BETA, COSMO_GAMMA,
                            MULT_2J1, MULT_2L1, MULT_2KM1, MULT_2KP1,
                            PHI0, MU, R_STAR, N_GRID, R_MIN, PHI_CMB)
from lib.vacuum_phi import solve, compute_I2, extract_phi2
from lib.dirac_formT import find_eigenvalues
from lib.utils import save_results

# === Planck 2018 observed peak positions ===
planck_peaks = np.array([220.0, 537.5, 810.8, 1120.9, 1444.2, 1776.0, 2081.0])
n_peaks = len(planck_peaks)

# === WKB2 corrections from Diophantine quantum numbers ===
q1 = MULT_2J1 / MULT_2L1**2          # (2j+1)/(2l+1)^2 = 2/9
q2 = -1.0 / (MULT_2J1**2 * MULT_2KM1)  # -1/((2j+1)^2*(2kmax-1)) = -1/20

print(f"WKB2 corrections from angular quantum numbers:")
print(f"  q1 = (2j+1)/(2l+1)^2 = {int(MULT_2J1)}/{int(MULT_2L1**2)} = {q1:.6f}")
print(f"  q2 = -1/((2j+1)^2*(2kmax-1)) = -1/{int(MULT_2J1**2 * MULT_2KM1)} = {q2:.6f}")

# === Cosmological polynomial ===

def phi_cosmo(r):
    return COSMO_K * r + COSMO_BETA * r**2 + COSMO_GAMMA * r**3

def phi_cosmo_prime(r):
    return COSMO_K + 2.0 * COSMO_BETA * r + 3.0 * COSMO_GAMMA * r**2

# === Cosmological boundary ===
r_star_cosmo = brentq(lambda r: phi_cosmo(r) - PHI_CMB, 0.1, 100.0, xtol=1e-12)
phi_prime_rstar = phi_cosmo_prime(r_star_cosmo)
eta_cosmo = 2.0 / r_star_cosmo - 2.0 * phi_prime_rstar
c_eff = np.arctan(eta_cosmo) / np.pi

print(f"\nCosmological boundary:")
print(f"  r*_cosmo = {r_star_cosmo:.4f} Gpc")
print(f"  phi'(r*) = {phi_prime_rstar:.6f}")
print(f"  eta = 2/r* - 2*phi' = {eta_cosmo:.4f}")
print(f"  c_eff = arctan(eta)/pi = {c_eff:.4f}")

# === WKB2 peak law ===

def wkb2_peaks(ell_A, eta, q1_val, q2_val, n_max=7):
    """Compute peak positions: ell_n = ell_A * (nu_n + q1/nu_n + q2/nu_n^3)."""
    c = np.arctan(eta) / np.pi
    peaks = []
    for n in range(1, n_max + 1):
        nu = n + c
        ell = ell_A * (nu + q1_val / nu + q2_val / nu**3)
        peaks.append(ell)
    return np.array(peaks), c

# === Method 1: WKB2 with Diophantine q, peak-1 anchored ===
# ell_A = 220 / F1 where F1 = nu_1 + q1/nu_1 + q2/nu_1^3
nu_1 = 1.0 + c_eff
F1 = nu_1 + q1 / nu_1 + q2 / nu_1**3
ell_A_wkb2 = planck_peaks[0] / F1

peaks_wkb2, _ = wkb2_peaks(ell_A_wkb2, eta_cosmo, q1, q2, n_max=n_peaks)
errs_wkb2 = (peaks_wkb2 - planck_peaks) / planck_peaks * 100
rms_wkb2 = np.sqrt(np.mean((errs_wkb2 / 100)**2)) * 100

print(f"\n--- Method 1: WKB2 + Diophantine q (peak 1 anchored) ---")
print(f"  ell_A = {ell_A_wkb2:.1f}")
print(f"  0 fitted physics parameters; 1 anchored scale (peak 1)")
for i in range(n_peaks):
    print(f"  peak {i+1}: pred={peaks_wkb2[i]:.1f}  obs={planck_peaks[i]:.1f}  err={errs_wkb2[i]:+.2f}%")
print(f"  RMS = {rms_wkb2:.2f}%")

# === Method 2: Simple WKB0 (q1=q2=0), peak-1 anchored ===
ell_A_wkb0 = planck_peaks[0] / (1.0 + c_eff)
peaks_wkb0, _ = wkb2_peaks(ell_A_wkb0, eta_cosmo, 0, 0, n_max=n_peaks)
rms_wkb0 = np.sqrt(np.mean(((peaks_wkb0 - planck_peaks) / planck_peaks)**2)) * 100

print(f"\n--- Method 2: WKB0 (q1=q2=0, peak 1 anchored) ---")
print(f"  ell_A = {ell_A_wkb0:.1f}")
print(f"  RMS = {rms_wkb0:.2f}%")

# === Method 3: Micro ell_A, no anchoring ===
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
    peaks_micro, _ = wkb2_peaks(ell_A_micro, eta_cosmo, q1, q2, n_max=n_peaks)
    rms_micro = np.sqrt(np.mean(((peaks_micro - planck_peaks) / planck_peaks)**2)) * 100
else:
    E1, E2 = np.nan, np.nan
    ell_A_micro = np.nan
    peaks_micro = np.full(n_peaks, np.nan)
    rms_micro = np.nan

print(f"\n--- Method 3: Micro ell_A + WKB2 (not anchored) ---")
print(f"  E1 = {E1:.6f}, E2 = {E2:.6f}, E2/E1 = {E2/E1:.4f}")
print(f"  I_2 = {I_2:.6f}")
print(f"  ell_A_micro = {ell_A_micro:.1f}")
if np.isfinite(rms_micro):
    for i in range(n_peaks):
        err = (peaks_micro[i] - planck_peaks[i]) / planck_peaks[i] * 100
        print(f"  peak {i+1}: pred={peaks_micro[i]:.1f}  obs={planck_peaks[i]:.1f}  err={err:+.2f}%")
    print(f"  RMS = {rms_micro:.2f}%")
else:
    print(f"  FAIL (eigenvalues not found)")

# === Save results ===
results = {
    'planck_peaks': planck_peaks.tolist(),
    'quantum_numbers': {
        'q1': float(q1), 'q1_formula': '(2j+1)/(2l+1)^2 = 2/9',
        'q2': float(q2), 'q2_formula': '-1/((2j+1)^2*(2kmax-1)) = -1/20',
    },
    'cosmological_boundary': {
        'r_star_cosmo_Gpc': float(r_star_cosmo),
        'phi_prime_rstar': float(phi_prime_rstar),
        'eta': float(eta_cosmo),
        'c_eff': float(c_eff),
        'mu_g': float(MU_G),
    },
    'method_1_wkb2_diophantine': {
        'ell_A': float(ell_A_wkb2),
        'peaks': peaks_wkb2.tolist(),
        'rms_pct': float(rms_wkb2),
        'note': '0 fitted physics params; q1,q2 from quantum numbers; peak 1 anchored',
    },
    'method_2_wkb0': {
        'ell_A': float(ell_A_wkb0),
        'rms_pct': float(rms_wkb0),
        'note': 'q1=q2=0, peak 1 anchored',
    },
    'method_3_micro': {
        'E1': float(E1), 'E2': float(E2),
        'E2_E1': float(E2 / E1) if E1 > 0 else None,
        'I_2': float(I_2),
        'ell_A_micro': float(ell_A_micro),
        'rms_pct': float(rms_micro) if np.isfinite(rms_micro) else None,
        'note': 'ell_A from microphysics + WKB2 Diophantine q, not anchored',
    },
}
save_results('14_cmb_peaks.json', results)

# === Summary ===
print(f"\nCMB peaks: {n_peaks} peaks, Planck 2018 comparison")
print(f"  WKB2 + Diophantine q (peak 1 anchored): {rms_wkb2:.2f}% RMS")
print(f"  WKB0 (peak 1 anchored): {rms_wkb0:.2f}% RMS")
if np.isfinite(rms_micro):
    print(f"  Micro ell_A + WKB2 (not anchored): {rms_micro:.2f}% RMS")
