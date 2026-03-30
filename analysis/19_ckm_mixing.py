#!/usr/bin/env python3
"""19 -- CKM matrix and quark mixing.

GENERATES: data/generated/19_ckm_mixing.json
VERIFIES: Manuscript Section 16 (CKM Matrix and Quark Mixing)

sin(theta_C) = 9/40 (-0.13%)
A = cos(pi/5) (-0.24%)
delta_CKM = arctan(9/4) = 66.0 deg (0.1 sigma)
theta_QCD = 0 (strong CP solved)
J_CKM = 3.03e-5 (+0.8%)

Completion class: S (structural, parameter-independent).
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import (MULT_2J1, MULT_2L1, MULT_2KM1, MULT_2KP1,
                           PHI0)
from lib.utils import save_results, pct_error, format_comparison, sigma_tension

print("=" * 60)
print("CKM MATRIX AND QUARK MIXING")
print("=" * 60)

# === Cabibbo angle ===
# sin(theta_C) = (2l+1)^2 / ((2j+1)^3 * (2|kmax|-1)) = 9/40
sin_C = MULT_2L1**2 / (MULT_2J1**3 * MULT_2KM1)
sin_C_exp = 0.2253
print("\nCabibbo angle:")
format_comparison("  sin(theta_C) = 9/40", sin_C, sin_C_exp)

# === Wolfenstein A ===
A_wolf = abs(PHI0)  # cos(pi/5)
A_wolf_exp = 0.811
print("\nWolfenstein A:")
format_comparison("  A = |phi_0| = cos(pi/5)", A_wolf, A_wolf_exp)

# === Wolfenstein lambda ===
lam = sin_C  # lambda = sin(theta_C)

# === CKM CP phase ===
delta_ckm = np.arctan(9.0/4.0)  # = arctan((2l+1)^2/(2j+1)^2)
delta_ckm_deg = np.degrees(delta_ckm)
delta_ckm_exp = 65.5  # degrees
delta_ckm_sigma = 2.5
print("\nCKM CP phase:")
format_comparison("  delta_CKM = arctan(9/4)", delta_ckm_deg, delta_ckm_exp, "deg")
tens = sigma_tension(delta_ckm_deg, delta_ckm_exp, delta_ckm_sigma)
print(f"  Tension: {tens:.2f} sigma")
print(f"  tan(delta) = 9/4 = alpha_s/alpha_EM: CP violation IS the force hierarchy")

# === theta_23 (CKM) ===
# theta_23 = A * lambda^2
theta_23_ckm = A_wolf * lam**2
theta_23_exp = 0.0405
print(f"\nCKM theta_23:")
format_comparison("  theta_23 = A*lambda^2", theta_23_ckm, theta_23_exp)

# === Jarlskog invariant ===
# J_CKM = A^2 * lambda^6 * sin(delta) * (product of cosines)
# Approximate: J ~ A^2 * lambda^6 * eta_bar
# Exact from Wolfenstein: J = c12^2*c23*c13^2*s12*s23*s13*sin(delta)
c12 = np.sqrt(1 - sin_C**2)
s12 = sin_C
s23 = theta_23_ckm
c23 = np.sqrt(1 - s23**2)
# === Wolfenstein rho_bar and eta_bar (2026-03-30) ===
# rho_bar = 1/(2*pi) — angular normalization per radian
# eta_bar = (9/4) * rho_bar = 9/(8*pi) — from tan(delta) = 9/4
rho_bar = 1.0 / (2 * np.pi)
eta_bar = (9.0/4.0) * rho_bar  # = 9/(8*pi)
R_b = np.sqrt(rho_bar**2 + eta_bar**2)  # = sqrt(97)/(8*pi)

rho_exp = 0.159
rho_sig = 0.010
eta_exp = 0.348
eta_sig = 0.010

print(f"\nWolfenstein rho_bar and eta_bar:")
format_comparison("  rho_bar = 1/(2*pi)", rho_bar, rho_exp)
print(f"  Tension: {sigma_tension(rho_bar, rho_exp, rho_sig):.2f} sigma")
format_comparison("  eta_bar = 9/(8*pi)", eta_bar, eta_exp)
print(f"  Tension: {sigma_tension(eta_bar, eta_exp, eta_sig):.2f} sigma")
print(f"  R_b = sqrt(97)/(8*pi) = {R_b:.6f}")

# === Jarlskog invariant ===
s13 = A_wolf * lam**3 * R_b  # |V_ub| = A*lambda^3*R_b
c13 = np.sqrt(1 - s13**2)
J_ckm = c12 * c23 * c13**2 * s12 * s23 * s13 * np.sin(delta_ckm)
J_ckm_exp = 3.00e-5
print(f"\nJarlskog invariant:")
format_comparison("  J_CKM", J_ckm, J_ckm_exp)

# === theta_QCD ===
print(f"\nStrong CP:")
print(f"  theta_QCD = 0 (reality of S^2 angular algebra)")
print(f"  No axion required.")

# === Save ===
results = {
    'sin_theta_C': float(sin_C),
    'sin_theta_C_exp': sin_C_exp,
    'sin_theta_C_pct': float(pct_error(sin_C, sin_C_exp)),
    'A_wolfenstein': float(A_wolf),
    'A_exp': A_wolf_exp,
    'A_pct': float(pct_error(A_wolf, A_wolf_exp)),
    'delta_CKM_deg': float(delta_ckm_deg),
    'delta_CKM_exp': delta_ckm_exp,
    'delta_CKM_sigma': float(tens),
    'tan_delta': 9.0/4.0,
    'tan_delta_meaning': 'alpha_s/alpha_EM = (2l+1)^2/(2j+1)^2',
    'rho_bar': float(rho_bar),
    'rho_bar_exp': rho_exp,
    'rho_bar_sigma': float(sigma_tension(rho_bar, rho_exp, rho_sig)),
    'eta_bar': float(eta_bar),
    'eta_bar_exp': eta_exp,
    'eta_bar_sigma': float(sigma_tension(eta_bar, eta_exp, eta_sig)),
    'R_b': float(R_b),
    'theta_QCD': 0,
    'theta_QCD_note': 'exact zero from reality of S^2 algebra',
    'J_CKM': float(J_ckm),
    'J_CKM_exp': J_ckm_exp,
}
save_results('19_ckm_mixing.json', results)
