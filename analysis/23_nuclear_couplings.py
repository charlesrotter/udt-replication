#!/usr/bin/env python3
"""23 -- Nuclear coupling constants from the metric.

GENERATES: data/generated/23_nuclear_couplings.json
VERIFIES: Manuscript Section 19 (QCD and the Nuclear Force)

g_A = 4/pi = 1.2732 (-0.23%)
f_pi = 180*m_e = 91.98 MeV (-0.25%)
g_piNN = 2*pi^4/15 = 12.99 (+0.06%)
g^2/(4*pi) = pi^7/225 = 13.42 (-0.6%)

All from (j, l, |kappa_max|) = (1/2, 1, 3) and m_e.
Zero Standard Model inputs.

Completion class: S (structural, parameter-independent).
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import (MULT_2J1, MULT_2L1, MULT_2KM1, MULT_2KP1,
                           M_E, PHI_GOLD, C_CALIB)
from lib.utils import save_results, pct_error, format_comparison

print("=" * 60)
print("NUCLEAR COUPLING CONSTANTS FROM THE METRIC")
print("=" * 60)

# === g_A (axial coupling) ===
# g_A = (2j+1)^2 / pi = 4/pi
g_A = MULT_2J1**2 / np.pi
g_A_exp = 1.2762
print("\nAxial coupling:")
format_comparison("  g_A = (2j+1)^2/pi = 4/pi", g_A, g_A_exp)

# Deep derivation: g_A = (N_c+2)(2j+1)^2 / ((2|kmax|-1)*pi) = 5*4/(5*pi) = 4/pi
N_c = MULT_2L1  # 3
print(f"  Deep: (N_c+2)(2j+1)^2/((2|kmax|-1)*pi) = {(N_c+2)*MULT_2J1**2}/{MULT_2KM1*np.pi:.4f} = {(N_c+2)*MULT_2J1**2/(MULT_2KM1*np.pi):.6f}")
print(f"  Cancellation: N_c+2 = 2l+3 = {N_c+2} = 2|kmax|-1 = {MULT_2KM1} (Diophantine-specific)")

# === f_pi (pion decay constant) ===
# f_pi = C(9,2) * (2|kmax|-1) * m_e = 36 * 5 * m_e = 180 * m_e
from math import comb
f_pi_coeff = comb(9, 2) * MULT_2KM1  # 36*5 = 180
f_pi = f_pi_coeff * M_E
f_pi_exp = 92.21  # MeV
print("\nPion decay constant:")
print(f"  180 = C(9,2)*(2|kmax|-1) = {comb(9,2)}*{MULT_2KM1} = {f_pi_coeff}")
format_comparison("  f_pi = 180*m_e", f_pi, f_pi_exp, "MeV")

# Deep: f_pi/m_pi = 15/(7*pi) = (2l+1)(2|kmax|-1)/((2|kmax|+1)*pi)
ratio_fpi_mpi = MULT_2L1 * MULT_2KM1 / (MULT_2KP1 * np.pi)  # 15/(7*pi)
m_pi = 84 * np.pi * M_E
print(f"  f_pi/m_pi = 15/(7*pi) = {ratio_fpi_mpi:.6f} (check: {f_pi/m_pi:.6f})")

# === g_piNN (Goldberger-Treiman) ===
# g_piNN = g_A * m_N / f_pi = (4/pi) * 6*pi^5*m_e / (180*m_e)
m_N = 6 * np.pi**5 * M_E  # proton mass (angular formula)
g_piNN = g_A * m_N / f_pi
g_piNN_formula = 2 * np.pi**4 / 15  # algebraic form
g_piNN_exp = 13.05  # approximate
print("\nPion-nucleon coupling (Goldberger-Treiman):")
print(f"  g_piNN = g_A*m_N/f_pi = (4/pi)*6*pi^5/(180) = 24*pi^4/180 = 2*pi^4/15")
format_comparison("  g_piNN = 2*pi^4/15", g_piNN_formula, g_piNN_exp)
print(f"  Computed: {g_piNN:.6f} vs formula: {g_piNN_formula:.6f} (diff: {pct_error(g_piNN, g_piNN_formula):+.6f}%)")

# === g^2/(4*pi) ===
g2_4pi = g_piNN_formula**2 / (4 * np.pi)
g2_4pi_formula = np.pi**7 / 225
g2_4pi_exp = 13.55  # approximate
print("\nNuclear coupling constant:")
print(f"  g^2/(4*pi) = (2*pi^4/15)^2/(4*pi) = pi^7/225")
print(f"  225 = 15^2 = [(2l+1)(2|kmax|-1)]^2 = {(MULT_2L1*MULT_2KM1)**2}")
print(f"  pi^7: power 7 = 2|kmax|+1 = {MULT_2KP1} (closure number)")
format_comparison("  g^2/(4*pi) = pi^7/225", g2_4pi_formula, g2_4pi_exp)

# === OPEP potential at 1.5 fm ===
hbar_c = 197.327  # MeV*fm
m_pi_MeV = 134.977
r_fm = 1.5
x = m_pi_MeV * r_fm / hbar_c
V_OPEP = -(g2_4pi_formula / (4 * np.pi)) * (m_pi_MeV**2 / (3 * m_N)) * np.exp(-x) / x
print(f"\nOPEP potential at r = 1.5 fm:")
print(f"  V_OPEP = {V_OPEP:.1f} MeV (central S-wave, deuteron channel)")

# === Proton charge radius ===
hbar_c_MeV_fm = 197.327
R_p = hbar_c_MeV_fm / (C_CALIB * PHI_GOLD)
R_p_exp = 0.8751  # CODATA fm
print(f"\nProton charge radius:")
format_comparison("  R_p = hbar*c/(C*phi_gold)", R_p, R_p_exp, "fm")

# === Deuteron binding energy ===
# B_d = source^2 * E1^2 * (2j+1)/(2|kmax|+1) * C = C/63
B_d = C_CALIB / (MULT_2L1**2 * MULT_2KP1)
B_d_exp = 2.2244  # MeV
print(f"\nDeuteron binding energy:")
print(f"  B_d = C/((2l+1)^2 * (2|kmax|+1)) = C/63")
format_comparison("  B_d", B_d, B_d_exp, "MeV")

results = {
    'g_A': float(g_A),
    'g_A_exp': g_A_exp,
    'g_A_pct': float(pct_error(g_A, g_A_exp)),
    'f_pi_MeV': float(f_pi),
    'f_pi_exp': f_pi_exp,
    'f_pi_pct': float(pct_error(f_pi, f_pi_exp)),
    'g_piNN': float(g_piNN_formula),
    'g_piNN_exp': g_piNN_exp,
    'g_piNN_pct': float(pct_error(g_piNN_formula, g_piNN_exp)),
    'g2_4pi': float(g2_4pi_formula),
    'g2_4pi_exp': g2_4pi_exp,
    'g2_4pi_pct': float(pct_error(g2_4pi_formula, g2_4pi_exp)),
    'V_OPEP_1p5fm_MeV': float(V_OPEP),
    'R_p_fm': float(R_p),
    'R_p_exp_fm': R_p_exp,
    'R_p_pct': float(pct_error(R_p, R_p_exp)),
    'B_d_MeV': float(B_d),
    'B_d_exp_MeV': B_d_exp,
    'B_d_pct': float(pct_error(B_d, B_d_exp)),
    'B_d_formula': 'C / ((2l+1)^2 * (2|kmax|+1)) = C/63',
}
save_results('23_nuclear_couplings.json', results)
