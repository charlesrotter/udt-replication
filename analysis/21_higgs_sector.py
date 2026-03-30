#!/usr/bin/env python3
"""21 -- Higgs sector: VEV, self-coupling, and boson mass.

GENERATES: data/generated/21_higgs_sector.json
VERIFIES: Manuscript Section 18 (The Higgs Mechanism)

v = 504*pi^6*m_e = 247.6 GeV (+0.6%)
  504 = (2j+1)^3 * (2l+1)^2 * (2|kmax|+1) = 8*9*7
lambda = pi/24 where 24 = (2j+1)^3*(2l+1) = 8*3
m_H = v*sqrt(2*lambda) = 126.7 GeV (+1.1%)

Completion class: S (structural, parameter-independent).
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import MULT_2J1, MULT_2L1, MULT_2KP1, M_E
from lib.utils import save_results, pct_error, format_comparison

print("=" * 60)
print("HIGGS SECTOR")
print("=" * 60)

# === Higgs VEV ===
mult_504 = MULT_2J1**3 * MULT_2L1**2 * MULT_2KP1  # 8*9*7 = 504
v_pred = mult_504 * np.pi**6 * M_E / 1000  # GeV
v_exp = 246.22  # GeV
print(f"\nHiggs VEV:")
print(f"  504 = (2j+1)^3*(2l+1)^2*(2|kmax|+1) = {MULT_2J1**3}*{MULT_2L1**2}*{MULT_2KP1} = {mult_504}")
format_comparison("  v = 504*pi^6*m_e", v_pred, v_exp, "GeV")

# === Self-coupling ===
mult_24 = MULT_2J1**3 * MULT_2L1  # 8*3 = 24
lam = np.pi / mult_24
lam_exp = 0.13  # approximate
print(f"\nSelf-coupling:")
print(f"  24 = (2j+1)^3*(2l+1) = {MULT_2J1**3}*{MULT_2L1} = {mult_24}")
format_comparison("  lambda = pi/24", lam, lam_exp)

# === Higgs mass ===
m_H = v_pred * np.sqrt(2 * lam)
m_H_exp = 125.25  # GeV
print(f"\nHiggs boson mass:")
format_comparison("  m_H = v*sqrt(2*lambda)", m_H, m_H_exp, "GeV")

# === v/m_t ratio ===
v_over_mt = MULT_2L1 / MULT_2J1  # 3/2
m_t_pred = v_exp / v_over_mt  # using experimental v
m_t_exp = 172.69
print(f"\nTop quark ratio:")
print(f"  v/m_t = (2l+1)/(2j+1) = 3/2 = {v_over_mt:.4f}")
format_comparison("  m_t = v/(3/2)", v_exp/v_over_mt, m_t_exp, "GeV")

# === Yukawa couplings as capacity fractions ===
print(f"\nYukawa couplings (capacity fractions):")
y_pi = 84 * np.pi / (504 * np.pi**6)
y_p = 6 * np.pi**5 / (504 * np.pi**6)
print(f"  y_pi = 84*pi / (504*pi^6) = {y_pi:.6e}")
print(f"  y_p  = 6*pi^5 / (504*pi^6) = {y_p:.6e}")

results = {
    'v_GeV': float(v_pred),
    'v_exp_GeV': v_exp,
    'v_pct': float(pct_error(v_pred, v_exp)),
    'multiplicity_504': int(mult_504),
    'lambda': float(lam),
    'lambda_exp': lam_exp,
    'multiplicity_24': int(mult_24),
    'm_H_GeV': float(m_H),
    'm_H_exp_GeV': m_H_exp,
    'm_H_pct': float(pct_error(m_H, m_H_exp)),
    'v_over_mt': float(v_over_mt),
}
save_results('21_higgs_sector.json', results)
