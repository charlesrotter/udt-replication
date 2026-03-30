#!/usr/bin/env python3
"""20 -- Quark masses and charges.

GENERATES: data/generated/20_quark_masses.json
VERIFIES: Manuscript Section 17 (Quark Masses and Charges)

N_c = (2l+1) = 3
Q_d = -1/3, Q_u = +2/3
Mass ratios: m_s/m_d = 2*pi^2 (-1.3%), etc.

Completion class: S (structural, parameter-independent).
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import MULT_2J1, MULT_2L1, M_E
from lib.utils import save_results, pct_error, format_comparison

print("=" * 60)
print("QUARK MASSES AND CHARGES")
print("=" * 60)

# === Color ===
N_c = MULT_2L1  # 3
print(f"\nNumber of colors: N_c = (2l+1) = {N_c}")

# === Charges ===
Q_d = -1.0 / N_c
Q_u = 2.0 / N_c
R_ud = N_c * (Q_u**2 + Q_d**2)
print(f"\nQuark charges:")
print(f"  Q_d = -1/N_c = {Q_d:.4f}")
print(f"  Q_u = +2/N_c = {Q_u:.4f}")
print(f"  R(u,d) = N_c*(Q_u^2+Q_d^2) = {R_ud:.4f} (SM: 5/3 = {5/3:.4f})")

# === Base masses ===
m_u_pred = MULT_2J1**2 * M_E  # 4*m_e
m_d_pred = MULT_2L1**2 * M_E  # 9*m_e
# PDG MS-bar masses at 2 GeV
m_u_exp = 2.16  # MeV
m_d_exp = 4.67  # MeV
print(f"\nBase quark masses:")
format_comparison("  m_u = (2j+1)^2 * m_e = 4*m_e", m_u_pred, m_u_exp, "MeV")
format_comparison("  m_d = (2l+1)^2 * m_e = 9*m_e", m_d_pred, m_d_exp, "MeV")
print("  (MS-bar scheme-dependent; approximate)")

# === Inter-generation mass ratios ===
ratios = [
    ("m_s/m_d", 2*np.pi**2, 93.4/4.67, "2*pi^2"),
    ("m_c/m_u", 6*np.pi**4, 1270/2.16, "6*pi^4"),
    ("m_b/m_s", 9*np.pi**2/2, 4180/93.4, "9*pi^2/2"),
    ("m_t/m_c", 14*np.pi**2, 172690/1270, "14*pi^2"),
]

print(f"\nInter-generation mass ratios:")
ratio_results = []
for name, pred, exp, formula in ratios:
    err = pct_error(pred, exp)
    print(f"  {name} = {formula} = {pred:.2f} (exp: {exp:.1f}, err: {err:+.2f}%)")
    ratio_results.append({
        'name': name,
        'formula': formula,
        'predicted': float(pred),
        'experimental': float(exp),
        'pct_error': float(err),
    })

# === Top quark via Higgs ===
# v/m_t = 3/2 = (2l+1)/(2j+1)
v_over_mt = MULT_2L1 / MULT_2J1  # 3/2
m_t_from_v = 246220 / v_over_mt  # MeV
m_t_exp = 172690  # MeV
print(f"\nTop quark from v/m_t = 3/2:")
format_comparison("  m_t = v/(3/2)", m_t_from_v/1000, m_t_exp/1000, "GeV")

results = {
    'N_c': int(N_c),
    'Q_d': float(Q_d),
    'Q_u': float(Q_u),
    'R_ud': float(R_ud),
    'm_u_pred_MeV': float(m_u_pred),
    'm_d_pred_MeV': float(m_d_pred),
    'mass_ratios': ratio_results,
    'm_t_from_v_GeV': float(m_t_from_v/1000),
}
save_results('20_quark_masses.json', results)
