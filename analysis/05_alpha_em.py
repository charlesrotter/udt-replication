#!/usr/bin/env python3
"""05 -- Fine structure constant from vacuum geometry.

SOURCE: data/generated/01_vacuum_profile.json
GENERATES: data/generated/05_alpha_em.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  I_2 = integral(e^{2*phi(r)} dr)
  1/alpha = 36*pi / I_2

Decomposition: 36 = (2j+1)^2 * (2l+1)^2 = 4 * 9 = C(9,2)
  j=1/2, l=1 from Diophantine selection.

p-n mass splitting: m_n - m_p = (5/4) * alpha_EM * C_CALIB
  Factor 5/4 = (2|kappa_max|-1) / (2j+1)^2

Completion class: A (uses A-class I_2 from 01).
"""
import os
import sys
import numpy as np
from math import comb

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import (ALPHA_EM_PDG, ALPHA_EM_INV_PDG, C_CALIB, M_E,
                           PDG_MASSES, MULT_2J1, MULT_2L1, MULT_2KM1)
from lib.utils import save_results, load_results, pct_error, format_comparison

# === Load I_2 from vacuum profile ===
vac = load_results('01_vacuum_profile.json')
I_2 = vac['I_2']

# === Compute alpha_EM ===
alpha_inv = 36 * np.pi / I_2
alpha_em = 1.0 / alpha_inv

# === Decomposition of the factor 36 ===
decomposition = {
    'factor_36': 36,
    '(2j+1)^2': MULT_2J1**2,       # = 4
    '(2l+1)^2': MULT_2L1**2,       # = 9
    'product': MULT_2J1**2 * MULT_2L1**2,  # = 36
    'C(9,2)': comb(9, 2),           # = 36
    'match': MULT_2J1**2 * MULT_2L1**2 == 36 == comb(9, 2),
}

# === p-n mass splitting ===
# m_n - m_p = (5/4) * alpha * C_CALIB
# factor 5/4 = (2|kappa_max|-1)/(2j+1)^2
pn_factor = MULT_2KM1 / MULT_2J1**2  # 5/4
delta_mn_mp_pred = pn_factor * alpha_em * C_CALIB
delta_mn_mp_pdg = PDG_MASSES['n'] - PDG_MASSES['p']

# === Print comparisons ===
print("Fine structure constant from vacuum geometry:")
format_comparison("1/alpha", alpha_inv, ALPHA_EM_INV_PDG)
format_comparison("alpha", alpha_em, ALPHA_EM_PDG)
print(f"\n  36 = (2j+1)^2*(2l+1)^2 = {MULT_2J1**2}*{MULT_2L1**2} = C(9,2) = {comb(9,2)}: {'MATCH' if decomposition['match'] else 'FAIL'}")
print(f"\np-n splitting:")
format_comparison("m_n - m_p", delta_mn_mp_pred, delta_mn_mp_pdg, "MeV")

# === Save results ===
results = {
    'I_2': I_2,
    'alpha_em': alpha_em,
    'alpha_em_inv': alpha_inv,
    'alpha_em_inv_pdg': ALPHA_EM_INV_PDG,
    'alpha_pct_error': pct_error(alpha_inv, ALPHA_EM_INV_PDG),
    'decomposition': decomposition,
    'pn_splitting': {
        'factor': pn_factor,
        'factor_label': '(2|kmax|-1)/(2j+1)^2 = 5/4',
        'predicted_MeV': delta_mn_mp_pred,
        'pdg_MeV': delta_mn_mp_pdg,
        'pct_error': pct_error(delta_mn_mp_pred, delta_mn_mp_pdg),
    },
}
save_results('05_alpha_em.json', results)
