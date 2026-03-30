#!/usr/bin/env python3
"""05 -- Fine structure constant: three-level derivation chain.

SOURCE: data/generated/01_vacuum_profile.json, data/generated/03_sources.json
GENERATES: data/generated/05_alpha_em.json
VERIFIES: Manuscript Section 12 (The Fine-Structure Constant)

Three-level chain:
  Level 1: Algebraic bridge 36*pi/I_2 = 137.43 (+0.29%)
  Level 2: Computed ODE alpha_0 = source^2*E_1^2*I_2/(2*pi) -> 1/137.14 (+0.074%)
  Level 3: One-loop alpha = alpha_0*(1 + alpha_0/pi^2) -> 1/137.036 (+0.0004%)

CRITICAL: alpha_0 uses computed ODE values, NOT algebraic bridge 36*pi/I_2.

Completion class: A (finite domain).
"""
import os, sys
import numpy as np
from math import comb

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import (ALPHA_EM_PDG, ALPHA_EM_INV_PDG, C_CALIB, M_E,
                           PDG_MASSES, MULT_2J1, MULT_2L1, MULT_2KM1)
from lib.utils import save_results, load_results, pct_error, format_comparison

# === Load computed values ===
vac = load_results('01_vacuum_profile.json')
I_2 = vac['I_2']

src = load_results('03_sources.json')
source_km1 = src['sources']['km1'][0]['source']
E1_km1 = src['sources']['km1'][0]['E']

print("=" * 60)
print("FINE STRUCTURE CONSTANT: THREE-LEVEL CHAIN")
print("=" * 60)

# === Level 1: Algebraic bridge ===
alpha_inv_alg = 36 * np.pi / I_2
alpha_alg = 1.0 / alpha_inv_alg
print("\nLevel 1 (algebraic bridge):")
format_comparison("  1/alpha_0^alg = 36*pi/I_2", alpha_inv_alg, ALPHA_EM_INV_PDG)

# === Level 2: Computed ODE values ===
source2_E12 = source_km1**2 * E1_km1**2
alpha_0_ode = source2_E12 * I_2 / (2 * np.pi)
alpha_inv_ode = 1.0 / alpha_0_ode
print("\nLevel 2 (computed ODE):")
print(f"  source(km1) = {source_km1:.8f}")
print(f"  E_1(km1) = {E1_km1:.8f}")
print(f"  source^2*E_1^2 = {source2_E12:.8f} (algebraic 1/18 = {1/18:.8f}, diff = {pct_error(source2_E12, 1/18):+.4f}%)")
format_comparison("  1/alpha_0^ode", alpha_inv_ode, ALPHA_EM_INV_PDG)

# === Level 3: One-loop angular correction ===
alpha_corr = alpha_0_ode * (1 + alpha_0_ode / np.pi**2)
alpha_inv_corr = 1.0 / alpha_corr
print("\nLevel 3 (one-loop correction):")
print(f"  alpha_0/pi^2 = {alpha_0_ode/np.pi**2:.8f}")
format_comparison("  1/alpha (corrected)", alpha_inv_corr, ALPHA_EM_INV_PDG)

# === Hydrogen binding energy ===
E_H = 0.5 * alpha_corr**2 * M_E * 1e6  # eV
E_H_exp = 13.605693  # eV
print(f"\nHydrogen binding energy:")
print(f"  E_H = alpha^2*m_e/2 = {E_H:.6f} eV (exp: {E_H_exp:.6f}, err: {pct_error(E_H, E_H_exp):+.6f}%)")

# === p-n mass splitting (using corrected alpha) ===
pn_factor = MULT_2KM1 / MULT_2J1**2  # 5/4
delta_mn_mp_pred = pn_factor * alpha_corr * C_CALIB
delta_mn_mp_pdg = PDG_MASSES['n'] - PDG_MASSES['p']
print(f"\np-n splitting:")
format_comparison("  m_n - m_p", delta_mn_mp_pred, delta_mn_mp_pdg, "MeV")

# === Decomposition of 36 ===
decomposition = {
    'factor_36': 36,
    '(2j+1)^2': int(MULT_2J1**2),
    '(2l+1)^2': int(MULT_2L1**2),
    'C(9,2)': comb(9, 2),
    'match': MULT_2J1**2 * MULT_2L1**2 == 36 == comb(9, 2),
}

# === Save results ===
results = {
    'I_2': I_2,
    'three_level_chain': {
        'level_1_algebraic': {
            'alpha_inv': float(alpha_inv_alg),
            'pct_error': float(pct_error(alpha_inv_alg, ALPHA_EM_INV_PDG)),
        },
        'level_2_computed': {
            'source_km1': float(source_km1),
            'E1_km1': float(E1_km1),
            'source2_E12': float(source2_E12),
            'alpha_0': float(alpha_0_ode),
            'alpha_inv': float(alpha_inv_ode),
            'pct_error': float(pct_error(alpha_inv_ode, ALPHA_EM_INV_PDG)),
        },
        'level_3_corrected': {
            'correction_factor': float(alpha_0_ode / np.pi**2),
            'alpha': float(alpha_corr),
            'alpha_inv': float(alpha_inv_corr),
            'pct_error': float(pct_error(alpha_inv_corr, ALPHA_EM_INV_PDG)),
        },
    },
    'hydrogen_binding_eV': float(E_H),
    'hydrogen_pct_error': float(pct_error(E_H, E_H_exp)),
    'decomposition': decomposition,
    'pn_splitting': {
        'factor': float(pn_factor),
        'predicted_MeV': float(delta_mn_mp_pred),
        'pdg_MeV': float(delta_mn_mp_pdg),
        'pct_error': float(pct_error(delta_mn_mp_pred, delta_mn_mp_pdg)),
    },
}
save_results('05_alpha_em.json', results)
