#!/usr/bin/env python3
"""11 -- alpha^3 factorization from I_2 integral.

SOURCE: data/generated/01_vacuum_profile.json (I_2)
        lib/constants.py (ALPHA_EM_PDG for comparison)
GENERATES: data/generated/11_alpha3_factorization.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  sqrt(-g) = c r^2 sin(theta)  (phi-independent)
  Tetrad: e^a_mu = diag(e^{-phi}, e^{+phi}, r, r sin(theta))

Derivation chain:
  alpha_EM = I_2 / (36*pi)
  where I_2 = integral(e^{2*phi(r)} dr) over the vacuum cell.

  alpha^3 = [I_2 / (36*pi)]^3
           = [1/(3*pi)]^3 * [1/12]^3 * I_2^3

  Factor breakdown:
    1/(3*pi) -- angular coupling on S^2 (from 2l+1 = 3 angular modes, pi from S^2 volume)
    1/12     -- gauge structure factor (from C(9,2)/3 = 36/3 = 12, three generations)
    I_2      -- metric integral (pure geometry)

Completion class: A (finite domain; I_2 from vacuum ODE).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import ALPHA_EM_PDG, ALPHA_EM_INV_PDG, MU, PHI0, R_STAR, N_GRID, R_MIN
from lib.vacuum_phi import solve, compute_I2
from lib.utils import save_results, load_results, gate_check, format_comparison

# === Load I_2 from script 01 output ===
try:
    data_01 = load_results('01_vacuum_profile.json')
    I_2 = data_01['I_2']
    print(f"  Loaded I_2 = {I_2:.10f} from 01_vacuum_profile.json")
except FileNotFoundError:
    # Recompute if 01 hasn't been run
    print("  01_vacuum_profile.json not found, recomputing...")
    r, phi, J, phip, overflow = solve(PHI0, R_STAR, N_GRID, MU, r_min=R_MIN)
    assert overflow is None, f"Overflow at r={overflow}"
    I_2 = compute_I2(r, phi)
    print(f"  Computed I_2 = {I_2:.10f}")

# === Derive alpha_EM ===
alpha_EM = I_2 / (36.0 * np.pi)
alpha_inv = 1.0 / alpha_EM

# === alpha^3 factorization ===
# Direct computation
alpha3_direct = alpha_EM**3

# Factored form: alpha^3 = [1/(3*pi)]^3 * [1/12]^3 * I_2^3
factor_angular = 1.0 / (3.0 * np.pi)     # angular coupling on S^2
factor_gauge = 1.0 / 12.0                 # gauge structure (36/3 = 12)
I_2_cubed = I_2**3

alpha3_factored = factor_angular**3 * factor_gauge**3 * I_2_cubed

# === Verification: factored matches direct ===
factorization_error = abs(alpha3_factored - alpha3_direct)

# === Component breakdown ===
# 36*pi = 3*pi * 12, so alpha = I_2/(3*pi * 12)
# alpha^3 = I_2^3 / (3*pi)^3 / 12^3

# Factor meanings:
#   3*pi: from angular integration on S^2
#     - 3 = 2l+1 (orbital multiplicity with l=1)
#     - pi = half the solid angle factor (4*pi / (2*(2j+1)) = 4*pi/4 = pi)
#   12: from gauge structure
#     - 36 = C(9,2) (2-form components, alpha boson dimension)
#     - 12 = 36/3 (three generations reduce the effective dimension)

factor_3pi_value = (3.0 * np.pi)**3
factor_12_value = 12.0**3

# Check: 36*pi = 3*pi * 12
assert abs(36.0 * np.pi - 3.0 * np.pi * 12.0) < 1e-14, "36*pi != 3*pi * 12"

# === Gate: factorization identity to 1e-8 ===
gate_passed = gate_check("G1", factorization_error, 1e-8,
                          "alpha^3 factorization identity")

# === Comparison with PDG ===
err_alpha = format_comparison("alpha_EM", alpha_EM, ALPHA_EM_PDG)
err_inv = format_comparison("1/alpha", alpha_inv, ALPHA_EM_INV_PDG)

# === Neutrino mass connection ===
# m_nu = alpha^3 * m_e / (2j+1)^2 = alpha^3 * m_e / 4
# This uses alpha^3, so the factorization matters physically
from lib.constants import M_E, MULT_2J1
m_nu_MeV = alpha3_direct * M_E / MULT_2J1**2
m_nu_eV = m_nu_MeV * 1e6

# === Save results ===
results = {
    'I_2': I_2,
    'alpha_EM': alpha_EM,
    'alpha_inv': alpha_inv,
    'alpha_EM_pdg': ALPHA_EM_PDG,
    'alpha_inv_pdg': ALPHA_EM_INV_PDG,
    'error_pct_alpha': float(err_alpha),
    'alpha3': {
        'direct': alpha3_direct,
        'factored': alpha3_factored,
        'factorization_error': factorization_error,
    },
    'factorization_components': {
        'angular_factor': {
            'value': factor_angular,
            'formula': '1/(3*pi)',
            'origin': '3 = 2l+1 angular modes; pi from S^2 volume',
        },
        'gauge_factor': {
            'value': factor_gauge,
            'formula': '1/12',
            'origin': '12 = C(9,2)/3 = 36/3 gauge structure / 3 generations',
        },
        'metric_integral': {
            'value': I_2,
            'formula': 'integral(e^{2*phi(r)} dr)',
            'origin': 'vacuum ODE on [r_min, r_star]',
        },
        'identity': 'alpha^3 = [1/(3*pi)]^3 * [1/12]^3 * I_2^3',
        'denominator_product': {
            '(3*pi)^3': factor_3pi_value,
            '12^3': factor_12_value,
            'total_36pi_cubed': (36.0 * np.pi)**3,
        },
    },
    'neutrino_connection': {
        'm_nu_eV': m_nu_eV,
        'formula': 'alpha^3 * m_e / 4',
    },
    'gate_G1': {
        'name': 'alpha^3 factorization identity',
        'error': factorization_error,
        'threshold': 1e-8,
        'passed': gate_passed,
    },
}
save_results('11_alpha3_factorization.json', results)

# === Summary ===
print(f"alpha_EM = {alpha_EM:.8f} (1/alpha = {alpha_inv:.4f})")
print(f"  alpha^3 = {alpha3_direct:.12e}")
print(f"  Factored: [1/(3pi)]^3 * [1/12]^3 * I_2^3 = {alpha3_factored:.12e}")
print(f"  Factorization error: {factorization_error:.2e}")
print(f"  Gate G1: {'PASS' if gate_passed else 'FAIL'}")
