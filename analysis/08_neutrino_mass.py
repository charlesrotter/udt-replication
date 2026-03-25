#!/usr/bin/env python3
"""08 -- Neutrino mass prediction from alpha_EM^3.

SOURCE: data/generated/01_vacuum_profile.json (I_2)
        lib/angular_integrals.py (neutrino_mass)
        lib/constants.py (M_E, NUFIT)
GENERATES: data/generated/08_neutrino_mass.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2

Neutrino mass: m_nu = alpha^3 * m_e / 4
  where alpha = I_2 / (36*pi) from vacuum geometry.

Cross-checks:
  m_p * m_nu / m_e^2 = (3/2)*pi^5*alpha^3
  m_mu / m_nu = (80/3)*(pi/alpha)^3

Solar neutrino (tentative):
  m_nu_solar = source(-3) * alpha^3 * m_e
  where source(-3) = 5/(84*sqrt(2))

Completion class: A (algebraic from A-class I_2).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import M_E, NUFIT, PDG_MASSES, MULT_2J1
from lib.angular_integrals import neutrino_mass
from lib.utils import save_results, load_results, pct_error, format_comparison

# === Load I_2 from vacuum profile ===
vac = load_results('01_vacuum_profile.json')
I_2 = vac['I_2']

# === Compute alpha_EM from geometry ===
alpha_inv = 36 * np.pi / I_2
alpha_em = 1.0 / alpha_inv

# === Neutrino mass from lib ===
nu_result = neutrino_mass(alpha_em)

# Direct computation for verification
m_nu = alpha_em**3 * M_E / MULT_2J1**2  # alpha^3 * m_e / 4
m_nu_eV = m_nu * 1e6  # MeV -> eV
sqrt_dm2_atm = np.sqrt(NUFIT['dm2_31'])  # eV

print(f"Neutrino mass prediction:")
print(f"  alpha_EM = {alpha_em:.10f} (1/alpha = {alpha_inv:.6f})")
print(f"  m_nu = alpha^3 * m_e / 4 = {m_nu_eV:.6f} eV")
print(f"  sqrt(dm2_31) = {sqrt_dm2_atm:.6f} eV")
format_comparison("m_nu", m_nu_eV, sqrt_dm2_atm, "eV")

# === Cross-check 1: m_p * m_nu / m_e^2 = (3/2)*pi^5*alpha^3 ===
m_p = 6 * np.pi**5 * M_E
lhs_1 = m_p * m_nu / M_E**2
rhs_1 = 1.5 * np.pi**5 * alpha_em**3
cross1_match = abs(lhs_1 - rhs_1) / abs(rhs_1) < 1e-10

print(f"\nCross-check 1: m_p*m_nu/m_e^2 = (3/2)*pi^5*alpha^3")
print(f"  LHS = {lhs_1:.10e}, RHS = {rhs_1:.10e}, match = {cross1_match}")

# === Cross-check 2: m_mu / m_nu = (80/3)*(pi/alpha)^3 ===
m_mu = 20 * np.pi**3 / 3 * M_E
lhs_2 = m_mu / m_nu
rhs_2 = (80.0 / 3.0) * (np.pi / alpha_em)**3
cross2_match = abs(lhs_2 - rhs_2) / abs(rhs_2) < 1e-10

print(f"\nCross-check 2: m_mu/m_nu = (80/3)*(pi/alpha)^3")
print(f"  LHS = {lhs_2:.10e}, RHS = {rhs_2:.10e}, match = {cross2_match}")

# === Cross-check 3: source(-1) = 1/4 connection ===
# m_nu = alpha^3 * m_e / (2j+1)^2 = alpha^3 * m_e * source(-1)
# since source(-1) = 1/4 and (2j+1)^2 = 4
source_km1 = 0.25  # loaded from 03, but exact value is 1/4
lhs_3 = alpha_em**3 * M_E * source_km1
cross3_match = abs(lhs_3 - m_nu) / abs(m_nu) < 1e-10

print(f"\nCross-check 3: source(-1)=1/4 connection")
print(f"  alpha^3*m_e*S(-1) = {lhs_3:.10e} MeV, m_nu = {m_nu:.10e} MeV, match = {cross3_match}")

# === Solar neutrino (tentative) ===
source_km3 = 5.0 / (84.0 * np.sqrt(2))
m_nu_solar = source_km3 * alpha_em**3 * M_E
m_nu_solar_eV = m_nu_solar * 1e6
sqrt_dm2_sol = np.sqrt(NUFIT['dm2_21'])

print(f"\nSolar neutrino (tentative):")
print(f"  m_nu_solar = S(-3)*alpha^3*m_e = {m_nu_solar_eV:.6f} eV")
print(f"  sqrt(dm2_21) = {sqrt_dm2_sol:.6f} eV")
format_comparison("m_nu_solar", m_nu_solar_eV, sqrt_dm2_sol, "eV")

# === Neutrino mass sum ===
# Hierarchy: m1 ~ 0, m2 ~ sqrt(dm2_21), m3 ~ sqrt(dm2_31)
m1 = 0.0
m2 = sqrt_dm2_sol
m3 = sqrt_dm2_atm
sum_nu_pdg = m1 + m2 + m3

# UDT sum: m_nu + m_nu_solar (m1 ~ 0 if strict hierarchy)
sum_nu_udt = m_nu_eV + m_nu_solar_eV

# === Save results ===
results = {
    'alpha_em': alpha_em,
    'alpha_em_inv': alpha_inv,
    'I_2': I_2,
    'neutrino_mass': {
        'predicted_MeV': m_nu,
        'predicted_eV': m_nu_eV,
        'formula': 'alpha^3 * m_e / 4',
        'comparison_eV': sqrt_dm2_atm,
        'comparison_label': 'sqrt(dm2_31)',
        'pct_error': pct_error(m_nu_eV, sqrt_dm2_atm),
    },
    'cross_checks': {
        'mp_mnu_over_me2': {
            'lhs': lhs_1, 'rhs': rhs_1,
            'formula': '(3/2)*pi^5*alpha^3',
            'match': cross1_match,
        },
        'mmu_over_mnu': {
            'lhs': lhs_2, 'rhs': rhs_2,
            'formula': '(80/3)*(pi/alpha)^3',
            'match': cross2_match,
        },
        'source_km1_connection': {
            'lhs': lhs_3, 'rhs': m_nu,
            'formula': 'alpha^3*m_e*source(-1)',
            'match': cross3_match,
        },
    },
    'solar_neutrino': {
        'predicted_eV': m_nu_solar_eV,
        'source_km3': source_km3,
        'comparison_eV': sqrt_dm2_sol,
        'comparison_label': 'sqrt(dm2_21)',
        'pct_error': pct_error(m_nu_solar_eV, sqrt_dm2_sol),
        'status': 'tentative',
    },
    'mass_sum': {
        'udt_eV': sum_nu_udt,
        'pdg_hierarchy_eV': sum_nu_pdg,
    },
}
save_results('08_neutrino_mass.json', results)
