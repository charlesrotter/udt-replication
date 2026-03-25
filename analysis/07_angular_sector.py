#!/usr/bin/env python3
"""07 -- Angular sector mass predictions (parameter-free).

SOURCE: lib/angular_integrals.py (angular_mass_formulas, diophantine_check)
        lib/constants.py (M_E, PDG_MASSES)
GENERATES: data/generated/07_angular_sector.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2

Angular sector masses derive from the spinor harmonic coupling
on S^2, with no radial integration needed:
  m_pi = C(9,3) * pi * m_e = 84 * pi * m_e
  m_mu = (2j+1)^2*(2|kmax|-1)/(2l+1) * pi^{2l+1} * m_e = 20*pi^3/3 * m_e
  m_p  = (2j+1)(2l+1) * pi^{2|kmax|-1} * m_e = 6*pi^5 * m_e

The Diophantine identity (2j+1)^2(2l+1)(2|kmax|+1) = C(2l+2|kmax|+1, 2l+1)
selects j=1/2, l=1, |kmax|=3 uniquely.

Completion class: A (algebraic from angular sector).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import M_E, PDG_MASSES
from lib.angular_integrals import angular_mass_formulas, diophantine_check
from lib.utils import save_results, format_comparison, gate_check

# === Diophantine verification ===
dioph = diophantine_check()
gate_check("D1", 0.0 if dioph['match'] else 1.0, 0.5,
           f"Diophantine: LHS={dioph['LHS']} == RHS={dioph['RHS']}")
print(f"  Quantum numbers: j={dioph['quantum_numbers']['j']}, "
      f"l={dioph['quantum_numbers']['l']}, "
      f"|kmax|={dioph['quantum_numbers']['kappa_max']}")

# === Angular sector masses ===
masses = angular_mass_formulas()

print("\nAngular sector mass predictions (parameter-free):")
for name, data in masses.items():
    if name == 'electron':
        print(f"  {name:8s}: {data['predicted']:.5f} MeV (anchor)")
    else:
        format_comparison(name, data['predicted'], data['pdg'], "MeV")

# === Explicit numerical checks ===
m_pi_check = 84 * np.pi * M_E
m_mu_check = 20 * np.pi**3 / 3 * M_E
m_p_check = 6 * np.pi**5 * M_E

assert abs(m_pi_check - masses['pion']['predicted']) < 1e-10
assert abs(m_mu_check - masses['muon']['predicted']) < 1e-10
assert abs(m_p_check - masses['proton']['predicted']) < 1e-10

# === Additional ratios ===
ratios = {
    'm_p/m_e': {
        'value': m_p_check / M_E,
        'formula': '6*pi^5',
        'exact': 6 * np.pi**5,
    },
    'm_mu/m_e': {
        'value': m_mu_check / M_E,
        'formula': '20*pi^3/3',
        'exact': 20 * np.pi**3 / 3,
    },
    'm_p/m_pi': {
        'value': m_p_check / m_pi_check,
        'formula': 'pi^4/14',
        'exact': np.pi**4 / 14,
    },
    'm_mu/m_pi': {
        'value': m_mu_check / m_pi_check,
        'formula': '(5/63)*pi^2',
        'exact': 5 * np.pi**2 / 63,
    },
}

print("\nMass ratios:")
for name, data in ratios.items():
    print(f"  {name}: {data['value']:.10f} = {data['formula']} = {data['exact']:.10f}")

# === Save results ===
results = {
    'diophantine': dioph,
    'masses': masses,
    'ratios': ratios,
    'formulas': {
        'm_pi': '84*pi*m_e',
        'm_mu': '20*pi^3/3 * m_e',
        'm_p': '6*pi^5 * m_e',
    },
}
save_results('07_angular_sector.json', results)
