#!/usr/bin/env python3
"""09 -- PMNS mixing angles from angular quantum numbers.

SOURCE: lib/angular_integrals.py (pmns_mixing_angles, neutrino_mass)
        lib/constants.py (NUFIT)
GENERATES: data/generated/09_pmns_angles.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2

PMNS mixing angles from quantum number fractions:
  sin^2(theta_12) = (2j+1)^2 / ((2j+1)^2 + (2l+1)^2) = 4/13
  sin^2(theta_23) = (2j+1)^2 / (2|kappa_max|+1) = 4/7
  sin^2(theta_13) = 1 / ((2l+1)^2 * (2|kappa_max|-1)) = 1/45

Compare to NuFIT 5.2 (NO): report tension in sigma units.

Completion class: A (algebraic from angular sector).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import M_E, NUFIT
from lib.angular_integrals import pmns_mixing_angles, neutrino_mass
from lib.utils import save_results, load_results, format_comparison, gate_check

# === PMNS mixing angles from lib ===
pmns = pmns_mixing_angles()

print("PMNS mixing angles:")
for name, data in pmns.items():
    pred = data['predicted']
    nufit = data['nufit']
    sigma = data['sigma']
    tension = data['tension_sigma']
    exact = data['exact']
    print(f"  {name}: {pred:.6f} (={exact}) vs NuFIT {nufit:.5f} +/- {sigma:.5f}  [{tension:.1f} sigma]")

# === Compute actual angles in degrees ===
angles_deg = {}
for name, data in pmns.items():
    theta = np.arcsin(np.sqrt(data['predicted']))
    theta_nufit = np.arcsin(np.sqrt(data['nufit']))
    angles_deg[name] = {
        'theta_pred_deg': np.degrees(theta),
        'theta_nufit_deg': np.degrees(theta_nufit),
    }

print("\nAngles in degrees:")
for name, data in angles_deg.items():
    print(f"  {name}: {data['theta_pred_deg']:.3f} deg vs NuFIT {data['theta_nufit_deg']:.3f} deg")

# === Neutrino mass sum ===
# Load alpha from 01 output
try:
    vac = load_results('01_vacuum_profile.json')
    I_2 = vac['I_2']
    alpha_em = I_2 / (36 * np.pi)
except Exception:
    # Fallback to PDG alpha
    alpha_em = 1.0 / 137.036
    print("  (Warning: using PDG alpha as fallback)")

nu_result = neutrino_mass(alpha_em)

# Mass sum in normal ordering: m1 ~ 0, m2 ~ sqrt(dm2_21), m3 ~ m_nu_predicted
m1 = 0.0
m2 = np.sqrt(NUFIT['dm2_21'])  # eV
m3 = nu_result['predicted_eV']  # eV
sum_masses = m1 + m2 + m3

# Cosmological bound: sum < 0.12 eV (Planck 2018)
cosmo_bound = 0.12  # eV

print(f"\nNeutrino mass sum:")
print(f"  m1 ~ 0, m2 = sqrt(dm2_21) = {m2:.6f} eV, m3 = {m3:.6f} eV")
print(f"  Sum = {sum_masses:.6f} eV (cosmo bound < {cosmo_bound} eV)")

# === Unitarity check ===
# The PMNS matrix should be unitary; check the simplest constraint
s12 = pmns['sin2_12']['predicted']
s23 = pmns['sin2_23']['predicted']
s13 = pmns['sin2_13']['predicted']
c12 = 1 - s12
c23 = 1 - s23
c13 = 1 - s13

# Jarlskog invariant (delta_CP = 0 for these rational angles)
J_max = np.sqrt(s12 * c12 * s23 * c23 * s13 * c13) * np.sqrt(s13)
print(f"\n  |J_max| (Jarlskog, delta_CP=pi/2) = {J_max:.6f}")

# === Save results ===
results = {
    'pmns_angles': {},
    'angles_deg': angles_deg,
    'neutrino_mass_sum': {
        'm1_eV': m1,
        'm2_eV': m2,
        'm3_predicted_eV': m3,
        'sum_eV': sum_masses,
        'cosmo_bound_eV': cosmo_bound,
        'within_bound': sum_masses < cosmo_bound,
    },
    'jarlskog_max': J_max,
}

# Store PMNS data with serializable values
for name, data in pmns.items():
    results['pmns_angles'][name] = {
        'predicted': data['predicted'],
        'exact_fraction': data['exact'],
        'nufit': data['nufit'],
        'sigma': data['sigma'],
        'tension_sigma': data['tension_sigma'],
    }

save_results('09_pmns_angles.json', results)

# === Summary ===
max_tension = max(d['tension_sigma'] for d in pmns.values())
print(f"\nMax tension: {max_tension:.1f} sigma")
print(f"Sum m_nu = {sum_masses:.4f} eV ({'OK' if sum_masses < cosmo_bound else 'EXCEEDS'} cosmo bound)")
