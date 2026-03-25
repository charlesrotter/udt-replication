#!/usr/bin/env python3
"""01 -- Vacuum phi(r) profile at locked parameters.

SOURCE: lib/vacuum_phi.py (solve, extract_phi2, compute_I2, compute_box_residual)
GENERATES: data/generated/01_vacuum_profile.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  sqrt(-g) = c r^2 sin(theta)  (phi-independent)
  KG equation: box_g phi - mu^2 phi = -S  (positive S deepens well)
  box_g phi = (1/r^2) d/dr(r^2 e^{-2phi} phi')

Completion class: A (finite domain [r_min, r_star]).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import PHI0, MU, MU2, R_STAR, N_GRID, R_MIN
from lib.vacuum_phi import solve, extract_phi2, compute_I2
from lib.utils import save_results, gate_check

# === Solve vacuum ODE ===
r, phi, J, phip, overflow = solve(PHI0, R_STAR, N_GRID, MU, r_min=R_MIN)
assert overflow is None, f"Overflow at r={overflow}"

# === Derived quantities ===
e2phi = np.exp(np.clip(2 * phi, -400, 400))
phi2 = extract_phi2(r, phi, PHI0)
I_2 = compute_I2(r, phi)
phi_min = float(np.min(phi))
phi_min_r = float(r[np.argmin(phi)])

# === Consistency gate: vacuum ODE residual ===
# The ODE is J'(r) = r^2 * mu^2 * phi (vacuum, S=0).
# J is the flux variable solved by RK4; differentiate J (smooth) and compare.
# Also verify phi' = J*e^{2phi}/r^2 consistency.
dr = r[1] - r[0]

# Gate 1: J' vs r^2 mu^2 phi (flux equation)
dJ_dr = np.gradient(J, dr)
rhs_J = r**2 * MU2 * phi
# Skip first/last few points (origin regularization, boundary effects)
interior = slice(50, -50)
res_J = np.abs(dJ_dr[interior] - rhs_J[interior])
max_res_J = float(np.max(res_J))

# Gate 2: phi' consistency: phip = J*e^{2phi}/r^2
phip_from_J = J * e2phi / r**2
phip_from_J[0] = 0.0  # regularize origin
res_phip = np.abs(phip[interior] - phip_from_J[interior])
max_res_phip = float(np.max(res_phip))

# The J-derivative gate uses finite differences on a smooth variable,
# so threshold scales with dr^2 (RK4 is 4th order, np.gradient is 2nd).
# At N=100000, dr ~ 7e-5, dr^2 ~ 5e-9. Use 1e-6 as conservative threshold.
gate_J = gate_check("G1a", max_res_J, 1e-6,
                     "J'(r) = r^2 mu^2 phi (flux ODE, interior)")
gate_phip = gate_check("G1b", max_res_phip, 1e-12,
                       "phi' = J e^{2phi}/r^2 (definition consistency)")
gate_passed = gate_J and gate_phip
max_res = max_res_J  # primary residual for reporting

# === Downsample for JSON output ===
stride = max(1, len(r) // 1000)
idx = np.arange(0, len(r), stride)

# === Save results ===
results = {
    'parameters': {
        'phi0': PHI0,
        'mu': MU,
        'mu2': MU2,
        'r_star': R_STAR,
        'n_grid': N_GRID,
        'r_min': R_MIN,
    },
    'profile': {
        'r': r[idx].tolist(),
        'phi': phi[idx].tolist(),
        'phip': phip[idx].tolist(),
        'e2phi': e2phi[idx].tolist(),
    },
    'I_2': I_2,
    'phi2': phi2,
    'phi_min': phi_min,
    'phi_min_r': phi_min_r,
    'phi_at_rstar': float(phi[-1]),
    'phip_at_rstar': float(phip[-1]),
    'gate_G1a': {
        'name': "J' = r^2 mu^2 phi (flux ODE)",
        'max_residual': max_res_J,
        'threshold': 1e-6,
        'passed': gate_J,
    },
    'gate_G1b': {
        'name': "phi' = J e^{2phi}/r^2 (consistency)",
        'max_residual': max_res_phip,
        'threshold': 1e-12,
        'passed': gate_phip,
    },
}
save_results('01_vacuum_profile.json', results)

# === Summary ===
print(f"Vacuum profile: phi0={PHI0:.5f}, mu={MU:.5f}, r*={R_STAR}")
print(f"  phi_min={phi_min:.6f} at r={phi_min_r:.4f}")
print(f"  I_2 = {I_2:.10f}")
print(f"  phi_2 = {phi2:.10f}")
print(f"  Gate G1: {'PASS' if gate_passed else 'FAIL'} (max |res|={max_res:.2e})")
