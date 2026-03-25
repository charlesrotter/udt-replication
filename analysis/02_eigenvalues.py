#!/usr/bin/env python3
"""02 -- Dirac eigenvalues for kappa in {-1,+1,-2,+2,-3,+3}.

SOURCE: lib/dirac_formT.py (find_eigenvalues, shoot), lib/vacuum_phi.py
GENERATES: data/generated/02_eigenvalues.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  Form-T: G' = (-kappa/r + phi')G + E e^{2phi} F
           F' = (+kappa/r + phi')F - E e^{2phi} G
  BC: G'(r*) = 0 (geometric Neumann)

Completion class: A (finite domain [r_min, r_star]).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import PHI0, MU, R_STAR, N_GRID, R_MIN, E1_EXACT
from lib.vacuum_phi import solve, extract_phi2
from lib.dirac_formT import find_eigenvalues, shoot
from lib.utils import save_results, load_results, gate_check, pct_error

# === Load or recompute vacuum ===
try:
    vac = load_results('01_vacuum_profile.json')
    r = np.linspace(R_MIN, R_STAR, N_GRID)
    # Recompute full-resolution arrays (JSON has downsampled)
    r, phi, J, phip, overflow = solve(PHI0, R_STAR, N_GRID, MU, r_min=R_MIN)
    phi2 = vac['phi2']
except Exception:
    r, phi, J, phip, overflow = solve(PHI0, R_STAR, N_GRID, MU, r_min=R_MIN)
    phi2 = extract_phi2(r, phi, PHI0)

e2phi = np.exp(np.clip(2 * phi, -400, 400))
ephi = np.exp(np.clip(phi, -200, 200))

# === Algebraic targets ===
algebraic_targets = {
    -1: {'E1': 2 * np.sqrt(2) / 3, 'label': '2*sqrt(2)/3'},
    +1: {'E1': None, 'label': '(no known closed form)'},
    -2: {'E1': 16.0 / 3.0, 'label': '16/3'},
    +2: {'E1': 35.0 / 9.0, 'label': '35/9'},
    -3: {'E1': 42.0 / 5.0, 'label': '42/5'},
    +3: {'E1': None, 'label': '(no known closed form)'},
}

# === Find eigenvalues for each kappa ===
kappas = [-1, +1, -2, +2, -3, +3]
all_evals = {}
all_gates = {}
algebraic_table = []

for kappa in kappas:
    # E_max depends on |kappa| (higher kappa -> higher eigenvalues)
    E_max = 40.0 * abs(kappa)
    evals = find_eigenvalues(kappa, r, phip, e2phi, phi0=PHI0, phi2=phi2,
                             E_min=0.01, E_max=E_max, n_scan=50000, n_modes=5)
    all_evals[str(kappa)] = evals

    # Gate: verify G'(r*) = 0 for each eigenvalue
    kappa_gates = []
    for n, E in enumerate(evals):
        bc, G, F = shoot(E, kappa, r, phip, e2phi, phi0=PHI0, phi2=phi2)
        gate_val = abs(bc)
        passed = gate_check(
            f"G1(kappa={kappa:+d},n={n+1})", gate_val, 1e-10,
            f"G'(r*)=0 for E={E:.10f}"
        )
        kappa_gates.append({
            'n': n + 1, 'E': E, 'bc_residual': gate_val, 'passed': passed
        })

    all_gates[str(kappa)] = kappa_gates

    # Algebraic comparison
    tgt = algebraic_targets[kappa]
    if tgt['E1'] is not None and len(evals) > 0:
        err = abs(evals[0] - tgt['E1'])
        err_pct = pct_error(evals[0], tgt['E1'])
        algebraic_table.append({
            'kappa': kappa,
            'E1_computed': evals[0],
            'E1_algebraic': tgt['E1'],
            'label': tgt['label'],
            'abs_error': err,
            'pct_error': err_pct,
        })

# === Specific verification: E1(kappa=-1) vs 2*sqrt(2)/3 ===
E1_km1 = all_evals['-1'][0]
e1_err = abs(E1_km1 - E1_EXACT)
print(f"\n  E1(kappa=-1) = {E1_km1:.15f}")
print(f"  2*sqrt(2)/3  = {E1_EXACT:.15f}")
print(f"  |difference| = {e1_err:.2e}")

# === Save results ===
results = {
    'eigenvalues': all_evals,
    'gates': all_gates,
    'algebraic_comparison': algebraic_table,
    'E1_km1_check': {
        'computed': E1_km1,
        'exact': E1_EXACT,
        'abs_error': e1_err,
    },
    'parameters': {
        'phi0': PHI0, 'mu': float(MU), 'r_star': R_STAR,
        'n_grid': N_GRID, 'phi2': phi2,
    },
}
save_results('02_eigenvalues.json', results)

# === Summary ===
print(f"\nEigenvalues found per kappa:")
for kappa in kappas:
    evals = all_evals[str(kappa)]
    e_str = ', '.join(f'{e:.6f}' for e in evals[:3])
    print(f"  kappa={kappa:+d}: [{e_str}] ({len(evals)} modes)")
