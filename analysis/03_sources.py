#!/usr/bin/env python3
"""03 -- Source integrals for each kappa channel.

SOURCE: lib/dirac_formT.py (wavefunction, source_integral), lib/vacuum_phi.py
GENERATES: data/generated/03_sources.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  source(kappa) = integral((G^2 - F^2) e^{2phi} r^2 dr)
  Normalization: integral(e^phi (G^2 + F^2) r^2 dr) = 1  (SL weight)

Completion class: A (finite domain [r_min, r_star]).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import PHI0, MU, R_STAR, N_GRID, R_MIN
from lib.vacuum_phi import solve, extract_phi2
from lib.dirac_formT import wavefunction, source_integral
from lib.utils import save_results, load_results, pct_error

# === Load eigenvalues from 02 ===
ev_data = load_results('02_eigenvalues.json')
eigenvalues = ev_data['eigenvalues']
phi2 = ev_data['parameters']['phi2']

# === Recompute vacuum (full resolution) ===
r, phi, J, phip, overflow = solve(PHI0, R_STAR, N_GRID, MU, r_min=R_MIN)
e2phi = np.exp(np.clip(2 * phi, -400, 400))
ephi = np.exp(np.clip(phi, -200, 200))

# === Algebraic targets for ground-state sources ===
algebraic_targets = {
    -1: {'value': 1.0 / 4.0, 'label': '1/4'},
    +1: {'value': None, 'label': '(no closed form)'},
    -2: {'value': None, 'label': '(no closed form)'},
    +2: {'value': 13.0 / 72.0, 'label': '13/72'},
    -3: {'value': 5.0 / (84.0 * np.sqrt(2)), 'label': '5/(84*sqrt(2))'},
    +3: {'value': None, 'label': '(no closed form)'},
}

# === Compute sources for each kappa, each mode ===
kappas = [-1, +1, -2, +2, -3, +3]
all_sources = {}
algebraic_table = []

for kappa in kappas:
    ks = str(kappa)
    evals = eigenvalues.get(ks, [])
    kappa_sources = []

    for n, E in enumerate(evals):
        G, F = wavefunction(E, kappa, r, phip, e2phi, ephi,
                            phi0=PHI0, phi2=phi2)

        # Verify normalization
        norm = np.trapezoid((G**2 + F**2) * ephi * r**2, r)

        # Compute source integral
        src = source_integral(G, F, r, e2phi)

        kappa_sources.append({
            'n': n + 1,
            'E': E,
            'source': src,
            'norm_check': norm,
        })

    all_sources[ks] = kappa_sources

    # Algebraic comparison (ground state only)
    # Note: for positive kappa, the leading Frobenius component is F,
    # so source = integral(G^2 - F^2) is negative. Compare magnitudes.
    tgt = algebraic_targets[kappa]
    if tgt['value'] is not None and len(kappa_sources) > 0:
        s_computed = kappa_sources[0]['source']
        s_mag = abs(s_computed)
        err = pct_error(s_mag, tgt['value'])
        algebraic_table.append({
            'kappa': kappa,
            'source_computed': s_computed,
            'source_magnitude': s_mag,
            'source_algebraic': tgt['value'],
            'label': tgt['label'],
            'pct_error': err,
        })

# === Save results ===
results = {
    'sources': all_sources,
    'algebraic_comparison': algebraic_table,
    'parameters': {
        'phi0': PHI0, 'mu': float(MU), 'r_star': R_STAR,
        'n_grid': N_GRID, 'phi2': phi2,
    },
}
save_results('03_sources.json', results)

# === Summary ===
print("Source integrals (ground state):")
for kappa in kappas:
    ks = str(kappa)
    srcs = all_sources[ks]
    if srcs:
        s = srcs[0]['source']
        tgt = algebraic_targets[kappa]
        if tgt['value'] is not None:
            err = pct_error(abs(s), tgt['value'])
            print(f"  kappa={kappa:+d}: |S|={abs(s):.10f} vs {tgt['label']}={tgt['value']:.10f} ({err:+.4f}%)")
        else:
            print(f"  kappa={kappa:+d}: S={s:.10f}")
print(f"  {len(algebraic_table)} algebraic comparisons performed")
