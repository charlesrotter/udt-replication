#!/usr/bin/env python3
"""06 -- Meson and baryon mass predictions from Dirac eigenvalues.

SOURCE: data/generated/02_eigenvalues.json, data/generated/03_sources.json
        lib/constants.py (C_CALIB, PHI_GOLD, R_STAR, PDG_MASSES)
GENERATES: data/generated/06_particle_masses.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  Mass formula: m = E_n(kappa) * C_CALIB
  C_CALIB = 4*pi^2 * m_e * r* = 140.95 MeV

Baryon sub-cavity: r_b = r* / phi_gold
  phi_gold = (1+sqrt(5))/2, the golden ratio
  Baryons are eigenvalues on [0, r_b] with same Neumann BC.

Density control: err/gap metric for spectral comparison.

Completion class: A (finite domain).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import (PHI0, MU, R_STAR, N_GRID, R_MIN, C_CALIB,
                           PHI_GOLD, PDG_MASSES, M_E)
from lib.vacuum_phi import solve, extract_phi2
from lib.dirac_formT import find_eigenvalues
from lib.utils import save_results, load_results, pct_error

# === Load meson eigenvalues from 02 ===
ev_data = load_results('02_eigenvalues.json')
meson_eigenvalues = ev_data['eigenvalues']
phi2 = ev_data['parameters']['phi2']

# === Meson mass assignments ===
# m = E_n(kappa) * C_CALIB
meson_assignments = [
    # (kappa, n_index, particle_name)
    (-1, 0, 'eta'),       # E1(kappa=-1)
    (-1, 1, 'rho'),       # E2(kappa=-1)
    (-1, 2, 'omega'),     # E3(kappa=-1) -- or separate kappa
    (+1, 0, 'K'),         # E1(kappa=+1)
    (+1, 1, 'phi_m'),     # E2(kappa=+1)
    (-2, 0, 'f2'),        # E1(kappa=-2)
    (+2, 0, 'a2'),        # E1(kappa=+2)
]

meson_table = []
for kappa, n_idx, name in meson_assignments:
    ks = str(kappa)
    evals = meson_eigenvalues.get(ks, [])
    if n_idx < len(evals):
        E = evals[n_idx]
        m_pred = E * C_CALIB
        m_pdg = PDG_MASSES.get(name, None)
        entry = {
            'particle': name,
            'kappa': kappa,
            'n': n_idx + 1,
            'E': E,
            'mass_predicted_MeV': m_pred,
        }
        if m_pdg is not None:
            entry['mass_pdg_MeV'] = m_pdg
            entry['pct_error'] = pct_error(m_pred, m_pdg)
        meson_table.append(entry)

# === Baryon sub-cavity ===
r_b = R_STAR / PHI_GOLD
print(f"Baryon sub-cavity: r_b = r*/phi_gold = {R_STAR}/{PHI_GOLD:.5f} = {r_b:.6f}")

# Solve vacuum on [r_min, r_b]
n_baryon = N_GRID
r_bar, phi_bar, J_bar, phip_bar, overflow_bar = solve(
    PHI0, r_b, n_baryon, MU, r_min=R_MIN
)
e2phi_bar = np.exp(np.clip(2 * phi_bar, -400, 400))
phi2_bar = extract_phi2(r_bar, phi_bar, PHI0)

baryon_kappas = [-1, +1, -2, +2, -3, +3]
baryon_eigenvalues = {}
for kappa in baryon_kappas:
    E_max = 50.0 * abs(kappa)
    evals = find_eigenvalues(kappa, r_bar, phip_bar, e2phi_bar,
                             phi0=PHI0, phi2=phi2_bar,
                             E_min=0.01, E_max=E_max, n_scan=50000, n_modes=5)
    baryon_eigenvalues[str(kappa)] = evals

# Baryon mass assignments
baryon_assignments = [
    (-1, 0, 'p'),        # proton (angular formula gives exact; radial as cross-check)
    (-1, 1, 'Lambda'),
    (-1, 2, 'Sigma'),
    (+1, 0, 'n'),        # neutron
    (+1, 1, 'Xi'),
    (-2, 0, 'Delta'),
    (+2, 0, 'Omega_b'),
]

baryon_table = []
for kappa, n_idx, name in baryon_assignments:
    ks = str(kappa)
    evals = baryon_eigenvalues.get(ks, [])
    if n_idx < len(evals):
        E = evals[n_idx]
        m_pred = E * C_CALIB
        m_pdg = PDG_MASSES.get(name, None)
        entry = {
            'particle': name,
            'kappa': kappa,
            'n': n_idx + 1,
            'E': E,
            'mass_predicted_MeV': m_pred,
            'domain': 'baryon_subcavity',
            'r_b': r_b,
        }
        if m_pdg is not None:
            entry['mass_pdg_MeV'] = m_pdg
            entry['pct_error'] = pct_error(m_pred, m_pdg)
        baryon_table.append(entry)

# === Density control: err/gap ===
# For mesons: compute average gap between adjacent eigenvalues
all_meson_E = sorted([e['E'] for e in meson_table])
if len(all_meson_E) > 1:
    gaps = [all_meson_E[i+1] - all_meson_E[i] for i in range(len(all_meson_E)-1)]
    avg_gap = np.mean(gaps)
    uniform_gap = (all_meson_E[-1] - all_meson_E[0]) / (len(all_meson_E) - 1)
    actual_uniform_ratio = avg_gap / uniform_gap if uniform_gap > 0 else float('inf')
else:
    avg_gap = uniform_gap = actual_uniform_ratio = 0.0

# Compute err/gap for each meson
for entry in meson_table:
    if 'pct_error' in entry and avg_gap > 0:
        entry['err_over_gap'] = abs(entry['pct_error']) / 100.0 * entry['E'] / avg_gap

# === Save results ===
results = {
    'C_CALIB_MeV': C_CALIB,
    'mesons': meson_table,
    'baryons': baryon_table,
    'baryon_subcavity': {
        'r_b': r_b,
        'phi_gold': PHI_GOLD,
        'eigenvalues': baryon_eigenvalues,
    },
    'density_control': {
        'n_meson_eigenvalues': len(all_meson_E),
        'avg_gap': avg_gap,
        'uniform_gap': uniform_gap,
        'actual_uniform_ratio': actual_uniform_ratio,
    },
}
save_results('06_particle_masses.json', results)

# === Summary ===
print(f"\nMeson masses (C = {C_CALIB:.2f} MeV):")
for m in meson_table:
    pdg_str = f" vs {m.get('mass_pdg_MeV', '?'):.1f} ({m.get('pct_error', 0):+.2f}%)" if 'mass_pdg_MeV' in m else ""
    print(f"  {m['particle']:8s}: {m['mass_predicted_MeV']:.1f} MeV{pdg_str}")
print(f"\nBaryon masses (sub-cavity r_b = {r_b:.4f}):")
for b in baryon_table:
    pdg_str = f" vs {b.get('mass_pdg_MeV', '?'):.1f} ({b.get('pct_error', 0):+.2f}%)" if 'mass_pdg_MeV' in b else ""
    print(f"  {b['particle']:8s}: {b['mass_predicted_MeV']:.1f} MeV{pdg_str}")
