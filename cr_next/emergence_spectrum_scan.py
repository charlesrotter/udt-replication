#!/usr/bin/env python3
"""Emergence spectrum: eigenvalue structure as phi_0 deepens from 0 to -0.809.

Maps the full re-emergence sequence: which particles appear first as the
scalar field reforms after dissolution at phi=0.

Scans phi_0 from -0.01 to -0.81 in 80 steps. For each phi_0, finds all
eigenvalues for kappa = -1, +1, -2, +2, -3, +3 up to E=30.

Output: cr_next/outputs/emergence_spectrum.json
"""
import os, sys, json, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import MU, R_STAR, N_GRID, R_MIN, C_CALIB, PHI_GOLD
from lib.vacuum_phi import solve, extract_phi2
from lib.dirac_formT import find_eigenvalues

pi = np.pi
me = 0.511  # MeV

# Scan parameters
phi0_values = np.concatenate([
    np.linspace(-0.01, -0.10, 10),
    np.linspace(-0.12, -0.30, 10),
    np.linspace(-0.32, -0.60, 15),
    np.linspace(-0.62, -0.81, 15),
])
kappas = [-1, +1, -2, +2, -3, +3]

print(f"Emergence spectrum scan: {len(phi0_values)} phi_0 values x {len(kappas)} kappa channels")
print(f"r* = {R_STAR}, mu^2 = pi/3, C = {C_CALIB:.2f} MeV")
print(f"Sub-cavity r_b = r*/phi_gold = {R_STAR/PHI_GOLD:.4f}")
print()

results = []
t0 = time.time()

for i, phi0_t in enumerate(phi0_values):
    entry = {'phi0': float(phi0_t), 'eigenvalues': {}}

    try:
        r, phi, J, phip, overflow = solve(phi0_t, R_STAR, N_GRID, MU, r_min=R_MIN)
        if overflow is not None:
            entry['status'] = f'overflow at r={overflow}'
            results.append(entry)
            continue

        e2phi = np.exp(np.clip(2 * phi, -400, 400))
        phi2 = extract_phi2(r, phi, phi0_t)
        entry['phi_min'] = float(phi[-1])

        for kappa in kappas:
            ev = find_eigenvalues(kappa, r, phip, e2phi, phi0=phi0_t, phi2=phi2,
                                  E_min=0.01, E_max=30.0, n_scan=5000, n_modes=5)
            masses = [float(e * C_CALIB) for e in ev]
            entry['eigenvalues'][str(kappa)] = {
                'E': [float(e) for e in ev],
                'mass_MeV': masses,
                'n_modes': len(ev),
            }

        entry['status'] = 'ok'

        # Print progress
        n_total = sum(entry['eigenvalues'][str(k)]['n_modes'] for k in kappas)
        elapsed = time.time() - t0
        print(f"[{i+1}/{len(phi0_values)}] phi0={phi0_t:.4f}: {n_total} total modes, "
              f"phi_min={phi[-1]:.4f}, elapsed={elapsed:.0f}s")
        sys.stdout.flush()

    except Exception as e:
        entry['status'] = f'error: {str(e)}'

    results.append(entry)

# Summary
print()
print("=" * 60)
print("EMERGENCE SUMMARY")
print("=" * 60)

# Find threshold phi_0 where first eigenvalue appears for each kappa
for kappa in kappas:
    threshold = None
    for entry in results:
        if entry.get('status') != 'ok':
            continue
        if entry['eigenvalues'][str(kappa)]['n_modes'] > 0:
            threshold = entry['phi0']
            first_mass = entry['eigenvalues'][str(kappa)]['mass_MeV'][0]
            break
    if threshold:
        print(f"  kappa={kappa:+d}: first appears at phi0 = {threshold:.4f}, mass = {first_mass:.1f} MeV")
    else:
        print(f"  kappa={kappa:+d}: no eigenvalues found in scan range")

# At locked phi_0, print full spectrum
print()
locked = [e for e in results if abs(e['phi0'] - (-0.809)) < 0.01]
if locked:
    entry = locked[-1]
    print(f"At phi0 = {entry['phi0']:.4f} (near locked):")
    for kappa in kappas:
        masses = entry['eigenvalues'][str(kappa)]['mass_MeV']
        print(f"  kappa={kappa:+d}: {[f'{m:.1f}' for m in masses[:3]]} MeV")

# Re-emergence fractions
print()
print("Re-emergence fractions (from sub-cavity structure):")
print(f"  He-4 (sub-cavity volume): (r_b/r*)^3 = 1/phi_gold^3 = sqrt(5)-2 = {np.sqrt(5)-2:.6f}")
print(f"  H (remainder):           3-sqrt(5) = {3-np.sqrt(5):.6f}")
print(f"  D/H (EM branching):      alpha/(9*pi^3) = {1/137.036/(9*pi**3):.4e}")

# Save
out_dir = os.path.join(os.path.dirname(__file__), 'outputs')
os.makedirs(out_dir, exist_ok=True)
outpath = os.path.join(out_dir, 'emergence_spectrum.json')
with open(outpath, 'w') as f:
    json.dump({'scan': results, 'parameters': {
        'r_star': R_STAR, 'mu2': float(MU**2), 'C_MeV': float(C_CALIB),
        'phi_gold': float(PHI_GOLD), 'n_grid': N_GRID,
    }}, f, indent=2)
print(f"\nSaved: {outpath}")
print(f"Total time: {time.time()-t0:.0f}s")
