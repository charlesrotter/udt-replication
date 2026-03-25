#!/usr/bin/env python3
"""04 -- Bridge identity products: source(kappa) * E1(kappa).

SOURCE: data/generated/02_eigenvalues.json, data/generated/03_sources.json
GENERATES: data/generated/04_bridge_identities.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2

This is pure data analysis -- no ODE integration needed.

Bridge products:
  |kappa|=1: source * E1 = sqrt(2)/6
  |kappa|=2: source * E1 = pi/6
  |kappa|=3: source * E1 = sqrt(2)/4

Squared products:
  |kappa|=1: source^2 * E1^2 = 1/18
  |kappa|=2: source^2 * E1^2 = pi^2/36
  |kappa|=3: source^2 * E1^2 = 1/8

Completion class: A (pure algebra from A-class inputs).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.utils import save_results, load_results, pct_error, gate_check

# === Load eigenvalues and sources ===
ev_data = load_results('02_eigenvalues.json')
src_data = load_results('03_sources.json')

eigenvalues = ev_data['eigenvalues']
sources = src_data['sources']

# === Algebraic targets for bridge products ===
bridge_targets = {
    1: {'se': np.sqrt(2) / 6, 'label': 'sqrt(2)/6',
        's2e2': 1.0 / 18.0, 's2e2_label': '1/18'},
    2: {'se': np.pi / 6, 'label': 'pi/6',
        's2e2': np.pi**2 / 36.0, 's2e2_label': 'pi^2/36'},
    3: {'se': np.sqrt(2) / 4, 'label': 'sqrt(2)/4',
        's2e2': 1.0 / 8.0, 's2e2_label': '1/8'},
}

# === Compute bridge products for each |kappa| ===
bridge_results = []

for ak in [1, 2, 3]:
    for kappa in [-ak, +ak]:
        ks = str(kappa)
        evals = eigenvalues.get(ks, [])
        srcs = sources.get(ks, [])

        if not evals or not srcs:
            continue

        E1 = evals[0]
        S1 = srcs[0]['source']

        # Bridge identity: |S| * E1 = target for each |kappa|.
        # For negative kappa, source is positive (G-dominant Frobenius).
        # For positive kappa, source is negative (F-dominant).
        se_product = abs(S1) * E1
        s2e2_product = S1**2 * E1**2

        tgt = bridge_targets[ak]
        se_err = pct_error(se_product, tgt['se'])
        s2e2_err = pct_error(s2e2_product, tgt['s2e2'])

        bridge_results.append({
            'kappa': kappa,
            'abs_kappa': ak,
            'E1': E1,
            'source': S1,
            'source_abs': abs(S1),
            'SE_product': se_product,
            'SE_target': tgt['se'],
            'SE_label': tgt['label'],
            'SE_pct_error': se_err,
            'S2E2_product': s2e2_product,
            'S2E2_target': tgt['s2e2'],
            'S2E2_label': tgt['s2e2_label'],
            'S2E2_pct_error': s2e2_err,
        })

# === Gate: bridge products match for canonical (kappa<0) channels ===
# The algebraic bridge identities are derived for the kappa<0 (G-dominant) sector.
# Positive kappa has different E1 values, so the same identities do not apply.
neg_kappa_se_pass = True
neg_kappa_s2e2_pass = True
for br in bridge_results:
    if br['kappa'] < 0:
        if abs(br['SE_pct_error']) > 1.0:
            neg_kappa_se_pass = False
        if abs(br['S2E2_pct_error']) > 1.0:
            neg_kappa_s2e2_pass = False

gate_check("G1", 0.0 if neg_kappa_se_pass else 1.0, 0.5,
           "|S|*E1 products (kappa<0) match algebraic targets to < 1%")
gate_check("G2", 0.0 if neg_kappa_s2e2_pass else 1.0, 0.5,
           "S^2*E1^2 products (kappa<0) match to < 1%")

# === Save results ===
results = {
    'bridge_products': bridge_results,
    'targets': {
        str(ak): {
            'SE': bridge_targets[ak]['se'],
            'SE_label': bridge_targets[ak]['label'],
            'S2E2': bridge_targets[ak]['s2e2'],
            'S2E2_label': bridge_targets[ak]['s2e2_label'],
        }
        for ak in [1, 2, 3]
    },
    'kappa_neg_SE_match': neg_kappa_se_pass,
    'kappa_neg_S2E2_match': neg_kappa_s2e2_pass,
}
save_results('04_bridge_identities.json', results)

# === Summary ===
print("Bridge identity verification:")
for br in bridge_results:
    k = br['kappa']
    print(f"  kappa={k:+d}: S*E1={br['SE_product']:.10f} vs {br['SE_label']}={br['SE_target']:.10f} ({br['SE_pct_error']:+.4f}%)")
print(f"Squared products:")
# Print one per |kappa| (signs should agree)
for ak in [1, 2, 3]:
    entries = [b for b in bridge_results if b['abs_kappa'] == ak]
    if entries:
        b = entries[0]
        print(f"  |kappa|={ak}: S^2*E1^2={b['S2E2_product']:.10f} vs {b['S2E2_label']}={b['S2E2_target']:.10f} ({b['S2E2_pct_error']:+.4f}%)")
