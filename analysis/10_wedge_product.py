#!/usr/bin/env python3
"""10 -- Exterior algebra structure: pion (grade 3) to neutrino (grade 6).

SOURCE: lib/constants.py (angular quantum numbers)
GENERATES: data/generated/10_wedge_product.json

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  sqrt(-g) = c r^2 sin(theta)  (phi-independent)
  Tetrad: e^a_mu = diag(e^{-phi}, e^{+phi}, r, r sin(theta))

The angular Diophantine selects a 9-dimensional base space (2l + 2|kappa_max| + 1 = 9).
The pion lives as a 3-form (Lambda^3): C(9,3) = 84 components.
The alpha boson lives as a 2-form (Lambda^2): C(9,2) = 36 components.
The neutrino lives as (Lambda^2)^{otimes 3}: 36^3 = 46656 components.
Hodge duality on Lambda^9: *Lambda^3 = Lambda^6 (grade 3 -> grade 6).

Completion class: A (algebraic, no ODE).
"""
import os
import sys
import numpy as np
from math import comb, factorial

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import (J_HALF, ELL, KAPPA_MAX,
                            MULT_2J1, MULT_2L1, MULT_2KM1, MULT_2KP1)
from lib.utils import save_results, gate_check

# === Base space dimension ===
# From Diophantine: n = 2l + 2|kappa_max| + 1
n_base = MULT_2L1 + MULT_2KP1  # 3 + 7 = 10? No: 2l+1=3, 2|kmax|+1=7 -> sum
# Actually: the base space is n = 2l + 2|kappa_max| + 1 = 2*1 + 2*3 + 1 = 9
n_base_direct = 2 * ELL + 2 * KAPPA_MAX + 1  # = 9
assert n_base_direct == 9, f"Expected base space dim 9, got {n_base_direct}"

n = n_base_direct  # 9

# === Grade-k exterior algebras ===
# Lambda^k on R^n has dimension C(n, k)

dim_Lambda1 = comb(n, 1)   # 9
dim_Lambda2 = comb(n, 2)   # 36
dim_Lambda3 = comb(n, 3)   # 84
dim_Lambda4 = comb(n, 4)   # 126
dim_Lambda5 = comb(n, 5)   # 126
dim_Lambda6 = comb(n, 6)   # 84
dim_Lambda7 = comb(n, 7)   # 36
dim_Lambda8 = comb(n, 8)   # 9
dim_Lambda9 = comb(n, 9)   # 1

# === Pion: grade 3 ===
# C(9,3) = 84 = Sym^3(R^7) dimension check
dim_pion = dim_Lambda3  # 84
sym3_R7 = comb(7 + 3 - 1, 3)  # C(9,3) = 84 (Sym^3 on R^7 = C(7+3-1,3))
assert dim_pion == sym3_R7 == 84, f"Pion dimension mismatch: {dim_pion} vs {sym3_R7}"

# === Alpha: grade 2 ===
dim_alpha = dim_Lambda2  # 36

# === Neutrino: (Lambda^2)^{otimes 3} ===
dim_neutrino_tensor = dim_Lambda2**3  # 36^3 = 46656

# === Hodge duality verification ===
# On Lambda^n, *: Lambda^k -> Lambda^{n-k}
# *Lambda^3 = Lambda^{9-3} = Lambda^6
# dim(Lambda^3) should equal dim(Lambda^6) by Hodge symmetry
hodge_grade_in = 3
hodge_grade_out = n - hodge_grade_in  # 6
hodge_dim_match = (dim_Lambda3 == dim_Lambda6)
assert hodge_dim_match, f"Hodge dims don't match: {dim_Lambda3} vs {dim_Lambda6}"

# === Wedge product computation ===
# Lambda^2 ^ Lambda^2 ^ Lambda^2 -> Lambda^6
# The angular determinant of this triple wedge is a specific combinatorial factor.
#
# For a 2-form alpha in Lambda^2(R^n), alpha ^ alpha ^ alpha in Lambda^6(R^n)
# is determined by the Pfaffian structure.
#
# The number of independent triple-wedge components:
# Given a basis {e_i ^ e_j} for Lambda^2, the triple wedge maps
# (Lambda^2)^{otimes 3} -> Lambda^6.
# The image is Sym^3(Lambda^2) restricted to Lambda^6.
#
# Counting: the angular determinant is the number of distinct 6-element subsets
# of {1,...,9} that can be partitioned into three disjoint pairs.
# That is: C(9,6) * (number of perfect matchings on 6 vertices) / normalization

# Perfect matchings on 6 vertices: (2k)!! / 2^k = 6! / (2^3 * 3!) = 15
n_perfect_matchings_6 = factorial(6) // (2**3 * factorial(3))  # = 15
assert n_perfect_matchings_6 == 15

# C(9,6) = 84 subsets, each with 15 matchings
# But the angular determinant is C(9,6) = 84 (same as Lambda^3 by Hodge)
angular_det = comb(n, 6)  # 84

# === Verification: Hodge star maps grade 3 -> grade 6 ===
# The Hodge dual of the 84-component 3-form gives an 84-component 6-form
# The triple wedge Lambda^2 ^ Lambda^2 ^ Lambda^2 also lands in Lambda^6
# This is the bridge: pion (grade 3) <-> neutrino (grade 6) via Hodge

# === Dimension chain verification ===
# Poincare duality: C(n,k) = C(n, n-k)
poincare_check = all(comb(n, k) == comb(n, n - k) for k in range(n + 1))

# Total exterior algebra dimension: sum_k C(n,k) = 2^n
total_dim = sum(comb(n, k) for k in range(n + 1))
assert total_dim == 2**n, f"Total dim {total_dim} != 2^{n} = {2**n}"

# Euler characteristic: sum_k (-1)^k C(n,k) = 0 for n > 0
euler = sum((-1)**k * comb(n, k) for k in range(n + 1))
assert euler == 0, f"Euler characteristic {euler} != 0"

# === Physical interpretation ===
# The 36^3 = 46656 tensor product counts the full (unconstrained) neutrino state space.
# The actual neutrino lives in the antisymmetric subspace Lambda^6, which has dim 84.
# The ratio 46656/84 = 555.43... gives the combinatorial reduction factor.
reduction_factor = dim_neutrino_tensor / dim_Lambda6

# Damping ratio from the manuscript: C(9,3)/C(9,4) = 84/126 = 2/3
damping_ratio = dim_Lambda3 / dim_Lambda4
assert abs(damping_ratio - 2.0 / 3.0) < 1e-14

# Even/odd ratio: C(9,2)/C(9,3) = 36/84 = 3/7
even_odd_ratio = dim_Lambda2 / dim_Lambda3
assert abs(even_odd_ratio - 3.0 / 7.0) < 1e-14

# === Gate: all combinatorial identities consistent ===
gate_passed = gate_check("G1", float(not hodge_dim_match), 0.5,
                          "Hodge dim(Lambda^3) = dim(Lambda^6)")
gate_poincare = gate_check("G2", float(not poincare_check), 0.5,
                            "Poincare duality for all grades")
gate_euler = gate_check("G3", float(abs(euler)), 0.5,
                         "Euler characteristic = 0")

# === Save results ===
results = {
    'base_space_dim': n,
    'base_space_derivation': '2l + 2|kappa_max| + 1 = 2*1 + 2*3 + 1 = 9',
    'exterior_algebra': {
        'total_dim': total_dim,
        'grade_dims': {str(k): comb(n, k) for k in range(n + 1)},
    },
    'pion': {
        'grade': 3,
        'dim_Lambda3': dim_Lambda3,
        'sym3_R7_check': sym3_R7,
        'interpretation': 'C(9,3) = 84 = dim(3-form on R^9)',
    },
    'alpha': {
        'grade': 2,
        'dim_Lambda2': dim_Lambda2,
        'interpretation': 'C(9,2) = 36 = dim(2-form on R^9)',
    },
    'neutrino': {
        'tensor_grade': '(Lambda^2)^{otimes 3}',
        'tensor_dim': dim_neutrino_tensor,
        'antisymmetric_dim': dim_Lambda6,
        'reduction_factor': reduction_factor,
        'interpretation': '36^3 = 46656 -> Lambda^6 = 84 via antisymmetrization',
    },
    'hodge_dual': {
        'map': '*Lambda^3 -> Lambda^{n-3} = Lambda^6',
        'dim_source': dim_Lambda3,
        'dim_target': dim_Lambda6,
        'dims_match': hodge_dim_match,
        'interpretation': 'grade 3 (pion) <-> grade 6 (neutrino) via Hodge star',
    },
    'angular_determinant': {
        'value': angular_det,
        'perfect_matchings_6': n_perfect_matchings_6,
        'interpretation': 'C(9,6) = 84 independent triple-wedge components',
    },
    'derived_ratios': {
        'damping_2_3': {
            'value': damping_ratio,
            'formula': 'C(9,3)/C(9,4) = 84/126 = 2/3',
        },
        'even_odd_3_7': {
            'value': even_odd_ratio,
            'formula': 'C(9,2)/C(9,3) = 36/84 = 3/7',
        },
    },
    'gates': {
        'G1_hodge_match': gate_passed,
        'G2_poincare': gate_poincare,
        'G3_euler': gate_euler,
    },
}
save_results('10_wedge_product.json', results)

# === Summary ===
print(f"Exterior algebra on R^{n}: total dim = {total_dim} = 2^{n}")
print(f"  Pion: Lambda^3 = {dim_Lambda3}, Alpha: Lambda^2 = {dim_Lambda2}")
print(f"  Neutrino: (Lambda^2)^3 = {dim_neutrino_tensor}, Hodge target Lambda^6 = {dim_Lambda6}")
print(f"  Damping = {damping_ratio:.4f} = 2/3, Even/odd = {even_odd_ratio:.4f} = 3/7")
print(f"  All gates: {'PASS' if all([gate_passed, gate_poincare, gate_euler]) else 'FAIL'}")
