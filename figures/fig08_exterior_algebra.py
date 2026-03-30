#!/usr/bin/env python3
# Generates: Figure 8 in manuscript section 14
# Exterior algebra grade structure: Lambda^2 -> Lambda^3 -> Lambda^6
"""
fig08 -- Exterior algebra grade structure schematic.

Shows the chain: Lambda^2 (alpha) -> Lambda^3 (pion) -> Lambda^6 (neutrino)
with Hodge dual arrow *Lambda^3 = Lambda^6.
Includes grade numbers, dimensions, and physical interpretations.

Saves: manuscript/figures/fig08_exterior_algebra.{pdf,png}
"""
import os
import sys
import numpy as np
from math import comb

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc

# === Style setup ===
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['CMU Serif', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'font.size': 9,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.3,
    'ytick.minor.width': 0.3,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'savefig.dpi': 300,
})

CB_BLUE   = '#0072B2'
CB_ORANGE = '#E69F00'
CB_GREEN  = '#009E73'
CB_RED    = '#D55E00'
CB_PURPLE = '#CC79A7'
CB_CYAN   = '#56B4E9'

n = 9  # base space dimension

# === Create figure ===
fig, ax = plt.subplots(figsize=(7.0, 3.5))
ax.set_xlim(-0.5, 10.5)
ax.set_ylim(-1.5, 4.5)
ax.set_aspect('equal')
ax.axis('off')

# --- Title ---
ax.text(5, 4.0, r'Exterior Algebra on $\mathbb{R}^9$:  '
        r'$\Lambda^k(\mathbb{R}^9)$,  dim $= \binom{9}{k}$',
        ha='center', va='center', fontsize=11, fontweight='bold')

# --- Grade boxes along bottom ---
grades = list(range(10))
dims = [comb(9, k) for k in grades]

# Physical assignments
phys = {
    0: r'$\mathbf{1}$ (scalar)',
    2: r'$\alpha_{\rm EM}$',
    3: r'$\pi$ (pion)',
    6: r'$\nu$ (neutrino)',
    9: r'$\mathbf{1}$ (volume)',
}

highlight_colors = {
    2: '#E3F2FD',
    3: '#FFF3E0',
    6: '#F3E5F5',
}

box_w = 0.85
box_h = 1.4
y_base = 0.5

for k in grades:
    x_c = k + 0.5
    fc = highlight_colors.get(k, '#F5F5F5')
    ec = '0.3' if k in phys else '0.7'
    lw = 1.2 if k in phys else 0.5

    box = FancyBboxPatch((x_c - box_w/2, y_base - box_h/2), box_w, box_h,
                          boxstyle="round,pad=0.08",
                          facecolor=fc, edgecolor=ec, linewidth=lw)
    ax.add_patch(box)

    # Grade label
    ax.text(x_c, y_base + 0.35, f'$k={k}$', ha='center', va='center',
            fontsize=7, fontweight='bold')

    # Dimension
    ax.text(x_c, y_base - 0.0, f'{dims[k]}', ha='center', va='center',
            fontsize=9, color=CB_BLUE if k in phys else '0.5')

    # Physical assignment
    if k in phys:
        ax.text(x_c, y_base - 0.45, phys[k], ha='center', va='center',
                fontsize=5.5, color='0.3')

# --- Wedge product arrow: Lambda^2 -> Lambda^3 ---
ax.annotate('', xy=(3.05, y_base + 0.8), xytext=(2.95, y_base + 0.8),
            arrowprops=dict(arrowstyle='->', color=CB_GREEN, lw=1.5))
# Broader arrow
ax.annotate(r'$\wedge$', xy=(3.5, y_base + 1.3), xytext=(2.5, y_base + 1.3),
            fontsize=10, color=CB_GREEN, ha='center', va='center',
            arrowprops=dict(arrowstyle='->', color=CB_GREEN, lw=1.5,
                            connectionstyle='arc3,rad=0'))

# --- Triple wedge arrow: Lambda^2^3 -> Lambda^6 ---
ax.annotate(r'$\alpha\wedge\alpha\wedge\alpha$',
            xy=(6.5, y_base + 1.3), xytext=(2.5, y_base + 2.3),
            fontsize=8, color=CB_PURPLE, ha='center', va='center',
            arrowprops=dict(arrowstyle='->', color=CB_PURPLE, lw=1.5,
                            connectionstyle='arc3,rad=-0.2'))

# --- Hodge dual arrow: *Lambda^3 = Lambda^6 ---
ax.annotate(r'Hodge: $\,\star\Lambda^3 = \Lambda^6$',
            xy=(6.5, y_base - 0.9), xytext=(3.5, y_base - 0.9),
            fontsize=8, color=CB_RED, ha='center', va='center',
            arrowprops=dict(arrowstyle='->', color=CB_RED, lw=1.5,
                            connectionstyle='arc3,rad=0.3'))

# --- Poincare duality brackets ---
pairs = [(0, 9), (1, 8), (2, 7), (3, 6), (4, 5)]
for k1, k2 in pairs:
    x1, x2 = k1 + 0.5, k2 + 0.5
    y_br = y_base - 0.85 - 0.2 * (k2 - k1 < 6)
    if k1 in phys or k2 in phys:
        ax.plot([x1, x1, x2, x2], [y_base - 0.72, y_br, y_br, y_base - 0.72],
                color='0.8', linewidth=0.5, zorder=0)

# --- Key relations annotation ---
ax.text(5, -1.1,
        r'$\binom{9}{2} = 36 \;\to\; 1/\alpha = 36\pi/I_2$'
        r'$\qquad\qquad$'
        r'$\binom{9}{3} = 84 \;\to\; m_\pi = 84\,\pi\, m_e$'
        r'$\qquad\qquad$'
        r'$36^3 = 46656 \;\to\; \nu$ tensor space',
        ha='center', va='center', fontsize=7,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFFDE7',
                  edgecolor='0.6', alpha=0.9))

fig.tight_layout()

# === Save ===
savefig(fig, 'fig08_exterior_algebra')
plt.close(fig)
print("fig08_exterior_algebra: done")
