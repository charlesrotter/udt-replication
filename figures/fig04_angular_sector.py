#!/usr/bin/env python3
# Generates: Figure 4 in manuscript section 10
# Schematic: angular quantum number assignments for the four particles
"""
fig04 -- Angular sector schematic.

Shows how the Diophantine solution (j=1/2, l=1, |kappa_max|=3)
maps to the four angular particles (e, pi, mu, p) via
multiplicity-weighted pi-power assignments.

This is a conceptual/schematic figure using text annotations.

Saves: manuscript/figures/fig04_angular_sector.{pdf,png}
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# === Style setup ===
plt.style.use('default')
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif', 'Computer Modern Roman'],
    'mathtext.fontset': 'cm',
    'font.size': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
})

CB_BLUE = '#0072B2'
CB_ORANGE = '#E69F00'
CB_GREEN = '#009E73'
CB_RED = '#D55E00'

# === Create figure ===
fig, ax = plt.subplots(figsize=(3.4, 4.5))
ax.set_xlim(0, 10)
ax.set_ylim(0, 12)
ax.set_aspect('equal')
ax.axis('off')

# --- Title box ---
title_box = FancyBboxPatch((0.5, 10.5), 9, 1.2, boxstyle="round,pad=0.2",
                            facecolor='#E3F2FD', edgecolor='0.3', linewidth=1)
ax.add_patch(title_box)
ax.text(5, 11.1, r'Diophantine: $\,(2j\!+\!1)^2(2l\!+\!1)(2|\kappa_{\max}|\!+\!1)'
        r' = \binom{2l+2|\kappa_{\max}|+1}{2l+1}$',
        ha='center', va='center', fontsize=7.5, fontweight='bold')

# --- Quantum number solution box ---
qn_box = FancyBboxPatch((1.5, 8.8), 7, 1.3, boxstyle="round,pad=0.2",
                          facecolor='#FFF3E0', edgecolor='0.3', linewidth=0.8)
ax.add_patch(qn_box)
ax.text(5, 9.75, r'Unique solution: $j=\frac{1}{2},\; l=1,\; |\kappa_{\max}|=3$',
        ha='center', va='center', fontsize=9, fontweight='bold')
ax.text(5, 9.2, r'Multiplicities: $2j\!+\!1=2,\;\; 2l\!+\!1=3,\;\;'
        r' 2|\kappa_{\max}|\!-\!1=5,\;\; 2|\kappa_{\max}|\!+\!1=7$',
        ha='center', va='center', fontsize=7)

# --- Arrow from QN to particles ---
ax.annotate('', xy=(5, 7.8), xytext=(5, 8.7),
            arrowprops=dict(arrowstyle='->', color='0.4', lw=1.2))
ax.text(5, 8.25, 'mass formulas', ha='center', va='center', fontsize=7,
        style='italic', color='0.4')

# --- Four particle boxes ---
particles = [
    {
        'name': r'$e$ (electron)', 'color': CB_BLUE,
        'mass': r'$m_e$  (anchor)',
        'assign': r'grade 0: $m_e$',
        'value': '0.511 MeV',
        'x': 1.3, 'y': 5.8,
    },
    {
        'name': r'$\pi^0$ (pion)', 'color': CB_ORANGE,
        'mass': r'$C(9,3)\cdot\pi\cdot m_e$',
        'assign': r'$j \to \pi$: grade 3',
        'value': '134.7 MeV',
        'x': 5.7, 'y': 5.8,
    },
    {
        'name': r'$\mu$ (muon)', 'color': CB_GREEN,
        'mass': r'$\frac{20\pi^3}{3} m_e$',
        'assign': r'$l \to \mu$: $\pi^{2l+1}$',
        'value': '106.5 MeV',
        'x': 1.3, 'y': 3.0,
    },
    {
        'name': r'$p$ (proton)', 'color': CB_RED,
        'mass': r'$6\pi^5\, m_e$',
        'assign': r'$|\kappa_{\max}| \to p$: $\pi^{2|\kappa_{\max}|-1}$',
        'value': '938.5 MeV',
        'x': 5.7, 'y': 3.0,
    },
]

box_w, box_h = 3.5, 2.2
for p in particles:
    box = FancyBboxPatch((p['x'] - 0.3, p['y'] - 0.3), box_w, box_h,
                          boxstyle="round,pad=0.15",
                          facecolor='white', edgecolor=p['color'],
                          linewidth=1.5)
    ax.add_patch(box)
    cx = p['x'] - 0.3 + box_w / 2
    ax.text(cx, p['y'] + 1.5, p['name'],
            ha='center', va='center', fontsize=9, fontweight='bold',
            color=p['color'])
    ax.text(cx, p['y'] + 0.85, p['mass'],
            ha='center', va='center', fontsize=8)
    ax.text(cx, p['y'] + 0.25, p['assign'],
            ha='center', va='center', fontsize=7, style='italic',
            color='0.4')
    ax.text(cx, p['y'] - 0.1, p['value'],
            ha='center', va='center', fontsize=7, color='0.5')

# --- Bottom: base space note ---
ax.text(5, 0.7, r'Base space $\mathbb{R}^9$:  $n = 2l + 2|\kappa_{\max}| + 1 = 9$',
        ha='center', va='center', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#F3E5F5',
                  edgecolor='0.6', alpha=0.9))
ax.text(5, 0.15, r'$\Lambda^k(\mathbb{R}^9)$: exterior algebra determines '
        r'$\alpha_{\rm EM}$, neutrino sector',
        ha='center', va='center', fontsize=6.5, color='0.5')

fig.tight_layout()

# === Save ===
savefig(fig, 'fig04_angular_sector')
plt.close(fig)
print("fig04_angular_sector: done")
