#!/usr/bin/env python3
# Generates: Figure 9 in manuscript section 15
# Unified connection map: metric -> all predictions
"""
fig09 -- Unified connection map.

Central node: ds^2 metric
Branching to sectors via the Diophantine primes (2,3,5,7).
Shows: masses, alpha, neutrino, PMNS, CMB, BAO as leaf predictions.

Uses manual matplotlib drawing (no networkx dependency required).

Saves: manuscript/figures/fig09_connection_map.{pdf,png}
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

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

CB_BLUE = '#0072B2'
CB_ORANGE = '#E69F00'
CB_GREEN = '#009E73'
CB_RED = '#D55E00'
CB_PURPLE = '#CC79A7'
CB_CYAN = '#56B4E9'

# === Create figure ===
fig, ax = plt.subplots(figsize=(7.0, 5.5))
ax.set_xlim(-5.5, 5.5)
ax.set_ylim(-4.5, 4.5)
ax.set_aspect('equal')
ax.axis('off')

def draw_node(ax, x, y, text, size=0.6, color='white', ec='0.3', fontsize=8,
              bold=False, shape='round'):
    """Draw a node (rounded box) at (x,y)."""
    w = size * 2.2
    h = size * 1.2
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                          boxstyle=f"round,pad=0.15",
                          facecolor=color, edgecolor=ec, linewidth=1.0)
    ax.add_patch(box)
    fw = 'bold' if bold else 'normal'
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
            fontweight=fw, zorder=10)

def draw_edge(ax, x1, y1, x2, y2, label='', color='0.5', lw=1.0):
    """Draw a connection line between two nodes."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                connectionstyle='arc3,rad=0.0'))
    if label:
        mx, my = 0.5*(x1+x2), 0.5*(y1+y2)
        ax.text(mx, my + 0.15, label, ha='center', va='bottom',
                fontsize=6, color=color, style='italic')

# === Central node: the metric ===
draw_node(ax, 0, 3, r'$ds^2 = -e^{-2\phi}c^2dt^2 + e^{2\phi}dr^2 + r^2d\Omega^2$',
          size=2.2, color='#E3F2FD', ec=CB_BLUE, fontsize=9, bold=True)

# === First level: ODE + Angular ===
draw_node(ax, -2.5, 1.2, 'Scalar ODE\n' + r'$\nabla^2_g\phi=\mu^2\phi$',
          size=1.0, color='#FFF3E0', ec=CB_ORANGE, fontsize=7)
draw_node(ax, 2.5, 1.2, 'Angular\n' + r'$S^2$ harmonics',
          size=1.0, color='#E8F5E9', ec=CB_GREEN, fontsize=7)

draw_edge(ax, -0.5, 2.3, -2.0, 1.8, r'$\phi(r)$', CB_ORANGE, 1.2)
draw_edge(ax, 0.5, 2.3, 2.0, 1.8, 'Diophantine', CB_GREEN, 1.2)

# === Second level: branching ===
# Left branch: micro + cosmo from ODE
draw_node(ax, -4.0, -0.5, 'Microphysics\nDirac on' + r' $\phi$',
          size=1.0, color='#FCE4EC', ec=CB_RED, fontsize=7)
draw_node(ax, -1.0, -0.5, 'Cosmology\n' + r'$\phi(r)$ polynomial',
          size=1.0, color='#F3E5F5', ec=CB_PURPLE, fontsize=7)

draw_edge(ax, -2.8, 0.6, -3.8, 0.1, r'$E_n$', CB_RED, 1.0)
draw_edge(ax, -2.2, 0.6, -1.2, 0.1, r'$\kappa,\beta,\gamma$', CB_PURPLE, 1.0)

# Right branch: masses + alpha from angular
draw_node(ax, 1.5, -0.5, r'$\alpha_{\rm EM}$' + '\n' + r'$36\pi/I_2$',
          size=0.8, color='#E3F2FD', ec=CB_BLUE, fontsize=7)
draw_node(ax, 4.0, -0.5, 'Mass formulas\n' + r'$m_e,\pi,\mu,p$',
          size=1.0, color='#E8F5E9', ec=CB_GREEN, fontsize=7)

draw_edge(ax, 2.0, 0.6, 1.5, 0.1, r'$I_2$', CB_BLUE, 1.0)
draw_edge(ax, 3.0, 0.6, 3.8, 0.1, r'$(j,l,|\kappa|)$', CB_GREEN, 1.0)

# === Third level: leaf predictions ===
# From microphysics
draw_node(ax, -5.0, -2.5, 'Mass ladder\n' + r'$m_n = E_n C$',
          size=0.9, color='white', fontsize=6.5)
draw_node(ax, -3.0, -2.5, 'Hierarchy\n' + r'$H \approx 40$',
          size=0.9, color='white', fontsize=6.5)

draw_edge(ax, -4.2, -1.1, -4.8, -1.9, '', CB_RED, 0.7)
draw_edge(ax, -3.8, -1.1, -3.2, -1.9, '', CB_RED, 0.7)

# From cosmology
draw_node(ax, -1.5, -2.5, 'SNe Ia\n0.166 mag',
          size=0.8, color='white', fontsize=6.5)
draw_node(ax, 0.3, -2.5, 'BAO\n' + r'$D_V/r_d$',
          size=0.7, color='white', fontsize=6.5)

draw_edge(ax, -1.0, -1.1, -1.3, -1.9, '', CB_PURPLE, 0.7)
draw_edge(ax, -0.8, -1.1, 0.1, -1.9, '', CB_PURPLE, 0.7)

# From alpha
draw_node(ax, 1.5, -2.5, r'$\nu$ mass' + '\n' + r'$\alpha^3 m_e/4$',
          size=0.8, color='white', fontsize=6.5)

draw_edge(ax, 1.5, -1.1, 1.5, -1.9, '', CB_BLUE, 0.7)

# From mass formulas
draw_node(ax, 3.2, -2.5, 'PMNS\n' + r'$\frac{4}{13},\frac{4}{7},\frac{1}{45}$',
          size=0.9, color='white', fontsize=6.5)
draw_node(ax, 5.0, -2.5, r'$\Lambda^k$ algebra' + '\n' + 'grades',
          size=0.8, color='white', fontsize=6.5)

draw_edge(ax, 3.8, -1.1, 3.2, -1.9, '', CB_GREEN, 0.7)
draw_edge(ax, 4.2, -1.1, 4.8, -1.9, '', CB_GREEN, 0.7)

# === Bottom: CMB from cosmo + angular ===
draw_node(ax, -0.5, -3.8, r'CMB:  $\phi(r_*) = \ln(1+z_{\rm CMB}) = 7.004$',
          size=2.0, color='#FFFDE7', ec='0.5', fontsize=7)

draw_edge(ax, -1.5, -3.1, -1.0, -3.3, '', '0.6', 0.5)
draw_edge(ax, 0.3, -3.1, -0.2, -3.3, '', '0.6', 0.5)

# === Count annotations ===
ax.text(0, -4.3, r'1 metric $\,\to\,$ 0 free parameters $\,\to\,$ '
        r'masses, $\alpha$, PMNS, CMB, BAO, hierarchy',
        ha='center', va='center', fontsize=8, fontweight='bold',
        color='0.3')

fig.tight_layout()

# === Save ===
savefig(fig, 'fig09_connection_map')
plt.close(fig)
print("fig09_connection_map: done")
