#!/usr/bin/env python3
# Generates: Figure 3 in manuscript section 9
# Predicted vs observed mass ratios for angular sector particles
"""
fig03 -- Angular sector mass ratios: m_predicted / m_observed.

Includes electron (anchor), pion, muon, proton, and tau (Koide Z3).
Horizontal line at 1.0 (perfect prediction).

Saves: manuscript/figures/fig03_mass_ratios.{pdf,png}
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import M_E, PDG_MASSES
from lib.angular_integrals import angular_mass_formulas
from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# === Style setup ===
plt.style.use('default')
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif', 'Computer Modern Roman'],
    'mathtext.fontset': 'cm',
    'font.size': 10,
    'axes.labelsize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'axes.linewidth': 0.8,
})

CB_BLUE = '#0072B2'
CB_ORANGE = '#E69F00'
CB_GREEN = '#009E73'
CB_RED = '#D55E00'
CB_PURPLE = '#CC79A7'

# === Compute mass ratios ===
masses = angular_mass_formulas()

particles = ['electron', 'pion', 'muon', 'proton', 'tau']
labels = [r'$e$', r'$\pi^0$', r'$\mu$', r'$p$', r'$\tau$']
formulas = [
    r'anchor',
    r'$84\pi\, m_e$',
    r'$\frac{20\pi^3}{3}\, m_e$',
    r'$6\pi^5\, m_e$',
    r'Koide $Z_3$',
]
colors = [CB_BLUE, CB_ORANGE, CB_GREEN, CB_RED, CB_PURPLE]

ratios = []
for p in particles:
    m = masses[p]
    ratios.append(m['predicted'] / m['pdg'])

ratios = np.array(ratios)
errors_pct = (ratios - 1.0) * 100

# === Create figure ===
fig, ax = plt.subplots(figsize=(3.8, 3.0))

x = np.arange(len(particles))

# Perfect prediction line
ax.axhline(1.0, color='0.7', linewidth=0.8, linestyle='--', zorder=1)

# 0.1% band
ax.axhspan(0.999, 1.001, color='#E8F5E9', alpha=0.6, zorder=0,
            label=r'$\pm 0.1\%$ band')

# Data points
ax.bar(x, ratios - 1.0, bottom=1.0, width=0.55,
       color=colors,
       edgecolor='0.3', linewidth=0.5, zorder=3)

# Add error text on each bar
for i, (r, ep) in enumerate(zip(ratios, errors_pct)):
    y_off = 0.0004 if r >= 1.0 else -0.0004
    va = 'bottom' if r >= 1.0 else 'top'
    ax.text(i, r + y_off, f'{ep:+.3f}%', ha='center', va=va, fontsize=6.5,
            fontweight='bold')

# Formula annotations below
for i, f in enumerate(formulas):
    ax.text(i, 0.9965, f, ha='center', va='top', fontsize=5.5,
            style='italic', rotation=0)

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel(r'$m_{\rm pred} \,/\, m_{\rm PDG}$')
ax.set_ylim(0.996, 1.004)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.legend(loc='upper right', frameon=True, fancybox=False,
          edgecolor='0.7', framealpha=0.9, fontsize=6.5)

ax.set_title('Angular sector mass predictions\n'
             r'$(j,\,\ell,\,|\kappa_{\max}|) = (1/2,\,1,\,3)$',
             fontsize=9, pad=6)

fig.tight_layout()

# === Save ===
savefig(fig, 'fig03_mass_ratios')
plt.close(fig)
print("fig03_mass_ratios: done")
