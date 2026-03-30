#!/usr/bin/env python3
"""
fig03 -- Angular sector mass predictions: deviation from experiment.

Five particles from (j, l, |kappa_max|) = (1/2, 1, 3) and m_e.
Horizontal lollipop plot showing (m_pred - m_PDG) / m_PDG in percent.

Saves: manuscript/figures/fig03_mass_ratios.{pdf,png}
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from lib.constants import M_E, PDG_MASSES
from lib.angular_integrals import angular_mass_formulas
from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# === Publication style ===
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['CMU Serif', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'font.size': 9,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 9,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'xtick.minor.width': 0.3,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'savefig.dpi': 300,
})

# Okabe-Ito colorblind-safe
CB_BLUE   = '#0072B2'
CB_ORANGE = '#E69F00'
CB_GREEN  = '#009E73'
CB_RED    = '#D55E00'
CB_PURPLE = '#CC79A7'

# === Compute ===
masses = angular_mass_formulas()

# Bottom to top
particles = ['tau', 'proton', 'muon', 'pion', 'electron']
ylabels = [
    r'$\tau$' + '\n' + r'{\fontsize{6}{6}\selectfont Koide}',
    r'$p$' + '\n' + r'$6\pi^5 m_e$',
    r'$\mu$',
    r'$\pi^0$' + '\n' + r'$84\pi m_e$',
    r'$e$' + '\n' + r'anchor',
]
# Use simple labels + separate formula text
particle_labels = [r'$\tau$', r'$p$', r'$\mu$', r'$\pi^0$', r'$e$']
formula_labels = [
    r'Koide $Z_3$',
    r'$6\pi^5 m_e$',
    r'$\frac{20}{3}\pi^3 m_e$',
    r'$84\pi\, m_e$',
    'anchor',
]
colors = [CB_PURPLE, CB_RED, CB_GREEN, CB_ORANGE, CB_BLUE]

errors_pct = []
for p in particles:
    m = masses[p]
    errors_pct.append((m['predicted'] / m['pdg'] - 1.0) * 100)

# === Create figure ===
fig, ax = plt.subplots(figsize=(3.4, 2.8))

y = np.arange(len(particles))

# ±0.1% band
ax.axvspan(-0.1, 0.1, color='#E8F4E8', alpha=0.6, zorder=0)
ax.axvline(0, color='0.55', linewidth=0.5, zorder=1)

# Lollipop: thin line + marker
for i, (err, col) in enumerate(zip(errors_pct, colors)):
    ax.plot([0, err], [i, i], '-', color=col, linewidth=0.9, zorder=3,
            solid_capstyle='round')
    ax.plot(err, i, 'o', color=col, markersize=6.5, markeredgecolor='0.15',
            markeredgewidth=0.5, zorder=5)

# Error + formula annotations above each row
for i, (err, form) in enumerate(zip(errors_pct, formula_labels)):
    # Error percentage: offset above marker
    ax.text(err, i + 0.32, f'{err:+.3f}%', fontsize=6.5,
            ha='center', va='bottom', color='0.15', fontweight='bold')
    # Formula: right-aligned near right edge
    ax.text(0.97, i, form, fontsize=6, ha='right', va='center',
            color='0.45', transform=ax.get_yaxis_transform())

ax.set_yticks(y)
ax.set_yticklabels(particle_labels)
ax.set_xlabel(r'$(m_{\rm pred} - m_{\rm PDG})\,/\,m_{\rm PDG}$ (%)')
ax.set_xlim(-0.13, 0.13)
ax.set_ylim(-0.6, len(particles) - 0.3)

# ±0.1% label
ax.text(0.098, len(particles) - 0.45, r'$\pm 0.1\%$', fontsize=5.5,
        ha='right', va='top', color='#388E3C', alpha=0.8)

# Title
ax.set_title(r'Angular sector:  $(j,\,\ell,\,|\kappa_{\max}|) = (1/2,\,1,\,3)$',
             fontsize=8.5, pad=6)

ax.tick_params(axis='y', length=0)
fig.tight_layout()

# === Save ===
savefig(fig, 'fig03_mass_ratios')
plt.close(fig)
print("fig03_mass_ratios: done")
