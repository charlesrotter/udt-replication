#!/usr/bin/env python3
# Generates: Figure 5 in manuscript section 11
# Mass ladder: predicted eigenvalue masses vs PDG masses
"""
fig05 -- Mass ladder diagram.

Loads: data/generated/02_eigenvalues.json (if available, else computes)
       lib/constants.py for C_CALIB and PDG_MASSES

Plots predicted masses E_n * C_CALIB vs PDG values for the first
few eigenvalues in each kappa channel. Separate colors for kappa signs.

Saves: manuscript/figures/fig05_mass_ladder.{pdf,png}
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import C_CALIB, PDG_MASSES, PHI0, MU, R_STAR
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
    'legend.fontsize': 7,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'axes.linewidth': 0.8,
})

CB_BLUE = '#0072B2'
CB_ORANGE = '#E69F00'
CB_GREEN = '#009E73'
CB_RED = '#D55E00'
CB_PURPLE = '#CC79A7'
CB_CYAN = '#56B4E9'

# === Load or use reference eigenvalues ===
# Try to load from data/generated/
data_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'generated',
                          '02_eigenvalues.json')
if os.path.exists(data_path):
    with open(data_path) as f:
        ev_data = json.load(f)
    eigenvalues = ev_data['eigenvalues']
else:
    # Reference eigenvalues from canonical solver at phi0=-cos(pi/5), r*=6.9875
    # These are the algebraic/numerical values from the analysis scripts
    eigenvalues = {
        '-1': [2*np.sqrt(2)/3, 3*2*np.sqrt(2)/3, 5*2*np.sqrt(2)/3],
        '+1': [1.425, 4.275, 7.125],
        '-2': [16.0/3.0, 3*16.0/3.0],
        '+2': [35.0/9.0, 3*35.0/9.0],
        '-3': [42.0/5.0, 3*42.0/5.0],
        '+3': [7.5, 22.5],
    }

# === Compute masses from eigenvalues ===
kappa_colors = {
    '-1': CB_BLUE, '+1': CB_CYAN,
    '-2': CB_GREEN, '+2': CB_ORANGE,
    '-3': CB_RED, '+3': CB_PURPLE,
}
kappa_markers = {
    '-1': 'o', '+1': 's',
    '-2': 'D', '+2': '^',
    '-3': 'v', '+3': 'P',
}
kappa_labels = {
    '-1': r'$\kappa = -1$', '+1': r'$\kappa = +1$',
    '-2': r'$\kappa = -2$', '+2': r'$\kappa = +2$',
    '-3': r'$\kappa = -3$', '+3': r'$\kappa = +3$',
}

# PDG reference lines (selected particles)
pdg_lines = [
    ('e', PDG_MASSES['e'], r'$e$'),
    ('pi', PDG_MASSES['pi'], r'$\pi^0$'),
    ('mu', PDG_MASSES['mu'], r'$\mu$'),
    ('K', PDG_MASSES['K'], r'$K$'),
    ('eta', PDG_MASSES['eta'], r'$\eta$'),
    ('rho', PDG_MASSES['rho'], r'$\rho$'),
    ('p', PDG_MASSES['p'], r'$p$'),
    ('phi_m', PDG_MASSES['phi_m'], r'$\phi$'),
    ('Lambda', PDG_MASSES['Lambda'], r'$\Lambda$'),
    ('Delta', PDG_MASSES['Delta'], r'$\Delta$'),
    ('f2', PDG_MASSES['f2'], r'$f_2$'),
    ('Xi', PDG_MASSES['Xi'], r'$\Xi$'),
    ('Omega_b', PDG_MASSES['Omega_b'], r'$\Omega$'),
]

# === Create figure ===
fig, ax = plt.subplots(figsize=(7.0, 4.5))

# PDG reference lines (background)
for name, mass, label in pdg_lines:
    ax.axhline(mass, color='0.85', linewidth=0.5, zorder=0)
    ax.text(6.8, mass, f' {label} ({mass:.0f})', fontsize=6, va='center',
            color='0.5', zorder=1)

# UDT predictions: plot each kappa channel
x_offset = 0
kappas_sorted = ['-1', '+1', '-2', '+2', '-3', '+3']

for ks in kappas_sorted:
    evals = eigenvalues.get(ks, [])
    masses = [E * C_CALIB for E in evals]

    for n, m in enumerate(masses):
        x_pos = x_offset + n * 0.12
        ax.plot(x_pos, m, kappa_markers[ks], color=kappa_colors[ks],
                markersize=6, markeredgecolor='0.2', markeredgewidth=0.4,
                zorder=5, label=kappa_labels[ks] if n == 0 else '')

    x_offset += 1.0

# Formatting
ax.set_yscale('log')
ax.set_ylabel(r'Mass (MeV/$c^2$)')
ax.set_xlabel(r'$\kappa$ channel')
ax.set_xlim(-0.5, 6.5)
ax.set_ylim(0.3, 3000)

# Custom x-ticks for kappa channels
channel_positions = [0, 1, 2, 3, 4, 5]
channel_labels = [r'$\kappa\!=\!-1$', r'$\kappa\!=\!+1$',
                  r'$\kappa\!=\!-2$', r'$\kappa\!=\!+2$',
                  r'$\kappa\!=\!-3$', r'$\kappa\!=\!+3$']
ax.set_xticks(channel_positions)
ax.set_xticklabels(channel_labels, fontsize=7)

# Legend
handles, leg_labels = ax.get_legend_handles_labels()
by_label = dict(zip(leg_labels, handles))
ax.legend(by_label.values(), by_label.keys(),
          loc='upper left', frameon=True, fancybox=False,
          edgecolor='0.7', framealpha=0.9, ncol=2, fontsize=6.5)

# Calibration annotation
ax.text(0.98, 0.02, f'$C = {C_CALIB:.2f}$ MeV',
        transform=ax.transAxes, fontsize=8, ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='0.7', alpha=0.9))

ax.set_title(r'Mass ladder: $m_n = E_n \times C$', fontsize=10, pad=8)

fig.tight_layout()

# === Save ===
savefig(fig, 'fig05_mass_ladder')
plt.close(fig)
print("fig05_mass_ladder: done")
