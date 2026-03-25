#!/usr/bin/env python3
# Generates: Figure 6 in manuscript section 13
# Neutrino mass scale comparison
"""
fig06 -- Neutrino mass scale.

Loads: lib/constants.py, lib/angular_integrals.py
UDT prediction: m_nu = alpha^3 * m_e / 4
Compare with sqrt(dm2_atm) and experimental bounds.

Shows hierarchy: m1 = 0, m2 = sqrt(dm2_21), m3 = sqrt(dm2_31)
Overlays KATRIN upper bound and Planck cosmological bound.

Saves: manuscript/figures/fig06_neutrino_scale.{pdf,png}
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import M_E, ALPHA_EM_PDG, NUFIT, MULT_2J1
from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
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

# === UDT prediction ===
alpha_em = ALPHA_EM_PDG  # Use PDG alpha for comparison (UDT derives ~same)
m_nu_MeV = alpha_em**3 * M_E / MULT_2J1**2  # alpha^3 * m_e / 4
m_nu_eV = m_nu_MeV * 1e6  # convert to eV

# === Experimental/oscillation values ===
dm2_21 = NUFIT['dm2_21']  # eV^2
dm2_31 = NUFIT['dm2_31']  # eV^2

# Normal ordering: m1 = 0, m2 = sqrt(dm2_21), m3 = sqrt(dm2_31)
m1 = 0.0
m2 = np.sqrt(dm2_21) * 1e3  # meV
m3 = np.sqrt(dm2_31) * 1e3  # meV
m_nu_meV = m_nu_eV * 1e3    # meV

# Experimental bounds (meV)
katrin_bound = 450.0    # 0.45 eV -> 450 meV (Aker et al. 2024, 90% CL)
planck_sum = 120.0      # sum_m_nu < 0.12 eV -> each < ~40 meV (if degenerate)
planck_per_nu = 40.0    # rough per-species from sum bound

# === Create figure ===
fig, ax = plt.subplots(figsize=(3.4, 3.2))

# --- Mass hierarchy bars ---
bar_y = [3, 2, 1, 0]
bar_labels = [r'$m_1$', r'$m_2$', r'$m_3$', r'UDT $m_\nu$']
bar_vals = [m1, m2, m3, m_nu_meV]
bar_colors = ['0.8', CB_BLUE, CB_GREEN, CB_RED]

bars = ax.barh(bar_y, bar_vals, height=0.5, color=bar_colors,
               edgecolor='0.3', linewidth=0.5, zorder=3)

# Value labels
for y, v, c in zip(bar_y, bar_vals, bar_colors):
    if v > 0:
        ax.text(v + 1.5, y, f'{v:.2f} meV', va='center', fontsize=8,
                fontweight='bold', color='0.2')
    else:
        ax.text(1.0, y, '0 (NO)', va='center', fontsize=8,
                color='0.5')

# --- Experimental bounds ---
# KATRIN
ax.axvline(katrin_bound, color=CB_PURPLE, linewidth=1.0, linestyle='--',
           zorder=2)
ax.text(katrin_bound + 5, 3.3, 'KATRIN\n(0.45 eV)', fontsize=6.5,
        color=CB_PURPLE, va='bottom')

# Planck cosmological bound (per species)
ax.axvline(planck_per_nu, color=CB_ORANGE, linewidth=1.0, linestyle=':',
           zorder=2)
ax.text(planck_per_nu + 3, 3.3, r'Planck $\Sigma m_\nu$', fontsize=6.5,
        color=CB_ORANGE, va='bottom')

# --- Formula annotation ---
ax.text(0.95, 0.05,
        r'UDT: $m_\nu = \frac{\alpha^3 m_e}{4} = '
        f'{m_nu_meV:.2f}$ meV',
        transform=ax.transAxes, fontsize=8, ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF3E0',
                  edgecolor='0.7', alpha=0.9))

# Compare with sqrt(dm2_atm)
sqrt_atm = np.sqrt(dm2_31) * 1e3  # meV
err_pct = (m_nu_meV - sqrt_atm) / sqrt_atm * 100
ax.text(0.95, 0.18,
        r'$\sqrt{\Delta m^2_{\rm atm}} = $' + f'{sqrt_atm:.2f} meV ({err_pct:+.1f}%)',
        transform=ax.transAxes, fontsize=7, ha='right', va='bottom',
        color='0.4')

ax.set_xlabel('Mass (meV)')
ax.set_yticks(bar_y)
ax.set_yticklabels(bar_labels)
ax.set_xlim(0, 520)
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.set_title('Neutrino mass scale', fontsize=10, pad=8)

fig.tight_layout()

# === Save ===
savefig(fig, 'fig06_neutrino_scale')
plt.close(fig)
print("fig06_neutrino_scale: done")
