#!/usr/bin/env python3
"""
fig06 -- Neutrino mass scale.

UDT prediction m_nu = alpha^3 * m_e / 4 vs oscillation data.
Zoomed to the physics (0-75 meV). KATRIN bound noted as text.

Saves: manuscript/figures/fig06_neutrino_scale.{pdf,png}
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from lib.constants import M_E, ALPHA_EM_PDG, NUFIT, MULT_2J1
from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# === Publication style ===
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
CB_GREEN  = '#009E73'
CB_RED    = '#D55E00'
CB_ORANGE = '#E69F00'
CB_PURPLE = '#CC79A7'

# === UDT prediction ===
alpha_em = ALPHA_EM_PDG
m_nu_MeV = alpha_em**3 * M_E / MULT_2J1**2
m_nu_eV = m_nu_MeV * 1e6
m_nu_meV = m_nu_eV * 1e3

# === Oscillation values ===
m1 = 0.0
m2 = np.sqrt(NUFIT['dm2_21']) * 1e3  # meV
m3 = np.sqrt(NUFIT['dm2_31']) * 1e3  # meV

err_pct = (m_nu_meV - m3) / m3 * 100

# === Create figure ===
fig, ax = plt.subplots(figsize=(3.4, 2.6))

bar_y = [3, 2, 1, 0]
bar_labels = [r'$m_1$', r'$m_2$', r'$m_3$',
              r'UDT $m_\nu$']
bar_vals = [m1, m2, m3, m_nu_meV]
bar_colors = ['0.82', CB_BLUE, CB_GREEN, CB_RED]

ax.barh(bar_y, bar_vals, height=0.55, color=bar_colors,
        edgecolor='0.25', linewidth=0.5, zorder=3)

# Value labels
for y, v, c in zip(bar_y, bar_vals, bar_colors):
    if v > 0:
        ax.text(v + 1.0, y, f'{v:.1f} meV', va='center', fontsize=7.5,
                fontweight='bold', color='0.15')
    else:
        ax.text(1.0, y, '0 (NO)', va='center', fontsize=7.5, color='0.5')

# Formula box — in the whitespace above m1
ax.text(0.97, 0.88,
        r'$m_\nu = \frac{\alpha^3 m_e}{4}$' + f' = {m_nu_meV:.1f} meV'
        + f'\n' + r'$\sqrt{\Delta m^2_{\rm atm}}$' + f' = {m3:.1f} meV ({err_pct:+.1f}%)',
        transform=ax.transAxes, fontsize=7, ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF3E0',
                  edgecolor='0.6', alpha=0.95, linewidth=0.4))

# Bounds as text annotations — below formula box
ax.text(0.97, 0.65, 'KATRIN: $< 450$ meV\nPlanck: $\\Sigma < 120$ meV',
        transform=ax.transAxes, fontsize=6, ha='right', va='top',
        color='0.5')

ax.set_xlabel('Mass (meV)')
ax.set_yticks(bar_y)
ax.set_yticklabels(bar_labels)
ax.set_xlim(0, 72)
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylim(-0.6, 3.8)

ax.tick_params(axis='y', length=0)
fig.tight_layout()

# === Save ===
savefig(fig, 'fig06_neutrino_scale')
plt.close(fig)
print("fig06_neutrino_scale: done")
