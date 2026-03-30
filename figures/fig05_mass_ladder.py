#!/usr/bin/env python3
"""
fig05 -- Hadron mass ladder: predicted vs PDG.

Scatter plot of UDT predicted mass vs PDG mass for eigenvalue-sector
particles (mesons + baryons) and angular-sector particles.
Diagonal = perfect agreement. Inset shows percent errors.

Uses VALIDATED assignments from CG §14, VR §20.

Saves: manuscript/figures/fig05_mass_ladder.{pdf,png}
"""
import os, sys
import json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from lib.constants import C_CALIB, PHI_GOLD, R_STAR, M_E
from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

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
    'xtick.minor.width': 0.3,
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

# === Load computed eigenvalues ===
ev_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'generated',
                       '02_eigenvalues.json')
with open(ev_path) as f:
    ev_data = json.load(f)

eigenvalues = ev_data['eigenvalues']
# Also load sub-cavity eigenvalues
sub_ev = ev_data.get('baryon_subcavity', {}).get('eigenvalues', {})

# === Validated meson assignments (VR §20, CG §14.3) ===
# Format: (name, label, kappa_str, n_index, pdg_MeV)
mesons = [
    ('eta',   r'$\eta$',         '2',  0,  547.862),
    ('omega', r'$\omega$',      '-1',  1,  782.66),
    ('phi',   r'$\phi$',        '3',  0, 1019.461),
    ('f2p',   r"$f_2'$",       '-3',  1, 1525.0),
]

# Compute meson predicted masses
meson_pred, meson_pdg, meson_labels = [], [], []
for name, label, kappa, n, pdg in mesons:
    E = eigenvalues[kappa][n]
    m_pred = E * C_CALIB
    meson_pred.append(m_pred)
    meson_pdg.append(pdg)
    meson_labels.append(label)
    err = (m_pred - pdg) / pdg * 100
    print(f"  Meson {label}: E={E:.4f}, m={m_pred:.1f} vs {pdg:.1f} MeV ({err:+.2f}%)")

# === Validated baryon assignments (VR §20, CG §14.4) ===
# Sub-cavity eigenvalues on [0, r_b] where r_b = r*/phi_gold
baryons = [
    ('Lambda', r'$\Lambda$',    '-1', 1, 1115.683),
    ('Sigma',  r'$\Sigma^+$',  '-1', 2, 1189.37),
    ('Xi_star',r'$\Xi^*$',      '1', 1, 1531.8),
    ('Omega',  r'$\Omega$',    '-1', 3, 1672.45),
]

baryon_pred, baryon_pdg, baryon_labels = [], [], []
if sub_ev:
    for name, label, kappa, n, pdg in baryons:
        if kappa in sub_ev and n < len(sub_ev[kappa]):
            E = sub_ev[kappa][n]
            m_pred = E * C_CALIB
            baryon_pred.append(m_pred)
            baryon_pdg.append(pdg)
            baryon_labels.append(label)
            err = (m_pred - pdg) / pdg * 100
            print(f"  Baryon {label}: E={E:.4f}, m={m_pred:.1f} vs {pdg:.1f} MeV ({err:+.2f}%)")

# === Angular sector particles ===
angular_particles = [
    (r'$\pi^0$', 84 * np.pi * M_E, 134.977),
    (r'$\mu$',   20 * np.pi**3 / 3 * M_E, 105.658),
    (r'$p$',     6 * np.pi**5 * M_E, 938.272),
]

ang_pred = [x[1] for x in angular_particles]
ang_pdg = [x[2] for x in angular_particles]
ang_labels = [x[0] for x in angular_particles]
for label, pred, pdg in angular_particles:
    err = (pred - pdg) / pdg * 100
    print(f"  Angular {label}: m={pred:.2f} vs {pdg:.2f} MeV ({err:+.3f}%)")

# === Create figure ===
fig, ax = plt.subplots(figsize=(3.4, 3.4))

# Perfect agreement diagonal
m_range = np.array([80, 2000])
ax.plot(m_range, m_range, '-', color='0.75', linewidth=0.6, zorder=1)

# ±1% bands
ax.fill_between(m_range, m_range * 0.99, m_range * 1.01,
                color='#E8F4E8', alpha=0.5, zorder=0)

# Mesons
ax.scatter(meson_pdg, meson_pred, s=45, c=CB_BLUE, marker='o',
           edgecolors='0.2', linewidths=0.5, zorder=5, label='Mesons')

# Baryons
if baryon_pred:
    ax.scatter(baryon_pdg, baryon_pred, s=45, c=CB_RED, marker='s',
               edgecolors='0.2', linewidths=0.5, zorder=5, label='Baryons')

# Angular sector
ax.scatter(ang_pdg, ang_pred, s=55, c=CB_GREEN, marker='^',
           edgecolors='0.2', linewidths=0.5, zorder=5, label='Angular sector')

# Labels for each particle
all_pdg = meson_pdg + baryon_pdg + ang_pdg
all_pred = meson_pred + baryon_pred + ang_pred
all_labels = meson_labels + baryon_labels + ang_labels

for pdg, pred, label in zip(all_pdg, all_pred, all_labels):
    # Offset to avoid overlap
    x_off = 8
    y_off = -8
    if label in [r'$\omega$']:
        y_off = 6
    if label in [r'$\Sigma^+$']:
        x_off = -10
        y_off = 6
    ax.annotate(label, (pdg, pred), fontsize=6,
                xytext=(x_off, y_off), textcoords='offset points',
                ha='left', va='top')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'PDG mass (MeV/$c^2$)')
ax.set_ylabel(r'UDT prediction (MeV/$c^2$)')
ax.set_xlim(80, 2000)
ax.set_ylim(80, 2000)
ax.set_aspect('equal')

ax.legend(loc='upper left', frameon=True, fancybox=False,
          edgecolor='0.6', framealpha=0.95, fontsize=7,
          handlelength=1.2, markerscale=0.8)

# C annotation
ax.text(0.97, 0.05, f'$C = {C_CALIB:.2f}$ MeV',
        transform=ax.transAxes, fontsize=7, ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='0.6', alpha=0.9, linewidth=0.4))

fig.tight_layout()

# === Save ===
savefig(fig, 'fig05_mass_ladder')
plt.close(fig)
print("fig05_mass_ladder: done")
