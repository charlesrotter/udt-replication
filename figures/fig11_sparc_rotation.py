#!/usr/bin/env python3
# Generates: Figure 11 in manuscript section 18
# Representative SPARC rotation curve
"""
fig11 -- SPARC rotation curve: baryons-only vs observed.

Loads real NGC 2403 data from SPARC database (Lelli, McGaugh & Schombert 2016).
Plots Vobs with error bars, Vbar = sqrt(Vgas^2 + Vdisk^2 + Vbul^2),
and shades the gap between Vobs and Vbar (the "missing mass" problem).

This is an honest "open problem" figure: the UDT remnant model is not yet
implemented, so we show only the baryonic shortfall.

Saves: manuscript/figures/fig11_sparc_rotation.{pdf,png}
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

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

# === Load real SPARC data ===
data_path = os.path.join(os.path.dirname(__file__), '..',
                         'data', 'external', 'sparc', 'Rotmod_LTG',
                         'NGC2403_rotmod.dat')
data = np.loadtxt(data_path)
# Columns: Rad(kpc) Vobs(km/s) errV Vgas Vdisk Vbul SBdisk SBbul
r_kpc = data[:, 0]
Vobs  = data[:, 1]
errV  = data[:, 2]
Vgas  = data[:, 3]
Vdisk = data[:, 4]
Vbul  = data[:, 5]

# Baryonic velocity: quadrature sum of gas + disk + bulge
# Note: SPARC velocities can be negative (sign encodes direction of
# contribution for decomposition), so we use signed squares.
Vbar = np.sqrt(np.abs(Vgas)*Vgas + np.abs(Vdisk)*Vdisk + np.abs(Vbul)*Vbul)

# === Create figure ===
fig, ax = plt.subplots(figsize=(3.4, 3.0))

# Shade the gap between Vobs and Vbar (the missing-mass region)
ax.fill_between(r_kpc, Vbar, Vobs, where=(Vobs > Vbar),
                color=CB_PURPLE, alpha=0.15, zorder=1,
                label='Missing mass gap')

# Baryonic curve
ax.plot(r_kpc, Vbar, '-', color=CB_ORANGE, linewidth=1.2,
        label=r'$V_{\rm bar}$ (gas + disk + bulge)', zorder=2)

# Observed data with error bars
ax.errorbar(r_kpc, Vobs, yerr=errV, fmt='o',
            color=CB_BLUE, markersize=2.5, elinewidth=0.5,
            capsize=1.2, capthick=0.4, label=r'$V_{\rm obs}$ (SPARC)',
            zorder=4)

ax.set_xlabel('Radius (kpc)')
ax.set_ylabel(r'$V_{\rm rot}$ (km s$^{-1}$)')
ax.set_xlim(0, r_kpc.max() * 1.05)
ax.set_ylim(0, max(Vobs.max(), Vbar.max()) * 1.2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax.legend(loc='lower right', frameon=True, fancybox=False,
          edgecolor='0.7', framealpha=0.9, fontsize=7)

# Galaxy name
ax.text(0.05, 0.95, 'NGC 2403',
        transform=ax.transAxes, fontsize=9, va='top',
        fontweight='bold', color='0.3')

# Open-problem annotation
ax.text(0.05, 0.84,
        'Baryons-only: gap = open problem\n'
        r'UDT: $g_{\rm tt}$ remnants (not yet computed)',
        transform=ax.transAxes, fontsize=6.5, va='top',
        color='0.45', linespacing=1.4)

fig.tight_layout()

# === Save ===
savefig(fig, 'fig11_sparc_rotation')
plt.close(fig)
print("fig11_sparc_rotation: done")
