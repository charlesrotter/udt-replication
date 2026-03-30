#!/usr/bin/env python3
"""
fig07 -- PMNS mixing angles and CP phase: UDT vs NuFIT 5.2.

Four panels: sin^2(theta_12), sin^2(theta_23), sin^2(theta_13), delta_CP.
NuFIT 1/2/3 sigma bands, UDT prediction as vertical line.
All within 1 sigma. Zero free parameters.

Saves: manuscript/figures/fig07_pmns_contours.{pdf,png}
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from lib.angular_integrals import pmns_mixing_angles
from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator

# === Publication style ===
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['CMU Serif', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'font.size': 9,
    'axes.labelsize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'xtick.minor.width': 0.3,
    'xtick.direction': 'in',
    'xtick.top': True,
    'savefig.dpi': 300,
})

CB_RED = '#D55E00'

# === Predictions ===
pmns = pmns_mixing_angles()
delta_cp_pred = np.pi + np.arcsin(4.0/13.0)
delta_cp_deg = np.degrees(delta_cp_pred)

# (key, xlabel, frac_str, predicted, central, sigma, bands)
panels = [
    ('sin2_12', r'$\sin^2\!\theta_{12}$', r'$\frac{4}{13}$',
     pmns['sin2_12']['predicted'], 0.304, 0.012,
     {'1s': (0.292, 0.316), '2s': (0.280, 0.328), '3s': (0.269, 0.343)}),
    ('sin2_23', r'$\sin^2\!\theta_{23}$', r'$\frac{4}{7}$',
     pmns['sin2_23']['predicted'], 0.573, 0.016,
     {'1s': (0.557, 0.589), '2s': (0.541, 0.605), '3s': (0.415, 0.616)}),
    ('sin2_13', r'$\sin^2\!\theta_{13}$', r'$\frac{1}{45}$',
     pmns['sin2_13']['predicted'], 0.02219, 0.00062,
     {'1s': (0.02157, 0.02281), '2s': (0.02096, 0.02344), '3s': (0.02034, 0.02404)}),
    ('delta_cp', r'$\delta_{\mathrm{CP}}$ (deg)',
     r'$\pi\!+\!\arcsin\!\frac{4}{13}$',
     delta_cp_deg, 197.0, 25.0,
     {'1s': (172.0, 222.0), '2s': (135.0, 250.0), '3s': (107.0, 297.0)}),
]

# === Colors ===
band_colors = {'1s': '#66BB6A', '2s': '#FFD54F', '3s': '#FFB74D'}
band_alphas = {'1s': 0.45, '2s': 0.30, '3s': 0.18}

# === Create figure ===
fig, axes = plt.subplots(1, 4, figsize=(7.0, 2.2))

for ax, (key, xlabel, frac_str, pred, central, sigma, bands) in zip(axes, panels):
    # Bands: 3s first (back), then 2s, 1s on top
    for sig_key in ['3s', '2s', '1s']:
        lo, hi = bands[sig_key]
        ax.axvspan(lo, hi, color=band_colors[sig_key],
                   alpha=band_alphas[sig_key], zorder=1)

    # UDT prediction
    ax.axvline(pred, color=CB_RED, linewidth=1.8, zorder=5)

    # Annotation box
    tension = abs(pred - central) / sigma
    if key == 'delta_cp':
        val_str = f'{pred:.1f}'
    else:
        val_str = f'{pred:.4f}'
    box_text = f'UDT $=$ {frac_str}\n$= {val_str}$\n$({tension:.2f}\\sigma)$'
    ax.text(0.95, 0.93, box_text, transform=ax.transAxes,
            fontsize=5.5, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='0.6', alpha=0.95, linewidth=0.4))

    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_yticks([])

    # x-range: pad beyond 3σ
    lo3, hi3 = bands['3s']
    pad = (hi3 - lo3) * 0.2
    ax.set_xlim(lo3 - pad, hi3 + pad)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_locator(MaxNLocator(5))

# Legend across top
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [
    Patch(facecolor=band_colors['1s'], alpha=band_alphas['1s'],
          edgecolor='none', label=r'NuFIT $1\sigma$'),
    Patch(facecolor=band_colors['2s'], alpha=band_alphas['2s'],
          edgecolor='none', label=r'NuFIT $2\sigma$'),
    Patch(facecolor=band_colors['3s'], alpha=band_alphas['3s'],
          edgecolor='none', label=r'NuFIT $3\sigma$'),
    Line2D([0], [0], color=CB_RED, linewidth=1.8, label='UDT'),
]
fig.legend(handles=legend_elements, loc='upper center', ncol=4,
           bbox_to_anchor=(0.5, 1.03), frameon=True, fancybox=False,
           edgecolor='0.6', framealpha=0.95, fontsize=7,
           handlelength=1.5, columnspacing=1.2)

fig.tight_layout(rect=[0, 0, 1, 0.90])

# === Save ===
savefig(fig, 'fig07_pmns_contours')
plt.close(fig)
print("fig07_pmns_contours: done")
