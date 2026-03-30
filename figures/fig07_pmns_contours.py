#!/usr/bin/env python3
# Generates: Figure 7 in manuscript section 15
# PMNS mixing angles + CP phase: UDT predictions vs NuFIT bands
"""
fig07 -- PMNS mixing angles and CP phase comparison.

Four panels: sin^2(theta_12), sin^2(theta_23), sin^2(theta_13), delta_CP.
Each shows NuFIT 1/2/3 sigma bands and UDT prediction as vertical line.

All UDT predictions within 1 sigma of NuFIT 5.2 (NO).
delta_CP = pi + arcsin(4/13) = 197.9 deg (0.04 sigma).

Saves: manuscript/figures/fig07_pmns_contours.{pdf,png}
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.angular_integrals import pmns_mixing_angles
from lib.constants import MULT_2J1, MULT_2L1
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

CB_RED = '#D55E00'

# === Get mixing angle predictions ===
pmns = pmns_mixing_angles()

# delta_CP prediction
delta_cp_pred = np.pi + np.arcsin(4.0/13.0)
delta_cp_deg = np.degrees(delta_cp_pred)

# Panel data: (key, xlabel, fraction_label, predicted, central, sigma, bands)
panels = [
    ('sin2_12', r'$\sin^2\theta_{12}$', r'$\frac{4}{13}$',
     pmns['sin2_12']['predicted'], 0.304, 0.012,
     {'1sig': (0.292, 0.316), '2sig': (0.280, 0.328), '3sig': (0.269, 0.343)}),
    ('sin2_23', r'$\sin^2\theta_{23}$', r'$\frac{4}{7}$',
     pmns['sin2_23']['predicted'], 0.573, 0.016,
     {'1sig': (0.557, 0.589), '2sig': (0.541, 0.605), '3sig': (0.415, 0.616)}),
    ('sin2_13', r'$\sin^2\theta_{13}$', r'$\frac{1}{45}$',
     pmns['sin2_13']['predicted'], 0.02219, 0.00062,
     {'1sig': (0.02157, 0.02281), '2sig': (0.02096, 0.02344), '3sig': (0.02034, 0.02404)}),
    ('delta_cp', r'$\delta_{\rm CP}$ [deg]',
     r'$\pi + \arcsin\!\frac{4}{13}$',
     delta_cp_deg, 197.0, 25.0,
     {'1sig': (172.0, 222.0), '2sig': (135.0, 250.0), '3sig': (107.0, 297.0)}),
]

# === Create figure ===
fig, axes = plt.subplots(1, 4, figsize=(7.5, 2.8))

sigma_colors = ['#4CAF50', '#FFC107', '#FF9800']
sigma_alphas = [0.4, 0.25, 0.15]
sigma_labels = [r'$1\sigma$', r'$2\sigma$', r'$3\sigma$']

for ax, (key, xlabel, frac, pred, central, sigma, bands) in zip(axes, panels):
    # Plot sigma bands (3sig first, then 2, then 1 on top)
    for sigma_i, (sig_label, sig_color, sig_alpha) in enumerate(
            zip(reversed(sigma_labels), reversed(sigma_colors),
                reversed(sigma_alphas))):
        sig_key = ['3sig', '2sig', '1sig'][sigma_i]
        lo, hi = bands[sig_key]
        ax.axvspan(lo, hi, color=sig_color, alpha=sig_alpha,
                   label=f'NuFIT {sig_label}' if ax == axes[0] else '',
                   zorder=1)

    # NuFIT central value
    ax.axvline(central, color='0.4', linewidth=0.8, linestyle=':',
               zorder=2)

    # UDT prediction
    ax.axvline(pred, color=CB_RED, linewidth=2.0, zorder=5,
               label='UDT' if ax == axes[0] else '')

    # Fraction label and tension
    tension = abs(pred - central) / sigma
    ax.text(0.95, 0.92, f'UDT = {frac}\n= {pred:.4f}\n({tension:.2f}$\\sigma$)',
            transform=ax.transAxes, fontsize=6, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='0.7', alpha=0.9))

    ax.set_xlabel(xlabel)
    ax.set_yticks([])

    # Set x-range
    lo3, hi3 = bands['3sig']
    pad = (hi3 - lo3) * 0.25
    ax.set_xlim(lo3 - pad, hi3 + pad)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

# Legend from first axis
handles, labels_leg = axes[0].get_legend_handles_labels()
fig.legend(handles, labels_leg, loc='upper center', ncol=4,
           bbox_to_anchor=(0.5, 1.02), frameon=True, fancybox=False,
           edgecolor='0.7', framealpha=0.9, fontsize=7)

fig.tight_layout(rect=[0, 0, 1, 0.92])

# === Save ===
savefig(fig, 'fig07_pmns_contours')
plt.close(fig)
print("fig07_pmns_contours: done")
