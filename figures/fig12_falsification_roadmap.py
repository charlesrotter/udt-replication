#!/usr/bin/env python3
# Generates: Figure 12 in manuscript section 20
# Falsification roadmap: experiments vs UDT predictions timeline
"""
fig12 -- Falsification roadmap.

Timeline of experiments that can test UDT predictions.
X-axis: time (2025-2035)
Y-axis: experiment categories
Color-coded bars showing when each experiment reaches decisive sensitivity.

Saves: manuscript/figures/fig12_falsification_roadmap.{pdf,png}
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.utils import savefig

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import MultipleLocator

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

# Color categories
CAT_COLORS = {
    'neutrino': '#0072B2',
    'cosmology': '#E69F00',
    'particle': '#009E73',
    'gravity': '#D55E00',
    'astro': '#CC79A7',
}

# === Experiment data ===
# (name, category, start_year, decisive_year, UDT_prediction, status)
experiments = [
    # Neutrino experiments
    ('JUNO', 'neutrino', 2025, 2028,
     r'NO hierarchy, $\Delta m^2_{31}$', 'commissioning'),
    ('DUNE', 'neutrino', 2028, 2033,
     r'$\theta_{23}$ octant, $\delta_{\rm CP}$', 'construction'),
    ('Hyper-K', 'neutrino', 2027, 2032,
     r'$\sin^2\theta_{23} = 4/7$', 'construction'),
    ('KATRIN (final)', 'neutrino', 2025, 2026,
     r'$m_\nu < 0.3$ eV', 'running'),
    ('Project 8', 'neutrino', 2028, 2032,
     r'$m_\nu \sim 48$ meV', 'R&D'),

    # Cosmology experiments
    ('DESI (full)', 'cosmology', 2025, 2028,
     r'$D_V/r_d$ at $z < 2.5$', 'running'),
    ('CMB-S4', 'cosmology', 2030, 2034,
     r'TT peaks, $\Sigma m_\nu$', 'planned'),
    ('Euclid', 'cosmology', 2025, 2029,
     r'BAO, weak lensing', 'running'),
    ('Rubin/LSST', 'cosmology', 2025, 2030,
     r'SNe Ia $\mu(z)$', 'commissioning'),

    # Particle physics
    ('LHCb Run 3', 'particle', 2025, 2028,
     r'$m_p/m_\pi$ precision', 'running'),
    ('CODATA 2026', 'particle', 2025, 2026,
     r'$\alpha_{\rm EM}^{-1} = 36\pi/I_2$', 'in progress'),

    # Gravitational
    ('LISA Pathfinder+', 'gravity', 2030, 2035,
     r'$G_{\rm eff}(r)$ profile', 'planned'),

    # Astrophysics
    ('JWST', 'astro', 2025, 2030,
     r'High-$z$ structure', 'running'),
    ('SKA Phase 1', 'astro', 2027, 2032,
     r'Rotation curves', 'construction'),
]

# Sort by start year within category
experiments.sort(key=lambda x: (list(CAT_COLORS.keys()).index(x[1]), x[2]))

# === Create figure ===
fig, ax = plt.subplots(figsize=(7.0, 5.0))

n_exp = len(experiments)
y_positions = list(range(n_exp))

# Status hatching patterns
status_hatch = {
    'running': '',
    'commissioning': '///',
    'construction': '...',
    'R&D': 'xxx',
    'planned': '\\\\\\',
    'in progress': '',
}

for i, (name, cat, start, decisive, prediction, status) in enumerate(experiments):
    y = n_exp - 1 - i
    color = CAT_COLORS[cat]
    hatch = status_hatch.get(status, '')

    # Bar from start to decisive year
    duration = decisive - start
    bar = ax.barh(y, duration, left=start, height=0.6,
                  color=color, alpha=0.7, edgecolor='0.3', linewidth=0.5,
                  hatch=hatch, zorder=3)

    # Decisive year marker
    ax.plot(decisive, y, 'D', color='0.2', markersize=4, zorder=5)

    # Experiment name
    ax.text(start - 0.15, y, name, ha='right', va='center', fontsize=7,
            fontweight='bold')

    # UDT prediction (to the right)
    ax.text(decisive + 0.2, y, prediction, ha='left', va='center',
            fontsize=5.5, color='0.4', style='italic')

# Formatting
ax.set_xlim(2024, 2036)
ax.set_ylim(-0.8, n_exp - 0.2)
ax.set_xlabel('Year')
ax.set_yticks([])
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(MultipleLocator(1))

# Vertical line for "now"
ax.axvline(2026, color='0.6', linewidth=0.8, linestyle=':', zorder=1)
ax.text(2026, n_exp - 0.5, 'now', ha='center', fontsize=7, color='0.5')

# Category legend
legend_elements = []
for cat_name, cat_color in CAT_COLORS.items():
    from matplotlib.patches import Patch
    legend_elements.append(Patch(facecolor=cat_color, alpha=0.7,
                                  edgecolor='0.3', linewidth=0.5,
                                  label=cat_name.capitalize()))

# Status legend
from matplotlib.patches import Patch
status_items = [
    Patch(facecolor='0.7', hatch='', edgecolor='0.3', label='Running'),
    Patch(facecolor='0.7', hatch='///', edgecolor='0.3', label='Commissioning'),
    Patch(facecolor='0.7', hatch='...', edgecolor='0.3', label='Construction'),
    Patch(facecolor='0.7', hatch='\\\\\\', edgecolor='0.3', label='Planned'),
]

leg1 = ax.legend(handles=legend_elements, loc='lower right',
                  title='Category', title_fontsize=7,
                  frameon=True, fancybox=False, edgecolor='0.7',
                  framealpha=0.9, fontsize=6)
ax.add_artist(leg1)

leg2 = ax.legend(handles=status_items, loc='lower left',
                  title='Status', title_fontsize=7,
                  frameon=True, fancybox=False, edgecolor='0.7',
                  framealpha=0.9, fontsize=6)

ax.set_title('Falsification roadmap: decisive experiments for UDT',
             fontsize=10, pad=10)

# Grid
ax.grid(axis='x', color='0.9', linewidth=0.3, zorder=0)

fig.tight_layout()

# === Save ===
savefig(fig, 'fig12_falsification_roadmap')
plt.close(fig)
print("fig12_falsification_roadmap: done")
