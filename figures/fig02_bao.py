#!/usr/bin/env python3
# Generates: Figure 2 in manuscript section 5
# BAO: D_V/r_d vs redshift with BOSS + DESI data
"""
fig02 -- BAO distance-redshift diagram.

Loads: lib/constants.py for cosmological polynomial coefficients
Computes D_V/r_d from UDT geometric polynomial.
Overlays BOSS and DESI data points with error bars.

Note: UDT does not use r_d = 147 Mpc (LCDM acoustic horizon).
r_d is a geometric coherence ruler from the metric.

Saves: manuscript/figures/fig02_bao.{pdf,png}
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import COSMO_K, COSMO_BETA, COSMO_GAMMA, MU_G
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

# === UDT cosmological model ===
def phi_cosmo(r):
    return COSMO_K * r + COSMO_BETA * r**2 + COSMO_GAMMA * r**3

def phi_prime_cosmo(r):
    return COSMO_K + 2 * COSMO_BETA * r + 3 * COSMO_GAMMA * r**2

def z_of_r(r):
    return np.exp(phi_cosmo(r)) - 1.0

def r_of_z(z):
    """Invert phi(r) = ln(1+z) by root finding."""
    from scipy.optimize import brentq
    target = np.log(1.0 + z)
    try:
        return brentq(lambda r: phi_cosmo(r) - target, 1e-6, 60.0)
    except ValueError:
        return np.nan

def H_udt(z):
    """UDT Hubble parameter: H(z) = c * dphi/dr / (1+z) evaluated at r(z).
    In UDT: dz/dr = (1+z)*phi'(r), so H(z) = phi'(r(z)) * c.
    We work in Gpc units so H is in km/s/Mpc via conversion.
    """
    r = r_of_z(z)
    pp = phi_prime_cosmo(r)
    # H = c * phi'(r) in units of Gpc^{-1}
    # Convert: c/Gpc = 3e5 km/s / 1e3 Mpc = 300 km/s/Mpc (roughly)
    return pp  # Keep in Gpc^{-1}

def D_V_over_rd(z, r_d_Gpc):
    """Volume-averaged distance D_V/r_d.

    D_V(z) = [z * D_M(z)^2 / H_Gpc(z)]^{1/3}
    where D_M = r(z), H_Gpc = phi'(r(z))
    """
    r_z = r_of_z(z)
    pp = phi_prime_cosmo(r_z)
    if pp <= 0 or r_z <= 0:
        return np.nan
    D_V = (z * r_z**2 / pp) ** (1.0 / 3.0)
    return D_V / r_d_Gpc

# === UDT geometric r_d ===
# r_d = 1000 / sqrt(2*pi*|k|*|beta|*(1+z_d))
# with z_d ~ 1060 (drag epoch)
z_drag = 1060.0
r_d_Gpc = 1000.0 / np.sqrt(2 * np.pi * abs(COSMO_K) * abs(COSMO_BETA)
                             * (1 + z_drag))

# === Generate UDT curve ===
z_curve = np.linspace(0.1, 2.5, 300)
dv_curve = np.array([D_V_over_rd(z, r_d_Gpc) for z in z_curve])

# === BAO data: BOSS + DESI ===
# BOSS DR12 consensus (Alam et al. 2017)
boss_z = np.array([0.38, 0.51, 0.61])
boss_dv = np.array([10.23, 13.36, 15.45])  # D_V/r_d
boss_err = np.array([0.17, 0.21, 0.25])
boss_labels = ['BOSS DR12'] * 3

# DESI Y1 (2024) -- representative values
desi_z = np.array([0.30, 0.51, 0.71, 0.93, 1.32, 1.49, 2.33])
desi_dv = np.array([7.93, 13.62, 16.85, 21.71, 27.79, 30.69, 39.71])
desi_err = np.array([0.15, 0.25, 0.32, 0.28, 0.69, 0.80, 0.95])

# Ly-alpha BOSS (z~2.34)
lya_z = np.array([2.34])
lya_dv = np.array([37.6])
lya_err = np.array([1.5])

# Compute UDT predictions at data redshifts for residual annotation
udt_at_boss = np.array([D_V_over_rd(z, r_d_Gpc) for z in boss_z])
udt_at_desi = np.array([D_V_over_rd(z, r_d_Gpc) for z in desi_z])

# === Create figure ===
fig, ax = plt.subplots(figsize=(3.4, 3.0))

# UDT curve
ax.plot(z_curve, dv_curve, '-', color=CB_RED, linewidth=1.2,
        label='UDT prediction', zorder=3)

# BOSS data
ax.errorbar(boss_z, boss_dv, yerr=boss_err, fmt='s',
            color=CB_BLUE, markersize=4, elinewidth=0.7,
            capsize=2, capthick=0.5, label='BOSS DR12', zorder=4)

# DESI data
ax.errorbar(desi_z, desi_dv, yerr=desi_err, fmt='D',
            color=CB_GREEN, markersize=3.5, elinewidth=0.7,
            capsize=2, capthick=0.5, label='DESI Y1', zorder=4)

# Ly-alpha
ax.errorbar(lya_z, lya_dv, yerr=lya_err, fmt='^',
            color=CB_ORANGE, markersize=4, elinewidth=0.7,
            capsize=2, capthick=0.5, label=r'Ly$\alpha$ forest', zorder=4)

ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$D_V / r_d$')
ax.set_xlim(0, 2.5)
ax.set_ylim(0, 45)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.legend(loc='upper left', frameon=True, fancybox=False,
          edgecolor='0.7', framealpha=0.9)

# Note on data contamination
ax.text(0.95, 0.08,
        r'$r_d$ geometric (no $\Lambda$CDM prior)',
        transform=ax.transAxes, fontsize=7,
        ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFFDE7',
                  edgecolor='0.7', alpha=0.9))

fig.tight_layout()

# === Save ===
savefig(fig, 'fig02_bao')
plt.close(fig)
print("fig02_bao: done")
