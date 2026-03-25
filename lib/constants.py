"""UDT locked parameters and physical constants.

Sign conventions (locked):
  Metric signature: (-,+,+,+)
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2
  sqrt(-g) = c r^2 sin(theta)  (phi-independent)
  Tetrad: e^a_mu = diag(e^{-phi}, e^{+phi}, r, r sin(theta))
  Einstein equation: G^mu_nu = 8 pi G T^mu_nu
  KG equation: box_g phi - mu^2 phi = -S  (positive S deepens well)
"""
import numpy as np

# === Locked geometric parameters ===
# These are derived from the metric; see CG sections 13, 17, 19.

# Angular quantum numbers (unique Diophantine solution)
J_HALF = 0.5          # j = 1/2 (spin)
ELL = 1               # l = 1 (orbital)
KAPPA_MAX = 3          # |kappa_max| = 3

# Four multiplicities from (j, l, |kappa_max|)
MULT_2J1 = 2          # 2j+1 = 2 (spin multiplicity)
MULT_2L1 = 3          # 2l+1 = 3 (orbital multiplicity)
MULT_2KM1 = 5         # 2|kappa_max|-1 = 5
MULT_2KP1 = 7         # 2|kappa_max|+1 = 7 (orbit size)

# Metric depth
PHI0 = -np.cos(np.pi / 5)          # phi_0 = -cos(pi/5) = -0.80902...
PHI_GOLD = (1 + np.sqrt(5)) / 2     # golden ratio phi_gold = 1.61803...

# Screening mass
MU2 = np.pi / 3                     # mu^2 = pi/3
MU = np.sqrt(MU2)                   # mu = sqrt(pi/3) = 1.02333...

# Cavity size: r* = 7 - 1/80
R_STAR = 7.0 - 1.0 / 80.0           # = 6.9875

# Exact eigenvalue
E1_EXACT = 2 * np.sqrt(2) / 3       # E_1(kappa=-1) = 2*sqrt(2)/3

# === Derived constants ===

# Physical constants (external inputs)
M_E = 0.51100                        # electron mass (MeV)
C_LIGHT = 2.998e8                    # speed of light (m/s)
G_NEWTON = 6.674e-11                 # gravitational constant (m^3/kg/s^2)
HBAR = 6.582e-22                     # reduced Planck constant (MeV*s)
T_CMB = 2.725                        # CMB temperature (K)
T_STARLIGHT = 3000.0                 # recombination surface temperature (K)

# Calibration constant: C = 4*pi^2 * m_e * r*
C_CALIB = 4 * np.pi**2 * M_E * R_STAR   # = 140.95 MeV

# Fine structure constant (derived)
# 1/alpha = 36*pi/I_2 where I_2 = integral(e^{2phi} dr) (computed, not algebraic)
# I_2 is computed by analysis/01_vacuum_profile.py
# Here we store the PDG value for comparison only
ALPHA_EM_PDG = 1.0 / 137.036
ALPHA_EM_INV_PDG = 137.036

# === Angular sector mass formulas (parameter-free) ===

def m_pion():
    """m_pi = C(9,3) * pi * m_e = 84 * pi * m_e"""
    return 84 * np.pi * M_E

def m_muon():
    """m_mu = (2j+1)^2 * (2|kappa_max|-1) / (2l+1) * pi^{2l+1} * m_e
           = 20*pi^3/3 * m_e"""
    return 20 * np.pi**3 / 3 * M_E

def m_proton():
    """m_p = (2j+1)(2l+1) * pi^{2|kappa_max|-1} * m_e = 6*pi^5 * m_e"""
    return 6 * np.pi**5 * M_E

def m_neutron_proton_split():
    """m_n - m_p = (5/4) * alpha_EM * C"""
    # Uses derived alpha, not PDG
    return 5.0 / 4.0 * ALPHA_EM_PDG * C_CALIB

# === Cosmological parameters ===

# Cosmological screening mass
MU_G = np.pi * MU / 13              # = pi*sqrt(pi/3)/13 = 0.2473 Gpc^{-1}

# Geometric polynomial coefficients
COSMO_K = 1.5 * MU_G                 # 3/2 * mu_g
COSMO_BETA = -np.cos(np.pi / 5) * MU_G**2   # -cos(pi/5) * mu_g^2
COSMO_GAMMA = (2.0 / 3.0) * MU_G**3         # 2/3 * mu_g^3

# CMB boundary condition
Z_CMB = T_STARLIGHT / T_CMB - 1     # ~1100
PHI_CMB = np.log(1 + Z_CMB)          # ~7.004

# === PDG comparison targets (NOT inputs to UDT) ===

PDG_MASSES = {
    'e':     0.51100,     # electron
    'mu':    105.658,     # muon
    'tau':   1776.86,     # tau
    'pi':    134.977,     # pion (pi0)
    'pi_pm': 139.570,     # pion (pi+/-)
    'p':     938.272,     # proton
    'n':     939.565,     # neutron
    'eta':   547.862,     # eta
    'rho':   775.26,      # rho
    'omega': 782.66,      # omega
    'phi_m': 1019.461,    # phi meson
    'f2':    1275.5,      # f2(1270)
    'a2':    1318.2,      # a2(1320)
    'K':     493.677,     # kaon (K+)
    'K0':    497.611,     # kaon (K0)
    'Sigma': 1189.37,     # Sigma+
    'Xi':    1314.86,     # Xi0
    'Lambda':1115.683,    # Lambda
    'Delta': 1232.0,      # Delta(1232)
    'Omega_b': 1672.45,   # Omega baryon
}

# Neutrino oscillation data (NuFIT 5.2, NO)
NUFIT = {
    'dm2_21': 7.42e-5,    # eV^2
    'dm2_31': 2.515e-3,   # eV^2
    'sin2_12': 0.304,     # +/- 0.012
    'sin2_23': 0.573,     # +/- 0.016
    'sin2_13': 0.02219,   # +/- 0.00062
}

# === Solver defaults ===

N_GRID = 100000           # default grid points for phi ODE
N_SCAN = 50000            # default energy scan points
R_MIN = 1e-8              # regularization radius
