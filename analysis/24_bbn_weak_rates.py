#!/usr/bin/env python3
"""24 -- BBN weak rates, detailed balance, and electroweak constants.

GENERATES: data/generated/24_bbn_weak_rates.json
VERIFIES: Manuscript BBN section

All constants from UDT metric. No ΛCDM inputs.

tau_n = 953.7 s (+8.6% -- accumulated parameter errors in G_F, g_A, Q)
Detailed balance: n/p = exp(-Q/T) exact to 4 decimal places
G_F = 1/(sqrt(2)*v^2) with v = 504*pi^6*m_e (-1.11%)
M_W = v*sqrt(pi/30) (-0.32%)
M_Z = M_W*sqrt(13/10) (+0.18%)

Completion class: S (structural, parameter-independent).
"""
import os, sys
import numpy as np
from scipy import integrate

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.constants import (MULT_2J1, MULT_2L1, MULT_2KM1, MULT_2KP1,
                           M_E, C_CALIB, ALPHA_EM_PDG)
from lib.utils import save_results, pct_error, format_comparison

pi = np.pi
me = M_E  # 0.511 MeV

# === UDT-derived electroweak constants ===
v_udt = 504.0 * pi**6 * me                  # Higgs VEV in MeV
GF_udt = 1.0 / (np.sqrt(2.0) * v_udt**2)   # Fermi constant MeV^-2
g_A = MULT_2J1**2 / pi                       # 4/pi = 1.2732
axial = 1.0 + 3.0 * g_A**2                   # 5.863
Q_udt = (5.0/4.0) * ALPHA_EM_PDG * C_CALIB  # n-p mass splitting
hbar = 6.582119e-22                           # MeV*s
G0 = GF_udt**2 / (2.0 * pi**3) / hbar       # rate prefactor s^-1

# Electroweak boson masses
MW = v_udt * np.sqrt(pi / 30.0)             # W mass
MZ = MW * np.sqrt(13.0 / 10.0)              # Z mass
sin2_thW = 3.0 / 13.0

print("=" * 60)
print("BBN WEAK RATES FROM THE METRIC")
print("=" * 60)

# === Print electroweak constants ===
print("\nElectroweak constants:")
format_comparison("  G_F", GF_udt, 1.1664e-11, "MeV^-2")
format_comparison("  g_A = 4/pi", g_A, 1.2762)
format_comparison("  Q = (5/4)*alpha*C", Q_udt, 1.29333, "MeV")
format_comparison("  M_W = v*sqrt(pi/30)", MW/1000, 80.379, "GeV")
format_comparison("  M_Z = M_W*sqrt(13/10)", MZ/1000, 91.188, "GeV")
format_comparison("  sin2(thW) = 3/13", sin2_thW, 0.2312)
print(f"  axial factor (1+3*g_A^2) = {axial:.4f}")
print(f"  G0 = {G0:.4e} s^-1")

# === Scalar-safe Fermi-Dirac ===
def fd(x):
    x = float(x)
    if x > 500: return 0.0
    if x < -500: return 1.0
    return 1.0 / (1.0 + np.exp(x))

def focc(E, T): return fd(E / T)
def fvac(E, T): return fd(-E / T)

# === Five weak-rate integrands ===
def I1_decay(Ee, T, Q):
    if Ee <= me or Ee >= Q: return 0.0
    pe = np.sqrt(Ee**2 - me**2)
    Enu = Q - Ee
    return Ee * pe * Enu**2 * fvac(Ee, T) * fvac(Enu, T)

def I2_poscap(Ee, T, Q):
    if Ee <= me: return 0.0
    pe = np.sqrt(Ee**2 - me**2)
    Enu = Q + Ee
    return Ee * pe * Enu**2 * focc(Ee, T) * fvac(Enu, T)

def I3_nucap(Enu, T, Q):
    if Enu <= 0.0: return 0.0
    Ee = Q + Enu
    if Ee <= me: return 0.0
    pe = np.sqrt(Ee**2 - me**2)
    return Ee * pe * Enu**2 * focc(Enu, T) * fvac(Ee, T)

def I4_eleccap(Ee, T, Q):
    if Ee <= Q: return 0.0
    pe = np.sqrt(Ee**2 - me**2)
    Enu = Ee - Q
    return Ee * pe * Enu**2 * focc(Ee, T) * fvac(Enu, T)

def I5_anucap(Enu, T, Q):
    if Enu <= Q + me: return 0.0
    Ee = Enu - Q
    if Ee <= me: return 0.0
    pe = np.sqrt(Ee**2 - me**2)
    return Ee * pe * Enu**2 * focc(Enu, T) * fvac(Ee, T)

def qint(f, a, b):
    try:
        v, _ = integrate.quad(f, a, b, limit=500, epsrel=1e-9, epsabs=1e-30)
        return max(v, 0.0)
    except Exception:
        return 0.0

def Gnp(T, Q):
    Ehi = max(30 * T, 10.0)
    i1 = qint(lambda E: I1_decay(E, T, Q), me + 1e-8, Q - 1e-8)
    i2 = qint(lambda E: I2_poscap(E, T, Q), me + 1e-8, Ehi)
    i3 = qint(lambda E: I3_nucap(E, T, Q), 1e-8, Ehi)
    return axial * G0 * (i1 + i2 + i3)

def Gpn(T, Q):
    Ehi = max(30 * T + Q + me, 10.0)
    i4 = qint(lambda E: I4_eleccap(E, T, Q), Q + 1e-8, Ehi)
    i5 = qint(lambda E: I5_anucap(E, T, Q), Q + me + 1e-8, Ehi)
    return axial * G0 * (i4 + i5)

# === Vacuum neutron lifetime ===
print("\n" + "=" * 60)
print("Vacuum neutron lifetime")
T0 = 0.001
gnp0 = Gnp(T0, Q_udt)
tau_n = 1.0 / gnp0 if gnp0 > 0 else np.inf
tau_n_pdg = 878.4
print(f"  tau_n = {tau_n:.1f} s")
format_comparison("  tau_n", tau_n, tau_n_pdg, "s")

# === Detailed balance ===
print("\n" + "=" * 60)
print("Detailed balance: n/p vs exp(-Q/T)")
print(f"{'T':>6} {'Gnp':>12} {'Gpn':>12} {'n/p_rate':>10} {'exp(-Q/T)':>10} {'ratio':>7}")
print("-" * 62)

balance_data = []
for T in [0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
    gnp = Gnp(T, Q_udt)
    gpn = Gpn(T, Q_udt)
    np_r = gpn / gnp if gnp > 0 else 0.0
    np_B = np.exp(-Q_udt / T)
    ratio = np_r / np_B if np_B > 0 else 0.0
    print(f"{T:>6.1f} {gnp:>12.4e} {gpn:>12.4e} {np_r:>10.5f} {np_B:>10.5f} {ratio:>7.4f}")
    balance_data.append({
        'T_MeV': T, 'Gnp': float(gnp), 'Gpn': float(gpn),
        'np_rate': float(np_r), 'np_boltzmann': float(np_B),
        'ratio': float(ratio),
    })

# === He-4 mass fraction vs T ===
print("\n" + "=" * 60)
print("He-4 mass fraction Y = 2*(n/p)/(1+n/p)")
Y_data = []
for T in [0.3, 0.5, 0.65, 0.7, 1.0, 1.5, 2.0, 3.0]:
    gnp = Gnp(T, Q_udt)
    gpn = Gpn(T, Q_udt)
    np_r = gpn / gnp if gnp > 0 else 0.0
    Y = 2 * np_r / (1 + np_r) if np_r > 0 else 0.0
    print(f"  T={T:.2f} MeV:  n/p={np_r:.4f}  Y={Y:.4f}")
    Y_data.append({'T_MeV': T, 'np': float(np_r), 'Y': float(Y)})
print(f"  Observed Y(He4) = 0.245 +/- 0.003")

# === Gate: detailed balance ===
worst_ratio = max(abs(d['ratio'] - 1.0) for d in balance_data if d['T_MeV'] >= 1.0)
gate_pass = worst_ratio < 1e-3
print(f"\n  Gate: detailed balance worst deviation (T>=1): {worst_ratio:.2e}")
print(f"  Gate: {'PASS' if gate_pass else 'FAIL'}")

results = {
    'electroweak': {
        'G_F_MeV2': float(GF_udt),
        'G_F_pdg': 1.1664e-11,
        'G_F_pct': float(pct_error(GF_udt, 1.1664e-11)),
        'g_A': float(g_A),
        'g_A_pdg': 1.2762,
        'Q_MeV': float(Q_udt),
        'Q_pdg': 1.29333,
        'MW_GeV': float(MW / 1000),
        'MZ_GeV': float(MZ / 1000),
        'sin2_thW': float(sin2_thW),
        'v_GeV': float(v_udt / 1000),
    },
    'tau_n_s': float(tau_n),
    'tau_n_pdg': tau_n_pdg,
    'tau_n_pct': float(pct_error(tau_n, tau_n_pdg)),
    'detailed_balance': balance_data,
    'Y_vs_T': Y_data,
    'gate_detailed_balance': bool(gate_pass),
}
save_results('24_bbn_weak_rates.json', results)
