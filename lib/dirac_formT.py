"""Form-T Dirac solver on UDT vacuum background.

Derives from the UDT metric via tetrad and spin connection:
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2

Canonical Form-T (massless m=0):
  G' = (-kappa/r + phi')G + E*e^{2phi}*F
  F' = (+kappa/r + phi')F - E*e^{2phi}*G

Boundary conditions:
  Origin: Frobenius regularity (kappa-dependent power-law start)
  r_star: Geometric Neumann G'(r*) = 0

Shooting residual (no finite differences):
  R(E) = (-kappa/r* + phi'*)G* + E*e^{2phi*}*F*
  Root R(E)=0 <=> G'(r*)=0.

Normalization: integral(e^phi (G^2 + F^2) r^2 dr) = 1  (SL weight)

Completion class: A (finite domain [r_min, r_star]).

# SOURCE: Adapted from UDT/lib/dirac_formT.py (canonical)
# GENERATES: All Dirac eigenvalues and wavefunctions for mass predictions
"""
import numpy as np
from scipy.optimize import brentq

try:
    import torch
    _HAS_GPU = torch.cuda.is_available()
    if _HAS_GPU:
        _device = torch.device('cuda')
except ImportError:
    _HAS_GPU = False


def _sexp(x, cap=200.0):
    """Safe exponential with overflow protection."""
    return np.exp(np.clip(x, -cap, cap))


# ----------------------------------------------------------------
# Frobenius starts
# ----------------------------------------------------------------

def frobenius_start(kappa, E, phi0, phi2, r_min):
    """Compute Frobenius-regular (G, F) at r = r_min.

    For kappa=-1: G ~ r, F ~ -(E0/3)r^2
    For kappa=+1: F ~ r, G ~ +(E0/3)r^2
    General |kappa|: power-law fallback.

    Parameters
    ----------
    kappa : int -- angular quantum number
    E : float -- trial eigenvalue
    phi0 : float -- phi(0)
    phi2 : float -- coefficient in phi(r) ~ phi_0 + (1/2)*phi_2*r^2
    r_min : float -- starting radius

    Returns
    -------
    G0, F0 : float
    """
    E0 = E * np.exp(2 * phi0)
    r, r2, r4 = r_min, r_min**2, r_min**4

    if kappa == -1:
        G0 = r * (1.0 + 0.5 * (phi2 - E0**2 / 3.0) * r2)
        F0 = -(E0 / 3.0) * r2 + (E0 / 30.0) * (E0**2 - 11.0 * phi2) * r4
    elif kappa == +1:
        F0 = r * (1.0 + 0.5 * (phi2 - E0**2 / 3.0) * r2)
        G0 = (E0 / 3.0) * r2 + (E0 / 30.0) * (11.0 * phi2 - E0**2) * r4
    else:
        ak = abs(kappa)
        if kappa > 0:
            G0 = r**ak
            F0 = r**(ak - 1)
        else:
            G0 = r**(ak - 1)
            F0 = r**ak

    return G0, F0


# ----------------------------------------------------------------
# CPU shooting
# ----------------------------------------------------------------

def shoot(E, kappa, r, phip, e2phi, phi0=None, phi2=0.0):
    """Shoot Form-T from r[0] to r[-1]. Returns (residual, G, F).

    The residual is G'(r*) computed from the ODE:
      R(E) = (-kappa/r* + phi'*)G* + E*e^{2phi*}*F*
    """
    n = len(r)
    dr = r[1] - r[0]
    r_min = r[0]

    if phi0 is not None:
        G0, F0 = frobenius_start(kappa, E, phi0, phi2, r_min)
    else:
        ak = abs(kappa)
        G0 = r_min**ak if kappa > 0 else r_min**(ak - 1)
        F0 = r_min**(ak - 1) if kappa > 0 else r_min**ak

    G = np.zeros(n)
    F = np.zeros(n)
    G[0], F[0] = G0, F0

    for i in range(n - 1):
        ri, ppi, e2pi = r[i], phip[i], e2phi[i]
        Gi, Fi = G[i], F[i]

        if i < n - 2:
            ri_m = 0.5 * (r[i] + r[i + 1])
            ppi_m = 0.5 * (phip[i] + phip[i + 1])
            e2pi_m = 0.5 * (e2phi[i] + e2phi[i + 1])
        else:
            ri_m, ppi_m, e2pi_m = ri, ppi, e2pi

        def dG(G_, F_, ri_, ppi_, e2pi_):
            return (-kappa / ri_ + ppi_) * G_ + E * e2pi_ * F_

        def dF(G_, F_, ri_, ppi_, e2pi_):
            return (kappa / ri_ + ppi_) * F_ - E * e2pi_ * G_

        ri1 = r[i + 1]
        ppi1 = phip[i + 1]
        e2pi1 = e2phi[i + 1]

        k1G = dG(Gi, Fi, ri, ppi, e2pi) * dr
        k1F = dF(Gi, Fi, ri, ppi, e2pi) * dr
        k2G = dG(Gi + k1G / 2, Fi + k1F / 2, ri_m, ppi_m, e2pi_m) * dr
        k2F = dF(Gi + k1G / 2, Fi + k1F / 2, ri_m, ppi_m, e2pi_m) * dr
        k3G = dG(Gi + k2G / 2, Fi + k2F / 2, ri_m, ppi_m, e2pi_m) * dr
        k3F = dF(Gi + k2G / 2, Fi + k2F / 2, ri_m, ppi_m, e2pi_m) * dr
        k4G = dG(Gi + k3G, Fi + k3F, ri1, ppi1, e2pi1) * dr
        k4F = dF(Gi + k3G, Fi + k3F, ri1, ppi1, e2pi1) * dr

        G[i + 1] = Gi + (k1G + 2 * k2G + 2 * k3G + k4G) / 6
        F[i + 1] = Fi + (k1F + 2 * k2F + 2 * k3F + k4F) / 6

    bc = (-kappa / r[-1] + phip[-1]) * G[-1] + E * e2phi[-1] * F[-1]
    return bc, G, F


# ----------------------------------------------------------------
# GPU batch scanning
# ----------------------------------------------------------------

def gpu_scan(E_arr, kappa, r, phip, e2phi, phi0=None, phi2=0.0):
    """GPU-vectorized shooting over many E values. Returns BC residuals."""
    if not _HAS_GPU:
        return np.array([shoot(E, kappa, r, phip, e2phi, phi0, phi2)[0]
                         for E in E_arr])

    n_E = len(E_arr)
    n_r = len(r)
    dr = r[1] - r[0]
    r_min = r[0]

    E_t = torch.tensor(E_arr, dtype=torch.float64, device=_device)
    e2phi0 = np.exp(2 * phi0) if phi0 is not None else 1.0

    if phi0 is not None and abs(kappa) == 1:
        E0_t = E_t * e2phi0
        r2 = r_min**2
        r4 = r_min**4
        if kappa == -1:
            G_t = r_min * (1.0 + 0.5 * (phi2 - E0_t**2 / 3.0) * r2)
            F_t = -(E0_t / 3.0) * r2 + (E0_t / 30.0) * (E0_t**2 - 11.0 * phi2) * r4
        else:
            F_t = r_min * (1.0 + 0.5 * (phi2 - E0_t**2 / 3.0) * r2)
            G_t = (E0_t / 3.0) * r2 + (E0_t / 30.0) * (11.0 * phi2 - E0_t**2) * r4
    else:
        ak = abs(kappa)
        g0 = r_min**ak if kappa > 0 else r_min**(ak - 1)
        f0 = r_min**(ak - 1) if kappa > 0 else r_min**ak
        G_t = torch.full((n_E,), g0, dtype=torch.float64, device=_device)
        F_t = torch.full((n_E,), f0, dtype=torch.float64, device=_device)

    r_t = torch.tensor(r, dtype=torch.float64, device=_device)
    pp_t = torch.tensor(phip, dtype=torch.float64, device=_device)
    e2p_t = torch.tensor(e2phi, dtype=torch.float64, device=_device)

    for i in range(n_r - 1):
        ri = r_t[i]
        ppi = pp_t[i]
        e2pi = e2p_t[i]
        if i < n_r - 2:
            ri_m = 0.5 * (r_t[i] + r_t[i + 1])
            ppi_m = 0.5 * (pp_t[i] + pp_t[i + 1])
            e2pi_m = 0.5 * (e2p_t[i] + e2p_t[i + 1])
        else:
            ri_m, ppi_m, e2pi_m = ri, ppi, e2pi

        ri1 = r_t[i + 1]
        ppi1 = pp_t[i + 1]
        e2pi1 = e2p_t[i + 1]

        def step(G_, F_, ri_, ppi_, e2pi_):
            dG = (-kappa / ri_ + ppi_) * G_ + E_t * e2pi_ * F_
            dF = (kappa / ri_ + ppi_) * F_ - E_t * e2pi_ * G_
            return dG, dF

        k1G, k1F = step(G_t, F_t, ri, ppi, e2pi)
        k1G *= dr; k1F *= dr
        k2G, k2F = step(G_t + k1G / 2, F_t + k1F / 2, ri_m, ppi_m, e2pi_m)
        k2G *= dr; k2F *= dr
        k3G, k3F = step(G_t + k2G / 2, F_t + k2F / 2, ri_m, ppi_m, e2pi_m)
        k3G *= dr; k3F *= dr
        k4G, k4F = step(G_t + k3G, F_t + k3F, ri1, ppi1, e2pi1)
        k4G *= dr; k4F *= dr

        G_t = G_t + (k1G + 2 * k2G + 2 * k3G + k4G) / 6
        F_t = F_t + (k1F + 2 * k2F + 2 * k3F + k4F) / 6

    bc = (-kappa / r_t[-1] + pp_t[-1]) * G_t + E_t * e2p_t[-1] * F_t
    return bc.cpu().numpy()


# ----------------------------------------------------------------
# Eigenvalue finder
# ----------------------------------------------------------------

def find_eigenvalues(kappa, r, phip, e2phi, phi0=None, phi2=0.0,
                     E_min=0.1, E_max=120.0, n_scan=50000, n_modes=10):
    """Find Dirac eigenvalues by bracket scanning + brentq refinement.

    Returns list of eigenvalues (sorted, ascending).
    """
    E_scan = np.linspace(E_min, E_max, n_scan)

    # Coarse scan (GPU if available)
    bc = gpu_scan(E_scan, kappa, r, phip, e2phi, phi0, phi2)

    # Find brackets and refine with brentq
    evals = []
    for i in range(len(bc) - 1):
        if np.isfinite(bc[i]) and np.isfinite(bc[i + 1]) and bc[i] * bc[i + 1] < 0:
            try:
                E_r = brentq(
                    lambda E: shoot(E, kappa, r, phip, e2phi, phi0, phi2)[0],
                    E_scan[i], E_scan[i + 1], xtol=1e-10
                )
                evals.append(E_r)
                if len(evals) >= n_modes:
                    break
            except Exception:
                pass

    return sorted(evals)


# ----------------------------------------------------------------
# Wavefunction (normalized)
# ----------------------------------------------------------------

def wavefunction(E, kappa, r, phip, e2phi, ephi, phi0=None, phi2=0.0):
    """Compute physically normalized wavefunction.

    Returns (G, F) with integral(e^phi (G^2 + F^2) r^2 dr) = 1 (SL weight).
    """
    _, G, F = shoot(E, kappa, r, phip, e2phi, phi0, phi2)
    norm2 = np.trapezoid((G**2 + F**2) * ephi * r**2, r)
    if norm2 > 0:
        s = 1.0 / np.sqrt(norm2)
        G = G * s
        F = F * s
    return G, F


# ----------------------------------------------------------------
# Source integral
# ----------------------------------------------------------------

def source_integral(G, F, r, e2phi):
    """Compute source(kappa) = integral((G^2 - F^2) e^{2phi} r^2 dr)."""
    return np.trapezoid((G**2 - F**2) * e2phi * r**2, r)
