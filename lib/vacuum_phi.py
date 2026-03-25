"""Vacuum scalar flux integrator for phi(r) background.

Derives from the UDT metric:
  ds^2 = -e^{-2phi}c^2 dt^2 + e^{2phi}dr^2 + r^2 dOmega^2

The covariant Klein-Gordon equation box_g phi = mu^2 phi in flux form:
  J'(r) = r^2 (mu^2 phi - S(r))
  phi'(r) = J e^{2phi} / r^2

where J = r^2 e^{-2phi} phi' is the conserved flux variable.
ICs: phi(0) = phi_0, J(0) = 0 (regularity at origin).

For vacuum (S=0), the nonlinear ODE is integrated with fixed-step RK4.

Sign conventions: see lib/constants.py header.

Completion class: A (finite domain [r_min, r_star]).

# SOURCE: Adapted from UDT/lib/scalar.py (canonical)
# GENERATES: Vacuum phi profile used by all analysis scripts
"""
import numpy as np


def _sexp(x, cap=200.0):
    """Safe exponential with overflow protection."""
    return np.exp(np.clip(x, -cap, cap))


def solve(phi0, r_star, n_steps, mu, r_min=1e-8, source=None):
    """Integrate scalar flux equation on [r_min, r_star].

    Parameters
    ----------
    phi0 : float
        Initial condition phi(0).
    r_star : float
        Domain boundary radius.
    n_steps : int
        Number of grid points.
    mu : float
        Screening mass (mu^2 enters ODE).
    r_min : float
        Starting radius (Frobenius regularization).
    source : array of length n_steps, or None
        Source S(r). None for vacuum (S=0).

    Returns
    -------
    r : ndarray of shape (n_steps,)
    phi : ndarray of shape (n_steps,)
    J : ndarray of shape (n_steps,)
    phip : ndarray of shape (n_steps,) -- phi'(r) = J*e^{2phi}/r^2
    overflow : float or None -- radius where overflow occurred, if any
    """
    r = np.linspace(r_min, r_star, n_steps)
    dr = r[1] - r[0]
    phi = np.zeros(n_steps)
    J = np.zeros(n_steps)
    phi[0] = phi0
    J[0] = 0.0

    mu2 = mu * mu
    overflow = None

    for i in range(n_steps - 1):
        ri, pi, Ji = r[i], phi[i], J[i]
        Si = source[i] if source is not None else 0.0

        def dJ(r_, p_, J_, S_):
            return r_**2 * (mu2 * p_ - S_)

        def dp(r_, p_, J_):
            return J_ * _sexp(2 * p_) / r_**2 if r_ > 1e-14 else 0.0

        # Source at midpoint and next point
        if source is not None:
            Sm = 0.5 * (source[min(i, n_steps - 2)] + source[min(i + 1, n_steps - 1)])
            S1 = source[min(i + 1, n_steps - 1)]
        else:
            Sm = S1 = 0.0

        # RK4 (frozen-coefficient midpoint -- ablation-verified as correct)
        k1J = dJ(ri, pi, Ji, Si) * dr
        k1p = dp(ri, pi, Ji) * dr
        rm = ri + dr / 2
        k2J = dJ(rm, pi + k1p / 2, Ji + k1J / 2, Sm) * dr
        k2p = dp(rm, pi + k1p / 2, Ji + k1J / 2) * dr
        k3J = dJ(rm, pi + k2p / 2, Ji + k2J / 2, Sm) * dr
        k3p = dp(rm, pi + k2p / 2, Ji + k2J / 2) * dr
        r1 = ri + dr
        k4J = dJ(r1, pi + k3p, Ji + k3J, S1) * dr
        k4p = dp(r1, pi + k3p, Ji + k3J) * dr

        J[i + 1] = Ji + (k1J + 2 * k2J + 2 * k3J + k4J) / 6
        phi[i + 1] = pi + (k1p + 2 * k2p + 2 * k3p + k4p) / 6

        if abs(phi[i + 1]) > 200 or abs(J[i + 1]) > 1e80:
            overflow = r[i + 1]
            phi[i + 2:] = phi[i + 1]
            J[i + 2:] = J[i + 1]
            break

    # Derived: phi'(r) = J * e^{2phi} / r^2
    phip = J * _sexp(2 * phi) / r**2
    phip[0] = 0.0  # regularity at origin

    return r, phi, J, phip, overflow


def extract_phi2(r, phi, phi0):
    """Extract phi_2 from phi(r) ~ phi_0 + (1/2)*phi_2*r^2 near origin.

    Uses least-squares fit to the first ~20 grid points.
    """
    n_fit = min(20, len(r) // 10)
    if n_fit < 3:
        return 0.0
    dphi = phi[:n_fit] - phi0
    r2 = r[:n_fit]**2
    phi2 = 2.0 * np.dot(r2, dphi) / np.dot(r2, r2)
    return phi2


def compute_I2(r, phi):
    """Compute I_2 = integral(e^{2*phi(r)} dr) over [r_min, r_star]."""
    return np.trapezoid(_sexp(2 * phi), r)


def compute_box_residual(r, phi, phip, mu):
    """Compute box_g phi - mu^2 phi residual for consistency gate.

    box_g phi = (1/r^2) d/dr(r^2 e^{-2phi} phi')
    For vacuum: should equal mu^2 phi everywhere.
    """
    e2phi = _sexp(2 * phi)
    em2phi = _sexp(-2 * phi)

    # box_g phi via flux derivative
    flux = r**2 * em2phi * phip
    dflux = np.gradient(flux, r)
    box = dflux / r**2
    box[0] = box[1]  # regularize origin

    residual = box - mu**2 * phi
    return residual
