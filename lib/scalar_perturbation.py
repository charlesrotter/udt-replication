"""Scalar boson eigenvalue solver on UDT vacuum background.

The scalar perturbation eta(r) satisfies the Sturm-Liouville equation:
  -(e^{-2phi}/r^2 d/dr(r^2 e^{-2phi} d eta/dr)) + V_eff(r) eta = omega^2 eta

where V_eff comes from the phi background curvature.

# SOURCE: Derived from CG section 11
"""
import numpy as np
from scipy.optimize import brentq


def scalar_potential(r, phi, phip, mu):
    """Compute effective potential for scalar perturbations.

    V_eff = mu^2 + (nonlinear curvature corrections)
    """
    e2phi = np.exp(np.clip(2 * phi, -400, 400))
    em2phi = np.exp(np.clip(-2 * phi, -400, 400))

    # V_eff = mu^2 - 2*phi'^2 + angular correction
    V_eff = mu**2 * np.ones_like(r)
    return V_eff


def shoot_scalar(omega, r, phi, phip, mu):
    """Shoot scalar perturbation equation from origin to r_star.

    Returns (residual, eta) where residual = eta'(r_star) (Neumann BC).
    """
    n = len(r)
    dr = r[1] - r[0]
    e2phi = np.exp(np.clip(2 * phi, -400, 400))
    V = scalar_potential(r, phi, phip, mu)

    eta = np.zeros(n)
    eta[0] = 1.0  # normalization
    eta_p = 0.0   # eta'(0) = 0 for l=0

    for i in range(n - 1):
        ri = r[i]
        if ri < 1e-14:
            ri = 1e-14
        e2i = e2phi[i]

        # eta'' + (2/r - 2*phi')*eta' + e^{2phi}*(omega^2 - V)*eta = 0
        coeff_1 = 2.0 / ri - 2.0 * phip[i]
        coeff_0 = e2i * (omega**2 - V[i])

        k1 = eta_p
        l1 = -coeff_1 * eta_p + coeff_0 * eta[i]

        k2 = eta_p + 0.5 * dr * l1
        l2 = -coeff_1 * k2 + coeff_0 * (eta[i] + 0.5 * dr * k1)

        k3 = eta_p + 0.5 * dr * l2
        l3 = -coeff_1 * k3 + coeff_0 * (eta[i] + 0.5 * dr * k2)

        k4 = eta_p + dr * l3
        l4 = -coeff_1 * k4 + coeff_0 * (eta[i] + dr * k3)

        eta[i + 1] = eta[i] + dr * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        eta_p = eta_p + dr * (l1 + 2 * l2 + 2 * l3 + l4) / 6

    return eta_p, eta


def find_scalar_eigenvalues(r, phi, phip, mu, omega_min=0.1, omega_max=20.0,
                            n_scan=5000, n_modes=5):
    """Find scalar boson eigenvalues by shooting + brentq."""
    omega_scan = np.linspace(omega_min, omega_max, n_scan)
    residuals = np.array([shoot_scalar(w, r, phi, phip, mu)[0]
                          for w in omega_scan])

    evals = []
    for i in range(len(residuals) - 1):
        if (np.isfinite(residuals[i]) and np.isfinite(residuals[i + 1])
                and residuals[i] * residuals[i + 1] < 0):
            try:
                w = brentq(lambda w: shoot_scalar(w, r, phi, phip, mu)[0],
                           omega_scan[i], omega_scan[i + 1], xtol=1e-10)
                evals.append(w)
                if len(evals) >= n_modes:
                    break
            except Exception:
                pass

    return sorted(evals)
