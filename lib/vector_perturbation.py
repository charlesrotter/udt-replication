"""EM vector mode solver on UDT vacuum background.

The vector perturbation A_mu satisfies Maxwell's equations on the UDT metric.
For a radial mode with angular momentum l=1:
  A''(r) + (2/r)*A'(r) + [omega^2*e^{2phi} - l(l+1)/r^2]*A = 0

The identity g^{tt}g^{rr} = -1/c^2 (phi-independent) gives exact Coulomb
from the metric algebraically.

# SOURCE: Derived from CG section 14
"""
import numpy as np
from scipy.optimize import brentq


def shoot_vector(omega, r, phi, phip, l_ang=1):
    """Shoot vector perturbation from origin to r_star.

    Returns (residual, A) where residual = A'(r_star).
    """
    n = len(r)
    dr = r[1] - r[0]
    e2phi = np.exp(np.clip(2 * phi, -400, 400))

    A = np.zeros(n)
    A[0] = r[0]**(l_ang)  # regular at origin
    Ap = l_ang * r[0]**(l_ang - 1)

    for i in range(n - 1):
        ri = max(r[i], 1e-14)
        e2i = e2phi[i]
        ll1 = l_ang * (l_ang + 1)

        def rhs(A_, Ap_, ri_, e2i_):
            return -(2.0 / ri_) * Ap_ + (ll1 / ri_**2 - omega**2 * e2i_) * A_

        k1 = Ap
        l1 = rhs(A[i], Ap, ri, e2i)

        if i < n - 2:
            rm = 0.5 * (r[i] + r[i + 1])
            e2m = 0.5 * (e2phi[i] + e2phi[i + 1])
        else:
            rm, e2m = ri, e2i

        k2 = Ap + 0.5 * dr * l1
        l2 = rhs(A[i] + 0.5 * dr * k1, k2, rm, e2m)

        k3 = Ap + 0.5 * dr * l2
        l3 = rhs(A[i] + 0.5 * dr * k2, k3, rm, e2m)

        k4 = Ap + dr * l3
        l4 = rhs(A[i] + dr * k3, k4, r[i + 1], e2phi[i + 1])

        A[i + 1] = A[i] + dr * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        Ap = Ap + dr * (l1 + 2 * l2 + 2 * l3 + l4) / 6

    return Ap, A


def find_vector_eigenvalues(r, phi, phip, omega_min=0.1, omega_max=20.0,
                            n_scan=5000, n_modes=5, l_ang=1):
    """Find vector boson eigenvalues by shooting + brentq."""
    omega_scan = np.linspace(omega_min, omega_max, n_scan)
    residuals = np.array([shoot_vector(w, r, phi, phip, l_ang)[0]
                          for w in omega_scan])

    evals = []
    for i in range(len(residuals) - 1):
        if (np.isfinite(residuals[i]) and np.isfinite(residuals[i + 1])
                and residuals[i] * residuals[i + 1] < 0):
            try:
                w = brentq(lambda w: shoot_vector(w, r, phi, phip, l_ang)[0],
                           omega_scan[i], omega_scan[i + 1], xtol=1e-10)
                evals.append(w)
                if len(evals) >= n_modes:
                    break
            except Exception:
                pass

    return sorted(evals)
