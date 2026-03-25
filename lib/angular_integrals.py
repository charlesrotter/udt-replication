"""Angular integrals on S^2 for spinor harmonics.

The Dirac equation on the UDT metric separates into radial (Form-T) and
angular (spinor harmonics Omega_kappa) parts. The angular coupling matrix
elements <Omega_kappa|sigma.r_hat|Omega_kappa'> determine selection rules
and coupling strengths.

The Diophantine identity:
  (2j+1)^2 (2l+1)(2|kappa_max|+1) = C(2l+2|kappa_max|+1, 2l+1)
selects j=1/2, l=1, |kappa_max|=3 uniquely.

# SOURCE: Derived from CG sections 7, 13, 19
"""
import numpy as np
from .constants import (J_HALF, ELL, KAPPA_MAX,
                         MULT_2J1, MULT_2L1, MULT_2KM1, MULT_2KP1)


def diophantine_check():
    """Verify the Diophantine selection rule.

    (2j+1)^2 (2l+1)(2|kappa_max|+1) = C(2l+2|kappa_max|+1, 2l+1)

    Returns dict with LHS, RHS, and match status.
    """
    from math import comb
    j2p1 = MULT_2J1
    l2p1 = MULT_2L1
    km2m1 = MULT_2KM1
    km2p1 = MULT_2KP1

    lhs = j2p1**2 * l2p1 * km2p1
    # C(2l + 2|kmax| + 1, 2l + 1) — note: NOT (2l+1)+(2|kmax|+1)
    rhs = comb(2 * ELL + 2 * KAPPA_MAX + 1, 2 * ELL + 1)

    return {
        'LHS': lhs,
        'RHS': rhs,
        'match': lhs == rhs,
        'quantum_numbers': {'j': J_HALF, 'l': ELL, 'kappa_max': KAPPA_MAX},
        'multiplicities': {'2j+1': j2p1, '2l+1': l2p1,
                           '2|kmax|-1': km2m1, '2|kmax|+1': km2p1},
    }


def angular_mass_formulas():
    """Compute the four angular-sector mass formulas.

    Returns dict with predicted masses and PDG comparison.
    """
    from .constants import M_E, PDG_MASSES

    m_pi = 84 * np.pi * M_E          # C(9,3)*pi*m_e
    m_mu = 20 * np.pi**3 / 3 * M_E   # (2j+1)^2*(2|kmax|-1)/(2l+1) * pi^{2l+1} * m_e
    m_p = 6 * np.pi**5 * M_E          # (2j+1)(2l+1) * pi^{2|kmax|-1} * m_e

    return {
        'electron': {'predicted': M_E, 'pdg': PDG_MASSES['e'],
                     'error_pct': 0.0, 'formula': 'anchor'},
        'pion': {'predicted': m_pi, 'pdg': PDG_MASSES['pi'],
                 'error_pct': (m_pi - PDG_MASSES['pi']) / PDG_MASSES['pi'] * 100,
                 'formula': 'C(9,3)*pi*m_e = 84*pi*m_e'},
        'muon': {'predicted': m_mu, 'pdg': PDG_MASSES['mu'],
                 'error_pct': (m_mu - PDG_MASSES['mu']) / PDG_MASSES['mu'] * 100,
                 'formula': '20*pi^3/3 * m_e'},
        'proton': {'predicted': m_p, 'pdg': PDG_MASSES['p'],
                   'error_pct': (m_p - PDG_MASSES['p']) / PDG_MASSES['p'] * 100,
                   'formula': '6*pi^5 * m_e'},
    }


def pmns_mixing_angles():
    """Compute PMNS mixing angles from quantum number fractions.

    sin^2(theta_12) = (2j+1)^2 / ((2j+1)^2 + (2l+1)^2) = 4/13
    sin^2(theta_23) = (2j+1)^2 / (2|kappa_max|+1) = 4/7
    sin^2(theta_13) = 1 / ((2l+1)^2 * (2|kappa_max|-1)) = 1/45

    Returns dict with predictions and NuFIT comparison.
    """
    from .constants import NUFIT

    s12 = MULT_2J1**2 / (MULT_2J1**2 + MULT_2L1**2)  # 4/13
    s23 = MULT_2J1**2 / MULT_2KP1                      # 4/7
    s13 = 1.0 / (MULT_2L1**2 * MULT_2KM1)              # 1/45

    return {
        'sin2_12': {
            'predicted': s12, 'exact': '4/13',
            'nufit': NUFIT['sin2_12'], 'sigma': 0.012,
            'tension_sigma': abs(s12 - NUFIT['sin2_12']) / 0.012,
        },
        'sin2_23': {
            'predicted': s23, 'exact': '4/7',
            'nufit': NUFIT['sin2_23'], 'sigma': 0.016,
            'tension_sigma': abs(s23 - NUFIT['sin2_23']) / 0.016,
        },
        'sin2_13': {
            'predicted': s13, 'exact': '1/45',
            'nufit': NUFIT['sin2_13'], 'sigma': 0.00062,
            'tension_sigma': abs(s13 - NUFIT['sin2_13']) / 0.00062,
        },
    }


def neutrino_mass(alpha_em, m_e=None):
    """Compute neutrino mass from alpha^3 * m_e / (2j+1)^2.

    Parameters
    ----------
    alpha_em : float -- the fine structure constant (derived from I_2)
    m_e : float, optional -- electron mass in MeV (default from constants)

    Returns dict with prediction and comparison.
    """
    from .constants import M_E, NUFIT
    if m_e is None:
        m_e = M_E

    m_nu = alpha_em**3 * m_e / MULT_2J1**2  # alpha^3 * m_e / 4
    m_nu_eV = m_nu * 1e6  # convert MeV to eV

    sqrt_dm2_atm = np.sqrt(NUFIT['dm2_31'])  # eV

    return {
        'predicted_eV': m_nu_eV,
        'sqrt_dm2_atm_eV': sqrt_dm2_atm,
        'error_pct': (m_nu_eV - sqrt_dm2_atm) / sqrt_dm2_atm * 100,
        'formula': 'alpha^3 * m_e / 4',
        'sum_m_nu_eV': m_nu_eV,  # m1=0, m2~sqrt(dm2_21), m3~m_nu
    }
