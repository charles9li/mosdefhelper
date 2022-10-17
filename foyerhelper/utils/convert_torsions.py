__all__ = []
import numpy as np


def opls_to_rb(opls_coefficients):
    """Converts OPLS proper torsion coefficients to Ryckaert-Bellemans (RB)
    torsion coefficients.

    The OPLS torsion has the following form:
        u(phi) = c0
                + c1 * [1 + cos(phi)]
                + c2 * [1 - cos(2*phi)]
                + c3 * [1 + cos(3*phi)]
    The RB torsion has the following form:
        u(psi) = sum_n Vn * cos^n(psi)
    where
        phi = psi - pi
    The following are conversion between cn and Vn:
        V0 = c0 + c1 + 2*c2 + c3
        V1 = -c1 + 3*c3
        V2 = -2*c2
        V3 = -4*c3
    """
    V0 = opls_coefficients[0] + opls_coefficients[1] + 2*opls_coefficients[2] + opls_coefficients[3]
    V1 = -opls_coefficients[1] + 3*opls_coefficients[3]
    V2 = -2*opls_coefficients[2]
    V3 = -4*opls_coefficients[3]
    return np.array([V0, V1, V2, V3])


if __name__ == '__main__':
    from scipy.constants import R
    coeffs = R / 1000 * opls_to_rb([0, 5.24, 13.0, 0.])
    for c in coeffs:
        print(c)
