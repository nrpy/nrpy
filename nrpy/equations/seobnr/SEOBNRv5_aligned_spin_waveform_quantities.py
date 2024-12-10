"""
Construct symbolic expression for the SEOBNRv5 aligned-spin gravitational-wave strain and flux.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
from typing import Any, Dict, List

import sympy as sp

from nrpy.equations.grhd.Min_Max_and_Piecewise_Expressions import coord_greater_bound

# The name of this module ("WaveEquation") is given by __name__:
thismodule = __name__


def complex_mult(z1: List[Any], z2: List[Any]) -> List[Any]:
    """
    Multiply two complex numbers given as list of real and imaginary parts.

    This functions takes two lists containing the real and imaginary part of a complex number
    and returns a list with the real and imaginary part of the resulting multiple.

    :param z1: Complex number 1 as list [Real(z1),Imag(z1)]
    :param z2: Complex number 2 as list [Real(z2),Imag(z2)]
    :return: Complex number z1 x z2 as list [Real(z1*z2),Imag(z1*z2)]

    >>> z1 = [1,2]
    >>> z2 = [3,5]
    >>> complex_mult(z1,z2)
    [-7, 11]

    >>> import sympy as sp
    >>> x1 , y1 , x2 , y2 = sp.symbols('x1 y1 x2 y2',real = True)
    >>> z1 = [x1,y1]
    >>> z2 = [x2,y2]
    >>> complex_mult(z1,z2)
    [x1*x2 - y1*y2, x1*y2 + x2*y1]
    """
    # complex multiplication
    # z1 = x1 + I*y1
    # z2 = x2 + I*y2
    # z1*z2 = x1*x2 - y1*y2 + I*(x1*y2 + x2*y1)
    RE = 0
    IM = 1
    return [z1[RE] * z2[RE] - z1[IM] * z2[IM], z1[RE] * z2[IM] + z1[IM] * z2[RE]]


def f2r(input_float: float) -> sp.Rational:
    """
    Convert a floating-point number to a high-precision rational number.

    This function takes a floating-point number, converts it to a string,
    and appends 60 zeros to increase the precision of the conversion to a rational number.

    :param input_float: The floating-point number to convert.
    :return: A sympy Rational number with high precision.

    >>> f2r(0.1)
    1/10
    >>> f2r(1.5)
    3/2
    >>> f2r(2.0)
    2
    """
    # Convert the input float to a string
    float_as_string = str(input_float)

    # Ensure the string has a decimal point
    if "." not in float_as_string:
        float_as_string = f"{float_as_string}."

    # Append 60 zeros after the decimal of the floating point number to increase precision
    return sp.Rational(float_as_string + "0" * 60)


class SEOBNRv5_aligned_spin_waveform_quantities:
    """Class for computing the SEOBNRv5 aligned-spin gravitational-wave strain and flux."""

    def __init__(self) -> None:
        """
        Compute the SEOBNRv5 aligned-spin Hamiltonian.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the aligned-spin SEOBNRv5 flux and strain. It
        initializes class variables like mass parameters, spin parameters, and
        various coefficients required for the waveforms's multipole modes.

        Inputs: 'm1', 'm2', 'r', 'phi', 'prstar', 'pphi', 'chi1', 'chi2', 'Hreal', 'Omega' and 'Omega_circ'
        Outputs: 'flux' and 'hlms'
        """
        (
            m1,
            m2,
            self.r,
            self.phi,
            self.prstar,
            self.pphi,
            chi1,
            chi2,
            self.Hreal,
            self.Omega,
            self.Omega_circ,
        ) = sp.symbols(
            "m1 m2 r phi prstar pphi chi1 chi2 Hreal Omega Omega_circ", real=True
        )
        self.rho = {}
        for l in range(2, 9):
            for m in range(1, l + 1):
                self.rho.update({f"({l} , {m})": 0})
        self.deltalm = {}
        self.fspin = {}
        self.fspin_limit = {}
        for l in range(2, 6):
            for m in range(l, min(l - 2, 4), -1):
                self.deltalm.update({f"({l} , {m})": 0})
                self.fspin.update({f"({l} , {m})": 0})
                self.fspin_limit.update({f"({l} , {m})": 0})
        M = m1 + m2
        self.nu = m1 * m2 / (M**2)
        self.delta = (m1 - m2) / M
        self.noneqcond = coord_greater_bound(self.delta, sp.sympify(1e-14)).subs(
            sp.Function("nrpyAbs"), sp.Abs
        )
        self.eqcond = 1 - self.noneqcond
        self.deltainvertible = self.delta * self.noneqcond + 1 * self.eqcond
        self.deltainv = 1 / self.deltainvertible
        self.chi_A = sp.Rational(1, 2) * (chi1 - chi2)
        self.chi_S = sp.Rational(1, 2) * (chi1 + chi2)
        self.vomega = self.Omega ** (sp.Rational(1, 3))
        self.vphi = self.Omega * (self.Omega_circ ** (-sp.Rational(2, 3)))
        self.vh3 = self.Omega * self.Hreal
        Heff = (self.Hreal**2 - 1) / (2 * self.nu) + 1
        EulerGamma = f2r(0.577215664901532860606512090082402431)
        eulerlog = sp.zeros(9)
        for m in range(1, 9):
            eulerlog[m] = EulerGamma + sp.log(2 * m * self.vomega)
        self.khat = sp.zeros(9)
        for m in range(9):
            self.khat[m] = m * self.Omega * self.Hreal
        self.deltalm["(5 , 5)"] = (
            (sp.Rational(96875, 131250) + sp.Rational(857528, 131250) * self.nu)
            * (self.vh3)
            / (1 - 2 * self.nu)
            + sp.Rational(3865, 429) * sp.pi * (self.vh3) ** 2
            + (
                (
                    sp.Rational(-7686949127, 31783752)
                    + sp.Rational(954500400, 31783752) * sp.pi**2
                )
                * (self.vh3) ** 3
            )
        )
        self.deltalm["(4 , 4)"] = (
            (sp.Rational(112, 120) + sp.Rational(219, 120) * self.nu) * (self.vh3)
            + sp.Rational(25136, 3465) * (self.vh3) ** 2
            + (sp.Rational(201088, 10395) * sp.pi**2 - sp.Rational(55144, 375))
            * (self.vh3) ** 3
        )
        self.deltalm["(4 , 3)"] = (
            (sp.Rational(4961, 810) * self.nu + sp.Rational(3, 5)) * (self.vh3)
        ) / (1 - 2 * self.nu) + sp.Rational(1571, 385) * sp.pi * (self.vh3) ** 2
        self.deltalm["(3 , 3)"] = (
            sp.Rational(13, 10) * (self.vh3)
            + sp.Rational(39, 7) * (self.vh3) ** 2
            + (-sp.Rational(227827, 3000) + sp.Rational(78, 7) * sp.pi**2)
            * (self.vh3) ** 3
            - sp.Rational(80897, 2430) * self.nu * self.vomega**5
        )
        self.deltalm["(3 , 2)"] = (
            ((sp.Rational(11, 5) * self.nu + sp.Rational(2, 3)) * (self.vh3))
            / (1 - 3 * self.nu)
            + sp.Rational(52, 21) * sp.pi * (self.vh3) ** 2
            + ((sp.Rational(208, 63) * sp.pi**2) - sp.Rational(9112, 405))
            * (self.vh3) ** 3
        )
        self.deltalm["(2 , 2)"] = (
            sp.Rational(7, 3) * self.vh3
            + (self.vh3) ** 2
            * (
                ((sp.Rational(8, 3) * self.nu - sp.Rational(4, 3)) * self.chi_S)
                - (sp.Rational(4, 3) * self.delta * self.chi_A)
                + (sp.Rational(428, 105) * sp.pi)
            )
            + (self.vh3) ** 3
            * ((sp.Rational(1712, 315) * sp.pi**2) - (sp.Rational(2203, 81)))
            - 24 * self.nu * self.vomega**5
        )
        self.deltalm["(2 , 1)"] = (
            sp.Rational(2, 3) * self.vh3
            + (sp.Rational(107, 105) * sp.pi) * (self.vh3) ** 2
            + ((sp.Rational(214, 315) * sp.pi**2) - (sp.Rational(272, 81)))
            * (self.vh3) ** 3
            - (sp.Rational(25, 2) * self.nu * self.vomega**5)
        )
        self.fspin["(5 , 5)"] = (
            (
                -70 * self.nu / (3 * (-1 + 2 * self.nu))
                + 110 * self.nu**2 / (3 * (-1 + 2 * self.nu))
                + 10 / (3 * (-1 + 2 * self.nu))
            )
            * (self.chi_A * self.deltainv)
            + (
                10 / (3 * (-1 + 2 * self.nu))
                - 10 * self.nu / (-1 + 2 * self.nu)
                + 10 * self.nu**2 / (-1 + 2 * self.nu)
            )
            * self.chi_S
        ) * self.vomega**3 + (
            sp.Rational(5, 2) * self.chi_S**2
            + (
                -5 / (-1 + 2 * self.nu)
                + 30 * self.nu / (-1 + 2 * self.nu)
                - 40 * self.nu**2 / (-1 + 2 * self.nu)
            )
            * self.chi_A
            * self.chi_S
            * self.deltainv
            + (
                -5 / (2 * (-1 + 2 * self.nu))
                + 15 * self.nu / (-1 + 2 * self.nu)
                - 20 * self.nu**2 / (-1 + 2 * self.nu)
            )
            * self.chi_A**2
        ) * self.vomega**4
        self.fspin_limit["(5 , 5)"] = (
            (
                sp.Rational(-70, 3) * self.nu / ((-1 + 2 * self.nu))
                + sp.Rational(110, 3) * self.nu**2 / ((-1 + 2 * self.nu))
                + sp.Rational(10, 3) / ((-1 + 2 * self.nu))
            )
            * (self.chi_A)
        ) * self.vomega**3 + (
            (
                -5 / (-1 + 2 * self.nu)
                + 30 * self.nu / (-1 + 2 * self.nu)
                - 40 * self.nu**2 / (-1 + 2 * self.nu)
            )
            * self.chi_A
            * self.chi_S
        ) * self.vomega**4
        self.fspin["(4 , 1)"] = (
            sp.Rational(5, 2)
            * self.nu
            * self.vomega
            / (1 - 2 * self.nu)
            * (self.chi_S - self.chi_A * self.deltainv)
        )
        self.fspin_limit["(4 , 1)"] = (
            sp.Rational(5, 2)
            * self.nu
            * self.vomega
            / (1 - 2 * self.nu)
            * (-self.chi_A)
        )
        self.fspin["(4 , 3)"] = (
            self.vomega
            / (1 - 2 * self.nu)
            * (
                sp.Rational(5, 2) * self.nu * self.chi_S
                - sp.Rational(5, 2) * self.nu * self.chi_A * self.deltainv
            )
            + self.vomega**3
            / (1 - 2 * self.nu)
            * (
                (sp.Rational(887, 44) * self.nu - sp.Rational(3143, 132) * self.nu**2)
                * self.chi_A
                * self.deltainv
                + (-sp.Rational(529, 132) * self.nu**2 - sp.Rational(667, 44) * self.nu)
                * self.chi_S
            )
            + self.vomega**4
            / (1 - 2 * self.nu)
            * (
                (12 * self.nu**2 - sp.Rational(37, 3) * self.nu + sp.Rational(3, 2))
                * self.chi_A**2
                + (sp.Rational(137, 6) * self.nu**2 - 18 * self.nu + 3)
                * self.chi_A
                * self.chi_S
                * self.deltainv
                + (
                    sp.Rational(35, 6) * self.nu**2
                    + sp.Rational(1, 3) * self.nu
                    + sp.Rational(3, 2)
                )
                * self.chi_S**2
            )
        )
        self.fspin_limit["(4 , 3)"] = (
            self.vomega
            / (1 - 2 * self.nu)
            * (-sp.Rational(5, 2) * self.nu * self.chi_A)
            + self.vomega**3
            / (1 - 2 * self.nu)
            * (
                (sp.Rational(887, 44) * self.nu - sp.Rational(3143, 132) * self.nu**2)
                * self.chi_A
            )
            + self.vomega**4
            / (1 - 2 * self.nu)
            * (
                (sp.Rational(137, 6) * self.nu**2 - 18 * self.nu + 3)
                * self.chi_A
                * self.chi_S
            )
        )
        self.fspin["(3 , 1)"] = (
            self.vomega**3
            * (
                self.chi_A * self.deltainv * (-4 + 11 * self.nu)
                + self.chi_S * (-4 + 13 * self.nu)
            )
            / (2)
        )
        self.fspin_limit["(3 , 1)"] = -self.chi_A * sp.Rational(5, 8) * self.vomega**3
        self.fspin["(3 , 3)"] = (
            self.vomega**3
            * (
                ((sp.Rational(19, 2) * self.nu - 2) * self.chi_A * self.deltainv)
                + ((sp.Rational(5, 2) * self.nu - 2) * self.chi_S)
            )
            + self.vomega**4
            * (
                (sp.Rational(3, 2) - 6 * self.nu) * self.chi_A**2
                + (3 - 12 * self.nu) * (self.chi_A * self.chi_S * self.deltainv)
                + sp.Rational(3, 2) * self.chi_S**2
            )
            + self.vomega**5
            * (
                (
                    (
                        sp.Rational(407, 30) * self.nu**2
                        - sp.Rational(593, 60) * self.nu
                        + sp.Rational(2, 3)
                    )
                    * self.chi_A
                    * self.deltainv
                )
                + (
                    (
                        sp.Rational(241, 30) * self.nu**2
                        + sp.Rational(11, 20) * self.nu
                        + sp.Rational(2, 3)
                    )
                    * self.chi_S
                )
            )
            + self.vomega**6
            * (
                (-12 * self.nu**2 + sp.Rational(11, 2) * self.nu - sp.Rational(7, 4))
                * self.chi_A**2
                + (44 * self.nu**2 - self.nu - sp.Rational(7, 2))
                * (self.chi_A * self.chi_S * self.deltainv)
                + (6 * self.nu**2 - sp.Rational(27, 2) * self.nu - sp.Rational(7, 4))
                * self.chi_S**2
            )
        )
        self.fspinimag = (self.vh3) ** 2 * (
            sp.Rational(7339, 540) * self.nu - sp.Rational(81, 20)
        ) * self.chi_A * self.deltainv + (
            sp.Rational(593, 108) * self.nu - sp.Rational(81, 20)
        ) * self.chi_S
        self.fspin_limit["(3 , 3)"] = (
            self.vomega**3 * ((sp.Rational(19, 2) * self.nu - 2) * self.chi_A)
            + self.vomega**4 * ((3 - 12 * self.nu) * (self.chi_A * self.chi_S))
            + self.vomega**5
            * (
                (
                    sp.Rational(407, 30) * self.nu**2
                    - sp.Rational(593, 60) * self.nu
                    + sp.Rational(2, 3)
                )
                * self.chi_A
            )
            + self.vomega**6
            * (
                (44 * self.nu**2 - self.nu - sp.Rational(7, 2))
                * (self.chi_A * self.chi_S)
            )
        )
        self.fspinimag_limit = (
            (self.vh3) ** 2
            * (sp.Rational(7339, 540) * self.nu - sp.Rational(81, 20))
            * self.chi_A
        )
        self.fspin["(2 , 1)"] = (
            -sp.Rational(3, 2) * self.vomega * (self.chi_A * self.deltainv + self.chi_S)
            + self.vomega**3
            * (
                (sp.Rational(131, 84) * self.nu + sp.Rational(61, 12))
                * (self.chi_A * self.deltainv)
                + (sp.Rational(79, 84) * self.nu + sp.Rational(61, 12)) * self.chi_S
            )
            + self.vomega**4
            * (
                (-2 * self.nu - 3) * self.chi_A**2
                + (sp.Rational(21, 2) * self.nu - 6)
                * (self.chi_A * self.chi_S * self.deltainv)
                + (sp.Rational(1, 2) * self.nu - 3) * self.chi_S**2
            )
            + self.vomega**5
            * (
                (
                    sp.Rational(-703, 112) * self.nu**2
                    + sp.Rational(8797, 1008) * self.nu
                    - sp.Rational(81, 16)
                )
                * (self.chi_A * self.deltainv)
                + (
                    sp.Rational(613, 1008) * self.nu**2
                    + sp.Rational(1709, 1008) * self.nu
                    - sp.Rational(81, 16)
                )
                * self.chi_S
                + (sp.Rational(3, 4) - 3 * self.nu) * (self.chi_A**3 * self.deltainv)
                + (sp.Rational(9, 4) - 6 * self.nu)
                * (self.chi_A * self.chi_S**2 * self.deltainv)
                + (sp.Rational(9, 4) - 3 * self.nu) * (self.chi_A**2 * self.chi_S)
                + sp.Rational(3, 4) * self.chi_S**3
            )
            + self.vomega**6
            * (
                (
                    sp.Rational(5, 7) * self.nu**2
                    - sp.Rational(9287, 1008) * self.nu
                    + sp.Rational(4163, 252)
                )
                * self.chi_A**2
                + (
                    sp.Rational(139, 72) * self.nu**2
                    - sp.Rational(2633, 1008) * self.nu
                    + sp.Rational(4163, 252)
                )
                * self.chi_S**2
                + (
                    sp.Rational(9487, 504) * self.nu**2
                    - sp.Rational(1636, 21) * self.nu
                    + sp.Rational(4163, 126)
                )
                * (self.chi_A * self.chi_S * self.deltainv)
            )
        )
        self.fspin_limit["(2 , 1)"] = (
            -sp.Rational(3, 2) * self.vomega * (self.chi_A)
            + self.vomega**3
            * ((sp.Rational(131, 84) * self.nu + sp.Rational(61, 12)) * (self.chi_A))
            + self.vomega**4
            * (+(sp.Rational(21, 2) * self.nu - 6) * (self.chi_A * self.chi_S))
            + self.vomega**5
            * (
                (
                    sp.Rational(-703, 112) * self.nu**2
                    + sp.Rational(8797, 1008) * self.nu
                    - sp.Rational(81, 16)
                )
                * (self.chi_A)
                + (sp.Rational(3, 4) - 3 * self.nu) * (self.chi_A**3)
                + (sp.Rational(9, 4) - 6 * self.nu) * (self.chi_A * self.chi_S**2)
            )
            + self.vomega**6
            * (
                (
                    sp.Rational(9487, 504) * self.nu**2
                    - sp.Rational(1636, 21) * self.nu
                    + sp.Rational(4163, 126)
                )
                * (self.chi_A * self.chi_S)
            )
        )
        self.rho["(8 , 8)"] = (
            1
            + (
                (
                    3482
                    - 26778 * self.nu
                    + 64659 * self.nu**2
                    - 53445 * self.nu**3
                    + 12243 * self.nu**4
                )
                / (2736 * (-1 + 7 * self.nu - 14 * self.nu**2 + 7 * self.nu**3))
            )
            * self.vomega**2
        )
        self.rho["(8 , 7)"] = (
            1
            + (
                (
                    23478
                    - 154099 * self.nu
                    + 309498 * self.nu**2
                    - 207550 * self.nu**3
                    + 38920 * self.nu**4
                )
                / (18240 * (-1 + 6 * self.nu - 10 * self.nu**2 + 4 * self.nu**3))
            )
            * self.vomega**2
        )
        self.rho["(8 , 6)"] = (
            1
            + (
                (
                    1002
                    - 7498 * self.nu
                    + 17269 * self.nu**2
                    - 13055 * self.nu**3
                    + 2653 * self.nu**4
                )
                / (912 * (-1 + 7 * self.nu - 14 * self.nu**2 + 7 * self.nu**3))
            )
            * self.vomega**2
        )
        self.rho["(8 , 5)"] = (
            1
            + (
                (
                    4350
                    - 28055 * self.nu
                    + 54642 * self.nu**2
                    - 34598 * self.nu**3
                    + 6056 * self.nu**4
                )
                / (3648 * (-1 + 6 * self.nu - 10 * self.nu**2 + 4 * self.nu**3))
            )
            * self.vomega**2
        )
        self.rho["(8 , 4)"] = (
            1
            + (
                (
                    2666
                    - 19434 * self.nu
                    + 42627 * self.nu**2
                    - 28965 * self.nu**3
                    + 4899 * self.nu**4
                )
                / (2736 * (-1 + 7 * self.nu - 14 * self.nu**2 + 7 * self.nu**3))
            )
            * self.vomega**2
        )
        self.rho["(8 , 3)"] = (
            1
            + (
                (
                    20598
                    - 131059 * self.nu
                    + 249018 * self.nu**2
                    - 149950 * self.nu**3
                    + 24520 * self.nu**4
                )
                / (18240 * (-1 + 6 * self.nu - 10 * self.nu**2 + 4 * self.nu**3))
            )
            * self.vomega**2
        )
        self.rho["(8 , 2)"] = (
            1
            + (
                (
                    2462
                    - 17598 * self.nu
                    + 37119 * self.nu**2
                    - 22845 * self.nu**3
                    + 3063 * self.nu**4
                )
                / (2736 * (-1 + 7 * self.nu - 14 * self.nu**2 + 7 * self.nu**3))
            )
            * self.vomega**2
        )
        self.rho["(8 , 1)"] = (
            1
            + (
                (
                    20022
                    - 126451 * self.nu
                    + 236922 * self.nu**2
                    - 138430 * self.nu**3
                    + 21640 * self.nu**4
                )
                / (18240 * (-1 + 6 * self.nu - 10 * self.nu**2 + 4 * self.nu**3))
            )
            * self.vomega**2
        )
        self.rho["(7 , 7)"] = (
            1
            + (
                (1380 * self.nu**3 - 4963 * self.nu**2 + 4246 * self.nu - 906)
                / (714 * (3 * self.nu**2 - 4 * self.nu + 1))
            )
            * self.vomega**2
        )
        self.rho["(7 , 6)"] = (
            1
            + (
                (
                    6104 * self.nu**4
                    - 29351 * self.nu**3
                    + 37828 * self.nu**2
                    - 16185 * self.nu
                    + 2144
                )
                / (1666 * (7 * self.nu**3 - 14 * self.nu**2 + 7 * self.nu - 1))
            )
            * self.vomega**2
        )
        self.rho["(7 , 5)"] = (
            1
            + (
                (804 * self.nu**3 - 3523 * self.nu**2 + 3382 * self.nu - 762)
                / (714 * (3 * self.nu**2 - 4 * self.nu + 1))
            )
            * self.vomega**2
        )
        self.rho["(7 , 4)"] = (
            1
            + (
                (
                    41076 * self.nu**4
                    - 217959 * self.nu**3
                    + 298872 * self.nu**2
                    - 131805 * self.nu
                    + 17756
                )
                / (14994 * (7 * self.nu**3 - 14 * self.nu**2 + 7 * self.nu - 1))
            )
            * self.vomega**2
        )
        self.rho["(7 , 3)"] = (
            1
            + (
                (420 * self.nu**3 - 2563 * self.nu**2 + 2806 * self.nu - 666)
                / (714 * (3 * self.nu**2 - 4 * self.nu + 1))
            )
            * self.vomega**2
        )
        self.rho["(7 , 2)"] = (
            1
            + (
                (
                    32760 * self.nu**4
                    - 190239 * self.nu**3
                    + 273924 * self.nu**2
                    - 123489 * self.nu
                    + 16832
                )
                / (14994 * (7 * self.nu**3 - 14 * self.nu**2 + 7 * self.nu - 1))
            )
            * self.vomega**2
        )
        self.rho["(7 , 1)"] = (
            1
            + (
                (228 * self.nu**3 - 2083 * self.nu**2 + 2518 * self.nu - 618)
                / (714 * (3 * self.nu**2 - 4 * self.nu + 1))
            )
            * self.vomega**2
        )
        self.rho["(6 , 6)"] = (
            1
            + (
                (273 * self.nu**3 - 861 * self.nu**2 + 602 * self.nu - 106)
                / (84 * (5 * self.nu**2 - 5 * self.nu + 1))
            )
            * self.vomega**2
            - sp.Rational(1025435, 659736) * self.vomega**4
        )
        self.rho["(6 , 5)"] = (
            1
            + (
                (220 * self.nu**3 - 910 * self.nu**2 + 838 * self.nu - 185)
                / (144 * (3 * self.nu**2 - 4 * self.nu + 1))
            )
            * self.vomega**2
        )
        self.rho["(6 , 4)"] = (
            1
            + self.vomega**2
            * (
                (133 * self.nu**3 - 581 * self.nu**2 + 462 * self.nu - 86)
                / (84 * (5 * self.nu**2 - 5 * self.nu + 1))
            )
            - sp.Rational(476887, 659736) * self.vomega**4
        )
        self.rho["(6 , 3)"] = (
            1
            + (
                (156 * self.nu**3 - 750 * self.nu**2 + 742 * self.nu - 169)
                / (144 * (3 * self.nu**2 - 4 * self.nu + 1))
            )
            * self.vomega**2
        )
        self.rho["(6 , 2)"] = (
            1
            + (
                (49 * self.nu**3 - 413 * self.nu**2 + 378 * self.nu - 74)
                / (84 * (5 * self.nu**2 - 5 * self.nu + 1))
            )
            * self.vomega**2
            - sp.Rational(817991, 3298680) * self.vomega**4
        )
        self.rho["(6 , 1)"] = (
            1
            + (
                (124 * self.nu**3 - 670 * self.nu**2 + 694 * self.nu - 161)
                / (144 * (3 * self.nu**2 - 4 * self.nu + 1))
            )
            * self.vomega**2
        )
        self.rho["(5 , 5)"] = (
            1
            + self.vomega**2
            * (
                sp.Rational(487, 390) / (-1 + 2 * self.nu)
                - sp.Rational(649, 195) * self.nu / (-1 + 2 * self.nu)
                + sp.Rational(256, 195) * self.nu**2 / (-1 + 2 * self.nu)
            )
            - sp.Rational(3353747, 2129400) * self.vomega**4
            + (
                sp.Rational(190606537999247, 11957879934000)
                - sp.Rational(1546, 429) * eulerlog[5]
            )
            * self.vomega**6
            + (
                -sp.Rational(1213641959949291437, 118143853747920000)
                + sp.Rational(376451, 83655) * eulerlog[5]
            )
            * self.vomega**8
            + (
                -sp.Rational(150082616449726042201261, 4837990810977324000000)
                + sp.Rational(2592446431, 456756300) * eulerlog[5]
            )
            * self.vomega**10
            + self.nu
            * (
                -f2r(2.61) * self.vomega**4
                + f2r(1.25) * self.vomega**6
                + -f2r(35.7) * self.vomega**8
            )
        )
        self.rho["(5 , 4)"] = (
            1
            + (
                (33320 * self.nu**3 - 127610 * self.nu**2 + 96019 * self.nu - 17448)
                / (13650 * (5 * self.nu**2 - 5 * self.nu + 1))
            )
            * self.vomega**2
            - sp.Rational(16213384, 15526875) * self.vomega**4
        )
        self.rho["(5 , 3)"] = (
            1
            + ((176 * self.nu**2 - 850 * self.nu + 375) / (390 * (2 * self.nu - 1)))
            * self.vomega**2
            - sp.Rational(410833, 709800) * self.vomega**4
        )
        self.rho["(5 , 2)"] = (
            1
            + (
                (21980 * self.nu**3 - 104930 * self.nu**2 + 84679 * self.nu - 15828)
                / (13650 * (5 * self.nu**2 - 5 * self.nu + 1))
            )
            * self.vomega**2
            - sp.Rational(7187914, 15526875) * self.vomega**4
        )
        self.rho["(5 , 1)"] = (
            1
            + ((8 * self.nu**2 - 626 * self.nu + 319) / (390 * (2 * self.nu - 1)))
            * self.vomega**2
            - sp.Rational(31877, 304200) * self.vomega**4
        )
        self.rho["(4 , 1)"] = (
            1
            + ((288 * self.nu**2 - 1385 * self.nu + 602) / (528 * (2 * self.nu - 1)))
            * self.vomega**2
            - (sp.Rational(7775491, 21141120)) * self.vomega**4
            + (
                sp.Rational(1227423222031, 1758095539200)
                - sp.Rational(1571, 6930) * eulerlog[1]
            )
            * self.vomega**6
        )
        self.rho["(4 , 2)"] = (
            1
            + self.vomega**2
            * ((1146 - 3530 * self.nu + 285 * self.nu**2) / (1320 * (-1 + 3 * self.nu)))
            + self.vomega**3
            * (
                (
                    self.chi_A * (10 - 21 * self.nu) * self.delta
                    + self.chi_S * (10 - 59 * self.nu + 78 * self.nu**2)
                )
                / (15 * (-1 + 3 * self.nu))
            )
            + self.vomega**4
            * (
                (
                    -114859044
                    + 295834536 * self.nu
                    + 1204388696 * self.nu**2
                    - 3047981160 * self.nu**3
                    - 379526805 * self.nu**4
                )
                / (317116800 * (-1 + 3 * self.nu) ** 2)
            )
            + self.vomega**6
            * (848238724511 / 219761942400 - (3142 / 3465) * eulerlog[2])
        )
        self.rho["(4 , 3)"] = (
            1
            + self.vomega**2
            / (1 - 2 * self.nu)
            * (
                -sp.Rational(10, 11) * self.nu**2
                + sp.Rational(547, 176) * self.nu
                - sp.Rational(111, 88)
            )
            - sp.Rational(6894273, 7047040) * self.vomega**4
            + self.vomega**6
            * (
                sp.Rational(1664224207351, 195343948800)
                - sp.Rational(1571, 770) * eulerlog[3]
            )
            + self.vomega**8
            * (
                -sp.Rational(2465107182496333, 460490801971200)
                + sp.Rational(174381, 67760) * eulerlog[3]
            )
            + self.nu
            * (
                -f2r(0.654) * self.vomega**4
                + -f2r(3.69) * self.vomega**6
                + f2r(18.5) * self.vomega**8
            )
        )
        self.rho["(4 , 4)"] = (
            1
            + self.vomega**2
            * (
                (1614 - 5870 * self.nu + 2625 * self.nu**2)
                / (1320 * (-1 + 3 * self.nu))
            )
            + self.vomega**3
            * (
                self.chi_A * (10 - 39 * self.nu) * self.delta
                + self.chi_S * (10 - 41 * self.nu + 42 * self.nu**2)
            )
            / (15 * (-1 + 3 * self.nu))
            + self.vomega**4
            * (
                (
                    -511573572
                    + 2338945704 * self.nu
                    - 313857376 * self.nu**2
                    - 6733146000 * self.nu**3
                    + 1252563795 * self.nu**4
                )
                / (317116800 * (-1 + 3 * self.nu) ** 2)
                + self.chi_S**2 / 2
                + self.delta * self.chi_S * self.chi_A
                + self.delta**2 * self.chi_A**2 / 2
            )
            + self.vomega**5
            * (
                self.chi_A
                * self.delta
                * (-8280 + 42716 * self.nu - 57990 * self.nu**2 + 8955 * self.nu**3)
                / (6600 * (-1 + 3 * self.nu) ** 2)
                + self.chi_S
                * (
                    -8280
                    + 66284 * self.nu
                    - 176418 * self.nu**2
                    + 128085 * self.nu**3
                    + 88650 * self.nu**2 * self.nu**2
                )
                / (6600 * (-1 + 3 * self.nu) ** 2)
            )
            + self.vomega**6
            * (
                sp.Rational(16600939332793, 1098809712000)
                - sp.Rational(12568, 3465) * eulerlog[4]
            )
            + self.vomega**8
            * (
                -sp.Rational(172066910136202271, 19426955708160000)
                + sp.Rational(845198, 190575) * eulerlog[4]
            )
            + self.vomega**10
            * (
                -sp.Rational(17154485653213713419357, 568432724020761600000)
                + sp.Rational(22324502267, 3815311500) * eulerlog[4]
            )
            + self.nu
            * (
                -f2r(3.56) * self.vomega**6
                + f2r(15.6) * self.vomega**8
                + -216 * self.vomega**10
            )
        )
        self.rho["(3 , 1)"] = (
            1
            - self.vomega**2 * (sp.Rational(2, 9) * self.nu + sp.Rational(13, 18))
            + self.vomega**4
            * (
                -sp.Rational(829, 1782) * self.nu**2
                - sp.Rational(1685, 1782) * self.nu
                + sp.Rational(101, 7128)
            )
            + self.vomega**6
            * (sp.Rational(11706720301, 6129723600) - sp.Rational(26, 63) * eulerlog[1])
            + self.vomega**8
            * (
                sp.Rational(169, 567) * eulerlog[1]
                + sp.Rational(2606097992581, 4854741091200)
            )
        )
        self.rho["(3 , 2)"] = (
            1
            + self.vomega * ((4 * self.nu * self.chi_S) / (3 * (1 - 3 * self.nu)))
            + self.vomega**2
            * (
                (
                    sp.Rational(-32, 27) * self.nu**2
                    + sp.Rational(223, 54) * self.nu
                    - sp.Rational(164, 135)
                )
                / (1 - 3 * self.nu)
                - (16 * self.nu**2 * self.chi_S**2) / (9 * (1 - 3 * self.nu) ** 2)
            )
            + self.vomega**3
            * (
                (sp.Rational(13, 9) * self.nu + sp.Rational(2, 9))
                * (self.delta * self.chi_A)
                / (1 - 3 * self.nu)
                + (
                    sp.Rational(607, 81) * self.nu**3
                    + sp.Rational(503, 81) * self.nu**2
                    - sp.Rational(1478, 405) * self.nu
                    + sp.Rational(2, 9)
                )
                * self.chi_S
                / (1 - 3 * self.nu) ** 2
                + (320 * self.nu**3 * self.chi_S**3) / (81 * (1 - 3 * self.nu) ** 3)
            )
            + self.vomega**4
            * (
                (
                    sp.Rational(77141, 40095) * self.nu**4
                    - sp.Rational(508474, 40095) * self.nu**3
                    - sp.Rational(945121, 320760) * self.nu**2
                    + sp.Rational(1610009, 320760) * self.nu
                    - sp.Rational(180566, 200475)
                )
                / (1 - 3 * self.nu) ** 2
                + (4 * self.nu**2 - 3 * self.nu + sp.Rational(1, 3))
                * (self.chi_A**2)
                / (1 - 3 * self.nu)
                + (
                    sp.Rational(-50, 27) * self.nu**2
                    - sp.Rational(88, 27) * self.nu
                    + sp.Rational(2, 3)
                )
                * (self.delta * self.chi_A * self.chi_S)
                / (1 - 3 * self.nu) ** 2
                + (
                    sp.Rational(-2452, 243) * self.nu**4
                    - sp.Rational(1997, 243) * self.nu**3
                    + sp.Rational(1435, 243) * self.nu**2
                    - sp.Rational(43, 27) * self.nu
                    + sp.Rational(1, 3)
                )
                * (self.chi_S**2)
                / ((1 - 3 * self.nu) ** 3)
            )
            + self.vomega**5
            * (
                (
                    sp.Rational(-1184225, 96228) * self.nu**5
                    - sp.Rational(40204523, 962280) * self.nu**4
                    + sp.Rational(101706029, 962280) * self.nu**3
                    - sp.Rational(14103833, 192456) * self.nu**2
                    + sp.Rational(20471053, 962280) * self.nu
                    - sp.Rational(2788, 1215)
                )
                * self.chi_S
                / (1 - 3 * self.nu) ** 3
                + (
                    sp.Rational(608, 81) * self.nu**3
                    + sp.Rational(736, 81) * self.nu**2
                    - sp.Rational(16, 9) * self.nu
                )
                * (self.delta * self.chi_A * self.chi_S**2)
                / (1 - 3 * self.nu) ** 3
                + (
                    sp.Rational(889673, 106920) * self.nu**3
                    - sp.Rational(75737, 5346) * self.nu**2
                    + sp.Rational(376177, 35640) * self.nu
                    - sp.Rational(2788, 1215)
                )
                * (self.delta * self.chi_A)
                / (1 - 3 * self.nu) ** 2
                + (
                    sp.Rational(96176, 2187) * self.nu**5
                    + sp.Rational(43528, 2187) * self.nu**4
                    - sp.Rational(40232, 2187) * self.nu**3
                    + sp.Rational(376, 81) * self.nu**2
                    - sp.Rational(8, 9) * self.nu
                )
                * (self.chi_S**3)
                / (1 - 3 * self.nu) ** 4
                + (
                    sp.Rational(-32, 3) * self.nu**3
                    + 8 * self.nu**2
                    - sp.Rational(8, 9) * self.nu
                )
                * (self.chi_A**2 * self.chi_S)
                / ((1 - 3 * self.nu) ** 2)
            )
            + self.vomega**6
            * (sp.Rational(5849948554, 940355325) - sp.Rational(104, 63) * eulerlog[2])
            + self.vomega**8
            * (
                sp.Rational(17056, 8505) * eulerlog[2]
                - sp.Rational(10607269449358, 3072140846775)
            )
            + self.vomega**10
            * (
                -sp.Rational(1312549797426453052, 176264081083715625)
                + sp.Rational(18778864, 12629925) * eulerlog[2]
            )
            + self.nu
            * (
                +f2r(0.333) * self.vomega**6
                - f2r(6.5) * self.vomega**8
                + 98 * self.vomega**10
            )
        )
        self.rho["(3 , 3)"] = (
            1
            + self.vomega**2 * (sp.Rational(2, 3) * self.nu - sp.Rational(7, 6))
            + self.vomega**4
            * (
                -sp.Rational(6719, 3960)
                - sp.Rational(1861, 990) * self.nu
                + sp.Rational(149, 330) * self.nu**2
            )
            + self.vomega**6
            * (
                sp.Rational(3203101567, 227026800)
                + (-sp.Rational(129509, 25740) + sp.Rational(41, 192) * sp.pi**2)
                * self.nu
                - sp.Rational(274621, 154440) * self.nu**2
                + sp.Rational(12011, 46332) * self.nu**3
                - sp.Rational(26, 7) * eulerlog[3]
            )
            + self.vomega**8
            * (-sp.Rational(57566572157, 8562153600) + sp.Rational(13, 3) * eulerlog[3])
            + self.vomega**10
            * (
                -sp.Rational(903823148417327, 30566888352000)
                + sp.Rational(87347, 13860) * eulerlog[3]
            )
            + self.nu * (12 * self.vomega**8 + -215 * self.vomega**10)
        )
        self.rho["(2 , 2)"] = (
            1
            + self.vomega**2 * (sp.Rational(55, 84) * self.nu - sp.Rational(43, 42))
            + self.vomega**3
            * ((-2 * (self.chi_S + self.chi_A * self.delta - self.chi_S * self.nu)) / 3)
            + self.vomega**4
            * (
                sp.Rational(19583, 42336) * self.nu**2
                - sp.Rational(33025, 21168) * self.nu
                - sp.Rational(20555, 10584)
                + (sp.Rational(1, 2) - 2 * self.nu) * self.chi_A**2
                + self.delta * self.chi_A * self.chi_S
                + sp.Rational(1, 2) * self.chi_S**2
            )
            + self.vomega**5
            * (
                self.delta
                * (-sp.Rational(19, 42) * self.nu - sp.Rational(34, 21))
                * self.chi_A
                + (
                    sp.Rational(209, 126) * self.nu**2
                    + sp.Rational(49, 18) * self.nu
                    - sp.Rational(34, 21)
                )
                * self.chi_S
            )
            + self.vomega**6
            * (
                sp.Rational(10620745, 39118464) * self.nu**3
                - sp.Rational(6292061, 3259872) * self.nu**2
                + sp.Rational(41, 192) * sp.pi**2 * self.nu
                - sp.Rational(48993925, 9779616) * self.nu
                - sp.Rational(428, 105) * eulerlog[2]
                + sp.Rational(1556919113, 122245200)
                + self.delta
                * (sp.Rational(89, 126) - sp.Rational(781, 252) * self.nu)
                * self.chi_A
                * self.chi_S
                + (
                    -sp.Rational(27, 14) * self.nu**2
                    - sp.Rational(457, 504) * self.nu
                    + sp.Rational(89, 252)
                )
                * self.chi_A**2
                + (
                    sp.Rational(10, 9) * self.nu**2
                    - sp.Rational(1817, 504) * self.nu
                    + sp.Rational(89, 252)
                )
                * self.chi_S**2
            )
            + self.vomega**7
            * (
                self.delta
                * (
                    sp.Rational(97865, 63504) * self.nu**2
                    + sp.Rational(50140, 3969) * self.nu
                    + sp.Rational(18733, 15876)
                )
                * self.chi_A
                + (
                    sp.Rational(50803, 63504) * self.nu**3
                    - sp.Rational(245717, 63504) * self.nu**2
                    + sp.Rational(74749, 5292) * self.nu
                    + sp.Rational(18733, 15876)
                )
                * self.chi_S
                + self.delta
                * self.chi_A**3
                * (sp.Rational(1, 3) - sp.Rational(4, 3) * self.nu)
                + self.delta * (2 * self.nu + 1) * self.chi_A * self.chi_S**2
                + (
                    sp.Rational(-4, 1) * self.nu**2
                    - sp.Rational(3, 1) * self.nu
                    + sp.Rational(1, 1)
                )
                * self.chi_A**2
                * self.chi_S
                + (self.nu + sp.Rational(1, 3)) * self.chi_S**3
            )
            + self.vomega**8
            * (
                sp.Rational(9202, 2205) * eulerlog[2]
                - sp.Rational(387216563023, 160190110080)
            )
            + self.vomega**10
            * (
                sp.Rational(439877, 55566) * eulerlog[2]
                - sp.Rational(16094530514677, 533967033600)
            )
            + self.nu * (f2r(21.2) * self.vomega**8 + -411 * self.vomega**10)
        )
        self.rho["(2 , 1)"] = (
            1
            + self.vomega**2 * (sp.Rational(23, 84) * self.nu - sp.Rational(59, 56))
            + self.vomega**4
            * (
                sp.Rational(617, 4704) * self.nu**2
                - sp.Rational(10993, 14112) * self.nu
                - sp.Rational(47009, 56448)
            )
            + self.vomega**6
            * (
                sp.Rational(7613184941, 2607897600)
                - sp.Rational(107, 105) * eulerlog[1]
            )
            + self.vomega**8
            * (
                -sp.Rational(1168617463883, 911303737344)
                + sp.Rational(6313, 5880) * eulerlog[1]
            )
            + self.vomega**10
            * (
                -sp.Rational(63735873771463, 16569158860800)
                + sp.Rational(5029963, 5927040) * eulerlog[1]
            )
            + self.nu
            * (
                f2r(1.65) * self.vomega**6
                + f2r(26.5) * self.vomega**8
                + 80 * self.vomega**10
            )
        )
        self.Y = [
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            [
                0,
                f2r(0.3454941494713355),
                0,
                f2r(0.3231801841141506),
                0,
                f2r(0.32028164857621516),
                0,
                f2r(0.31937046138540076),
                0,
            ],
            [
                0,
                0,
                f2r(0.3862742020231896),
                0,
                f2r(0.33452327177864466),
                0,
                f2r(0.32569524293385776),
                0,
                f2r(0.32254835519288305),
            ],
            [
                0,
                0,
                0,
                f2r(0.4172238236327842),
                0,
                f2r(0.34594371914684025),
                0,
                f2r(0.331899519333737),
                0,
            ],
            [
                0,
                0,
                0,
                0,
                f2r(0.4425326924449826),
                0,
                f2r(0.3567812628539981),
                0,
                f2r(0.3382915688890245),
            ],
            [
                0,
                0,
                0,
                0,
                0,
                f2r(0.46413220344085826),
                0,
                f2r(0.3669287245764378),
                0,
            ],
            [0, 0, 0, 0, 0, 0, f2r(0.48308411358006625), 0, f2r(0.3764161087284946)],
            [0, 0, 0, 0, 0, 0, 0, f2r(0.5000395635705508), 0],
            [0, 0, 0, 0, 0, 0, 0, 0, f2r(0.5154289843972844)],
        ]
        self.c = sp.zeros(10)
        for k in range(2, 10):
            self.c[k] = (m1 / M) ** (k - 1) + ((-1) ** k) * ((m2 / M) ** (k - 1))
        self.c[3] *= self.noneqcond
        self.c[3] += -1 * self.eqcond
        self.c[5] *= self.noneqcond
        self.c[5] += sp.Rational(-1, 2) * self.eqcond
        self.r0 = 2 / sp.sqrt(sp.E)
        self.Tlm = sp.zeros(9, 9)
        tlmprod_fac = 1
        for m in range(1, 9):
            hathatk = self.khat[m]
            hathatksq4 = 4 * hathatk * hathatk
            hathatk4pi = 4 * sp.pi * hathatk
            tlmprod_fac = 1
            tlmprefac = sp.sqrt(hathatk4pi / (1 - sp.exp(-hathatk4pi)))
            for j in range(1, 9):
                z2 = sp.factorial(j)
                tlmprod_fac *= hathatksq4 + j**2
                if m > j:
                    continue
                self.Tlm[j, m] = tlmprefac * sp.sqrt(tlmprod_fac) / z2
        self.effective_source = [Heff, self.vomega * self.pphi]

    def flux(self) -> sp.Mul:
        """
        Compute the SEOBNRv5 aligned-spin flux.

        This function returns the symbolic expression for the factorized resummed flux
        in the SEOBNRv5 aligned spin approximant.

        :return: The symbolic expression representing the flux.
        """
        factorized_flux = sp.sympify(0)
        for l in range(2, 9):
            for m in range(1, l + 1):
                if not (l + m) % 2:
                    n_complex = (
                        8 * sp.pi * (sp.I * m) ** l / sp.factorial2(2 * l + 1)
                    ) * sp.sqrt(((l + 1) * (l + 2)) / ((l) * (l - 1)))
                    n_abs = sp.Abs(n_complex)
                else:
                    n_complex = (
                        -16 * sp.I * sp.pi * (sp.I * m) ** l / sp.factorial2(2 * l + 1)
                    ) * sp.sqrt(
                        ((2 * l + 1) * (l + 2) * (l**2 - m**2))
                        / ((2 * l - 1) * (l + 1) * (l) * (l - 1))
                    )
                    n_abs = sp.Abs(n_complex)
                pn_contribution_f = 0
                if not m % 2:
                    pn_contribution_f += self.rho[f"({l} , {m})"] ** l
                if m % 2:
                    pn_contribution_f += self.rho[f"({l} , {m})"] ** l
                    if (l < 5) or (l == 5 and m == 5):
                        pn_contribution_f = (
                            self.noneqcond
                            * (
                                self.rho[f"({l} , {m})"] ** l
                                + self.fspin[f"({l} , {m})"]
                            )
                            + self.eqcond * self.fspin_limit[f"({l} , {m})"]
                        )

                tail_term = self.Tlm[l, m]
                newtonian_strain_amplitude = (
                    self.nu
                    * n_abs
                    * self.c[l + (l + m) % 2]
                    * (self.vphi ** (l + (l + m) % 2))
                    * self.Y[m][l - (l + m) % 2]
                )
                strain_amplitude = (
                    newtonian_strain_amplitude
                    * self.effective_source[(l + m) % 2]
                    * tail_term
                    * pn_contribution_f
                )
                factorized_flux += m * m * strain_amplitude**2
        factorized_flux *= -(sp.Rational(1, 8) * self.Omega**2 / sp.pi)
        return factorized_flux / self.nu

    def strain(self) -> Dict[str, sp.core.Mul]:
        """
        Compute the SEOBNRv5 aligned-spin gravitational-wave strain.

        This function returns the symbolic expression for the factorized resummed strains
        in the SEOBNRv5 aligned spin approximant.

        :return: dictionary containing the symbolic expressions representing each mode of the strain
        """
        hlms = {}
        # modes = [(2, 2)]
        l, m = 2, 2
        pn_contribution_f = self.rho["(2 , 2)"] ** l
        pn_contribution_delta = sp.exp(sp.I * self.deltalm["(2 , 2)"])
        pn_contribution = pn_contribution_f * pn_contribution_delta
        gamma_term = sp.Symbol(f"gamma_{l}{m}")
        khat2 = self.khat[2]
        tail_prefactor = (
            sp.exp(sp.pi * khat2)
            * (sp.exp(2 * sp.I * khat2 * sp.log(2 * 2 * self.Omega * self.r0)))
            / sp.factorial(2)
        )
        tail_term = gamma_term * tail_prefactor
        non_newtonian_contribution = tail_term * pn_contribution
        n = (8 * sp.I * sp.pi * (sp.I * 2) ** 2 / sp.factorial2(2 * 2 + 1)) * sp.sqrt(
            ((2 + 1) * (2 + 2)) / ((2) * (2 - 1))
        )
        newtonian_strain = (
            self.nu
            * n
            * self.c[2]
            * self.vphi**2
            * self.Y[2][2]
            * sp.exp(-2 * sp.I * self.phi)
        )
        hlms_no_source = newtonian_strain * non_newtonian_contribution
        hlms[f"({l} , {m})"] = hlms_no_source * self.effective_source[0]
        return hlms


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    results_dict = ve.process_dictionary_of_expressions(
        SEOBNRv5_aligned_spin_waveform_quantities().__dict__,
        fixed_mpfs_for_free_symbols=True,
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )
