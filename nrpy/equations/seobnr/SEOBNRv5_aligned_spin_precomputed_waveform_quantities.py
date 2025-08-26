"""
Construct symbolic expression for the SEOBNRv5 aligned-spin gravitational-wave strain and flux when precomputing the waveform coefficients.
The coefficients entering the expressions for the strain and flux are
functions of mass ratio and spin. Since the spins do not evolve when they are aligned,
we can compute these coefficients once at the start of the code to accelerate
the computation of fluxes and strain.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
from typing import Dict

import sympy as sp

from nrpy.equations.grhd.Min_Max_and_Piecewise_Expressions import coord_greater_bound
from nrpy.helpers.float_to_rational import f2r

# The name of this module ("WaveEquation") is given by __name__:
thismodule = __name__


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
        self.max_vh3_order = 3
        self.max_vomega_order = 10
        self.max_l = 8
        self.max_l_f_modes = 5
        self.max_l_delta_modes = 5
        self.rho = sp.MutableDenseNDimArray.zeros(
            self.max_l + 1, self.max_l + 1, self.max_vomega_order + 1, self.max_l + 1
        )
        self.deltalm = sp.MutableDenseNDimArray.zeros(
            self.max_l_delta_modes + 1,
            self.max_l_delta_modes + 1,
            self.max_vh3_order + 1,
            self.max_vomega_order + 1,
        )
        self.fspin = sp.MutableDenseNDimArray.zeros(
            self.max_l_f_modes + 1,
            self.max_l_delta_modes + 1,
            self.max_vomega_order + 1,
        )
        self.fspin_limit = sp.MutableDenseNDimArray.zeros(
            self.max_l_f_modes + 1, self.max_l_f_modes + 1, self.max_vomega_order + 1
        )
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
        self.eulerlog = sp.zeros(9)
        self.eulerlog[0] = 1
        for m in range(1, 9):
            self.eulerlog[m] = EulerGamma + sp.log(2 * m * self.vomega)
        self.khat = sp.zeros(9)
        for m in range(9):
            self.khat[m] = m * self.Omega * self.Hreal

        self.deltalm[5, 5, 1, 0] = (
            sp.Rational(96875, 131250) + sp.Rational(857528, 131250) * self.nu
        ) / (1 - 2 * self.nu)

        self.deltalm[5, 5, 2, 0] = sp.Rational(3865, 429) * sp.pi * (self.vh3) ** 2
        self.deltalm[5, 5, 3, 0] = (
            sp.Rational(-7686949127, 31783752)
            + sp.Rational(954500400, 31783752) * sp.pi**2
        )

        self.deltalm[4, 4, 1, 0] = (
            sp.Rational(112, 120) + sp.Rational(219, 120) * self.nu
        )
        self.deltalm[4, 4, 2, 0] = sp.Rational(25136, 3465)
        self.deltalm[4, 4, 3, 0] = sp.Rational(
            201088, 10395
        ) * sp.pi**2 - sp.Rational(55144, 375)

        self.deltalm[4, 3, 1, 0] = (
            (sp.Rational(4961, 810) * self.nu + sp.Rational(3, 5))
        ) / (1 - 2 * self.nu)
        self.deltalm[4, 3, 2, 0] = sp.Rational(1571, 385) * sp.pi

        self.deltalm[3, 3, 1, 0] = sp.Rational(13, 10)
        self.deltalm[3, 3, 2, 0] = sp.Rational(39, 7)
        self.deltalm[3, 3, 3, 0] = (
            -sp.Rational(227827, 3000) + sp.Rational(78, 7) * sp.pi**2
        )
        self.deltalm[3, 3, 0, 5] = -sp.Rational(80897, 2430) * self.nu

        self.deltalm[3, 2, 1, 0] = (
            (sp.Rational(11, 5) * self.nu + sp.Rational(2, 3))
        ) / (1 - 3 * self.nu)
        self.deltalm[3, 2, 2, 0] = sp.Rational(52, 21) * sp.pi
        self.deltalm[3, 2, 3, 0] = (sp.Rational(208, 63) * sp.pi**2) - sp.Rational(
            9112, 405
        )

        self.deltalm[2, 2, 1, 0] = sp.Rational(7, 3)
        self.deltalm[2, 2, 2, 0] = (
            ((sp.Rational(8, 3) * self.nu - sp.Rational(4, 3)) * self.chi_S)
            - (sp.Rational(4, 3) * self.delta * self.chi_A)
            + (sp.Rational(428, 105) * sp.pi)
        )
        self.deltalm[2, 2, 3, 0] = (sp.Rational(1712, 315) * sp.pi**2) - (
            sp.Rational(2203, 81)
        )
        self.deltalm[2, 2, 0, 5] = -24 * self.nu

        self.deltalm[2, 1, 1, 0] = sp.Rational(2, 3)
        self.deltalm[2, 1, 2, 0] = sp.Rational(107, 105) * sp.pi
        self.deltalm[2, 1, 1, 0] = (sp.Rational(214, 315) * sp.pi**2) - (
            sp.Rational(272, 81)
        )
        self.deltalm[2, 1, 1, 0] = -(sp.Rational(25, 2) * self.nu)

        self.fspin[5, 5, 3] = (
            -70 * self.nu / (3 * (-1 + 2 * self.nu))
            + 110 * self.nu**2 / (3 * (-1 + 2 * self.nu))
            + 10 / (3 * (-1 + 2 * self.nu))
        ) * (self.chi_A * self.deltainv) + (
            10 / (3 * (-1 + 2 * self.nu))
            - 10 * self.nu / (-1 + 2 * self.nu)
            + 10 * self.nu**2 / (-1 + 2 * self.nu)
        ) * self.chi_S
        self.fspin[5, 5, 4] = (
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
        )

        self.fspin_limit[5, 5, 3] = (
            sp.Rational(-70, 3) * self.nu / ((-1 + 2 * self.nu))
            + sp.Rational(110, 3) * self.nu**2 / ((-1 + 2 * self.nu))
            + sp.Rational(10, 3) / ((-1 + 2 * self.nu))
        ) * (self.chi_A)
        self.fspin_limit[5, 5, 4] = (
            (
                -5 / (-1 + 2 * self.nu)
                + 30 * self.nu / (-1 + 2 * self.nu)
                - 40 * self.nu**2 / (-1 + 2 * self.nu)
            )
            * self.chi_A
            * self.chi_S
        )

        self.fspin[4, 1, 1] = (
            sp.Rational(5, 2)
            * self.nu
            / (1 - 2 * self.nu)
            * (self.chi_S - self.chi_A * self.deltainv)
        )

        self.fspin_limit[4, 1, 1] = (
            sp.Rational(5, 2) * self.nu / (1 - 2 * self.nu) * (-self.chi_A)
        )

        self.fspin[4, 3, 1] = (
            1
            / (1 - 2 * self.nu)
            * (
                sp.Rational(5, 2) * self.nu * self.chi_S
                - sp.Rational(5, 2) * self.nu * self.chi_A * self.deltainv
            )
        )
        self.fspin[4, 3, 3] = (
            1
            / (1 - 2 * self.nu)
            * (
                (sp.Rational(887, 44) * self.nu - sp.Rational(3143, 132) * self.nu**2)
                * self.chi_A
                * self.deltainv
                + (-sp.Rational(529, 132) * self.nu**2 - sp.Rational(667, 44) * self.nu)
                * self.chi_S
            )
        )
        self.fspin[4, 3, 4] = (
            1
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

        self.fspin_limit[4, 3, 1] = (
            1 / (1 - 2 * self.nu) * (-sp.Rational(5, 2) * self.nu * self.chi_A)
        )
        self.fspin_limit[4, 3, 3] = (
            1
            / (1 - 2 * self.nu)
            * (
                (sp.Rational(887, 44) * self.nu - sp.Rational(3143, 132) * self.nu**2)
                * self.chi_A
            )
        )
        self.fspin_limit[4, 3, 4] = (
            1
            / (1 - 2 * self.nu)
            * (
                (sp.Rational(137, 6) * self.nu**2 - 18 * self.nu + 3)
                * self.chi_A
                * self.chi_S
            )
        )

        self.fspin[3, 1, 3] = (
            self.chi_A * self.deltainv * (-4 + 11 * self.nu)
            + self.chi_S * (-4 + 13 * self.nu)
        ) / (2)

        self.fspin_limit[3, 1, 3] = -self.chi_A * sp.Rational(5, 8)

        self.fspin[3, 3, 3] = (
            (sp.Rational(19, 2) * self.nu - 2) * self.chi_A * self.deltainv
        ) + ((sp.Rational(5, 2) * self.nu - 2) * self.chi_S)
        self.fspin[3, 3, 4] = (
            (sp.Rational(3, 2) - 6 * self.nu) * self.chi_A**2
            + (3 - 12 * self.nu) * (self.chi_A * self.chi_S * self.deltainv)
            + sp.Rational(3, 2) * self.chi_S**2
        )
        self.fspin[3, 3, 5] = (
            (
                sp.Rational(407, 30) * self.nu**2
                - sp.Rational(593, 60) * self.nu
                + sp.Rational(2, 3)
            )
            * self.chi_A
            * self.deltainv
        ) + (
            (
                sp.Rational(241, 30) * self.nu**2
                + sp.Rational(11, 20) * self.nu
                + sp.Rational(2, 3)
            )
            * self.chi_S
        )
        self.fspin[3, 3, 6] = (
            (-12 * self.nu**2 + sp.Rational(11, 2) * self.nu - sp.Rational(7, 4))
            * self.chi_A**2
            + (44 * self.nu**2 - self.nu - sp.Rational(7, 2))
            * (self.chi_A * self.chi_S * self.deltainv)
            + (6 * self.nu**2 - sp.Rational(27, 2) * self.nu - sp.Rational(7, 4))
            * self.chi_S**2
        )

        self.fspinimag = (self.vh3) ** 2 * (
            sp.Rational(7339, 540) * self.nu - sp.Rational(81, 20)
        ) * self.chi_A * self.deltainv + (
            sp.Rational(593, 108) * self.nu - sp.Rational(81, 20)
        ) * self.chi_S

        self.fspin_limit[3, 3, 3] = (sp.Rational(19, 2) * self.nu - 2) * self.chi_A
        self.fspin_limit[3, 3, 4] = (3 - 12 * self.nu) * (self.chi_A * self.chi_S)
        self.fspin_limit[3, 3, 5] = (
            sp.Rational(407, 30) * self.nu**2
            - sp.Rational(593, 60) * self.nu
            + sp.Rational(2, 3)
        ) * self.chi_A
        self.fspin_limit[3, 3, 6] = (
            44 * self.nu**2 - self.nu - sp.Rational(7, 2)
        ) * (self.chi_A * self.chi_S)

        self.fspinimag_limit = (
            (self.vh3) ** 2
            * (sp.Rational(7339, 540) * self.nu - sp.Rational(81, 20))
            * self.chi_A
        )

        self.fspin[2, 1, 1] = -sp.Rational(3, 2) * (
            self.chi_A * self.deltainv + self.chi_S
        )
        self.fspin[2, 1, 3] = (sp.Rational(131, 84) * self.nu + sp.Rational(61, 12)) * (
            self.chi_A * self.deltainv
        ) + (sp.Rational(79, 84) * self.nu + sp.Rational(61, 12)) * self.chi_S
        self.fspin[2, 1, 4] = (
            (-2 * self.nu - 3) * self.chi_A**2
            + (sp.Rational(21, 2) * self.nu - 6)
            * (self.chi_A * self.chi_S * self.deltainv)
            + (sp.Rational(1, 2) * self.nu - 3) * self.chi_S**2
        )
        self.fspin[2, 1, 5] = (
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
        self.fspin[2, 1, 6] = (
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

        self.fspin_limit[2, 1, 1] = -sp.Rational(3, 2) * (self.chi_A)
        self.fspin_limit[2, 1, 3] = (
            sp.Rational(131, 84) * self.nu + sp.Rational(61, 12)
        ) * (self.chi_A)
        self.fspin_limit[2, 1, 4] = (sp.Rational(21, 2) * self.nu - 6) * (
            self.chi_A * self.chi_S
        )
        self.fspin_limit[2, 1, 5] = (
            (
                sp.Rational(-703, 112) * self.nu**2
                + sp.Rational(8797, 1008) * self.nu
                - sp.Rational(81, 16)
            )
            * (self.chi_A)
            + (sp.Rational(3, 4) - 3 * self.nu) * (self.chi_A**3)
            + (sp.Rational(9, 4) - 6 * self.nu) * (self.chi_A * self.chi_S**2)
        )
        self.fspin_limit[2, 1, 6] = (
            sp.Rational(9487, 504) * self.nu**2
            - sp.Rational(1636, 21) * self.nu
            + sp.Rational(4163, 126)
        ) * (self.chi_A * self.chi_S)

        self.rho[8, 8, 0, 0] = 1
        self.rho[8, 8, 2, 0] = (
            3482
            - 26778 * self.nu
            + 64659 * self.nu**2
            - 53445 * self.nu**3
            + 12243 * self.nu**4
        ) / (2736 * (-1 + 7 * self.nu - 14 * self.nu**2 + 7 * self.nu**3))

        self.rho[8, 7, 0, 0] = 1
        self.rho[8, 7, 2, 0] = (
            23478
            - 154099 * self.nu
            + 309498 * self.nu**2
            - 207550 * self.nu**3
            + 38920 * self.nu**4
        ) / (18240 * (-1 + 6 * self.nu - 10 * self.nu**2 + 4 * self.nu**3))

        self.rho[8, 6, 0, 0] = 1
        self.rho[8, 6, 2, 0] = (
            1002
            - 7498 * self.nu
            + 17269 * self.nu**2
            - 13055 * self.nu**3
            + 2653 * self.nu**4
        ) / (912 * (-1 + 7 * self.nu - 14 * self.nu**2 + 7 * self.nu**3))

        self.rho[8, 5, 0, 0] = 1
        self.rho[8, 5, 2, 0] = (
            4350
            - 28055 * self.nu
            + 54642 * self.nu**2
            - 34598 * self.nu**3
            + 6056 * self.nu**4
        ) / (3648 * (-1 + 6 * self.nu - 10 * self.nu**2 + 4 * self.nu**3))

        self.rho[8, 4, 0, 0] = 1
        self.rho[8, 4, 2, 0] = (
            2666
            - 19434 * self.nu
            + 42627 * self.nu**2
            - 28965 * self.nu**3
            + 4899 * self.nu**4
        ) / (2736 * (-1 + 7 * self.nu - 14 * self.nu**2 + 7 * self.nu**3))

        self.rho[8, 3, 0, 0] = 1
        self.rho[8, 3, 2, 0] = (
            20598
            - 131059 * self.nu
            + 249018 * self.nu**2
            - 149950 * self.nu**3
            + 24520 * self.nu**4
        ) / (18240 * (-1 + 6 * self.nu - 10 * self.nu**2 + 4 * self.nu**3))

        self.rho[8, 2, 0, 0] = 1
        self.rho[8, 2, 2, 0] = (
            2462
            - 17598 * self.nu
            + 37119 * self.nu**2
            - 22845 * self.nu**3
            + 3063 * self.nu**4
        ) / (2736 * (-1 + 7 * self.nu - 14 * self.nu**2 + 7 * self.nu**3))

        self.rho[8, 1, 0, 0] = 1
        self.rho[8, 1, 2, 0] = (
            20022
            - 126451 * self.nu
            + 236922 * self.nu**2
            - 138430 * self.nu**3
            + 21640 * self.nu**4
        ) / (18240 * (-1 + 6 * self.nu - 10 * self.nu**2 + 4 * self.nu**3))

        self.rho[7, 7, 0, 0] = 1
        self.rho[7, 7, 2, 0] = (
            1380 * self.nu**3 - 4963 * self.nu**2 + 4246 * self.nu - 906
        ) / (714 * (3 * self.nu**2 - 4 * self.nu + 1))

        self.rho[7, 6, 0, 0] = 1
        self.rho[7, 6, 2, 0] = (
            6104 * self.nu**4
            - 29351 * self.nu**3
            + 37828 * self.nu**2
            - 16185 * self.nu
            + 2144
        ) / (1666 * (7 * self.nu**3 - 14 * self.nu**2 + 7 * self.nu - 1))

        self.rho[7, 5, 0, 0] = 1
        self.rho[7, 5, 2, 0] = (
            804 * self.nu**3 - 3523 * self.nu**2 + 3382 * self.nu - 762
        ) / (714 * (3 * self.nu**2 - 4 * self.nu + 1))

        self.rho[7, 4, 0, 0] = 1
        self.rho[7, 4, 2, 0] = (
            41076 * self.nu**4
            - 217959 * self.nu**3
            + 298872 * self.nu**2
            - 131805 * self.nu
            + 17756
        ) / (14994 * (7 * self.nu**3 - 14 * self.nu**2 + 7 * self.nu - 1))

        self.rho[7, 3, 0, 0] = 1
        self.rho[7, 3, 2, 0] = (
            420 * self.nu**3 - 2563 * self.nu**2 + 2806 * self.nu - 666
        ) / (714 * (3 * self.nu**2 - 4 * self.nu + 1))

        self.rho[7, 2, 0, 0] = 1
        self.rho[7, 2, 2, 0] = (
            32760 * self.nu**4
            - 190239 * self.nu**3
            + 273924 * self.nu**2
            - 123489 * self.nu
            + 16832
        ) / (14994 * (7 * self.nu**3 - 14 * self.nu**2 + 7 * self.nu - 1))

        self.rho[7, 1, 0, 0] = 1
        self.rho[7, 1, 2, 0] = (
            228 * self.nu**3 - 2083 * self.nu**2 + 2518 * self.nu - 618
        ) / (714 * (3 * self.nu**2 - 4 * self.nu + 1))

        self.rho[6, 6, 0, 0] = 1
        self.rho[6, 6, 2, 0] = (
            273 * self.nu**3 - 861 * self.nu**2 + 602 * self.nu - 106
        ) / (84 * (5 * self.nu**2 - 5 * self.nu + 1))
        self.rho[6, 6, 4, 0] = -sp.Rational(1025435, 659736)

        self.rho[6, 5, 0, 0] = 1
        self.rho[6, 5, 2, 0] = (
            220 * self.nu**3 - 910 * self.nu**2 + 838 * self.nu - 185
        ) / (144 * (3 * self.nu**2 - 4 * self.nu + 1))

        self.rho[6, 4, 0, 0] = 1
        self.rho[6, 4, 2, 0] = (
            133 * self.nu**3 - 581 * self.nu**2 + 462 * self.nu - 86
        ) / (84 * (5 * self.nu**2 - 5 * self.nu + 1))
        self.rho[6, 4, 4, 0] = -sp.Rational(476887, 659736)

        self.rho[6, 3, 0, 0] = 1
        self.rho[6, 3, 2, 0] = (
            156 * self.nu**3 - 750 * self.nu**2 + 742 * self.nu - 169
        ) / (144 * (3 * self.nu**2 - 4 * self.nu + 1))

        self.rho[6, 2, 0, 0] = 1
        self.rho[6, 2, 2, 0] = (
            49 * self.nu**3 - 413 * self.nu**2 + 378 * self.nu - 74
        ) / (84 * (5 * self.nu**2 - 5 * self.nu + 1))
        self.rho[6, 2, 4, 0] = -sp.Rational(817991, 3298680)

        self.rho[6, 1, 0, 0] = 1
        self.rho[6, 1, 2, 0] = (
            124 * self.nu**3 - 670 * self.nu**2 + 694 * self.nu - 161
        ) / (144 * (3 * self.nu**2 - 4 * self.nu + 1))

        self.rho[5, 5, 0, 0] = 1
        self.rho[5, 5, 2, 0] = (
            sp.Rational(487, 390) / (-1 + 2 * self.nu)
            - sp.Rational(649, 195) * self.nu / (-1 + 2 * self.nu)
            + sp.Rational(256, 195) * self.nu**2 / (-1 + 2 * self.nu)
        )
        self.rho[5, 5, 4, 0] = -sp.Rational(3353747, 2129400) - f2r(2.61) * self.nu
        self.rho[5, 5, 6, 0] = (
            sp.Rational(190606537999247, 11957879934000) + f2r(1.25) * self.nu
        )
        self.rho[5, 5, 6, 5] = -sp.Rational(1546, 429)
        self.rho[5, 5, 8, 0] = (
            -sp.Rational(1213641959949291437, 118143853747920000) - f2r(35.7) * self.nu
        )
        self.rho[5, 5, 8, 5] = sp.Rational(376451, 83655)
        self.rho[5, 5, 10, 0] = -sp.Rational(
            150082616449726042201261, 4837990810977324000000
        )
        self.rho[5, 5, 10, 5] = sp.Rational(2592446431, 456756300)

        self.rho[5, 4, 1, 0] = 1
        self.rho[5, 4, 2, 0] = (
            33320 * self.nu**3 - 127610 * self.nu**2 + 96019 * self.nu - 17448
        ) / (13650 * (5 * self.nu**2 - 5 * self.nu + 1))
        self.rho[5, 4, 4, 0] = -sp.Rational(16213384, 15526875)

        self.rho[5, 3, 0, 0] = 1
        self.rho[5, 3, 2, 0] = (176 * self.nu**2 - 850 * self.nu + 375) / (
            390 * (2 * self.nu - 1)
        )
        self.rho[5, 3, 4, 0] = -sp.Rational(410833, 709800)

        self.rho[5, 2, 0, 0] = 1
        self.rho[5, 2, 2, 0] = (
            21980 * self.nu**3 - 104930 * self.nu**2 + 84679 * self.nu - 15828
        ) / (13650 * (5 * self.nu**2 - 5 * self.nu + 1))
        self.rho[5, 2, 4, 0] = -sp.Rational(7187914, 15526875)

        self.rho[5, 1, 0, 0] = 1
        self.rho[5, 1, 2, 0] = (8 * self.nu**2 - 626 * self.nu + 319) / (
            390 * (2 * self.nu - 1)
        )
        self.rho[5, 1, 4, 0] = -sp.Rational(31877, 304200)

        self.rho[4, 1, 0, 0] = 1
        self.rho[4, 1, 2, 0] = (288 * self.nu**2 - 1385 * self.nu + 602) / (
            528 * (2 * self.nu - 1)
        )
        self.rho[4, 1, 4, 0] = -(sp.Rational(7775491, 21141120))
        self.rho[4, 1, 6, 0] = sp.Rational(1227423222031, 1758095539200)
        self.rho[4, 1, 6, 1] = -sp.Rational(1571, 6930)

        self.rho[4, 2, 0, 0] = 1
        self.rho[4, 2, 2, 0] = (1146 - 3530 * self.nu + 285 * self.nu**2) / (
            1320 * (-1 + 3 * self.nu)
        )
        self.rho[4, 2, 3, 0] = (
            self.chi_A * (10 - 21 * self.nu) * self.delta
            + self.chi_S * (10 - 59 * self.nu + 78 * self.nu**2)
        ) / (15 * (-1 + 3 * self.nu))
        self.rho[4, 2, 4, 0] = (
            -114859044
            + 295834536 * self.nu
            + 1204388696 * self.nu**2
            - 3047981160 * self.nu**3
            - 379526805 * self.nu**4
        ) / (317116800 * (-1 + 3 * self.nu) ** 2)
        self.rho[4, 2, 6, 0] = sp.Rational(848238724511, 219761942400)
        self.rho[4, 2, 6, 2] = -sp.Rational(3142, 3465)

        self.rho[4, 3, 0, 0] = 1
        self.rho[4, 3, 2, 0] = (
            -sp.Rational(10, 11) * self.nu**2
            + sp.Rational(547, 176) * self.nu
            - sp.Rational(111, 88)
        ) / (1 - 2 * self.nu)
        self.rho[4, 3, 4, 0] = -sp.Rational(6894273, 7047040) - f2r(0.654) * self.nu
        self.rho[4, 3, 6, 0] = (
            sp.Rational(1664224207351, 195343948800) - f2r(3.69) * self.nu
        )
        self.rho[4, 3, 6, 3] = -sp.Rational(1571, 770)
        self.rho[4, 3, 8, 0] = (
            -sp.Rational(2465107182496333, 460490801971200) + f2r(18.5) * self.nu
        )
        self.rho[4, 3, 8, 3] = sp.Rational(174381, 67760)

        self.rho[4, 4, 1, 0] = 1
        self.rho[4, 4, 2, 0] = (1614 - 5870 * self.nu + 2625 * self.nu**2) / (
            1320 * (-1 + 3 * self.nu)
        )
        self.rho[4, 4, 3, 0] = (
            self.chi_A * (10 - 39 * self.nu) * self.delta
            + self.chi_S * (10 - 41 * self.nu + 42 * self.nu**2)
        ) / (15 * (-1 + 3 * self.nu))
        self.rho[4, 4, 4, 0] = (
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
        self.rho[4, 4, 5, 0] = self.chi_A * self.delta * (
            -8280 + 42716 * self.nu - 57990 * self.nu**2 + 8955 * self.nu**3
        ) / (6600 * (-1 + 3 * self.nu) ** 2) + self.chi_S * (
            -8280
            + 66284 * self.nu
            - 176418 * self.nu**2
            + 128085 * self.nu**3
            + 88650 * self.nu**2 * self.nu**2
        ) / (
            6600 * (-1 + 3 * self.nu) ** 2
        )
        self.rho[4, 4, 6, 0] = (
            sp.Rational(16600939332793, 1098809712000) - f2r(3.56) * self.nu
        )
        self.rho[4, 4, 6, 4] = -sp.Rational(12568, 3465)
        self.rho[4, 4, 8, 0] = (
            -sp.Rational(172066910136202271, 19426955708160000) + f2r(15.6) * self.nu
        )
        self.rho[4, 4, 8, 4] = sp.Rational(845198, 190575)
        self.rho[4, 4, 10, 0] = (
            -sp.Rational(17154485653213713419357, 568432724020761600000) - 216 * self.nu
        )
        self.rho[4, 4, 10, 4] = sp.Rational(22324502267, 3815311500)

        self.rho[3, 1, 0, 0] = 1
        self.rho[3, 1, 2, 0] = -(sp.Rational(2, 9) * self.nu + sp.Rational(13, 18))
        self.rho[3, 1, 4, 0] = (
            -sp.Rational(829, 1782) * self.nu**2
            - sp.Rational(1685, 1782) * self.nu
            + sp.Rational(101, 7128)
        )
        self.rho[3, 1, 6, 0] = sp.Rational(11706720301, 6129723600)
        self.rho[3, 1, 6, 1] = -sp.Rational(26, 63)
        self.rho[3, 1, 8, 1] = sp.Rational(169, 567)
        self.rho[3, 1, 8, 0] = sp.Rational(2606097992581, 4854741091200)

        self.rho[3, 2, 0, 0] = 1
        self.rho[3, 2, 1, 0] = (4 * self.nu * self.chi_S) / (3 * (1 - 3 * self.nu))
        self.rho[3, 2, 2, 0] = (
            sp.Rational(-32, 27) * self.nu**2
            + sp.Rational(223, 54) * self.nu
            - sp.Rational(164, 135)
        ) / (1 - 3 * self.nu) - (16 * self.nu**2 * self.chi_S**2) / (
            9 * (1 - 3 * self.nu) ** 2
        )
        self.rho[3, 2, 3, 0] = (
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
        self.rho[3, 2, 4, 0] = (
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
        self.rho[3, 2, 5, 0] = (
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
        self.rho[3, 2, 6, 0] = sp.Rational(5849948554, 940355325) + f2r(0.333) * self.nu
        self.rho[3, 2, 6, 2] = -sp.Rational(104, 63)
        self.rho[3, 2, 8, 2] = sp.Rational(17056, 8505)
        self.rho[3, 2, 8, 0] = (
            -sp.Rational(10607269449358, 3072140846775) - f2r(6.5) * self.nu
        )
        self.rho[3, 2, 10, 0] = (
            -sp.Rational(1312549797426453052, 176264081083715625) + 98 * self.nu
        )
        self.rho[3, 2, 10, 2] = sp.Rational(18778864, 12629925)

        self.rho[3, 3, 0, 0] = 1
        self.rho[3, 3, 2, 0] = sp.Rational(2, 3) * self.nu - sp.Rational(7, 6)
        self.rho[3, 3, 4, 0] = (
            -sp.Rational(6719, 3960)
            - sp.Rational(1861, 990) * self.nu
            + sp.Rational(149, 330) * self.nu**2
        )
        self.rho[3, 3, 6, 0] = (
            sp.Rational(3203101567, 227026800)
            + (-sp.Rational(129509, 25740) + sp.Rational(41, 192) * sp.pi**2) * self.nu
            - sp.Rational(274621, 154440) * self.nu**2
            + sp.Rational(12011, 46332) * self.nu**3
        )
        self.rho[3, 3, 6, 3] = -sp.Rational(26, 7)
        self.rho[3, 3, 8, 0] = -sp.Rational(57566572157, 8562153600) + self.nu * 12
        self.rho[3, 3, 8, 3] = sp.Rational(13, 3)
        self.rho[3, 3, 10, 0] = (
            -sp.Rational(903823148417327, 30566888352000) - 215 * self.nu
        )
        self.rho[3, 3, 10, 3] = sp.Rational(87347, 13860)

        self.rho[2, 2, 0, 0] = 1
        self.rho[2, 2, 2, 0] = sp.Rational(55, 84) * self.nu - sp.Rational(43, 42)
        self.rho[2, 2, 3, 0] = (
            -2 * (self.chi_S + self.chi_A * self.delta - self.chi_S * self.nu)
        ) / 3
        self.rho[2, 2, 4, 0] = (
            sp.Rational(19583, 42336) * self.nu**2
            - sp.Rational(33025, 21168) * self.nu
            - sp.Rational(20555, 10584)
            + (sp.Rational(1, 2) - 2 * self.nu) * self.chi_A**2
            + self.delta * self.chi_A * self.chi_S
            + sp.Rational(1, 2) * self.chi_S**2
        )
        self.rho[2, 2, 5, 0] = (
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
        self.rho[2, 2, 6, 2] = -sp.Rational(428, 105)
        self.rho[2, 2, 6, 0] = (
            sp.Rational(10620745, 39118464) * self.nu**3
            - sp.Rational(6292061, 3259872) * self.nu**2
            + sp.Rational(41, 192) * sp.pi**2 * self.nu
            - sp.Rational(48993925, 9779616) * self.nu
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
        self.rho[2, 2, 7, 0] = (
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
        self.rho[2, 2, 8, 2] = sp.Rational(9202, 2205)
        self.rho[2, 2, 8, 0] = -sp.Rational(387216563023, 160190110080) + self.nu * f2r(
            21.2
        )
        self.rho[2, 2, 10, 2] = sp.Rational(439877, 55566)
        self.rho[2, 2, 10, 0] = (
            sp.Rational(16094530514677, 533967033600) - 411 * self.nu
        )

        self.rho[2, 1, 0, 0] = 1
        self.rho[2, 1, 2, 0] = sp.Rational(23, 84) * self.nu - sp.Rational(59, 56)
        self.rho[2, 1, 4, 0] = (
            sp.Rational(617, 4704) * self.nu**2
            - sp.Rational(10993, 14112) * self.nu
            - sp.Rational(47009, 56448)
        )
        self.rho[2, 1, 6, 0] = sp.Rational(7613184941, 2607897600) + f2r(1.65) * self.nu
        self.rho[2, 1, 6, 1] = -sp.Rational(107, 105)
        self.rho[2, 1, 8, 0] = (
            -sp.Rational(1168617463883, 911303737344) + f2r(26.5) * self.nu
        )
        self.rho[2, 1, 8, 1] = sp.Rational(6313, 5880)
        self.rho[2, 1, 10, 0] = (
            -sp.Rational(63735873771463, 16569158860800) + 80 * self.nu
        )
        self.rho[2, 1, 10, 1] = sp.Rational(5029963, 5927040)

        self.Y = [
            [f2r(0), f2r(0), f2r(0), f2r(0), f2r(0), f2r(0), f2r(0), f2r(0), f2r(0)],
            [
                f2r(0),
                f2r(0.3454941494713355),
                f2r(0),
                f2r(0.3231801841141506),
                f2r(0),
                f2r(0.32028164857621516),
                f2r(0),
                f2r(0.31937046138540076),
                f2r(0),
            ],
            [
                f2r(0),
                f2r(0),
                f2r(0.3862742020231896),
                f2r(0),
                f2r(0.33452327177864466),
                f2r(0),
                f2r(0.32569524293385776),
                f2r(0),
                f2r(0.32254835519288305),
            ],
            [
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0.4172238236327842),
                f2r(0),
                f2r(0.34594371914684025),
                f2r(0),
                f2r(0.331899519333737),
                f2r(0),
            ],
            [
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0.4425326924449826),
                f2r(0),
                f2r(0.3567812628539981),
                f2r(0),
                f2r(0.3382915688890245),
            ],
            [
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0.46413220344085826),
                f2r(0),
                f2r(0.3669287245764378),
                f2r(0),
            ],
            [
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0.48308411358006625),
                f2r(0),
                f2r(0.3764161087284946),
            ],
            [
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0.5000395635705508),
                f2r(0),
            ],
            [
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0),
                f2r(0.5154289843972844),
            ],
        ]
        self.c = sp.zeros(10)
        for k in range(2, 10):
            self.c[k] = (m1 / M) ** (k - 1) + ((-1) ** k) * ((m2 / M) ** (k - 1))
        self.c[3] *= self.noneqcond
        self.c[3] += -1 * self.eqcond
        self.c[5] *= self.noneqcond
        self.c[5] += sp.Rational(-1, 2) * self.eqcond
        self.r0 = 2 / sp.sqrt(sp.E)
        self.n_complex = sp.zeros(9, 9)
        for l in range(2, 9):
            for m in range(1, l + 1):
                if not (l + m) % 2:
                    self.n_complex[l, m] = (
                        8 * sp.pi * (sp.I * m) ** l / sp.factorial2(2 * l + 1)
                    ) * sp.sqrt(((l + 1) * (l + 2)) / ((l) * (l - 1)))
                else:
                    self.n_complex[l, m] = (
                        -16 * sp.I * sp.pi * (sp.I * m) ** l / sp.factorial2(2 * l + 1)
                    ) * sp.sqrt(
                        ((2 * l + 1) * (l + 2) * (l**2 - m**2))
                        / ((2 * l - 1) * (l + 1) * (l) * (l - 1))
                    )
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
        # self.max_vomega_order = 10
        # self.max_l = 8
        # self.max_l_f_modes = 5

        factorized_flux = sp.sympify(0)
        for l in range(2, 9):
            for m in range(1, l + 1):
                n_abs = sp.Symbol(f"commondata->n_abs_{l}_{m}", real=True)
                pn_contribution_f = 0
                rholm = 0
                for j in range(self.max_vomega_order + 1):
                    for k in range(m + 1):
                        if sp.sympify(self.rho[l, m, j, k]).is_nonzero is not False:
                            rholm += (
                                sp.Symbol(f"commondata->rho_{l}_{m}_{j}_{k}")
                                * self.vomega**j
                                * self.eulerlog[k]
                            )
                if not m % 2:
                    pn_contribution_f += rholm**l
                if m % 2:
                    pn_contribution_f += rholm**l
                    fspinlm = 0
                    if (l < 5) or (l == 5 and m == 5):
                        for j in range(self.max_vomega_order + 1):
                            if sp.sympify(self.fspin[l, m, j]).is_nonzero is not False:
                                fspinlm += (
                                    sp.Symbol(
                                        f"commondata->fspin_{l}_{m}_{j}", real=True
                                    )
                                    * self.vomega**j
                                )
                        pn_contribution_f = (
                            self.noneqcond * (fspinlm + rholm**l)
                            + self.eqcond * fspinlm
                        )

                tail_term = self.Tlm[l, m]
                newtonian_strain_amplitude = (
                    self.nu
                    * n_abs
                    * sp.Symbol(f"commondata->c{l + (l + m) % 2}", real=True)
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
        # self.max_vomega_order = 10

        l, m = 2, 2
        rholm = 0
        deltalm = 0
        for j in range(self.max_vomega_order + 1):
            for k in range(m + 1):
                if sp.sympify(self.rho[l, m, j, k]).is_nonzero is not False:
                    rholm += (
                        sp.Symbol(f"commondata->rho_{l}_{m}_{j}_{k}", real=True)
                        * self.vomega**j
                        * self.eulerlog[k]
                    )
            for k in range(self.max_vh3_order + 1):
                if sp.sympify(self.deltalm[l, m, k, j]).is_nonzero is not False:
                    deltalm += (
                        sp.Symbol(f"commondata->deltalm_{l}_{m}_{k}_{j}", real=True)
                        * self.vh3**k
                        * self.vomega**j
                    )
        pn_contribution_f = rholm**l
        pn_contribution_delta = sp.exp(sp.I * deltalm)
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
        n = sp.Symbol(f"commondata->n_complex_{l}_{m}", real=True)
        newtonian_strain = (
            self.nu
            * n
            * sp.Symbol(f"commondata->c{l + (l + m) % 2}", real=True)
            * self.vphi**2
            * self.Y[m][l - (l + m) % 2]
            * sp.exp(-m * sp.I * self.phi)
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
