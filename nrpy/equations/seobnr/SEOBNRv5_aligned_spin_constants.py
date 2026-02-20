"""
Construct the SEOBNRv5 NR-derived constants.

Authors: Siddharth Mahesh
sm0193 at mix dot wvu dot edu
Zachariah B. Etienne
zachetie at gmail *dot com

The SEOBNRv5 model relies on NR-calibrated or inferred constants
for improving the accuracy of the EOB Hamiltonian and merger-ringdown waveform.
The Backwards-One Body (BOB) formalism relies on the properties of the remnant black hole
for the merger-ringdown waveform and NQC corrections.
The NR constants are expressed in terms of the binary masses (m1, m2), and spins (chi1, chi2).

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
from typing import cast

import sympy as sp

from nrpy.equations.grhd.Min_Max_and_Piecewise_Expressions import (
    coord_greater_bound,
    coord_less_bound,
)
from nrpy.helpers.float_to_rational import f2r

# The name of this module is given by __name__:
thismodule = __name__


class SEOBNR_aligned_spin_constants:
    """Class for computing the SEOBNR aligned-spin constants."""

    def __init__(
        self,
        calibration_no_spin: bool = False,
        calibration_spin: bool = False,
    ) -> None:
        """
        Compute the SEOBNR aligned-spin constants.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the SEOBNR aligned-spin constants. It
        initializes class variables like mass parameters, spin parameters, and
        computes the NR calibration as well as the remnant properties of the black hole.
        The key outputs of the SEOBNRv5_aligned_spin_constants class are:
            - 'a_6' : the pseudo-5PN non-spinning Hamiltonian
                            calibration parameter for the SEOBNRv5 model.
                            Equation 78 of https://arxiv.org/pdf/2303.18039.
                            This is only computed if calibration_no_spin is False.
            - 'Delta_t' : the time delay between the peak of the (l=2,m=2) mode
                            and the time when the EOB perturber crosses
                            the innermost stable circular orbit (ISCO) of the remnant.
                            Always defined as Delta_t_NS + Delta_t_S.
                            In production mode (calibration_no_spin and calibration_spin are False),
                            it is evaluated from Equations 79 & 80 of https://arxiv.org/pdf/2303.18039.
                            In calibration modes either Delta_t_NS or Delta_t_S is
                            treated as an external calibration input.
            - 'd_SO' : the spin-orbit calibration parameter for the SEOBNRv5 model.
                            Equation 81 of https://arxiv.org/pdf/2303.18039.
                            This is only computed if calibration_spin is False.
            - 'a_f' : the final spin of the remnant black hole
                            using https://arxiv.org/pdf/1605.01938 and implemented in
                            https://lscsoft.docs.ligo.org/lalsuite/lalinference/nrutils_8py_source.html#l01020
                            with inputs "M3J4" for M=3 and J=4.
            - 'M_f' : the final mass of the remnant black hole
                            using https://arxiv.org/abs/1611.00332 and implemented in
                            bbh_final_mass_non_precessing_UIB2016() in
                            https://arxiv.org/src/1611.00332v2/anc/FinalStateUIB2016.py
                            with version="v2" as input.
            - 'rISCO' : the ISCO radius of the remnant black hole.
                            using Equation 2.21 of Bardeen, Press, and Teukolsky,
                            https://ui.adsabs.harvard.edu/abs/1972ApJ...178..347B/abstract
            - 'rstop' : the radius of the remnant black hole at the time of ISCO crossing.

        :param calibration_no_spin: Flag to enable/disable calibration of the non-spinning parameters
        :param calibration_spin: Flag to enable/disable calibration of the spinning parameters
        :raises ValueError: If both calibration_no_spin and calibration_spin are True
        :return None:
        """
        # The calibration process for the SEOBNRv5 is done in two steps:
        # 1. Calibration of the non-spinning coefficients
        # 2. Calibration of the spin-dependent coefficients
        # Therefore, the C code can only be generated for one of the above calibration options.
        if calibration_no_spin and calibration_spin:
            raise ValueError(
                "calibration_no_spin and calibration_spin cannot both be True."
            )
        self.m1, self.m2, self.chi1, self.chi2 = sp.symbols(
            "m1 m2 chi1 chi2", real=True
        )
        # compute calibration parameters
        if calibration_no_spin:
            # This is the first (non-spinning) calibration stage so we have no precalculated values
            self.a6, self.Delta_t_NS = sp.symbols("a6 Delta_t_NS", real=True)
            # In non-spinning calibration, spin-dependent calibration terms are disabled:
            # we set dSO and Delta_t_S to 0. The calibration workflow should supply chi1=chi2=0.
            self.dSO = sp.sympify(0)
            self.Delta_t_S = sp.sympify(0)
        elif calibration_spin:
            # This is the second (spinning) calibration stage where we have precalculated values for a6 and Delta_t_NS
            self.compute_calibration_params()
            # overwrite Delta_t_S and dSO to symbols
            self.dSO, self.Delta_t_S = sp.symbols("dSO Delta_t_S", real=True)
        else:
            # This is the post-calibration stage and all values are precalculated
            self.compute_calibration_params()
        # For all stages Delta_t is defined as the sum of non-spinning and spinning parameters
        self.Delta_t = self.Delta_t_NS + self.Delta_t_S
        # Final mass and spin computation
        self.final_spin_non_precessing_HBR2016()
        self.final_mass_non_precessing_UIB2016()
        self.rISCO = self.Kerr_ISCO_radius(self.a_f)
        self.rstop = -1 * coord_less_bound(self.Delta_t, sp.sympify(0)).subs(
            sp.Function("nrpyAbs"), sp.Abs
        ) + f2r(0.98) * self.rISCO * coord_greater_bound(
            self.Delta_t, sp.sympify(0)
        ).subs(
            sp.Function("nrpyAbs"), sp.Abs
        )

    def Kerr_ISCO_radius(self, a: sp.Expr) -> sp.Expr:
        """
        Compute the radius of the innermost stable circular orbit (ISCO) of a Kerr black hole.

        :param a: dimensionless spin parameter of the Kerr black hole.
        :return: radius of the ISCO in Boyer-Lindquist coordinates.
        """
        # restrict a to 1.
        a_ceil_one = (a + 1 - sp.Abs(a - 1)) / 2
        z1 = 1 + (1 - a_ceil_one**2) ** (sp.Rational(1, 3)) * (
            (1 + a_ceil_one) ** (sp.Rational(1, 3))
            + (1 - a_ceil_one) ** (sp.Rational(1, 3))
        )
        z2 = sp.sqrt(3 * a_ceil_one**2 + z1**2)
        a_sign = sp.sign(a_ceil_one)
        return cast(sp.Expr, 3 + z2 - sp.sqrt((3 - z1) * (3 + z1 + 2 * z2)) * a_sign)

    def compute_calibration_params(
        self,
    ) -> None:
        """
        Compute the calibration parameters for the SEOBNRv5 aligned-spin model.

        :return None:
        """
        m1 = self.m1
        m2 = self.m2
        chi1 = self.chi1
        chi2 = self.chi2
        M = m1 + m2
        nu = m1 * m2 / M**2
        ap = (m1 * chi1 + m2 * chi2) / M**2
        am = (m1 * chi1 - m2 * chi2) / M**2
        par_a6 = [
            f2r(4.17877875e01),
            f2r(-3.02193382e03),
            f2r(3.34144394e04),
            f2r(-1.69019140e05),
            f2r(3.29523262e05),
        ]
        self.pyseobnr_a6 = (
            par_a6[0]
            + par_a6[1] * nu
            + par_a6[2] * nu**2
            + par_a6[3] * nu**3
            + par_a6[4] * nu**4
        )
        self.pyseobnr_dSO = (
            -f2r(7.71251231383957) * am**3
            - f2r(17.2294679794015) * am**2 * ap
            - f2r(238.430383378296) * am**2 * nu
            + f2r(69.5461667822545) * am**2
            - f2r(10.5225438990315) * am * ap**2
            + f2r(362.767393298729) * am * ap * nu
            - f2r(85.8036338010274) * am * ap
            - f2r(1254.66845939312) * am * nu**2
            + f2r(472.431937787377) * am * nu
            - f2r(39.742317057316) * am
            - f2r(7.58458103577458) * ap**3
            - f2r(42.7601129678844) * ap**2 * nu
            + f2r(18.1783435552183) * ap**2
            - f2r(201.905934468847) * ap * nu**2
            - f2r(90.5790079104259) * ap * nu
            + f2r(49.6299175121658) * ap
            + f2r(478.546231305475) * nu**3
            + f2r(679.521769948995) * nu**2
            - f2r(177.334831768076) * nu
            - f2r(37.6897780220529)
        )
        self.Delta_t_S = nu ** (sp.Rational(-1, 5) + 0 * nu) * (
            f2r(8.39238879807543) * am**2 * ap
            - f2r(16.9056858928167) * am**2 * nu
            + f2r(7.23410583477034) * am**2
            + f2r(6.38975598319936) * am * ap**2
            + f2r(179.569824846781) * am * ap * nu
            - f2r(40.6063653476775) * am * ap
            + f2r(144.253395844761) * am * nu**2
            - f2r(90.1929138487509) * am * nu
            + f2r(14.2203101910927) * am
            - f2r(6.78913884987037) * ap**4
            + f2r(5.39962303470497) * ap**3
            - f2r(132.224950777226) * ap**2 * nu
            + f2r(49.8016443361381) * ap**2
            + f2r(384.201018794943) * ap * nu**2
            - f2r(141.253181790353) * ap * nu
            + f2r(17.5710132409988) * ap
        )
        par_dtns = [
            f2r(1.00513217e01),
            -f2r(5.96231800e01),
            -f2r(1.05687385e03),
            -f2r(9.79317619e03),
            f2r(5.55652392e04),
        ]
        self.Delta_t_NS = nu ** (sp.Rational(-1, 5) + par_dtns[0] * nu) * (
            par_dtns[1] + par_dtns[2] * nu + par_dtns[3] * nu**2 + par_dtns[4] * nu**3
        )

    def final_spin_non_precessing_HBR2016(
        self,
    ) -> None:
        """
        Compute the final spin for the SEOBNRv5 aligned-spin model.
        The spin is calculated using the non-precessing HBR2016 fits with version "M3J4"
        as outlined in https://lscsoft.docs.ligo.org/lalsuite/lalinference/nrutils_8py_source.html

        :return None:
        """
        m1 = self.m1
        m2 = self.m2
        chi1 = self.chi1
        chi2 = self.chi2
        nM = 3
        nJ = 4
        k = sp.zeros(4, 5)
        k[0, 0] = -f2r(5.97723)
        k[0, 1] = f2r(3.39221)
        k[0, 2] = f2r(4.48865)
        k[0, 3] = -f2r(5.77101)
        k[0, 4] = -f2r(13.0459)
        k[1, 0] = f2r(35.1278)
        k[1, 1] = -f2r(72.9336)
        k[1, 2] = -f2r(86.0036)
        k[1, 3] = f2r(93.7371)
        k[1, 4] = f2r(200.975)
        k[2, 0] = -f2r(146.822)
        k[2, 1] = f2r(387.184)
        k[2, 2] = f2r(447.009)
        k[2, 3] = -f2r(467.383)
        k[2, 4] = -f2r(884.339)
        k[3, 0] = f2r(223.911)
        k[3, 1] = -f2r(648.502)
        k[3, 2] = -f2r(697.177)
        k[3, 3] = f2r(753.738)
        k[3, 4] = f2r(1166.89)
        xi = f2r(0.474046)
        q = m2 / m1
        chi1z = chi1
        chi2z = chi2
        nu = m1 * m2 / (m1 + m2) ** 2
        atot = (chi1z + chi2z * q * q) / ((1 + q) * (1 + q))
        aeff = atot + xi * nu * (chi1z + chi2z)
        rISCOeff = self.Kerr_ISCO_radius(aeff)
        LISCOeff = (sp.Rational(2, 3) / sp.sqrt(3)) * (
            1 + 2 * sp.sqrt(3 * rISCOeff - 2)
        )
        EISCOeff = sp.sqrt(1 - 2 / (3 * rISCOeff))
        aeff_j = [1, aeff, aeff**2, aeff**3, aeff**4]
        nu_i = [1, nu, nu**2, nu**3]
        ksum = 0
        for i in range(nM + 1):
            for j in range(nJ + 1):
                ksum += k[i, j] * nu_i[i] * aeff_j[j]

        ell = sp.Abs(LISCOeff - 2 * atot * (EISCOeff - 1) + nu * ksum)
        self.a_f = atot + ell / (1 / q + 2 + q)

    def final_mass_non_precessing_UIB2016(
        self,
    ) -> None:
        """
        Compute the final mass for the SEOBNRv5 aligned-spin model.
        The mass is calculated using the non-precessing UIB2016 fits with version "v2"
        as outlined in https://lscsoft.docs.ligo.org/lalsuite/lalinference/nrutils_8py_source.html

        :return None:
        """
        m1 = self.m1
        m2 = self.m2
        chi1 = self.chi1
        chi2 = self.chi2
        Shat = (m1 * m1 * chi1 + m2 * m2 * chi2) / (m1 * m1 + m2 * m2)
        nu = m1 * m2 / (m1 + m2) ** 2
        Shat2 = Shat * Shat
        Shat3 = Shat2 * Shat
        nu2 = nu * nu
        nu3 = nu2 * nu
        nu4 = nu3 * nu
        sqrt1m4nu = sp.sqrt(1 - 4 * nu)
        Deltachi = chi1 - chi2
        Deltachi2 = Deltachi * Deltachi
        a2 = f2r(0.5609904135313374)
        a3 = -f2r(0.84667563764404)
        a4 = f2r(3.145145224278187)
        Erad_nu_0 = (
            a4 * nu4 + a3 * nu3 + a2 * nu2 + (1 - sp.Rational(2, 3) * sp.sqrt(2)) * nu
        )
        b1 = -f2r(0.2091189048177395)
        b2 = -f2r(0.19709136361080587)
        b3 = -f2r(0.1588185739358418)
        b5 = f2r(2.9852925538232014)
        f20 = f2r(4.271313308472851)
        f30 = f2r(31.08987570280556)
        f50 = f2r(1.5673498395263061)
        f10 = f2r(1.8083565298668276)
        f21 = 0
        d10 = -f2r(0.09803730445895877)
        d11 = -f2r(3.2283713377939134)
        d20 = f2r(0.01118530335431078)
        d30 = -f2r(0.01978238971523653)
        d31 = -f2r(4.91667749015812)
        f11 = f2r(15.738082204419655)
        f31 = -f2r(243.6299258830685)
        f51 = -f2r(0.5808669012986468)
        bfin1 = b1 * (f10 + f11 * nu + (16 - 16 * f10 - 4 * f11) * nu2)
        bfin2 = b2 * (f20 + f21 * nu + (16 - 16 * f20 - 4 * f21) * nu2)
        bfin3 = b3 * (f30 + f31 * nu + (16 - 16 * f30 - 4 * f31) * nu2)
        bfin5 = b5 * (f50 + f51 * nu + (16 - 16 * f50 - 4 * f51) * nu2)
        Erad_eq_Shat = (
            f2r(0.128) * bfin3 * Shat3
            + f2r(0.211) * bfin2 * Shat2
            + f2r(0.346) * bfin1 * Shat
            + 1
        ) / (1 - f2r(0.212) * bfin5 * Shat)
        Erad_nu_Shat = Erad_nu_0 * Erad_eq_Shat
        d10 = -f2r(0.09803730445895877)
        d11 = -f2r(3.2283713377939134)
        d20 = f2r(0.01118530335431078)
        d30 = -f2r(0.01978238971523653)
        d31 = -f2r(4.91667749015812)
        A_1 = d10 * sqrt1m4nu * nu2 * (d11 * nu + 1)
        A_2 = d20 * nu3
        A_3 = d30 * sqrt1m4nu * nu * (d31 * nu + 1)
        DeltaErad_nu_Shat_Deltachi = (
            A_1 * Deltachi + A_2 * Deltachi2 + A_3 * Deltachi * Shat
        )
        self.M_f = 1 - (Erad_nu_Shat + DeltaErad_nu_Shat_Deltachi)


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
        SEOBNR_aligned_spin_constants().__dict__,
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
