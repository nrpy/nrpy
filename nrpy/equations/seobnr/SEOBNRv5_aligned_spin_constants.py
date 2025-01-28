"""
Construct symbolic expression for the SEOBNRv5 aligned-spin constants.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
from typing import Any, List

import sympy as sp

from nrpy.equations.grhd.Min_Max_and_Piecewise_Expressions import (
    coord_greater_bound,
    coord_less_bound,
)

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

    return [z1[0] * z2[0] - z1[1] * z2[1], z1[0] * z2[1] + z1[1] * z2[0]]


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


class SEOBNR_aligned_spin_constants:
    """Class for computing the BOB aligned-spin gravitational-wave strain and NQC corrections."""

    def __init__(self) -> None:
        """
        Compute the BOB aligned-spin waveform.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the aligned-spin BOB strain. It
        initializes class variables like mass parameters, spin parameters, and
        various coefficients required for the waveforms's amplitude and phase.
        :return None:
        """
        (m1, m2, chi1, chi2) = sp.symbols("m1 m2 chi1 chi2", real=True)

        # Delta_t computation
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
        Delta_t_S = nu ** (sp.Rational(-1, 5) + 0 * nu) * (
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
        Delta_t_NS = nu ** (sp.Rational(-1, 5) + par_dtns[0] * nu) * (
            par_dtns[1] + par_dtns[2] * nu + par_dtns[3] * nu**2 + par_dtns[4] * nu**3
        )
        self.Delta_t = Delta_t_NS + Delta_t_S

        # Final mass and spin computation
        atot = m1 * m1 * chi1 + m2 * m2 * chi2
        aeff = atot + f2r(0.474046) * nu * (chi1 + chi2)
        aeff_max1 = (aeff + 1 - sp.Abs(aeff - 1)) / 2
        Z1eff = 1 + (sp.cbrt(1 - aeff_max1 * aeff_max1)) * (
            sp.cbrt(1 + aeff_max1) + sp.cbrt(1 - aeff_max1)
        )
        Z2eff = sp.sqrt(3 * aeff_max1 * aeff_max1 + Z1eff * Z1eff)
        rISCOeff = (
            3 + Z2eff - sp.sign(aeff_max1) * sp.sqrt((3 - Z1eff) * (3 + Z1eff + Z2eff))
        )
        LISCOeff = (sp.Rational(2, 3) / sp.sqrt(3)) * (
            1 + 2 * sp.sqrt(3 * rISCOeff - 2)
        )
        EISCOeff = sp.sqrt(1 - 2 / (3 * rISCOeff))
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

        NRfactor = 0
        nu_i = 1
        for i in range(4):
            aeff_j = 1
            for j in range(5):
                NRfactor += k[i, j] * nu_i * aeff_j
                aeff_j *= aeff
            nu_i *= nu
        ell = sp.Abs(LISCOeff - 2 * atot * (EISCOeff - 1) + nu * NRfactor)
        self.a_f = atot + nu * ell
        Z1f = 1 + (sp.cbrt(1 - self.a_f * self.a_f)) * (
            sp.cbrt(1 + self.a_f) + sp.cbrt(1 - self.a_f)
        )
        Z2f = sp.sqrt(3 * self.a_f * self.a_f + Z1f * Z1f)
        self.rISCO = 3 + Z2f - sp.sign(self.a_f) * sp.sqrt((3 - Z1f) * (3 + Z1f + Z2f))
        self.rstop = -1 * coord_less_bound(self.Delta_t, 0).subs(
            sp.Function("nrpyAbs"), sp.Abs
        ) + f2r(0.98) * self.rISCO * coord_greater_bound(self.Delta_t, 0).subs(
            sp.Function("nrpyAbs"), sp.Abs
        )
        Shat = (m1 * m1 * chi1 + m2 * m2 * chi2) / (m1 * m1 + m2 * m2)
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
