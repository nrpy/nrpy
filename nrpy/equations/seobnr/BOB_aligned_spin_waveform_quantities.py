"""
Construct symbolic expression for the BOB aligned-spin gravitational-wave strain and NQC corrections.

Author: Siddharth Mahesh

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


class BOB_aligned_spin_waveform_quantities:
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
        (m1, m2, chi1, chi2, omega_qnm, tau_qnm, t_0, t) = sp.symbols(
            "m1 m2 chi1 chi2 omega_qnm tau_qnm t_0 t", real=True
        )

        # Delta_t computation
        M = m1 + m2
        nu = m1 * m2 / M**2
        ap = (m1 * chi1 + m2 * chi2) / M**2
        am = (m2 * chi1 - m2 * chi2) / M**2
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
        par = [
            f2r(1.00513217e01),
            -f2r(5.96231800e01),
            -f2r(1.05687385e03),
            -f2r(9.79317619e03),
            f2r(5.55652392e04),
        ]
        Delta_t_NS = nu ** (sp.Rational(-1, 5) + par[0] * nu) * (
            par[1] + par[2] * nu + par[3] * nu**2 + par[4] * nu**3
        )
        self.Delta_t = Delta_t_NS + Delta_t_S

        # Final mass and spin computation
        atot = m1 * m1 * chi1 + m2 * m2 * chi2
        aeff = atot + f2r(0.474046) * nu * (chi1 + chi2)
        aeff_max1 = aeff + 1 - sp.Abs(aeff - 1)
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
        Shat4 = Shat3 * Shat
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
        h22NR = nu * sp.Abs(
            f2r(71.97969776036882194603) * nu4
            - f2r(13.35761402231352157344) * nu3 * Shat
            - f2r(46.87585958426210908101) * nu3
            + f2r(0.61988944517825661507) * nu2 * Shat2
            + f2r(7.19426416189229733789) * nu2 * Shat
            + f2r(12.44040490932310127903) * nu2
            + f2r(0.43014673069078152023) * nu * Shat3
            - f2r(1.74313546783413597652) * nu * Shat
            - f2r(0.86828935763242798274) * nu
            - f2r(0.08493901280736430859) * Shat3
            - f2r(0.02082621429567295401) * Shat2
            + f2r(0.18693991146784910695) * Shat
            + f2r(1.46709663479911811557)
        )
        omega22NR = -(
            f2r(5.89352329617707670906) * nu4
            + f2r(3.75145580491965446868) * nu3 * Shat
            - f2r(3.34930536209472551334) * nu3
            - f2r(0.97140932625194231775) * nu2 * Shat2
            - f2r(1.69734302394369973577) * nu2 * Shat
            + f2r(0.28539204856044564362) * nu2
            + f2r(0.2419483723662931296) * nu * Shat3
            + f2r(0.51801427018052081941) * nu * Shat2
            + f2r(0.25096450064948544467) * nu * Shat
            - f2r(0.31709602351033533418) * nu
            - f2r(0.01525897668158244028) * Shat4
            - f2r(0.06692658483513916345) * Shat3
            - f2r(0.08715176045684569495) * Shat2
            - f2r(0.09133931944098934441) * Shat
            - f2r(0.2685414392185025978)
        )
        Omega_0 = omega22NR / 2
        Omega_qnm = omega_qnm / 2
        t_p = t_0 - 2 * tau_qnm * sp.log(Omega_0 / Omega_qnm)
        Ap = h22NR * (omega22NR**2) * sp.cosh((t_0 - t_p) / tau_qnm)
        k = (Omega_qnm**4 - Omega_0**4) / (1 - sp.tanh((t_0 - t_p) / tau_qnm))
        kappa_m = Omega_0 * Omega_0 / Omega_qnm
        kappa_p = Omega_qnm
        Omega = (
            Omega_0**4
            + k * (sp.tanh((t - t_p) / tau_qnm) - sp.tanh((t_0 - t_p) / tau_qnm))
        ) ** (1 / 4)
        Omega_0_over_kappa_p = Omega_0 / kappa_p
        Omega_0_over_kappa_m = Omega_0 / kappa_m
        Omega_over_kappa_p = Omega / kappa_p
        Omega_over_kappa_m = Omega / kappa_m
        arctanh_m = (
            0.5
            * kappa_m
            * tau_qnm
            * sp.log(
                (1 + Omega_over_kappa_m)
                * (1 - Omega_0_over_kappa_m)
                / ((1 - Omega_over_kappa_m) * (1 + Omega_0_over_kappa_m))
            )
        )
        arctanh_p = (
            0.5
            * kappa_p
            * tau_qnm
            * sp.log(
                (1 + Omega_over_kappa_p)
                * (1 - Omega_0_over_kappa_p)
                / ((1 - Omega_over_kappa_p) * (1 + Omega_0_over_kappa_p))
            )
        )
        arctan_m = (
            kappa_m * tau_qnm * (sp.atan2(Omega, kappa_m) - sp.atan2(Omega_0, kappa_m))
        )
        arctan_p = (
            kappa_p * tau_qnm * (sp.atan2(Omega, kappa_p) - sp.atan2(Omega_0, kappa_p))
        )
        h = (Ap / 4 / (Omega**2)) * (1 / sp.cosh((t - t_p) / tau_qnm))
        Phi = arctan_p + arctanh_p - arctan_m - arctanh_m
        phi = 2 * Phi
        self.hplus = h * sp.cos(phi)
        self.hcross = h * sp.cos(-phi)
        # mostly trivial
        t_attach = t_0
        self.h_t_attach = h22NR
        self.hdot_t_attach = sp.sympify(0)
        self.hddot_t_attach = sp.diff(sp.diff(h, t), t).subs(t, t_attach)
        self.w_t_attach = omega22NR
        tanht0 = (Omega_0**4 - Omega_qnm**4) / (Omega_0**4 + Omega_qnm**4)
        self.wdot_t_attach = -omega22NR * tanht0 / (2 * tau_qnm)


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
        BOB_aligned_spin_waveform_quantities().__dict__,
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
