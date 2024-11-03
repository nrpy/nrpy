"""
Construct symbolic expression for the BOB aligned-spin gravitational-wave strain and NQC corrections.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
from typing import Any, List

import sympy as sp

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
        Shat = (m1 * m1 * chi1 + m2 * m2 * chi2) / (m1 * m1 + m2 * m2)
        Shat2 = Shat * Shat
        Shat3 = Shat2 * Shat
        Shat4 = Shat3 * Shat
        nu2 = nu * nu
        nu3 = nu2 * nu
        nu4 = nu3 * nu
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
        self.h = (Ap / 4 / (Omega**2)) * (1 / sp.cosh((t - t_p) / tau_qnm))
        Phi = arctan_p + arctanh_p - arctan_m - arctanh_m
        self.phi = 2 * Phi
        # mostly trivial
        # t_attach = t_0
        self.h_t_attach = h22NR
        self.hdot_t_attach = sp.sympify(0)
        hddot = sp.diff(sp.diff(self.h, t), t)
        self.hddot_t_attach = hddot.subs(t, t_0)
        self.w_t_attach = omega22NR
        Omegadot = sp.diff(Omega, t)
        self.wdot_t_attach = 2 * Omegadot.subs(t, t_0)


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
