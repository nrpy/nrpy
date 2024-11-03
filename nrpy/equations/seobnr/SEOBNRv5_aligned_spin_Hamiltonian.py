"""
Construct symbolic expression for the SEOBNRv5 aligned-spin Hamiltonian and associated quantities.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
import sympy as sp

# The name of this module ("WaveEquation") is given by __name__:
thismodule = __name__


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


class SEOBNRv5_aligned_spin_Hamiltonian_quantities:
    """Class for computing the SEOBNRv5 aligned-spin Hamiltonian and associated quantities."""

    def __init__(self) -> None:
        """
        Compute the SEOBNRv5 aligned-spin Hamiltonian.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the aligned-spin SEOBNRv5 Hamiltonian. It
        initializes class variables like mass parameters, spin parameters, and
        various coefficients required for the Hamiltonian's effective potential
        calculations.

        Inputs: 'm1', 'm2', 'r', 'prstar', 'pphi', 'chi1', 'chi2', 'a6', and 'dSO'
        Outputs: 'xi' and 'Hreal'
        """
        m1, m2, self.r, self.prstar, self.pphi, chi1, chi2, a6, dSO = sp.symbols(
            "m1 m2 r prstar pphi chi1 chi2 a6 dSO", real=True
        )
        u = 1 / self.r
        M = m1 + m2
        delta = (m1 - m2) / M
        self.nu = m1 * m2 / (M**2)
        gap = (
            sp.Rational(7, 4)
            + (self.pphi**2 / self.r**2)
            * (-sp.Rational(45, 32) * self.nu - sp.Rational(15, 32))
            + (1 / self.r) * (sp.Rational(23, 32) * self.nu - sp.Rational(3, 32))
            + (self.pphi**4 / self.r**4)
            * (
                sp.Rational(345, 256) * self.nu**2
                + sp.Rational(75, 128) * self.nu
                + sp.Rational(105, 256)
            )
            + (self.pphi**2 / self.r**3)
            * (
                -sp.Rational(1591, 768) * self.nu**2
                - sp.Rational(267, 128) * self.nu
                + sp.Rational(59, 256)
            )
            + (1 / self.r**2)
            * (
                sp.Rational(109, 192) * self.nu**2
                - sp.Rational(177, 32) * self.nu
                - sp.Rational(5, 64)
            )
        )
        gam = (
            sp.Rational(1, 4)
            + (self.pphi**2 / self.r**2)
            * (sp.Rational(15, 32) - sp.Rational(9, 32) * self.nu)
            + (1 / self.r) * (sp.Rational(11, 32) * self.nu + sp.Rational(3, 32))
            + (self.pphi**4 / self.r**4)
            * (
                sp.Rational(75, 256) * self.nu**2
                - sp.Rational(45, 128) * self.nu
                - sp.Rational(105, 256)
            )
            + (self.pphi**2 / self.r**3)
            * (
                -sp.Rational(613, 768) * self.nu**2
                - sp.Rational(35, 128) * self.nu
                - sp.Rational(59, 256)
            )
            + (1 / self.r**2)
            * (
                sp.Rational(103, 192) * self.nu**2
                - sp.Rational(1, 32) * self.nu
                + sp.Rational(5, 64)
            )
        )
        ap = (m1 * chi1 + m2 * chi2) / M
        am = (m1 * chi1 - m2 * chi2) / M
        Qnos = (
            (
                f2r(0.121954868780449) * self.nu * self.prstar**8 / self.r
                + self.prstar**6
                * (
                    6 * self.nu**3
                    - f2r(5.4) * self.nu**2
                    - f2r(2.78300763695006) * self.nu
                )
                / self.r**2
                + self.prstar**4
                * (10 * self.nu**3 - 131 * self.nu**2 + f2r(92.7110442849544) * self.nu)
                / self.r**3
            )
            + (
                self.prstar**8
                * (
                    -6 * self.nu**4
                    + f2r(3.42857142857143) * self.nu**3
                    + f2r(3.33842023648322) * self.nu**2
                    + f2r(1.38977750996128) * self.nu
                )
                / self.r**2
                + self.prstar**6
                * (
                    -14 * self.nu**4
                    + 188 * self.nu**3
                    - f2r(89.5298327361234) * self.nu**2
                    - f2r(33.9782122170436) * self.nu
                )
                / self.r**3
                + self.prstar**4
                * (
                    f2r(602.318540416564) * self.nu**3
                    + self.nu**2 * (118.4 * sp.log(self.r) - f2r(1796.13660498019))
                    + self.nu
                    * (f2r(452.542166996693) - f2r(51.6952380952381) * sp.log(self.r))
                )
                / self.r**4
            )
            + (
                f2r(1.48275342024365)
                * self.nu
                * self.prstar**8
                / self.r ** sp.Rational(5, 2)
                - f2r(11.3175085791863)
                * self.nu
                * self.prstar**6
                / self.r ** sp.Rational(7, 2)
                + f2r(147.443752990146)
                * self.nu
                * self.prstar**4
                / self.r ** sp.Rational(9, 2)
            )
            + self.prstar**4 * (-6 * self.nu**2 + 8 * self.nu) / self.r**2
        )
        d5 = 0
        Dnons = (
            self.r
            * (
                f2r(6730497718123.02) * self.nu**3
                + 22295347200 * self.nu**2 * d5
                + 133772083200 * self.nu**2 * self.r**2
                + f2r(1822680546449.21) * self.nu**2 * self.r
                + f2r(80059249540278.2) * self.nu**2
                + 22295347200 * self.nu * d5 * self.r
                - 193226342400 * self.nu * d5
                + f2r(2589101062873.81) * self.nu * self.r**2
                + f2r(10611661054566.2) * self.nu * self.r
                - f2r(12049908701745.2) * self.nu
                + f2r(5107745331375.71) * self.r**2
                - f2r(326837426.241486) * self.r * (14700 * self.nu + 42911)
                - f2r(39476764256925.6) * self.r
                - (
                    -5041721180160 * self.nu**2
                    - f2r(25392914995744.3) * self.nu
                    - 879923036160 * self.r**2
                    - 283115520 * self.r * (14700 * self.nu + 42911)
                    + 104186110149937
                )
                * sp.log(self.r)
                + 5787938193408 * sp.log(self.r) ** 2
                + 275059053208689
            )
            / (
                55296
                * self.nu
                * (
                    14515200 * self.nu**3
                    - f2r(42636451.6032331) * self.nu**2
                    - 7680 * self.nu * (315 * d5 + f2r(890888.810272497))
                    + f2r(4331361844.61149) * self.nu
                    + f2r(1002013764.01019)
                )
                - 967680
                * self.r**3
                * (
                    -138240 * self.nu**2
                    - f2r(2675575.66847905) * self.nu
                    - f2r(5278341.3229329)
                )
                - 9216
                * self.r**2
                * (
                    -f2r(197773496.793534) * self.nu**2
                    - 7680 * self.nu * (315 * d5 + f2r(405152.309729121))
                    + f2r(2481453539.84635) * self.nu
                    + f2r(5805304367.87913)
                )
                + self.r
                * (
                    f2r(5927865218923.02) * self.nu**3
                    + 70778880 * self.nu**2 * (315 * d5 + f2r(2561145.80918574))
                    - 138141470005001 * self.nu**2
                    - 4718592 * self.nu * (40950 * d5 + f2r(86207832.4415642))
                    + 450172889755120 * self.nu
                    + f2r(86618264430493.3) * (1 - 0.496948781616935 * self.nu) ** 2
                    + 188440788778196
                )
                + 5787938193408 * self.r * sp.log(self.r) ** 2
                + (
                    -1698693120 * self.nu * (11592 * self.nu + 69847)
                    + 879923036160 * self.r**3
                    + 283115520 * self.r**2 * (14700 * self.nu + 42911)
                    + 49152
                    * self.r
                    * (
                        102574080 * self.nu**2
                        + f2r(409207698.136075) * self.nu
                        - f2r(2119671837.36038)
                    )
                )
                * sp.log(self.r)
            )
        )
        Anons = (
            7680
            * self.r**4
            * (
                -f2r(5416406.59541186) * self.nu**2
                + 28 * self.nu * (1920 * a6 + f2r(733955.307463037))
                + 2048 * self.nu * (756 * self.nu + 336 * self.r + 407) * sp.log(self.r)
                - 7
                * self.r
                * (
                    -f2r(185763.092693281) * self.nu**2
                    + f2r(938918.400156317) * self.nu
                    - 245760
                )
                - 3440640
            )
            / (
                f2r(241555486248.807) * self.nu**4
                + 1120
                * self.nu**3
                * (
                    -f2r(17833256.898555) * self.r**2
                    - f2r(163683964.822551) * self.r
                    - f2r(1188987459.03162)
                )
                + 7
                * self.nu**2
                * (
                    -39321600 * a6 * (3 * self.r + 59)
                    + f2r(745857848.115604) * a6
                    + f2r(1426660551.8844) * self.r**5
                    - f2r(3089250703.76879) * self.r**4
                    - f2r(6178501407.53758) * self.r**3
                    + f2r(2064783811.32587) * self.r**2
                    + f2r(122635399361.987) * self.r
                    + f2r(276057889687.011)
                )
                + 67645734912 * self.nu**2 * sp.log(self.r) ** 2
                + 53760
                * self.nu
                * (
                    7680
                    * a6
                    * (self.r**4 + 2 * self.r**3 + 4 * self.r**2 + 8 * self.r + 16)
                    + 128
                    * self.r
                    * (
                        -f2r(6852.34813868015) * self.r**4
                        + f2r(4264.6962773603) * self.r**3
                        + f2r(8529.39255472061) * self.r**2
                        + f2r(13218.7851094412) * self.r
                        - f2r(33722.4297811176)
                    )
                    + f2r(113485.217444961)
                    * self.r
                    * (-self.r**4 + 2 * self.r**3 + 4 * self.r**2 + 8 * self.r + 16)
                    + f2r(148.04406601634)
                    * self.r
                    * (
                        349 * self.r**4
                        + 1926 * self.r**3
                        + 3852 * self.r**2
                        + 7704 * self.r
                        + 36400
                    )
                )
                + 32768
                * self.nu
                * (
                    -f2r(1882456.23663972) * self.nu**2
                    - f2r(38842241.4769507) * self.nu
                    + 161280 * self.r**5
                    + 480 * self.r**4 * (756 * self.nu + 1079)
                    + 960 * self.r**3 * (756 * self.nu + 1079)
                    + 1920 * self.r**2 * (588 * self.nu + 1079)
                    + 240
                    * self.r
                    * (-3024 * self.nu**2 - f2r(7466.27061066206) * self.nu + 17264)
                    + 13447680
                )
                * sp.log(self.r)
                + 13212057600 * self.r**5
            )
        )
        self.xi = sp.sqrt(Dnons) * (Anons + ap**2 * u * u) / (1 + ap**2 * u * u)
        pr = self.prstar / self.xi
        QalignSS = ((pr**4) / (self.r**3)) * (
            ap**2
            * (
                -5 * self.nu * self.nu
                + self.nu * sp.Rational(165, 32)
                + sp.Rational(25, 32)
            )
            + delta * ap * am * (self.nu * sp.Rational(45, 8) - sp.Rational(5, 16))
            + am**2
            * (
                -self.nu * self.nu * sp.Rational(15, 8)
                + self.nu * sp.Rational(75, 32)
                - sp.Rational(15, 32)
            )
        )
        BnpalignSS = (1 / self.r**3) * (
            ap**2 * (3 * self.nu + sp.Rational(45, 16))
            - delta * ap * am * sp.Rational(21, 8)
            + am**2 * (self.nu * sp.Rational(3, 4) - sp.Rational(3, 16))
        ) + (1 / self.r**4) * (
            ap**2 * (-self.nu * sp.Rational(1171, 64) - sp.Rational(861, 64))
            + delta * ap * am * (self.nu * sp.Rational(13, 16) + sp.Rational(449, 32))
            + am**2
            * (
                self.nu * self.nu * sp.Rational(1, 16)
                + self.nu * sp.Rational(115, 64)
                - sp.Rational(37, 64)
            )
        )
        AalignSS = (1 / self.r**4) * (
            ap**2 * sp.Rational(9, 8)
            - delta * ap * am * sp.Rational(5, 4)
            + am**2 * (self.nu * sp.Rational(1, 2) + sp.Rational(1, 8))
        ) + (1 / self.r**5) * (
            ap**2 * (-self.nu * sp.Rational(175, 64) - sp.Rational(225, 64))
            + delta * ap * am * (-self.nu * sp.Rational(39, 16) + sp.Rational(117, 32))
            + am**2
            * (
                self.nu * self.nu * sp.Rational(21, 16)
                - self.nu * sp.Rational(81, 64)
                - sp.Rational(9, 64)
            )
        )
        Qalign = Qnos + QalignSS
        Balignnp = -1 + ap**2 * u * u + Anons * Dnons + BnpalignSS
        Bkerreqnp = -(1 + 2 / self.r) / (self.r**2 + ap**2 * (1 + 2 / self.r))
        Aalign = (ap**2 * u * u + Anons + AalignSS) / (
            1 + ap**2 * (1 + 2 / self.r) / (self.r**2)
        )
        Galigna3 = self.pphi * (delta * am * ap**2 - ap**3) / (4 * self.r**2)
        SOcalib = self.nu * dSO * ap * self.pphi * (u**3)
        Heven = sp.sqrt(
            Aalign
            * (
                1
                + self.pphi * self.pphi / self.r**2
                + (1 + Balignnp) * pr * pr
                + Bkerreqnp * self.pphi * self.pphi * ap**2 / self.r**2
                + Qalign
            )
        )
        Hodd = (self.pphi * (gap * ap + delta * gam * am) + SOcalib + Galigna3) / (
            self.r**3 + ap**2 * (self.r + 2)
        )
        Heff = Hodd + Heven
        self.Hreal = sp.sqrt(1 + 2 * self.nu * (Heff - 1))
        self.dHreal_dr = sp.diff(self.Hreal, self.r) / self.nu
        self.dHreal_dprstar = sp.diff(self.Hreal, self.prstar) / self.nu
        self.dHreal_dpphi = sp.diff(self.Hreal, self.pphi) / self.nu
        self.dHreal_dr_dr = sp.diff(self.dHreal_dr, self.r) / self.nu
        self.dHreal_dr_dpphi = sp.diff(self.dHreal_dr, self.pphi) / self.nu
        Hreal_circ = self.Hreal.subs(self.prstar, 0)
        self.dHreal_dr_circ = sp.diff(Hreal_circ, self.r) / self.nu
        self.dHreal_dpphi_circ = sp.diff(Hreal_circ, self.pphi) / self.nu
        self.dHreal_dr_dr_circ = sp.diff(self.dHreal_dr_circ, self.r) / self.nu
        self.dHreal_dr_dpphi_circ = sp.diff(self.dHreal_dr_circ, self.pphi) / self.nu
        self.dHreal_dpphi_dpphi_circ = (
            sp.diff(self.dHreal_dpphi_circ, self.pphi) / self.nu
        )


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
        SEOBNRv5_aligned_spin_Hamiltonian_quantities().__dict__,
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
