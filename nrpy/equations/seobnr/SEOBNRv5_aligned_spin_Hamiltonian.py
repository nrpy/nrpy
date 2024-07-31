"""
Construct symbolic expression for the SEOBNRv5 aligned-spin Hamiltonian and associated quantities.

Author: Siddharth Mahesh

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
import sympy as sp
from mpmath import mpf

# The name of this module ("WaveEquation") is given by __name__:
thismodule = __name__


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
                mpf(0.121954868780449) * self.nu * self.prstar**8 / self.r
                + self.prstar**6
                * (
                    mpf(6.0) * self.nu**3
                    - mpf(5.4) * self.nu**2
                    - mpf(2.78300763695006) * self.nu
                )
                / self.r**2
                + self.prstar**4
                * (
                    mpf(10.0) * self.nu**3
                    - mpf(131.0) * self.nu**2
                    + mpf(92.7110442849544) * self.nu
                )
                / self.r**3
            )
            + (
                self.prstar**8
                * (
                    -mpf(6.0) * self.nu**4
                    + mpf(3.42857142857143) * self.nu**3
                    + mpf(3.33842023648322) * self.nu**2
                    + mpf(1.38977750996128) * self.nu
                )
                / self.r**2
                + self.prstar**6
                * (
                    -mpf(14.0) * self.nu**4
                    + mpf(188.0) * self.nu**3
                    - mpf(89.5298327361234) * self.nu**2
                    - mpf(33.9782122170436) * self.nu
                )
                / self.r**3
                + self.prstar**4
                * (
                    mpf(602.318540416564) * self.nu**3
                    + self.nu**2 * (mpf(118.4) * sp.log(self.r) - mpf(1796.13660498019))
                    + self.nu
                    * (mpf(452.542166996693) - mpf(51.6952380952381) * sp.log(self.r))
                )
                / self.r**4
            )
            + (
                mpf(1.48275342024365) * self.nu * self.prstar**8 / self.r**2.5
                - mpf(11.3175085791863) * self.nu * self.prstar**6 / self.r**3.5
                + mpf(147.443752990146) * self.nu * self.prstar**4 / self.r**4.5
            )
            + self.prstar**4 * (-mpf(6.0) * self.nu**2 + mpf(8.0) * self.nu) / self.r**2
        )
        d5 = 0
        Dnons = (
            self.r
            * (
                mpf(6730497718123.02) * self.nu**3
                + mpf(22295347200.0) * self.nu**2 * d5
                + mpf(133772083200.0) * self.nu**2 * self.r**2
                + mpf(1822680546449.21) * self.nu**2 * self.r
                + mpf(80059249540278.2) * self.nu**2
                + mpf(22295347200.0) * self.nu * d5 * self.r
                - mpf(193226342400.0) * self.nu * d5
                + mpf(2589101062873.81) * self.nu * self.r**2
                + mpf(10611661054566.2) * self.nu * self.r
                - mpf(12049908701745.2) * self.nu
                + mpf(5107745331375.71) * self.r**2
                - mpf(326837426.241486) * self.r * (14700.0 * self.nu + 42911.0)
                - mpf(39476764256925.6) * self.r
                - (
                    -mpf(5041721180160.0) * self.nu**2
                    - mpf(25392914995744.3) * self.nu
                    - mpf(879923036160.0) * self.r**2
                    - mpf(283115520.0) * self.r * (14700.0 * self.nu + 42911.0)
                    + mpf(104186110149937.0)
                )
                * sp.log(self.r)
                + mpf(5787938193408.0) * sp.log(self.r) ** 2
                + mpf(275059053208689.0)
            )
            / (
                mpf(55296.0)
                * self.nu
                * (
                    mpf(14515200.0) * self.nu**3
                    - mpf(42636451.6032331) * self.nu**2
                    - mpf(7680.0) * self.nu * (315.0 * d5 + mpf(890888.810272497))
                    + mpf(4331361844.61149) * self.nu
                    + mpf(1002013764.01019)
                )
                - mpf(967680.0)
                * self.r**3
                * (
                    -mpf(138240.0) * self.nu**2
                    - mpf(2675575.66847905) * self.nu
                    - mpf(5278341.3229329)
                )
                - mpf(9216.0)
                * self.r**2
                * (
                    -mpf(197773496.793534) * self.nu**2
                    - mpf(7680.0) * self.nu * (315.0 * d5 + mpf(405152.309729121))
                    + mpf(2481453539.84635) * self.nu
                    + mpf(5805304367.87913)
                )
                + self.r
                * (
                    mpf(5927865218923.02) * self.nu**3
                    + mpf(70778880.0)
                    * self.nu**2
                    * (315.0 * d5 + mpf(2561145.80918574))
                    - mpf(138141470005001.0) * self.nu**2
                    - mpf(4718592.0) * self.nu * (40950.0 * d5 + mpf(86207832.4415642))
                    + mpf(450172889755120.0) * self.nu
                    + mpf(86618264430493.3)
                    * (1 - mpf(0.496948781616935) * self.nu) ** 2
                    + mpf(188440788778196.0)
                )
                + mpf(5787938193408.0) * self.r * sp.log(self.r) ** 2
                + (
                    -mpf(1698693120.0) * self.nu * (11592.0 * self.nu + 69847.0)
                    + mpf(879923036160.0) * self.r**3
                    + mpf(283115520.0) * self.r**2 * (14700.0 * self.nu + 42911.0)
                    + mpf(49152.0)
                    * self.r
                    * (
                        mpf(102574080.0) * self.nu**2
                        + mpf(409207698.136075) * self.nu
                        - mpf(2119671837.36038)
                    )
                )
                * sp.log(self.r)
            )
        )
        Anons = (
            7680.0
            * self.r**4
            * (
                -mpf(5416406.59541186) * self.nu**2
                + 28.0 * self.nu * (1920.0 * a6 + mpf(733955.307463037))
                + 2048.0
                * self.nu
                * (756.0 * self.nu + 336.0 * self.r + 407.0)
                * sp.log(self.r)
                - 7.0
                * self.r
                * (
                    -mpf(185763.092693281) * self.nu**2
                    + mpf(938918.400156317) * self.nu
                    - 245760.0
                )
                - 3440640.0
            )
            / (
                mpf(241555486248.807) * self.nu**4
                + 1120.0
                * self.nu**3
                * (
                    -mpf(17833256.898555) * self.r**2
                    - mpf(163683964.822551) * self.r
                    - mpf(1188987459.03162)
                )
                + 7.0
                * self.nu**2
                * (
                    -39321600.0 * a6 * (3.0 * self.r + 59.0)
                    + mpf(745857848.115604) * a6
                    + mpf(1426660551.8844) * self.r**5
                    - mpf(3089250703.76879) * self.r**4
                    - mpf(6178501407.53758) * self.r**3
                    + mpf(2064783811.32587) * self.r**2
                    + mpf(122635399361.987) * self.r
                    + mpf(276057889687.011)
                )
                + 67645734912.0 * self.nu**2 * sp.log(self.r) ** 2
                + 53760.0
                * self.nu
                * (
                    7680.0
                    * a6
                    * (
                        self.r**4
                        + 2.0 * self.r**3
                        + 4.0 * self.r**2
                        + 8.0 * self.r
                        + 16.0
                    )
                    + 128.0
                    * self.r
                    * (
                        -mpf(6852.34813868015) * self.r**4
                        + mpf(4264.6962773603) * self.r**3
                        + mpf(8529.39255472061) * self.r**2
                        + mpf(13218.7851094412) * self.r
                        - mpf(33722.4297811176)
                    )
                    + mpf(113485.217444961)
                    * self.r
                    * (
                        -self.r**4
                        + 2.0 * self.r**3
                        + 4.0 * self.r**2
                        + 8.0 * self.r
                        + 16.0
                    )
                    + mpf(148.04406601634)
                    * self.r
                    * (
                        349.0 * self.r**4
                        + 1926.0 * self.r**3
                        + 3852.0 * self.r**2
                        + 7704.0 * self.r
                        + 36400.0
                    )
                )
                + 32768.0
                * self.nu
                * (
                    -mpf(1882456.23663972) * self.nu**2
                    - mpf(38842241.4769507) * self.nu
                    + mpf(161280.0) * self.r**5
                    + 480.0 * self.r**4 * (756.0 * self.nu + 1079.0)
                    + 960.0 * self.r**3 * (756.0 * self.nu + 1079.0)
                    + 1920.0 * self.r**2 * (588.0 * self.nu + 1079.0)
                    + 240.0
                    * self.r
                    * (-3024.0 * self.nu**2 - mpf(7466.27061066206) * self.nu + 17264.0)
                    + 13447680.0
                )
                * sp.log(self.r)
                + 13212057600.0 * self.r**5
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

        # define the derivatives
        self.dHreal_dr = sp.diff(self.Hreal, self.r) / self.nu
        self.dHreal_dprstar = sp.diff(self.Hreal, self.prstar) / self.nu
        self.dHreal_dpphi = sp.diff(self.Hreal, self.pphi) / self.nu
        self.dHreal_dr_dr = sp.diff(self.dHreal_dr, self.r) / self.nu
        self.dHreal_dr_dpphi = sp.diff(self.dHreal_dr, self.pphi) / self.nu
        Hreal_circ = self.Hreal.subs(self.prstar, 0)
        self.dHreal_dr_circ = sp.diff(Hreal_circ, self.r) / self.nu
        self.dHreal_dpphi_circ = sp.diff(Hreal_circ, self.pphi) / self.nu


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
