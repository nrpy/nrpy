"""
Construct the SEOBNRv5 Hamiltonian for spin-precessing binaries in the partially-precessing (quasi-precessing) limit.

Authors: Siddharth Mahesh
sm0193 at mix dot wvu dot edu
Zachariah B. Etienne
zachetie at gmail *dot com

The Hamiltonian is central to the Effective-One-Body (EOB) formalism, which maps the
conservative part of the two-body problem to the motion of an effective particle in a deformed metric.
This approach combines post-Newtonian (PN) theory with numerical relativity (NR)
calibrations to model the inspiral phase of a binary black hole merger.

The Hamiltonian is expressed in terms of the binary masses (m1, m2),
vectorial spins (chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z),
their projections onto the orbital angular momentum (chi1_l, chi2_l),
their projections onto the Newtonian angular momentum (chi1_ln, chi2_ln),
and canonical variables (r, phi, pr, pphi). This implementation is based on
the SEOBNRv5HM model (see Section III.A of https://arxiv.org/pdf/2303.18143 for the full list of terms)
and is used in the SEBOBv2 formalism described in Mahesh, McWilliams, and Etienne,
"Spinning Effective-to-Backwards-One Body".
License: BSD 2-Clause
"""

# Step P1: Import needed modules:
import sympy as sp

from nrpy.helpers.float_to_rational import f2r

# The name of this module ("WaveEquation") is given by __name__:
thismodule = __name__


class SEOBNRv5_quasi_precessing_spin_Hamiltonian_quantities:
    """Class for computing the SEOBNRv5 generic-spin Hamiltonian and associated quantities."""

    def __init__(self) -> None:
        """
        Compute the SEOBNRv5 generic-spin Hamiltonian.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the generic-spin SEOBNRv5 Hamiltonian. It
        initializes class variables like mass parameters, spin parameters, and
        various coefficients required for the Hamiltonian's effective potential
        calculations.

        Care needs to be taken in documenting the precessing v/s partially precessing
        equations as the orbit averaging and assumptions of circularity result in
        different expressions for the metric potential. Without the correct equation numbers,
        it can be hard to visually differentiate between the two cases and can lead
        to false positives in bug analyses.

        All potentials entering the effective Hamiltonian are computed under the assumption
        that the system of interest is a black hole binary. That is all multipolar coefficients
        hat{C}_{*} in Equations 4 and 5 of https://arxiv.org/pdf/2303.18143 are set to zero,
        and their contributions to the potentials are neglected in the code.

        The key outputs of the SEOBNRv5_generic_spin_Hamiltonian_quantities class are:
            - 'xi': The tortoise parameter
                        (Equation 44 of https://arxiv.org/pdf/2303.18143).
            - 'Hreal': The Hamiltonian
                        (Equation 7 of https://arxiv.org/pdf/2303.18143).
            - First derivatives of the Hamiltonian
                        needed for ODE integration.
            - First and second derivatives of the Hamiltonian
                        needed for initial conditions.
            - The instantaneous and circular frequency
                        needed for waveform and flux calculations.


        Inputs: 'm1', 'm2', 'r', 'phi', 'prstar', 'pphi', 'chi1_x', 'chi1_y', 'chi1_z', 'chi1_l', 'chi1_ln', 'chi2_x', 'chi2_y', 'chi2_z', 'chi2_l', 'chi2_ln', 'a6', and 'dSO'
        Outputs: 'xi' and 'Hreal'
        """
        (
            m1,
            m2,
            self.r,
            self.phi,
            self.prstar,
            self.pphi,
            self.chi1_x,  # x-component of the spin of the first black hole
            self.chi1_y,  # y-component of the spin of the first black hole
            self.chi1_z,  # z-component of the spin of the first black hole
            self.chi1_l,  # projection of the spin of the first black hole onto the orbital angular momentum
            self.chi1_ln,  # projection of the spin of the first black hole onto the Newtonian orbital angular momentum (needed for orbit averaging)
            self.chi2_x,  # x-component of the spin of the second black hole
            self.chi2_y,  # y-component of the spin of the second black hole
            self.chi2_z,  # z-component of the spin of the second black hole
            self.chi2_l,  # projection of the spin of the second black hole onto the orbital angular momentum
            self.chi2_ln,  # projection of the spin of the second black hole onto the Newtonian angular momentum (needed for orbit averaging)
            a6,
            dSO,
        ) = sp.symbols(
            "m1 m2 r phi prstar pphi chi1_x chi1_y chi1_z chi1_l chi1_ln chi2_x chi2_y chi2_z chi2_l chi2_ln a6 dSO",
            real=True,
        )
        # mass combinations
        M = m1 + m2
        m1_norm = m1 / M
        m2_norm = m2 / M
        delta = m1_norm - m2_norm
        self.nu = m1_norm * m2_norm
        # spin combinations
        ap_x = m1_norm * self.chi1_x + m2_norm * self.chi2_x
        ap_y = m1_norm * self.chi1_y + m2_norm * self.chi2_y
        ap_z = m1_norm * self.chi1_z + m2_norm * self.chi2_z
        apsq = ap_x**2 + ap_y**2 + ap_z**2
        am_x = m1_norm * self.chi1_x - m2_norm * self.chi2_x
        am_y = m1_norm * self.chi1_y - m2_norm * self.chi2_y
        am_z = m1_norm * self.chi1_z - m2_norm * self.chi2_z
        amsq = am_x**2 + am_y**2 + am_z**2
        ap_dot_am = ap_x * am_x + ap_y * am_y + ap_z * am_z
        ap_dot_ln = m1_norm * self.chi1_ln + m2_norm * self.chi2_ln
        am_dot_ln = m1_norm * self.chi1_ln - m2_norm * self.chi2_ln
        ap_dot_ln_times_am_dot_ln = ap_dot_ln * am_dot_ln
        ap_l = m1_norm * self.chi1_l + m2_norm * self.chi2_l
        am_l = m1_norm * self.chi1_l - m2_norm * self.chi2_l
        ap_dot_l = self.pphi * ap_l
        am_dot_l = self.pphi * am_l
        # orbit averaged version of the spin combinations
        # see 54 of https://arxiv.org/pdf/2303.18143
        n_dot_ap_2_avg = sp.Rational(1, 2) * (apsq - ap_dot_ln**2)
        n_dot_am_2_avg = sp.Rational(1, 2) * (amsq - am_dot_ln**2)
        n_dot_ap_times_n_dot_am_avg = sp.Rational(1, 2) * (
            ap_dot_am - ap_dot_ln_times_am_dot_ln
        )
        # define inverse radius
        u = 1 / self.r
        # non-spinning metric potentials
        # Equation 25 of https://arxiv.org/pdf/2303.18143
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
        # below two terms are given as Taylor expanded potentials
        # in
        # however, the Pade approximant is what we want.
        # The Pade approximated versions are taken from
        # https://git.ligo.org/waveforms/software/pyseobnr/-/blob/2aebc3347dfade6d1df282ae6eedd34083a980b0/pyseobnr/eob/hamiltonian/Ham_AvgS2precess_simple_cython_PA_AD.pyx
        # Dnons = Dbpm (line 111 ONLY in the above permalink)
        # Anons = Apm (line 113 ONLY in the above permalink)
        # d5 = 0 (Appendix A of https://arxiv.org/pdf/2303.18039)
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
        # tortoise factor xi (Equation 44 of https://arxiv.org/pdf/2303.18143)
        self.xi = sp.sqrt(Dnons) * (Anons + apsq * u * u) / (1 + apsq * u * u)
        # radial momentum pr in terms of tortoise momentum prstar
        # (Equation 42 of https://arxiv.org/pdf/2303.18143)
        pr = self.prstar / self.xi
        # odd-in-spin potentials
        # gyro-gravitomagnetic factors
        # v5PHM only uses up to 3.5 PN for these terms
        # Equation 28a of https://arxiv.org/pdf/2303.18143
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
        # Equation 28b of https://arxiv.org/pdf/2303.18143
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
        # SO calibration factor
        # Equation 29 of https://arxiv.org/pdf/2303.18143
        SOcalib = self.nu * dSO * ap_dot_l * (u**3)
        # cubic-in-spin SO terms
        # Equation 50 of https://arxiv.org/pdf/2303.18143
        # same as Equation 37 but with with p_r = 0
        # and some values changed due to the
        # orbit averaging procedure
        Gppreca3_avg = ap_dot_l * (
            self.pphi**2
            * u**3
            * (
                sp.Rational(1, 2) * delta * n_dot_ap_times_n_dot_am_avg
                - sp.Rational(1, 4) * n_dot_ap_2_avg
            )
            + u**2
            * (
                -sp.Rational(1, 4) * apsq
                - sp.Rational(3, 4) * n_dot_ap_2_avg
                + sp.Rational(5, 24) * delta * ap_dot_am
                - sp.Rational(5, 3) * delta * n_dot_ap_times_n_dot_am_avg
            )
        ) + am_dot_l * (
            -sp.Rational(1, 4) * self.pphi**2 * u**3 * delta * n_dot_ap_2_avg
            + u**2
            * (
                sp.Rational(1, 24) * delta * apsq
                + sp.Rational(5, 12) * delta * n_dot_ap_2_avg
            )
        )
        # quadratic-in-spin metric potentials
        # in-plane spin components
        # Equation 49a of https://arxiv.org/pdf/2303.18143
        AinplaneSS_avg = (
            2 * u**3 * n_dot_ap_2_avg
            + u**4
            * (
                sp.Rational(33, 8) * delta * n_dot_ap_times_n_dot_am_avg
                + (-sp.Rational(1, 2) * self.nu - sp.Rational(3, 8)) * n_dot_am_2_avg
                + (sp.Rational(7, 4) * self.nu - sp.Rational(31, 4)) * n_dot_ap_2_avg
            )
            + u**5
            * (
                +(17 * self.nu - sp.Rational(1, 4))
                * delta
                * n_dot_ap_times_n_dot_am_avg
                + (
                    -sp.Rational(41, 8) * self.nu**2
                    + sp.Rational(583, 2) * self.nu
                    - sp.Rational(171, 64)
                )
                * n_dot_am_2_avg
                + (
                    -sp.Rational(11, 8) * self.nu**2
                    + sp.Rational(1435, 96) * self.nu
                    + sp.Rational(187, 64)
                )
                * n_dot_ap_2_avg
            )
        )
        # Equation 49b of https://arxiv.org/pdf/2303.18143
        BinplanepSS_avg = (
            -(u**2) * n_dot_ap_2_avg
            + u**3
            * (
                sp.Rational(3, 4) * delta * n_dot_ap_times_n_dot_am_avg
                + (-sp.Rational(3, 4) * self.nu + sp.Rational(3, 16)) * n_dot_am_2_avg
                + (-sp.Rational(7, 4) * self.nu - sp.Rational(15, 16)) * n_dot_ap_2_avg
            )
            + u**4
            * (
                +(-sp.Rational(49, 4) * self.nu - sp.Rational(43, 8))
                * delta
                * n_dot_ap_times_n_dot_am_avg
                + (
                    sp.Rational(19, 8) * self.nu**2
                    - sp.Rational(545, 32) * self.nu
                    + sp.Rational(219, 64)
                )
                * n_dot_am_2_avg
                + (
                    sp.Rational(11, 8) * self.nu**2
                    - sp.Rational(805, 96) * self.nu
                    + sp.Rational(125, 64)
                )
                * n_dot_ap_2_avg
            )
        )
        # out-of-plane spin components, same as aligned for quasi-precessing
        # binaries with values replaced by dot products.
        # Equation 39f of https://arxiv.org/pdf/2303.18143
        QalignSS = ((pr**4) / (self.r**3)) * (
            apsq
            * (
                -5 * self.nu * self.nu
                + self.nu * sp.Rational(165, 32)
                + sp.Rational(25, 32)
            )
            + delta * ap_dot_am * (self.nu * sp.Rational(45, 8) - sp.Rational(5, 16))
            + amsq
            * (
                -self.nu * self.nu * sp.Rational(15, 8)
                + self.nu * sp.Rational(75, 32)
                - sp.Rational(15, 32)
            )
        )
        # Equation 39d of https://arxiv.org/pdf/2303.18143
        BnpalignSS = (1 / self.r**3) * (
            apsq * (3 * self.nu + sp.Rational(45, 16))
            - delta * ap_dot_am * sp.Rational(21, 8)
            + amsq * (self.nu * sp.Rational(3, 4) - sp.Rational(3, 16))
        ) + (1 / self.r**4) * (
            apsq * (-self.nu * sp.Rational(1171, 64) - sp.Rational(861, 64))
            + delta * ap_dot_am * (self.nu * sp.Rational(13, 16) + sp.Rational(449, 32))
            + amsq
            * (
                self.nu * self.nu * sp.Rational(1, 16)
                + self.nu * sp.Rational(115, 64)
                - sp.Rational(37, 64)
            )
        )
        # Equation 39a of https://arxiv.org/pdf/2303.18143
        AalignSS = (1 / self.r**4) * (
            apsq * sp.Rational(9, 8)
            - delta * ap_dot_am * sp.Rational(5, 4)
            + amsq * (self.nu * sp.Rational(1, 2) + sp.Rational(1, 8))
        ) + (1 / self.r**5) * (
            apsq * (-self.nu * sp.Rational(175, 64) - sp.Rational(225, 64))
            + delta
            * ap_dot_am
            * (-self.nu * sp.Rational(39, 16) + sp.Rational(117, 32))
            + amsq
            * (
                self.nu * self.nu * sp.Rational(21, 16)
                - self.nu * sp.Rational(81, 64)
                - sp.Rational(9, 64)
            )
        )
        # metric potentials
        # Equation 48 of https://arxiv.org/pdf/2303.18143
        Qpprec = Qnos + QalignSS
        Bpprecp = 1 + BinplanepSS_avg
        Bpprecnp = -1 + apsq * u * u + Anons * Dnons + BnpalignSS
        # Equation 17 of https://arxiv.org/pdf/2303.18143
        Bkerreqnpa = -(1 + 2 / self.r) / (self.r**2 + apsq * (1 + 2 / self.r))
        Apprec = (apsq * u * u + Anons + AalignSS + AinplaneSS_avg) / (
            1 + apsq * (1 + 2 / self.r) / (self.r**2)
        )
        # even-in-spin Hamiltonian
        # Equation 47 of https://arxiv.org/pdf/2303.18143
        # broken into two parts.
        Heven = sp.sqrt(
            Apprec
            * (
                1
                + Bpprecp * self.pphi * self.pphi / self.r**2
                + (1 + Bpprecnp) * pr * pr
                + Bkerreqnpa * ap_dot_l**2 / self.r**2
                + Qpprec
            )
        )
        # odd-in-spin Hamiltonian
        Hodd = ((gap * ap_dot_l + delta * gam * am_dot_l) + SOcalib + Gppreca3_avg) / (
            self.r**3 + apsq * (self.r + 2)
        )
        Heff = Hodd + Heven
        # Equation 7 of https://arxiv.org/pdf/2303.18143
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
        SEOBNRv5_quasi_precessing_spin_Hamiltonian_quantities().__dict__,
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
