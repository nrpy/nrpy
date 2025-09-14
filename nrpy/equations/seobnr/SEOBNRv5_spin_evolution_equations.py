"""
Construct the SEOBNRv5 spin evolution equations for spin-precessing binaries.

Authors: Siddharth Mahesh
sm0193 at mix dot wvu dot edu
Zachariah B. Etienne
zachetie at gmail *dot com

In the 5th version of the Spinning Effective One Body (SEOBNRv5) model,
the orbital quantities and the spins are evolved separately, with the spin evolution equations are informed
by Taylor expanded post-Newtonian expressions in the quasi-circular limit.

The spin evolution equations are a set of equations that depend on
the binary masses (m1,m2), vectorial spins (chi1_x, chi1_y, chi1_z, chi2_x, chi2_y, chi2_z),
the normal vector (ln_x, ln_y, ln_z) (the unit vector in the direction of the Newtonian angular momentum),
and the orbital frequency (omega).

The equations describe the time derivatives of the spin vectors, the orbital plane and the orbital frequency.
They also specify the orbital angular momentum in terms of the dependent variables.
They are used to evolve the spins and angular momenta in the SEOBNRv5PHM model
(see Section III.C-G of https://arxiv.org/pdf/2303.18143 for the full list of terms)
and are used in the SEBOBv2 formalism described in Mahesh, McWilliams, and Etienne,
"Spinning Effective-to-Backwards-One Body".
License: BSD 2-Clause
"""

# Step P1: Import needed modules:
import sympy as sp

from nrpy.helpers.float_to_rational import f2r

# The name of this module ("WaveEquation") is given by __name__:
thismodule = __name__


class SEOBNRv5_spin_evolution_equations:
    """Class for computing the SEOBNRv5 spin evolution equations."""

    def __init__(self) -> None:
        """
        Compute the SEOBNRv5 spin evolution equations.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the SEOBNRv5 spin evolution equations. It
        initializes class variables like mass parameters, spin parameters, and
        various coefficients required for the dynamics of the spins, angular momenta
        and PN velocity parameter v.

        All potentials entering the effective Hamiltonian are computed under the assumption
        that the system of interest is a black hole binary. That is all multipolar coefficients
        hat{C}_{*} in Equations 4 and 5 of https://arxiv.org/pdf/2303.18143 are set to zero,
        and their contributions to the potentials are neglected in the code.

        The key outputs of the SEOBNRv5_generic_spin_Hamiltonian_quantities class are:
            - 'S1dot': The time derivative of the spin vector of the first black hole
                        (Equation 66 of https://arxiv.org/pdf/2303.18143).
            - 'S2dot': The time derivative of the spin vector of the second black hole
                        (Equation 66 of https://arxiv.org/pdf/2303.18143).
            - 'lndot': The time derivative of the Newtonian angular momentum unit vector
                        (Equation 71 of https://arxiv.org/pdf/2303.18143).
            - 'omegadot': The time derivative of the orbital frequency
                        (Equation 69 of https://arxiv.org/pdf/2303.18143).
            - 'L'   : The orbital angular momentum
                        (Equation 65 of https://arxiv.org/pdf/2303.18143).
        """
        (
            self.m1,
            self.m2,
            self.omega,
            self.chi1_x,
            self.chi1_y,
            self.chi1_z,
            self.chi2_x,
            self.chi2_y,
            self.chi2_z,
            self.ln_x,
            self.ln_y,
            self.ln_z,
        ) = sp.symbols(
            "m1 m2 omega chi1_x chi1_y chi1_z chi2_x chi2_y chi2_z ln_x ln_y ln_z",
            real=True,
        )
        EulerGamma = f2r(0.577215664901532860606512090082402431)
        # use sympy.Matrix to define vectors
        chi1 = sp.Matrix([self.chi1_x, self.chi1_y, self.chi1_z])
        chi2 = sp.Matrix([self.chi2_x, self.chi2_y, self.chi2_z])
        ln = sp.Matrix([self.ln_x, self.ln_y, self.ln_z])
        # mass combinations
        M = self.m1 + self.m2
        m1_norm = self.m1 / M
        m2_norm = self.m2 / M
        nu = m1_norm * m2_norm
        delta = m1_norm - m2_norm
        # spin combinations
        S1 = chi1 * self.m1**2
        S2 = chi2 * self.m2**2
        ln_dot_S1 = ln.dot(S1)
        ln_dot_S2 = ln.dot(S2)
        S1_dot_S2 = S1.dot(S2)
        S1sq = S1.dot(S1)
        S2sq = S2.dot(S2)
        ln_cross_S1 = ln.cross(S1)
        ln_cross_S2 = ln.cross(S2)
        S1_cross_S2 = S1.cross(S2)
        S2_cross_S1 = S2.cross(S1)
        ln_dot_S1_cross_S2 = ln.dot(S1_cross_S2)
        # PN expansion parameter v
        v = self.omega ** sp.Rational(1, 3)
        # S1_dot: Equation 66 of https://arxiv.org/pdf/2303.18143
        S1_dot = (
            ln_cross_S1
            * (
                v**5 * (sp.Rational(3, 2) * m2_norm + sp.Rational(1, 2) * nu)
                + v**7
                * (
                    (sp.Rational(9, 8) - sp.Rational(5, 4) * nu) * m2_norm
                    - sp.Rational(1, 24) * nu**2
                    + sp.Rational(15, 8) * nu
                )
                + v**9
                * (
                    (
                        sp.Rational(5, 16) * nu**2
                        - sp.Rational(39, 4) * nu
                        + sp.Rational(27, 16)
                    )
                    * m2_norm
                    - sp.Rational(1, 48) * nu**3
                    - sp.Rational(55, 16) * nu**2
                    + sp.Rational(81, 16) * nu
                )
            )
            + v**6
            * (
                ln_cross_S1
                * (
                    ln_dot_S1 * (sp.Rational(3, 2) * nu - sp.Rational(3, 2) * m2_norm)
                    - sp.Rational(3, 2) * nu * ln_dot_S2
                )
                + sp.Rational(1, 2) * nu * S2_cross_S1
            )
            + v**8
            * (
                ln_cross_S1
                * (
                    ln_dot_S1
                    * (
                        -sp.Rational(17, 12) * nu**2
                        - sp.Rational(9, 4) * nu
                        + (sp.Rational(9, 4) - sp.Rational(15, 4)) * m2_norm
                    )
                    + ln_dot_S2 * (sp.Rational(1, 12) * nu**2 - sp.Rational(1, 2) * nu)
                )
                - sp.Rational(1, 4) * nu**2 * S2_cross_S1
            )
            + v**10
            * (
                ln_cross_S1
                * (
                    ln_dot_S1
                    * (
                        sp.Rational(121, 144) * nu**3
                        - sp.Rational(91, 16) * nu**2
                        - sp.Rational(27, 16) * nu
                        + (
                            sp.Rational(385, 48) * nu**2
                            - sp.Rational(97, 16) * nu
                            + sp.Rational(27, 16)
                        )
                        * m2_norm
                    )
                    + ln_dot_S2
                    * (
                        sp.Rational(103, 144) * nu**3
                        + sp.Rational(139, 48) * nu**2
                        - sp.Rational(9, 4) * nu
                    )
                )
                + S2_cross_S1
                * (
                    sp.Rational(1, 48) * nu**3
                    + sp.Rational(49, 16) * nu**2
                    + sp.Rational(3, 8) * nu
                )
            )
        )
        S2_dot = (
            ln_cross_S2
            * (
                v**5 * (sp.Rational(3, 2) * m1_norm + sp.Rational(1, 2) * nu)
                + v**7
                * (
                    (sp.Rational(9, 8) - sp.Rational(5, 4) * nu) * m1_norm
                    - sp.Rational(1, 24) * nu**2
                    + sp.Rational(15, 8) * nu
                )
                + v**9
                * (
                    (
                        sp.Rational(5, 16) * nu**2
                        - sp.Rational(39, 4) * nu
                        + sp.Rational(27, 16)
                    )
                    * m1_norm
                    - sp.Rational(1, 48) * nu**3
                    - sp.Rational(55, 16) * nu**2
                    + sp.Rational(81, 16) * nu
                )
            )
            + v**6
            * (
                ln_cross_S2
                * (
                    ln_dot_S2 * (sp.Rational(3, 2) * nu - sp.Rational(3, 2) * m1_norm)
                    - sp.Rational(3, 2) * nu * ln_dot_S1
                )
                + sp.Rational(1, 2) * nu * S1_cross_S2
            )
            + v**8
            * (
                ln_cross_S2
                * (
                    ln_dot_S2
                    * (
                        -sp.Rational(17, 12) * nu**2
                        - sp.Rational(9, 4) * nu
                        + (sp.Rational(9, 4) - sp.Rational(15, 4)) * m1_norm
                    )
                    + ln_dot_S1 * (sp.Rational(1, 12) * nu**2 - sp.Rational(1, 2) * nu)
                )
                - sp.Rational(1, 4) * nu**2 * S1_cross_S2
            )
            + v**10
            * (
                ln_cross_S2
                * (
                    ln_dot_S2
                    * (
                        sp.Rational(121, 144) * nu**3
                        - sp.Rational(91, 16) * nu**2
                        - sp.Rational(27, 16) * nu
                        + (
                            sp.Rational(385, 48) * nu**2
                            - sp.Rational(97, 16) * nu
                            + sp.Rational(27, 16)
                        )
                        * m1_norm
                    )
                    + ln_dot_S1
                    * (
                        sp.Rational(103, 144) * nu**3
                        + sp.Rational(139, 48) * nu**2
                        - sp.Rational(9, 4) * nu
                    )
                )
                + S1_cross_S2
                * (
                    sp.Rational(1, 48) * nu**3
                    + sp.Rational(49, 16) * nu**2
                    + sp.Rational(3, 8) * nu
                )
            )
        )
        # v_dot: Equation 69 of https://arxiv.org/pdf/2303.18143
        # it is given as a combination of
        # nospin + spin-orbit + spin1-spin2 + spin-squared
        # Equation 69b of https://arxiv.org/pdf/2303.18143
        v_dot_nospin = (
            1
            + v**2 * (-sp.Rational(11, 4) * nu - sp.Rational(743, 336))
            + 4 * sp.pi * v**3
            + v**4
            * (
                sp.Rational(59, 18) * nu**2
                + sp.Rational(13661, 2016) * nu
                + sp.Rational(34103, 18144)
            )
            + sp.pi * v**5 * (-sp.Rational(189, 8) * nu - sp.Rational(4159, 672))
            + v**6
            * (
                sp.Rational(541, 896) * nu**2
                - sp.Rational(5605, 2592) * nu**3
                - sp.Rational(56198689, 217728) * nu
                + sp.pi**2 * (sp.Rational(451, 48) * nu + sp.Rational(16, 3))
                - sp.Rational(1712, 105) * sp.log(v)
                - sp.Rational(1712, 105) * EulerGamma
                + sp.Rational(16447322263, 139708800)
                - sp.Rational(3424, 105) * sp.log(2)
            )
            + sp.pi
            * v**7
            * (
                sp.Rational(91495, 1512) * nu**2
                + sp.Rational(358675, 6048) * nu
                - sp.Rational(4415, 4032)
            )
        )
        # Equation 69c of https://arxiv.org/pdf/2303.18143
        v_dot_spin_orbit = (
            # S1 part
            ln_dot_S1
            * (
                v**3 * (-sp.Rational(19, 6) * nu - sp.Rational(25, 4) * m2_norm)
                + v**5
                * (
                    sp.Rational(79, 6) * nu**2
                    - sp.Rational(21611, 1008) * nu
                    + (sp.Rational(281, 8) * nu - sp.Rational(809, 84)) * m2_norm
                )
                + sp.pi
                * v**6
                * (-sp.Rational(37, 3) * nu - sp.Rational(151, 6) * m2_norm)
                + v**7
                * (
                    -sp.Rational(10819, 432) * nu**3
                    + sp.Rational(40289, 288) * nu**2
                    - sp.Rational(1932041, 18144) * nu
                    + (
                        -sp.Rational(2903, 32) * nu**2
                        + sp.Rational(257023, 1008) * nu
                        - sp.Rational(1195759, 18144)
                    )
                    * m2_norm
                )
                + sp.pi
                * v**8
                * (
                    sp.Rational(34303, 336) * nu**2
                    - sp.Rational(46957, 504) * nu
                    + (sp.Rational(50483, 224) * nu - sp.Rational(1665, 28)) * m2_norm
                )
            )
            # 1 <-> 2 part
            + ln_dot_S2
            * (
                v**3 * (-sp.Rational(19, 6) * nu - sp.Rational(25, 4) * m1_norm)
                + v**5
                * (
                    sp.Rational(79, 6) * nu**2
                    - sp.Rational(21611, 1008) * nu
                    + (sp.Rational(281, 8) * nu - sp.Rational(809, 84)) * m1_norm
                )
                + sp.pi
                * v**6
                * (-sp.Rational(37, 3) * nu - sp.Rational(151, 6) * m1_norm)
                + v**7
                * (
                    -sp.Rational(10819, 432) * nu**3
                    + sp.Rational(40289, 288) * nu**2
                    - sp.Rational(1932041, 18144) * nu
                    + (
                        -sp.Rational(2903, 32) * nu**2
                        + sp.Rational(257023, 1008) * nu
                        - sp.Rational(1195759, 18144)
                    )
                    * m1_norm
                )
                + sp.pi
                * v**8
                * (
                    sp.Rational(34303, 336) * nu**2
                    - sp.Rational(46957, 504) * nu
                    + (sp.Rational(50483, 224) * nu - sp.Rational(1665, 28)) * m1_norm
                )
            )
        )
        # Equation 69d of https://arxiv.org/pdf/2303.18143
        v_dot_spin1_spin2 = nu * (
            v**4
            * (
                sp.Rational(721, 48) * ln_dot_S1 * ln_dot_S2
                - sp.Rational(247, 48) * S1_dot_S2
            )
            + v**6
            * (
                (sp.Rational(14433, 224) - sp.Rational(11779, 288) * nu)
                * ln_dot_S1
                * ln_dot_S2
                + (sp.Rational(6373, 288) * nu + sp.Rational(16255, 672)) * S1_dot_S2
            )
            + sp.pi
            * v**7
            * (sp.Rational(207, 4) * ln_dot_S1 * ln_dot_S2 - 12 * S1_dot_S2)
            + v**8
            * (
                (
                    -sp.Rational(162541, 3456) * nu**2
                    - sp.Rational(195697, 896) * nu
                    - sp.Rational(9355721, 72576)
                )
                * S1_dot_S2
                + (
                    sp.Rational(33163, 3456) * nu**2
                    - sp.Rational(10150387, 24192) * nu
                    + sp.Rational(21001565, 24192)
                )
                * ln_dot_S1
                * ln_dot_S2
            )
        )
        # Equation 69e of https://arxiv.org/pdf/2303.18143
        v_dot_spin_squared = (
            # S1 part
            (
                v**4
                * (
                    ln_dot_S1**2
                    * (sp.Rational(719, 96) * m2_norm - sp.Rational(719, 96) * nu)
                    + S1sq
                    * (sp.Rational(233, 96) * nu - sp.Rational(233, 96) * m2_norm)
                )
                + v**6
                * (
                    ln_dot_S1**2
                    * (
                        sp.Rational(25373, 576) * nu**2
                        + sp.Rational(2185, 448) * nu
                        + (sp.Rational(19423, 576) * nu - sp.Rational(2185, 448))
                        * m2_norm
                    )
                    + S1sq
                    * (
                        -sp.Rational(6011, 576) * nu**2
                        - sp.Rational(8503, 448) * nu
                        + (sp.Rational(8503, 448) - sp.Rational(1177, 576) * nu)
                        * m2_norm
                    )
                )
                + sp.pi
                * v**7
                * (
                    ln_dot_S1**2
                    * (sp.Rational(209, 8) * m2_norm - sp.Rational(209, 8) * nu)
                    + S1sq * (6 * nu - 6 * m2_norm)
                )
                + v**8
                * (
                    ln_dot_S1**2
                    * (
                        -sp.Rational(869429, 6912) * nu**3
                        + sp.Rational(14283281, 48384) * nu**2
                        - sp.Rational(11888267, 48384) * nu
                        + (
                            sp.Rational(11888267, 48384)
                            - sp.Rational(2392243, 6912) * nu**2
                            + sp.Rational(4063301, 16128) * nu
                        )
                        * m2_norm
                    )
                    + S1sq
                    * (
                        sp.Rational(138323, 6912) * nu**3
                        + sp.Rational(711521, 5376) * nu**2
                        + sp.Rational(8207303, 145152) * nu
                        + (
                            sp.Rational(250693, 6912) * nu**2
                            - sp.Rational(812353, 5376) * nu
                            - sp.Rational(8207303, 145152)
                        )
                        * m2_norm
                    )
                )
            )
            # 1 <-> 2 part
            + (
                v**4
                * (
                    ln_dot_S2**2
                    * (sp.Rational(719, 96) * m1_norm - sp.Rational(719, 96) * nu)
                    + S2sq
                    * (sp.Rational(233, 96) * nu - sp.Rational(233, 96) * m1_norm)
                )
                + v**6
                * (
                    ln_dot_S2**2
                    * (
                        sp.Rational(25373, 576) * nu**2
                        + sp.Rational(2185, 448) * nu
                        + (sp.Rational(19423, 576) * nu - sp.Rational(2185, 448))
                        * m1_norm
                    )
                    + S2sq
                    * (
                        -sp.Rational(6011, 576) * nu**2
                        - sp.Rational(8503, 448) * nu
                        + (sp.Rational(8503, 448) - sp.Rational(1177, 576) * nu)
                        * m1_norm
                    )
                )
                + sp.pi
                * v**7
                * (
                    ln_dot_S2**2
                    * (sp.Rational(209, 8) * m1_norm - sp.Rational(209, 8) * nu)
                    + S2sq * (6 * nu - 6 * m1_norm)
                )
                + v**8
                * (
                    ln_dot_S2**2
                    * (
                        -sp.Rational(869429, 6912) * nu**3
                        + sp.Rational(14283281, 48384) * nu**2
                        - sp.Rational(11888267, 48384) * nu
                        + (
                            sp.Rational(11888267, 48384)
                            - sp.Rational(2392243, 6912) * nu**2
                            + sp.Rational(4063301, 16128) * nu
                        )
                        * m1_norm
                    )
                    + S2sq
                    * (
                        sp.Rational(138323, 6912) * nu**3
                        + sp.Rational(711521, 5376) * nu**2
                        + sp.Rational(8207303, 145152) * nu
                        + (
                            sp.Rational(250693, 6912) * nu**2
                            - sp.Rational(812353, 5376) * nu
                            - sp.Rational(8207303, 145152)
                        )
                        * m1_norm
                    )
                )
            )
        )
        # Equation 69a of https://arxiv.org/pdf/2303.18143
        v_dot = (
            sp.Rational(32, 5)
            * nu
            * v**9
            * (v_dot_nospin + v_dot_spin_orbit + v_dot_spin1_spin2 + v_dot_spin_squared)
        )
        # L: Equation 65 of https://arxiv.org/pdf/2303.18143
        # is given as a combination of
        # nospin + spin-orbit + spin1-spin2 + spin-squared
        # Equation 65b of https://arxiv.org/pdf/2303.18143
        L_nospin = ln * (
            1
            + (sp.Rational(1, 6) * nu + sp.Rational(3, 2)) * v**2
            + (
                sp.Rational(1, 24) * nu**2
                - sp.Rational(19, 8) * nu
                + sp.Rational(27, 8)
            )
            * v**4
            + (
                sp.Rational(7, 1296) * nu**3
                + sp.Rational(31, 24) * nu**2
                + (sp.Rational(41, 24) * sp.pi**2 - sp.Rational(6889, 144)) * nu
                + sp.Rational(135, 16)
            )
            * v**6
            + (
                -sp.Rational(55, 31104) * nu**4
                - sp.Rational(215, 1728) * nu**3
                + (sp.Rational(356035, 3456) - sp.Rational(2255, 576) * sp.pi**2)
                * nu**2
                + nu
                * (
                    sp.Rational(98869, 5760)
                    - sp.Rational(128, 3) * EulerGamma
                    - sp.Rational(6455, 1536) * sp.pi**2
                    - sp.Rational(256, 3) * sp.log(2)
                    - sp.Rational(128, 3) * sp.log(v)
                )
                + sp.Rational(2835, 128)
            )
            * v**8
        )
        # Equation 65c of https://arxiv.org/pdf/2303.18143
        L_spin_orbit = (
            # S1 part
            (
                v**3
                * (
                    ln
                    * ln_dot_S1
                    * (-sp.Rational(7, 12) * nu - sp.Rational(7, 4) * m2_norm)
                    + S1 * (-sp.Rational(1, 4) * nu - sp.Rational(3, 4) * m2_norm)
                )
                + v**5
                * (
                    ln
                    * ln_dot_S1
                    * (
                        sp.Rational(11, 144) * nu**2
                        - sp.Rational(55, 16) * nu
                        + (sp.Rational(55, 24) * nu - sp.Rational(33, 16)) * m2_norm
                    )
                    + S1
                    * (
                        sp.Rational(1, 48) * nu**2
                        - sp.Rational(15, 16) * nu
                        + (sp.Rational(5, 8) * nu - sp.Rational(9, 16)) * m2_norm
                    )
                )
                + v**7
                * (
                    ln
                    * ln_dot_S1
                    * (
                        sp.Rational(5, 96) * nu**3
                        + sp.Rational(275, 32) * nu**2
                        - sp.Rational(405, 32) * nu
                        + (
                            -sp.Rational(25, 32) * nu**2
                            + sp.Rational(195, 8) * nu
                            - sp.Rational(135, 32)
                        )
                        * m2_norm
                    )
                    + S1
                    * (
                        sp.Rational(1, 96) * nu**3
                        + sp.Rational(55, 32) * nu**2
                        - sp.Rational(81, 32) * nu
                        + (
                            -sp.Rational(5, 32) * nu**2
                            + sp.Rational(39, 8) * nu
                            - sp.Rational(27, 32)
                        )
                        * m2_norm
                    )
                )
            )
            # 1 <-> 2 part
            + (
                v**3
                * (
                    ln
                    * ln_dot_S2
                    * (-sp.Rational(7, 12) * nu - sp.Rational(7, 4) * m1_norm)
                    + S2 * (-sp.Rational(1, 4) * nu - sp.Rational(3, 4) * m1_norm)
                )
                + v**5
                * (
                    ln
                    * ln_dot_S2
                    * (
                        sp.Rational(11, 144) * nu**2
                        - sp.Rational(55, 16) * nu
                        + (sp.Rational(55, 24) * nu - sp.Rational(33, 16)) * m1_norm
                    )
                    + S2
                    * (
                        sp.Rational(1, 48) * nu**2
                        - sp.Rational(15, 16) * nu
                        + (sp.Rational(5, 8) * nu - sp.Rational(9, 16)) * m1_norm
                    )
                )
                + v**7
                * (
                    ln
                    * ln_dot_S2
                    * (
                        sp.Rational(5, 96) * nu**3
                        + sp.Rational(275, 32) * nu**2
                        - sp.Rational(405, 32) * nu
                        + (
                            -sp.Rational(25, 32) * nu**2
                            + sp.Rational(195, 8) * nu
                            - sp.Rational(135, 32)
                        )
                        * m1_norm
                    )
                    + S2
                    * (
                        sp.Rational(1, 96) * nu**3
                        + sp.Rational(55, 32) * nu**2
                        - sp.Rational(81, 32) * nu
                        + (
                            -sp.Rational(5, 32) * nu**2
                            + sp.Rational(39, 8) * nu
                            - sp.Rational(27, 32)
                        )
                        * m1_norm
                    )
                )
            )
        )
        # Equation 65d of https://arxiv.org/pdf/2303.18143
        L_spin1_spin2 = nu * (
            v**4
            * (
                ln * (2 * ln_dot_S1 * ln_dot_S2 - S1_dot_S2)
                + sp.Rational(1, 2) * ln_dot_S1 * S2
                + sp.Rational(1, 2) * ln_dot_S2 * S1
            )
            + v**6
            * (
                ln
                * (
                    ln_dot_S1
                    * ln_dot_S2
                    * (sp.Rational(13, 36) * nu - sp.Rational(7, 6))
                    + sp.Rational(2, 3) * nu * S1_dot_S2
                )
                + S2 * ln_dot_S1 * (sp.Rational(5, 4) - sp.Rational(7, 24) * nu)
                + S1 * ln_dot_S2 * (sp.Rational(5, 4) - sp.Rational(7, 24) * nu)
            )
            + v**8
            * (
                ln
                * (
                    ln_dot_S1
                    * ln_dot_S2
                    * (
                        -sp.Rational(361, 432) * nu**2
                        + sp.Rational(361, 288) * nu
                        + sp.Rational(15, 4)
                    )
                    + (
                        -sp.Rational(5, 72) * nu**2
                        - sp.Rational(245, 24) * nu
                        - sp.Rational(5, 4)
                    )
                    * S1_dot_S2
                )
                + S2
                * ln_dot_S1
                * (
                    -sp.Rational(223, 288) * nu**2
                    - sp.Rational(349, 64) * nu
                    + sp.Rational(15, 8)
                )
                + S1
                * ln_dot_S2
                * (
                    -sp.Rational(223, 288) * nu**2
                    - sp.Rational(349, 64) * nu
                    + sp.Rational(15, 8)
                )
            )
        )
        # Equation 65e of https://arxiv.org/pdf/2303.18143
        L_spin_squared = (
            # S1 part
            (
                v**4
                * (
                    ln
                    * (
                        ln_dot_S1**2 * (m2_norm - nu)
                        + S1sq * (sp.Rational(1, 2) * nu - sp.Rational(1, 2) * m2_norm)
                    )
                    + ln_dot_S1
                    * S1
                    * (sp.Rational(1, 2) * m2_norm - sp.Rational(1, 2) * nu)
                )
                + v**6
                * (
                    ln
                    * (
                        ln_dot_S1**2
                        * (
                            sp.Rational(121, 72) * nu**2
                            + sp.Rational(35, 8) * nu
                            + (sp.Rational(11, 2) * nu - sp.Rational(35, 8)) * m2_norm
                        )
                        + S1sq
                        * (
                            -sp.Rational(1, 3) * nu**2
                            - nu
                            + (1 - sp.Rational(5, 3) * nu) * m2_norm
                        )
                    )
                    + S1
                    * ln_dot_S1
                    * (
                        sp.Rational(5, 24) * nu**2
                        - sp.Rational(11, 8) * nu
                        + (sp.Rational(11, 8) - sp.Rational(1, 2) * nu) * m2_norm
                    )
                )
                + v**8
                * (
                    ln
                    * ln_dot_S1**2
                    * (
                        -sp.Rational(505, 864) * nu**3
                        + sp.Rational(347, 96) * nu**2
                        + sp.Rational(111, 32) * nu
                        + (
                            -sp.Rational(2833, 288) * nu**2
                            + sp.Rational(199, 16) * nu
                            - sp.Rational(111, 32)
                        )
                        * m2_norm
                    )
                    + ln
                    * S1sq
                    * (
                        sp.Rational(5, 144) * nu**3
                        + sp.Rational(275, 48) * nu**2
                        - sp.Rational(15, 16) * nu
                        + (
                            sp.Rational(295, 144) * nu**2
                            - sp.Rational(455, 48) * nu
                            + sp.Rational(15, 16)
                        )
                        * m2_norm
                    )
                    + S1
                    * ln_dot_S1
                    * (
                        -sp.Rational(235, 288) * nu**3
                        + sp.Rational(563, 96) * nu**2
                        - sp.Rational(21, 32) * nu
                        + (
                            -sp.Rational(113, 32) * nu**2
                            - sp.Rational(7, 3) * nu
                            + sp.Rational(21, 32)
                        )
                        * m2_norm
                    )
                )
            )
            # 1 <-> 2 part
            + (
                v**4
                * (
                    ln
                    * (
                        ln_dot_S2**2 * (m1_norm - nu)
                        + S2sq * (sp.Rational(1, 2) * nu - sp.Rational(1, 2) * m1_norm)
                    )
                    + ln_dot_S2
                    * S2
                    * (sp.Rational(1, 2) * m1_norm - sp.Rational(1, 2) * nu)
                )
                + v**6
                * (
                    ln
                    * (
                        ln_dot_S2**2
                        * (
                            sp.Rational(121, 72) * nu**2
                            + sp.Rational(35, 8) * nu
                            + (sp.Rational(11, 2) * nu - sp.Rational(35, 8)) * m1_norm
                        )
                        + S2sq
                        * (
                            -sp.Rational(1, 3) * nu**2
                            - nu
                            + (1 - sp.Rational(5, 3) * nu) * m1_norm
                        )
                    )
                    + S2
                    * ln_dot_S2
                    * (
                        sp.Rational(5, 24) * nu**2
                        - sp.Rational(11, 8) * nu
                        + (sp.Rational(11, 8) - sp.Rational(1, 2) * nu) * m1_norm
                    )
                )
                + v**8
                * (
                    ln
                    * ln_dot_S2**2
                    * (
                        -sp.Rational(505, 864) * nu**3
                        + sp.Rational(347, 96) * nu**2
                        + sp.Rational(111, 32) * nu
                        + (
                            -sp.Rational(2833, 288) * nu**2
                            + sp.Rational(199, 16) * nu
                            - sp.Rational(111, 32)
                        )
                        * m1_norm
                    )
                    + ln
                    * S2sq
                    * (
                        sp.Rational(5, 144) * nu**3
                        + sp.Rational(275, 48) * nu**2
                        - sp.Rational(15, 16) * nu
                        + (
                            sp.Rational(295, 144) * nu**2
                            - sp.Rational(455, 48) * nu
                            + sp.Rational(15, 16)
                        )
                        * m1_norm
                    )
                    + S2
                    * ln_dot_S2
                    * (
                        -sp.Rational(235, 288) * nu**3
                        + sp.Rational(563, 96) * nu**2
                        - sp.Rational(21, 32) * nu
                        + (
                            -sp.Rational(113, 32) * nu**2
                            - sp.Rational(7, 3) * nu
                            + sp.Rational(21, 32)
                        )
                        * m1_norm
                    )
                )
            )
        )
        # Equation 65a of https://arxiv.org/pdf/2303.18143
        L = (L_nospin + L_spin_orbit + L_spin1_spin2 + L_spin_squared) / v
        # ln_dot: Equation 71 of https://arxiv.org/pdf/2303.18143
        # is given as a combination of
        # spin-orbit + spin1-spin2 + spin-squared
        # Equation 71b of https://arxiv.org/pdf/2303.18143
        ln_dot_spin_orbit = (
            # S1 part
            ln_cross_S1
            * (
                v**6 * (-sp.Rational(1, 2) * nu - sp.Rational(3, 2) * m2_norm)
                + v**8
                * (
                    sp.Rational(1, 4) * nu**2
                    - sp.Rational(9, 4) * nu
                    + (sp.Rational(9, 4) * nu + sp.Rational(9, 4)) * m2_norm
                )
                + v**10
                * (
                    -sp.Rational(1, 48) * nu**3
                    + sp.Rational(81, 16) * nu**2
                    - sp.Rational(27, 16) * nu
                    + (
                        -sp.Rational(21, 16) * nu**2
                        + sp.Rational(63, 16) * nu
                        + sp.Rational(27, 16)
                    )
                    * m2_norm
                )
            )
            # 1 <-> 2 part
            + ln_cross_S2
            * (
                v**6 * (-sp.Rational(1, 2) * nu - sp.Rational(3, 2) * m1_norm)
                + v**8
                * (
                    sp.Rational(1, 4) * nu**2
                    - sp.Rational(9, 4) * nu
                    + (sp.Rational(9, 4) * nu + sp.Rational(9, 4)) * m1_norm
                )
                + v**10
                * (
                    -sp.Rational(1, 48) * nu**3
                    + sp.Rational(81, 16) * nu**2
                    - sp.Rational(27, 16) * nu
                    + (
                        -sp.Rational(21, 16) * nu**2
                        + sp.Rational(63, 16) * nu
                        + sp.Rational(27, 16)
                    )
                    * m1_norm
                )
            )
        )
        # Equation 71c of https://arxiv.org/pdf/2303.18143
        ln_dot_spin1_spin2 = nu * (
            sp.Rational(3, 2)
            * v**7
            * (ln_cross_S1 * ln_dot_S2 + ln_cross_S2 * ln_dot_S1)
            + v**9
            * (
                ln_cross_S1
                * ln_dot_S2
                * (
                    -sp.Rational(5, 4) * nu
                    - sp.Rational(15, 8) * m2_norm
                    - sp.Rational(21, 4)
                )
                + ln_cross_S2
                * ln_dot_S1
                * (
                    -sp.Rational(5, 4) * nu
                    + sp.Rational(15, 8) * m2_norm
                    - sp.Rational(57, 8)
                )
                - sp.Rational(5, 8) * delta * ln * ln_dot_S1_cross_S2
                + sp.Rational(3, 8) * delta * S1_cross_S2
            )
            + v**11
            * (
                ln_cross_S1
                * ln_dot_S2
                * (
                    -sp.Rational(1, 6) * nu**2
                    + sp.Rational(25, 4) * nu
                    + (sp.Rational(71, 32) * nu + sp.Rational(9, 32)) * m2_norm
                    + sp.Rational(15, 16)
                )
                + ln_cross_S2
                * ln_dot_S1
                * (
                    -sp.Rational(1, 6) * nu**2
                    + sp.Rational(271, 32) * nu
                    + (-sp.Rational(71, 32) * nu - sp.Rational(9, 32)) * m2_norm
                    + sp.Rational(39, 32)
                )
                + sp.Rational(1, 96) * delta * (89 * nu - 27) * ln * ln_dot_S1_cross_S2
                - sp.Rational(9, 32) * delta * (2 * nu + 1) * S1_cross_S2
            )
        )
        # Equation 71d of https://arxiv.org/pdf/2303.18143
        ln_dot_spin_squared = (
            # S1 part
            (ln_cross_S1 * ln_dot_S1)
            * (
                v**7 * (sp.Rational(3, 2) * m2_norm - sp.Rational(3, 2) * nu)
                + v**9 * (2 * nu**2 + 9 * nu + (3 * nu - 9) * m2_norm)
                + v**11
                * (
                    -sp.Rational(23, 16) * nu**3
                    - sp.Rational(157, 16) * nu**2
                    - sp.Rational(93, 16) * nu
                    + (-sp.Rational(439, 48) * nu**2 + 4 * nu + sp.Rational(93, 16))
                    * m2_norm
                )
            )
            # 1 <-> 2 part
            + (ln_cross_S2 * ln_dot_S2)
            * (
                v**7 * (sp.Rational(3, 2) * m1_norm - sp.Rational(3, 2) * nu)
                + v**9 * (2 * nu**2 + 9 * nu + (3 * nu - 9) * m1_norm)
                + v**11
                * (
                    -sp.Rational(23, 16) * nu**3
                    - sp.Rational(157, 16) * nu**2
                    - sp.Rational(93, 16) * nu
                    + (-sp.Rational(439, 48) * nu**2 + 4 * nu + sp.Rational(93, 16))
                    * m1_norm
                )
            )
        )
        # Equation 71a of https://arxiv.org/pdf/2303.18143
        ln_dot = ln_dot_spin_orbit + ln_dot_spin1_spin2 + ln_dot_spin_squared
        # Collect all vector components
        self.S1_dot_x = S1_dot[0]
        self.S1_dot_y = S1_dot[1]
        self.S1_dot_z = S1_dot[2]
        self.S2_dot_x = S2_dot[0]
        self.S2_dot_y = S2_dot[1]
        self.S2_dot_z = S2_dot[2]
        self.ln_dot_x = ln_dot[0]
        self.ln_dot_y = ln_dot[1]
        self.ln_dot_z = ln_dot[2]
        self.L_x = L[0]
        self.L_y = L[1]
        self.L_z = L[2]
        # Use chain rule on omega = v**3 to get omega_dot = 3 v**2 v_dot
        self.omega_dot = 3 * v**2 * v_dot


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
        SEOBNRv5_spin_evolution_equations().__dict__,
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
