"""
Construct the co-precessing frame transformations and inertial GW polarizations for SEOBNRv5PHM.

Authors:
Suchindram Dasgupta
sd00113 at mix dot wvu dot edu

Zachariah B. Etienne
zachetie at gmail *dot com

This module implements the transformation from the co-precessing frame to the
inertial observer frame, computing the final GW polarizations (h_+, h_x).
It encapsulates:
1. The J_f-frame triad construction from the remnant angular momentum.
2. The extraction of the inertial Euler angles.
3. The Wigner small-d rotations to project the co-precessing modes into the
   final inertial polarizations.

This implementation is based on Equations 15, 17, and 26-28 of
https://arxiv.org/pdf/2303.18046.pdf (SEOBNRv5PHM paper).

License: BSD 2-Clause
"""

from functools import lru_cache
from typing import List, Optional, Tuple, cast

import sympy as sp

thismodule = __name__


class SEOBNRv5_Coprecessing_Rotations:
    """
    Compute the SEOBNRv5 co-precessing frame transformations and polarizations.

    This class sets up the symbolic expressions for converting the co-precessing
    GW modes (h^P_{lm}) into the inertial frame polarizations (h_+, h_x), chaining
    the J_f-frame, the co-precessing Euler angles, and the Wigner rotations.
    """

    _WIGNER_BETA = sp.Symbol("_wigner_beta", real=True)

    @staticmethod
    @lru_cache(maxsize=None)
    def wigner_d_small_template(l: int, m1: int, m2: int) -> sp.Expr:
        """
        Return a cached Wigner small-d template in terms of `_WIGNER_BETA`.

        This optimization prevents SymPy from recomputing factorials and trigonometric
        expansions repeatedly, cutting generation time to milliseconds.

        :param l: Degree index of the Wigner d-matrix element.
        :param m1: First magnetic index.
        :param m2: Second magnetic index.
        :return: Symbolic template expression d^l_{m1,m2}(_WIGNER_BETA).
        """
        d_val = sp.sympify(0)
        k_min = max(0, m2 - m1)
        k_max = min(l - m1, l + m2)

        for k in range(k_min, k_max + 1):
            num = ((-1) ** (k + m1 - m2)) * sp.sqrt(
                sp.factorial(l + m1)
                * sp.factorial(l - m1)
                * sp.factorial(l + m2)
                * sp.factorial(l - m2)
            )
            den = sp.sympify(
                sp.factorial(l - m1 - k)
                * sp.factorial(l + m2 - k)
                * sp.factorial(k)
                * sp.factorial(k + m1 - m2)
            )

            term = (
                (sp.S(num) / den)
                * (
                    sp.cos(SEOBNRv5_Coprecessing_Rotations._WIGNER_BETA / 2)
                    ** (2 * l - 2 * k - m1 + m2)
                )
                * (
                    sp.sin(SEOBNRv5_Coprecessing_Rotations._WIGNER_BETA / 2)
                    ** (2 * k + m1 - m2)
                )
            )
            d_val += term

        return d_val

    @classmethod
    def wigner_d_small(cls, l: int, m1: int, m2: int, beta: sp.Expr) -> sp.Expr:
        """
        Compute the Wigner small-d matrix element d^l_{m1,m2}(beta).

        :param l: Degree index of the Wigner d-matrix element.
        :param m1: First magnetic index.
        :param m2: Second magnetic index.
        :param beta: Euler-angle argument used in the rotation.
        :return: Symbolic value of d^l_{m1,m2}(beta).
        :raises ValueError: If ``l < 0`` or if ``|m1| > l`` or ``|m2| > l``.
        """
        if l < 0:
            raise ValueError(
                f"Wigner-d matrix index l must be non-negative. Received l={l}"
            )
        if abs(m1) > l or abs(m2) > l:
            raise ValueError(
                f"Wigner-d matrix indices m1 and m2 must satisfy |m| <= l. "
                f"Received m1={m1}, m2={m2}, l={l}"
            )

        return cast(
            sp.Expr,
            cls.wigner_d_small_template(l, m1, m2).xreplace({cls._WIGNER_BETA: beta}),
        )

    def _polarizations_from_angles(
        self,
        alpha: sp.Expr,
        beta: sp.Expr,
        gamma: sp.Expr,
        modes: List[Tuple[int, int]],
    ) -> Tuple[sp.Expr, sp.Expr]:
        """
        Generate polarizations for a specific Euler-angle branch.

        :param alpha: First Euler angle for the inertial-frame rotation.
        :param beta: Second Euler angle for the inertial-frame rotation.
        :param gamma: Third Euler angle for the inertial-frame rotation.
        :param modes: List of ``(l, m)`` co-precessing modes with ``m >= 0``.
        :return: Tuple ``(h_plus_expr, h_cross_expr)`` of symbolic polarizations.
        """
        h_plus_expr = sp.sympify(0)
        h_cross_expr = sp.sympify(0)

        for l, m in modes:
            hR = sp.Symbol(f"hP_l{l}m{m}_Re", real=True)
            hI = sp.Symbol(f"hP_l{l}m{m}_Im", real=True)

            pref = sp.sqrt(sp.Rational(2 * l + 1, 4) / sp.pi)

            d_mat_pos = self.wigner_d_small(l, m, 2, beta)
            A_pos = pref * d_mat_pos
            phase_pos = 2 * alpha + m * gamma

            Re_Cp = A_pos * (hR * sp.cos(phase_pos) - hI * sp.sin(phase_pos))
            Im_Cp = A_pos * (hR * sp.sin(phase_pos) + hI * sp.cos(phase_pos))

            if m == 0:
                h_plus_expr += Re_Cp
                h_cross_expr -= Im_Cp
                continue

            hR_neg = ((-1) ** l) * hR
            hI_neg = ((-1) ** l) * (-hI)

            d_mat_neg = self.wigner_d_small(l, -m, 2, beta)
            A_neg = pref * d_mat_neg
            phase_neg = 2 * alpha - m * gamma

            Re_Cm = A_neg * (hR_neg * sp.cos(phase_neg) - hI_neg * sp.sin(phase_neg))
            Im_Cm = A_neg * (hR_neg * sp.sin(phase_neg) + hI_neg * sp.cos(phase_neg))

            h_plus_expr += Re_Cp + Re_Cm
            h_cross_expr -= Im_Cp + Im_Cm

        return h_plus_expr, h_cross_expr

    def __init__(self, modes: Optional[List[Tuple[int, int]]] = None) -> None:
        """
        Compute symbolic expressions for SEOBNRv5 co-precessing rotations.

        :param modes: Optional list of ``(l, m)`` modes with nonnegative ``m``.
            If ``None``, use the default mode set.
        :raises ValueError: If ``modes`` is empty or contains any mode with ``m < 0``.
        """
        if modes is None:
            modes = [(2, 2), (2, 1), (3, 3), (3, 2), (4, 4), (4, 3), (5, 5)]

        if not modes:
            raise ValueError(
                "The modes list cannot be empty. Must provide at least one (l, m) mode."
            )

        for l, m in modes:
            if m < 0:
                raise ValueError(
                    f"Only m >= 0 modes should be specified. Received m={m} for l={l}. "
                    "Negative m symmetries are handled internally."
                )

        # Step 1: Inputs for J_f frame construction (Eq. 15)
        Jfx, Jfy, Jfz = sp.symbols("J_f_x J_f_y J_f_z", real=True)

        J_norm = sp.sqrt(sp.Max(Jfx**2 + Jfy**2 + Jfz**2, sp.sympify("1e-30")))
        e3_J_x = Jfx / J_norm
        e3_J_y = Jfy / J_norm
        e3_J_z = Jfz / J_norm

        # Continuous blending for Gram-Schmidt.
        abs_e3_x = sp.Abs(e3_J_x)
        weight_x = sp.Max(
            sp.sympify(0),
            sp.Min(
                sp.sympify(1),
                (sp.sympify(1) - abs_e3_x - sp.sympify("1e-5")) / sp.sympify("9e-5"),
            ),
        )
        weight_y = sp.sympify(1) - weight_x

        # Preserve J-frame continuity by flipping the fallback y-projection sign.
        # We use an algebraic sign formulation (x / |x|) from throwing "undefined reference to 'sign'".
        fallback_sign = -e3_J_x / sp.Max(abs_e3_x, sp.sympify("1e-30"))

        v_x_x = 1 - e3_J_x**2
        v_x_y = -e3_J_x * e3_J_y
        v_x_z = -e3_J_x * e3_J_z

        v_y_x = fallback_sign * (-e3_J_y * e3_J_x)
        v_y_y = fallback_sign * (1 - e3_J_y**2)
        v_y_z = fallback_sign * (-e3_J_y * e3_J_z)

        # Normalize projections individually before blending
        norm_vx = sp.sqrt(sp.Max(v_x_x**2 + v_x_y**2 + v_x_z**2, sp.sympify("1e-30")))
        norm_vy = sp.sqrt(sp.Max(v_y_x**2 + v_y_y**2 + v_y_z**2, sp.sympify("1e-30")))

        e1_J_unnorm_x = weight_x * v_x_x / norm_vx + weight_y * v_y_x / norm_vy
        e1_J_unnorm_y = weight_x * v_x_y / norm_vx + weight_y * v_y_y / norm_vy
        e1_J_unnorm_z = weight_x * v_x_z / norm_vx + weight_y * v_y_z / norm_vy

        norm_e1 = sp.sqrt(
            sp.Max(
                e1_J_unnorm_x**2 + e1_J_unnorm_y**2 + e1_J_unnorm_z**2,
                sp.sympify("1e-30"),
            )
        )

        e1_J_x = e1_J_unnorm_x / norm_e1
        e1_J_y = e1_J_unnorm_y / norm_e1
        e1_J_z = e1_J_unnorm_z / norm_e1

        e2_J_x = e3_J_y * e1_J_z - e3_J_z * e1_J_y
        e2_J_y = e3_J_z * e1_J_x - e3_J_x * e1_J_z
        e2_J_z = e3_J_x * e1_J_y - e3_J_y * e1_J_x

        R_Eq15 = sp.Matrix(
            [
                [e1_J_x, e1_J_y, e1_J_z],
                [e2_J_x, e2_J_y, e2_J_z],
                [e3_J_x, e3_J_y, e3_J_z],
            ]
        )

        for i in range(3):
            for j in range(3):
                setattr(self, f"R_Eq15_{i}{j}", R_Eq15[i, j])

        # Step 2: Combine matrices for the inertial Euler angles (Eqs. 17 and 26).
        alpha_JP, beta_JP, gamma_JP = sp.symbols("alpha_JP beta_JP gamma_JP", real=True)
        iota, phi0 = sp.symbols("iota varphi_0", real=True)

        R_obs_z = sp.Matrix(
            [
                [sp.cos(phi0), -sp.sin(phi0), sp.sympify(0)],
                [sp.sin(phi0), sp.cos(phi0), sp.sympify(0)],
                [sp.sympify(0), sp.sympify(0), sp.sympify(1)],
            ]
        )
        R_obs_y = sp.Matrix(
            [
                [sp.cos(iota), sp.sympify(0), sp.sin(iota)],
                [sp.sympify(0), sp.sympify(1), sp.sympify(0)],
                [-sp.sin(iota), sp.sympify(0), sp.cos(iota)],
            ]
        )
        R_obs = R_obs_y * R_obs_z
        R_JI = R_Eq15.T

        R_z_alpha = sp.Matrix(
            [
                [sp.cos(alpha_JP), -sp.sin(alpha_JP), sp.sympify(0)],
                [sp.sin(alpha_JP), sp.cos(alpha_JP), sp.sympify(0)],
                [sp.sympify(0), sp.sympify(0), sp.sympify(1)],
            ]
        )
        R_y_beta = sp.Matrix(
            [
                [sp.cos(beta_JP), sp.sympify(0), sp.sin(beta_JP)],
                [sp.sympify(0), sp.sympify(1), sp.sympify(0)],
                [-sp.sin(beta_JP), sp.sympify(0), sp.cos(beta_JP)],
            ]
        )
        R_z_gamma = sp.Matrix(
            [
                [sp.cos(gamma_JP), -sp.sin(gamma_JP), sp.sympify(0)],
                [sp.sin(gamma_JP), sp.cos(gamma_JP), sp.sympify(0)],
                [sp.sympify(0), sp.sympify(0), sp.sympify(1)],
            ]
        )

        R_PJ = R_z_gamma.T * R_y_beta.T * R_z_alpha.T
        M_PI = R_obs * R_JI * R_PJ

        # Store M33 so the C-infrastructure can generate the if/else pole branches
        self.M33 = M_PI[2, 2]
        M13 = M_PI[0, 2]
        M23 = M_PI[1, 2]
        M31 = M_PI[2, 0]
        M32 = M_PI[2, 1]
        M11 = M_PI[0, 0]
        M21 = M_PI[1, 0]

        M33_clamped = sp.Max(sp.sympify(-1), sp.Min(sp.sympify(1), self.M33))
        self.beta_PI = sp.acos(M33_clamped)

        # Output explicit branches to bypass Piecewise Boolean CSE bugs
        self.alpha_PI_generic = sp.atan2(M23, M13)
        self.gamma_PI_generic = sp.atan2(M32, -M31)

        self.alpha_PI_pole_pos = sp.atan2(M21, M11)
        self.gamma_PI_pole_pos = sp.sympify(0)

        self.alpha_PI_pole_neg = sp.atan2(-M21, -M11)
        self.gamma_PI_pole_neg = sp.sympify(0)

        # Step 3: Project co-precessing modes to inertial polarizations (Eq. 28).
        self.h_plus_I_generic, self.h_cross_I_generic = self._polarizations_from_angles(
            self.alpha_PI_generic, self.beta_PI, self.gamma_PI_generic, modes
        )

        self.h_plus_I_pole_pos, self.h_cross_I_pole_pos = (
            self._polarizations_from_angles(
                self.alpha_PI_pole_pos, self.beta_PI, self.gamma_PI_pole_pos, modes
            )
        )

        self.h_plus_I_pole_neg, self.h_cross_I_pole_neg = (
            self._polarizations_from_angles(
                self.alpha_PI_pole_neg, self.beta_PI, self.gamma_PI_pole_neg, modes
            )
        )

        # Aliases for API compatibility and docstring Doctests
        self.alpha_PI = self.alpha_PI_generic
        self.gamma_PI = self.gamma_PI_generic
        self.h_plus_I = self.h_plus_I_generic
        self.h_cross_I = self.h_cross_I_generic


if __name__ == "__main__":
    import doctest
    import os
    import sys

    project_root = os.path.join(os.path.dirname(__file__), "..", "..", "..")
    sys.path.append(project_root)

    try:
        import nrpy.validate_expressions.validate_expressions as ve
    except ImportError as exc:
        raise ImportError(
            "Could not import NRPy validation module. Please ensure "
            "the script is in the correct directory structure "
            "and __init__.py files are present."
        ) from exc

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    print("Running NRPy expression validation...")

    rotations_class = SEOBNRv5_Coprecessing_Rotations()

    results_dict = ve.process_dictionary_of_expressions(
        rotations_class.__dict__,
        fixed_mpfs_for_free_symbols=True,
    )

    module_name = os.path.splitext(os.path.basename(__file__))[0]

    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        module_name,
        results_dict,
    )
