"""
Construct the co-precessing frame transformations and inertial GW polarizations for SEOBNRv5PHM.

Authors:
Suchindram Dasgupta
sd00113 at mix dot wvu dot edu

Zachariah B. Etienne
zachetie at gmail *dot com

This module implements the transformation from the co-precessing frame to the
inertial observer frame, computing the final GW polarizations (h_+, h_x).
It includes:
1. The J_f-frame triad construction from the remnant angular momentum.
2. The extraction of the inertial Euler angles.
3. The Wigner small-d rotations to project the co-precessing modes into the
   final inertial polarizations.

This implementation is based on Equations 15, 17, and 26-28 of
https://arxiv.org/pdf/2303.18046.pdf (SEOBNRv5PHM paper).

License: BSD 2-Clause
"""

import os
import sys
from typing import List, Optional, Tuple

import sympy as sp

thismodule = __name__


def wigner_d_small(l: int, m1: int, m2: int, beta: sp.Expr) -> sp.Expr:
    """
    Computes the Wigner small-d matrix element d^l_{m1, m2}(beta) using the
    explicit polynomial sum. This ensures pure real algebraic SymPy expressions,
    avoiding complex branches in generated C-code.

    Inputs:
        l (int): The l index of the d-matrix.
        m1 (int): The first m index.
        m2 (int): The second m index.
        beta (sp.Expr): The angle (SymPy expression) to evaluate.

    Outputs:
        sp.Expr: The exact symbolic expression for the d-matrix element.

    Raises:
        ValueError: If l < 0, or if |m1| > l or |m2| > l.

    >>> import sympy as sp
    >>> beta = sp.Symbol("beta", real=True)
    >>> wigner_d_small(2, 2, 2, beta)
    cos(beta/2)**4
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

    d_val = sp.sympify(0)
    k_min = max(0, m2 - m1)
    k_max = min(l - m1, l + m2)

    for k in range(k_min, k_max + 1):
        num = ((-1) ** k) * sp.sqrt(
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
            * (sp.cos(beta / 2) ** (2 * l - 2 * k - m1 + m2))
            * (sp.sin(beta / 2) ** (2 * k + m1 - m2))
        )
        d_val += term

    return sp.simplify(d_val)


class SEOBNRv5_Coprecessing_Rotations:
    """
    Compute the SEOBNRv5 co-precessing frame transformations and polarizations.

    This class sets up the symbolic expressions for converting the co-precessing
    GW modes (h^P_{lm}) into the inertial frame polarizations (h_+, h_x), chaining
    the J_f-frame, the co-precessing Euler angles, and the Wigner rotations.
    """

    def __init__(self, modes: Optional[List[Tuple[int, int]]] = None) -> None:
        """
        Compute the symbolic expressions for the SEOBNRv5 co-precessing rotations.

        Inputs:
        -------
        - modes: A list of (l, m) tuples specifying which modes to include in the
          polarization sum. Only m >= 0 modes should be passed; the -m symmetry
          is handled automatically.

        Outputs (as class attributes):
        ------------------------------
        - self.R_Eq15_*: Matrix elements for the J_f-frame transformation (Eq. 15).
        - self.alpha_PI, self.beta_PI, self.gamma_PI: The combined P -> I Euler angles.
        - self.h_plus_I, self.h_cross_I: The final inertial GW polarizations (Eq. 28).

        Raises:
            ValueError: If the `modes` list is empty or invalid.

        >>> rot = SEOBNRv5_Coprecessing_Rotations(modes=[(2, 2)])
        >>> type(rot.h_plus_I)
        <class 'sympy.core.add.Add'>
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

        J_norm = sp.sqrt(sp.Max(Jfx**2 + Jfy**2 + Jfz**2, 1e-30))
        e3_J_x = Jfx / J_norm
        e3_J_y = Jfy / J_norm
        e3_J_z = Jfz / J_norm

        # Continuous blending for Gram-Schmidt (prevents MCMC gradient discontinuities).
        # Fallback to the y-axis (0, 1, 0) is required if J_f is nearly parallel
        # to the x-axis, which would cause the orthogonal projection to collapse
        # to zero (NaNs) during normalization.
        abs_e3_x = sp.Abs(e3_J_x)
        weight_x = sp.Piecewise(
            (sp.sympify(0), 1 - abs_e3_x <= 1e-5),
            (sp.sympify(1), 1 - abs_e3_x >= 1e-4),
            ((1 - abs_e3_x - 1e-5) / 9e-5, True),
        )
        weight_y = 1 - weight_x

        v_x_x = 1 - e3_J_x**2
        v_x_y = -e3_J_x * e3_J_y
        v_x_z = -e3_J_x * e3_J_z

        v_y_x = -e3_J_y * e3_J_x
        v_y_y = 1 - e3_J_y**2
        v_y_z = -e3_J_y * e3_J_z

        e1_J_unnorm_x = weight_x * v_x_x + weight_y * v_y_x
        e1_J_unnorm_y = weight_x * v_x_y + weight_y * v_y_y
        e1_J_unnorm_z = weight_x * v_x_z + weight_y * v_y_z

        norm_e1 = sp.sqrt(
            sp.Max(e1_J_unnorm_x**2 + e1_J_unnorm_y**2 + e1_J_unnorm_z**2, 1e-30)
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
                setattr(self, f"R_Eq15_{i}{j}", sp.simplify(R_Eq15[i, j]))

        # Step 2: Combine matrices for the Inertial Euler Angles (Eqs. 17, 26)

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
        R_obs = R_obs_z * R_obs_y

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
        R_PJ = R_z_alpha * R_y_beta * R_z_gamma

        M_PI = R_obs * R_JI * R_PJ

        M33 = M_PI[2, 2]
        M13 = M_PI[0, 2]
        M23 = M_PI[1, 2]
        M31 = M_PI[2, 0]
        M32 = M_PI[2, 1]
        M11 = M_PI[0, 0]
        M21 = M_PI[1, 0]

        M33_clamped = sp.Max(-1, sp.Min(1, M33))
        TOL = 1e-14

        self.beta_PI = sp.Piecewise(
            (0, M33 > 1 - TOL), (sp.pi, M33 < -1 + TOL), (sp.acos(M33_clamped), True)
        )
        self.alpha_PI = sp.Piecewise(
            (sp.atan2(M21, M11), M33 > 1 - TOL),
            (sp.atan2(-M21, -M11), M33 < -1 + TOL),
            (sp.atan2(M23, M13), True),
        )
        self.gamma_PI = sp.Piecewise(
            (0, M33 > 1 - TOL), (0, M33 < -1 + TOL), (sp.atan2(M32, -M31), True)
        )

        # Step 3: Project Co-Precessing modes to Inertial Polarizations (Eq. 28)

        h_plus_expr = sp.sympify(0)
        h_cross_expr = sp.sympify(0)

        for l, m in modes:
            hR = sp.Symbol(f"hP_l{l}m{m}_Re", real=True)
            hI = sp.Symbol(f"hP_l{l}m{m}_Im", real=True)

            pref = sp.sqrt((2 * l + 1) / (4 * sp.pi))

            d_mat_pos = wigner_d_small(l, m, 2, self.beta_PI)
            A_pos = pref * d_mat_pos
            phase_pos = 2 * self.alpha_PI + m * self.gamma_PI

            Re_Cp = A_pos * (hR * sp.cos(phase_pos) - hI * sp.sin(phase_pos))
            Im_Cp = A_pos * (hR * sp.sin(phase_pos) + hI * sp.cos(phase_pos))

            if m == 0:
                h_plus_expr += Re_Cp
                h_cross_expr -= Im_Cp
                continue

            hR_neg = ((-1) ** l) * hR
            hI_neg = ((-1) ** l) * (-hI)

            d_mat_neg = wigner_d_small(l, -m, 2, self.beta_PI)
            A_neg = pref * d_mat_neg
            phase_neg = 2 * self.alpha_PI - m * self.gamma_PI

            Re_Cm = A_neg * (hR_neg * sp.cos(phase_neg) - hI_neg * sp.sin(phase_neg))
            Im_Cm = A_neg * (hR_neg * sp.sin(phase_neg) + hI_neg * sp.cos(phase_neg))

            h_plus_expr += Re_Cp + Re_Cm
            h_cross_expr -= Im_Cp + Im_Cm

        self.h_plus_I = h_plus_expr
        self.h_cross_I = h_cross_expr


if __name__ == "__main__":
    import doctest

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
