"""
Compute GRMHD signal-speed bounds for a reconstructed face.

The approximate MHD dispersion relation is Duez et al., Phys. Rev. D 72,
024028 (2005), Eqs. (49)-(50), https://arxiv.org/abs/astro-ph/0503420v2.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple

import sympy as sp

from nrpy.equations.grhd.characteristic_speeds import (
    find_cmax_cmin as find_grhd_cmax_cmin,
)
from nrpy.equations.grmhd.GRMHD_equations import compute_smallb2, compute_smallb4U


def find_cmax_cmin(
    flux_dirn: int,
    gamma_faceDD: List[List[sp.Expr]],
    beta_faceU: List[sp.Expr],
    alpha_face: sp.Expr,
    u4U_r: List[sp.Expr],
    u4U_l: List[sp.Expr],
    BmagU_r: List[sp.Expr],
    BmagU_l: List[sp.Expr],
    rho_b_r: sp.Expr,
    rho_b_l: sp.Expr,
    h_r: sp.Expr,
    h_l: sp.Expr,
    cs2_r: sp.Expr,
    cs2_l: sp.Expr,
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Compute the nonnegative HLL speed bounds from both GRMHD face states.

    The squared Alfvén speed v_A^2 and v_0^2 are defined beside Duez et al. (2005),
    Eq. (50). The bounds use the derived roots in the GRHD
    characteristic-speed helper and Duez et al., Eq. (48). All face inputs,
    including BmagU_r and BmagU_l, must use the basis of gamma_faceDD; for
    reference-metric evolutions this is the rescaled basis, and the returned
    bounds are rescaled speeds.

    :param flux_dirn: Spatial flux direction, 0 through 2.
    :param gamma_faceDD: Covariant face metric.
    :param beta_faceU: Face shift vector.
    :param alpha_face: Face lapse.
    :param u4U_r: Right fluid four-velocity.
    :param u4U_l: Left fluid four-velocity.
    :param BmagU_r: Right Eulerian magnetic field scaled by 1/sqrt(4 pi), in the basis of gamma_faceDD.
    :param BmagU_l: Left Eulerian magnetic field scaled by 1/sqrt(4 pi), in the basis of gamma_faceDD.
    :param rho_b_r: Right baryon density.
    :param rho_b_l: Left baryon density.
    :param h_r: Right specific enthalpy.
    :param h_l: Left specific enthalpy.
    :param cs2_r: Right sound speed squared.
    :param cs2_l: Left sound speed squared.
    :return: Left-going and right-going nonnegative wave-speed bounds.
    """
    smallb4U_r = compute_smallb4U(gamma_faceDD, beta_faceU, alpha_face, u4U_r, BmagU_r)
    smallb4U_l = compute_smallb4U(gamma_faceDD, beta_faceU, alpha_face, u4U_l, BmagU_l)
    smallb2_r = compute_smallb2(gamma_faceDD, beta_faceU, alpha_face, smallb4U_r)
    smallb2_l = compute_smallb2(gamma_faceDD, beta_faceU, alpha_face, smallb4U_l)
    vA2_r = smallb2_r / (rho_b_r * h_r + smallb2_r)
    vA2_l = smallb2_l / (rho_b_l * h_l + smallb2_l)
    v02_r = vA2_r + cs2_r * (1 - vA2_r)
    v02_l = vA2_l + cs2_l * (1 - vA2_l)
    return find_grhd_cmax_cmin(
        flux_dirn,
        gamma_faceDD,
        beta_faceU,
        alpha_face,
        u4U_r,
        u4U_l,
        v02_r,
        v02_l,
    )


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.indexedexp as ixp
    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Step 1: Compare B=0 bounds with GRHD for arbitrary face states in every direction.
    gamma_faceDD_sym = ixp.declarerank2("gamma_faceDD", symmetry="sym01", dimension=3)
    beta_faceU_sym = ixp.declarerank1("beta_faceU", dimension=3)
    alpha_face_sym = sp.symbols("alpha_face", real=True)
    u4U_r_sym = ixp.declarerank1("u4U_r", dimension=4)
    u4U_l_sym = ixp.declarerank1("u4U_l", dimension=4)
    rho_b_r_sym, rho_b_l_sym, h_r_sym, h_l_sym, cs2_r_sym, cs2_l_sym = sp.symbols(
        "rho_b_r rho_b_l h_r h_l cs2_r cs2_l", real=True
    )
    zero_BmagU = ixp.zerorank1(dimension=3)
    for zero_field_dirn in range(3):
        grmhd_zero_speeds = find_cmax_cmin(
            zero_field_dirn,
            gamma_faceDD_sym,
            beta_faceU_sym,
            alpha_face_sym,
            u4U_r_sym,
            u4U_l_sym,
            zero_BmagU,
            zero_BmagU,
            rho_b_r_sym,
            rho_b_l_sym,
            h_r_sym,
            h_l_sym,
            cs2_r_sym,
            cs2_l_sym,
        )
        grhd_speeds = find_grhd_cmax_cmin(
            zero_field_dirn,
            gamma_faceDD_sym,
            beta_faceU_sym,
            alpha_face_sym,
            u4U_r_sym,
            u4U_l_sym,
            cs2_r_sym,
            cs2_l_sym,
        )
        if grmhd_zero_speeds != grhd_speeds:
            raise AssertionError(
                f"B=0 speed bounds differ from GRHD in flux direction {zero_field_dirn}"
            )

    # Step 2: Compare magnetic wave speeds with trusted values. The right state
    #         moves in +y and sets cmax; the left state moves in -y and sets cmin.
    gamma_faceDD_test = sp.diag(1, 4, 9).tolist()
    beta_faceU_test = [sp.Rational(1, 10), sp.Rational(1, 20), -sp.Rational(1, 30)]
    alpha_face_test = sp.Rational(5, 4)
    u4U_r_test = [
        sp.sympify(1),
        -sp.Rational(1, 10),
        sp.Rational(13, 40),
        sp.Rational(1, 30),
    ]
    u4U_l_test = [
        sp.Rational(13, 15),
        -sp.Rational(13, 150),
        -sp.Rational(151, 600),
        sp.Rational(13, 450),
    ]
    rho_b_r_test = sp.sympify(1)
    rho_b_l_test = sp.Rational(4, 5)
    h_r_test = sp.Rational(3, 2)
    h_l_test = sp.Rational(7, 5)
    cs2_r_test = sp.Rational(1, 4)
    cs2_l_test = sp.Rational(1, 5)
    BmagU_r_test = [sp.symbols("Bmag_r", real=True), sp.Rational(1, 3), 0]
    BmagU_l_test = [sp.symbols("Bmag_l", real=True), 0, sp.Rational(1, 4)]
    cmin, cmax = find_cmax_cmin(
        1,
        gamma_faceDD_test,
        beta_faceU_test,
        alpha_face_test,
        u4U_r_test,
        u4U_l_test,
        BmagU_r_test,
        BmagU_l_test,
        rho_b_r_test,
        rho_b_l_test,
        h_r_test,
        h_l_test,
        cs2_r_test,
        cs2_l_test,
    )
    expressions = {
        "cmin": cmin.subs(sp.Function("nrpyAbs"), sp.Abs),
        "cmax": cmax.subs(sp.Function("nrpyAbs"), sp.Abs),
    }
    sampled_results = ve.process_dictionary_of_expressions(
        expressions, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__), os.getcwd(), "characteristic_speeds", sampled_results
    )
