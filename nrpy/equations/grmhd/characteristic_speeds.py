"""
Construct GRMHD characteristic speeds at left and right interfaces of grid cells.
These are to be fed into an approximate Riemann solver.

The approximate MHD dispersion relation is Duez et al., Phys. Rev. D 72,
024028 (2005), Eqs. (49)-(50), https://arxiv.org/abs/astro-ph/0503420v2.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Step 1.a: import all needed modules from NRPy:
from typing import Dict, List, Tuple, cast

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends

from nrpy.equations.grhd.characteristic_speeds import (
    find_cmax_cmin as find_grhd_cmax_cmin,
)
from nrpy.equations.grmhd.GRMHD_equations import compute_smallb2, compute_smallb4U


def compute_v02(
    gamma_faceDD: List[List[sp.Expr]],
    beta_faceU: List[sp.Expr],
    alpha_face: sp.Expr,
    u4U: List[sp.Expr],
    BmagU: List[sp.Expr],
    rho_b: sp.Expr,
    h: sp.Expr,
    cs2: sp.Expr,
) -> sp.Expr:
    """
    Compute the squared fast magnetosonic speed estimate v_0^2 for one face state.

    The squared Alfvén speed v_A^2 = b^2/(rho_b h + b^2) and
    v_0^2 = v_A^2 + c_s^2 (1 - v_A^2) are defined beside Duez et al. (2005),
    Eq. (50). All inputs, including BmagU, must use the basis of gamma_faceDD;
    for reference-metric evolutions this is the rescaled basis.

    :param gamma_faceDD: spatial metric at the cell interface.
    :param beta_faceU: shift vector at cell interface
    :param alpha_face: lapse function at cell interface
    :param u4U: Four-velocity u^mu in the reconstructed face state.
    :param BmagU: Eulerian magnetic field scaled by 1/sqrt(4 pi), in the basis of gamma_faceDD.
    :param rho_b: Baryon density in the reconstructed face state.
    :param h: Specific enthalpy in the reconstructed face state.
    :param cs2: Sound speed squared in the reconstructed face state.

    :return: Squared fluid-frame signal-speed estimate v_0^2.
    """
    smallb4U = compute_smallb4U(gamma_faceDD, beta_faceU, alpha_face, u4U, BmagU)
    smallb2 = compute_smallb2(gamma_faceDD, beta_faceU, alpha_face, smallb4U)
    vA2 = smallb2 / (rho_b * h + smallb2)
    return vA2 + cs2 * (1 - vA2)


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
    Compute nonnegative HLL speed bounds from two reconstructed GRMHD face states.

    The max/min prescription is Duez et al. (2005), Eq. (48) and the
    definitions immediately preceding it. All face inputs, including
    BmagU_r and BmagU_l, must use the basis of gamma_faceDD; for
    reference-metric evolutions this is the rescaled basis, and the returned
    bounds are rescaled speeds.

    :param flux_dirn: Direction for flux calculation.
    :param gamma_faceDD: spatial metric at the cell interface.
    :param beta_faceU: shift vector at cell interface
    :param alpha_face: lapse function at cell interface
    :param u4U_r: Four-velocity u^mu in the right reconstructed face state.
    :param u4U_l: Four-velocity u^mu in the left reconstructed face state.
    :param BmagU_r: Right Eulerian magnetic field scaled by 1/sqrt(4 pi), in the basis of gamma_faceDD.
    :param BmagU_l: Left Eulerian magnetic field scaled by 1/sqrt(4 pi), in the basis of gamma_faceDD.
    :param rho_b_r: Baryon density in the right reconstructed face state.
    :param rho_b_l: Baryon density in the left reconstructed face state.
    :param h_r: Specific enthalpy in the right reconstructed face state.
    :param h_l: Specific enthalpy in the left reconstructed face state.
    :param cs2_r: Sound speed squared in the right reconstructed face state.
    :param cs2_l: Sound speed squared in the left reconstructed face state.

    :return: Nonnegative HLL bounds (cmin, cmax).

    The face metric is shared by the right and left reconstructed states.

    """
    # First, we need to find the signal-speed estimates on each face
    v02_r = compute_v02(
        gamma_faceDD, beta_faceU, alpha_face, u4U_r, BmagU_r, rho_b_r, h_r, cs2_r
    )
    v02_l = compute_v02(
        gamma_faceDD, beta_faceU, alpha_face, u4U_l, BmagU_l, rho_b_l, h_l, cs2_l
    )

    return find_grhd_cmax_cmin(
        flux_dirn, gamma_faceDD, beta_faceU, alpha_face, u4U_r, u4U_l, v02_r, v02_l
    )


def _nrpyAbs_to_Abs(expr: sp.Expr) -> sp.Expr:
    """
    Replace nrpyAbs with unevaluated sp.Abs so trusted values can be computed.

    The traversal is memoized by object, so shared subexpressions are rebuilt
    once, and changed nodes are rebuilt with evaluate=False to avoid SymPy's
    costly absolute-value and canonicalization logic.

    :param expr: Expression containing nrpyAbs calls.
    :return: The same expression with sp.Abs in place of nrpyAbs.
    """
    nrpyAbs = sp.Function("nrpyAbs")
    rebuilt: Dict[int, sp.Basic] = {}
    stack: List[Tuple[sp.Basic, bool]] = [(expr, False)]
    while stack:
        node, expanded = stack.pop()
        if id(node) in rebuilt:
            continue
        if not node.args:
            rebuilt[id(node)] = node
            continue
        if not expanded:
            stack.append((node, True))
            stack.extend((arg, False) for arg in node.args if id(arg) not in rebuilt)
            continue
        args = [rebuilt[id(arg)] for arg in node.args]
        if node.func == nrpyAbs:
            rebuilt[id(node)] = sp.Abs(args[0], evaluate=False)
        elif all(new is old for new, old in zip(args, node.args)):
            rebuilt[id(node)] = node
        else:
            rebuilt[id(node)] = node.func(*args, evaluate=False)
    return cast(sp.Expr, rebuilt[id(expr)])


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

    exprs_dict: Dict[str, sp.Expr] = {}

    import nrpy.indexedexp as ixp  # NRPy: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
    import nrpy.reference_metric as refmetric  # NRPy: Reference metric support
    from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
    from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
    from nrpy.equations.general_relativity.g4munu_conversions import ADM_to_g4UU
    from nrpy.equations.grhd.characteristic_speeds import find_cp_cm

    rfm = refmetric.reference_metric["Cartesian"]

    alpha_face_test = sp.symbols("alpha_face", real=True)
    cf_face = sp.symbols("cf_face", real=True)
    h_faceDD = ixp.declarerank2("h_faceDD", symmetry="sym01", dimension=3)
    vet_faceU = ixp.declarerank1("vet_faceU", dimension=3)

    # We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
    # on the right and left faces
    rescaledvrU = ixp.declarerank1("rescaledvrU", dimension=3)
    rescaledvlU = ixp.declarerank1("rescaledvlU", dimension=3)
    BmagrU = ixp.declarerank1("BmagrU", dimension=3)
    BmaglU = ixp.declarerank1("BmaglU", dimension=3)

    VrU = ixp.zerorank1()
    VlU = ixp.zerorank1()

    u4rUt = sp.symbols("u4rUt", real=True)
    u4lUt = sp.symbols("u4lUt", real=True)
    u4rU = ixp.zerorank1(dimension=4)
    u4rU[0] = u4rUt
    u4lU = ixp.zerorank1(dimension=4)
    u4lU[0] = u4lUt

    beta_faceU_test = ixp.zerorank1()
    for i in range(3):
        VrU[i] = rescaledvrU[i] * rfm.ReU[i]
        VlU[i] = rescaledvlU[i] * rfm.ReU[i]
        u4rU[i + 1] = VrU[i] * u4rU[0]
        u4lU[i + 1] = VlU[i] * u4lU[0]
        beta_faceU_test[i] = vet_faceU[i] * rfm.ReU[i]

    rho_b_r_test = sp.symbols("rhob_r", real=True)
    rho_b_l_test = sp.symbols("rhob_l", real=True)

    h_r_test = sp.symbols("h_r", real=True)
    h_l_test = sp.symbols("h_l", real=True)

    cs2_r_test = sp.symbols("cs2_r", real=True)
    cs2_l_test = sp.symbols("cs2_l", real=True)

    # ADM in terms of BSSN
    AitoB = BSSN_to_ADM(CoordSystem="Cartesian")

    Bq = BSSN_quantities["Cartesian"]

    gamma_faceDD_test = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            gamma_faceDD_test[i][j] = (
                AitoB.gammaDD[i][j]
                .subs(Bq.hDD[i][j], h_faceDD[i][j])
                .subs(Bq.cf, cf_face)
            )

    # Step 1: With B^i = 0, the bounds must equal the GRHD bounds in every direction.
    zero_BmagU = ixp.zerorank1(dimension=3)
    for flux_dirn_test in range(3):
        if find_cmax_cmin(
            flux_dirn_test,
            gamma_faceDD_test,
            beta_faceU_test,
            alpha_face_test,
            u4rU,
            u4lU,
            zero_BmagU,
            zero_BmagU,
            rho_b_r_test,
            rho_b_l_test,
            h_r_test,
            h_l_test,
            cs2_r_test,
            cs2_l_test,
        ) != find_grhd_cmax_cmin(
            flux_dirn_test,
            gamma_faceDD_test,
            beta_faceU_test,
            alpha_face_test,
            u4rU,
            u4lU,
            cs2_r_test,
            cs2_l_test,
        ):
            raise AssertionError(
                f"B=0 speed bounds differ from GRHD in flux direction {flux_dirn_test}"
            )

    # Step 2: Each face state's v_0^2 must reach the GRHD bounds in every direction.
    #         A flat face metric keeps these exact comparisons small.
    flat_gammaDD = ixp.zerorank2(dimension=3)
    for i in range(3):
        flat_gammaDD[i][i] = sp.sympify(1)
    zero_betaU = ixp.zerorank1(dimension=3)
    unit_alpha = sp.sympify(1)
    flat_v02_r = compute_v02(
        flat_gammaDD,
        zero_betaU,
        unit_alpha,
        u4rU,
        BmagrU,
        rho_b_r_test,
        h_r_test,
        cs2_r_test,
    )
    flat_v02_l = compute_v02(
        flat_gammaDD,
        zero_betaU,
        unit_alpha,
        u4lU,
        BmaglU,
        rho_b_l_test,
        h_l_test,
        cs2_l_test,
    )
    for flux_dirn_test in range(3):
        if find_cmax_cmin(
            flux_dirn_test,
            flat_gammaDD,
            zero_betaU,
            unit_alpha,
            u4rU,
            u4lU,
            BmagrU,
            BmaglU,
            rho_b_r_test,
            rho_b_l_test,
            h_r_test,
            h_l_test,
            cs2_r_test,
            cs2_l_test,
        ) != find_grhd_cmax_cmin(
            flux_dirn_test,
            flat_gammaDD,
            zero_betaU,
            unit_alpha,
            u4rU,
            u4lU,
            flat_v02_r,
            flat_v02_l,
        ):
            raise AssertionError(
                f"Face-state v_0^2 values are misassigned in flux direction {flux_dirn_test}"
            )

    # Step 3: Compare magnetic signal speeds with trusted values.
    v02_l_test = compute_v02(
        gamma_faceDD_test,
        beta_faceU_test,
        alpha_face_test,
        u4lU,
        BmaglU,
        rho_b_l_test,
        h_l_test,
        cs2_l_test,
    )
    g4_faceUU = ADM_to_g4UU(gamma_faceDD_test, beta_faceU_test, alpha_face_test)

    (
        cminus_tmp_test,
        cplus_tmp_test,
    ) = find_cp_cm(1, g4_faceUU, u4lU, v02_l_test)

    exprs_dict["cminus"] = _nrpyAbs_to_Abs(cminus_tmp_test)
    exprs_dict["cplus"] = _nrpyAbs_to_Abs(cplus_tmp_test)

    (
        cmin_tmp_test,
        cmax_tmp_test,
    ) = find_cmax_cmin(
        2,
        gamma_faceDD_test,
        beta_faceU_test,
        alpha_face_test,
        u4rU,
        u4lU,
        BmagrU,
        BmaglU,
        rho_b_r_test,
        rho_b_l_test,
        h_r_test,
        h_l_test,
        cs2_r_test,
        cs2_l_test,
    )

    exprs_dict["cmin"] = _nrpyAbs_to_Abs(cmin_tmp_test)
    exprs_dict["cmax"] = _nrpyAbs_to_Abs(cmax_tmp_test)

    results_dict = ve.process_dictionary_of_expressions(
        exprs_dict, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )
