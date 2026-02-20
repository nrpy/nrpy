"""
Construct characteristic speeds at left and right interfaces of grid cells.
These are to be fed into an approximate Riemann solver.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

# Step Import needed modules:
from typing import Any, Dict, List, Tuple

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends

import nrpy.equations.grhd.Min_Max_and_Piecewise_Expressions as noif
from nrpy.equations.general_relativity.g4munu_conversions import ADM_to_g4UU


def find_cp_cm(
    flux_dirn: int,
    g4UU: List[List[sp.Expr]],
    u4U: List[sp.Expr],
    cs2: sp.Expr,
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Compute the in-going and out-going characteristic speeds c_+ and c_-.

    :param flux_dirn: direction for flux calculation (int)
    :param g4UU: contravariant form of the four-dimensional spacetime tensor g^{mu nu}
    :param u4U: four-velocity u^mu
    :param cs2: sound speed squared

    :return: symbolic expressions of characteristic speeds c_+ and c_-, cminus and cplus

    ::note: Here we actually need to ensure that we're using rescaled quantities
    for the speeds and spatial metric. However, as discussed in the fluxes module,
    a successful rescaling amounts to the same symbolic expressions as just using
    Cartesian coordinates.

    """
    v02 = cs2
    a = (1 - v02) * (u4U[0] ** 2) - v02 * g4UU[0][0]
    b = 2 * v02 * g4UU[flux_dirn + 1][0] - 2 * u4U[flux_dirn + 1] * u4U[0] * (1 - v02)
    c = (1 - v02) * (u4U[flux_dirn + 1] ** 2) - v02 * g4UU[flux_dirn + 1][flux_dirn + 1]

    # Now, we are free to solve the quadratic equation as usual. We take care to avoid passing a
    # negative value to the sqrt function.
    detm = b * b - sp.sympify(4) * a * c

    detm = sp.sqrt(noif.max_noif(sp.sympify(0), detm))

    # note that these correspond to a single interface, left or right
    cplus_tmp = sp.Rational(1, 2) * (-b / a + detm / a)
    cminus_tmp = sp.Rational(1, 2) * -(b / a + detm / a)

    cminus = noif.min_noif(cplus_tmp, cminus_tmp)
    cplus = noif.max_noif(cplus_tmp, cminus_tmp)

    # the above in C code
    # if (cplus < cminus) {
    # CCTK_REAL cp = cminus;
    # cminus = cplus;
    # cplus = cp;

    return cminus, cplus


def find_cmax_cmin(
    flux_dirn: int,
    gamma_faceDD: List[List[sp.Expr]],
    beta_faceU: List[sp.Expr],
    alpha_face: sp.Expr,
    u4U_r: List[sp.Expr],
    u4U_l: List[sp.Expr],
    cs2_r: sp.Expr,
    cs2_l: sp.Expr,
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Compute the maximum and minimum characteristic speeds c_max and c_min.

    :param flux_dirn: direction for flux calculation (int)
    :param gamma_faceDD: spatial metric at the cell interface.
    :param beta_faceU: shift vector at cell interface
    :param alpha_face: lapse function at cell interface
    :param u4U_r: four-velocity u^mu at the right interface of a grid cell
    :param u4U_l: four-velocity u^mu at the left interface of a grid cell
    :param cs2_r: sound speed squared at the right interface of a grid cell
    :param cs2_l: sound speed squared at the left interface of a grid cell

    :return: symbolic expressions of the maximum and minimum
             characteristic speeds c_max and c_min

    ::note: Note that we do not distinguish metric quantites between left and right
            interfaces. This is because these quantities are usually smooth.

    """
    # First, we need to find the characteristic speeds on each face
    g4UU = ADM_to_g4UU(gamma_faceDD, beta_faceU, alpha_face)

    # Original needed for GRMHD
    cmr, cpr = find_cp_cm(flux_dirn, g4UU, u4U_r, cs2_r)

    cml, cpl = find_cp_cm(flux_dirn, g4UU, u4U_l, cs2_l)

    # The following algorithms have been verified with random floats:

    #   // Then compute cmax, cmin. This is required for the HLL flux.
    #   original C code
    #   CCTK_REAL cmaxL =  MAX(0.0,MAX(cplusl,cplusr));
    #   CCTK_REAL cminL = -MIN(0.0,MIN(cminusl,cminusr));

    # Now, we need to set cmax to the larger of cpr,cpl, and 0
    cmax = noif.max_noif(sp.sympify(0.0), noif.max_noif(cpl, cpr))
    # And then, set cmin to the smaller of cmr,cml, and 0
    cmin = -noif.min_noif(sp.sympify(0.0), noif.min_noif(cml, cmr))

    # save the rescaled char. speeds
    return cmin, cmax


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

    exprs_dict: Dict[str, Any] = {}

    import nrpy.indexedexp as ixp  # NRPy: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
    import nrpy.reference_metric as refmetric  # NRPy: Reference metric support
    from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
    from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM

    rfm = refmetric.reference_metric["Cartesian"]

    alpha_face_test = sp.symbols("alpha_face", real=True)
    cf_face = sp.symbols("cf_face", real=True)
    h_faceDD = ixp.declarerank2("h_faceDD", symmetry="sym01", dimension=3)
    vet_faceU = ixp.declarerank1("vet_faceU", dimension=3)

    # We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
    # on the right and left faces
    rescaledvrU = ixp.declarerank1("rescaledvrU", dimension=3)
    rescaledvlU = ixp.declarerank1("rescaledvlU", dimension=3)

    VU = ixp.zerorank1()
    VrU = ixp.zerorank1()
    VlU = ixp.zerorank1()

    u4rUt = sp.symbols("u4rUt", real=True)
    u4lUt = sp.symbols("u4lUt", real=True)
    u4rU = ixp.zerorank1(dimension=4)
    u4rU[0] = u4rUt
    u4lU = ixp.zerorank1(dimension=4)
    u4lU[0] = u4lUt

    betaU = ixp.zerorank1()
    beta_faceU_test = ixp.zerorank1()
    for i in range(3):
        VrU[i] = rescaledvrU[i] * rfm.ReU[i]
        VlU[i] = rescaledvlU[i] * rfm.ReU[i]
        u4rU[i + 1] = VrU[i] * u4rU[0]
        u4lU[i + 1] = VlU[i] * u4lU[0]
        beta_faceU_test[i] = vet_faceU[i] * rfm.ReU[i]

    rho_b_r = sp.symbols("rhob_r", real=True)
    rho_b_l = sp.symbols("rhob_l", real=True)

    h_r = sp.symbols("h_r", real=True)
    h_l = sp.symbols("h_l", real=True)

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

    g4_faceUU = ADM_to_g4UU(gamma_faceDD_test, beta_faceU_test, alpha_face_test)

    (
        cminus_tmp_test,
        cplus_tmp_test,
    ) = find_cp_cm(1, g4_faceUU, u4lU, cs2_r_test)

    exprs_dict["cminus"] = cminus_tmp_test.subs(sp.Function("nrpyAbs"), sp.Abs)
    exprs_dict["cplus"] = cplus_tmp_test.subs(sp.Function("nrpyAbs"), sp.Abs)

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
        cs2_r_test,
        cs2_l_test,
    )

    exprs_dict["cmin"] = cmin_tmp_test.subs(sp.Function("nrpyAbs"), sp.Abs)
    exprs_dict["cmax"] = cmax_tmp_test.subs(sp.Function("nrpyAbs"), sp.Abs)

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
