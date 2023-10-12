"""
Construct expressions for BSSN stress-energy source terms, in terms of elements of T^{mu nu}.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Sequence, Optional, Tuple, Dict, Any
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
import nrpy.equations.general_relativity.g4munu_conversions as g4conv  # NRPy+: ADM/BSSN <-> 4-metric conversions


def T4UU_and_ADM_to_SDD_SD_S_rho(
    gammaDD: Optional[Sequence[Sequence[sp.Expr]]] = None,
    betaU: Optional[Sequence[sp.Expr]] = None,
) -> Tuple[Sequence[Sequence[sp.Expr]], Sequence[sp.Expr], sp.Expr, sp.Expr]:
    """
    Define BSSN source terms in terms of T^{mu nu} and ADM variables.

    :param gammaDD: 3x3 spatial part of the ADM metric tensor, gamma_{ij}
    :return: Returns a tuple containing the BSSN source terms SDD, SD, S, and rho.
    """
    # Step 1: Define gamma4DD[mu][nu] = g_{mu nu} + n_{mu} n_{nu}
    alpha = sp.symbols("alpha", real=True)
    if betaU is None:
        betaU = ixp.declarerank1("betaU")
    if gammaDD is None:
        gammaDD = ixp.declarerank2("gammaDD", symmetry="sym01")

    zero = sp.sympify(0)
    negone = sp.sympify(-1)
    n4D = [negone * alpha, zero, zero, zero]
    g4DD = g4conv.ADM_to_g4DD(gammaDD, betaU, alpha)

    gamma4DD = ixp.zerorank2(dimension=4)
    for mu in range(4):
        for nu in range(4):
            gamma4DD[mu][nu] = g4DD[mu][nu] + n4D[mu] * n4D[nu]

    # Step 2: Declare T4UU here
    T4UU = ixp.declarerank2("T4UU", symmetry="sym01", dimension=4)

    # Step 3: Define BSSN source terms SDD, SD, S, and rho
    # Step 3.a: S_{ij} = gamma_{i mu} gamma_{j nu} T^{mu nu}
    SDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for mu in range(4):
                for nu in range(4):
                    SDD[i][j] += (
                        gamma4DD[i + 1][mu] * gamma4DD[j + 1][nu] * T4UU[mu][nu]
                    )

    # Step 3.b: S_{i} = -gamma_{i mu} n_{nu} T^{mu nu}
    SD = ixp.zerorank1()
    for i in range(3):
        for mu in range(4):
            for nu in range(4):
                SD[i] += -gamma4DD[i + 1][mu] * n4D[nu] * T4UU[mu][nu]

    # Step 3.c: S = gamma^{ij} S_{ij}
    gammaUU, _ = ixp.symm_matrix_inverter3x3(gammaDD)
    S = zero
    for i in range(3):
        for j in range(3):
            S += gammaUU[i][j] * SDD[i][j]

    # Step 3.d: rho = n_{mu} n_{nu} T^{mu nu}
    rho = zero
    for mu in range(4):
        for nu in range(4):
            rho += n4D[mu] * n4D[nu] * T4UU[mu][nu]

    return SDD, SD, S, rho


def T4UU_and_BSSN_to_SDD_SD_S_rho(
    CoordSystem: str = "Cartesian", enable_rfm_precompute: bool = False
) -> Tuple[Sequence[Sequence[sp.Expr]], Sequence[sp.Expr], sp.Expr, sp.Expr]:
    """
    Define BSSN source terms in terms of T^{mu nu} and BSSN variables.

    :param CoordSystem: Coordinate system to use, defaults to "Cartesian"
    :param enable_rfm_precompute: Flag for enabling reference metric precomputation, defaults to False
    :return: Returns a tuple containing the BSSN source terms SDD, SD, S, and rho.
    """
    BtoA = BSSN_to_ADM(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )
    return T4UU_and_ADM_to_SDD_SD_S_rho(gammaDD=BtoA.gammaDD, betaU=BtoA.betaU)


# Step 3: Add BSSN stress-energy source terms to BSSN RHSs
def BSSN_RHSs_T4UU_source_terms(
    CoordSystem: str = "Cartesian",
    enable_rfm_precompute: bool = False,
) -> Tuple[sp.Expr, Sequence[Sequence[sp.Expr]], Sequence[sp.Expr]]:
    """
    Compute BSSN source terms from T^{mu nu} and contribute them to BSSN RHSs.

    :param CoordSystem: Coordinate system to use, defaults to "Cartesian"
    :param enable_rfm_precompute: Flag for enabling reference metric precomputation, defaults to False
    :return: Returns a tuple containing the BSSN source terms for trK_rhs, Abar_rhsDD, and Lambdabar_rhsU.
    """
    SDD, SD, S, rho = T4UU_and_BSSN_to_SDD_SD_S_rho(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )
    PI = par.register_CodeParameter(
        "REAL",
        __name__,
        "PI",
        "3.14159265358979323846264338327950288",
        add_to_glb_code_params_dict=True,
        add_to_set_CodeParameters_h=True,
        add_to_parfile=False,
    )
    alpha = sp.symbols("alpha", real=True)
    Bq = BSSN_quantities[
        f"{CoordSystem}_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    sourceterm_trK_rhs = 4 * PI * alpha * (rho + S)

    gammabarUU, _ = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)
    tracefree_SDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            tracefree_SDD[i][j] = SDD[i][j]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for m in range(3):
                    tracefree_SDD[i][j] += (
                        -sp.Rational(1, 3)
                        * Bq.gammabarDD[i][j]
                        * gammabarUU[k][m]
                        * SDD[k][m]
                    )

    sourceterm_Abar_rhsDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            sourceterm_Abar_rhsDD[i][j] = (
                -8 * PI * alpha * Bq.exp_m4phi * tracefree_SDD[i][j]
            )

    sourceterm_Lambdabar_rhsU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            sourceterm_Lambdabar_rhsU[i] += -16 * PI * alpha * gammabarUU[i][j] * SD[j]

    return (
        sourceterm_trK_rhs,
        sourceterm_Abar_rhsDD,
        sourceterm_Lambdabar_rhsU,
    )


# Add BSSN stress-energy source terms to BSSN constraints
def BSSN_constraints_T4UU_source_terms(
    CoordSystem: str = "Cartesian", enable_rfm_precompute: bool = False
) -> Tuple[sp.Expr, Sequence[sp.Expr]]:
    """
    Add BSSN stress-energy source terms to BSSN constraints.

    :param CoordSystem: Coordinate system, default is "Cartesian"
    :param enable_rfm_precompute: Enable reference metric precomputation, default is False

    :return: Tuple containing sourceterm_H and sourceterm_MU
    """
    # Step 1: Call BSSN_source_terms_ito_T4UU to get SDD, SD, S, & rho
    _SDD, SD, _S, rho = T4UU_and_BSSN_to_SDD_SD_S_rho(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )  # _SDD, _S unused

    PI = par.register_CodeParameter(
        "REAL",
        __name__,
        "PI",
        "3.14159265358979323846264338327950288",
        add_to_glb_code_params_dict=True,
        add_to_set_CodeParameters_h=True,
        add_to_parfile=False,
    )

    # Step 2: Add source term to the Hamiltonian constraint H
    sourceterm_H = -16 * PI * rho

    # Step 3: Add source term to the momentum constraint M^i

    # Step 3.a: Compute gammaUU in terms of BSSN quantities
    BtoA = BSSN_to_ADM(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )

    # Step 3.b: Raise S_i
    SU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            SU[i] += BtoA.gammaUU[i][j] * SD[j]

    # Step 3.c: Add source term to momentum constraint
    sourceterm_MU = ixp.zerorank1()
    for i in range(3):
        sourceterm_MU[i] = -8 * PI * SU[i]

    return sourceterm_H, sourceterm_MU


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

    # Notify BSSN_quantities that T4munu is enabled.
    par.set_parval_from_str("enable_T4munu", True)
    exprs_dict: Dict[str, Any] = {}
    (
        exprs_dict["ADMSDD"],
        exprs_dict["ADMSD"],
        exprs_dict["ADMS"],
        exprs_dict["ADMrho"],
    ) = T4UU_and_ADM_to_SDD_SD_S_rho()

    for Coord in [
        "Spherical",
        "SinhSpherical",
        "SinhSpherical_rfm_precompute",
        "Cartesian",
        "SinhCartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        enable_rfm_pre = "rfm_precompute" in Coord
        (
            exprs_dict[f"BSSNSDD_{Coord}"],
            exprs_dict[f"BSSNSD_{Coord}"],
            exprs_dict[f"BSSNS_{Coord}"],
            exprs_dict[f"BSSNrho_{Coord}"],
        ) = T4UU_and_BSSN_to_SDD_SD_S_rho(
            CoordSystem=Coord.replace("_rfm_precompute", ""),
            enable_rfm_precompute=enable_rfm_pre,
        )
        (
            exprs_dict[f"BSSNsourceterm_trK_rhs_{Coord}"],
            exprs_dict[f"BSSNsourceterm_a_rhs_{Coord}"],
            exprs_dict[f"BSSNsourceterm_lambda_rhs_{Coord}"],
        ) = BSSN_RHSs_T4UU_source_terms(
            CoordSystem=Coord.replace("_rfm_precompute", ""),
            enable_rfm_precompute=enable_rfm_pre,
        )
        (
            exprs_dict[f"BSSNsourceterm_H_{Coord}"],
            exprs_dict[f"BSSNsourceterm_MU_{Coord}"],
        ) = BSSN_constraints_T4UU_source_terms(
            CoordSystem=Coord.replace("_rfm_precompute", ""),
            enable_rfm_precompute=enable_rfm_pre,
        )
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
