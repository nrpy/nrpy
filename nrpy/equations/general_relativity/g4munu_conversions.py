"""
Construct expressions for ADM or BSSN quantities in terms of the 4-metric g4DD, and g4DD/g4UU in terms of ADM/BSSN quantities.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""


from typing import Optional, Sequence, Tuple, Dict, Any
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.helpers.cached_functions import cached_simplify


# g_{mu nu} in terms of BSSN (if inputvars=="BSSN") or ADM (if inputvars=="ADM") variables.
def ADM_to_g4DD(
    gammaDD: Sequence[Sequence[sp.Expr]],
    betaU: Sequence[sp.Expr],
    alpha: sp.Expr,
) -> Sequence[Sequence[sp.Expr]]:
    """
    Convert the 3+1 decomposed ADM variables to the 4-metric g_{mu nu}.

    :param gammaDD: 3x3 spatial metric tensor, gamma_{ij}.
    :param betaU: Shift vector, beta^i.
    :param alpha: Lapse function, alpha.
    :return: 4x4 metric tensor, g_{mu nu}.
    """
    # Step 0: Initialize 4x4 tensor to store g_{mu nu}
    g4DD = ixp.zerorank2(dimension=4)

    # Step 1: Compute beta_i via Eq. 2.121 in B&S
    betaD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            betaD[i] += gammaDD[i][j] * betaU[j]

    # Step 2: Compute beta_i beta^i, the beta contraction
    beta2 = sp.sympify(0)
    for i in range(3):
        beta2 += betaU[i] * betaD[i]

    # Step 3: Construct g4DD via Eq. 2.122 in B&S
    g4DD[0][0] = -(alpha**2) + beta2
    for mu in range(1, 4):
        g4DD[mu][0] = g4DD[0][mu] = betaD[mu - 1]
    for mu in range(1, 4):
        for nu in range(1, 4):
            g4DD[mu][nu] = gammaDD[mu - 1][nu - 1]

    return g4DD


def BSSN_to_g4DD(
    CoordSystem: str = "Cartesian", enable_rfm_precompute: bool = False
) -> Sequence[Sequence[sp.Expr]]:
    """
    Convert BSSN variables to the 4-metric g_{mu nu}.

    :param CoordSystem: Coordinate system, default is "Cartesian".
    :param enable_rfm_precompute: Enable reference metric pre-computation, default is False.
    :return: 4x4 metric tensor, g_{mu nu}.
    """
    # Set gammaDD, betaU, and alpha in terms of generic BSSN variables
    BtoA = BSSN_to_ADM(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )
    gammaDD = BtoA.gammaDD
    Bq = BSSN_quantities[
        f"{CoordSystem}_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    betaU = Bq.betaU
    alpha = sp.symbols("alpha", real=True)

    return ADM_to_g4DD(gammaDD, betaU, alpha)


def ADM_to_g4UU(
    gammaDD: Sequence[Sequence[sp.Expr]],
    betaU: Sequence[sp.Expr],
    alpha: sp.Expr,
    gammaUU: Optional[Sequence[Sequence[sp.Expr]]] = None,
) -> Sequence[Sequence[sp.Expr]]:
    """
    Construct the contravariant 4-metric tensor, g^{mu nu}, using ADM variables.

    :param gammaDD: 3x3 spatial metric tensor, gamma_{ij}.
    :param betaU: Shift vector, beta^i.
    :param alpha: Lapse function, alpha.
    :param gammaUU: Optional 3x3 spatial metric tensor, gamma^{ij}. Default is None.
    :return: 4x4 contravariant metric tensor, g^{mu nu}.
    """
    # Compute gammaUU based on provided gammaDD:
    if gammaUU is None:
        gammaUU, _ = ixp.symm_matrix_inverter3x3(gammaDD)

    # Then evaluate g4UU:
    g4UU = ixp.zerorank2(dimension=4)
    g4UU[0][0] = -1 / alpha**2
    for mu in range(1, 4):
        g4UU[0][mu] = g4UU[mu][0] = betaU[mu - 1] / alpha**2
    for mu in range(1, 4):
        for nu in range(1, 4):
            g4UU[mu][nu] = (
                gammaUU[mu - 1][nu - 1] - betaU[mu - 1] * betaU[nu - 1] / alpha**2
            )

    return g4UU


def BSSN_to_g4UU(
    CoordSystem: str = "Cartesian", enable_rfm_precompute: bool = False
) -> Sequence[Sequence[sp.Expr]]:
    """
    Convert BSSN variables to the contravariant 4-metric g^{mu nu}.

    :param CoordSystem: Coordinate system, default is "Cartesian".
    :param enable_rfm_precompute: Enable reference metric pre-computation, default is False.
    :return: 4x4 contravariant metric tensor, g^{mu nu}.
    """
    BtoA = BSSN_to_ADM(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )
    gammaDD = BtoA.gammaDD
    Bq = BSSN_quantities[
        f"{CoordSystem}_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    betaU = Bq.betaU
    alpha = sp.symbols("alpha", real=True)

    return ADM_to_g4UU(gammaDD, betaU, alpha)


# g_{mu nu} -> ADM variables
def g4DD_to_ADM(
    g4DD: Optional[Sequence[Sequence[sp.Expr]]] = None,
) -> Tuple[Sequence[Sequence[sp.Expr]], sp.Expr, Sequence[sp.Expr]]:
    """
    Convert a 4-metric tensor (g_{mu nu}) to ADM variables (gamma_{ij}), (alpha), (beta^i).

    :param g4DD: The input 4-metric tensor (g_{mu nu}) represented as a 4x4 matrix.
    :return: Tuple containing (gamma_{ij}), (alpha), and (beta^i) in the same order.
    """
    # Step 1: Declare g4DD as a 4-metric tensor:
    g4DD_is_input_into_this_function = True
    if g4DD is None:
        g4DD_is_input_into_this_function = False
        g4DD = ixp.declarerank2("g4DD", symmetry="sym01", dimension=4)

    # Step 2: Compute gammaDD & betaD from the 4-metric tensor
    betaD = ixp.zerorank1()
    gammaDD = ixp.zerorank2()
    for i in range(3):
        betaD[i] = g4DD[0][i]
        for j in range(3):
            gammaDD[i][j] = g4DD[i + 1][j + 1]

    # Step 3: Compute betaU based on the 4-metric tensor
    # Step 3.a: Compute gammaUU based on provided gammaDD
    gammaUU, _gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)  # _gammaDET is unused.

    # Step 3.b: Use gammaUU to raise betaU
    betaU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            betaU[i] += gammaUU[i][j] * betaD[j]

    # Step 4: Compute alpha based on the 4-metric tensor
    # Step 4.a: Compute beta^2 = beta^k beta_k:
    beta_squared = sp.sympify(0)
    for k in range(3):
        beta_squared += betaU[k] * betaD[k]

    # Step 4.b: alpha = sqrt(beta^2 - g_{00}):
    alpha = (
        sp.sqrt(cached_simplify(beta_squared) - g4DD[0][0])
        if not g4DD_is_input_into_this_function
        else sp.sqrt(beta_squared - g4DD[0][0])
    )

    return gammaDD, alpha, betaU


def g4DD_to_BSSN(
    g4DD: Optional[Sequence[Sequence[sp.Expr]]] = None,
    CoordSystem: str = "Cartesian",
    enable_rfm_precompute: bool = False,
) -> Tuple[Sequence[Sequence[sp.Expr]], sp.Expr, Sequence[sp.Expr], sp.Expr]:
    """
    Convert a 4-metric tensor (g_{mu nu}) to BSSN variables.

    :param g4DD: The input 4-metric tensor (g_{mu nu}).
    :param CoordSystem: The coordinate system to be used. Default is "Cartesian".
    :param enable_rfm_precompute: Whether to enable precomputation of reference metric quantities.
    :return: Tuple containing BSSN variables (hDD), (cf), (vetU), and (alpha).
    """
    gammaDD, alpha, betaU = g4DD_to_ADM(g4DD)

    # Convert ADM to BSSN based on the 4-metric tensor
    dummyBU = ixp.zerorank1()
    dummyKDD = ixp.zerorank2()
    AtoB = ADM_to_BSSN(
        gammaDD,
        dummyKDD,
        betaU,
        dummyBU,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
    )
    return AtoB.hDD, AtoB.cf, AtoB.vetU, alpha


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
    ingammaDD = ixp.declarerank2("gammaDD", symmetry="sym01")
    inbetaU = ixp.declarerank1("betaU")
    inalpha = sp.Symbol("alpha", real=True)
    exprs_dict["ADMg4DD"] = ADM_to_g4DD(ingammaDD, inbetaU, inalpha)
    exprs_dict["ADMg4UU"] = ADM_to_g4UU(ingammaDD, inbetaU, inalpha)
    (
        exprs_dict["ADMgammaDD_ito_g4DD"],
        exprs_dict["ADMalpha_ito_g4DD"],
        exprs_dict["ADMbetaU_ito_g4DD"],
    ) = g4DD_to_ADM()
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
        exprs_dict[f"BSSNg4DD_{Coord}"] = BSSN_to_g4DD(
            CoordSystem=Coord.replace("_rfm_precompute", ""),
            enable_rfm_precompute=enable_rfm_pre,
        )
        exprs_dict[f"BSSNg4UU_{Coord}"] = BSSN_to_g4UU(
            CoordSystem=Coord.replace("_rfm_precompute", ""),
            enable_rfm_precompute=enable_rfm_pre,
        )
        (
            exprs_dict[f"BSSNhDD_ito_g4DD_{Coord}"],
            exprs_dict[f"BSSNcf_ito_g4DD_{Coord}"],
            exprs_dict[f"BSSNvetU_ito_g4DD_{Coord}"],
            exprs_dict[f"BSSNalpha_ito_g4DD_{Coord}"],
        ) = g4DD_to_BSSN(
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
