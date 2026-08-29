"""
Construct lapse and shift right-hand sides for fCCZ4.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple

import sympy as sp

import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.BSSN_gauge_RHSs import BSSN_gauge_RHSs
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.fCCZ4_RHSs import fCCZ4_RHSs


def fCCZ4_gauge_RHSs(
    CoordSystem: str = "Cartesian",
    enable_rfm_precompute: bool = False,
    enable_T4munu: bool = False,
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    enable_YBS_Gamma_constraint_adjustment: bool = False,
) -> Tuple[sp.Expr, List[sp.Expr], List[sp.Expr]]:
    """
    Build supported fCCZ4 moving-puncture gauge equations.

    :param CoordSystem: Reference-metric coordinate system.
    :param enable_rfm_precompute: Use precomputed reference-metric symbols.
    :param enable_T4munu: Include matter in connection evolution.
    :param LapseEvolutionOption: Any lapse condition supported by
        ``BSSN_gauge_RHSs``.
    :param ShiftEvolutionOption: Any shift condition supported by
        ``BSSN_gauge_RHSs``.
    :param enable_YBS_Gamma_constraint_adjustment: Enable the YBS
        connection-constraint adjustment.
    :return: Lapse RHS and rescaled shift and driver RHS vectors.

    """
    alpha_rhs, vet_rhsU, bssn_bet_rhsU = BSSN_gauge_RHSs(
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_T4munu=enable_T4munu,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
        enable_YBS_Gamma_constraint_adjustment=(
            enable_YBS_Gamma_constraint_adjustment
        ),
    )
    suffix = CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    Bq = BSSN_quantities[suffix]
    rhs = fCCZ4_RHSs.get_rhs(
        suffix + ("_T4munu" if enable_T4munu else ""),
        enable_YBS_Gamma_constraint_adjustment=(
            enable_YBS_Gamma_constraint_adjustment
        ),
    )
    if LapseEvolutionOption == "OnePlusLog":
        alpha_rhs += 4 * Bq.alpha * rhs.Theta

    bet_rhsU = [bssn_bet_rhsU[i] for i in range(3)]
    if ShiftEvolutionOption in (
        "GammaDriving2ndOrder_NoCovariant",
        "GammaDriving2ndOrder_Covariant",
        "GammaDriving2ndOrder_Covariant__Hatted",
        "NonAdvectingGammaDriving",
    ):
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        for i in range(3):
            bet_rhsU[i] += (
                sp.Rational(3, 4) * rhs.Lambdatilde_rhsU_delta[i] / rfm.ReU[i]
            )
    return alpha_rhs, vet_rhsU, bet_rhsU


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

    lapse_options = (
        "OnePlusLog",
        "BHSHarmonicSlicing",
        "Frozen",
        "OnePlusLogAlt",
    )
    shift_options = (
        "Frozen",
        "GammaDriving2ndOrder_NoCovariant",
        "GammaDriving2ndOrder_Covariant",
        "GammaDriving2ndOrder_Covariant__Hatted",
        "GammaDriving1stOrder_Covariant",
        "GammaDriving1stOrder_Covariant__Hatted",
        "NonAdvectingGammaDriving",
    )
    cases = [
        (coord, False, False, lapse, shift, coord)
        for coord in ("SinhCartesian", "SinhSpherical")
        for lapse in lapse_options
        for shift in shift_options
    ]
    cases.extend(
        [
            (
                "Cartesian",
                False,
                False,
                "OnePlusLog",
                "GammaDriving2ndOrder_Covariant__Hatted",
                "Cartesian",
            ),
            (
                "SinhSpherical",
                True,
                True,
                "OnePlusLog",
                "GammaDriving2ndOrder_Covariant__Hatted",
                "SinhSpherical_rfm_precompute_T4munu",
            ),
        ]
    )
    for (
        case_coord,
        case_enable_rfm_precompute,
        case_enable_T4munu,
        case_lapse,
        case_shift,
        trusted_coord_suffix,
    ) in cases:
        case_alpha_rhs, case_vet_rhsU, case_bet_rhsU = fCCZ4_gauge_RHSs(
            CoordSystem=case_coord,
            enable_rfm_precompute=case_enable_rfm_precompute,
            enable_T4munu=case_enable_T4munu,
            LapseEvolutionOption=case_lapse,
            ShiftEvolutionOption=case_shift,
        )
        expression_dict = {
            "alpha_rhs": case_alpha_rhs,
            "vet_rhsU0": case_vet_rhsU[0],
            "vet_rhsU1": case_vet_rhsU[1],
            "vet_rhsU2": case_vet_rhsU[2],
            "bet_rhsU0": case_bet_rhsU[0],
            "bet_rhsU1": case_bet_rhsU[1],
            "bet_rhsU2": case_bet_rhsU[2],
        }
        processed = ve.process_dictionary_of_expressions(
            expression_dict, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{case_lapse}_{case_shift}_{trusted_coord_suffix}",
            processed,
        )
