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
) -> Tuple[sp.Expr, List[sp.Expr], List[sp.Expr]]:
    """
    Build supported fCCZ4 moving-puncture gauge equations.

    :param CoordSystem: Reference-metric coordinate system.
    :param enable_rfm_precompute: Use precomputed reference-metric symbols.
    :param enable_T4munu: Include matter in connection evolution.
    :param LapseEvolutionOption: ``OnePlusLog`` or ``Frozen``.
    :param ShiftEvolutionOption: Hatted covariant second-order Gamma driver,
        nonadvecting Gamma driver, or ``Frozen``.
    :return: Lapse RHS and rescaled shift and driver RHS vectors.
    :raises ValueError: If a lapse or shift option is unsupported.

    """
    if LapseEvolutionOption not in ("OnePlusLog", "Frozen"):
        raise ValueError(
            f"Error: LapseEvolutionOption == {LapseEvolutionOption} not supported!"
        )
    supported_shifts = (
        "GammaDriving2ndOrder_Covariant__Hatted",
        "NonAdvectingGammaDriving",
        "Frozen",
    )
    if ShiftEvolutionOption not in supported_shifts:
        raise ValueError(
            f"Error: ShiftEvolutionOption == {ShiftEvolutionOption} not supported!"
        )

    alpha_rhs, vet_rhsU, bssn_bet_rhsU = BSSN_gauge_RHSs(
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_T4munu=enable_T4munu,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    suffix = CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    Bq = BSSN_quantities[suffix]
    rhs = fCCZ4_RHSs[suffix + ("_T4munu" if enable_T4munu else "")]
    if LapseEvolutionOption == "OnePlusLog":
        alpha_rhs += 4 * Bq.alpha * rhs.Theta

    bet_rhsU = [bssn_bet_rhsU[i] for i in range(3)]
    if ShiftEvolutionOption in (
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

    unsupported_lapse_residual = sp.Integer(1)
    try:
        fCCZ4_gauge_RHSs(LapseEvolutionOption="BHSHarmonicSlicing")
    except ValueError as error:
        unsupported_lapse_residual = sp.Integer(
            0
            if str(error)
            == "Error: LapseEvolutionOption == BHSHarmonicSlicing not supported!"
            else 1
        )

    unsupported_shift_residual = sp.Integer(1)
    try:
        fCCZ4_gauge_RHSs(ShiftEvolutionOption="GammaDriving2ndOrder_Covariant")
    except ValueError as error:
        unsupported_shift_residual = sp.Integer(
            0
            if str(error)
            == "Error: ShiftEvolutionOption == GammaDriving2ndOrder_Covariant not supported!"
            else 1
        )

    cases = (
        (
            "Cartesian",
            False,
            False,
            "OnePlusLog",
            "GammaDriving2ndOrder_Covariant__Hatted",
            "Cartesian",
        ),
        (
            "Cartesian",
            False,
            False,
            "OnePlusLog",
            "NonAdvectingGammaDriving",
            "Cartesian",
        ),
        (
            "Cartesian",
            False,
            False,
            "OnePlusLog",
            "Frozen",
            "Cartesian",
        ),
        (
            "Cartesian",
            False,
            False,
            "Frozen",
            "GammaDriving2ndOrder_Covariant__Hatted",
            "Cartesian",
        ),
        (
            "Cartesian",
            False,
            False,
            "Frozen",
            "NonAdvectingGammaDriving",
            "Cartesian",
        ),
        (
            "Cartesian",
            False,
            False,
            "Frozen",
            "Frozen",
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
        (
            "SinhSpherical",
            False,
            False,
            "OnePlusLog",
            "NonAdvectingGammaDriving",
            "SinhSpherical",
        ),
        (
            "SinhSpherical",
            False,
            False,
            "OnePlusLog",
            "Frozen",
            "SinhSpherical",
        ),
        (
            "SinhSpherical",
            False,
            False,
            "Frozen",
            "GammaDriving2ndOrder_Covariant__Hatted",
            "SinhSpherical",
        ),
        (
            "SinhSpherical",
            False,
            False,
            "Frozen",
            "NonAdvectingGammaDriving",
            "SinhSpherical",
        ),
        (
            "SinhSpherical",
            False,
            False,
            "Frozen",
            "Frozen",
            "SinhSpherical",
        ),
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
        (
            validation_bssn_alpha_rhs,
            validation_bssn_vet_rhsU,
            validation_bssn_bet_rhsU,
        ) = BSSN_gauge_RHSs(
            CoordSystem=case_coord,
            enable_rfm_precompute=case_enable_rfm_precompute,
            enable_T4munu=case_enable_T4munu,
            LapseEvolutionOption=case_lapse,
            ShiftEvolutionOption=case_shift,
        )
        validation_suffix = case_coord + (
            "_rfm_precompute" if case_enable_rfm_precompute else ""
        )
        validation_Bq = BSSN_quantities[validation_suffix]
        validation_rhs = fCCZ4_RHSs[
            validation_suffix + ("_T4munu" if case_enable_T4munu else "")
        ]
        expected_alpha_rhs = validation_bssn_alpha_rhs
        if case_lapse == "OnePlusLog":
            expected_alpha_rhs += 4 * validation_Bq.alpha * validation_rhs.Theta
        expected_bet_rhsU = [validation_bssn_bet_rhsU[i] for i in range(3)]
        if case_shift in (
            "GammaDriving2ndOrder_Covariant__Hatted",
            "NonAdvectingGammaDriving",
        ):
            validation_rfm = refmetric.reference_metric[
                case_coord + ("_rfm_precompute" if case_enable_rfm_precompute else "")
            ]
            for validation_i in range(3):
                expected_bet_rhsU[validation_i] += (
                    sp.Rational(3, 4)
                    * validation_rhs.Lambdatilde_rhsU_delta[validation_i]
                    / validation_rfm.ReU[validation_i]
                )
        expression_dict = {
            "alpha_rhs": case_alpha_rhs,
            "bet_rhsU": case_bet_rhsU,
            "vet_rhsU": case_vet_rhsU,
            "alpha_correction_residual": case_alpha_rhs - expected_alpha_rhs,
            "bet_correction_residual": [
                case_bet_rhsU[i] - expected_bet_rhsU[i] for i in range(3)
            ],
            "vet_correction_residual": [
                case_vet_rhsU[i] - validation_bssn_vet_rhsU[i] for i in range(3)
            ],
        }
        if (
            case_coord == "Cartesian"
            and case_lapse == "OnePlusLog"
            and case_shift == "Frozen"
        ):
            expression_dict["unsupported_lapse_residual"] = unsupported_lapse_residual
            expression_dict["unsupported_shift_residual"] = unsupported_shift_residual
        processed = ve.process_dictionary_of_expressions(
            expression_dict, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{case_lapse}_{case_shift}_{trusted_coord_suffix}",
            processed,
        )
