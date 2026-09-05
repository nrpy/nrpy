# nrpy/equations/general_relativity/kreiss_oliger_terms.py
"""
Shared Kreiss-Oliger dissipation terms for the BSSN and fCCZ4 RHS systems.

The two formulations add the same dissipation terms to the same evolved
quantities; fCCZ4 adds one more, for its ``Theta_fCCZ4`` field.  Keeping one
implementation here means a change to the dissipation (a strength, a
curvature-aware factor, a new field) is made once instead of twice.

The registering module is a parameter, so each caller keeps registering the
strengths under its own module name and no caller's generated output moves.

Author: NRPy Dendro fCCZ4 infrastructure (PR 4)
"""

from typing import Dict

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric


def add_KreissOliger_dissipation_terms(
    rhs_by_symbol_name: Dict[str, sp.Expr],
    *,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    registering_module: str,
    ShiftEvolutionOption: str,
    KreissOliger_strength_gauge: float,
    KreissOliger_strength_nongauge: float,
    enable_CAKO: bool,
    W: sp.Expr,
    include_Theta_fCCZ4: bool,
) -> None:
    """
    Add Kreiss-Oliger dissipation terms to an RHS dictionary, in place.

    Gauge quantities (lapse, shift, shift driver) receive the gauge strength;
    geometric and Z4 quantities receive the non-gauge strength.  With
    ``enable_CAKO`` both strengths are multiplied by the conformal factor W.

    :param rhs_by_symbol_name: RHS dictionary, keyed by NRPy RHS symbol name;
        mutated in place.
    :param CoordSystem: Reference-metric coordinate system.
    :param enable_rfm_precompute: Use precomputed reference-metric symbols.
    :param registering_module: Module name under which the two strength
        CodeParameters are registered (the caller's ``__name__``, so generated
        parameter metadata is unchanged by this extraction).
    :param ShiftEvolutionOption: Shift evolution option; the second-order
        shift driver has its own evolved field to dissipate.
    :param KreissOliger_strength_gauge: Default gauge dissipation strength.
    :param KreissOliger_strength_nongauge: Default non-gauge strength.
    :param enable_CAKO: Multiply both strengths by the conformal factor W.
    :param W: The conformal factor, in the caller's chosen representation.
    :param include_Theta_fCCZ4: Also dissipate the fCCZ4 ``Theta_fCCZ4``
        field, which BSSN does not have.
    """
    diss_strength_gauge, diss_strength_nongauge = par.register_CodeParameters(
        "REAL",
        registering_module,
        ["KreissOliger_strength_gauge", "KreissOliger_strength_nongauge"],
        [KreissOliger_strength_gauge, KreissOliger_strength_nongauge],
        commondata=True,
    )

    # vvv BEGIN CAKO vvv
    if enable_CAKO:
        diss_strength_gauge *= W
        diss_strength_nongauge *= W
    # ^^^ END CAKO ^^^

    rfm = refmetric.reference_metric[
        (CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem)
    ]
    alpha_dKOD = ixp.declarerank1("alpha_dKOD")
    cf_dKOD = ixp.declarerank1("cf_dKOD")
    trK_dKOD = ixp.declarerank1("trK_dKOD")
    betU_dKOD = ixp.declarerank2("betU_dKOD", symmetry="nosym")
    vetU_dKOD = ixp.declarerank2("vetU_dKOD", symmetry="nosym")
    lambdaU_dKOD = ixp.declarerank2("lambdaU_dKOD", symmetry="nosym")
    aDD_dKOD = ixp.declarerank3("aDD_dKOD", symmetry="sym01")
    hDD_dKOD = ixp.declarerank3("hDD_dKOD", symmetry="sym01")
    if include_Theta_fCCZ4:
        Theta_fCCZ4_dKOD = ixp.declarerank1("Theta_fCCZ4_dKOD")
    for k in range(3):
        rhs_by_symbol_name["alpha_rhs"] += (
            diss_strength_gauge * alpha_dKOD[k] * rfm.ReU[k]
        )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
        rhs_by_symbol_name["cf_rhs"] += (
            diss_strength_nongauge * cf_dKOD[k] * rfm.ReU[k]
        )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
        rhs_by_symbol_name["trK_rhs"] += (
            diss_strength_nongauge * trK_dKOD[k] * rfm.ReU[k]
        )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
        if include_Theta_fCCZ4:
            rhs_by_symbol_name["Theta_fCCZ4_rhs"] += (
                diss_strength_nongauge * Theta_fCCZ4_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
        for i in range(3):
            if "2ndOrder" in ShiftEvolutionOption:
                rhs_by_symbol_name[f"bet_rhsU{i}"] += (
                    diss_strength_gauge * betU_dKOD[i][k] * rfm.ReU[k]
                )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs_by_symbol_name[f"vet_rhsU{i}"] += (
                diss_strength_gauge * vetU_dKOD[i][k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs_by_symbol_name[f"lambda_rhsU{i}"] += (
                diss_strength_nongauge * lambdaU_dKOD[i][k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            for j in range(i, 3):
                rhs_by_symbol_name[f"a_rhsDD{i}{j}"] += (
                    diss_strength_nongauge * aDD_dKOD[i][j][k] * rfm.ReU[k]
                )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                rhs_by_symbol_name[f"h_rhsDD{i}{j}"] += (
                    diss_strength_nongauge * hDD_dKOD[i][j][k] * rfm.ReU[k]
                )  # ReU[k] = 1/scalefactor_orthog_funcform[k]


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
