# nrpy/equations/general_relativity/fCCZ4_system.py
"""
Shared, infrastructure-neutral fCCZ4 expression factory.

This module assembles the complete fCCZ4 right-hand-side system -- the
non-gauge fCCZ4 equations, the moving-puncture gauge equations, and the
optional Kreiss-Oliger, CAHD, and SSL terms -- into one infrastructure-
neutral expression bundle.  Both the BHaH and Dendro infrastructure
registrars consume this bundle, so there is exactly one fCCZ4 expression
source.

The factory returns expressions and semantic relationships only: it does not
emit a state order, infrastructure-specific names, parameter default table,
stencil plan, loop, or CFunction.

Author: NRPy Dendro fCCZ4 infrastructure (PR 4)
"""

from collections import OrderedDict as ODict
from dataclasses import dataclass
from typing import Dict, Mapping, Tuple

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.fCCZ4_constraints import FCCZ4Constraints
from nrpy.equations.general_relativity.fCCZ4_gauge_RHSs import fCCZ4_gauge_RHSs
from nrpy.equations.general_relativity.fCCZ4_RHSs import fCCZ4_RHSs
from nrpy.equations.general_relativity.kreiss_oliger_terms import (
    add_KreissOliger_dissipation_terms,
)


@dataclass(frozen=True)
class FCCZ4ExpressionBundle:
    """
    Immutable fCCZ4 expression bundle for one formulation profile.

    :param rhs_by_symbol_name: Mapping from NRPy RHS symbol name
        (e.g., ``"a_rhsDD00"``, ``"alpha_rhs"``) to its full symbolic
        right-hand side, gauge and optional terms included.
    :param upwind_control_vec: Upwind control vector (beta^i) used by NRPy
        finite-difference upwinding.
    :param projection_by_gridfunction_name: Algebraic projection expressions
        keyed by EVOL gridfunction name.  Empty until the PR8/9 projection
        profile qualifies its own closure, fingerprint, hashes, and tests
        (whitepaper section 10.4: each feature is a separate profile).
    :param diagnostics_by_name: Diagnostic expressions keyed by diagnostic
        name (``H_Z4`` and ``Z4constraintU0..2`` from the registered
        constraint factory).  Empty unless ``enable_diagnostics`` was set.
    :param provenance_files: Basenames of the equation modules the bundle was
        assembled from, in deterministic order, derived from the actual
        enabled inputs.
    """

    rhs_by_symbol_name: Mapping[str, sp.Expr]
    upwind_control_vec: Tuple[sp.Expr, sp.Expr, sp.Expr]
    projection_by_gridfunction_name: Mapping[str, sp.Expr]
    diagnostics_by_name: Mapping[str, sp.Expr]
    provenance_files: Tuple[str, ...]


def build_fccz4_expression_bundle(
    *,
    CoordSystem: str = "Cartesian",
    enable_rfm_precompute: bool = False,
    enable_RbarDD_gridfunctions: bool = False,
    enable_T4munu: bool = False,
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    enable_KreissOliger_dissipation: bool = False,
    KreissOliger_strength_gauge: float = 0.3,
    KreissOliger_strength_nongauge: float = 0.3,
    enable_CAKO: bool = False,
    enable_CAHD: bool = False,
    enable_SSL: bool = False,
    enable_YBS_Gamma_constraint_adjustment: bool = False,
    enable_YBS_momentum_constraint_adjustment: bool = False,
    enable_diagnostics: bool = False,
) -> FCCZ4ExpressionBundle:
    """
    Assemble the complete fCCZ4 RHS expression system for one profile.

    The assembly sequence is fixed:

    1. Build (or fetch from cache) the non-gauge fCCZ4 RHSs.
    2. Copy the returned non-gauge dictionary before adding terms.
    3. Build the selected lapse and shift gauge RHSs and add
       ``alpha_rhs``, ``vet_rhsU0..2``, and ``bet_rhsU0..2``.
    4. Add Kreiss-Oliger, CAHD, and SSL terms when enabled.
    5. Build the upwind control vector (beta^i).

    :param CoordSystem: Reference-metric coordinate system.
    :param enable_rfm_precompute: Use precomputed reference-metric symbols.
    :param enable_RbarDD_gridfunctions: Read the BSSN-shaped conformal Ricci
        components from auxiliary gridfunctions.
    :param enable_T4munu: Include stress-energy source terms.
    :param LapseEvolutionOption: Lapse evolution equation choice.
    :param ShiftEvolutionOption: Shift evolution equation choice.
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :param KreissOliger_strength_gauge: KO strength for gauge fields.
    :param KreissOliger_strength_nongauge: KO strength for geometry/Z4 fields.
    :param enable_CAKO: Multiply KO strengths by the conformal factor W.
    :param enable_CAHD: Enable curvature-aware Hamiltonian-constraint damping.
    :param enable_SSL: Enable slow-start lapse.
    :param enable_YBS_Gamma_constraint_adjustment: Enable the YBS
        connection-constraint adjustment.
    :param enable_YBS_momentum_constraint_adjustment: Enable the
        Yo--Lin--Cao momentum-constraint adjustment.
    :param enable_diagnostics: Also build the constraint diagnostics
        (``H_Z4`` and ``Z4constraintU0..2``).  Off by default: constructing
        the constraint factory registers auxiliary gridfunctions in some
        profiles, and no RHS consumer needs the result.
    :return: The assembled fCCZ4 expression bundle.
    :raises ValueError: If EvolvedConformalFactor_cf is not a supported value,
        or if an optional feature combination is unsupported (whitepaper
        section 10.2 step 8: curvature-aware KO requires KO dissipation).

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.grid as gri
    >>> par.set_parval_from_str("Infrastructure", "BHaH")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     bundle = build_fccz4_expression_bundle()

    Whitepaper section 5.5: the fCCZ4 state is exactly 25 evolved fields, and
    the Z4 scalar is one of them.

    >>> evol, _auxevol, _diag, _aux = gri.GridFunction.gridfunction_lists()
    >>> len(evol)
    25
    >>> "Theta_fCCZ4" in evol
    True
    >>> len(bundle.rhs_by_symbol_name)
    25

    Section 8.5, single Kreiss-Oliger ownership: dissipation terms appear if
    and only if they were requested.

    >>> def _has_ko(b):
    ...     return any(
    ...         "dKOD" in str(symbol)
    ...         for expr in b.rhs_by_symbol_name.values()
    ...         for symbol in expr.free_symbols
    ...     )
    >>> _has_ko(bundle)
    False
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     ko_bundle = build_fccz4_expression_bundle(
    ...         enable_KreissOliger_dissipation=True
    ...     )
    >>> _has_ko(ko_bundle)
    True

    Curvature-aware KO requires KO dissipation (section 10.2 step 8).

    >>> try:
    ...     build_fccz4_expression_bundle(enable_CAKO=True)
    ... except ValueError:
    ...     print("Rejected without KO dissipation. Good.")
    Rejected without KO dissipation. Good.

    Diagnostics are opt-in (section 10.4): the RHS consumers do not need them.

    >>> bundle.diagnostics_by_name
    {}
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     diag_bundle = build_fccz4_expression_bundle(enable_diagnostics=True)
    >>> sorted(diag_bundle.diagnostics_by_name)
    ['H_Z4', 'Z4constraintU0', 'Z4constraintU1', 'Z4constraintU2']
    """
    if enable_CAKO and not enable_KreissOliger_dissipation:
        raise ValueError(
            "enable_CAKO requires enable_KreissOliger_dissipation: "
            "curvature-aware KO multiplies the KO strengths."
        )
    rhs_cache_key = (
        CoordSystem
        + ("_rfm_precompute" if enable_rfm_precompute else "")
        + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
        + ("_T4munu" if enable_T4munu else "")
    )

    # Step 1+2: Copy the cached non-gauge fCCZ4 RHS dictionary before adding
    # terms; gauge and optional KO/CAHD/SSL terms belong only to this bundle.
    fccz4_rhs = fCCZ4_RHSs.get_rhs(
        rhs_cache_key,
        enable_YBS_Gamma_constraint_adjustment=(enable_YBS_Gamma_constraint_adjustment),
        enable_YBS_momentum_constraint_adjustment=(
            enable_YBS_momentum_constraint_adjustment
        ),
    )
    local_RHSs_varname_to_expr_dict = fccz4_rhs.fCCZ4_RHSs_varname_to_expr_dict.copy()

    # Step 3: Build the selected lapse and shift gauge RHSs.
    alpha_rhs, vet_rhsU, bet_rhsU = fCCZ4_gauge_RHSs(
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_T4munu=enable_T4munu,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
        enable_YBS_Gamma_constraint_adjustment=(enable_YBS_Gamma_constraint_adjustment),
    )
    local_RHSs_varname_to_expr_dict["alpha_rhs"] = alpha_rhs
    for i in range(3):
        local_RHSs_varname_to_expr_dict[f"vet_rhsU{i}"] = vet_rhsU[i]
        local_RHSs_varname_to_expr_dict[f"bet_rhsU{i}"] = bet_rhsU[i]

    local_RHSs_varname_to_expr_dict = ODict(
        sorted(local_RHSs_varname_to_expr_dict.items())
    )

    # Step 4.a: Define conformal factor W.
    Bq = BSSN_quantities[
        CoordSystem
        + ("_rfm_precompute" if enable_rfm_precompute else "")
        + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
    ]
    EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
    if EvolvedConformalFactor_cf == "W":
        W = Bq.cf
    elif EvolvedConformalFactor_cf == "chi":
        W = sp.sqrt(Bq.cf)
    elif EvolvedConformalFactor_cf == "phi":
        W = sp.exp(-2 * Bq.cf)
    else:
        raise ValueError(
            "Error: only EvolvedConformalFactor_cf = (W or chi or phi) supported."
        )

    # Step 4.b: Add Kreiss-Oliger dissipation to the fCCZ4 RHSs.  The terms
    # are shared with the BSSN assembly; fCCZ4 additionally dissipates
    # Theta_fCCZ4 (whitepaper section 10.2 step 7).
    if enable_KreissOliger_dissipation:
        add_KreissOliger_dissipation_terms(
            local_RHSs_varname_to_expr_dict,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=enable_rfm_precompute,
            registering_module=__name__,
            ShiftEvolutionOption=ShiftEvolutionOption,
            KreissOliger_strength_gauge=KreissOliger_strength_gauge,
            KreissOliger_strength_nongauge=KreissOliger_strength_nongauge,
            enable_CAKO=enable_CAKO,
            W=W,
            include_Theta_fCCZ4=True,
        )

    # Step 4.c: Add the curvature-aware Hamiltonian-constraint damping.
    if enable_CAHD:
        Bcon = BSSN_constraints[
            CoordSystem
            + ("_rfm_precompute" if enable_rfm_precompute else "")
            + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
            + ("_T4munu" if enable_T4munu else "")
        ]
        C_CAHD = par.register_CodeParameter(
            "REAL", __name__, "C_CAHD", 0.15, commondata=True, add_to_parfile=True
        )
        # Initialize CAHD_term assuming phi is the evolved conformal factor.
        # CFL_FACTOR is defined in MoL; dsmin stores raw physical grid spacing.
        CAHD_term = (
            -C_CAHD
            * sp.Symbol("CFL_FACTOR", real=True)
            * sp.Symbol("dsmin", real=True)
            * Bcon.H
        )
        if EvolvedConformalFactor_cf == "phi":
            pass  # CAHD_term already assumes phi is the evolved conformal factor.
        elif EvolvedConformalFactor_cf == "W":
            # \partial_t W = \partial_t e^{-2 phi} = -2 W \partial_t phi
            CAHD_term *= -2 * Bq.cf
        elif EvolvedConformalFactor_cf == "chi":
            # \partial_t chi = \partial_t e^{-4 phi} = -4 chi \partial_t phi
            CAHD_term *= -4 * Bq.cf
        else:
            raise ValueError(
                "Error: only EvolvedConformalFactor_cf = (W or chi or phi) supported."
            )
        local_RHSs_varname_to_expr_dict["cf_rhs"] += CAHD_term

    # Step 4.d: Add slow-start lapse.
    if enable_SSL:
        SSL_Gaussian_prefactor = par.register_CodeParameter(
            "REAL",
            __name__,
            "SSL_Gaussian_prefactor",
            1.0,
            commondata=True,
            add_to_parfile=False,
        )
        _SSL_h, _SSL_sigma = par.register_CodeParameters(
            "REAL",
            __name__,
            ["SSL_h", "SSL_sigma"],
            [0.6, 20.0],
            commondata=True,
            add_to_parfile=True,
        )
        local_RHSs_varname_to_expr_dict["alpha_rhs"] -= (
            W * SSL_Gaussian_prefactor * (Bq.alpha - W)
        )

    # Step 5: Build the upwind control vector (beta^i).
    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    betaU = ixp.zerorank1()
    vetU = ixp.declarerank1("vetU")
    for i in range(3):
        # self.lambda_rhsU[i] = self.Lambdabar_rhsU[i] / rfm.ReU[i]
        betaU[i] = vetU[i] * rfm.ReU[i]

    # Step 6: Diagnostics from the registered constraint factory (whitepaper
    # section 10.1: H_Z4 and Z4constraintU0..2 first).  Building them is
    # opt-in: constructing the constraint factory is a side effect on the
    # caller's registries (with enable_T4munu it registers auxiliary
    # gridfunctions), and the RHS consumers do not use the result.  The
    # diagnostics profile qualifies separately in PR 9 (section 10.4).
    diagnostics: Dict[str, sp.Expr] = {}
    if enable_diagnostics:
        constraints = FCCZ4Constraints(
            CoordSystem=CoordSystem,
            enable_rfm_precompute=enable_rfm_precompute,
            enable_RbarDD_gridfunctions=enable_RbarDD_gridfunctions,
            enable_T4munu=enable_T4munu,
        )
        diagnostics["H_Z4"] = constraints.H_Z4
        for i in range(3):
            diagnostics[f"Z4constraintU{i}"] = constraints.Z4constraintU[i]

    # Provenance lists the equation modules this bundle was actually
    # assembled from, so it tracks the enabled inputs rather than a fixed
    # list.
    provenance = [
        "BSSN_quantities.py",
        "fCCZ4_RHSs.py",
        "fCCZ4_gauge_RHSs.py",
    ]
    if enable_KreissOliger_dissipation:
        provenance.append("kreiss_oliger_terms.py")
    if enable_diagnostics:
        provenance.append("fCCZ4_constraints.py")
    if enable_CAHD:
        provenance.append("BSSN_constraints.py")
    if enable_T4munu:
        provenance.append("T4munu.py")

    return FCCZ4ExpressionBundle(
        rhs_by_symbol_name=local_RHSs_varname_to_expr_dict,
        upwind_control_vec=(betaU[0], betaU[1], betaU[2]),
        projection_by_gridfunction_name={},
        diagnostics_by_name=diagnostics,
        provenance_files=tuple(sorted(provenance)),
    )


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Trusted-expression validation of the assembled bundle, over the option
    # axes the qualified profile varies: Kreiss-Oliger dissipation off and on.
    # This pins the assembled right-hand sides themselves, so a change in any
    # contributing equation module is caught here at the symbolic layer.
    for case_enable_KO, trusted_suffix in (
        (False, "Cartesian"),
        (True, "Cartesian_KO"),
    ):
        case_bundle = build_fccz4_expression_bundle(
            CoordSystem="Cartesian",
            enable_KreissOliger_dissipation=case_enable_KO,
        )
        processed = ve.process_dictionary_of_expressions(
            dict(case_bundle.rhs_by_symbol_name), fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{trusted_suffix}",
            processed,
        )
