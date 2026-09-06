# nrpy/infrastructures/Dendro/general_relativity/BSSN/rhs_eval.py
"""
Direct finite-difference BSSN right-hand side for the Dendro infrastructure.

This module assembles the BSSN evolution system the way ETLegacy and BHaH do —
:class:`nrpy.equations.general_relativity.BSSN_RHSs.BSSNRHSs` for the non-gauge
equations, :func:`nrpy.equations.general_relativity.BSSN_gauge_RHSs.BSSN_gauge_RHSs`
for the lapse and shift, and the shared Kreiss-Oliger helper — and then lowers
it through the same formulation-agnostic Dendro helpers the fCCZ4 builder uses.

It is the second formulation lowered through this infrastructure, and it exists
partly as the test of whether the generic layer really is generic.  Adding it
did require generic-layer work — the formulation-agnostic lowering was
extracted into ``kernel_lowering`` and ``tensor_family_of`` moved into
``naming`` — but no formulation-specific content entered the generic layer and
no existing emitter changed behaviour: the fCCZ4 output is byte-identical
across every profile.

:func:`build_bssn_rhs` registers no CFunction; it records the padding and the
upwind-control set it derived into the Dendro area of the NRPy parameter
registry, and :func:`register_CFunctions_rhs_eval` registers the bodies it
returns.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from collections import OrderedDict
from dataclasses import dataclass
from typing import Dict, Tuple

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_gauge_RHSs import BSSN_gauge_RHSs
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs
from nrpy.equations.general_relativity.kreiss_oliger_terms import (
    add_KreissOliger_dissipation_terms,
)
from nrpy.infrastructures.Dendro import generation_parameters  # noqa: F401
from nrpy.infrastructures.Dendro import kernel_lowering as kl
from nrpy.infrastructures.Dendro import naming
from nrpy.infrastructures.Dendro import registration as reg
from nrpy.infrastructures.Dendro.block_loop import block_loop
from nrpy.infrastructures.Dendro.simple_loop import require_serial_parallelization

# Dendro names solver sources for the formulation (its own BSSN solver ships
# bssnCtx.cpp and bssneqs.cpp), so the emitted CFunctions carry the same stem.
RHS_BLOCK_CFUNCTION = "bssn_rhs_block"
RHS_GLOBAL_CFUNCTION = "bssn_rhs"
FLAT_BLOCK_CFUNCTION = "bssn_rhs_flat_block"

# The BSSN evolved state: hDD (6), aDD (6), cf, trK, lambdaU (3), alpha,
# vetU (3), betU (3).
BSSN_EVOL_COUNT = 24


@dataclass(frozen=True)
class BSSNRHSBuild:
    """
    Immutable result of building the direct-FD BSSN RHS for one profile.

    :param evol_order: The 24 EVOL gridfunction names, in registry order.
    :param upwind_control_fields: The EVOL names the emitted kernel upwinds on.
    :param lvalues: The 24 output lvalues (``rhs_<name>[pp]``).
    :param padding: Ghost points (px, py, pz) the emitted operators reach.
    :param block_body: The per-block CFunction body (point loop + bindings).
    :param block_params: The per-block CFunction parameter list.
    :param global_body: The all-block CFunction body (NRPy block loop).
    :param global_params: The all-block CFunction parameter list.
    :param flat_block_body: The flat-block adapter CFunction body.
    :param flat_block_params: The flat-block adapter CFunction parameter list.
    """

    evol_order: Tuple[str, ...]
    upwind_control_fields: Tuple[str, ...]
    lvalues: Tuple[str, ...]
    padding: Tuple[int, int, int]
    block_body: str
    block_params: str
    global_body: str
    global_params: str
    flat_block_body: str
    flat_block_params: str


def bssn_rhs_expressions(
    *,
    CoordSystem: str,
    LapseEvolutionOption: str,
    ShiftEvolutionOption: str,
    enable_KreissOliger_dissipation: bool,
) -> Tuple[Dict[str, sp.Expr], Tuple[sp.Expr, ...]]:
    """
    Assemble the BSSN RHS expression set and its upwind control vector.

    The assembly order matches ETLegacy's and BHaH's: the non-gauge RHSs come
    from the cached ``BSSN_RHSs`` object, the gauge RHSs are added to a copy of
    its dictionary, Kreiss-Oliger terms are applied through the shared helper,
    and the upwind control vector is the rescaled shift.

    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :return: (rhs_by_symbol_name, upwind_control_vec).
    """
    rhs = BSSN_RHSs[CoordSystem]
    alpha_rhs, vet_rhsU, bet_rhsU = BSSN_gauge_RHSs(
        CoordSystem=CoordSystem,
        enable_rfm_precompute=False,
        enable_T4munu=False,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    # Copy before adding: BSSN_RHSs caches its dictionary, and mutating it
    # would leak this profile's gauge and dissipation choices into every later
    # consumer in the same process.
    rhs_by_symbol_name: Dict[str, sp.Expr] = OrderedDict(
        sorted(rhs.BSSN_RHSs_varname_to_expr_dict.items())
    )
    rhs_by_symbol_name["alpha_rhs"] = alpha_rhs
    for i in range(3):
        rhs_by_symbol_name[f"vet_rhsU{i}"] = vet_rhsU[i]
        rhs_by_symbol_name[f"bet_rhsU{i}"] = bet_rhsU[i]
    rhs_by_symbol_name = OrderedDict(sorted(rhs_by_symbol_name.items()))

    if enable_KreissOliger_dissipation:
        add_KreissOliger_dissipation_terms(
            rhs_by_symbol_name,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=False,
            registering_module=__name__,
            ShiftEvolutionOption=ShiftEvolutionOption,
            KreissOliger_strength_gauge=0.3,
            KreissOliger_strength_nongauge=0.3,
            enable_CAKO=False,
            # W is consumed only under enable_CAKO, which this profile does
            # not offer, so no conformal factor needs deriving here.
            W=sp.sympify(1),
            # BSSN has no Z4 scalar; the fCCZ4 profile is the only one that
            # dissipates Theta_fCCZ4.
            include_Theta_fCCZ4=False,
        )

    # Upwind control vector is the rescaled shift, exactly as ETLegacy builds
    # it: betaU[i] = vetU[i] * ReU[i].
    rfm = refmetric.reference_metric[CoordSystem]
    vetU = ixp.declarerank1("vetU")
    betaU = ixp.zerorank1()
    for i in range(3):
        betaU[i] = vetU[i] * rfm.ReU[i]
    return rhs_by_symbol_name, tuple(betaU)


def build_bssn_rhs(
    *,
    fd_order: int,
    enable_KreissOliger_dissipation: bool,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
) -> BSSNRHSBuild:
    """
    Build the direct-FD BSSN RHS for one profile.

    The caller must set ``Infrastructure`` to ``Dendro`` and ``parallelization``
    to ``"none"`` before calling; both are asserted here rather than
    overwritten, because silently discarding a caller's request would produce
    an unqualified configuration.

    :param fd_order: The finite-difference order (2, 4, or 6).  Order 8
        reaches five ghost points, above the ``max_proven_padding`` of 4 in
        ``dendrolib_capabilities.json``, so it stays capability-gated.
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :return: The immutable :class:`BSSNRHSBuild` result.
    :raises ValueError: If ``Infrastructure`` is not Dendro, if ``fd_order`` is
        outside the qualified set, if the RHS symbols do not map bijectively
        onto the registered EVOL fields, if the EVOL count is not exactly 24,
        or if the kernel's dKOD presence does not match
        ``enable_KreissOliger_dissipation``.

    Doctests:
    >>> import contextlib, io
    >>> import nrpy.infrastructures.Dendro.generation_parameters  # noqa: F401
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fp_type", "double")
    >>> par.set_parval_from_str("Dendro_scalar_type", "DendroScalar")
    >>> par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    >>> try:
    ...     build_bssn_rhs(fd_order=8, enable_KreissOliger_dissipation=False)
    ... except ValueError as error:
    ...     print(str(error).splitlines()[0])
    Unsupported fd_order=8; allowed: (2, 4, 6). fd_order 8 reaches five ghost points, above the max_proven_padding of 4 recorded in dendrolib_capabilities.json, so it stays capability-gated until Dendrolib is pinned.
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _build = build_bssn_rhs(fd_order=4, enable_KreissOliger_dissipation=False)
    >>> len(_build.evol_order), "Theta_fCCZ4" in _build.evol_order
    (24, False)
    >>> sorted(_build.lvalues)[:2]
    ['rhs_aDD00[pp]', 'rhs_aDD01[pp]']
    >>> _build.padding
    (3, 3, 3)
    >>> _build.upwind_control_fields
    ('vetU0', 'vetU1', 'vetU2')
    >>> reg.upwind_control_fields() == _build.upwind_control_fields
    True
    >>> reg.required_padding()
    (3, 3, 3)
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro BSSN RHS, got "
            f"{par.parval_from_str('Infrastructure')!r}."
        )
    if fd_order not in (2, 4, 6):
        raise ValueError(
            f"Unsupported fd_order={fd_order!r}; allowed: (2, 4, 6). "
            "fd_order 8 reaches five ghost points, above the "
            "max_proven_padding of 4 recorded in dendrolib_capabilities.json, "
            "so it stays capability-gated until Dendrolib is pinned."
        )
    par.set_parval_from_str(
        "Dendro_enable_KreissOliger_dissipation", enable_KreissOliger_dissipation
    )
    generation_parameters.validate_generation_parameters()
    require_serial_parallelization()

    rhs_by_symbol_name, upwind_control_vec = bssn_rhs_expressions(
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    )
    rhs_symbols = tuple(rhs_by_symbol_name)
    lvalues = tuple(
        f"rhs_{naming.rhs_symbol_to_gridfunction_name(symbol)}[pp]"
        for symbol in rhs_symbols
    )
    evol_order = reg.registered_evol_order()
    if len(evol_order) != BSSN_EVOL_COUNT:
        raise ValueError(
            f"BSSN EVOL registry must hold exactly {BSSN_EVOL_COUNT} fields, "
            f"found {len(evol_order)}."
        )
    mapped = {naming.rhs_symbol_to_gridfunction_name(symbol) for symbol in rhs_symbols}
    if len(lvalues) != len(set(lvalues)) or mapped != set(evol_order):
        raise ValueError(
            "BSSN RHS symbols do not map bijectively onto the registered EVOL "
            f"fields: missing={sorted(set(evol_order) - mapped)} "
            f"extra={sorted(mapped - set(evol_order))}"
        )
    upwind_control_fields = kl.upwind_control_fields(upwind_control_vec, evol_order)

    par.set_parval_from_str("fd_order", fd_order)
    fp_type = str(par.parval_from_str("fp_type"))
    fp_type_alias = str(par.parval_from_str("Dendro_scalar_type"))
    kernel = c_codegen(
        list(rhs_by_symbol_name.values()),
        list(lvalues),
        enable_fd_codegen=True,
        enable_fd_functions=False,
        enable_simd=False,
        fp_type=fp_type,
        fp_type_alias=fp_type_alias,
        mem_alloc_style="210",
        rational_const_alias="static const",
        verbose=False,
        upwind_control_vec=list(upwind_control_vec),
    )
    operator_records = kl.emitted_operators(kernel, fd_order)
    padding = kl.padding_from_operators(operator_records)
    if (
        any(record["operator"].startswith("dKOD") for record in operator_records)
        != enable_KreissOliger_dissipation
    ):
        raise ValueError(
            "dKOD-operator presence does not match "
            f"enable_KreissOliger_dissipation={enable_KreissOliger_dissipation!r}."
        )

    used_codeparameters = kl.used_codeparameters(rhs_by_symbol_name.values())
    cparam_args = kl.cparam_declarations(used_codeparameters, fp_type_alias)
    cparam_values = kl.cparam_arguments(used_codeparameters)
    tail_args = f", {cparam_args}" if cparam_args else ""
    tail_values = f", {cparam_values}" if cparam_values else ""

    block_body = (
        kl.block_pointer_bindings(evol_order, fp_type_alias)
        + "\n"
        + kl.point_loop(kernel)
    )
    block_params = (
        f"const BlockGeometry& geom, const {fp_type_alias}* const* in_gfs, "
        f"{fp_type_alias}* const* rhs_gfs" + tail_args
    )
    global_body = block_loop(
        f"{RHS_BLOCK_CFUNCTION}(world.geom[blk], in_gfs, rhs_gfs{tail_values});",
        num_blocks="world.num_blocks",
    )
    global_params = (
        f"const MockWorld& world, const {fp_type_alias}* const* in_gfs, "
        f"{fp_type_alias}* const* rhs_gfs" + tail_args
    )
    # Single-body rule: the flat adapter derives flat-layout pointers, packs
    # the per-component pointer arrays, and calls the registered block kernel.
    # There is exactly one numerical body, so the two paths cannot diverge.
    flat_block_body = (
        kl.flat_block_pointer_bindings(evol_order, fp_type_alias)
        + f"\nconst {fp_type_alias}* const in_gfs_call[] = {{"
        + ", ".join(naming.input_pointer(name) for name in evol_order)
        + f"}};\n{fp_type_alias}* rhs_gfs_call[] = {{"
        + ", ".join(naming.rhs_pointer(name) for name in evol_order)
        + "};\n"
        + f"{RHS_BLOCK_CFUNCTION}(geom, in_gfs_call, rhs_gfs_call{tail_values});\n"
    )
    flat_block_params = (
        f"const BlockGeometry& geom, const {fp_type_alias}* const in_gfs_flat, "
        f"{fp_type_alias}* const rhs_gfs_flat" + tail_args
    )

    reg.set_upwind_control_fields(upwind_control_fields)
    reg.set_required_padding(padding)
    return BSSNRHSBuild(
        evol_order=evol_order,
        upwind_control_fields=upwind_control_fields,
        lvalues=lvalues,
        padding=padding,
        block_body=block_body,
        block_params=block_params,
        global_body=global_body,
        global_params=global_params,
        flat_block_body=flat_block_body,
        flat_block_params=flat_block_params,
    )


def register_CFunctions_rhs_eval(
    *,
    fd_order: int,
    enable_KreissOliger_dissipation: bool,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
) -> None:
    """
    Register the per-block, all-block, and flat-block BSSN RHS CFunctions.

    :param fd_order: The finite-difference order (2, 4, or 6).
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    """
    build = build_bssn_rhs(
        fd_order=fd_order,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    subdirectory = "generated/src/rhs"
    reg.register_Dendro_CFunction(
        role="rhs_block",
        name=RHS_BLOCK_CFUNCTION,
        desc="Per-block direct-FD BSSN RHS (24 fields).",
        subdirectory=subdirectory,
        params=build.block_params,
        body=build.block_body,
    )
    reg.register_Dendro_CFunction(
        role="rhs",
        name=RHS_GLOBAL_CFUNCTION,
        desc="All-block direct-FD BSSN RHS (NRPy block loop).",
        subdirectory=subdirectory,
        params=build.global_params,
        body=build.global_body,
    )
    reg.register_Dendro_CFunction(
        role="rhs_flat_block",
        name=FLAT_BLOCK_CFUNCTION,
        desc="LTS flat-block adapter (same numerical body, flat layout).",
        subdirectory=subdirectory,
        params=build.flat_block_params,
        body=build.flat_block_body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
