"""
Direct finite-difference BSSN or fCCZ4 right-hand sides for Dendro.

One module, one entry point, and the formulation chosen by the
``enable_fCCZ4`` argument, as ``BHaH/general_relativity/rhs_eval.py`` does for
the same two formulations.  The two expression-assembly routes stay separate
below: fCCZ4 comes from the shared ``fCCZ4_system`` factory, BSSN from
``BSSN_RHSs`` plus the gauge terms and optional Kreiss-Oliger dissipation, and
BHaH's own divergence guard records that unifying those branches is not wanted.

Everything downstream of the expressions is shared: the point kernel is lowered
by ``c_codegen`` with direct finite differences, wrapped in the Dendro interior
point loop, and registered as a per-block, an all-block and a flat-block
CFunction.  This module authors no field name, no physics default and no
finite-difference coefficient.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from collections import OrderedDict
from dataclasses import dataclass
from typing import Dict, Mapping, Tuple, Union

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_gauge_RHSs import BSSN_gauge_RHSs
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs
from nrpy.equations.general_relativity.fCCZ4_system import (
    build_fccz4_expression_bundle,
)
from nrpy.equations.general_relativity.kreiss_oliger_terms import (
    add_KreissOliger_dissipation_terms,
)
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import block_kernel_helpers as bkh
from nrpy.infrastructures.Dendro import generation_parameters
from nrpy.infrastructures.Dendro import gridfunction_name_decorations as gf_names
from nrpy.infrastructures.Dendro.simple_loop import (
    block_loop,
    require_serial_parallelization,
)

# The per-block CFunction name (the Dendro scheduling role key).
FCCZ4_RHS_EVAL_BLOCK_CFUNCTION = "fccz4_rhs_eval_block"

# All-block CFunction name.
FCCZ4_RHS_EVAL_ALL_BLOCKS_CFUNCTION = "fccz4_rhs_eval"

# LTS flat-block adapter CFunction name (same numerical body, flat layout).
FCCZ4_RHS_EVAL_FLAT_BLOCK_CFUNCTION = "fccz4_rhs_eval_flat_block"


@dataclass(frozen=True)
class FCCZ4RHSBuild:
    """
    Immutable result of building the direct-FD fCCZ4 RHS for one profile.

    :param evol_order: The 25 EVOL gridfunction names, in registry order.
    :param upwind_control_fields: The EVOL names the emitted kernel upwinds on.
    :param lvalues: The 25 output lvalues (``rhs_<name>[pp]``).
    :param padding: Ghost points the emitted operators reach, widest axis.
    :param block_body: The per-block CFunction body (point loop + bindings).
    :param block_params: The per-block CFunction parameter list.
    :param all_blocks_body: The all-block CFunction body (NRPy block loop).
    :param all_blocks_params: The all-block CFunction parameter list.
    :param flat_block_body: The flat-block adapter CFunction body.
    :param flat_block_params: The flat-block adapter CFunction parameter list.
    :param used_codeparameters: The CodeParameters the signatures forward, in
        signature order, computed from the expression free symbols.
    :param rhs_by_symbol_name: The assembled symbolic right-hand sides, kept so
        the module's ``__main__`` can pin them against trusted values without
        reassembling them.
    """

    evol_order: Tuple[str, ...]
    upwind_control_fields: Tuple[str, ...]
    lvalues: Tuple[str, ...]
    padding: int
    block_body: str
    block_params: str
    all_blocks_body: str
    all_blocks_params: str
    flat_block_body: str
    flat_block_params: str
    used_codeparameters: Tuple[str, ...]
    rhs_by_symbol_name: Mapping[str, sp.Expr]


def _build_rhs_eval_fCCZ4(
    *,
    fd_order: int,
    enable_KreissOliger_dissipation: bool,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
) -> FCCZ4RHSBuild:
    """
    Build the direct-FD fCCZ4 RHS for one profile.

    This is a pure builder (no CFunction registration); the caller registers
    the returned bodies with :func:`nrpy.c_function.register_CFunction`
    and records each scheduling role with
    :func:`nrpy.infrastructures.Dendro.CFunction_roles.set_CFunction_role`.

    The caller must set ``Infrastructure`` to ``Dendro`` and
    ``parallelization`` to ``"none"`` before calling; both are asserted here
    rather than overwritten, because silently discarding a caller's request
    would produce an unqualified configuration.  The registered Dendro
    generation parameters are validated before anything is lowered.

    :param fd_order: The finite-difference order (2, 4, or 6).  Order 8
        reaches five ghost points, which the pinned Dendrolib supports at
        element order 10; it stays outside this builder's qualified set.
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation in the shared factory.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :return: The immutable :class:`FCCZ4RHSBuild` result.
    :raises ValueError: If ``fd_order`` is outside the qualified set, if the
        fCCZ4 RHS symbols do not map bijectively onto the registered EVOL
        fields, if the EVOL count is not exactly 25, or if the kernel's dKOD
        presence does not match ``enable_KreissOliger_dissipation``.

    Doctests:
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fp_type", "double")
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    >>> par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    >>> import nrpy.grid as _gri
    >>> _gri.glb_gridfcs_dict.clear()
    >>> par.glb_extras_dict.pop("Dendro", None) and None
    >>> try:
    ...     build_rhs_eval(enable_fCCZ4=True, fd_order=8, enable_KreissOliger_dissipation=False)
    ... except ValueError as error:
    ...     print(str(error).splitlines()[0])
    Unsupported fd_order=8; allowed: (2, 4, 6). fd_order 8 reaches five ghost points, which the pinned Dendrolib proves at element order 10; it is outside this builder's qualified set rather than host-gated.
    >>> import contextlib, io
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _build = build_rhs_eval(
    ...         enable_fCCZ4=True, fd_order=4, enable_KreissOliger_dissipation=False
    ...     )
    >>> len(_build.evol_order), "Theta_fCCZ4" in _build.evol_order
    (25, True)
    >>> sorted(_build.lvalues)[:2]
    ['rhs_Theta_fCCZ4[pp]', 'rhs_aDD00[pp]']
    >>> _build.padding
    3
    >>> _build.upwind_control_fields
    ('vetU0', 'vetU1', 'vetU2')
    >>> roles.upwind_control_fields()
    ('vetU0', 'vetU1', 'vetU2')
    >>> roles.required_padding()
    3
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro fCCZ4 RHS, got "
            f"{par.parval_from_str('Infrastructure')!r}."
        )
    if fd_order not in (2, 4, 6):
        raise ValueError(
            f"Unsupported fd_order={fd_order!r}; allowed: (2, 4, 6). "
            "fd_order 8 reaches five ghost points, which the pinned "
            "Dendrolib proves at element order 10; it is outside this "
            "builder's qualified set rather than host-gated."
        )
    # Single-registry authority: the builder
    # profile is written into the registered Dendro generation parameters
    # and validated immediately, so generation parameters, equation hash,
    # and kernel can never skew.
    generation_parameters.validate_generation_parameters()
    # The qualified (non-nested) threading profile is serial.  The
    # point kernel runs inside Dendro's own block traversal, so an inner OpenMP
    # pragma would nest parallelism.  Assert rather than overwrite: silently
    # discarding a caller's request would produce an unqualified configuration
    # whose manifest disagrees with the invocation.
    require_serial_parallelization()
    bundle = build_fccz4_expression_bundle(
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    )
    rhs_symbols = tuple(bundle.rhs_by_symbol_name)
    lvalues = tuple(
        f"rhs_{gf_names.rhs_symbol_to_gridfunction_name(symbol)}[pp]"
        for symbol in rhs_symbols
    )
    evol_order = roles.registered_evol_order()
    # The qualified profile is exactly 25 EVOL fields.
    if len(evol_order) != 25:
        raise ValueError(
            "fCCZ4 EVOL registry must hold exactly 25 fields, found "
            f"{len(evol_order)}."
        )
    # The RHS symbols must map bijectively onto the registered EVOL fields
    # (25 EVOL fields, under their exact registered names).
    mapped = {
        gf_names.rhs_symbol_to_gridfunction_name(symbol) for symbol in rhs_symbols
    }
    if len(lvalues) != len(set(lvalues)) or mapped != set(evol_order):
        raise ValueError(
            "fCCZ4 RHS symbols do not map bijectively onto the registered "
            f"EVOL fields: missing={sorted(set(evol_order) - mapped)} "
            f"extra={sorted(mapped - set(evol_order))}"
        )
    # The upwind control fields are the EVOL gridfunctions that appear in the
    # shared factory's upwind control vector (e.g. vetU0/1/2 for the
    # canonical fCCZ4 profile).  Derived, not hardcoded, so a Gate 4 harness
    # can drive positive/negative/zero control on exactly these fields.
    upwind_control_fields = bkh.upwind_control_fields_from_control_vec(
        bundle.upwind_control_vec, evol_order
    )
    par.set_parval_from_str("fd_order", fd_order)
    # Both scalar spellings come from the registries,
    # and every codegen option is passed explicitly.
    fp_type = str(par.parval_from_str("fp_type"))
    scalar_type = gri.DENDRO_SCALAR_TYPE
    kernel = c_codegen(
        list(bundle.rhs_by_symbol_name.values()),
        list(lvalues),
        enable_fd_codegen=True,
        enable_fd_functions=False,
        enable_simd=False,
        fp_type=fp_type,
        fp_type_alias=scalar_type,
        mem_alloc_style="210",
        rational_const_alias="static const",
        verbose=False,
        upwind_control_vec=list(bundle.upwind_control_vec),
    )
    # Padding is the widest reach of the derivative operators the emitted
    # kernel actually contains, taken per axis from the same coefficient
    # source the kernel was lowered with.  It is NOT fd_order // 2: the
    # upwinded and Kreiss-Oliger families reach one point further than the
    # centred ones (at fd_order 4, dupD reaches 3 while dD reaches 2), so a
    # radius-derived padding reads past the end of a Dendro block.
    # The consumed CodeParameters are the expression free symbols that are
    # registered CodeParameters.  Reading the symbols the equations actually
    # contain is exact; scanning the emitted C text for names is not.  Sorted
    # for a deterministic CFunction parameter order (the caller forwards the
    # values in this order).
    used_codeparameters = bkh.used_codeparameters(bundle.rhs_by_symbol_name.values())
    point_loop_body = bkh.point_loop(kernel)
    block_body = (
        bkh.block_pointer_bindings(evol_order, scalar_type) + "\n" + point_loop_body
    )
    all_blocks_body = block_loop(
        f"{FCCZ4_RHS_EVAL_BLOCK_CFUNCTION}(mesh.geom[blk], in_gfs, rhs_gfs"
        + (
            f", {bkh.cparam_arguments(used_codeparameters)}"
            if used_codeparameters
            else ""
        )
        + ");",
        num_blocks="mesh.num_blocks",
    )
    # Single-body rule: the flat adapter derives flat-layout
    # pointers, packs the per-component pointer arrays, and calls the
    # registered block kernel.  There is exactly one numerical body, so the
    # two paths cannot diverge.
    flat_call_args = (
        f", {bkh.cparam_arguments(used_codeparameters)}" if used_codeparameters else ""
    )
    flat_block_body = (
        bkh.flat_block_pointer_bindings(evol_order, scalar_type)
        + f"\nconst {scalar_type}* const in_gfs_call[] = {{"
        + ", ".join(gf_names.input_pointer(name) for name in evol_order)
        + f"}};\n{scalar_type}* rhs_gfs_call[] = {{"
        + ", ".join(gf_names.rhs_pointer(name) for name in evol_order)
        + "};\n"
        + f"{FCCZ4_RHS_EVAL_BLOCK_CFUNCTION}(geom, in_gfs_call, rhs_gfs_call{flat_call_args});\n"
    )
    cparam_args = bkh.cparam_declarations(used_codeparameters)
    block_params = (
        f"const BlockGeometry& geom, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* rhs_gfs" + (f", {cparam_args}" if cparam_args else "")
    )
    # ``StandaloneHostMesh`` and ``BlockGeometry`` are the standalone host's
    # stand-ins for the real host's types.  The production
    # ``ot::Block``/``RuntimeGeometry`` signatures are settled only when the
    # generated solver is built against a real Dendro-GR checkout, so these
    # prototypes must not be treated as the module ABI.
    all_blocks_params = (
        f"const StandaloneHostMesh& mesh, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* rhs_gfs" + (f", {cparam_args}" if cparam_args else "")
    )
    flat_block_params = (
        f"const BlockGeometry& geom, const {scalar_type}* const in_gfs_flat, "
        f"{scalar_type}* const rhs_gfs_flat"
        + (f", {cparam_args}" if cparam_args else "")
    )
    kernel_expressions = list(bundle.rhs_by_symbol_name.values())
    operators = bkh.emitted_derivative_operators(
        kernel_expressions, list(bundle.upwind_control_vec)
    )
    padding = bkh.padding_from_derivative_operators(
        kernel_expressions, list(bundle.upwind_control_vec), fd_order
    )
    # Single KO ownership: dKOD operators are present in the emitted kernel if
    # and only if Kreiss-Oliger dissipation was requested.
    if (
        any(operator.startswith("dKOD") for operator in operators)
        != enable_KreissOliger_dissipation
    ):
        raise ValueError(
            "dKOD-operator presence does not match "
            f"enable_KreissOliger_dissipation={enable_KreissOliger_dissipation!r}."
        )
    # Recorded beside the padding so the state header renders the positions
    # the emitted kernel actually upwinds on, rather than an empty table.
    roles.set_upwind_control_fields(upwind_control_fields)
    roles.set_required_padding(padding)
    return FCCZ4RHSBuild(
        evol_order=evol_order,
        upwind_control_fields=upwind_control_fields,
        lvalues=lvalues,
        padding=padding,
        block_body=block_body,
        block_params=block_params,
        all_blocks_body=all_blocks_body,
        all_blocks_params=all_blocks_params,
        flat_block_body=flat_block_body,
        flat_block_params=flat_block_params,
        used_codeparameters=tuple(used_codeparameters),
        rhs_by_symbol_name=dict(bundle.rhs_by_symbol_name),
    )


def _register_CFunctions_rhs_eval_fCCZ4(
    *,
    fd_order: int,
    enable_KreissOliger_dissipation: bool,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    solver_stem: str,
) -> None:
    """
    Register the per-block, all-block, and flat-block fCCZ4 RHS CFunctions.

    The ghost points the emitted operators reach are recorded through
    :func:`nrpy.infrastructures.Dendro.CFunction_roles.set_required_padding`, so
    the state header and the parameter file read them from the registry.

    :param fd_order: The finite-difference order (2, 4, or 6).  Order 8
        reaches five ghost points; padding 5 is proven on the pinned
        Dendrolib, so order 8 is a generator limit rather than a host limit.
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :param solver_stem: Lowercase formulation stem for the emitted include.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    """
    build = _build_rhs_eval_fCCZ4(
        fd_order=fd_order,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    subdirectory = "generated/src/rhs_eval"
    cfc.register_CFunction(
        name=FCCZ4_RHS_EVAL_BLOCK_CFUNCTION,
        desc="Per-block direct-FD fCCZ4 RHS (25 fields).",
        subdirectory=subdirectory,
        params=build.block_params,
        body=build.block_body,
        includes=[f"{solver_stem}_defines.h"],
    )
    roles.set_CFunction_role(FCCZ4_RHS_EVAL_BLOCK_CFUNCTION, "rhs_eval_block")
    roles.set_CFunction_codeparameters(
        FCCZ4_RHS_EVAL_BLOCK_CFUNCTION, build.used_codeparameters
    )
    cfc.register_CFunction(
        name=FCCZ4_RHS_EVAL_ALL_BLOCKS_CFUNCTION,
        desc="All-block direct-FD fCCZ4 RHS (NRPy block loop).",
        subdirectory=subdirectory,
        params=build.all_blocks_params,
        body=build.all_blocks_body,
        includes=[f"{solver_stem}_defines.h"],
    )
    roles.set_CFunction_role(FCCZ4_RHS_EVAL_ALL_BLOCKS_CFUNCTION, "rhs_eval")
    roles.set_CFunction_codeparameters(
        FCCZ4_RHS_EVAL_ALL_BLOCKS_CFUNCTION, build.used_codeparameters
    )
    cfc.register_CFunction(
        name=FCCZ4_RHS_EVAL_FLAT_BLOCK_CFUNCTION,
        desc="LTS flat-block adapter (same numerical body, flat layout).",
        subdirectory=subdirectory,
        params=build.flat_block_params,
        body=build.flat_block_body,
        includes=[f"{solver_stem}_defines.h"],
    )
    roles.set_CFunction_role(FCCZ4_RHS_EVAL_FLAT_BLOCK_CFUNCTION, "rhs_eval_flat_block")
    roles.set_CFunction_codeparameters(
        FCCZ4_RHS_EVAL_FLAT_BLOCK_CFUNCTION, build.used_codeparameters
    )


BSSN_RHS_EVAL_BLOCK_CFUNCTION = "bssn_rhs_eval_block"
BSSN_RHS_EVAL_ALL_BLOCKS_CFUNCTION = "bssn_rhs_eval"
BSSN_RHS_EVAL_FLAT_BLOCK_CFUNCTION = "bssn_rhs_eval_flat_block"

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
    :param padding: Ghost points the emitted operators reach, widest axis.
    :param block_body: The per-block CFunction body (point loop + bindings).
    :param block_params: The per-block CFunction parameter list.
    :param all_blocks_body: The all-block CFunction body (NRPy block loop).
    :param all_blocks_params: The all-block CFunction parameter list.
    :param flat_block_body: The flat-block adapter CFunction body.
    :param flat_block_params: The flat-block adapter CFunction parameter list.
    :param used_codeparameters: The CodeParameters the signatures forward, in
        signature order, computed from the expression free symbols.
    :param rhs_by_symbol_name: The assembled symbolic right-hand sides, kept so
        the module's ``__main__`` can pin them against trusted values without
        reassembling them.
    """

    evol_order: Tuple[str, ...]
    upwind_control_fields: Tuple[str, ...]
    lvalues: Tuple[str, ...]
    padding: int
    block_body: str
    block_params: str
    all_blocks_body: str
    all_blocks_params: str
    flat_block_body: str
    flat_block_params: str
    used_codeparameters: Tuple[str, ...]
    rhs_by_symbol_name: Mapping[str, sp.Expr]


def BSSN_rhs_expressions(
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


def _build_rhs_eval_BSSN(
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
        reaches five ghost points, which the pinned Dendrolib proves at element
        order 10; the limit recorded in
        the pinned Dendrolib's proven axes, so it stays capability-gated.
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
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fp_type", "double")
    >>> par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    >>> import nrpy.grid as _gri
    >>> from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities as _bq
    >>> from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs as _brhs
    >>> _gri.glb_gridfcs_dict.clear()
    >>> par.glb_extras_dict.pop("Dendro", None) and None
    >>> # The equations layer memoizes per CoordSystem and registers the evolved
    >>> # state in its constructor, so the memo must go with the registry or the
    >>> # rebuild finds nothing registered.
    >>> _bq.clear() or _brhs.clear()
    >>> try:
    ...     _build_rhs_eval_BSSN(fd_order=8, enable_KreissOliger_dissipation=False)
    ... except ValueError as error:
    ...     print(str(error).splitlines()[0])
    Unsupported fd_order=8; allowed: (2, 4, 6). fd_order 8 reaches five ghost points, which the pinned Dendrolib proves at element order 10; it is outside this builder's qualified set rather than host-gated.
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _build = _build_rhs_eval_BSSN(fd_order=4, enable_KreissOliger_dissipation=False)
    >>> len(_build.evol_order), "Theta_fCCZ4" in _build.evol_order
    (24, False)
    >>> sorted(_build.lvalues)[:2]
    ['rhs_aDD00[pp]', 'rhs_aDD01[pp]']
    >>> _build.padding
    3
    >>> _build.upwind_control_fields
    ('vetU0', 'vetU1', 'vetU2')
    >>> roles.upwind_control_fields() == _build.upwind_control_fields
    True
    >>> roles.required_padding()
    3
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro BSSN RHS, got "
            f"{par.parval_from_str('Infrastructure')!r}."
        )
    if fd_order not in (2, 4, 6):
        raise ValueError(
            f"Unsupported fd_order={fd_order!r}; allowed: (2, 4, 6). "
            "fd_order 8 reaches five ghost points, which the pinned "
            "Dendrolib proves at element order 10; it is outside this "
            "builder's qualified set rather than host-gated."
        )
    generation_parameters.validate_generation_parameters()
    require_serial_parallelization()

    rhs_by_symbol_name, upwind_control_vec = BSSN_rhs_expressions(
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    )
    rhs_symbols = tuple(rhs_by_symbol_name)
    lvalues = tuple(
        f"rhs_{gf_names.rhs_symbol_to_gridfunction_name(symbol)}[pp]"
        for symbol in rhs_symbols
    )
    evol_order = roles.registered_evol_order()
    if len(evol_order) != BSSN_EVOL_COUNT:
        raise ValueError(
            f"BSSN EVOL registry must hold exactly {BSSN_EVOL_COUNT} fields, "
            f"found {len(evol_order)}."
        )
    mapped = {
        gf_names.rhs_symbol_to_gridfunction_name(symbol) for symbol in rhs_symbols
    }
    if len(lvalues) != len(set(lvalues)) or mapped != set(evol_order):
        raise ValueError(
            "BSSN RHS symbols do not map bijectively onto the registered EVOL "
            f"fields: missing={sorted(set(evol_order) - mapped)} "
            f"extra={sorted(mapped - set(evol_order))}"
        )
    upwind_control_fields = bkh.upwind_control_fields_from_control_vec(
        upwind_control_vec, evol_order
    )

    par.set_parval_from_str("fd_order", fd_order)
    fp_type = str(par.parval_from_str("fp_type"))
    scalar_type = gri.DENDRO_SCALAR_TYPE
    kernel = c_codegen(
        list(rhs_by_symbol_name.values()),
        list(lvalues),
        enable_fd_codegen=True,
        enable_fd_functions=False,
        enable_simd=False,
        fp_type=fp_type,
        fp_type_alias=scalar_type,
        mem_alloc_style="210",
        rational_const_alias="static const",
        verbose=False,
        upwind_control_vec=list(upwind_control_vec),
    )
    kernel_expressions = list(rhs_by_symbol_name.values())
    operators = bkh.emitted_derivative_operators(
        kernel_expressions, list(upwind_control_vec)
    )
    padding = bkh.padding_from_derivative_operators(
        kernel_expressions, list(upwind_control_vec), fd_order
    )
    if (
        any(operator.startswith("dKOD") for operator in operators)
        != enable_KreissOliger_dissipation
    ):
        raise ValueError(
            "dKOD-operator presence does not match "
            f"enable_KreissOliger_dissipation={enable_KreissOliger_dissipation!r}."
        )

    used_codeparameters = bkh.used_codeparameters(rhs_by_symbol_name.values())
    cparam_args = bkh.cparam_declarations(used_codeparameters)
    cparam_values = bkh.cparam_arguments(used_codeparameters)
    tail_args = f", {cparam_args}" if cparam_args else ""
    tail_values = f", {cparam_values}" if cparam_values else ""

    block_body = (
        bkh.block_pointer_bindings(evol_order, scalar_type)
        + "\n"
        + bkh.point_loop(kernel)
    )
    block_params = (
        f"const BlockGeometry& geom, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* rhs_gfs" + tail_args
    )
    all_blocks_body = block_loop(
        f"{BSSN_RHS_EVAL_BLOCK_CFUNCTION}(mesh.geom[blk], in_gfs, rhs_gfs{tail_values});",
        num_blocks="mesh.num_blocks",
    )
    all_blocks_params = (
        f"const StandaloneHostMesh& mesh, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* rhs_gfs" + tail_args
    )
    # Single-body rule: the flat adapter derives flat-layout pointers, packs
    # the per-component pointer arrays, and calls the registered block kernel.
    # There is exactly one numerical body, so the two paths cannot diverge.
    flat_block_body = (
        bkh.flat_block_pointer_bindings(evol_order, scalar_type)
        + f"\nconst {scalar_type}* const in_gfs_call[] = {{"
        + ", ".join(gf_names.input_pointer(name) for name in evol_order)
        + f"}};\n{scalar_type}* rhs_gfs_call[] = {{"
        + ", ".join(gf_names.rhs_pointer(name) for name in evol_order)
        + "};\n"
        + f"{BSSN_RHS_EVAL_BLOCK_CFUNCTION}(geom, in_gfs_call, rhs_gfs_call{tail_values});\n"
    )
    flat_block_params = (
        f"const BlockGeometry& geom, const {scalar_type}* const in_gfs_flat, "
        f"{scalar_type}* const rhs_gfs_flat" + tail_args
    )

    roles.set_upwind_control_fields(upwind_control_fields)
    roles.set_required_padding(padding)
    return BSSNRHSBuild(
        evol_order=evol_order,
        upwind_control_fields=upwind_control_fields,
        lvalues=lvalues,
        padding=padding,
        block_body=block_body,
        block_params=block_params,
        all_blocks_body=all_blocks_body,
        all_blocks_params=all_blocks_params,
        flat_block_body=flat_block_body,
        flat_block_params=flat_block_params,
        used_codeparameters=tuple(used_codeparameters),
        rhs_by_symbol_name=dict(rhs_by_symbol_name),
    )


def _register_CFunctions_rhs_eval_BSSN(
    *,
    fd_order: int,
    enable_KreissOliger_dissipation: bool,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    solver_stem: str,
) -> None:
    """
    Register the per-block, all-block, and flat-block BSSN RHS CFunctions.

    :param fd_order: The finite-difference order (2, 4, or 6).
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :param solver_stem: Lowercase formulation stem for the emitted include.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    """
    build = _build_rhs_eval_BSSN(
        fd_order=fd_order,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    subdirectory = "generated/src/rhs_eval"
    cfc.register_CFunction(
        name=BSSN_RHS_EVAL_BLOCK_CFUNCTION,
        desc="Per-block direct-FD BSSN RHS (24 fields).",
        subdirectory=subdirectory,
        params=build.block_params,
        body=build.block_body,
        includes=[f"{solver_stem}_defines.h"],
    )
    roles.set_CFunction_role(BSSN_RHS_EVAL_BLOCK_CFUNCTION, "rhs_eval_block")
    roles.set_CFunction_codeparameters(
        BSSN_RHS_EVAL_BLOCK_CFUNCTION, build.used_codeparameters
    )
    cfc.register_CFunction(
        name=BSSN_RHS_EVAL_ALL_BLOCKS_CFUNCTION,
        desc="All-block direct-FD BSSN RHS (NRPy block loop).",
        subdirectory=subdirectory,
        params=build.all_blocks_params,
        body=build.all_blocks_body,
        includes=[f"{solver_stem}_defines.h"],
    )
    roles.set_CFunction_role(BSSN_RHS_EVAL_ALL_BLOCKS_CFUNCTION, "rhs_eval")
    roles.set_CFunction_codeparameters(
        BSSN_RHS_EVAL_ALL_BLOCKS_CFUNCTION, build.used_codeparameters
    )
    cfc.register_CFunction(
        name=BSSN_RHS_EVAL_FLAT_BLOCK_CFUNCTION,
        desc="LTS flat-block adapter (same numerical body, flat layout).",
        subdirectory=subdirectory,
        params=build.flat_block_params,
        body=build.flat_block_body,
        includes=[f"{solver_stem}_defines.h"],
    )
    roles.set_CFunction_role(BSSN_RHS_EVAL_FLAT_BLOCK_CFUNCTION, "rhs_eval_flat_block")
    roles.set_CFunction_codeparameters(
        BSSN_RHS_EVAL_FLAT_BLOCK_CFUNCTION, build.used_codeparameters
    )


def build_rhs_eval(
    *,
    enable_fCCZ4: bool = False,
    fd_order: int = 4,
    enable_KreissOliger_dissipation: bool = True,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
) -> Union[BSSNRHSBuild, FCCZ4RHSBuild]:
    """
    Build the Dendro right-hand-side kernels for one formulation.

    One entry point for both formulations, taking the formulation as an
    argument, exactly as ``BHaH/general_relativity/rhs_eval.py`` does.  The two
    branches stay separate below: they assemble their expression sets by
    different routes, and BHaH's own divergence guard records that unifying the
    branches is not wanted.

    :param enable_fCCZ4: Build fCCZ4 instead of BSSN.
    :param fd_order: The finite-difference order (2, 4, or 6).
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :return: The formulation's build record.
    """
    if enable_fCCZ4:
        return _build_rhs_eval_fCCZ4(
            fd_order=fd_order,
            enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
            CoordSystem=CoordSystem,
            LapseEvolutionOption=LapseEvolutionOption,
            ShiftEvolutionOption=ShiftEvolutionOption,
        )
    return _build_rhs_eval_BSSN(
        fd_order=fd_order,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )


def register_CFunctions_rhs_eval(
    *,
    enable_fCCZ4: bool = False,
    fd_order: int = 4,
    enable_KreissOliger_dissipation: bool = True,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    solver_stem: str,
) -> None:
    """
    Register the right-hand-side CFunctions for one formulation.

    :param enable_fCCZ4: Register fCCZ4 instead of BSSN.
    :param fd_order: The finite-difference order (2, 4, or 6).
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :param solver_stem: Lowercase formulation stem for the emitted include.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    """
    if enable_fCCZ4:
        _register_CFunctions_rhs_eval_fCCZ4(
            solver_stem=solver_stem,
            fd_order=fd_order,
            enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
            CoordSystem=CoordSystem,
            LapseEvolutionOption=LapseEvolutionOption,
            ShiftEvolutionOption=ShiftEvolutionOption,
        )
        return
    _register_CFunctions_rhs_eval_BSSN(
        solver_stem=solver_stem,
        fd_order=fd_order,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Symbolic pinning of the two shipped right-hand sides.  coding_style.md
    # forbids a golden C file for a kernel this size and names this route
    # instead, so the trusted values are the assembled expressions, evaluated
    # at fixed mpf substitutions exactly as
    # nrpy/infrastructures/CarpetX/general_relativity/rhs_eval.py pins its own.
    # The profiles are the ones the two examples ship: Cartesian, chi for
    # fCCZ4 and W for BSSN, Kreiss-Oliger off.
    import os

    import nrpy.validate_expressions.validate_expressions as ve
    from nrpy.infrastructures.Dendro.general_relativity import trusted_capture

    par.set_parval_from_str("Infrastructure", "Dendro")
    par.set_parval_from_str("parallelization", "none")
    par.set_parval_from_str("fp_type", "double")
    par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    par.set_parval_from_str("fd_order", 4)
    for sweep_fCCZ4, sweep_cf in trusted_capture.SHIPPED_PROFILES:
        trusted_capture.reset_generation_state()
        par.set_parval_from_str("EvolvedConformalFactor_cf", sweep_cf)
        sweep_build = build_rhs_eval(
            enable_fCCZ4=sweep_fCCZ4,
            fd_order=4,
            enable_KreissOliger_dissipation=False,
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}"
            f"_{trusted_capture.SHIPPED_GAUGE}"
            f"_Cartesian_{sweep_cf}_fCCZ4{sweep_fCCZ4}_KOFalse",
            ve.process_dictionary_of_expressions(
                dict(sweep_build.rhs_by_symbol_name), fixed_mpfs_for_free_symbols=True
            ),
        )
