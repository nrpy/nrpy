"""
Direct finite-difference fCCZ4 right-hand side for the Dendro infrastructure.

This module lowers the shared fCCZ4 expression factory
(:func:`nrpy.equations.general_relativity.fCCZ4_system.build_fccz4_expression_bundle`)
into a Dendro per-block direct-FD RHS kernel:

* map each RHS symbol to its registered EVOL gridfunction name algorithmically
  (:func:`nrpy.infrastructures.Dendro.gridfunction_name_decorations.rhs_symbol_to_gridfunction_name`)
  and assert the bijection against the registry;
* derive the output lvalues and the ``in_<name>`` / ``rhs_<name>`` pointer
  bindings from the NRPy gridfunction registry (no field name is hardcoded);
* run :func:`nrpy.c_codegen.c_codegen` with direct finite differences, using
  the shared factory's upwind control vector and the ``DendroScalar`` alias;
* emit the point and block loops through the NRPy Dendro loop helpers;
* provide the three registered CFunction bodies (per-block, all-block, and the
  local-time-stepping flat-block adapter, which reuses the same numerical
  body); and
* expose the emitted operator records and the ghost-point padding those
  operators reach.

:func:`build_rhs_eval` registers no CFunction; it records the padding and the
upwind-control set it derived into the Dendro area of the NRPy parameter
registry, and :func:`register_CFunctions_rhs_eval` registers the bodies it
returns.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from dataclasses import dataclass
from typing import Tuple

import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.fCCZ4_system import (
    build_fccz4_expression_bundle,
)
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import block_kernel_helpers as bkh
from nrpy.infrastructures.Dendro import generation_parameters  # noqa: F401
from nrpy.infrastructures.Dendro import gridfunction_name_decorations as gf_names
from nrpy.infrastructures.Dendro.block_loop import block_loop
from nrpy.infrastructures.Dendro.simple_loop import require_serial_parallelization

# The per-block CFunction name (the Dendro scheduling role key).
RHS_EVAL_BLOCK_CFUNCTION = "fccz4_rhs_eval_block"

# All-block CFunction name.
RHS_EVAL_GLOBAL_CFUNCTION = "fccz4_rhs_eval"

# LTS flat-block adapter CFunction name (same numerical body, flat layout).
RHS_EVAL_FLAT_BLOCK_CFUNCTION = "fccz4_rhs_eval_flat_block"


@dataclass(frozen=True)
class FCCZ4RHSBuild:
    """
    Immutable result of building the direct-FD fCCZ4 RHS for one profile.

    :param evol_order: The 25 EVOL gridfunction names, in registry order.
    :param upwind_control_fields: The EVOL names the emitted kernel upwinds on.
    :param lvalues: The 25 output lvalues (``rhs_<name>[pp]``).
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


def build_rhs_eval(
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
    the returned bodies through :func:`register_Dendro_CFunction`.

    The caller must set ``Infrastructure`` to ``Dendro`` and
    ``parallelization`` to ``"none"`` before calling; both are asserted here
    rather than overwritten, because silently discarding a caller's request
    would produce an unqualified configuration.  The registered Dendro
    generation parameters are validated before anything is lowered.

    :param fd_order: The finite-difference order (2, 4, or 6).  Order 8
        reaches five ghost points and stays capability-gated until Dendrolib
        is pinned.
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
    >>> try:
    ...     build_rhs_eval(fd_order=8, enable_KreissOliger_dissipation=False)
    ... except ValueError as error:
    ...     print(str(error).splitlines()[0])
    Unsupported fd_order=8; allowed: (2, 4, 6). fd_order 8 reaches five ghost points, above the max_proven_padding of 4 recorded in dendrolib_capabilities.json, so it stays capability-gated until Dendrolib is pinned.
    >>> import contextlib, io
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _build = build_rhs_eval(fd_order=4, enable_KreissOliger_dissipation=False)
    >>> len(_build.evol_order), "Theta_fCCZ4" in _build.evol_order
    (25, True)
    >>> sorted(_build.lvalues)[:2]
    ['rhs_Theta_fCCZ4[pp]', 'rhs_aDD00[pp]']
    >>> _build.padding
    (3, 3, 3)
    >>> _build.upwind_control_fields
    ('vetU0', 'vetU1', 'vetU2')
    >>> roles.upwind_control_fields()
    ('vetU0', 'vetU1', 'vetU2')
    >>> roles.required_padding()
    (3, 3, 3)
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro fCCZ4 RHS, got "
            f"{par.parval_from_str('Infrastructure')!r}."
        )
    if fd_order not in (2, 4, 6):
        raise ValueError(
            f"Unsupported fd_order={fd_order!r}; allowed: (2, 4, 6). "
            "fd_order 8 reaches five ghost points, above the "
            "max_proven_padding of 4 recorded in dendrolib_capabilities.json, "
            "so it stays capability-gated until Dendrolib is pinned."
        )
    # Single-registry authority: the builder
    # profile is written into the registered Dendro generation parameters
    # and validated immediately, so generation parameters, equation hash,
    # and kernel can never skew.
    par.set_parval_from_str(
        "Dendro_enable_KreissOliger_dissipation", enable_KreissOliger_dissipation
    )
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
    upwind_control_fields = bkh.upwind_control_fields(
        bundle.upwind_control_vec, evol_order
    )
    par.set_parval_from_str("fd_order", fd_order)
    # Both scalar spellings come from the registries,
    # and every codegen option is passed explicitly.
    fp_type = str(par.parval_from_str("fp_type"))
    fp_type_alias = str(par.parval_from_str("Dendro_scalar_type"))
    kernel = c_codegen(
        list(bundle.rhs_by_symbol_name.values()),
        list(lvalues),
        enable_fd_codegen=True,
        enable_fd_functions=False,
        enable_simd=False,
        fp_type=fp_type,
        fp_type_alias=fp_type_alias,
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
        bkh.block_pointer_bindings(evol_order, fp_type_alias) + "\n" + point_loop_body
    )
    global_body = block_loop(
        f"{RHS_EVAL_BLOCK_CFUNCTION}(world.geom[blk], in_gfs, rhs_gfs"
        + (
            f", {bkh.cparam_arguments(used_codeparameters)}"
            if used_codeparameters
            else ""
        )
        + ");",
        num_blocks="world.num_blocks",
    )
    # Single-body rule: the flat adapter derives flat-layout
    # pointers, packs the per-component pointer arrays, and calls the
    # registered block kernel.  There is exactly one numerical body, so the
    # two paths cannot diverge.
    flat_call_args = (
        f", {bkh.cparam_arguments(used_codeparameters)}" if used_codeparameters else ""
    )
    flat_block_body = (
        bkh.flat_block_pointer_bindings(evol_order, fp_type_alias)
        + f"\nconst {fp_type_alias}* const in_gfs_call[] = {{"
        + ", ".join(gf_names.input_pointer(name) for name in evol_order)
        + f"}};\n{fp_type_alias}* rhs_gfs_call[] = {{"
        + ", ".join(gf_names.rhs_pointer(name) for name in evol_order)
        + "};\n"
        + f"{RHS_EVAL_BLOCK_CFUNCTION}(geom, in_gfs_call, rhs_gfs_call{flat_call_args});\n"
    )
    cparam_args = bkh.cparam_declarations(used_codeparameters, fp_type_alias)
    block_params = (
        f"const BlockGeometry& geom, const {fp_type_alias}* const* in_gfs, "
        f"{fp_type_alias}* const* rhs_gfs" + (f", {cparam_args}" if cparam_args else "")
    )
    # Mock-vehicle signatures: ``MockWorld``/``BlockGeometry``
    # are test doubles for the mock-host compile.  The production
    # ``ot::Block``/``RuntimeGeometry`` signatures are settled only once
    # Dendrolib is pinned, so these mock prototypes must not be treated as
    # the module ABI.
    global_params = (
        f"const MockWorld& world, const {fp_type_alias}* const* in_gfs, "
        f"{fp_type_alias}* const* rhs_gfs" + (f", {cparam_args}" if cparam_args else "")
    )
    flat_block_params = (
        f"const BlockGeometry& geom, const {fp_type_alias}* const in_gfs_flat, "
        f"{fp_type_alias}* const rhs_gfs_flat"
        + (f", {cparam_args}" if cparam_args else "")
    )
    operator_records = bkh.emitted_operators(kernel, fd_order)
    padding = bkh.padding_from_operators(operator_records)
    # Single KO ownership: dKOD operators are present in the emitted kernel if
    # and only if Kreiss-Oliger dissipation was requested.
    if (
        any(record["operator"].startswith("dKOD") for record in operator_records)
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
    Register the per-block, all-block, and flat-block fCCZ4 RHS CFunctions.

    The ghost points the emitted operators reach are recorded through
    :func:`nrpy.infrastructures.Dendro.CFunction_roles.set_required_padding`, so
    the state header and the parameter file read them from the registry.

    :param fd_order: The finite-difference order (2, 4, or 6).  Order 8
        reaches five ghost points, above the ``max_proven_padding`` of 4 in
        ``dendrolib_capabilities.json``, so it stays capability-gated.
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    """
    build = build_rhs_eval(
        fd_order=fd_order,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    subdirectory = "generated/src/rhs_eval"
    roles.register_Dendro_CFunction(
        role="rhs_eval_block",
        name=RHS_EVAL_BLOCK_CFUNCTION,
        desc="Per-block direct-FD fCCZ4 RHS (25 fields).",
        subdirectory=subdirectory,
        params=build.block_params,
        body=build.block_body,
    )
    roles.register_Dendro_CFunction(
        role="rhs_eval",
        name=RHS_EVAL_GLOBAL_CFUNCTION,
        desc="All-block direct-FD fCCZ4 RHS (NRPy block loop).",
        subdirectory=subdirectory,
        params=build.global_params,
        body=build.global_body,
    )
    roles.register_Dendro_CFunction(
        role="rhs_eval_flat_block",
        name=RHS_EVAL_FLAT_BLOCK_CFUNCTION,
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
