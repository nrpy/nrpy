"""
NRPy PR 5: direct finite-difference fCCZ4 RHS for the Dendro backend.

This module lowers the shared fCCZ4 expression factory (:func:
``nrpy.equations.general_relativity.fCCZ4_system.build_fccz4_expression_bundle``)
to a Dendro per-block direct-FD RHS kernel, exactly as the whitepaper mandates:

* map each RHS symbol to its registered EVOL gridfunction name algorithmically
  (:func:`nrpy.infrastructures.Dendro.naming.rhs_symbol_to_gridfunction_name`)
  and assert the bijection against the registry;
* derive the output lvalues and the ``in_<name>`` / ``rhs_<name>`` pointer
  bindings from the NRPy gridfunction registry (no field name is hardcoded);
* run :func:`nrpy.c_codegen.c_codegen` with direct finite differences inside a
  scoped Dendro access-capture context, using the shared factory's upwind
  control vector and the ``DendroScalar`` alias;
* emit the point and block loops through the NRPy Dendro loop helpers;
* provide the three registered CFunction bodies (per-block, all-block, and the
  LTS flat-block adapter, which reuses the same numerical body); and
* expose the operator and access manifests and the access-derived padding.

The module is a *builder*: it does not register anything by itself.  The caller
(usually the PR 5 exit test) registers the returned CFunctions through
:func:`nrpy.infrastructures.Dendro.registration.register_dendro_CFunction`,
freezes the environment, and runs the Gate 3 mock-compile and Gate 4
fixed-block BHaH equivalence checks.
"""

import re
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import sympy as sp

import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.fCCZ4_system import (
    build_fccz4_expression_bundle,
)
from nrpy.finite_difference import compute_fdcoeffs_fdstencl
from nrpy.infrastructures.Dendro import access_capture as cap
from nrpy.infrastructures.Dendro import generation_parameters  # noqa: F401
from nrpy.infrastructures.Dendro import gridfunction_output as gfo
from nrpy.infrastructures.Dendro import naming
from nrpy.infrastructures.Dendro import registration as reg
from nrpy.infrastructures.Dendro.block_loop import block_loop
from nrpy.infrastructures.Dendro.freeze import (
    body_uses_symbol,
    canonical_expression_digest,
)
from nrpy.infrastructures.Dendro.simple_loop import (
    interior_loop,
    require_serial_parallelization,
)

# The per-block CFunction name (the Dendro scheduling role key).
BLOCK_CFUNCTION = "fccz4_rhs_block"

# All-block CFunction name.
GLOBAL_CFUNCTION = "fccz4_rhs"

# LTS flat-block adapter CFunction name (same numerical body, flat layout).
FLAT_BLOCK_CFUNCTION = "fccz4_rhs_flat_block"

# Finite-difference operator families emitted by c_codegen.  Field names use
# uppercase ``DD`` / ``U`` components, so these lowercase tokens identify
# derivative operators unambiguously.  The digit run that follows is
# ``<tensor component indices><derivative direction indices>`` -- e.g.
# ``hDD_dDD0112`` is the mixed second derivative in directions (1, 2) of the
# component (0, 1) -- so the run must be matched in full (``\d+``) and the
# *trailing* index/indices taken as the direction.  Matching only one or two
# digits made every derivative of a rank-1 or rank-2 field invisible and
# recorded a rank-1 first derivative under its component index instead of its
# direction (whitepaper sections 8.6/8.7/15.4).
_OPERATOR_RE = re.compile(r"_(dD|dDD|dupD|ddnD|dKOD|dfullupD|dfulldnD)(\d+)\b")

# Families carrying a single direction index; the rest carry two.  The
# classification matches nrpy/finite_difference.py, where `dfullupD` and
# `dfulldnD` are first-derivative (single-direction) operators.
_SINGLE_DIRECTION_FAMILIES = frozenset(
    {"dD", "dupD", "ddnD", "dKOD", "dfullupD", "dfulldnD"}
)


@dataclass(frozen=True)
class FCCZ4RHSBuild:
    """
    Immutable result of building the direct-FD fCCZ4 RHS for one profile.

    :param fd_order: The finite-difference order used.
    :param enable_KO: Whether Kreiss-Oliger dissipation was enabled.
    :param kernel: The c_codegen point-kernel (direct FD, DendroScalar alias).
    :param evol_order: The 25 EVOL gridfunction names, in registry order.
    :param rhs_symbols: The 25 RHS symbol names, in rhs_by_symbol_name order.
    :param lvalues: The 25 output lvalues (``rhs_<name>[pp]``).
    :param padding: The access-capture-derived (px, py, pz).
    :param used_codeparameters: CodeParameter names referenced by the kernel.
    :param block_body: The per-block CFunction body (point loop + bindings).
    :param block_params: The per-block CFunction parameter list.
    :param global_body: The all-block CFunction body (NRPy block loop).
    :param global_params: The all-block CFunction parameter list.
    :param flat_block_body: The flat-block adapter CFunction body.
    :param flat_block_params: The flat-block adapter CFunction parameter list.
    :param operator_manifest: Canonical operator manifest mapping.
    :param access_manifest: The canonical captured access tuples.
    :param canonical_expression_digests: (RHS symbol, canonical digest) pairs
        fingerprinting the symbolic right-hand sides themselves (section
        9.10); freeze hashes these into ``equation_hash``.
    """

    fd_order: int
    enable_KO: bool
    kernel: str
    evol_order: Tuple[str, ...]
    upwind_control_fields: Tuple[str, ...]
    rhs_symbols: Tuple[str, ...]
    lvalues: Tuple[str, ...]
    padding: Tuple[int, int, int]
    used_codeparameters: Tuple[str, ...]
    block_body: str
    block_params: str
    global_body: str
    global_params: str
    flat_block_body: str
    flat_block_params: str
    operator_manifest: Dict[str, Any]
    access_manifest: Tuple[Tuple[str, int, int, int], ...]
    canonical_expression_digests: Tuple[Tuple[str, str], ...]


def _cparam_args(used_codeparameters: Tuple[str, ...], scalar_type: str) -> str:
    """
    Build the trailing ``const <scalar> <cp>, ...`` parameter list.

    :param used_codeparameters: CodeParameter names to add.
    :param scalar_type: The registered Dendro scalar alias (section 6.6).
    :return: Comma-joined parameter declarations (no leading comma).
    """
    return ", ".join(f"const {scalar_type} {name}" for name in used_codeparameters)


def _cparam_values(used_codeparameters: Tuple[str, ...]) -> str:
    """
    Build the argument-forwarding list for the used CodeParameters.

    :param used_codeparameters: CodeParameter names to forward.
    :return: Comma-joined name list (no trailing comma).
    """
    return ", ".join(used_codeparameters)


def _by_position(_name: str, position: int) -> str:
    """
    Return the frozen registry position as the component index expression.

    The mock-vehicle translation units compile the kernels without the
    generated state header, so the bindings use the integer registry position
    rather than ``to_index(EvolVar::…)``; both orderings come from the same
    NRPy list (whitepaper section 5.1).

    :param _name: Exact gridfunction name (unused; the index is positional).
    :param position: Position in the frozen EVOL order.
    :return: The component index expression.
    """
    return str(position)


def _block_pointer_bindings(evol_order: Tuple[str, ...], scalar_type: str) -> str:
    """
    Emit the per-field input and RHS pointer bindings for the block layout.

    The bindings are rendered by the single shared emitter
    (:func:`nrpy.infrastructures.Dendro.gridfunction_output.render_component_bindings`)
    from the registry order, so no field name is hardcoded and the roles and
    per-component base offset cannot drift from the frozen-snapshot renderer.
    Every binding adds ``geom.component_offset`` (whitepaper section 8.2): the
    pointer arrays are allocation-relative, so a nonzero per-component base
    must be applied or multi-block layouts read the wrong component.

    :param evol_order: The 25 EVOL names, in registry order.
    :param scalar_type: The registered Dendro scalar alias.
    :return: The binding statements.
    """
    return (
        gfo.render_component_bindings(
            evol_order,
            scalar_type,
            array="in_gfs",
            role=naming.input_pointer,
            const_pointee=True,
            index_expression=_by_position,
        )
        + "\n"
        + gfo.render_component_bindings(
            evol_order,
            scalar_type,
            array="rhs_gfs",
            role=naming.rhs_pointer,
            const_pointee=False,
            index_expression=_by_position,
        )
    )


def _flat_block_pointer_bindings(evol_order: Tuple[str, ...], scalar_type: str) -> str:
    """
    Emit the per-field bindings for the LTS flat-block layout.

    In this layout field ``f`` occupies
    ``in_gfs_flat + f * (nx * ny * nz)``.  Extents are hoisted into a
    ``ptrdiff_t`` local (Appendix A pattern) and the per-component base
    ``geom.component_offset`` is applied exactly as in the block layout; the
    same shared emitter renders both layouts.

    :param evol_order: The 25 EVOL names, in registry order.
    :param scalar_type: The registered Dendro scalar alias.
    :return: The binding statements.
    """
    return "\n".join(
        [
            "const std::ptrdiff_t vol = static_cast<std::ptrdiff_t>(geom.nx)"
            " * geom.ny * geom.nz;",
            gfo.render_component_bindings(
                evol_order,
                scalar_type,
                array="in_gfs_flat",
                role=naming.input_pointer,
                const_pointee=True,
                index_expression=_by_position,
                flat_stride="vol",
            ),
            gfo.render_component_bindings(
                evol_order,
                scalar_type,
                array="rhs_gfs_flat",
                role=naming.rhs_pointer,
                const_pointee=False,
                index_expression=_by_position,
                flat_stride="vol",
            ),
        ]
    )


def _point_loop(kernel: str) -> str:
    """
    Wrap the point kernel in the NRPy Dendro interior point loop.

    :param kernel: The c_codegen point kernel.
    :return: The interior-loop-wrapped body.
    """
    return interior_loop(
        kernel,
        nx="geom.nx",
        ny="geom.ny",
        nz="geom.nz",
        padding="geom.padding",
        pmin_padded="geom.pmin_padded",
        dx="geom.dx",
    )


def _operator_manifest(
    kernel: str,
    fd_order: int,
    enable_KO: bool,
    *,
    emitted_fields: Tuple[str, ...],
    access_count: int,
) -> Dict[str, Any]:
    """
    Derive the operator manifest from the emitted kernel.

    The set of derivative-operator families actually emitted is scanned from
    the kernel (single source of truth: the generated code), and each
    operator's exact rational stencil is obtained from the existing
    ``compute_fdcoeffs_fdstencl`` (one coefficient source, not a second
    stencil generator).  Each operator record carries the section 15.4
    per-operator fields (name, FD order, exact rational coefficient strings,
    signed offsets, per-axis and total reach); the consuming CFunction, the
    emitted-field set and the captured-access count are kernel-level facts,
    recorded once at the top level rather than repeated in every record,
    because per-operator field attribution is not recoverable after CSE.
    The operator-derived maximum offset per axis is the cross-check target
    for the access-capture-derived padding, which
    :func:`~nrpy.infrastructures.Dendro.freeze.freeze_nrpy_dendro_environment`
    enforces.

    :param kernel: The emitted point kernel.
    :param fd_order: The finite-difference order.
    :param enable_KO: Whether Kreiss-Oliger dissipation was enabled.
    :param emitted_fields: Exact EVOL names the consuming CFunction writes.
    :param access_count: Number of captured accesses of the consuming CFunction.
    :return: The operator manifest mapping.
    :raises ValueError: If a derivative token in the kernel is malformed.
    """
    operators: List[Tuple[str, str, str]] = []
    for match in _OPERATOR_RE.finditer(kernel):
        base = match.group(1)
        idx = match.group(2)
        directions = idx[-1] if base in _SINGLE_DIRECTION_FAMILIES else idx[-2:]
        if len(directions) != (1 if base in _SINGLE_DIRECTION_FAMILIES else 2):
            raise ValueError(
                f"Malformed derivative token {base}{idx} in the emitted kernel."
            )
        derivstring = f"{base}{directions}"
        if derivstring not in {op[0] for op in operators}:
            operators.append((derivstring, base, directions))
    operator_records = []
    for derivstring, _base, _idx in operators:
        # Pass the raw fd_order: compute_fdcoeffs_fdstencl applies the
        # +2 for dKOD internally (exactly as c_codegen lowers the kernel),
        # so this reproduces the operator offsets one coefficient source
        # (not a second stencil generator).
        coeffs, stencils = compute_fdcoeffs_fdstencl(derivstring, fd_order)
        max_offset = 0
        per_axis = [0, 0, 0]
        for stencil in stencils:
            for axis, step in enumerate(stencil):
                per_axis[axis] = max(per_axis[axis], abs(step))
                max_offset = max(max_offset, abs(step))
        operator_records.append(
            {
                "operator": derivstring,
                "fd_order": fd_order,
                "coefficients": [str(coeff) for coeff in coeffs],
                "offsets": [list(stencil) for stencil in stencils],
                "max_offset_per_axis": per_axis,
                "max_offset": max_offset,
            }
        )
    # Field attribution is kernel-level, so it is recorded once here rather
    # than repeated inside every operator record (whitepaper section 15.4).
    return {
        "derivative_backend": "full_stencil",
        "ko_owner": "nrpy_full_stencil",
        "ko_enabled": bool(enable_KO),
        "fd_order": fd_order,
        "consuming_cfunctions": [BLOCK_CFUNCTION],
        "emitted_fields": list(emitted_fields),
        "access_count": access_count,
        "operators": operator_records,
    }


def build_fccz4_rhs(
    *,
    fd_order: int,
    enable_KO: bool,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
) -> FCCZ4RHSBuild:
    """
    Build the direct-FD fCCZ4 RHS for one profile.

    This is a pure builder (no CFunction registration); the caller registers
    the returned bodies through :func:`register_dendro_CFunction`.

    The caller must set ``Infrastructure`` to ``Dendro`` before calling
    (asserted here).  The builder pins the qualified profile into the
    registries itself: ``parallelization`` is set to ``"none"`` (section
    7.2: the point kernel runs inside Dendro's own block traversal, so an
    inner OpenMP pragma would nest parallelism; the ``openmp`` profile is
    deferred until the Dendro thread model is measured) and the Dendro
    generation parameters are synced from the profile arguments, so the
    frozen snapshot describes the built kernel (sections 4.1/6.1/9.3).

    :param fd_order: The finite-difference order (2, 4, 6, or 8; 8 is
        capability-gated at freeze).
    :param enable_KO: Enable Kreiss-Oliger dissipation in the shared factory.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :return: The immutable :class:`FCCZ4RHSBuild` result.
    :raises ValueError: If the fCCZ4 RHS symbols do not map bijectively onto
        the registered EVOL fields, if the EVOL count is not exactly 25, or
        if the kernel's dKOD presence does not match ``enable_KO``.
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro fCCZ4 RHS, got "
            f"{par.parval_from_str('Infrastructure')!r}."
        )
    if fd_order not in (2, 4, 6, 8):
        raise ValueError(f"Unsupported fd_order={fd_order!r}; allowed: (2, 4, 6, 8).")
    # Single-registry authority (sections 4.1/6.1/9.3/11.3): the builder
    # profile is written into the registered Dendro generation parameters
    # and validated immediately, so generation parameters, equation hash,
    # and kernel can never skew.
    par.set_parval_from_str("Dendro_enable_KO", enable_KO)
    par.set_parval_from_str("Dendro_fccz4_CoordSystem", CoordSystem)
    par.set_parval_from_str("Dendro_fccz4_LapseEvolutionOption", LapseEvolutionOption)
    par.set_parval_from_str("Dendro_fccz4_ShiftEvolutionOption", ShiftEvolutionOption)
    generation_parameters.validate_generation_parameters()
    # Section 7.2: the qualified (non-nested) threading profile is serial.  The
    # point kernel runs inside Dendro's own block traversal, so an inner OpenMP
    # pragma would nest parallelism.  Assert rather than overwrite: silently
    # discarding a caller's request would produce an unqualified configuration
    # whose manifest disagrees with the invocation (section 9.3).
    require_serial_parallelization()
    bundle = build_fccz4_expression_bundle(
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
        enable_KreissOliger_dissipation=enable_KO,
    )
    rhs_symbols = tuple(bundle.rhs_by_symbol_name)
    lvalues = tuple(
        f"rhs_{naming.rhs_symbol_to_gridfunction_name(symbol)}[pp]"
        for symbol in rhs_symbols
    )
    evol_order = reg.registered_evol_order()
    # Section 11.4 pins the qualified profile at exactly 25 EVOL fields.
    if len(evol_order) != 25:
        raise ValueError(
            "fCCZ4 EVOL registry must hold exactly 25 fields, found "
            f"{len(evol_order)}."
        )
    # The RHS symbols must map bijectively onto the registered EVOL fields
    # (whitepaper section 5.5: 25 EVOL fields; section 7.1 exact-name rule).
    mapped = {naming.rhs_symbol_to_gridfunction_name(symbol) for symbol in rhs_symbols}
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
    _control_symbols = set()
    for _component in bundle.upwind_control_vec:
        _control_symbols |= set(sp.sympify(str(_component)).free_symbols)
    upwind_control_fields = tuple(
        name for name in evol_order if sp.Symbol(name) in _control_symbols
    )
    par.set_parval_from_str("fd_order", fd_order)
    # Both scalar spellings come from the registries (sections 4.1/6.6),
    # and every section 8.4 codegen option is passed explicitly.
    fp_type = str(par.parval_from_str("fp_type"))
    fp_type_alias = str(par.parval_from_str("Dendro_scalar_type"))
    with cap.capture_gridfunction_accesses(BLOCK_CFUNCTION):
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
    padding = cap.required_padding(BLOCK_CFUNCTION)
    # Sort for a deterministic CFunction parameter order (the caller forwards
    # the values in this order).  Membership is token-aware (the same helper
    # freeze uses for the CodeParameter closure): a substring test would
    # falsely match e.g. ``eta`` inside ``Theta_fCCZ4`` (section 6.2).
    used_codeparameters = tuple(
        sorted(
            name for name in par.glb_code_params_dict if body_uses_symbol(kernel, name)
        )
    )
    point_loop_body = _point_loop(kernel)
    block_body = (
        _block_pointer_bindings(evol_order, fp_type_alias) + "\n" + point_loop_body
    )
    global_body = block_loop(
        f"{BLOCK_CFUNCTION}(world.geom[blk_id], in_gfs, rhs_gfs"
        + (f", {_cparam_values(used_codeparameters)}" if used_codeparameters else "")
        + ");",
        count="world.num_blocks",
    )
    # Section 13.6 single-body rule: the flat adapter derives flat-layout
    # pointers, packs the per-component pointer arrays, and calls the
    # registered block kernel.  There is exactly one numerical body, so the
    # recorded ``calls`` edge is true and the two paths cannot diverge.
    flat_call_args = (
        f", {_cparam_values(used_codeparameters)}" if used_codeparameters else ""
    )
    flat_block_body = (
        _flat_block_pointer_bindings(evol_order, fp_type_alias)
        + f"\nconst {fp_type_alias}* const in_gfs_call[] = {{"
        + ", ".join(naming.input_pointer(name) for name in evol_order)
        + f"}};\n{fp_type_alias}* rhs_gfs_call[] = {{"
        + ", ".join(naming.rhs_pointer(name) for name in evol_order)
        + "};\n"
        + f"{BLOCK_CFUNCTION}(geom, in_gfs_call, rhs_gfs_call{flat_call_args});\n"
    )
    cparam_args = _cparam_args(used_codeparameters, fp_type_alias)
    block_params = (
        f"const BlockGeometry& geom, const {fp_type_alias}* const* in_gfs, "
        f"{fp_type_alias}* const* rhs_gfs" + (f", {cparam_args}" if cparam_args else "")
    )
    # Mock-vehicle signatures (section 13.3): ``MockWorld``/``BlockGeometry``
    # are test doubles for the Gate 3/4 mock compile.  Production
    # ``ot::Block``/``RuntimeGeometry`` signatures freeze only after the I0-1
    # Dendrolib gates; PR 6/7 snapshot consumption must gate on the real
    # types and must not calcify these mock prototypes into the ABI.
    global_params = (
        f"const MockWorld& world, const {fp_type_alias}* const* in_gfs, "
        f"{fp_type_alias}* const* rhs_gfs" + (f", {cparam_args}" if cparam_args else "")
    )
    flat_block_params = (
        f"const BlockGeometry& geom, const {fp_type_alias}* const in_gfs_flat, "
        f"{fp_type_alias}* const rhs_gfs_flat"
        + (f", {cparam_args}" if cparam_args else "")
    )
    # Physics provenance (section 9.10): fingerprint the symbolic expressions
    # themselves, before any CSE or C rendering, so the equation hash cannot
    # be fooled by a formatting change or fooled into agreement by one.
    canonical_expression_digests = tuple(
        (symbol, canonical_expression_digest(bundle.rhs_by_symbol_name[symbol]))
        for symbol in sorted(rhs_symbols, key=lambda s: (s.lower(), s))
    )
    accesses = cap.get_captured_accesses(BLOCK_CFUNCTION)
    operator_manifest = _operator_manifest(
        kernel,
        fd_order,
        enable_KO,
        emitted_fields=evol_order,
        access_count=len(accesses),
    )
    # Section 8.5 NRPy-once KO ownership: dKOD operators are present if and
    # only if KO dissipation was requested.
    if (
        any(
            record["operator"].startswith("dKOD")
            for record in operator_manifest["operators"]
        )
        != enable_KO
    ):
        raise ValueError(
            "dKOD-operator presence does not match enable_KO="
            f"{enable_KO!r} (sections 8.5/18: single KO ownership)."
        )
    return FCCZ4RHSBuild(
        fd_order=fd_order,
        enable_KO=enable_KO,
        kernel=kernel,
        evol_order=evol_order,
        upwind_control_fields=upwind_control_fields,
        rhs_symbols=rhs_symbols,
        lvalues=lvalues,
        padding=padding,
        used_codeparameters=used_codeparameters,
        block_body=block_body,
        block_params=block_params,
        global_body=global_body,
        global_params=global_params,
        flat_block_body=flat_block_body,
        flat_block_params=flat_block_params,
        operator_manifest=operator_manifest,
        access_manifest=accesses,
        canonical_expression_digests=canonical_expression_digests,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
