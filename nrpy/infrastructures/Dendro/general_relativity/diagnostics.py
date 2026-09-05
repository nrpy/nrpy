# nrpy/infrastructures/Dendro/general_relativity/diagnostics.py
"""
PR 9: fCCZ4 constraint diagnostics for the Dendro backend.

The first diagnostic set is the Hamiltonian constraint of the fCCZ4 system and
the spatial Z4 connection constraint (whitepaper section 14.6).  Both come from
the shared expression factory
(:func:`nrpy.equations.general_relativity.fCCZ4_system.build_fccz4_expression_bundle`
with ``enable_diagnostics=True``), so the Dendro kernel and any other
infrastructure lower the same expressions.

The diagnostic gridfunctions are registered in the DIAG group: they are
recomputed from the evolved state and are never authoritative checkpoint state
(section 14.6).  The kernel uses the same finite-difference order and the same
memory-access mechanism as the RHS, so no separate diagnostic profile exists.

Every name in this module is read back from the shared factory or from the
registry; none is written down here.

Author: NRPy Dendro fCCZ4 infrastructure (PR 9)
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import sympy as sp

import nrpy.grid as gri
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.fCCZ4_system import (
    build_fccz4_expression_bundle,
)
from nrpy.infrastructures.Dendro import access_capture as cap
from nrpy.infrastructures.Dendro import generation_parameters
from nrpy.infrastructures.Dendro import gridfunction_output as gfo
from nrpy.infrastructures.Dendro import naming
from nrpy.infrastructures.Dendro import registration as reg
from nrpy.infrastructures.Dendro.block_loop import block_loop
from nrpy.infrastructures.Dendro.simple_loop import (
    interior_loop,
    require_serial_parallelization,
)

# CFunction names (Dendro scheduling role keys, whitepaper section 7.1).
CONSTRAINTS_BLOCK_CFUNCTION = "fccz4_constraints_block"
CONSTRAINTS_GLOBAL_CFUNCTION = "fccz4_constraints"

# NRPy variance letters: an exact tensor-component name ends with a run of
# these followed by one index digit per letter.
_VARIANCE_LETTERS = frozenset({"U", "D"})
_INDEX_DIGITS = "0123456789"


@dataclass(frozen=True)
class FCCZ4DiagnosticsBuild:
    """
    Immutable result of building the diagnostic kernel for one profile.

    :param block_body: The per-block CFunction body (bindings + point loop).
    :param block_params: The per-block CFunction parameter list.
    :param global_body: The all-block CFunction body (NRPy block loop).
    :param global_params: The all-block CFunction parameter list.
    """

    block_body: str
    block_params: str
    global_body: str
    global_params: str


def tensor_family_of(name: str) -> Optional[Tuple[str, int]]:
    """
    Split an exact NRPy name into its tensor family and rank, if it has one.

    A component name ends with a variance run (``U``/``D``) followed by one
    index digit per variance letter, so ``Z4constraintU0`` is component 0 of the
    rank-1 family ``Z4constraintU`` while ``H_Z4`` is a scalar whose name merely
    ends in a digit.  Plain string methods decide this, so no regular
    expression is needed.

    :param name: Exact registered or factory-supplied name.
    :return: ``(family, rank)`` for a tensor component, or None for a scalar.

    Doctests:
    >>> tensor_family_of("Z4constraintU0")
    ('Z4constraintU', 1)
    >>> tensor_family_of("hDD01")
    ('hDD', 2)
    >>> print(tensor_family_of("H_Z4"))
    None
    >>> print(tensor_family_of("alpha"))
    None
    """
    index_run = len(name) - len(name.rstrip(_INDEX_DIGITS))
    if index_run == 0:
        return None
    family = name[: len(name) - index_run]
    if len(family) <= index_run:
        return None
    if not set(family[-index_run:]) <= _VARIANCE_LETTERS:
        return None
    return family, index_run


def build_diagnostics() -> FCCZ4DiagnosticsBuild:
    """
    Build the per-block and all-block constraint-diagnostic CFunction bodies.

    The formulation profile is read from the registered Dendro generation
    parameters, never from arguments: the registry is the single authority for
    it (whitepaper sections 4.1/6.1/9.3), and a builder that accepted its own
    copy could lower a kernel for a different reference metric than the frozen
    snapshot, the manifests and every sibling CFunction declare.

    :return: The immutable :class:`FCCZ4DiagnosticsBuild` result.
    :raises ValueError: If Infrastructure is not Dendro, if the shared factory
        supplies no diagnostics, if a diagnostic has no registered DIAG
        gridfunction, or if the kernel reads outside the evolved state.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fd_order", 4)
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    >>> par.set_parval_from_str("Dendro_enable_KO", False)
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _build = build_diagnostics()

    Section 14.6: the first diagnostic set is registered as DIAG, and the
    kernel writes exactly those gridfunctions.

    >>> reg.registered_diag_order()
    ('H_Z4', 'Z4constraintU0', 'Z4constraintU1', 'Z4constraintU2')
    >>> [
    ...     name
    ...     for name in reg.registered_diag_order()
    ...     if naming.diag_pointer(name) + "[pp] =" in _build.block_body
    ... ]
    ['H_Z4', 'Z4constraintU0', 'Z4constraintU1', 'Z4constraintU2']

    The connection constraint keeps its rank-1 metadata, and the scalar whose
    name ends in a digit is registered as an exact (non-base) name.

    >>> gri.glb_gridfcs_dict["Z4constraintU1"].rank
    1
    >>> gri.glb_gridfcs_dict["H_Z4"].is_basename
    False

    The emitted kernel is compared against the trusted source stored beside
    this module, so a change in the lowered text is caught here rather than by
    a standalone harness.

    >>> validate_strings(_build.block_body, "constraints_block", file_ext="cpp")
    >>> validate_strings(_build.global_body, "constraints_allblock", file_ext="cpp")
    """
    # Step 1: Require the qualified Dendro profile, and validate the registered
    # generation parameters before any expression is built.
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro diagnostics, "
            f"got {par.parval_from_str('Infrastructure')!r}."
        )
    require_serial_parallelization()
    generation_parameters.validate_generation_parameters()
    scalar_type = str(par.parval_from_str("Dendro_scalar_type"))
    fp_type = str(par.parval_from_str("fp_type"))

    # Step 2: Take the diagnostic expressions from the shared factory, so the
    # Dendro kernel and every other infrastructure lower the same expressions.
    bundle = build_fccz4_expression_bundle(
        CoordSystem=str(par.parval_from_str("Dendro_fccz4_CoordSystem")),
        LapseEvolutionOption=str(
            par.parval_from_str("Dendro_fccz4_LapseEvolutionOption")
        ),
        ShiftEvolutionOption=str(
            par.parval_from_str("Dendro_fccz4_ShiftEvolutionOption")
        ),
        enable_KreissOliger_dissipation=bool(par.parval_from_str("Dendro_enable_KO")),
        enable_diagnostics=True,
    )
    expressions: Dict[str, sp.Expr] = dict(bundle.diagnostics_by_name)
    if not expressions:
        raise ValueError(
            "The shared fCCZ4 factory supplied no diagnostics; PR 9 requires "
            "enable_diagnostics=True to produce H_Z4 and the connection "
            "constraint (section 14.6)."
        )

    # Step 3: Register the diagnostics as DIAG gridfunctions.  Components of
    # one tensor family register together so their rank is recorded correctly;
    # a name that is not a tensor component registers as an exact (non-base)
    # name, because a trailing digit is not a valid NRPy basename.
    # Registration is idempotent.
    families: Dict[str, int] = {}
    scalars: List[str] = []
    for name in sorted(expressions):
        family = tensor_family_of(name)
        if family is None:
            scalars.append(name)
            continue
        base, rank = family
        families[base] = rank
    for name in scalars:
        if name not in gri.glb_gridfcs_dict:
            gri.register_gridfunctions(name, group="DIAG", is_basename=False)
    for base, rank in sorted(families.items()):
        if f"{base}{'0' * rank}" not in gri.glb_gridfcs_dict:
            gri.register_gridfunctions_for_single_rankN(base, rank=rank, group="DIAG")
    diag_order = reg.registered_diag_order()
    missing = sorted(set(expressions) - set(diag_order))
    if missing:
        raise ValueError(
            f"Diagnostic expressions {missing} have no registered DIAG "
            "gridfunction; the exact-name rule (section 5.2) is violated."
        )
    written = tuple(name for name in diag_order if name in expressions)
    evol_order = reg.registered_evol_order()

    # Step 4: Lower the diagnostics inside an access capture, so the padding
    # comes from the reads the kernel actually emits.
    with cap.capture_gridfunction_accesses(CONSTRAINTS_BLOCK_CFUNCTION):
        kernel = c_codegen(
            [expressions[name] for name in written],
            [f"{naming.diag_pointer(name)}[pp]" for name in written],
            include_braces=False,
            enable_fd_codegen=True,
            enable_fd_functions=False,
            enable_simd=False,
            fp_type=fp_type,
            fp_type_alias=scalar_type,
            verbose=False,
        )
    accessed = set(cap.accessed_gridfunction_names(CONSTRAINTS_BLOCK_CFUNCTION))
    unexpected = sorted(accessed - set(evol_order))
    if unexpected:
        raise ValueError(
            "The diagnostic kernel may read only evolved state, but it also "
            f"read {unexpected}."
        )
    read_names = tuple(name for name in evol_order if name in accessed)

    # Step 5: Bind exactly the fields the capture recorded, plus the diagnostic
    # write targets, then wrap the kernel in the NRPy point and block loops.
    block_body = gfo.render_component_bindings(
        read_names,
        scalar_type,
        array="in_gfs",
        role=naming.input_pointer,
        const_pointee=True,
        index_expression=lambda name, _position: str(evol_order.index(name)),
    )
    block_body += "\n"
    block_body += gfo.render_component_bindings(
        written,
        scalar_type,
        array="diag_gfs",
        role=naming.diag_pointer,
        const_pointee=False,
        index_expression=lambda name, _position: str(diag_order.index(name)),
    )
    block_body += "\n"
    block_body += interior_loop(
        kernel,
        nx="geom.nx",
        ny="geom.ny",
        nz="geom.nz",
        padding="geom.padding",
        pmin_padded="geom.pmin_padded",
        dx="geom.dx",
    )
    block_params = (
        f"const BlockGeometry& geom, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* diag_gfs"
    )
    global_params = (
        f"const MockWorld& world, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* diag_gfs"
    )
    global_body = block_loop(
        f"{CONSTRAINTS_BLOCK_CFUNCTION}(world.geom[blk_id], in_gfs, diag_gfs);",
        count="world.num_blocks",
    )
    return FCCZ4DiagnosticsBuild(
        block_body=block_body,
        block_params=block_params,
        global_body=global_body,
        global_params=global_params,
    )


def register_diagnostics_CFunctions() -> None:
    """Register the per-block and all-block diagnostic CFunctions."""
    build = build_diagnostics()
    block_desc = (
        "Per-block fCCZ4 constraint diagnostics: the Hamiltonian constraint "
        "and the spatial Z4 connection constraint (recomputed, never "
        "checkpoint state)."
    )
    global_desc = "All-block fCCZ4 constraint diagnostics (NRPy block loop)."
    subdirectory = "generated/src/diagnostics"
    reg.register_dendro_CFunction(
        role="diagnostics_block",
        # Invoked by the all-block entry point's NRPy block loop, not by the
        # host, so it is not itself a Dendro entry point (section 4.5).
        entry_point=False,
        lifecycle_hook="diagnostic",
        name=CONSTRAINTS_BLOCK_CFUNCTION,
        desc=block_desc,
        subdirectory=subdirectory,
        params=build.block_params,
        body=build.block_body,
    )
    reg.register_dendro_CFunction(
        role="diagnostics",
        entry_point=True,
        calls=(CONSTRAINTS_BLOCK_CFUNCTION,),
        lifecycle_hook="diagnostic",
        name=CONSTRAINTS_GLOBAL_CFUNCTION,
        desc=global_desc,
        subdirectory=subdirectory,
        params=build.global_params,
        body=build.global_body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
