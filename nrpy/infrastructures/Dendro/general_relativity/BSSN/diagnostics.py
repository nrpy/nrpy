# nrpy/infrastructures/Dendro/general_relativity/BSSN/diagnostics.py
"""
BSSN constraint diagnostics for the Dendro infrastructure.

The diagnostic set is the Hamiltonian constraint and the three momentum
constraint components, taken from the established NRPy projector
:class:`nrpy.equations.general_relativity.BSSN_constraints.BSSNconstraints`, so
this module contributes no new formulation content: it registers the DIAG
gridfunctions and lowers those expressions into a Dendro point loop.

The diagnostics are recomputed from the evolved state and are never
authoritative checkpoint state, so they are registered in the DIAG group.  The
kernel uses the same finite-difference order and the same memory-access
mechanism as the right-hand side, so no separate diagnostic profile exists.

Every name is read back from the projector or from the registry; none is
written down here.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from dataclasses import dataclass
from typing import Dict, List

import sympy as sp

import nrpy.grid as gri
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.infrastructures.Dendro import (  # noqa: F401
    Dendro_state_h,
    generation_parameters,
)
from nrpy.infrastructures.Dendro import kernel_lowering as kl
from nrpy.infrastructures.Dendro import naming
from nrpy.infrastructures.Dendro import registration as reg
from nrpy.infrastructures.Dendro.block_loop import block_loop
from nrpy.infrastructures.Dendro.naming import tensor_family_of
from nrpy.infrastructures.Dendro.simple_loop import (
    require_serial_parallelization,
    simple_loop,
)

# Dendro names solver sources for the formulation; its own BSSN solver ships
# bssn_constraints.cpp, so the emitted CFunctions carry the same stem.
CONSTRAINTS_BLOCK_CFUNCTION = "bssn_constraints_block"
CONSTRAINTS_GLOBAL_CFUNCTION = "bssn_constraints"


@dataclass(frozen=True)
class BSSNDiagnosticsBuild:
    """
    Immutable result of building the BSSN constraint diagnostics.

    :param block_body: The per-block CFunction body.
    :param block_params: The per-block CFunction parameter list.
    :param global_body: The all-block CFunction body (NRPy block loop).
    :param global_params: The all-block CFunction parameter list.
    """

    block_body: str
    block_params: str
    global_body: str
    global_params: str


def build_diagnostics(*, CoordSystem: str = "Cartesian") -> BSSNDiagnosticsBuild:
    """
    Build the per-block and all-block BSSN constraint CFunction bodies.

    :param CoordSystem: Reference-metric coordinate system.
    :return: The immutable :class:`BSSNDiagnosticsBuild` result.
    :raises ValueError: If Infrastructure is not Dendro, if a diagnostic
        expression has no registered DIAG gridfunction, or if the kernel reads
        anything other than evolved state.
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro BSSN "
            f"diagnostics, got {par.parval_from_str('Infrastructure')!r}."
        )
    generation_parameters.validate_generation_parameters()
    require_serial_parallelization()
    scalar_type = str(par.parval_from_str("Dendro_scalar_type"))
    fp_type = str(par.parval_from_str("fp_type"))

    # Step 1: Register the constraint diagnostics as DIAG gridfunctions BEFORE
    # constructing the projector.  DIAG is the settled infrastructure group for
    # diagnostics (BHaH registers 31 of them across its wave-equation, elliptic
    # and GR diagnostics), whereas BSSN_constraints files H, M,
    # LAMBDA_CONSTRAINT and MU under AUX in the equations layer.  Every one of
    # those registrations is guarded by ``if <name> not in glb_gridfcs_dict``,
    # so registering the names this kernel writes first leaves the projector's
    # own registration a no-op and keeps the Dendro diagnostics contract
    # intact without touching the shared equations module.
    written_names = ("H",) + tuple(f"MU{i}" for i in range(3))
    families: Dict[str, int] = {}
    scalars: List[str] = []
    for name in sorted(written_names):
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

    # Step 2: Take the constraint expressions from the established NRPy
    # projector, which now finds its diagnostic names already registered.
    constraints = BSSN_constraints[CoordSystem]
    expressions: Dict[str, sp.Expr] = {"H": constraints.H}
    for i in range(3):
        expressions[f"MU{i}"] = constraints.MU[i]
    diag_order = reg.registered_diag_order()
    missing = sorted(set(expressions) - set(diag_order))
    if missing:
        raise ValueError(
            f"Diagnostic expressions {missing} have no registered DIAG "
            "gridfunction; the exact-name rule is violated."
        )
    written = tuple(name for name in diag_order if name in expressions)
    evol_order = reg.registered_evol_order()

    # Step 3: Lower the diagnostics.
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
    accessed = kl.accessed_gridfunctions(expressions[name] for name in written)
    unexpected = sorted(accessed - set(evol_order))
    if unexpected:
        raise ValueError(
            "The diagnostic kernel may read only evolved state, but it also "
            f"read {unexpected}."
        )
    read_names = tuple(name for name in evol_order if name in accessed)

    # Step 4: Bind exactly the fields the kernel reads, plus the diagnostic
    # write targets, then wrap the kernel in the NRPy point and block loops.
    block_body = Dendro_state_h.output_component_bindings(
        read_names,
        scalar_type,
        array="in_gfs",
        role=naming.input_pointer,
        const_pointee=True,
        index_expression=lambda name, _position: str(evol_order.index(name)),
    )
    block_body += "\n"
    block_body += Dendro_state_h.output_component_bindings(
        written,
        scalar_type,
        array="diagnostic_gfs",
        role=naming.diag_pointer,
        const_pointee=False,
        index_expression=lambda name, _position: str(diag_order.index(name)),
    )
    block_body += "\n"
    block_body += simple_loop(
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
        f"{scalar_type}* const* diagnostic_gfs"
    )
    global_params = (
        f"const MockWorld& world, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* diagnostic_gfs"
    )
    global_body = block_loop(
        f"{CONSTRAINTS_BLOCK_CFUNCTION}(world.geom[blk], in_gfs, diagnostic_gfs);",
        num_blocks="world.num_blocks",
    )
    return BSSNDiagnosticsBuild(
        block_body=block_body,
        block_params=block_params,
        global_body=global_body,
        global_params=global_params,
    )


def register_CFunctions_diagnostics(*, CoordSystem: str = "Cartesian") -> None:
    """
    Register the per-block and all-block BSSN diagnostic CFunctions.

    :param CoordSystem: Reference-metric coordinate system.
    """
    build = build_diagnostics(CoordSystem=CoordSystem)
    block_desc = (
        "Per-block BSSN constraint diagnostics: the Hamiltonian constraint and "
        "the three momentum constraint components (recomputed, never "
        "checkpoint state)."
    )
    global_desc = "All-block BSSN constraint diagnostics (NRPy block loop)."
    subdirectory = "generated/src/diagnostics"
    reg.register_Dendro_CFunction(
        role="diagnostics_block",
        name=CONSTRAINTS_BLOCK_CFUNCTION,
        desc=block_desc,
        subdirectory=subdirectory,
        params=build.block_params,
        body=build.block_body,
    )
    reg.register_Dendro_CFunction(
        role="diagnostics",
        name=CONSTRAINTS_GLOBAL_CFUNCTION,
        desc=global_desc,
        subdirectory=subdirectory,
        params=build.global_params,
        body=build.global_body,
    )
