# nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py
"""
BSSN or fCCZ4 constraint diagnostics for Dendro.

One module, one entry point, and the formulation chosen by the
``enable_fCCZ4`` argument.  The fCCZ4 branch takes H_Z4 and the Z4 connection
constraint from the shared ``fCCZ4_system`` factory; the BSSN branch takes the
Hamiltonian constraint and the momentum-constraint components from the
established ``BSSN_constraints`` factory.  Both register what they write as
DIAG gridfunctions, which is the Dendro diagnostics contract: they are
recomputed from the evolved state and are never checkpoint state.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from dataclasses import dataclass
from typing import Dict, List, Mapping, Union

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.equations.general_relativity.fCCZ4_system import (
    build_fccz4_expression_bundle,
)
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import block_kernel_helpers as bkh
from nrpy.infrastructures.Dendro import generation_parameters
from nrpy.infrastructures.Dendro import gridfunction_name_decorations as gf_names
from nrpy.infrastructures.Dendro import state_h
from nrpy.infrastructures.Dendro.gridfunction_name_decorations import tensor_family_of
from nrpy.infrastructures.Dendro.simple_loop import (
    block_loop,
    require_serial_parallelization,
)

# The emitted CFunction names are the established NRPy operation name,
# prefixed with the solver stem the caller supplies, as the sibling modules
# spell theirs.  The formulation is carried by the stem, so one suffix serves
# both formulations.
CONSTRAINTS_EVAL_BLOCK_SUFFIX = "constraints_eval_block"
CONSTRAINTS_EVAL_ALL_BLOCKS_SUFFIX = "constraints_eval"


@dataclass(frozen=True)
class FCCZ4ConstraintsEvalBuild:
    """
    Immutable result of building the diagnostic kernel for one profile.

    :param block_body: The per-block CFunction body (bindings + point loop).
    :param block_params: The per-block CFunction parameter list.
    :param all_blocks_body: The all-block CFunction body (NRPy block loop).
    :param all_blocks_params: The all-block CFunction parameter list.
    :param diagnostics_by_name: The assembled symbolic diagnostics, kept
        so the module's ``__main__`` can pin them against trusted values
        without reassembling them.
    """

    block_body: str
    block_params: str
    all_blocks_body: str
    all_blocks_params: str
    diagnostics_by_name: Mapping[str, sp.Expr]


def _build_constraints_eval_fCCZ4(
    solver_stem: str,
    *,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    enable_KreissOliger_dissipation: bool = False,
) -> FCCZ4ConstraintsEvalBuild:
    """
    Build the per-block and all-block constraint-diagnostic CFunction bodies.

    The formulation profile arrives as arguments, exactly as BHaH threads
    ``CoordSystem`` and ``LapseEvolutionOption`` into its own registration
    functions, so the diagnostics and the RHS lower the same profile because
    the caller passes the same values to both.

    :param solver_stem: Lowercase formulation stem, prefixed onto every
        emitted CFunction name.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :param enable_KreissOliger_dissipation: Passed to the shared factory so the
        diagnostics are built from the same expression set as the right-hand
        side; a per-call argument, as in BHaH and ETLegacy.
    :return: The immutable :class:`FCCZ4ConstraintsEvalBuild` result.
    :raises ValueError: If Infrastructure is not Dendro, if the shared factory
        supplies no diagnostics, if a diagnostic has no registered DIAG
        gridfunction, or if the kernel reads outside the evolved state.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fd_order", 4)
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    >>> from nrpy.equations.general_relativity.BSSN_constraints import (
    ...     BSSN_constraints as _bc,
    ... )
    >>> gri.glb_gridfcs_dict.clear()
    >>> par.glb_extras_dict.pop("Dendro", None) and None
    >>> _bc.clear()
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _build = build_constraints_eval("fccz4", enable_fCCZ4=True)

    The first diagnostic set is registered as DIAG, and the
    kernel writes exactly those gridfunctions.

    >>> roles.registered_diag_order()
    ('H_Z4', 'Z4constraintU0', 'Z4constraintU1', 'Z4constraintU2')
    >>> [
    ...     name
    ...     for name in roles.registered_diag_order()
    ...     if gf_names.diag_pointer(name) + "[pp] =" in _build.block_body
    ... ]
    ['H_Z4', 'Z4constraintU0', 'Z4constraintU1', 'Z4constraintU2']

    The connection constraint keeps its rank-1 metadata, and the scalar whose
    name ends in a digit is registered as an exact (non-base) name.

    >>> gri.glb_gridfcs_dict["Z4constraintU1"].rank
    1
    >>> gri.glb_gridfcs_dict["H_Z4"].is_basename
    False

    The emitted kernel is not captured as a trusted baseline: it is a large
    SymPy-lowered kernel, which ``coding_style.md`` excludes from golden-output
    files.  What pins it is the shared expression factory these diagnostics read
    from, which carries its own trusted expression dictionaries.
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
    scalar_type = gri.DENDRO_SCALAR_TYPE
    fp_type = str(par.parval_from_str("fp_type"))

    # Step 2: Take the diagnostic expressions from the shared factory, so the
    # Dendro kernel and every other infrastructure lower the same expressions.
    bundle = build_fccz4_expression_bundle(
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        enable_diagnostics=True,
    )
    expressions: Dict[str, sp.Expr] = dict(bundle.diagnostics_by_name)
    if not expressions:
        raise ValueError(
            "The shared fCCZ4 factory supplied no diagnostics; this builder needs "
            "enable_diagnostics=True to produce H_Z4 and the connection "
            "constraint."
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
    diag_order = roles.registered_diag_order()
    missing = sorted(set(expressions) - set(diag_order))
    if missing:
        raise ValueError(
            f"Diagnostic expressions {missing} have no registered DIAG "
            "gridfunction; the exact-name rule is violated."
        )
    written = tuple(name for name in diag_order if name in expressions)
    evol_order = roles.registered_evol_order()

    # Step 4: Lower the diagnostics.  The fields the kernel reads are the
    # expression free symbols that are registered gridfunctions.
    kernel = c_codegen(
        [expressions[name] for name in written],
        [f"{gf_names.diag_pointer(name)}[pp]" for name in written],
        include_braces=False,
        enable_fd_codegen=True,
        enable_fd_functions=False,
        enable_simd=False,
        fp_type=fp_type,
        fp_type_alias=scalar_type,
        verbose=False,
    )
    accessed = bkh.accessed_gridfunctions(expressions[name] for name in written)
    unexpected = sorted(accessed - set(evol_order))
    if unexpected:
        raise ValueError(
            "The diagnostic kernel may read only evolved state, but it also "
            f"read {unexpected}."
        )
    read_names = tuple(name for name in evol_order if name in accessed)

    # Step 5: Bind exactly the fields the kernel reads, plus the diagnostic
    # write targets, then wrap the kernel in the NRPy point and block loops.
    block_body = state_h.output_component_bindings(
        read_names,
        scalar_type,
        array="in_gfs",
        role=gf_names.input_pointer,
        const_pointee=True,
        index_expression=lambda name, _position: str(evol_order.index(name)),
    )
    block_body += "\n"
    block_body += state_h.output_component_bindings(
        written,
        scalar_type,
        array="diagnostic_gfs",
        role=gf_names.diag_pointer,
        const_pointee=False,
        index_expression=lambda name, _position: str(diag_order.index(name)),
    )
    block_body += "\n"
    block_body += bkh.point_loop(kernel)
    block_params = (
        f"const BlockGeometry& geom, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* diagnostic_gfs"
    )
    all_blocks_params = (
        f"const StandaloneHostMesh& mesh, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* diagnostic_gfs"
    )
    all_blocks_body = block_loop(
        f"{solver_stem}_{CONSTRAINTS_EVAL_BLOCK_SUFFIX}"
        f"(mesh.geom[blk], in_gfs, diagnostic_gfs);",
        num_blocks="mesh.num_blocks",
    )
    return FCCZ4ConstraintsEvalBuild(
        block_body=block_body,
        block_params=block_params,
        all_blocks_body=all_blocks_body,
        all_blocks_params=all_blocks_params,
        diagnostics_by_name=dict(expressions),
    )


def _register_CFunctions_constraints_eval_fCCZ4(
    solver_stem: str,
    *,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    enable_KreissOliger_dissipation: bool = False,
) -> None:
    """
    Register the per-block and all-block diagnostic CFunctions.

    :param solver_stem: Lowercase formulation stem, prefixed onto every
        emitted CFunction name and onto the emitted include.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :param enable_KreissOliger_dissipation: Forwarded to the builder so the
        diagnostics come from the same expression set as the right-hand side.
    """
    build = _build_constraints_eval_fCCZ4(
        solver_stem,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    subdirectory = "generated/src/diagnostics"
    includes = [f"{solver_stem}_defines.h"]
    cfunc_type = "void"
    block_name = f"{solver_stem}_{CONSTRAINTS_EVAL_BLOCK_SUFFIX}"
    block_desc = (
        "Per-block fCCZ4 constraint diagnostics: the Hamiltonian constraint "
        "and the spatial Z4 connection constraint (recomputed, never "
        "checkpoint state)."
    )
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=block_desc,
        cfunc_type=cfunc_type,
        name=block_name,
        params=build.block_params,
        body=build.block_body,
    )
    roles.set_CFunction_role(block_name, "constraints_eval_block")
    all_blocks_name = f"{solver_stem}_{CONSTRAINTS_EVAL_ALL_BLOCKS_SUFFIX}"
    all_blocks_desc = "All-block fCCZ4 constraint diagnostics (NRPy block loop)."
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=all_blocks_desc,
        cfunc_type=cfunc_type,
        name=all_blocks_name,
        params=build.all_blocks_params,
        body=build.all_blocks_body,
    )
    roles.set_CFunction_role(all_blocks_name, "constraints_eval")


@dataclass(frozen=True)
class BSSNConstraintsEvalBuild:
    """
    Immutable result of building the BSSN constraint diagnostics.

    :param block_body: The per-block CFunction body.
    :param block_params: The per-block CFunction parameter list.
    :param all_blocks_body: The all-block CFunction body (NRPy block loop).
    :param all_blocks_params: The all-block CFunction parameter list.
    :param diagnostics_by_name: The assembled symbolic diagnostics, kept
        so the module's ``__main__`` can pin them against trusted values
        without reassembling them.
    """

    block_body: str
    block_params: str
    all_blocks_body: str
    all_blocks_params: str
    diagnostics_by_name: Mapping[str, sp.Expr]


def _build_constraints_eval_BSSN(
    solver_stem: str, *, CoordSystem: str = "Cartesian"
) -> BSSNConstraintsEvalBuild:
    """
    Build the per-block and all-block BSSN constraint CFunction bodies.

    Either call order works: the factory's construction registers the
    evolved state through ``BSSN_quantities`` if the right-hand-side builder
    has not already done so, and the two AUX names this kernel does not write
    are suppressed at the source by the gate below rather than deleted
    afterwards, so nothing pre-existing is disturbed.

    :param solver_stem: Lowercase formulation stem, prefixed onto every
        emitted CFunction name.
    :param CoordSystem: Reference-metric coordinate system.
    :return: The immutable :class:`BSSNConstraintsEvalBuild` result.
    :raises ValueError: If Infrastructure is not Dendro, if a diagnostic
        expression has no registered DIAG gridfunction, or if the kernel reads
        anything other than evolved state.

    Doctests:
    >>> import contextlib, io
    >>> from nrpy.infrastructures.Dendro.general_relativity import rhs_eval
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
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _ = rhs_eval.build_rhs_eval(
    ...         "bssn", fd_order=4, enable_KreissOliger_dissipation=False)
    ...     _build = _build_constraints_eval_BSSN("bssn")

    The diagnostics register in DIAG, and the two suppressed AUX names were
    never registered:

    >>> roles.registered_diag_order()
    ('H', 'MU0', 'MU1', 'MU2')
    >>> [n for n in ("M", "LAMBDA_CONSTRAINT") if n in gri.glb_gridfcs_dict]
    []

    The kernel writes exactly those four and reads only evolved state:

    >>> sorted({line.split("[pp]")[0][len("diag_"):]
    ...         for line in _build.block_body.splitlines()
    ...         if line.strip().startswith("diag_")})
    ['H', 'MU0', 'MU1', 'MU2']
    >>> "in_trK" in _build.block_body, "aux_" in _build.block_body
    (True, False)
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro BSSN "
            f"diagnostics, got {par.parval_from_str('Infrastructure')!r}."
        )
    generation_parameters.validate_generation_parameters()
    require_serial_parallelization()
    scalar_type = gri.DENDRO_SCALAR_TYPE
    fp_type = str(par.parval_from_str("fp_type"))

    # Step 1: Register the constraint diagnostics as DIAG gridfunctions BEFORE
    # constructing the factory.  DIAG is the settled infrastructure group for
    # diagnostics (BHaH registers 31 of them across its wave-equation, elliptic
    # and GR diagnostics), whereas BSSN_constraints files H, M and
    # LAMBDA_CONSTRAINT under AUX in the equations layer, each guarded by
    # ``if <name> not in glb_gridfcs_dict``.  Registering H first therefore
    # makes the factory's registration of that name a no-op and keeps it in
    # DIAG.  MU is different: the factory registers it only under the
    # register_MU_gridfunctions CodeParameter, which defaults to False and
    # which no Dendro module sets, so nothing competes for those names and
    # this builder is their sole owner.  Neither case touches the shared
    # equations module.
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
    # factory.  It would otherwise also register M and LAMBDA_CONSTRAINT, which
    # this kernel does not compute and which would make the generated state
    # header advertise two variables no kernel writes and no vector backs.  The
    # core gate suppresses those two registrations rather than deleting them
    # from the registry afterwards; it is restored so no other caller inherits
    # this builder's choice.
    gate = "register_M_and_LAMBDA_CONSTRAINT_gridfunctions"
    previous_gate = par.parval_from_str(gate)
    par.set_parval_from_str(gate, False)
    try:
        constraints = BSSN_constraints[CoordSystem]
    finally:
        par.set_parval_from_str(gate, previous_gate)
    expressions: Dict[str, sp.Expr] = {"H": constraints.H}
    for i in range(3):
        expressions[f"MU{i}"] = constraints.MU[i]
    diag_order = roles.registered_diag_order()
    missing = sorted(set(expressions) - set(diag_order))
    if missing:
        raise ValueError(
            f"Diagnostic expressions {missing} have no registered DIAG "
            "gridfunction; the exact-name rule is violated."
        )
    written = tuple(name for name in diag_order if name in expressions)
    evol_order = roles.registered_evol_order()

    # Step 3: Lower the diagnostics.
    kernel = c_codegen(
        [expressions[name] for name in written],
        [f"{gf_names.diag_pointer(name)}[pp]" for name in written],
        include_braces=False,
        enable_fd_codegen=True,
        enable_fd_functions=False,
        enable_simd=False,
        fp_type=fp_type,
        fp_type_alias=scalar_type,
        verbose=False,
    )
    accessed = bkh.accessed_gridfunctions(expressions[name] for name in written)
    unexpected = sorted(accessed - set(evol_order))
    if unexpected:
        raise ValueError(
            "The diagnostic kernel may read only evolved state, but it also "
            f"read {unexpected}."
        )
    read_names = tuple(name for name in evol_order if name in accessed)

    # Step 4: Bind exactly the fields the kernel reads, plus the diagnostic
    # write targets, then wrap the kernel in the NRPy point and block loops.
    block_body = state_h.output_component_bindings(
        read_names,
        scalar_type,
        array="in_gfs",
        role=gf_names.input_pointer,
        const_pointee=True,
        index_expression=lambda name, _position: str(evol_order.index(name)),
    )
    block_body += "\n"
    block_body += state_h.output_component_bindings(
        written,
        scalar_type,
        array="diagnostic_gfs",
        role=gf_names.diag_pointer,
        const_pointee=False,
        index_expression=lambda name, _position: str(diag_order.index(name)),
    )
    block_body += "\n"
    block_body += bkh.point_loop(kernel)
    block_params = (
        f"const BlockGeometry& geom, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* diagnostic_gfs"
    )
    all_blocks_params = (
        f"const StandaloneHostMesh& mesh, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* diagnostic_gfs"
    )
    all_blocks_body = block_loop(
        f"{solver_stem}_{CONSTRAINTS_EVAL_BLOCK_SUFFIX}"
        f"(mesh.geom[blk], in_gfs, diagnostic_gfs);",
        num_blocks="mesh.num_blocks",
    )
    return BSSNConstraintsEvalBuild(
        block_body=block_body,
        block_params=block_params,
        all_blocks_body=all_blocks_body,
        all_blocks_params=all_blocks_params,
        diagnostics_by_name=dict(expressions),
    )


def _register_CFunctions_constraints_eval_BSSN(
    solver_stem: str,
    *,
    CoordSystem: str = "Cartesian",
) -> None:
    """
    Register the per-block and all-block BSSN diagnostic CFunctions.

    :param solver_stem: Lowercase formulation stem, prefixed onto every
        emitted CFunction name and onto the emitted include.
    :param CoordSystem: Reference-metric coordinate system.
    """
    build = _build_constraints_eval_BSSN(solver_stem, CoordSystem=CoordSystem)
    subdirectory = "generated/src/diagnostics"
    includes = [f"{solver_stem}_defines.h"]
    cfunc_type = "void"
    block_name = f"{solver_stem}_{CONSTRAINTS_EVAL_BLOCK_SUFFIX}"
    block_desc = (
        "Per-block BSSN constraint diagnostics: the Hamiltonian constraint and "
        "the three momentum constraint components (recomputed, never "
        "checkpoint state)."
    )
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=block_desc,
        cfunc_type=cfunc_type,
        name=block_name,
        params=build.block_params,
        body=build.block_body,
    )
    roles.set_CFunction_role(block_name, "constraints_eval_block")
    all_blocks_name = f"{solver_stem}_{CONSTRAINTS_EVAL_ALL_BLOCKS_SUFFIX}"
    all_blocks_desc = "All-block BSSN constraint diagnostics (NRPy block loop)."
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=all_blocks_desc,
        cfunc_type=cfunc_type,
        name=all_blocks_name,
        params=build.all_blocks_params,
        body=build.all_blocks_body,
    )
    roles.set_CFunction_role(all_blocks_name, "constraints_eval")


def build_constraints_eval(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    enable_KreissOliger_dissipation: bool = False,
) -> Union[BSSNConstraintsEvalBuild, FCCZ4ConstraintsEvalBuild]:
    """
    Build the Dendro constraint diagnostics for one formulation.

    One entry point taking the formulation as an argument, as
    ``BHaH/general_relativity/rhs_eval.py`` does for its own two formulations.

    :param solver_stem: Lowercase formulation stem, prefixed onto every
        emitted CFunction name.
    :param enable_fCCZ4: Build the fCCZ4 diagnostics instead of the BSSN ones.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option; ignored when
        ``enable_fCCZ4`` is False, where the expressions come from
        ``BSSN_constraints`` and carry no gauge choice.
    :param ShiftEvolutionOption: Shift evolution option; ignored on the same
        branch and for the same reason.
    :param enable_KreissOliger_dissipation: Forwarded to the fCCZ4 factory so
        the diagnostics come from the same expression set as the right-hand
        side; ignored when ``enable_fCCZ4`` is False, where the BSSN branch
        takes its expressions from ``BSSN_constraints``.
    :return: The formulation's build record.
    """
    if enable_fCCZ4:
        return _build_constraints_eval_fCCZ4(
            solver_stem,
            CoordSystem=CoordSystem,
            LapseEvolutionOption=LapseEvolutionOption,
            ShiftEvolutionOption=ShiftEvolutionOption,
            enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        )
    return _build_constraints_eval_BSSN(solver_stem, CoordSystem=CoordSystem)


def register_CFunctions_constraints_eval(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    enable_KreissOliger_dissipation: bool = False,
) -> None:
    """
    Register the constraint-diagnostic CFunctions for one formulation.

    :param solver_stem: Lowercase formulation stem, prefixed onto every
        emitted CFunction name and onto the emitted include.
    :param enable_fCCZ4: Register the fCCZ4 diagnostics instead of the BSSN ones.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :param enable_KreissOliger_dissipation: Forwarded to the fCCZ4 builder.
    """
    if enable_fCCZ4:
        _register_CFunctions_constraints_eval_fCCZ4(
            solver_stem,
            CoordSystem=CoordSystem,
            LapseEvolutionOption=LapseEvolutionOption,
            ShiftEvolutionOption=ShiftEvolutionOption,
            enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        )
        return
    _register_CFunctions_constraints_eval_BSSN(solver_stem, CoordSystem=CoordSystem)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Symbolic pinning of the two shipped diagnostic sets; rhs_eval.py's sweep
    # carries the rationale.  What is local here: the evolved state these
    # diagnostics read is registered by the right-hand-side builder, exactly as
    # in the generated projects, so that runs first.
    import os

    import nrpy.validate_expressions.validate_expressions as ve
    from nrpy.infrastructures.Dendro.general_relativity import (
        rhs_eval as sweep_rhs_eval,
    )
    from nrpy.infrastructures.Dendro.general_relativity import trusted_capture

    par.set_parval_from_str("Infrastructure", "Dendro")
    par.set_parval_from_str("parallelization", "none")
    par.set_parval_from_str("fp_type", "double")
    par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    par.set_parval_from_str("fd_order", 4)
    for sweep_fCCZ4, sweep_cf in trusted_capture.SHIPPED_PROFILES:
        trusted_capture.reset_generation_state()
        par.set_parval_from_str("EvolvedConformalFactor_cf", sweep_cf)
        _ = sweep_rhs_eval.build_rhs_eval(
            "fccz4" if sweep_fCCZ4 else "bssn",
            enable_fCCZ4=sweep_fCCZ4,
            fd_order=4,
            enable_KreissOliger_dissipation=False,
        )
        sweep_build = build_constraints_eval(
            "fccz4" if sweep_fCCZ4 else "bssn",
            enable_fCCZ4=sweep_fCCZ4,
            enable_KreissOliger_dissipation=False,
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}"
            f"_{trusted_capture.SHIPPED_GAUGE}"
            f"_Cartesian_{sweep_cf}_fCCZ4{sweep_fCCZ4}",
            ve.process_dictionary_of_expressions(
                dict(sweep_build.diagnostics_by_name), fixed_mpfs_for_free_symbols=True
            ),
        )
