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
from typing import Dict, List, Mapping

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
from nrpy.infrastructures.Dendro import gridfunction_name_decorations as gf_names
from nrpy.infrastructures.Dendro.general_relativity import generation_parameters
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
class ConstraintsEvalBuild:
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


def build_constraints_eval(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    enable_KreissOliger_dissipation: bool = False,
) -> ConstraintsEvalBuild:
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
    :raises ValueError: If the registered fields do not match the formulation.

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

    >>> tuple(gri.GridFunction.gridfunction_lists()[2])
    ('H_Z4', 'Z4constraintU0', 'Z4constraintU1', 'Z4constraintU2')
    >>> [
    ...     name
    ...     for name in gri.GridFunction.gridfunction_lists()[2]
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
    ...     _build = build_constraints_eval("bssn")

    The diagnostics register in DIAG, and the two suppressed AUX names were
    never registered:

    >>> tuple(gri.GridFunction.gridfunction_lists()[2])
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
            "Infrastructure must be 'Dendro' to build the Dendro diagnostics, "
            f"got {par.parval_from_str('Infrastructure')!r}."
        )
    require_serial_parallelization()
    generation_parameters.validate_generation_parameters()
    scalar_type = gri.DENDRO_SCALAR_TYPE
    fp_type = str(par.parval_from_str("fp_type"))

    if enable_fCCZ4:
        bundle = build_fccz4_expression_bundle(
            CoordSystem=CoordSystem,
            LapseEvolutionOption=LapseEvolutionOption,
            ShiftEvolutionOption=ShiftEvolutionOption,
            enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
            enable_diagnostics=True,
        )
        expressions: Dict[str, sp.Expr] = dict(bundle.diagnostics_by_name)
        if not expressions:
            raise ValueError("The shared fCCZ4 factory supplied no diagnostics.")
        diagnostic_names = tuple(expressions)
    else:
        diagnostic_names = ("H",) + tuple(f"MU{i}" for i in range(3))
    # BSSN requires H to exist as DIAG before its factory registers AUX fields.
    families: Dict[str, int] = {}
    scalars: List[str] = []
    for name in sorted(diagnostic_names):
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
    if not enable_fCCZ4:
        gate = "register_M_and_LAMBDA_CONSTRAINT_gridfunctions"
        previous_gate = par.parval_from_str(gate)
        par.set_parval_from_str(gate, False)
        try:
            constraints = BSSN_constraints[CoordSystem]
        finally:
            par.set_parval_from_str(gate, previous_gate)
        expressions = {"H": constraints.H}
        for i in range(3):
            expressions[f"MU{i}"] = constraints.MU[i]
    _evol, _auxevol, diag, _aux = gri.GridFunction.gridfunction_lists()
    diag_order = tuple(diag)
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
    block_body = bkh.output_component_bindings(
        read_names,
        scalar_type,
        array="in_gfs",
        role=gf_names.input_pointer,
        const_pointee=True,
        index_expression=lambda name, _position: str(evol_order.index(name)),
    )
    block_body += "\n"
    block_body += bkh.output_component_bindings(
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
        f"const block_geometry_struct& geom, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* diagnostic_gfs"
    )
    all_blocks_params = (
        f"const standalone_host_mesh_struct& mesh, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* diagnostic_gfs"
    )
    all_blocks_body = block_loop(
        f"{solver_stem}_{CONSTRAINTS_EVAL_BLOCK_SUFFIX}"
        f"(mesh.geom[blk], in_gfs, diagnostic_gfs);",
        num_blocks="mesh.num_blocks",
    )
    return ConstraintsEvalBuild(
        block_body=block_body,
        block_params=block_params,
        all_blocks_body=all_blocks_body,
        all_blocks_params=all_blocks_params,
        diagnostics_by_name=dict(expressions),
    )


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
    build = build_constraints_eval(
        solver_stem,
        enable_fCCZ4=enable_fCCZ4,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    subdirectory = "generated/src/diagnostics"
    includes = [f"{solver_stem}_defines.h"]
    cfunc_type = "void"
    block_name = f"{solver_stem}_{CONSTRAINTS_EVAL_BLOCK_SUFFIX}"
    if enable_fCCZ4:
        block_desc = (
            "Per-block fCCZ4 constraint diagnostics: the Hamiltonian constraint "
            "and the spatial Z4 connection constraint (recomputed, never "
            "checkpoint state)."
        )
    else:
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
    formulation = "fCCZ4" if enable_fCCZ4 else "BSSN"
    all_blocks_desc = (
        f"All-block {formulation} constraint diagnostics (NRPy block loop)."
    )
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
    from nrpy.equations.general_relativity.BSSN_quantities import (
        BSSN_quantities as SweepBSSNQuantities,
    )
    from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs as SweepBSSNRHSs
    from nrpy.equations.general_relativity.fCCZ4_constraints import (
        fCCZ4_constraints as SweepFCCZ4Constraints,
    )
    from nrpy.equations.general_relativity.fCCZ4_RHSs import (
        fCCZ4_RHSs as SweepFCCZ4RHSs,
    )
    from nrpy.infrastructures.Dendro.general_relativity import (
        rhs_eval as sweep_rhs_eval,
    )

    par.set_parval_from_str("Infrastructure", "Dendro")
    par.set_parval_from_str("parallelization", "none")
    par.set_parval_from_str("fp_type", "double")
    par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    par.set_parval_from_str("fd_order", 4)
    shipped_gauge = "OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted"
    for sweep_fCCZ4, sweep_cf in ((True, "chi"), (False, "W")):
        cfc.CFunction_dict.clear()
        gri.glb_gridfcs_dict.clear()
        par.glb_extras_dict.pop("Dendro", None)
        for factory in (
            SweepBSSNQuantities,
            SweepBSSNRHSs,
            BSSN_constraints,
            SweepFCCZ4RHSs,
            SweepFCCZ4Constraints,
        ):
            factory.clear()
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
            f"_{shipped_gauge}"
            f"_Cartesian_{sweep_cf}_fCCZ4{sweep_fCCZ4}",
            ve.process_dictionary_of_expressions(
                dict(sweep_build.diagnostics_by_name), fixed_mpfs_for_free_symbols=True
            ),
        )
