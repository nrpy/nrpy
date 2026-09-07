# nrpy/infrastructures/Dendro/general_relativity/enforce_detgbar_equals_detghat_trAzero.py
r"""
Enforce det(gammabar) = det(gammahat) and tr(Abar) = 0 for the Dendro infrastructure.

The kernel restores the two algebraic constraints of the conformal
decomposition at every point of a block:

.. math::

    \frac{\det\bar\gamma}{\det\hat\gamma} = 1,
    \qquad
    \bar\gamma^{ij}\widetilde{A}_{ij} = 0.

The projected values come from the established NRPy module
:func:`nrpy.equations.general_relativity.BSSN_algebraic_constraints.BSSN_algebraic_constraints`,
so this module contributes no new formulation content: it lowers those
expressions into a Dendro point loop, adds a structured status record the
generated host lifecycle consumes, and never calls ``exit()``.

Every written field name is read back from the registered BSSN quantities
(``Bq.hDD[i][j]`` and ``Bq.aDD[i][j]`` are the gridfunction symbols themselves),
so no field name is hardcoded here.  ``lambdaU`` and ``Theta_fCCZ4`` are never
written: the enforcement is purely algebraic in ``hDD`` and ``aDD``.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from dataclasses import dataclass
from typing import List

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_algebraic_constraints import (
    BSSN_algebraic_constraints,
)
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import block_kernel_helpers as bkh
from nrpy.infrastructures.Dendro import generation_parameters
from nrpy.infrastructures.Dendro import gridfunction_name_decorations as gf_names
from nrpy.infrastructures.Dendro import state_h
from nrpy.infrastructures.Dendro.simple_loop import (
    block_loop,
    require_serial_parallelization,
)

# CFunction name suffixes; the solver stem is threaded from the caller.
ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK_SUFFIX = (
    "enforce_detgbar_equals_detghat_trAzero_block"
)
ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_ALL_BLOCKS_SUFFIX = (
    "enforce_detgbar_equals_detghat_trAzero"
)

# Generated status record: formulation-neutral (determinant, trace residual,
# nonfinite counts, first failing field/index, rank-local failure), so
# it belongs to the generated scalar contract emitted by types_h rather
# than to a physics builder.  The namespace is threaded from the caller, as
# every other emitted identifier is.
STATUS_RECORD = "generated::detgtrazero_status_struct"


@dataclass(frozen=True)
class DetgtrazeroBuild:
    """
    Immutable result of building the constraint enforcement for one profile.

    :param block_body: The per-block CFunction body (bindings + point loop).
    :param block_params: The per-block CFunction parameter list.
    :param all_blocks_body: The all-block CFunction body (NRPy block loop).
    :param all_blocks_params: The all-block CFunction parameter list.
    """

    block_body: str
    block_params: str
    all_blocks_body: str
    all_blocks_params: str


def build_enforce_detgbar_equals_detghat_trAzero(
    solver_stem: str, solver_namespace: str, *, CoordSystem: str = "Cartesian"
) -> DetgtrazeroBuild:
    """
    Build the per-block and all-block constraint-enforcement CFunction bodies.

    The point body evaluates the determinant ratio and the conformal trace of
    ``Atilde`` first, refuses the point when the determinant ratio is not
    positive or either quantity is nonfinite, and otherwise computes all
    twelve projected values into locals
    before writing any of them.  Computing into locals is what makes the
    in-place enforcement safe: the input and output pointers alias the same
    block arrays, so a value written early must not be able to perturb a value
    read late.

    :param solver_stem: Lowercase stem for the emitted CFunction names.
    :param solver_namespace: Solver namespace, following Dendro's lowercase
        formulation habit (``namespace bssn``).
    :param CoordSystem: Reference-metric coordinate system.
    :return: The immutable :class:`DetgtrazeroBuild` result.
    :raises ValueError: If Infrastructure is not Dendro, or if a projected
        field is not a registered EVOL gridfunction.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.grid as gri
    >>> from nrpy.equations.general_relativity.fCCZ4_system import (
    ...     build_fccz4_expression_bundle,
    ... )
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fd_order", 4)
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _bundle = build_fccz4_expression_bundle()
    >>> _build = build_enforce_detgbar_equals_detghat_trAzero(solver_stem="fccz4", solver_namespace="fccz4")

    The kernel writes exactly the rescaled conformal metric
    and traceless-curvature components, and nothing else.

    >>> sorted(
    ...     name
    ...     for name in gri.glb_gridfcs_dict
    ...     if gf_names.out_pointer(name) + "[pp] =" in _build.block_body
    ... )
    ['aDD00', 'aDD01', 'aDD02', 'aDD11', 'aDD12', 'aDD22', 'hDD00', 'hDD01', 'hDD02', 'hDD11', 'hDD12', 'hDD22']

    The connection and the Z4 scalar are left alone.

    >>> any(
    ...     gf_names.out_pointer(name) in _build.block_body
    ...     for name in ("lambdaU0", "lambdaU1", "lambdaU2", "Theta_fCCZ4")
    ... )
    False

    The failure branch never terminates the process.

    >>> "exit(" in _build.block_body
    False

    The emitted kernel itself is compared against the trusted baselines this
    module's ``__main__`` sweep captures, one file per profile.
    """
    # Step 1: Require the qualified Dendro profile, and validate the registered
    # generation parameters before any expression is built.
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro constraint enforcement, "
            f"got {par.parval_from_str('Infrastructure')!r}."
        )
    require_serial_parallelization()
    generation_parameters.validate_generation_parameters()
    scalar_type = gri.DENDRO_SCALAR_TYPE
    fp_type = str(par.parval_from_str("fp_type"))
    evol_order = roles.registered_evol_order()

    # Step 2: Collect the projected fields from BSSN_algebraic_constraints.
    # The target names are the registered BSSN gridfunction symbols themselves,
    # so the exact-name rule holds by construction.
    Bq = BSSN_quantities[CoordSystem]
    rfm = refmetric.reference_metric[CoordSystem]
    hprimeDD, aprimeDD = BSSN_algebraic_constraints(CoordSystem, False)
    projected_names: List[str] = []
    projected_exprs: List[sp.Expr] = []
    for tensor, projected in ((Bq.hDD, hprimeDD), (Bq.aDD, aprimeDD)):
        for i in range(3):
            for j in range(i, 3):
                projected_names.append(str(tensor[i][j]))
                projected_exprs.append(projected[i][j])
    unknown = sorted(set(projected_names) - set(evol_order))
    if unknown:
        raise ValueError(
            f"Projected fields {unknown} are not registered EVOL gridfunctions; "
            "the enforcement must write exact registered names."
        )

    # Step 3: Build the two reported residuals from the same registered BSSN
    # quantities the enforcement uses, so they describe the state it saw.  The
    # determinant is recomputed rather than read from ``Bq.detgammabar``, which
    # is the *assumed* value and would make the residual identically zero.
    _gammabarUU, detgammabar = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)
    det_ratio_expr = detgammabar / rfm.detgammahat
    trace_expr = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            trace_expr += Bq.gammabarUU[i][j] * Bq.AbarDD[i][j]

    # Step 4: Lower the residuals and the twelve projected values in ONE
    # c_codegen call.  Everything lands in locals,
    # so common subexpressions -- the determinant above all -- are shared
    # instead of evaluated twice per point.
    lvalues = [f"const {scalar_type} det_ratio", f"const {scalar_type} trace_residual"]
    lvalues += [f"const {scalar_type} projected_{name}" for name in projected_names]
    kernel = c_codegen(
        [det_ratio_expr, trace_expr] + projected_exprs,
        lvalues,
        include_braces=False,
        enable_fd_codegen=True,
        enable_fd_functions=False,
        enable_simd=False,
        fp_type=fp_type,
        fp_type_alias=scalar_type,
        verbose=False,
    )
    # Raw free symbols would miss a field read only through a derivative, so
    # this goes through the canonical derivative-aware reader.
    accessed = bkh.accessed_gridfunctions(
        [det_ratio_expr, trace_expr] + projected_exprs
    )
    read_names = tuple(name for name in evol_order if name in accessed)

    # Step 5: Assemble the point body top-to-bottom, matching the order the
    # generated C executes in: compute everything into locals, refuse the point
    # if the state is inadmissible, then store.  Storing last is what makes the
    # in-place enforcement safe, since the input and output pointers alias.
    point_body = kernel
    # Step 5.a: The structured-failure branch records the determinant, the
    # trace residual, the nonfinite flags, and the first failing field and
    # index, then leaves the point untouched.  It never calls exit(): the
    # rank-local count goes to the host, which owns the global reduction.  The
    # reported field is an index into the generated EVOL name array, so the
    # host can name it without a second table.  The projected values computed
    # above are simply not stored for a refused point.
    first_failing_field_scan = "".join(
        f"    if (status->first_failing_field < 0 && "
        f"!std::isfinite({gf_names.input_pointer(name)}[pp])) "
        f"status->first_failing_field = {evol_order.index(name)};\n"
        for name in read_names
    )
    point_body += r"""if (!(det_ratio > 0) || !std::isfinite(det_ratio) ||
    !std::isfinite(trace_residual)) {
  if (!std::isfinite(det_ratio) || !std::isfinite(trace_residual)) {
    status->nonfinite_points += 1;
  } // END IF: nonfinite determinant or trace
  if (status->failed_points == 0) {
    status->first_failing_index = static_cast<long long>(pp);
"""
    point_body += first_failing_field_scan
    point_body += r"""  } // END IF: first refused point
  status->failed_points += 1;
  continue;
} // END IF: nonpositive determinant or nonfinite value
"""
    point_body += r"""status->projected_points += 1;
status->max_abs_det_minus_one = std::fmax(
    status->max_abs_det_minus_one,
    static_cast<double>(std::fabs(det_ratio - 1)));
status->max_abs_trace_residual = std::fmax(
    status->max_abs_trace_residual,
    static_cast<double>(std::fabs(trace_residual)));
"""
    point_body += "".join(
        f"{gf_names.out_pointer(name)}[pp] = projected_{name};\n"
        for name in projected_names
    )

    # Step 6: Bind exactly the fields the kernel reads, plus the twelve
    # write targets.  Binding an unread pointer would trip -Wunused-variable.
    bindings = state_h.output_component_bindings(
        read_names,
        scalar_type,
        array="in_gfs",
        role=gf_names.input_pointer,
        const_pointee=True,
        index_expression=lambda name, _position: str(evol_order.index(name)),
    )
    bindings += "\n"
    bindings += state_h.output_component_bindings(
        projected_names,
        scalar_type,
        array="in_gfs",
        role=gf_names.out_pointer,
        const_pointee=False,
        index_expression=lambda name, _position: str(evol_order.index(name)),
    )

    # Step 7: Wrap the point body in the NRPy point loop and the NRPy block
    # loop.  Padding is zero: the enforcement is algebraic, so every cell of the
    # padded block is projected, ghost cells included.
    block_body = bindings + "\n"
    block_body += bkh.point_loop(point_body, padding="0")
    block_params = (
        f"const BlockGeometry& geom, {scalar_type}* const* in_gfs, "
        f"{solver_namespace}::{STATUS_RECORD}* const status"
    )
    all_blocks_body = block_loop(
        f"{solver_stem}_{ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK_SUFFIX}(mesh.geom[blk], in_gfs, status);",
        num_blocks="mesh.num_blocks",
    )
    all_blocks_params = (
        f"const StandaloneHostMesh& mesh, {scalar_type}* const* in_gfs, "
        f"{solver_namespace}::{STATUS_RECORD}* const status"
    )
    return DetgtrazeroBuild(
        block_body=block_body,
        block_params=block_params,
        all_blocks_body=all_blocks_body,
        all_blocks_params=all_blocks_params,
    )


def register_CFunctions_enforce_detgbar_equals_detghat_trAzero(
    solver_stem: str, solver_namespace: str, *, CoordSystem: str = "Cartesian"
) -> None:
    """
    Register the per-block and all-block constraint-enforcement CFunctions.

    :param solver_stem: Lowercase stem for the emitted CFunction names.
    :param solver_namespace: Solver namespace, following Dendro's lowercase
        formulation habit (``namespace bssn``).
    :param CoordSystem: Reference-metric coordinate system.
    """
    build = build_enforce_detgbar_equals_detghat_trAzero(
        solver_stem, solver_namespace, CoordSystem=CoordSystem
    )
    subdirectory = "generated/src/enforce_detgbar_equals_detghat_trAzero"
    includes = [f"{solver_stem}_defines.h"]
    cfunc_type = "void"
    block_name = f"{solver_stem}_{ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_BLOCK_SUFFIX}"
    block_desc = (
        "Enforce det(gammabar) = det(gammahat) and tr(Abar) = 0 per block: rescale the conformal metric to unit "
        "determinant ratio and remove the conformal trace of Atilde "
        "(structured status, no exit())."
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
    roles.set_CFunction_role(block_name, "enforce_detgbar_equals_detghat_trAzero_block")
    all_blocks_name = (
        f"{solver_stem}_{ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_ALL_BLOCKS_SUFFIX}"
    )
    all_blocks_desc = "Enforce det(gammabar) = det(gammahat) and tr(Abar) = 0 over all blocks (NRPy block loop)."
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=all_blocks_desc,
        cfunc_type=cfunc_type,
        name=all_blocks_name,
        params=build.all_blocks_params,
        body=build.all_blocks_body,
    )
    roles.set_CFunction_role(all_blocks_name, "enforce_detgbar_equals_detghat_trAzero")


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Trusted baseline for the emitted enforcement kernel, one file per shipped
    # profile.  general_relativity/initial_data.py's sweep carries the full
    # rationale for the capture object and the axes.
    from nrpy.helpers.generic import clang_format, validate_strings
    from nrpy.infrastructures.Dendro.general_relativity import trusted_capture

    par.set_parval_from_str("Infrastructure", "Dendro")
    par.set_parval_from_str("parallelization", "none")
    par.set_parval_from_str("fp_type", "double")
    par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    TRUSTED_FD_ORDER = 4
    par.set_parval_from_str("fd_order", TRUSTED_FD_ORDER)

    from nrpy.infrastructures.Dendro.general_relativity import (
        rhs_eval as sweep_rhs_eval,
    )

    for sweep_fCCZ4, sweep_cf in trusted_capture.SHIPPED_PROFILES:
        trusted_capture.reset_generation_state()
        sweep_stem = "fccz4" if sweep_fCCZ4 else "bssn"
        par.set_parval_from_str("EvolvedConformalFactor_cf", sweep_cf)
        # The right-hand-side registrar is what registers the evolved state
        # this enforcement rescales.
        sweep_rhs_eval.register_CFunctions_rhs_eval(
            enable_fCCZ4=sweep_fCCZ4,
            fd_order=TRUSTED_FD_ORDER,
            enable_KreissOliger_dissipation=False,
            solver_stem=sweep_stem,
        )
        register_CFunctions_enforce_detgbar_equals_detghat_trAzero(
            solver_stem=sweep_stem, solver_namespace=sweep_stem
        )
        trusted_kernel = cfc.CFunction_dict[
            roles.CFunction_name_for_role(
                "enforce_detgbar_equals_detghat_trAzero_block"
            )
        ]
        validate_strings(
            clang_format(trusted_kernel.full_function),
            f"Cartesian_{sweep_cf}_fCCZ4{sweep_fCCZ4}" f"_fdorder{TRUSTED_FD_ORDER}",
            file_ext="cpp",
        )
