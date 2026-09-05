# nrpy/infrastructures/Dendro/general_relativity/projection.py
r"""
PR 8: algebraic determinant/trace-free projection for the Dendro backend.

The projection restores the two algebraic constraints of the conformal
decomposition at every point of a block (whitepaper section 14.5):

.. math::

    \frac{\det\bar\gamma}{\det\hat\gamma} = 1,
    \qquad
    \bar\gamma^{ij}\widetilde{A}_{ij} = 0.

The projected values come from the established NRPy projector
:func:`nrpy.equations.general_relativity.BSSN_algebraic_constraints.BSSN_algebraic_constraints`,
so this module contributes no new formulation content: it lowers those
expressions into a Dendro point loop, adds the structured status record the
whitepaper requires, and never calls ``exit()``.

Every written field name is read back from the registered BSSN quantities
(``Bq.hDD[i][j]`` and ``Bq.aDD[i][j]`` are the gridfunction symbols themselves),
so no field name is hardcoded here.  ``lambdaU`` and ``Theta_fCCZ4`` are never
written: the projection is purely algebraic in ``hDD`` and ``aDD``.

Author: NRPy Dendro fCCZ4 infrastructure (PR 8)
"""

from dataclasses import dataclass
from typing import List

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_algebraic_constraints import (
    BSSN_algebraic_constraints,
)
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
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
PROJECTION_BLOCK_CFUNCTION = "fccz4_project_block"
PROJECTION_GLOBAL_CFUNCTION = "fccz4_project"

# Generated status type (fixed template ``fccz4_types_hpp.in``): the record is
# formulation-neutral (determinant, trace residual, nonfinite counts, floors,
# first failing field/index, rank-local failure), so it belongs to the
# generated scalar contract rather than to a physics builder.
STATUS_TYPE = "fccz4::generated::ProjectionStatus"


@dataclass(frozen=True)
class FCCZ4ProjectionBuild:
    """
    Immutable result of building the algebraic projection for one profile.

    :param block_body: The per-block CFunction body (bindings + point loop).
    :param block_params: The per-block CFunction parameter list.
    :param global_body: The all-block CFunction body (NRPy block loop).
    :param global_params: The all-block CFunction parameter list.
    """

    block_body: str
    block_params: str
    global_body: str
    global_params: str


def build_projection() -> FCCZ4ProjectionBuild:
    """
    Build the per-block and all-block algebraic projection CFunction bodies.

    The point body evaluates the determinant ratio and the conformal trace of
    ``Atilde`` first, refuses the point when either is nonpositive or
    nonfinite, and otherwise computes all twelve projected values into locals
    before writing any of them.  Computing into locals is what makes the
    in-place projection safe: the input and output pointers alias the same
    block arrays, so a value written early must not be able to perturb a value
    read late.

    :return: The immutable :class:`FCCZ4ProjectionBuild` result.
    :raises ValueError: If Infrastructure is not Dendro, or if a projected
        field is not a registered EVOL gridfunction.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.grid as gri
    >>> from nrpy.helpers.generic import validate_strings
    >>> from nrpy.equations.general_relativity.fCCZ4_system import (
    ...     build_fccz4_expression_bundle,
    ... )
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fd_order", 4)
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _bundle = build_fccz4_expression_bundle()
    >>> _build = build_projection()

    Section 14.5: the projection writes exactly the rescaled conformal metric
    and traceless-curvature components, and nothing else.

    >>> sorted(
    ...     name
    ...     for name in gri.glb_gridfcs_dict
    ...     if naming.out_pointer(name) + "[pp] =" in _build.block_body
    ... )
    ['aDD00', 'aDD01', 'aDD02', 'aDD11', 'aDD12', 'aDD22', 'hDD00', 'hDD01', 'hDD02', 'hDD11', 'hDD12', 'hDD22']

    Section 14.5 step 10: the connection and the Z4 scalar are left alone.

    >>> any(
    ...     naming.out_pointer(name) in _build.block_body
    ...     for name in ("lambdaU0", "lambdaU1", "lambdaU2", "Theta_fCCZ4")
    ... )
    False

    The failure branch never terminates the process (section 14.5 step 3).

    >>> "exit(" in _build.block_body
    False

    The emitted kernel is compared against the trusted source stored beside
    this module, so a change in the lowered text is caught here rather than by
    a standalone harness.

    >>> validate_strings(_build.block_body, "projection_block", file_ext="cpp")
    >>> validate_strings(_build.global_body, "projection_allblock", file_ext="cpp")
    """
    # Step 1: Require the qualified Dendro profile, and validate the registered
    # generation parameters before any expression is built.
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro projection, "
            f"got {par.parval_from_str('Infrastructure')!r}."
        )
    require_serial_parallelization()
    generation_parameters.validate_generation_parameters()
    scalar_type = str(par.parval_from_str("Dendro_scalar_type"))
    fp_type = str(par.parval_from_str("fp_type"))
    CoordSystem = str(par.parval_from_str("Dendro_fccz4_CoordSystem"))
    evol_order = reg.registered_evol_order()

    # Step 2: Collect the projected fields from the established NRPy projector.
    # The target names are the registered BSSN gridfunction symbols themselves,
    # so the section 5.2 exact-name rule holds by construction.
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
            "the projection must write exact registered names (section 5.2)."
        )

    # Step 3: Build the two reported residuals from the same registered BSSN
    # quantities the projector uses, so they describe the state it saw.  The
    # determinant is recomputed rather than read from ``Bq.detgammabar``, which
    # is the *assumed* value and would make the residual identically zero.
    _gammabarUU, detgammabar = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)
    det_ratio_expr = detgammabar / rfm.detgammahat
    trace_expr = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            trace_expr += Bq.gammabarUU[i][j] * Bq.AbarDD[i][j]

    # Step 4: Lower the residuals and the twelve projected values in ONE
    # c_codegen call, inside one access capture.  Everything lands in locals,
    # so common subexpressions -- the determinant above all -- are shared
    # instead of evaluated twice per point.
    lvalues = [f"const {scalar_type} det_ratio", f"const {scalar_type} trace_residual"]
    lvalues += [f"const {scalar_type} projected_{name}" for name in projected_names]
    with cap.capture_gridfunction_accesses(PROJECTION_BLOCK_CFUNCTION):
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
    accessed = set(cap.accessed_gridfunction_names(PROJECTION_BLOCK_CFUNCTION))
    read_names = tuple(name for name in evol_order if name in accessed)

    # Step 5: Assemble the point body top-to-bottom, matching the order the
    # generated C executes in: compute everything into locals, refuse the point
    # if the state is inadmissible, then store.  Storing last is what makes the
    # in-place projection safe, since the input and output pointers alias.
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
        f"!std::isfinite({naming.input_pointer(name)}[pp])) "
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
} // END IF: nonpositive or nonfinite determinant
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
        f"{naming.out_pointer(name)}[pp] = projected_{name};\n"
        for name in projected_names
    )

    # Step 6: Bind exactly the fields the capture recorded, plus the twelve
    # write targets.  Binding an unread pointer would trip -Wunused-variable.
    bindings = gfo.render_component_bindings(
        read_names,
        scalar_type,
        array="state_gfs",
        role=naming.input_pointer,
        const_pointee=True,
        index_expression=lambda name, _position: str(evol_order.index(name)),
    )
    bindings += "\n"
    bindings += gfo.render_component_bindings(
        projected_names,
        scalar_type,
        array="state_gfs",
        role=naming.out_pointer,
        const_pointee=False,
        index_expression=lambda name, _position: str(evol_order.index(name)),
    )

    # Step 7: Wrap the point body in the NRPy point loop and the NRPy block
    # loop.  Padding is zero: the projection is algebraic, so every cell of the
    # padded block is projected, ghost cells included.
    block_body = bindings + "\n"
    block_body += interior_loop(
        point_body,
        nx="geom.nx",
        ny="geom.ny",
        nz="geom.nz",
        padding="0",
        pmin_padded="geom.pmin_padded",
        dx="geom.dx",
    )
    block_params = (
        f"const BlockGeometry& geom, {scalar_type}* const* state_gfs, "
        f"{STATUS_TYPE}* const status"
    )
    global_body = block_loop(
        f"{PROJECTION_BLOCK_CFUNCTION}(world.geom[blk_id], state_gfs, status);",
        count="world.num_blocks",
    )
    global_params = (
        f"const MockWorld& world, {scalar_type}* const* state_gfs, "
        f"{STATUS_TYPE}* const status"
    )
    return FCCZ4ProjectionBuild(
        block_body=block_body,
        block_params=block_params,
        global_body=global_body,
        global_params=global_params,
    )


def register_projection_CFunctions() -> None:
    """Register the per-block and all-block projection CFunctions."""
    build = build_projection()
    block_desc = (
        "Per-block algebraic projection: rescale the conformal metric to unit "
        "determinant ratio and remove the conformal trace of Atilde "
        "(structured status, no exit())."
    )
    global_desc = "All-block algebraic projection (NRPy block loop)."
    subdirectory = "generated/src/projection"
    reg.register_dendro_CFunction(
        role="projection_block",
        # Invoked by the all-block entry point's NRPy block loop, not by the
        # host, so it is not itself a Dendro entry point (section 4.5).
        entry_point=False,
        lifecycle_hook="post_timestep",
        name=PROJECTION_BLOCK_CFUNCTION,
        desc=block_desc,
        subdirectory=subdirectory,
        params=build.block_params,
        body=build.block_body,
    )
    reg.register_dendro_CFunction(
        role="projection",
        entry_point=True,
        calls=(PROJECTION_BLOCK_CFUNCTION,),
        lifecycle_hook="post_timestep",
        name=PROJECTION_GLOBAL_CFUNCTION,
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
