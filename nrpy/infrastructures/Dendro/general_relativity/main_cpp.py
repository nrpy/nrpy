# nrpy/infrastructures/Dendro/general_relativity/main_cpp.py
"""GR lifecycle policy for generated Dendro entry points and CTest gates."""

from typing import Tuple

from nrpy.infrastructures.Dendro import main_cpp as generic_main
from nrpy.infrastructures.Dendro.general_relativity.solver_context import (
    substitute_application_identifiers,
)

_STANDALONE_INITIALIZATION = r"""  if (ctx.minkowski_initial_data()) {
    MPI_Finalize();
    return 1;
  }  // END IF: GR initial data failed"""

_STANDALONE_BEFORE_STEPS = r"""  if (global_max(ctx.enforce_detgbar_equals_detghat_trAzero_all_blocks() != 0 ? 1.0 : 0.0) > 0.0) {
    if (rank == 0)
      std::fprintf(stderr, "FAIL: GR projection refused initial data\n");
    MPI_Finalize();
    return 1;
  }  // END IF: initial projection failed
  const double projection_residual = global_max(std::fmax(
      ctx.last_detgtrazero_status.max_abs_det_minus_one,
      ctx.last_detgtrazero_status.max_abs_trace_residual));
  if (rank == 0) std::printf("DETGTRAZERO_RESIDUAL %.3e\n", projection_residual);
  if (projection_residual > 1e-13) {
    if (rank == 0)
      std::fprintf(stderr, "FAIL: initial det/trace residual exceeds 1e-13\n");
    MPI_Finalize();
    return 1;
  }  // END IF: initial projection residual is too large
  const double constraint = global_max(ctx.max_constraint_violation());
  if (rank == 0) std::printf("MAXCONSTRAINT %.3e\n", constraint);
  if (constraint > 1e-12) {
    if (rank == 0)
      std::fprintf(stderr, "FAIL: initial constraint exceeds 1e-12\n");
    MPI_Finalize();
    return 1;
  }  // END IF: initial constraint is too large
  const double rhs = global_max(ctx.max_interior_rhs());
  if (rank == 0) std::printf("MINKOWSKIRHS %.3e\n", rhs);
  if (rhs > 1e-13) {
    if (rank == 0) std::fprintf(stderr, "FAIL: initial RHS exceeds 1e-13\n");
    MPI_Finalize();
    return 1;
  }  // END IF: Minkowski RHS is too large
  ctx.perturb_state();
  const double adapter_difference = global_max(ctx.flat_adapter_max_difference());
  if (rank == 0) std::printf("FLATADAPTER %.3e\n", adapter_difference);
  if (adapter_difference != 0.0) {
    if (rank == 0)
      std::fprintf(stderr, "FAIL: flat adapter disagrees with block RHS\n");
    MPI_Finalize();
    return 1;
  }  // END IF: flat adapter disagrees
  const double perturbed_rhs = global_max(ctx.max_interior_rhs());
  if (rank == 0) std::printf("PERTURBEDRHS %.3e\n", perturbed_rhs);
  if (!(perturbed_rhs > 1e-12)) {
    if (rank == 0)
      std::fprintf(stderr, "FAIL: perturbed RHS does not discriminate the kernel\n");
    MPI_Finalize();
    return 1;
  }  // END IF: perturbed RHS vanished
  const double order = $NAMESPACE::observed_convergence_order(dx, ctx.params);
  if (rank == 0) std::printf("ORDER %.3f\n", order);
  if (!(order >= static_cast<double>($NAMESPACE::generated::FD_ORDER) - 0.5)) {
    if (rank == 0)
      std::fprintf(stderr, "FAIL: observed convergence order %.3f is too low\n", order);
    MPI_Finalize();
    return 1;
  }  // END IF: convergence order is too low
  if (ctx.minkowski_initial_data() ||
      global_max(ctx.enforce_detgbar_equals_detghat_trAzero_all_blocks() != 0 ? 1.0 : 0.0) > 0.0 ||
      ctx.snapshot_state()) {
    MPI_Finalize();
    return 1;
  }  // END IF: reset before evolution failed"""

_STANDALONE_AFTER_STEP = r"""    if (global_max(ctx.enforce_detgbar_equals_detghat_trAzero_all_blocks() != 0 ? 1.0 : 0.0) > 0.0) {
      if (rank == 0)
        std::fprintf(stderr, "FAIL: GR projection refused an evolved point\n");
      MPI_Finalize();
      return 1;
    }  // END IF: accepted-step projection failed"""

_STANDALONE_FINAL_CHECKS = r"""  const double drift = global_max(ctx.max_drift_from_snapshot());
  if (rank == 0) std::printf("DRIFT100 %.3e\n", drift);
  if (drift > 1e-11) {
    if (rank == 0) std::fprintf(stderr, "FAIL: 100-step drift exceeds 1e-11\n");
    MPI_Finalize();
    return 1;
  }  // END IF: GR state drift is too large
  const int initial_data_constructions = 2;
  if (rank == 0)
    std::printf("DETGTRAZERO_PASSES %llu STEPS %d INITIALDATA %d\n",
                ctx.detgtrazero_passes, nsteps, initial_data_constructions);
  if (ctx.detgtrazero_passes !=
      static_cast<unsigned long long>(nsteps + initial_data_constructions)) {
    if (rank == 0)
      std::fprintf(stderr, "FAIL: GR projection ran the wrong number of times\n");
    MPI_Finalize();
    return 1;
  }  // END IF: projection schedule is wrong
  if (rank == 0)
    std::printf("MINKOWSKI_OK blocks=%d extent=%d ranks=%d\n",
                n_blocks, extent, size);"""

_REAL_FINAL_CHECKS = r"""      const double rhs = context.max_rhs();
      const double constraints = context.max_constraints();
      const double drift = context.max_drift();
      double residual = 0.0;
      MPI_Allreduce(&context.projection_residual, &residual, 1, MPI_DOUBLE,
                    MPI_MAX, MPI_COMM_WORLD);
      double local_spacing = std::numeric_limits<double>::max();
      if (mesh->isActive())
        for (const auto& block : mesh->getLocalBlockList()) {
          const auto geometry = $NAMESPACE::block_geometry(
              *mesh, block, minimum, maximum);
          for (double spacing : geometry.dx)
            local_spacing = std::min(local_spacing, spacing);
        }  // END LOOP: find smallest physical spacing
      double spacing = 0.0;
      MPI_Allreduce(&local_spacing, &spacing, 1, MPI_DOUBLE, MPI_MIN,
                    MPI_COMM_WORLD);
      const double derivative_tolerance =
          256 * std::numeric_limits<double>::epsilon() / (spacing * spacing);
      int local_ok =
          (!mesh->isActive() || context.projection_passes == 1 + 5ULL * steps) &&
          stepper.curr_step() == steps &&
          std::abs(stepper.curr_time() - dt * steps) <=
              1e-11 * std::max(1.0, dt * steps) &&
          rhs <= derivative_tolerance && constraints <= derivative_tolerance &&
          drift <= 1e-11 && residual <= 1e-13;
      int global_ok = 0;
      MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      int active = mesh->isActive(), active_ranks = 0;
      MPI_Allreduce(&active, &active_ranks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (rank == 0)
        std::printf("REAL_MINKOWSKI %s active_ranks=%d steps=%u "
                    "time=%.17g rhs=%.17g constraints=%.17g drift=%.17g "
                    "projection=%.17g hmin=%.17g derivative_tolerance=%.17g\n",
                    global_ok ? "PASS" : "FAIL", active_ranks, steps,
                    stepper.curr_time(), rhs, constraints, drift, residual,
                    spacing, derivative_tolerance);
      if (!global_ok)
        throw std::runtime_error("fixed-mesh Minkowski check failed");"""


def output_main_cpp(
    solver_stem: str, solver_namespace: str, exec_or_library_name: str
) -> str:
    """
    Emit the generic process shell with explicit GR lifecycle policy.

    :param solver_stem: Lowercase formulation stem used in emitted names.
    :param solver_namespace: Namespace containing the generated solver.
    :param exec_or_library_name: Generated executable target name.
    :return: Complete generated C++ entry-point source.
    """
    text = generic_main.output_main_cpp(
        solver_stem,
        solver_namespace,
        exec_or_library_name,
        _STANDALONE_INITIALIZATION,
        _STANDALONE_BEFORE_STEPS,
        _STANDALONE_AFTER_STEP,
        _STANDALONE_FINAL_CHECKS,
        "100",
        "0.5 * dx",
        _REAL_FINAL_CHECKS,
    )
    return substitute_application_identifiers(text)


def standalone_ctest_statements(
    solver_stem: str, exec_or_library_name: str
) -> Tuple[str, ...]:
    """
    Return the standalone GR lifecycle registration.

    :param solver_stem: Lowercase formulation stem used in test names.
    :param exec_or_library_name: Generated executable target name.
    :return: CMake statements registering the standalone lifecycle test.
    """
    return (
        "# The GR Minkowski lifecycle checks every printed acceptance gate.",
        f"add_test(NAME {solver_stem}_minkowski_lifecycle",
        f"         COMMAND {exec_or_library_name} -b 2 -n 25 -d 0.25)",
        f"set_tests_properties({solver_stem}_minkowski_lifecycle PROPERTIES TIMEOUT 300)",
    )


def real_ctest_statements(
    solver_stem: str, exec_or_library_name: str
) -> Tuple[str, ...]:
    """
    Return the real-host GR lifecycle registration.

    :param solver_stem: Lowercase formulation stem used in test names.
    :param exec_or_library_name: Generated executable target name.
    :return: CMake statements registering the real-host lifecycle test.
    """
    return (
        f"add_test(NAME {solver_stem}_real_minkowski COMMAND ${{MPIEXEC_EXECUTABLE}} ${{MPIEXEC_NUMPROC_FLAG}} 2 ${{MPIEXEC_PREFLAGS}} $<TARGET_FILE:{exec_or_library_name}> ${{MPIEXEC_POSTFLAGS}})",
        f"set_tests_properties({solver_stem}_real_minkowski PROPERTIES TIMEOUT 300)",
    )


if __name__ == "__main__":
    import doctest

    raise SystemExit(doctest.testmod().failed)
