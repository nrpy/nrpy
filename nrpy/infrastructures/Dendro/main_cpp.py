# nrpy/infrastructures/Dendro/main_cpp.py
"""
Emit the standalone-host entry point for a generated Dendro solver.

The executable drives the registered generated CFunctions through the
NRPy-supplied standalone host declarations, so the Minkowski lifecycle is
exercised end to end.  Gates, in order:

    DETGTRAZERO_RESIDUAL  max det/trace residual after initial data        (<= 1e-13)
    MAXCONSTRAINT max |constraint diagnostic| after initial data   (<= 1e-12)
    MINKOWSKIRHS  max |RHS| at the Minkowski fixed point           (<= 1e-13)
    FLATADAPTER   |per-block RHS - flat-block adapter RHS|         (== 0)
    PERTURBEDRHS  max |RHS| on a smooth perturbed state            (> 0)
    ORDER         observed convergence order under h -> h/2 -> h/4 (>= N-0.5)
    DRIFT100      max drift after 100 CFL-limited steps            (<= 1e-11)
    DETGTRAZERO_PASSES    one enforcement per initial-data construction plus
                  one per accepted step                            (exact)

Every gate but ORDER and DETGTRAZERO_PASSES is an ``MPI_Allreduce(MAX)`` over ranks,
and each rank owns a disjoint subdomain carrying a different piece of the
analytic profile, so a rank-dependent fault cannot hide behind a rank-0 print.
ORDER is a single-block refinement study and is rank-independent by
construction.  DETGTRAZERO_PASSES compares the rank-local counter, because every rank
runs the same schedule and an exact per-rank count is the stronger check.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.solver_context import (
    substitute_solver_identifiers,
)

BANNER = generated_file_banner()

_MAIN = """// Standalone-host entry point.  The real Dendro-GR main (parameter file, Dendro
// mesh, MPI decomposition) lands when the generated solver is built against a
// real Dendro-GR checkout.
//
// Usage: $EXEC [-b n_blocks] [-n extent] [-d dx] [-t parfile] [-r name]...
//    -b number of local blocks per rank
//    -n padded block extent per axis (>= 2*REQUIRED_PADDING + 1)
//    -d block spacing
//    -t parameter-file path (rejected while this profile has no binding)
//    -r exact NRPy variable name to select for output/refinement (repeatable;
//       an unknown name is fatal and every valid generated name is listed)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

#include "$STEMCtx.h"

namespace {

double global_max(double local) {
  double global = local;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return global;
}  // END FUNCTION: global_max

// clang-format off
}  // END NAMESPACE: internal linkage
// clang-format on

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int n_blocks = 1, extent = 25;
  double dx = 0.25;
  const char* parfile = nullptr;
  // Exact NRPy names selected for output/refinement.  The names arrive from
  // the host command line; none is written down in this emitter.
  const int MAX_SELECTED = 16;
  const char* selected[16];
  unsigned n_selected = 0;
  for (int i = 1; i < argc; ++i) {
    if (!std::strcmp(argv[i], "-b") && i + 1 < argc) {
      n_blocks = std::atoi(argv[++i]);
    } else if (!std::strcmp(argv[i], "-n") && i + 1 < argc) {
      extent = std::atoi(argv[++i]);
    } else if (!std::strcmp(argv[i], "-d") && i + 1 < argc) {
      dx = std::atof(argv[++i]);
    } else if (!std::strcmp(argv[i], "-t") && i + 1 < argc) {
      parfile = argv[++i];
    } else if (!std::strcmp(argv[i], "-r") && i + 1 < argc) {
      if (n_selected < MAX_SELECTED) {
        selected[n_selected++] = argv[++i];
      } else {
        if (rank == 0) {
          std::fprintf(stderr, "too many -r selections (max %d)\\n",
                       MAX_SELECTED);
        }  // END IF: rank 0 reports it
        MPI_Finalize();
        return 2;
      }  // END ELSE: selection table full
    } else {  // END ELSE IF: -r name selection
      if (rank == 0) {
        std::fprintf(stderr,
                     "usage: $EXEC [-b n_blocks] [-n extent] [-d dx] "
                     "[-t parfile] [-r exact_name]...\\n");
      }  // END IF: rank 0 reports it
      MPI_Finalize();
      return 2;
    }  // END ELSE: usage and exit
  }  // END LOOP: for i over command-line arguments

  $NAMESPACE::Ctx ctx;
  if (ctx.initialize_mesh(n_blocks, extent, dx, rank, parfile)) {
    MPI_Finalize();
    return 1;
  }  // END IF: mesh setup failed
  if (ctx.startup_checks(rank)) {
    MPI_Finalize();
    return 1;
  }  // END IF: startup checks failed
  if (ctx.minkowski_initial_data()) {
    MPI_Finalize();
    return 1;
  }  // END IF: initial data failed
  // Resolve any host-selected exact names before evolving, so a misspelled
  // output or refinement variable fails immediately.
  if (n_selected > 0 && ctx.select_variables(selected, n_selected) != 0) {
    MPI_Finalize();
    return 1;
  }  // END IF: name selection failed

  // Project after initial-data construction.  The failure decision is
  // collective: an enforcement failure is data dependent, so a rank-local
  // `return` here would leave the other ranks waiting inside the next
  // reduction.  Every rank reduces the same flag and then agrees.
  if (global_max(ctx.enforce_detgbar_equals_detghat_trAzero_all_blocks() != 0 ? 1.0 : 0.0) > 0.0) {
    if (rank == 0) {
      std::fprintf(stderr,
                   "FAIL: constraint enforcement refused a point after initial data\\n");
    }  // END IF: rank 0 reports the refusal
    MPI_Finalize();
    return 1;
  }  // END IF: enforcement refused after initial data
  // Measured on flat initial data, so it is identically zero whatever the
  // kernel computes.  The discriminating enforcement evidence is the
  // detgtrazero self-test and the owner doctests that pin the exact
  // projected write set, not this line.
  const double detgtrazero_residual = global_max(
      std::fmax(ctx.last_detgtrazero_status.max_abs_det_minus_one,
                ctx.last_detgtrazero_status.max_abs_trace_residual));
  if (rank == 0) std::printf("DETGTRAZERO_RESIDUAL %.3e\\n", detgtrazero_residual);
  if (detgtrazero_residual > 1e-13) {
    if (rank == 0) {
      std::fprintf(stderr, "FAIL: initial det/trace residual exceeds 1e-13\\n");
    }  // END IF: rank 0 reports the residual
    MPI_Finalize();
    return 1;
  }  // END IF: initial residual too large

  // The constraint diagnostics of an exact solution vanish.  A kernel that
  // computed nothing would also report zero, so this gate is a necessary
  // condition only, and it establishes no pointwise value: the equations
  // layer's trusted dictionaries pin a different construction profile than
  // this solver lowers, so none of these expressions is pinned there.
  const double max_constraint = global_max(ctx.max_constraint_violation());
  if (rank == 0) std::printf("MAXCONSTRAINT %.3e\\n", max_constraint);
  if (max_constraint > 1e-12) {
    if (rank == 0) {
      std::fprintf(stderr, "FAIL: initial constraint violation exceeds 1e-12\\n");
    }  // END IF: rank 0 reports the violation
    MPI_Finalize();
    return 1;
  }  // END IF: initial constraint violation too large

  // Gate 1: the Minkowski state is a fixed point.
  const double rhs_max = global_max(ctx.max_interior_rhs());
  if (rank == 0) std::printf("MINKOWSKIRHS %.3e\\n", rhs_max);
  if (rhs_max > 1e-13) {
    if (rank == 0) std::fprintf(stderr, "FAIL: initial RHS exceeds 1e-13\\n");
    MPI_Finalize();
    return 1;
  }  // END IF: initial RHS too large

  // Gate 2: the LTS flat-block adapter and the per-block entry point are one
  // numerical body, so they must agree exactly.  Checked on a perturbed state,
  // where the RHS is not identically zero.
  // The perturbation runs on the effective parameters the context printed at
  // startup, so what the gates measure is what the run reports.
  ctx.perturb_state();
  const double adapter_difference =
      global_max(ctx.flat_adapter_max_difference());
  if (rank == 0) std::printf("FLATADAPTER %.3e\\n", adapter_difference);
  if (adapter_difference != 0.0) {
    if (rank == 0) {
      std::fprintf(stderr,
                   "FAIL: flat-block adapter disagrees with the per-block "
                   "entry point\\n");
    }  // END IF: rank 0 reports it
    MPI_Finalize();
    return 1;
  }  // END IF: flat adapter disagrees

  // Gate 3: on a spatially varying state the generated RHS must not vanish.  A
  // Minkowski state alone cannot distinguish a working kernel from one whose
  // derivative terms are absent, because every stencil difference of a
  // constant field is exactly zero whatever the coefficients.
  const double perturbed_rhs = global_max(ctx.max_interior_rhs());
  if (rank == 0) std::printf("PERTURBEDRHS %.3e\\n", perturbed_rhs);
  if (!(perturbed_rhs > 1e-12)) {
    if (rank == 0) {
      std::fprintf(stderr,
                   "FAIL: the RHS vanishes on a perturbed state; the "
                   "generated derivative terms are not being evaluated\\n");
    }  // END IF: rank 0 reports it
    MPI_Finalize();
    return 1;
  }  // END IF: perturbed RHS vanished

  // Gate 4: the RHS converges under refinement at the requested
  // finite-difference order.  This catches stencils that are mis-shaped or
  // applied at the wrong offsets.  It does NOT catch a uniform sign or scale
  // error in the coefficients: that kernel approximates a different continuum
  // operator and converges just as well.
  const double order = $NAMESPACE::observed_convergence_order(dx, ctx.params);
  if (rank == 0) std::printf("ORDER %.3f\\n", order);
  const double expected_order =
      static_cast<double>($NAMESPACE::generated::FD_ORDER) - 0.5;
  if (!(order >= expected_order)) {
    if (rank == 0) {
      std::fprintf(stderr,
                   "FAIL: observed convergence order %.3f is below the "
                   "requested %u - 0.5\\n",
                   order, $NAMESPACE::generated::FD_ORDER);
    }  // END IF: rank 0 reports it
    MPI_Finalize();
    return 1;
  }  // END IF: convergence order too low

  // Gate 5: 100 CFL-limited steps from the Minkowski fixed point must not
  // drift.
  if (ctx.minkowski_initial_data()) {
    MPI_Finalize();
    return 1;
  }  // END IF: initial data failed
  // The same scheduling rule applies to this second initial-data
  // construction, so the pass count below stays exactly one per construction
  // plus one per accepted step.
  if (global_max(ctx.enforce_detgbar_equals_detghat_trAzero_all_blocks() != 0 ? 1.0 : 0.0) > 0.0) {
    if (rank == 0) {
      std::fprintf(stderr,
                   "FAIL: constraint enforcement refused a point after initial data\\n");
    }  // END IF: rank 0 reports it
    MPI_Finalize();
    return 1;
  }  // END IF: enforcement refused a point
  const int initial_data_constructions = 2;
  ctx.snapshot_state();
  const int nsteps = 100;
  const double dt = 0.5 * dx;  // CFL-limited for unit wavespeed.
  for (int s = 0; s < nsteps; ++s) {
    if (ctx.euler_step(dt)) {
      MPI_Finalize();
      return 1;
    }  // END IF: timestep failed
    // Project after every accepted timestep, which is what the post_timestep
    // hook does on the real host.  The failure decision is collective, for the
    // reason given above.
    if (global_max(ctx.enforce_detgbar_equals_detghat_trAzero_all_blocks() != 0 ? 1.0 : 0.0) > 0.0) {
      if (rank == 0) {
        std::fprintf(stderr,
                     "FAIL: constraint enforcement refused a point during evolution\\n");
      }  // END IF: rank 0 reports the refusal
      MPI_Finalize();
      return 1;
    }  // END IF: enforcement refused during evolution
  }  // END LOOP: for s over CFL-limited steps
  const double drift = global_max(ctx.max_drift_from_snapshot());
  if (rank == 0) std::printf("DRIFT100 %.3e\\n", drift);
  if (drift > 1e-11) {
    if (rank == 0) std::fprintf(stderr, "FAIL: 100-step drift exceeds 1e-11\\n");
    MPI_Finalize();
    return 1;
  }  // END IF: 100-step drift too large
  // The enforcement hooks ran exactly as configured -- once per initial-data
  // construction and once per accepted timestep.  The count is incremented by
  // the context when the generated CFunction actually runs.
  if (rank == 0) {
    std::printf("DETGTRAZERO_PASSES %llu STEPS %d INITIALDATA %d\\n",
                ctx.detgtrazero_passes, nsteps, initial_data_constructions);
  }  // END IF: rank 0 reports pass count
  if (ctx.detgtrazero_passes !=
      static_cast<unsigned long long>(nsteps + initial_data_constructions)) {
    if (rank == 0) {
      std::fprintf(stderr, "FAIL: constraint enforcement ran %llu times, expected %d\\n",
                   ctx.detgtrazero_passes,
                   nsteps + initial_data_constructions);
    }  // END IF: rank 0 reports the count
    MPI_Finalize();
    return 1;
  }  // END IF: wrong enforcement pass count
  if (rank == 0) {
    std::printf("MINKOWSKI_OK blocks=%d extent=%d ranks=%d\\n", n_blocks,
                extent, size);
  }  // END IF: rank 0 reports it
  MPI_Finalize();
  return 0;
}  // END FUNCTION: main
"""


def output_main_cpp(
    solver_stem: str,
    solver_namespace: str,
    exec_or_library_name: str,
) -> str:
    """
    Emit the standalone-host entry point source.

    :param solver_stem: Lowercase formulation stem for emitted file names.
    :param solver_namespace: Solver namespace.
    :param exec_or_library_name: Name of the solver executable target, used in
        the usage message.
    :return: The complete C++ source text.

    Doctests:
    >>> from nrpy.infrastructures.Dendro.clang_format_guards import (
    ...     unguarded_end_namespace_markers,
    ... )
    >>> unguarded_end_namespace_markers(_MAIN)
    []
    >>> _MAIN.count("}  // END NAMESPACE:")
    1
    """
    text = _MAIN.replace("$EXEC", exec_or_library_name)
    return BANNER + substitute_solver_identifiers(text, solver_stem, solver_namespace)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
