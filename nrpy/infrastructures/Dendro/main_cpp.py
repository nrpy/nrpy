# nrpy/infrastructures/Dendro/main_cpp.py
"""
Emit standalone and real-host entry points for a generated Dendro solver.

The executable drives the registered generated CFunctions through the
NRPy-supplied standalone host declarations or actual Dendro runtime types.
The standalone branch checks these gates, in order:

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

The real branch checks fixed-mesh Minkowski evolution through Dendro ETS RK4,
with runtime TOML binding and a roundoff-scaled derivative tolerance. Its
transport and parameter-response oracles live in tests_infra.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from nrpy.infrastructures.Dendro.CodeParameters import output_toml_bindings
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.solver_context import (
    substitute_solver_identifiers,
)

BANNER = generated_file_banner()

_MAIN = """// Standalone-host entry point. The real Dendro-GR build selects the separate
// entry point below, with parameter files and actual distributed mesh types.
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


_REAL_MAIN = r"""#include "$STEMCtx.h"
#include "ets.h"
#include "meshUtils.h"
#include "octUtils.h"
#include <toml.hpp>
#include <fstream>
#include <sstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
/**
 * Run a fixed-mesh Minkowski qualification through the pinned real ETS host.
 *
 * @param argc Number of command-line arguments.
 * @param[in,out] argv Argument vector parsed by MPI and this entry point.
 * @return 0 on success, 1 if MPI_Abort unexpectedly returns after a failure.
 *
 * @note Detected failures abort MPI_COMM_WORLD, including inactive ranks.
 */
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    unsigned steps = 100;
    double dt = 0.001;
    std::string filename;
    for (int a = 1; a < argc; ++a) {
      const std::string arg(argv[a]);
      if (a + 1 >= argc) throw std::runtime_error("expected --steps N, --dt T, or -t FILE");
      const std::string value(argv[++a]);
      std::size_t used = 0;
      if (arg == "--steps") {
        const auto parsed = std::stoul(value, &used);
        if (used != value.size() || parsed < 1 || parsed > 1000000) throw std::runtime_error("invalid step count");
        steps = static_cast<unsigned>(parsed);
      } // END IF: parse step count
else if (arg == "--dt") {
        dt = std::stod(value, &used);
        if (used != value.size()) throw std::runtime_error("invalid timestep");
      } // END ELSE IF: parse timestep
else if (arg == "-t") filename = value;
      else throw std::runtime_error("unknown command-line option");
    } // END LOOP: parse runtime arguments
    if (!(dt > 0.0) || !std::isfinite(dt) || !std::isfinite(dt*steps)) throw std::runtime_error("invalid timestep");
    m_uiMaxDepth = 8;
    _InitializeHcurve(m_uiDim);
    std::vector<ot::TreeNode> octree;
    std::function<double(double,double,double)> refine = [](double x,double y,double z) {
      return std::exp(-(x*x+y*y+z*z)/0.5);
    }; // END LAMBDA: choose initial octree refinement
    // Build once; no remeshing during this qualified fixed-mesh run.
    const unsigned order = 2 * $NAMESPACE::generated::REQUIRED_PADDING;
    function2Octree(refine, octree, 5, 1e-3, order, MPI_COMM_WORLD);
    std::unique_ptr<ot::Mesh> mesh(ot::createMesh(octree.data(), octree.size(), order, MPI_COMM_WORLD, 0, ot::SM_TYPE::FDM));
    if (!mesh) throw std::runtime_error("mesh construction failed");
    const Point minimum(-1.0,-2.0,-4.0), maximum(3.0,2.0,4.0);
    mesh->setDomainBounds(minimum, maximum);
    {
      $NAMESPACE::Ctx context(mesh.get(), minimum, maximum, dt);
      if (!filename.empty()) {
        std::string contents;
        if (rank == 0) {
          std::ifstream input(filename);
          if (!input) throw std::runtime_error("cannot open parameter file");
          contents.assign(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
          if (contents.size() > 1048576) throw std::runtime_error("parameter file exceeds 1 MiB");
        } // END IF: read file on root
        int length = static_cast<int>(contents.size());
        MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
        contents.resize(length);
        MPI_Bcast(contents.data(), length, MPI_CHAR, 0, MPI_COMM_WORLD);
        std::istringstream input(contents);
        const auto document = toml::parse(input, filename);
        for (const auto& item : document.as_table())
          if (item.first != "params" && item.first != "$STEM") throw std::runtime_error("unknown parameter table");
        if (document.contains("$STEM")) {
          const auto& app = document.at("$STEM");
          if (app.as_table().size() != 1 || !app.contains("profile")) throw std::runtime_error("unknown solver table");
          const auto& profile = app.at("profile");
          for (const auto& item : profile.as_table())
            if (item.first != "name" && item.first != "fd_order" && item.first != "required_padding" && item.first != "ko_enabled") throw std::runtime_error("unknown profile key");
          if (toml::find<unsigned>(profile,"fd_order") != $NAMESPACE::generated::FD_ORDER ||
              toml::find<unsigned>(profile,"required_padding") != $NAMESPACE::generated::REQUIRED_PADDING ||
              toml::find<bool>(profile,"ko_enabled") != $NAMESPACE::generated::KO_ENABLED)
            throw std::runtime_error("parameter profile does not match generated kernels");
        } // END IF: validate generated profile
        if (document.contains("params")) {
          const auto& table = document.at("params");
          auto& params = context.params;
          for (const auto& item : table.as_table()) {
$PARAMETER_BINDINGS
            throw std::runtime_error("unknown runtime parameter: " + item.first);
          } // END LOOP: bind registered runtime parameters
        } // END IF: apply parameter table
      } // END IF: read requested TOML file
      if (!$VALIDATE(context.params)) throw std::runtime_error("invalid runtime parameters");
      if (rank == 0) $PRINT_EFFECTIVE(context.params);
      ts::ETS<DendroScalar, $NAMESPACE::Ctx> stepper(&context);
      stepper.set_ets_coefficients(ts::ETSType::RK4);
      stepper.init();
      for (unsigned step = 0; step < steps; ++step) stepper.evolve();
      const double rhs = context.max_rhs(), constraints = context.max_constraints(), drift = context.max_drift();
      double residual = 0.0;
      MPI_Allreduce(&context.projection_residual, &residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      double local_spacing = std::numeric_limits<double>::max();
      if (mesh->isActive())
        for (const auto& block : mesh->getLocalBlockList()) {
          const auto geometry = $NAMESPACE::block_geometry(*mesh, block, minimum, maximum);
          for (double spacing : geometry.dx) local_spacing = std::min(local_spacing, spacing);
        } // END LOOP: find smallest physical spacing
      double spacing = 0.0;
      MPI_Allreduce(&local_spacing, &spacing, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      // Coarse/fine interpolation introduces roundoff even for constant data.
      // Second derivatives amplify it by h^-2; bound the normalized residual
      // by 256 double-precision ulps. The state-drift bound stays independent.
      const double derivative_tolerance = 256*std::numeric_limits<double>::epsilon()/(spacing*spacing);
      int local_ok = (!mesh->isActive() || context.projection_passes == 1 + 5ULL*steps) &&
        stepper.curr_step() == steps && std::abs(stepper.curr_time() - dt*steps) <= 1e-11*std::max(1.0, dt*steps) &&
        rhs <= derivative_tolerance && constraints <= derivative_tolerance && drift <= 1e-11 && residual <= 1e-13;
      int global_ok = 0;
      MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      int active = mesh->isActive(), active_ranks = 0;
      MPI_Allreduce(&active, &active_ranks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (rank == 0) std::printf("REAL_MINKOWSKI %s active_ranks=%d steps=%u time=%.17g rhs=%.17g constraints=%.17g drift=%.17g projection=%.17g hmin=%.17g derivative_tolerance=%.17g\n", global_ok ? "PASS" : "FAIL", active_ranks, steps, stepper.curr_time(), rhs, constraints, drift, residual, spacing, derivative_tolerance);
      if (!global_ok) throw std::runtime_error("fixed-mesh Minkowski check failed");
    } // END BLOCK: own context before mesh destruction
  } // END TRY: run fixed mesh qualification
catch (const std::exception& error) {
    std::fprintf(stderr, "rank %d: %s\n", rank, error.what());
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  } // END CATCH: abort all parent ranks
  MPI_Finalize();
  return 0;
} // END FUNCTION: run real host solver
"""


def output_main_cpp(
    solver_stem: str,
    solver_namespace: str,
    exec_or_library_name: str,
) -> str:
    """
    Emit standalone and checked real-host entry points.

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
    text = (
        "#if defined(NRPY_DENDRO_STANDALONE_HOST)\n"
        + _MAIN.replace("$EXEC", exec_or_library_name)
        + "\n#else\n"
        + _REAL_MAIN.replace("$PARAMETER_BINDINGS", output_toml_bindings())
        + "\n#endif\n"
    )
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
