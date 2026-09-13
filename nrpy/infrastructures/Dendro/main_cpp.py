# nrpy/infrastructures/Dendro/main_cpp.py
"""
Emit generic standalone and real-host Dendro process shells.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from nrpy.infrastructures.Dendro.CodeParameters import output_toml_bindings
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.solver_context import substitute_solver_identifiers

BANNER = generated_file_banner()

_STANDALONE_MAIN = r"""// Standalone-host entry point.
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
        if (rank == 0)
          std::fprintf(stderr, "too many -r selections (max %d)\n", MAX_SELECTED);
        MPI_Finalize();
        return 2;
      }  // END ELSE: selection table full
    }  // END ELSE IF: select requested result
    else {
      if (rank == 0)
        std::fprintf(stderr,
                     "usage: $EXEC [-b n_blocks] [-n extent] [-d dx] "
                     "[-t parfile] [-r exact_name]...\n");
      MPI_Finalize();
      return 2;
    }  // END ELSE: reject unknown argument
  }  // END LOOP: parse standalone arguments

  $NAMESPACE::Ctx ctx;
  if (ctx.initialize_mesh(n_blocks, extent, dx, rank, parfile) ||
      ctx.startup_checks(rank)) {
    MPI_Finalize();
    return 1;
  }  // END IF: generic context setup failed
$STANDALONE_APPLICATION_INITIALIZATION
  if (n_selected > 0 && ctx.select_variables(selected, n_selected) != 0) {
    MPI_Finalize();
    return 1;
  }  // END IF: selected name is unknown
$STANDALONE_APPLICATION_BEFORE_STEPS
  const int nsteps = $STANDALONE_STEP_COUNT;
  const double dt = $STANDALONE_TIMESTEP;
  for (int step = 0; step < nsteps; ++step) {
    if (ctx.euler_step(dt)) {
      MPI_Finalize();
      return 1;
    }  // END IF: standalone timestep failed
$STANDALONE_APPLICATION_AFTER_STEP
  }  // END LOOP: evolve standalone host
$STANDALONE_APPLICATION_FINAL_CHECKS
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
#include <limits>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    unsigned steps = 100;
    double dt = 0.001;
    std::string filename;
    for (int argument = 1; argument < argc; ++argument) {
      const std::string option(argv[argument]);
      if (argument + 1 >= argc)
        throw std::runtime_error("expected --steps N, --dt T, or -t FILE");
      const std::string value(argv[++argument]);
      std::size_t used = 0;
      if (option == "--steps") {
        const auto parsed = std::stoul(value, &used);
        if (used != value.size() || parsed < 1 || parsed > 1000000)
          throw std::runtime_error("invalid step count");
        steps = static_cast<unsigned>(parsed);
      } else if (option == "--dt") {
        dt = std::stod(value, &used);
        if (used != value.size()) throw std::runtime_error("invalid timestep");
      } else if (option == "-t") {
        filename = value;
      } else {
        throw std::runtime_error("unknown command-line option");
      }  // END ELSE: reject unknown real-host argument
    }  // END LOOP: parse real-host arguments
    if (!(dt > 0.0) || !std::isfinite(dt) || !std::isfinite(dt * steps))
      throw std::runtime_error("invalid timestep");
    m_uiMaxDepth = 8;
    _InitializeHcurve(m_uiDim);
    std::vector<ot::TreeNode> octree;
    std::function<double(double, double, double)> refine =
        [](double x, double y, double z) {
          return std::exp(-(x*x + y*y + z*z) / 0.5);
        };  // END LAMBDA: choose initial octree refinement
    const unsigned order = 2 * $NAMESPACE::generated::REQUIRED_PADDING;
    function2Octree(refine, octree, 5, 1e-3, order, MPI_COMM_WORLD);
    std::unique_ptr<ot::Mesh> mesh(ot::createMesh(
        octree.data(), octree.size(), order, MPI_COMM_WORLD, 0,
        ot::SM_TYPE::FDM));
    if (!mesh) throw std::runtime_error("mesh construction failed");
    const Point minimum(-1.0, -2.0, -4.0), maximum(3.0, 2.0, 4.0);
    mesh->setDomainBounds(minimum, maximum);
    {
      $NAMESPACE::Ctx context(mesh.get(), minimum, maximum, dt);
      if (!filename.empty()) {
        std::string contents;
        if (rank == 0) {
          std::ifstream input(filename);
          if (!input) throw std::runtime_error("cannot open parameter file");
          contents.assign(std::istreambuf_iterator<char>(input),
                          std::istreambuf_iterator<char>());
          if (contents.size() > 1048576)
            throw std::runtime_error("parameter file exceeds 1 MiB");
        }  // END IF: root reads parameter file
        int length = static_cast<int>(contents.size());
        MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
        contents.resize(length);
        MPI_Bcast(contents.data(), length, MPI_CHAR, 0, MPI_COMM_WORLD);
        std::istringstream input(contents);
        const auto document = toml::parse(input, filename);
        for (const auto& item : document.as_table())
          if (item.first != "params" && item.first != "$STEM")
            throw std::runtime_error("unknown parameter table");
        if (document.contains("$STEM")) {
          const auto& app = document.at("$STEM");
          if (app.as_table().size() != 1 || !app.contains("profile"))
            throw std::runtime_error("unknown solver table");
          const auto& profile = app.at("profile");
          for (const auto& item : profile.as_table())
            if (item.first != "name" && item.first != "fd_order" &&
                item.first != "required_padding" && item.first != "ko_enabled")
              throw std::runtime_error("unknown profile key");
          if (toml::find<unsigned>(profile, "fd_order") != $NAMESPACE::generated::FD_ORDER ||
              toml::find<unsigned>(profile, "required_padding") != $NAMESPACE::generated::REQUIRED_PADDING ||
              toml::find<bool>(profile, "ko_enabled") != $NAMESPACE::generated::KO_ENABLED)
            throw std::runtime_error("parameter profile does not match generated kernels");
        }  // END IF: validate generated profile
        if (document.contains("params")) {
          const auto& table = document.at("params");
          auto& params = context.params;
          for (const auto& item : table.as_table()) {
$PARAMETER_BINDINGS
            throw std::runtime_error("unknown runtime parameter: " + item.first);
          }  // END LOOP: bind registered runtime parameters
        }  // END IF: apply runtime parameters
      }  // END IF: read requested TOML file
      if (!$VALIDATE(context.params))
        throw std::runtime_error("invalid runtime parameters");
      if (rank == 0) $PRINT_EFFECTIVE(context.params);
      ts::ETS<DendroScalar, $NAMESPACE::Ctx> stepper(&context);
      stepper.set_ets_coefficients(ts::ETSType::RK4);
      stepper.init();
      for (unsigned step = 0; step < steps; ++step) stepper.evolve();
$REAL_APPLICATION_FINAL_CHECKS
    }  // END BLOCK: destroy context before borrowed mesh
  }  // END TRY: run real host solver
  catch (const std::exception& error) {
    std::fprintf(stderr, "rank %d: %s\n", rank, error.what());
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }  // END CATCH: abort all parent ranks
  MPI_Finalize();
  return 0;
}  // END FUNCTION: main
"""


def output_main_cpp(
    solver_stem: str,
    solver_namespace: str,
    exec_or_library_name: str,
    standalone_application_initialization: str,
    standalone_application_before_steps: str,
    standalone_application_after_step: str,
    standalone_application_final_checks: str,
    standalone_step_count: str,
    standalone_timestep: str,
    real_application_final_checks: str,
) -> str:
    """
    Emit process shells using explicit application lifecycle statements.

    :param solver_stem: Lowercase formulation stem used in emitted names.
    :param solver_namespace: Namespace containing the generated solver.
    :param exec_or_library_name: Executable name used in standalone help.
    :param standalone_application_initialization: Application initialization
        statements after generic context setup.
    :param standalone_application_before_steps: Application checks and setup
        before the generic standalone step loop.
    :param standalone_application_after_step: Application policy after each
        accepted standalone step.
    :param standalone_application_final_checks: Application acceptance checks
        after standalone stepping.
    :param standalone_step_count: C++ expression for the standalone step count.
    :param standalone_timestep: C++ expression for the standalone timestep.
    :param real_application_final_checks: Application acceptance checks after
        real-host ETS stepping.
    :return: Complete generated C++ source.
    :raises ValueError: If a known application insertion remains unresolved.

    Doctests:
    >>> wave = output_main_cpp(
    ...     "wave", "wave", "waveSolver",
    ...     "  if (ctx.initialize_scalar_vector()) return 1;",
    ...     "  const double wave_norm = ctx.max_wave_rhs();",
    ...     "    ctx.apply_wave_boundary();",
    ...     "  if (!std::isfinite(wave_norm)) return 1;",
    ...     "8", "0.25 * dx",
    ...     "      if (!std::isfinite(context.max_wave_rhs())) return 1;",
    ... )
    >>> all(token in wave for token in (
    ...     "initialize_scalar_vector", "apply_wave_boundary", "max_wave_rhs"
    ... ))
    True
    >>> all(token not in wave for token in (
    ...     "minkowski", "detgtrazero", "max_constraints"
    ... ))
    True
    """
    text = (
        "#if defined(NRPY_DENDRO_STANDALONE_HOST)\n"
        + _STANDALONE_MAIN.replace("$EXEC", exec_or_library_name)
        + "\n#else\n"
        + _REAL_MAIN.replace("$PARAMETER_BINDINGS", output_toml_bindings())
        + "\n#endif\n"
    )
    replacements = (
        (
            "$STANDALONE_APPLICATION_INITIALIZATION",
            standalone_application_initialization,
        ),
        ("$STANDALONE_APPLICATION_BEFORE_STEPS", standalone_application_before_steps),
        ("$STANDALONE_APPLICATION_AFTER_STEP", standalone_application_after_step),
        ("$STANDALONE_APPLICATION_FINAL_CHECKS", standalone_application_final_checks),
        ("$STANDALONE_STEP_COUNT", standalone_step_count),
        ("$STANDALONE_TIMESTEP", standalone_timestep),
        ("$REAL_APPLICATION_FINAL_CHECKS", real_application_final_checks),
    )
    for token, value in replacements:
        text = text.replace(token, value)
    unresolved = tuple(token for token, _value in replacements if token in text)
    if unresolved:
        raise ValueError(
            f"Application lifecycle insertions were not resolved: {unresolved}"
        )
    return BANNER + substitute_solver_identifiers(text, solver_stem, solver_namespace)


if __name__ == "__main__":
    import doctest

    raise SystemExit(doctest.testmod().failed)
