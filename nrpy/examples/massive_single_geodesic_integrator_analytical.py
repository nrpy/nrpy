# nrpy/examples/massive_single_geodesic_integrator_analytical.py
r"""
Defines a C project for integrating massive geodesics in curved spacetime.

The script provides a standalone C application that evolves the trajectory
of a massive test particle. It coordinates the registration of spacetime-specific
physics kernels and links them with the GNU Scientific Library for time integration.
Structure of Arrays (SoA) layouts provide universal memory compatibility for shared
physics kernels. Struct mapping establishes the layout required by downstream conserved
quantities kernels. The initial metric is evaluated at the starting position so that
u^t can be computed from the massive-particle normalization condition. The final
normalization residual is reported at the last trajectory point.

The simulation solves the geodesic equation for a particle with non-zero mass:
    $d(u^\mu)/d(\tau) = -\Gamma^\mu_{\alpha \beta} u^\alpha u^\beta$
subject to the normalization constraint $u^\mu u_\mu = -1$.

Numerical fidelity is validated by monitoring constants of motion
associated with the spacetime's symmetries, including Killing vectors and tensors.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import argparse
import os
import shutil
from typing import List, Optional

import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.params as par
from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
    Analytic_Spacetimes,
)
from nrpy.equations.general_relativity.geodesics.geodesics import Geodesic_Equations
from nrpy.infrastructures.BHaH import cmdline_input_and_parfiles
from nrpy.infrastructures.BHaH.general_relativity.geodesics.connections import (
    connections,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.conserved_quantities import (
    conserved_quantities,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.g4DD_metric import (
    g4DD_metric,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive.calculate_ode_rhs_massive import (
    calculate_ode_rhs_massive,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive.main_single import (
    main_single,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive.ode_gsl_wrapper_massive import (
    ode_gsl_wrapper_massive,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive.single_integrator_analytical import (
    single_integrator_analytical,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive.u0_massive import (
    u0_massive,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.normalization_constraint import (
    normalization_constraint,
)


def _build_parser() -> argparse.ArgumentParser:
    """
    Build the massive single-geodesic generator argument parser.

    :return: Configured command-line argument parser.
    """
    parser = argparse.ArgumentParser(
        description="Generate a standalone analytical massive geodesic project."
    )
    parser.add_argument(
        "--initial-time",
        dest="initial_t",
        type=float,
        default=0.0,
        help="Initial coordinate time.",
    )
    parser.add_argument(
        "--initial-position",
        dest="initial_position",
        nargs=3,
        type=float,
        default=[4.0123, 0.0, 0.0],
        metavar=("X", "Y", "Z"),
        help="Initial Cartesian position.",
    )
    parser.add_argument(
        "--initial-momentum",
        dest="initial_momentum",
        nargs=3,
        type=float,
        default=[-0.5641, 0.0, 0.0],
        metavar=("PX", "PY", "PZ"),
        help="Initial spatial momentum components.",
    )
    parser.add_argument(
        "--r-escape",
        dest="r_escape",
        type=float,
        default=150.0,
        help="Escape radius used to terminate the trajectory.",
    )
    return parser


def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """
    Parse massive single-geodesic generator arguments.

    :param argv: Optional argument list; defaults to ``sys.argv`` when omitted.
    :return: Parsed command-line arguments.
    """
    return _build_parser().parse_args(argv)


if __name__ == "__main__":
    args = _parse_args()

    # Step P1: Set code-generation parameters and register the geodesic kernels.
    enable_parallel_codegen = True
    if enable_parallel_codegen:
        pcg.do_parallel_codegen()

    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")

    project_name = "massive_single_geodesic_integrator_analytical"
    project_dir = os.path.join("project", project_name)

    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "massive"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"

    if os.path.exists(project_dir):
        shutil.rmtree(project_dir)
    os.makedirs(project_dir, exist_ok=True)

    print(f"Acquiring symbolic data for {GEO_KEY}...")
    metric_data = Analytic_Spacetimes[SPACETIME]
    geodesic_data = Geodesic_Equations[GEO_KEY]

    print("Registering Physics Kernels...")
    g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)
    connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)
    calculate_ode_rhs_massive(geodesic_data.geodesic_rhs, metric_data.xx)

    if geodesic_data.u0_massive is None:
        raise ValueError(f"u0_massive is None for {GEO_KEY}")
    u0_massive(geodesic_data.u0_massive)

    conserved_quantities(SPACETIME, PARTICLE)
    normalization_constraint(geodesic_data.norm_constraint_expr, PARTICLE)
    ode_gsl_wrapper_massive(SPACETIME)

    single_integrator_analytical(SPACETIME, PARTICLE)
    main_single("single_integrator_analytical")

    # Step P2: Apply command-line initial-state and escape-radius overrides.
    for name, value in zip(
        ["initial_x", "initial_y", "initial_z"], args.initial_position
    ):
        par.adjust_CodeParam_default(name, value)
    for name, value in zip(
        ["initial_p_x", "initial_p_y", "initial_p_z"], args.initial_momentum
    ):
        par.adjust_CodeParam_default(name, value)
    par.adjust_CodeParam_default("initial_t", args.initial_t)
    par.adjust_CodeParam_default("r_escape", args.r_escape)

    print("Generating header files and Makefile...")

    CPs.write_CodeParameters_h_files(set_commondata_only=True, project_dir=project_dir)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()

    cmdline_input_and_parfiles.generate_default_parfile(
        project_dir=project_dir, project_name=project_name
    )
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
        project_name=project_name
    )

    macro_defs = r"""
    // Provide standalone indexing helpers used by the shared geodesic kernels.
    #ifndef IDX_LOCAL
    #define IDX_LOCAL(component, batch_id, batch_size) ((component) * (batch_size) + (batch_id))
    #endif

    #ifndef IDX_GLOBAL
    #define IDX_GLOBAL(component, ray_id, num_rays) ((component) * (num_rays) + (ray_id))
    #endif

    // Single-particle GSL driver: set BUNDLE_CAPACITY to 1 for shared memory-layout macros.
    #ifndef BUNDLE_CAPACITY
    #define BUNDLE_CAPACITY 1
    #endif
    """
    Bdefines_h.register_BHaH_defines("gpu_batch_macros", macro_defs)

    cpu_macros = {
        "BHAH_MALLOC_DEVICE(a, sz)": "#define BHAH_MALLOC_DEVICE(a, sz) BHAH_MALLOC((a), (sz))",
        "BHAH_FREE_DEVICE(a)": "#define BHAH_FREE_DEVICE(a) BHAH_FREE((a))",
        "BHAH_HD_INLINE": "#define BHAH_HD_INLINE",
    }

    additional_includes = [
        "gsl/gsl_vector.h",
        "gsl/gsl_matrix.h",
        "gsl/gsl_odeiv2.h",
        "gsl/gsl_errno.h",
        "gsl/gsl_math.h",
    ]

    Bdefines_h.output_BHaH_defines_h(
        project_dir=project_dir,
        additional_includes=additional_includes,
        enable_rfm_precompute=False,
        supplemental_defines_dict=cpu_macros,
    )

    addl_cflags = ["$(shell gsl-config --cflags)", "-Wno-stringop-truncation"]
    addl_libs = ["$(shell gsl-config --libs)"]

    Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=project_name,
        compiler_opt_option="fast",
        addl_CFLAGS=addl_cflags,
        addl_libraries=addl_libs,
        CC="gcc",
        src_code_file_ext="c",
    )

    # Step P2: Copy the trajectory visualizer into the generated project.
    vis_dir = os.path.join("nrpy", "examples", "geodesic_visualizations")
    vis_script_src = os.path.join(vis_dir, "visualize_trajectory.py")

    if os.path.exists(vis_script_src):
        shutil.copy(vis_script_src, project_dir)
    else:
        print(
            f"Warning: Visualization script not found at {vis_script_src}; trajectory plotting script was not copied."
        )

    print(
        f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
    )
    print(f"    Parameter file can be found in {project_name}.par\n")
    print(
        "    To generate the trajectory plot after running the C executable, ensure you have the required Python packages:"
    )
    print("    pip install matplotlib numpy\n")
    print(
        "    Then, execute the visualization script directly from the project directory:"
    )
    print("    python3 visualize_trajectory.py --particle_type Massive\n")
