# nrpy/examples/photon_single_geodesic_integrator_analytical.py
r"""
Generate a standalone C project for integrating a single photon geodesic.

The generated project evolves one massless test particle in an analytic
spacetime using the split RKF45 photon pipeline. It writes trajectory samples
and reports normalization and conserved-quantity diagnostics while preserving
the Structure of Arrays layout expected by the shared geodesic kernels.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import argparse
import os
import shutil
from typing import List, Optional

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
    Analytic_Spacetimes,
)
from nrpy.equations.general_relativity.geodesics.geodesics import Geodesic_Equations
from nrpy.infrastructures.BHaH import (
    BHaH_defines_h,
)
from nrpy.infrastructures.BHaH import CodeParameters as CPs
from nrpy.infrastructures.BHaH import Makefile_helpers as Makefile
from nrpy.infrastructures.BHaH import (
    cmdline_input_and_parfiles,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics import (
    connections,
    conserved_quantities,
    g4DD_metric,
    normalization_constraint,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
    calculate_ode_rhs_kernel,
    interpolation_kernel,
    main_single,
    p0_reverse_kernel,
    rkf45_finalize_and_control_kernel,
    rkf45_stage_update,
    single_integrator_analytical,
)


def _build_parser() -> argparse.ArgumentParser:
    """
    Build the analytical photon single-geodesic generator argument parser.

    :return: Configured command-line argument parser.
    """
    parser = argparse.ArgumentParser(
        description="Generate a standalone analytical photon geodesic project."
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
    Parse analytical photon single-geodesic generator arguments.

    :param argv: Optional argument list; defaults to ``sys.argv`` when omitted.
    :return: Parsed command-line arguments.
    """
    return _build_parser().parse_args(argv)


if __name__ == "__main__":
    args = _parse_args()

    # Step 1: Select the BHaH OpenMP backend for this example.
    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")

    # Step 2: Define project-level constants and output paths.
    project_name = "photon_single_geodesic_integrator_analytical"
    project_dir = os.path.abspath(os.path.join("project", project_name))
    parfile_path = os.path.join(project_dir, f"{project_name}.par")

    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"

    # Step 3: Recreate the generated project directory.
    shutil.rmtree(project_dir, ignore_errors=True)
    os.makedirs(project_dir, exist_ok=True)

    # Step 4: Acquire symbolic data for the selected spacetime.
    print(f"Acquiring symbolic data for {GEO_KEY}...")
    metric_data = Analytic_Spacetimes[SPACETIME]
    geodesic_data = Geodesic_Equations[GEO_KEY]

    # Step 5: Register the shared physics and diagnostics kernels.
    print("Registering Split-Pipeline Physics Kernels...")
    g4DD_metric.g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)
    connections.connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)
    conserved_quantities.conserved_quantities(SPACETIME, PARTICLE)
    normalization_constraint.normalization_constraint(
        geodesic_data.norm_constraint_expr, PARTICLE
    )

    if geodesic_data.p0_photon is None:
        raise ValueError(f"p0_photon is None for {GEO_KEY}")
    p0_reverse_kernel.p0_reverse_kernel(geodesic_data.p0_photon)
    interpolation_kernel.interpolation_kernel(SPACETIME)
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(
        geodesic_data.geodesic_eom_rhs_photon_christoffel(),
        geodesic_data.xx,
        use_metric_derivative_rhs=False,
        normalized_eom=False,
    )
    rkf45_stage_update.rkf45_stage_update()
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel(
        normalized_eom=False
    )

    # Step 5.a: Register the single-ray C main function.
    single_integrator_analytical.single_integrator_analytical(
        SPACETIME, PARTICLE, normalized_eom=False
    )
    main_single.main_single("single_integrator_analytical")

    # Step 5.b: Remove helper registrations emitted only through other kernels.
    # The generated main() calls the metric and connection logic through shared
    # kernels, so emitting standalone wrappers here would add redundant source
    # files and prototypes.
    for internal_func in [f"g4DD_metric_{SPACETIME}", f"connections_{SPACETIME}"]:
        cfc.CFunction_dict.pop(internal_func, None)

    # Step 6: Set relevant CodeParameter defaults.
    par.adjust_CodeParam_default("rkf45_absolute_error_tolerance", 1e-17)
    par.adjust_CodeParam_default("rkf45_error_tolerance", 1e-17)
    par.adjust_CodeParam_default("rkf45_h_max", 10.0)
    par.adjust_CodeParam_default("rkf45_h_min", 1e-20)
    par.adjust_CodeParam_default("rkf45_max_retries", 15)
    for name, value in zip(
        ["initial_x", "initial_y", "initial_z"], args.initial_position
    ):
        par.adjust_CodeParam_default(name, value)
    for name, value in zip(
        ["initial_p_x", "initial_p_y", "initial_p_z"], args.initial_momentum
    ):
        par.adjust_CodeParam_default(name, value)
    par.adjust_CodeParam_default("initial_t", args.initial_t)
    par.adjust_CodeParam_default("initial_integration_param", 0.0)
    par.adjust_CodeParam_default("initial_eulerian_distance", 0.0)
    par.adjust_CodeParam_default("initial_h", 0.1)
    par.adjust_CodeParam_default("r_escape", args.r_escape)
    par.adjust_CodeParam_default("evolution_measure_max", 1000.0)
    par.adjust_CodeParam_default("max_steps", 200000)

    # Step 7: Generate headers, default parameters, and the Makefile.
    print("Generating header files and Makefile...")
    CPs.write_CodeParameters_h_files(set_commondata_only=True, project_dir=project_dir)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()

    cmdline_input_and_parfiles.generate_default_parfile(
        project_dir=project_dir, project_name=project_name
    )
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
        project_name=project_name
    )

    # BHAH_MALLOC expects an assignable lvalue, so keep direct forwarding in
    # the compatibility wrappers below.
    cpu_macros = {
        "ReadCUDA(ptr)": "#define ReadCUDA(ptr) (*(ptr))\n",
        "WriteCUDA(ptr, val)": "#define WriteCUDA(ptr, val) (*(ptr) = (val))\n",
        "BHAH_MALLOC_DEVICE(a, sz)": (
            "#define BHAH_MALLOC_DEVICE(a, sz) BHAH_MALLOC(a, sz)\n"
        ),
        "BHAH_FREE_DEVICE(a)": "#define BHAH_FREE_DEVICE(a) BHAH_FREE(a)\n",
        "MulCUDA(a, b)": "#define MulCUDA(a, b) ((a) * (b))\n",
        "DivCUDA(a, b)": "#define DivCUDA(a, b) ((a) / (b))\n",
        "AddCUDA(a, b)": "#define AddCUDA(a, b) ((a) + (b))\n",
        "FusedMulAddCUDA(a, b, c)": (
            "#define FusedMulAddCUDA(a, b, c) ((a) * (b) + (c))\n"
        ),
        "AbsCUDA(val)": "#define AbsCUDA(val) fabs(val)\n",
        "SqrtCUDA(val)": "#define SqrtCUDA(val) sqrt(val)\n",
        "PowCUDA(a, b)": "#define PowCUDA(a, b) pow(a, b)\n",
        "BHAH_HD_FUNC": "#define BHAH_HD_FUNC\n",
        "BHAH_HD_INLINE": "#define BHAH_HD_INLINE static inline\n",
        "BHAH_WARP_ATOMIC_ADD(ptr, val)": (
            "#define BHAH_WARP_ATOMIC_ADD(ptr, val) "
            '_Pragma("omp atomic") *(ptr) += (val)\n'
        ),
        "GLOBAL_COMMONDATA_EXTERN": (
            "// Commondata is passed by reference in this OpenMP example.\n"
        ),
        "BHAH_DEVICE_SYNC()": "#define BHAH_DEVICE_SYNC() do {} while(0)\n",
    }

    BHaH_defines_h.output_BHaH_defines_h(
        project_dir=project_dir,
        enable_rfm_precompute=False,
        supplemental_defines_dict=cpu_macros,
    )

    addl_cflags = [
        "-fopenmp",
        "-O3",
        "-DDEBUG",
        "-Wno-stringop-truncation",
    ]
    addl_libs = ["-lm"]
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

    # Step 8: Copy the trajectory-visualization helper and print usage.
    vis_dir = os.path.join("nrpy", "examples", "geodesic_visualizations")
    vis_script_src = os.path.join(vis_dir, "visualize_trajectory.py")

    if os.path.exists(vis_script_src):
        shutil.copy(vis_script_src, project_dir)
    else:
        print(
            f"Warning: Visualization script not found at {vis_script_src}. "
            "Please ensure it exists."
        )

    print(
        f"Finished! Now go into {project_dir} and type `make` to build, "
        f"then ./{project_name} to run."
    )
    print(f"    Parameter file can be found at {parfile_path}\n")
    print(
        "    To generate the trajectory plot after running the C executable, "
        "ensure you have the required Python packages:"
    )
    print("    pip install matplotlib numpy\n")
    print(
        "    Then, execute the visualization script directly from the "
        "project directory:"
    )
    print("    python3 visualize_trajectory.py --particle_type Photon\n")
