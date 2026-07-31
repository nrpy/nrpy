# nrpy/examples/photon_single_geodesic_integrator_analytical.py
r"""
Generate a standalone C project for integrating a single photon geodesic.

The generated project evolves one massless test particle in an analytic
spacetime using the split RKF45 photon pipeline. Its one-ray initialization
uses the observer metric tetrad and the shared one-tile sampling contract. It
writes trajectory samples and reports normalization and conserved-quantity
diagnostics while preserving the Structure of Arrays layout expected by the
shared geodesic kernels.

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
    event_detection_manager_kernel,
    interpolation_kernel,
    main_single,
    normal_observer_log_energy,
    rkf45_finalize_and_control_kernel,
    rkf45_stage_update,
    set_initial_conditions_kernel,
    single_integrator_analytical,
)


def _build_parser() -> argparse.ArgumentParser:
    """
    Build the analytical photon single-geodesic generator argument parser.

    :return: Configured command-line argument parser.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Generate a standalone analytical photon geodesic project using "
            "the shared observer-tetrad initializer."
        )
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="Parent directory for the generated C project.",
    )
    parser.add_argument(
        "--observer-position",
        nargs=3,
        type=float,
        required=True,
        metavar=("X", "Y", "Z"),
        help="Observer coordinate position.",
    )
    parser.add_argument(
        "--observer-look-forward",
        nargs=3,
        type=float,
        required=True,
        metavar=("X", "Y", "Z"),
        help="Observer look-forward direction seed.",
    )
    parser.add_argument(
        "--observer-up",
        nargs=3,
        type=float,
        required=True,
        metavar=("UX", "UY", "UZ"),
        help="Observer up/roll direction seed.",
    )
    parser.add_argument(
        "--observer-fov",
        nargs=2,
        type=float,
        default=[1.0, 1.0],
        metavar=("ALPHA_W", "ALPHA_H"),
        help="Observer width and height fields of view, in radians.",
    )
    parser.add_argument(
        "--t-start",
        type=float,
        default=0.0,
        help="Initial observer-event coordinate time; default: 0.0.",
    )
    parser.add_argument(
        "--escape-radius",
        type=float,
        required=True,
        metavar="R",
        help="Coordinate-radius termination threshold.",
    )
    parser.add_argument(
        "--terminal-plane-center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="Terminal-plane center; requires all terminal-plane options.",
    )
    parser.add_argument(
        "--terminal-plane-normal",
        nargs=3,
        type=float,
        metavar=("NX", "NY", "NZ"),
        help="Terminal-plane normal; requires all terminal-plane options.",
    )
    parser.add_argument(
        "--terminal-plane-up",
        nargs=3,
        type=float,
        metavar=("UX", "UY", "UZ"),
        help="Terminal-plane up direction; requires all terminal-plane options.",
    )
    parser.add_argument(
        "--terminal-plane-radius",
        nargs=2,
        type=float,
        metavar=("MIN_RADIUS", "MAX_RADIUS"),
        help="Terminal-plane accepted coordinate-radius range.",
    )
    parser.add_argument(
        "--non-terminal-plane-center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="Nonterminal-plane center; requires all nonterminal-plane options.",
    )
    parser.add_argument(
        "--non-terminal-plane-normal",
        nargs=3,
        type=float,
        metavar=("NX", "NY", "NZ"),
        help="Nonterminal-plane normal; requires all nonterminal-plane options.",
    )
    parser.add_argument(
        "--non-terminal-plane-up",
        nargs=3,
        type=float,
        metavar=("UX", "UY", "UZ"),
        help="Nonterminal-plane up direction; requires all nonterminal-plane options.",
    )
    parser.add_argument(
        "--initial-step",
        type=float,
        metavar="H",
        help="Initial RKF45 step magnitude; default: current analytical value.",
    )
    parser.add_argument(
        "--rkf45-tolerances",
        nargs=2,
        type=float,
        metavar=("ABSOLUTE", "RELATIVE"),
        help="RKF45 absolute and relative error tolerances.",
    )
    parser.add_argument(
        "--rkf45-step-range",
        nargs=2,
        type=float,
        metavar=("H_MIN", "H_MAX"),
        help="RKF45 minimum and maximum step magnitudes.",
    )
    parser.add_argument(
        "--rkf45-max-retries",
        type=int,
        metavar="N",
        help="Maximum RKF45 rejected-step retries; default: current value.",
    )
    parser.add_argument(
        "--evolution-measure-max",
        type=float,
        metavar="VALUE",
        help="Upper limit on the common normal-observer log-energy measure ln(abs(alpha*p^0)).",
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
    import sys

    if len(sys.argv) == 1:
        _build_parser().print_help()
        sys.exit(0)
    args = _parse_args()

    if any(value <= 0.0 for value in args.observer_fov):
        raise ValueError("--observer-fov values must be positive.")
    if args.escape_radius <= 0.0:
        raise ValueError("--escape-radius must be positive.")
    terminal_plane_values = (
        args.terminal_plane_center,
        args.terminal_plane_normal,
        args.terminal_plane_up,
        args.terminal_plane_radius,
    )
    if any(value is not None for value in terminal_plane_values) and not all(
        value is not None for value in terminal_plane_values
    ):
        raise ValueError(
            "Terminal-plane options must be supplied as one complete group."
        )
    non_terminal_plane_values = (
        args.non_terminal_plane_center,
        args.non_terminal_plane_normal,
        args.non_terminal_plane_up,
    )
    if any(value is not None for value in non_terminal_plane_values) and not all(
        value is not None for value in non_terminal_plane_values
    ):
        raise ValueError(
            "Nonterminal-plane options must be supplied as one complete group."
        )

    # Step 1: Select the BHaH OpenMP backend for this example.
    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")

    # Step 2: Define project-level constants and output paths.
    project_name = "photon_single_geodesic_integrator_analytical"
    project_dir = os.path.abspath(os.path.join(args.outdir, project_name))
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
    u_expr, _ = geodesic_data.photon_momentum_to_normalized_quantities()
    normal_observer_log_energy.normal_observer_log_energy(u_expr)

    set_initial_conditions_kernel.register_photon_batch_structs()
    set_initial_conditions_kernel.set_initial_conditions_kernel(normalized_eom=False)
    event_detection_manager_kernel.register_event_plane_parameters()
    interpolation_kernel.interpolation_kernel(SPACETIME)
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(
        geodesic_data.geodesic_eom_rhs_photon_christoffel(),
        geodesic_data.xx,
        rhs_uses_metric_derivatives=False,
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
        ["observer_x", "observer_y", "observer_z"], args.observer_position
    ):
        par.adjust_CodeParam_default(name, value)
    for name, value in zip(
        [
            "observer_look_forward_x",
            "observer_look_forward_y",
            "observer_look_forward_z",
        ],
        args.observer_look_forward,
    ):
        par.adjust_CodeParam_default(name, value)
    for name, value in zip(
        ["observer_up_x", "observer_up_y", "observer_up_z"], args.observer_up
    ):
        par.adjust_CodeParam_default(name, value)
    par.adjust_CodeParam_default("alpha_w", args.observer_fov[0])
    par.adjust_CodeParam_default("alpha_h", args.observer_fov[1])
    par.adjust_CodeParam_default("t_start", args.t_start)
    par.adjust_CodeParam_default("initial_h", 0.1)
    par.adjust_CodeParam_default("r_escape", args.escape_radius)
    par.adjust_CodeParam_default("evolution_measure_max", 3.0)

    terminal_defaults = {
        "terminal_plane_center_x": -1.0e4,
        "terminal_plane_center_y": 0.0,
        "terminal_plane_center_z": 0.0,
        "terminal_plane_normal_x": 1.0,
        "terminal_plane_normal_y": 0.0,
        "terminal_plane_normal_z": 0.0,
        "terminal_plane_up_x": 0.0,
        "terminal_plane_up_y": 0.0,
        "terminal_plane_up_z": 1.0,
        "terminal_plane_min_coord_radius": 0.0,
        "terminal_plane_max_coord_radius": 1.0,
        "terminal_plane_enabled": False,
    }
    if args.terminal_plane_center is not None:
        terminal_defaults.update(
            {
                "terminal_plane_center_x": args.terminal_plane_center[0],
                "terminal_plane_center_y": args.terminal_plane_center[1],
                "terminal_plane_center_z": args.terminal_plane_center[2],
                "terminal_plane_normal_x": args.terminal_plane_normal[0],
                "terminal_plane_normal_y": args.terminal_plane_normal[1],
                "terminal_plane_normal_z": args.terminal_plane_normal[2],
                "terminal_plane_up_x": args.terminal_plane_up[0],
                "terminal_plane_up_y": args.terminal_plane_up[1],
                "terminal_plane_up_z": args.terminal_plane_up[2],
                "terminal_plane_min_coord_radius": args.terminal_plane_radius[0],
                "terminal_plane_max_coord_radius": args.terminal_plane_radius[1],
                "terminal_plane_enabled": True,
            }
        )
    non_terminal_defaults = {
        "non_terminal_plane_center_x": -1.0e4,
        "non_terminal_plane_center_y": 0.0,
        "non_terminal_plane_center_z": 0.0,
        "non_terminal_plane_normal_x": 1.0,
        "non_terminal_plane_normal_y": 0.0,
        "non_terminal_plane_normal_z": 0.0,
        "non_terminal_plane_up_x": 0.0,
        "non_terminal_plane_up_y": 0.0,
        "non_terminal_plane_up_z": 1.0,
        "non_terminal_plane_enabled": False,
    }
    if args.non_terminal_plane_center is not None:
        non_terminal_defaults.update(
            {
                "non_terminal_plane_center_x": args.non_terminal_plane_center[0],
                "non_terminal_plane_center_y": args.non_terminal_plane_center[1],
                "non_terminal_plane_center_z": args.non_terminal_plane_center[2],
                "non_terminal_plane_normal_x": args.non_terminal_plane_normal[0],
                "non_terminal_plane_normal_y": args.non_terminal_plane_normal[1],
                "non_terminal_plane_normal_z": args.non_terminal_plane_normal[2],
                "non_terminal_plane_up_x": args.non_terminal_plane_up[0],
                "non_terminal_plane_up_y": args.non_terminal_plane_up[1],
                "non_terminal_plane_up_z": args.non_terminal_plane_up[2],
                "non_terminal_plane_enabled": True,
            }
        )
    for name, value in {**terminal_defaults, **non_terminal_defaults}.items():
        par.adjust_CodeParam_default(name, value)

    if args.initial_step is not None:
        if args.initial_step <= 0.0:
            raise ValueError("--initial-step must be positive.")
        par.adjust_CodeParam_default("initial_h", args.initial_step)
    if args.rkf45_tolerances is not None:
        absolute_tolerance, relative_tolerance = args.rkf45_tolerances
        if absolute_tolerance <= 0.0 or relative_tolerance <= 0.0:
            raise ValueError("RKF45 tolerances must be positive.")
        par.adjust_CodeParam_default(
            "rkf45_absolute_error_tolerance", absolute_tolerance
        )
        par.adjust_CodeParam_default("rkf45_error_tolerance", relative_tolerance)
    if args.rkf45_step_range is not None:
        h_min, h_max = args.rkf45_step_range
        if h_min <= 0.0 or h_max < h_min:
            raise ValueError("Invalid RKF45 step range.")
        par.adjust_CodeParam_default("rkf45_h_min", h_min)
        par.adjust_CodeParam_default("rkf45_h_max", h_max)
    if args.rkf45_max_retries is not None:
        if args.rkf45_max_retries < 0:
            raise ValueError("--rkf45-max-retries must be nonnegative.")
        par.adjust_CodeParam_default("rkf45_max_retries", args.rkf45_max_retries)
    if args.evolution_measure_max is not None:
        if args.evolution_measure_max <= 0.0:
            raise ValueError("--evolution-measure-max must be positive.")
        par.adjust_CodeParam_default(
            "evolution_measure_max", args.evolution_measure_max
        )

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
    # these CPU-only macro definitions.
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
    vis_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "geodesic_visualizations"
    )
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
