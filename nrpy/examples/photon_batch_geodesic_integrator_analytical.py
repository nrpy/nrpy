r"""
Orchestrator for the Split-Pipeline Photon Geodesic Integrator.

This module constructs a high-performance C project for simulating photon trajectories
in curved spacetimes.

The generated code employs a Structure of Arrays (SoA) memory layout and an
adaptive RKF45 integration scheme managed by a lock-free TimeSlotManager system.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com

"""

import argparse
import os
import shutil

# NRPy core and helper modules for C code generation
import nrpy.c_function as cfc
import nrpy.helpers.generic as gh
import nrpy.params as par

# Physics/Math Generators (Symbolic definitions of spacetimes and geodesics)
from nrpy.equations.general_relativity.geodesics import analytic_spacetimes as anasp
from nrpy.equations.general_relativity.geodesics import geodesics as geo

# NRPy BlackHoles@Home (BHaH) infrastructure modules for C project management
from nrpy.infrastructures.BHaH import BHaH_defines_h, BHaH_device_defines_h
from nrpy.infrastructures.BHaH import CodeParameters as CPs
from nrpy.infrastructures.BHaH import Makefile_helpers as Makefile
from nrpy.infrastructures.BHaH import cmdline_input_and_parfiles

# C-Code Builder Functions (Registers specific physical and numerical routines)
from nrpy.infrastructures.BHaH.general_relativity.geodesics import (
    connections,
    conserved_quantities,
    g4DD_metric,
    normalization_constraint,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
    batch_integrator_analytical,
    calculate_and_fill_blueprint_data_universal,
    calculate_ode_rhs_kernel,
    event_detection_manager_kernel,
    find_event_time_and_state,
    handle_non_terminal_plane_intersection,
    handle_terminal_plane_intersection,
    interpolation_kernel,
    main_batch,
    normal_observer_log_energy,
    rkf45_finalize_and_control_kernel,
    rkf45_stage_update,
    set_initial_conditions_kernel,
    time_slot_manager_helpers,
)

if __name__ == "__main__":
    import sys

    # Step 1: Configure command-line arguments for the generation pipeline
    parser = argparse.ArgumentParser(
        description=(
            "Generate a tiled analytical photon raytracing project. "
            "The project emits v6 light-blueprint records with normalized "
            "image-sample coordinates."
        )
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="Parent directory for the generated C project.",
    )
    parser.add_argument(
        "--cuda",
        action="store_true",
        help="Generate CUDA code; omit this option for the OpenMP project.",
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
        metavar=("X", "Y", "Z"),
        help="Observer up/roll direction seed.",
    )
    parser.add_argument(
        "--observer-fov",
        nargs=2,
        type=float,
        required=True,
        metavar=("ALPHA_W", "ALPHA_H"),
        help="Observer width and height fields of view, in radians.",
    )
    parser.add_argument(
        "--scan-density",
        type=int,
        required=True,
        metavar="SAMPLES_PER_TILE_WIDTH",
        help="Width-side ray samples per tile; height-side count is derived.",
    )
    parser.add_argument(
        "--tile-counts",
        nargs=2,
        type=int,
        default=[1, 1],
        metavar=("TILES_WIDTH", "TILES_HEIGHT"),
        help="Number of angular image tiles; default: 1 1.",
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
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    if any(value <= 0.0 for value in args.observer_fov):
        parser.error("--observer-fov values must be positive.")
    if args.scan_density < 1:
        parser.error("--scan-density must be a positive integer.")
    if any(value < 1 for value in args.tile_counts):
        parser.error("--tile-counts values must be positive integers.")
    if args.escape_radius <= 0.0:
        parser.error("--escape-radius must be positive.")

    terminal_plane_values = (
        args.terminal_plane_center,
        args.terminal_plane_normal,
        args.terminal_plane_up,
        args.terminal_plane_radius,
    )
    if any(value is not None for value in terminal_plane_values) and not all(
        value is not None for value in terminal_plane_values
    ):
        parser.error("Terminal-plane options must be supplied as one complete group.")
    non_terminal_plane_values = (
        args.non_terminal_plane_center,
        args.non_terminal_plane_normal,
        args.non_terminal_plane_up,
    )
    if any(value is not None for value in non_terminal_plane_values) and not all(
        value is not None for value in non_terminal_plane_values
    ):
        parser.error(
            "Nonterminal-plane options must be supplied as one complete group."
        )

    # Step 2: Define strict project constants and simulation targets
    project_name = "photon_batch_geodesic_integrator_analytical"
    exec_name = "photon_batch_geodesic_integrator_analytical"
    project_dir = os.path.abspath(os.path.join(args.outdir, project_name))

    # Step 2.a: Define the target spacetime and particle type.
    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"

    # Step 3: Initialize the project directory and select the infrastructure backend
    print(f"Initializing project: {project_name}")
    shutil.rmtree(project_dir, ignore_errors=True)
    os.makedirs(project_dir, exist_ok=True)

    # Step 3.a: Select the BHaH infrastructure for code generation.
    par.set_parval_from_str("Infrastructure", "BHaH")

    # Step 3.b: Map the CUDA flag onto the NRPy parallelization parameter.
    parallelization_mode = "cuda" if args.cuda else "openmp"
    par.set_parval_from_str("parallelization", parallelization_mode)

    # Step 4: Acquire symbolic physics expressions.
    print(f" -> Acquiring symbolic data for {GEO_KEY}...")
    metric_data = anasp.Analytic_Spacetimes[SPACETIME]
    geodesic_data = geo.Geodesic_Equations[GEO_KEY]

    # Step 5: Register C functions in split-pipeline order.
    print(" -> Registering C functions and local CodeParameters...")

    # Step 5.a: Register batch state definitions and initialization kernels.
    set_initial_conditions_kernel.register_photon_batch_structs()
    set_initial_conditions_kernel.set_initial_conditions_kernel(normalized_eom=False)

    u_expr, _ = geo.GeodesicEquations.photon_momentum_to_normalized_quantities()
    normal_observer_log_energy.normal_observer_log_energy(u_expr)

    # Step 5.b: Register fundamental tensor and diagnostic kernels.
    g4DD_metric.g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)
    connections.connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)
    conserved_quantities.conserved_quantities(SPACETIME, PARTICLE)
    normalization_constraint.normalization_constraint(
        geodesic_data.norm_constraint_expr, PARTICLE
    )

    # Step 5.c: Register RKF45 evolution kernels.
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

    # Step 5.d: Register event-detection and boundary-intersection kernels.
    find_event_time_and_state.find_event_time_and_state()
    handle_terminal_plane_intersection.handle_terminal_plane_intersection()
    handle_non_terminal_plane_intersection.handle_non_terminal_plane_intersection()
    event_detection_manager_kernel.event_detection_manager_kernel(normalized_eom=False)
    calculate_and_fill_blueprint_data_universal.calculate_and_fill_blueprint_data_universal(
        normalized_eom=False
    )

    # Step 5.e: Register project-level orchestration helpers.
    time_slot_manager_helpers.time_slot_manager_helpers()
    batch_integrator_analytical.batch_integrator_analytical(SPACETIME)
    main_batch.main(SPACETIME, normalized_eom=False)

    # Step 5.f: Remove helper registrations emitted only through other kernels.
    # The event manager emits its local event helpers through prefunc, so keeping
    # standalone registrations would add unused source files and prototypes. The
    # metric and connection helpers are likewise inlined where they are consumed
    # in this generated project. Removing these entries prevents redundant emitted
    # wrappers in both OpenMP and CUDA builds, and also avoids duplicate CUDA
    # definitions during device linking.
    for internal_func in [
        "find_event_time_and_state",
        "handle_terminal_plane_intersection",
        "handle_non_terminal_plane_intersection",
        f"g4DD_metric_{SPACETIME}",
        f"connections_{SPACETIME}",
    ]:
        cfc.CFunction_dict.pop(internal_func, None)

    # Step 6: Override CodeParameter defaults before parfile generation.
    print(" -> Overriding desired CodeParameters before .par generation...")

    # Step 6.a: Set analytic spacetime defaults in units M_scale=1.
    # a_spin is dimensional Kerr a=J/M, so a_spin=0.9 means a/M=0.9 here.
    par.adjust_CodeParam_default("M_scale", 1.0)
    par.adjust_CodeParam_default("a_spin", 0.9)

    # Step 6.b: Set batch-integrator and numerical-limit defaults.
    par.adjust_CodeParam_default("evolution_measure_max", 3.0)
    par.adjust_CodeParam_default("perform_conservation_check", True)
    par.adjust_CodeParam_default("r_escape", args.escape_radius)
    par.adjust_CodeParam_default("slot_manager_delta_t", 100.0)
    par.adjust_CodeParam_default("slot_manager_t_min", -1000.0)

    # Step 6.c: Set observer and independent plane geometry from the CLI.
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

    # Step 6.e: Set angular tile-grid defaults for the selected backend.
    par.adjust_CodeParam_default("tiles_width", args.tile_counts[0])
    par.adjust_CodeParam_default("tiles_height", args.tile_counts[1])

    # Step 6.f: Set RKF45 controller defaults.
    par.adjust_CodeParam_default("initial_h", 0.1)
    par.adjust_CodeParam_default("rkf45_absolute_error_tolerance", 1e-10)
    par.adjust_CodeParam_default("rkf45_error_tolerance", 1e-10)
    par.adjust_CodeParam_default("rkf45_h_max", 10.0)
    par.adjust_CodeParam_default("rkf45_h_min", 1e-15)

    # Step 6.g: Set initial-condition defaults.
    par.adjust_CodeParam_default("t_start", args.t_start)

    # Step 6.h: Set width-side ray-sample density for the selected backend.
    # The initializer derives the height-side density from the fields of view
    # and tile-count aspect ratio. Final image resolution remains a renderer
    # setting, not a photon-grid parameter.
    par.adjust_CodeParam_default("scan_density", args.scan_density)

    if args.initial_step is not None:
        if args.initial_step <= 0.0:
            parser.error("--initial-step must be positive.")
        par.adjust_CodeParam_default("initial_h", args.initial_step)
    if args.rkf45_tolerances is not None:
        absolute_tolerance, relative_tolerance = args.rkf45_tolerances
        if absolute_tolerance <= 0.0 or relative_tolerance <= 0.0:
            parser.error("RKF45 tolerances must be positive.")
        par.adjust_CodeParam_default(
            "rkf45_absolute_error_tolerance", absolute_tolerance
        )
        par.adjust_CodeParam_default("rkf45_error_tolerance", relative_tolerance)
    if args.rkf45_step_range is not None:
        h_min, h_max = args.rkf45_step_range
        if h_min <= 0.0 or h_max < h_min:
            parser.error("RKF45 step range requires 0 < H_MIN <= H_MAX.")
        par.adjust_CodeParam_default("rkf45_h_min", h_min)
        par.adjust_CodeParam_default("rkf45_h_max", h_max)
    if args.rkf45_max_retries is not None:
        if args.rkf45_max_retries < 0:
            parser.error("--rkf45-max-retries must be nonnegative.")
        par.adjust_CodeParam_default("rkf45_max_retries", args.rkf45_max_retries)
    if args.evolution_measure_max is not None:
        if args.evolution_measure_max <= 0.0:
            parser.error("--evolution-measure-max must be positive.")
        par.adjust_CodeParam_default(
            "evolution_measure_max", args.evolution_measure_max
        )

    # Step 6.i: Generate C code for parameter handling.
    print(" -> Generating parameter handling code...")
    CPs.write_CodeParameters_h_files(project_dir=project_dir, set_commondata_only=True)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()
    cmdline_input_and_parfiles.generate_default_parfile(
        project_dir=project_dir, project_name=project_name
    )
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
        project_name=project_name,
        cmdline_inputs=[
            "terminal_plane_min_coord_radius",
            "terminal_plane_max_coord_radius",
            "tiles_width",
            "tiles_height",
            "scan_density",
        ],
    )

    # Step 7: Assemble the generated C project.
    print(" -> Assembling C project on disk...")

    if parallelization_mode == "cuda":
        if "DEVICE_THREAD_MACROS" not in par.glb_extras_dict:
            par.glb_extras_dict["DEVICE_THREAD_MACROS"] = {}

        # Use larger default CUDA thread blocks for the register-heavy kernels
        # emitted by this project.
        par.glb_extras_dict["DEVICE_THREAD_MACROS"].update(
            {
                "BHAH_THREADS_IN_X_DIR_DEFAULT": 256,
                "BHAH_THREADS_IN_Y_DIR_DEFAULT": 1,
                "BHAH_THREADS_IN_Z_DIR_DEFAULT": 1,
            }
        )

        cuda_macros = {
            "BHAH_HD_FUNC": "#define BHAH_HD_FUNC __device__\n",
            "BHAH_HD_INLINE": "#define BHAH_HD_INLINE __device__ __inline__\n",
            "BHAH_WARP_ATOMIC_ADD(ptr, val)": (
                "#define BHAH_WARP_ATOMIC_ADD(ptr, val) atomicAdd(ptr, val)\n"
            ),
            "GLOBAL_COMMONDATA_EXTERN": "extern __constant__ commondata_struct d_commondata;\n",
        }
        BHaH_defines_h.output_BHaH_defines_h(
            project_dir=project_dir,
            enable_rfm_precompute=False,
            supplemental_defines_dict=cuda_macros,
        )
        BHaH_device_defines_h.output_device_headers(project_dir=project_dir)

        print(" -> Copying hardware intrinsics...")
        gh.copy_files(
            package="nrpy.helpers",
            filenames_list=["cuda_intrinsics.h"],
            project_dir=project_dir,
            subdirectory=".",
        )

        compiler = "nvcc"
        cflags = ["-lcudart", "-DUSE_GPU", "-rdc=true", "-DDEBUG"]
        libs = ["-lm", "-lcudart"]
        ext = "cu"
    else:
        # Step 7.a: Map the shared CUDA-style helpers onto OpenMP.
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
                "// CPU passes commondata by reference, no global needed.\n"
            ),
            "BHAH_DEVICE_SYNC()": "#define BHAH_DEVICE_SYNC() do {} while(0)\n",
        }

        BHaH_defines_h.output_BHaH_defines_h(
            project_dir=project_dir,
            enable_rfm_precompute=False,
            supplemental_defines_dict=cpu_macros,
        )

        compiler = "gcc"
        cflags = ["-fopenmp", "-O3", "-DDEBUG", "-Wno-stringop-truncation"]
        libs = ["-lm"]
        ext = "c"

    print(" -> Generating Makefile")

    # Step 7.b: Select the Makefile optimization profile.
    opt_option = "nvcc" if parallelization_mode == "cuda" else "fast"

    Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=exec_name,
        compiler_opt_option=opt_option,
        addl_CFLAGS=cflags,
        addl_libraries=libs,
        CC=compiler,
        src_code_file_ext=ext,
    )
    # Step 8: Copy the v6 blueprint schema/reader/renderer/diagnostic helpers
    # and print usage instructions.
    vis_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "geodesic_visualizations"
    )
    vis_script_src = os.path.join(vis_dir, "visualize_lensed_image.py")
    config_src = os.path.join(vis_dir, "blueprint_config_and_schema.py")
    render_src = os.path.join(vis_dir, "render_lensed_image.py")
    blueprint_io_src = os.path.join(vis_dir, "blueprint_io.py")
    blueprint_analysis_src = os.path.join(vis_dir, "blueprint_analysis.py")

    for script_src in (
        vis_script_src,
        config_src,
        render_src,
        blueprint_io_src,
        blueprint_analysis_src,
    ):
        shutil.copy(script_src, project_dir)

    # Step 8.a: Build visualization command arguments from CodeParameter defaults.
    c_r_min = float(
        par.glb_code_params_dict["terminal_plane_min_coord_radius"].defaultvalue
    )
    c_r_max = float(
        par.glb_code_params_dict["terminal_plane_max_coord_radius"].defaultvalue
    )
    c_scan_density = int(par.glb_code_params_dict["scan_density"].defaultvalue)
    # Visualization output width; angular ray-sample density comes from
    # CodeParameters and the renderer maps normalized fractions to output pixels.
    c_pixel_width = 600
    parfile_path = os.path.join(project_dir, f"{project_name}.par")

    vis_command_parts = ["python3 visualize_lensed_image.py"]
    if args.terminal_plane_radius is not None:
        vis_command_parts.append(f"--terminal-plane-radius {c_r_min} {c_r_max}")
    vis_command_parts.append(f"--pixel-width {c_pixel_width}")
    vis_command = " ".join(vis_command_parts)

    print(
        f"Finished! Now go into {project_dir} and type `make` to build, "
        f"then ./{exec_name} to run."
    )
    print(f"    Parameter file can be found at {parfile_path}\n")
    print(
        "    To generate the lensed image after running the C executable, "
        "ensure you have the required Python packages:"
    )
    print("    pip install matplotlib numpy numba Pillow\n")
    print(
        "    Then, execute the visualization script directly from the project directory:"
    )
    print(f"    {vis_command}\n")
    print(
        "    Ray sampling uses "
        f"{c_scan_density} width-side samples per tile; height-side sampling "
        "is derived from the fields of view and tile-grid aspect ratio.\n"
    )

    blueprint_command = "python3 blueprint_analysis.py"

    print("    To run the blueprint diagnostic:")
    print(f"    {blueprint_command}\n")
