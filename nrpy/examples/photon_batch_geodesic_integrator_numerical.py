r"""
Orchestrator for the CPU numerical-spacetime photon geodesic integrator.

This module constructs the C project that evolves batched photon trajectories
through a numerical spacetime sourced from the combined raytracing ``.bin``
generated from ``two_blackholes_collide.py --raytracing-outputs``.

The generated code keeps the existing Structure-of-Arrays photon pipeline and
adaptive RKF45 stepping, but all metric and Christoffel evaluations now come
from the numerical interpolation helpers instead of analytic spacetime kernels.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

# ##############################################################################
# PART 0: IMPORTS AND PATH SETUP
# ##############################################################################

import argparse
import os
import shutil
import sys

import sympy as sp

# NRPy core and helper modules for C code generation
import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.reference_metric as refmetric

# Physics/Math Generators (Symbolic definitions of geodesics)
from nrpy.equations.general_relativity.geodesics import geodesics as geo

# NRPy BlackHoles@Home (BHaH) infrastructure modules for C project management
from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.infrastructures.BHaH import CodeParameters as CPs
from nrpy.infrastructures.BHaH import Makefile_helpers as Makefile
from nrpy.infrastructures.BHaH import cmdline_input_and_parfiles

# C-Code Builder Functions (Registers specific physical and numerical routines)
from nrpy.infrastructures.BHaH.general_relativity.geodesics import (
    normalization_constraint,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.azimuthal_symmetry_spatial_lagrange_interpolation import (
    register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.numerical_interpolation import (
    register_CFunction_numerical_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.temporal_lagrange_interpolation import (
    register_CFunction_temporal_lagrange_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
    batch_integrator_numerical,
    calculate_and_fill_blueprint_data_universal,
    calculate_ode_rhs_kernel,
    event_detection_manager_kernel,
    find_event_time_and_state,
    handle_source_plane_intersection,
    handle_window_plane_intersection,
    main,
    p0_reverse_kernel,
    rkf45_finalize_and_control_kernel,
    rkf45_stage_update,
    set_initial_conditions_kernel,
)

# Establish the script directory for reliable relative path resolution
script_dir = os.path.dirname(os.path.abspath(__file__))


def _require(condition: bool, message: str) -> None:
    """
    Raise ``ValueError`` if a derived example constraint is violated.

    :param condition: Boolean condition that must hold.
    :param message: Error message describing the violated constraint.
    :raises ValueError: If ``condition`` is false.
    """
    if not condition:
        raise ValueError(message)


# ##############################################################################
# PART 1: MAIN CONFIGURATION
# ##############################################################################

if __name__ == "__main__":
    # Step 1: Configure command-line arguments for the generation pipeline
    parser = argparse.ArgumentParser(
        description=(
            "Generate the Split-Pipeline Photon Geodesic Integrator "
            "(Numerical). Requires a numerical spacetime .bin file generated "
            "by two_blackholes_collide.py."
        ),
        epilog=(
            "To generate the desired numerical spacetime data, run:\n"
            "python3 two_blackholes_collide.py "
            "--raytracing-spacetime T_FINAL GRID_PHYSICAL_SIZE "
            "DIAGNOSTICS_OUTPUT_EVERY "
            "--raytracing-coord-system CoordSystem "
            "--raytracing-Nxx NXX0 NXX1 NXX2\n\n"
            "Then rerun this photon script using the .bin filename printed "
            "to the terminal by two_blackholes_collide.py, e.g.:\n"
            "python3 photon_batch_geodesic_integrator_numerical.py "
            "--bin_name two_blackholes_collide_7p5_7p5_0p25_SinhSpherical_72_12_2.bin "
            "--dataset-coord-system SinhSpherical "
            "--dataset-grid-physical-size 7.5 "
            "--dataset-sinhw 0.4"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bin_name",
        type=str,
        required=True,
        help=(
            "Name of the numerical spacetime .bin file inside "
            "project/raytracing_data."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="The parent directory where the generated C project will reside.",
    )
    parser.add_argument(
        "--dataset-coord-system",
        type=str,
        required=True,
        choices=refmetric.supported_CoordSystems,
        help=(
            "Coordinate system used when generating the numerical spacetime " "dataset."
        ),
    )
    parser.add_argument(
        "--dataset-grid-physical-size",
        type=float,
        required=True,
        help=(
            "Physical grid size used when generating the numerical spacetime "
            "dataset."
        ),
    )
    parser.add_argument(
        "--dataset-sinhw",
        type=float,
        default=None,
        help=(
            "SinhSpherical SINHW value used when generating the numerical "
            "spacetime dataset. Used only with --dataset-coord-system "
            "SinhSpherical."
        ),
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    _require(
        os.path.basename(args.bin_name) == args.bin_name,
        "--bin_name must be a filename only, not a path.",
    )
    _require(
        args.bin_name.endswith(".bin"),
        "--bin_name must name a numerical spacetime .bin file.",
    )
    _require(
        args.dataset_grid_physical_size > 0.0,
        "--dataset-grid-physical-size must be positive.",
    )
    if args.dataset_coord_system == "SinhSpherical":
        _require(
            args.dataset_sinhw is not None,
            "--dataset-sinhw is required for --dataset-coord-system SinhSpherical.",
        )
    else:
        _require(
            args.dataset_sinhw is None,
            "--dataset-sinhw is supported only for --dataset-coord-system SinhSpherical.",
        )

    # Step 2: Define strict project constants and simulation targets
    project_name = "photon_batch_geodesic_integrator_numerical"
    exec_name = "photon_batch_geodesic_integrator_numerical"
    project_dir = os.path.abspath(os.path.join(args.outdir, project_name))

    # Define the physical regime: Target spacetime and particle classification
    SPACETIME = "Numerical"
    integrator_mode = "Numerical"
    PARTICLE = "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"
    dataset_coord_system = args.dataset_coord_system
    enable_simd = False

    # Step 3: Initialize the project directory and select the infrastructure backend
    print(f"Initializing project: {project_name}")
    shutil.rmtree(project_dir, ignore_errors=True)
    os.makedirs(project_dir, exist_ok=True)

    # Instruct NRPy to emit the CPU/OpenMP BHaH pipeline used by the numerical integrator.
    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")
    par.set_parval_from_str(
        "CoordSystem_to_register_CodeParameters", dataset_coord_system
    )

    # Step 4: Build the generic symbolic photon equations consumed by the
    # runtime numerical metric and Christoffel interpolation pipeline.
    print(f" -> Assembling symbolic data for {GEO_KEY}...")
    generic_geodesic_equations = geo.GeodesicEquations.__new__(geo.GeodesicEquations)
    coordinate_symbols = list(sp.symbols("t x y z", real=True))
    geodesic_rhs = generic_geodesic_equations.geodesic_eom_rhs_photon()
    p0_photon = generic_geodesic_equations.hamiltonian_constraint_photon()
    normalization_constraint_expr = (
        generic_geodesic_equations.normalization_constraint()
    )

    # ##########################################################################
    # Step 5: Execute Modules and Register C Functions (Split-Pipeline Order)
    # ##########################################################################
    print(" -> Registering C functions and local CodeParameters...")

    # --- Initialization Phase Kernels ---
    set_initial_conditions_kernel.set_initial_conditions_kernel(SPACETIME)

    p0_reverse_kernel.p0_reverse_kernel(p0_photon)

    # --- Fundamental Constraint Evaluation ---
    normalization_constraint.normalization_constraint(
        normalization_constraint_expr, PARTICLE
    )

    # --- Numerical Interpolation Pipeline ---
    register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation(
        dataset_coord_system, enable_simd=enable_simd, project_dir=project_dir
    )
    register_CFunction_temporal_lagrange_interpolation(
        enable_simd=enable_simd, project_dir=project_dir
    )
    register_CFunction_numerical_interpolation(
        dataset_coord_system, enable_simd=enable_simd, project_dir=project_dir
    )

    # --- Core RKF45 and Geodesic Update Kernels ---
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(geodesic_rhs, coordinate_symbols)
    rkf45_stage_update.rkf45_stage_update()
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel(
        enable_numerical_time_window_step_cap=True
    )

    # --- Event Detection & Boundary Management ---
    find_event_time_and_state.find_event_time_and_state()
    handle_source_plane_intersection.handle_source_plane_intersection()
    handle_window_plane_intersection.handle_window_plane_intersection()
    event_detection_manager_kernel.event_detection_manager_kernel()
    calculate_and_fill_blueprint_data_universal.calculate_and_fill_blueprint_data_universal()

    # --- Infrastructure Helpers ---
    # The numerical interpolation and RKF45 finalization registrations above
    # already pull in the shared slot/time-window helpers on demand. Calling
    # them again here would append duplicate definitions into BHaH_defines.h.
    batch_integrator_numerical.batch_integrator_numerical(SPACETIME)
    main.main(SPACETIME, integrator_mode)

    # --- Native NRPy Cleanup ---
    # Remove helper registrations that the event manager inlines through its
    # prefunc so the project does not emit redundant standalone C files or
    # extra prototypes for functions that are no longer called directly.
    for internal_func in [
        "find_event_time_and_state",
        "handle_source_plane_intersection",
        "handle_window_plane_intersection",
    ]:
        cfc.CFunction_dict.pop(internal_func, None)

    # ##########################################################################
    # Step 5.5: OVERRIDE DEFAULT CODE PARAMETERS
    # ##########################################################################

    print(" -> Overriding desired photon CodeParameters before .par generation...")

    # --------------------------------------------------------------------------
    # Numerical Spacetime Data
    # --------------------------------------------------------------------------
    # The photon executable only needs to know where the validated combined
    # numerical spacetime .bin file lives. Generation/validation of this file is
    # handled outside this photon script.
    numerical_spacetime_bin_path = os.path.abspath(
        os.path.join(
            "project",
            "raytracing_data",
            args.bin_name,
        )
    )

    par.glb_code_params_dict["numerical_spacetime_bin_path"].defaultvalue = (
        numerical_spacetime_bin_path
    )

    # Step 5.5.a: Match the dataset reference-metric defaults used by the
    # combined numerical spacetime file before metadata overrides Nxx/dxx/xxmin/xxmax.
    par.adjust_CodeParam_default("grid_physical_size", args.dataset_grid_physical_size)
    rfm = refmetric.reference_metric[dataset_coord_system]
    for param_name, grid_size_mapping in rfm.grid_physical_size_dict.items():
        if grid_size_mapping == "grid_physical_size":
            par.adjust_CodeParam_default(param_name, args.dataset_grid_physical_size)
        elif grid_size_mapping == "-grid_physical_size":
            par.adjust_CodeParam_default(param_name, -args.dataset_grid_physical_size)
        else:
            raise ValueError(
                f"Unsupported grid_physical_size mapping '{grid_size_mapping}' "
                f"for {dataset_coord_system}:{param_name}"
            )
    if dataset_coord_system == "SinhSpherical":
        par.adjust_CodeParam_default("SINHW", args.dataset_sinhw)

    # --------------------------------------------------------------------------
    #  Lagrange interpolation orders used by the numerical spacetime interpolator
    # --------------------------------------------------------------------------
    par.glb_code_params_dict[
        "numerical_spacetime_spatial_interp_order"
    ].defaultvalue = 3
    par.glb_code_params_dict[
        "numerical_spacetime_temporal_interp_order"
    ].defaultvalue = 3

    # --------------------------------------------------------------------------
    # Execution Initial Conditions
    # --------------------------------------------------------------------------
    par.glb_code_params_dict["t_start"].defaultvalue = 100.0
    par.glb_code_params_dict["scan_density"].defaultvalue = 100

    # --------------------------------------------------------------------------
    # Batch Integrator & Numerical Limits
    # --------------------------------------------------------------------------
    par.glb_code_params_dict["p_t_max"].defaultvalue = 1000.0
    par.glb_code_params_dict["perform_normalization_check"].defaultvalue = True
    par.glb_code_params_dict["r_escape"].defaultvalue = 40.0

    # Maximum coordinate-time step allowed by the RKF45 controller when using
    # numerical spacetime data. This protects the time-window manager from
    # accepting steps that jump across too much numerical data at once.
    par.glb_code_params_dict["rkf45_max_delta_t"].defaultvalue = 0.5

    # Time-slot manager parameters. These must be chosen consistently with the
    # time range covered by numerical_spacetime_bin_path.
    par.glb_code_params_dict["slot_manager_delta_t"].defaultvalue = 2.0
    par.glb_code_params_dict["slot_manager_t_min"].defaultvalue = 5.0

    # --------------------------------------------------------------------------
    # Source Plane Geometric Mapping
    # --------------------------------------------------------------------------
    par.glb_code_params_dict["source_plane_center_x"].defaultvalue = -100.0
    par.glb_code_params_dict["source_plane_center_y"].defaultvalue = 0.0
    par.glb_code_params_dict["source_plane_center_z"].defaultvalue = 0.0

    par.glb_code_params_dict["source_plane_normal_x"].defaultvalue = 1.0
    par.glb_code_params_dict["source_plane_normal_y"].defaultvalue = 0.0
    par.glb_code_params_dict["source_plane_normal_z"].defaultvalue = 0.0

    par.glb_code_params_dict["source_r_max"].defaultvalue = 30.0
    par.glb_code_params_dict["source_r_min"].defaultvalue = 0.0

    par.glb_code_params_dict["source_up_vec_x"].defaultvalue = 0.0
    par.glb_code_params_dict["source_up_vec_y"].defaultvalue = 1.0
    par.glb_code_params_dict["source_up_vec_z"].defaultvalue = 0.0

    # --------------------------------------------------------------------------
    # Camera Window Geometric Mapping
    # --------------------------------------------------------------------------
    par.glb_code_params_dict["camera_pos_x"].defaultvalue = 11.0
    par.glb_code_params_dict["camera_pos_y"].defaultvalue = 0.0
    par.glb_code_params_dict["camera_pos_z"].defaultvalue = 0.0

    par.glb_code_params_dict["original_window_center_x"].defaultvalue = 10.0
    par.glb_code_params_dict["original_window_center_y"].defaultvalue = 0.0
    par.glb_code_params_dict["original_window_center_z"].defaultvalue = 0.0

    par.glb_code_params_dict["window_height"].defaultvalue = 1.0
    par.glb_code_params_dict["window_width"].defaultvalue = 1.0

    par.glb_code_params_dict["window_up_vec_x"].defaultvalue = 0.0
    par.glb_code_params_dict["window_up_vec_y"].defaultvalue = 0.0
    par.glb_code_params_dict["window_up_vec_z"].defaultvalue = 1.0

    # Numerical photon script is CPU/OpenMP-only, so use fixed CPU tiling.
    par.glb_code_params_dict["window_tiles_width"].defaultvalue = 1
    par.glb_code_params_dict["window_tiles_height"].defaultvalue = 1

    # --------------------------------------------------------------------------
    # RKF45 Adaptive Control Tolerances
    # --------------------------------------------------------------------------
    par.glb_code_params_dict["numerical_initial_h"].defaultvalue = 0.05
    par.glb_code_params_dict["rkf45_absolute_error_tolerance"].defaultvalue = 1.0e-8
    par.glb_code_params_dict["rkf45_error_tolerance"].defaultvalue = 1.0e-8
    par.glb_code_params_dict["rkf45_h_max"].defaultvalue = 10.0
    par.glb_code_params_dict["rkf45_h_min"].defaultvalue = 1.0e-15

    print(f" -> Numerical spacetime .bin path: {numerical_spacetime_bin_path}")
    print(f" -> Dataset coordinate system: {dataset_coord_system}")
    print(f" -> Dataset grid physical size: {args.dataset_grid_physical_size}")
    if args.dataset_sinhw is not None:
        print(f" -> Dataset SINHW: {args.dataset_sinhw}")

    # Step 6: Generate C Code for Parameter Handling
    print(" -> Generating parameter handling code...")
    CPs.write_CodeParameters_h_files(project_dir=project_dir, set_commondata_only=True)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()
    cmdline_input_and_parfiles.generate_default_parfile(
        project_dir=project_dir, project_name=project_name
    )
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
        project_name=project_name,
        cmdline_inputs=[
            "source_r_min",
            "source_r_max",
            "window_width",
            "window_height",
            "window_tiles_width",
            "window_tiles_height",
            "scan_density",
        ],
    )

    # ##########################################################################
    # Step 7: Assemble Final C Project Infrastructure
    # ##########################################################################
    print(" -> Assembling C project on disk...")

    # Map the shared CUDA-flavored kernel abstraction macros onto CPU/OpenMP.
    cpu_macros = {
        # --- Memory Access Redirection ---
        "ReadCUDA(ptr)": "#define ReadCUDA(ptr) (*(ptr))\n",
        "WriteCUDA(ptr, val)": "#define WriteCUDA(ptr, val) (*(ptr) = (val))\n",
        # --- Device Memory Fallbacks ---
        "BHAH_MALLOC_DEVICE(a, sz)": "#define BHAH_MALLOC_DEVICE(a, sz) BHAH_MALLOC(a, sz)\n",
        "BHAH_FREE_DEVICE(a)": "#define BHAH_FREE_DEVICE(a) BHAH_FREE(a)\n",
        # --- Basic Arithmetic Intrinsics (Required by RKF45 Kernels) ---
        "MulCUDA(a, b)": "#define MulCUDA(a, b) ((a) * (b))\n",
        "DivCUDA(a, b)": "#define DivCUDA(a, b) ((a) / (b))\n",
        "AddCUDA(a, b)": "#define AddCUDA(a, b) ((a) + (b))\n",
        "FusedMulAddCUDA(a, b, c)": "#define FusedMulAddCUDA(a, b, c) ((a) * (b) + (c))\n",
        # --- Hardware-Specific Math Redirection ---
        "AbsCUDA(val)": "#define AbsCUDA(val) fabs(val)\n",
        "SqrtCUDA(val)": "#define SqrtCUDA(val) sqrt(val)\n",
        "PowCUDA(a, b)": "#define PowCUDA(a, b) pow(a, b)\n",
        # --- Function Decorators ---
        "BHAH_HD_FUNC": "#define BHAH_HD_FUNC\n",
        "BHAH_HD_INLINE": "#define BHAH_HD_INLINE static inline\n",
        # --- Parallelization & Scope Helpers ---
        "BHAH_WARP_ATOMIC_ADD(ptr, val)": '#define BHAH_WARP_ATOMIC_ADD(ptr, val) _Pragma("omp atomic") *(ptr) += (val)\n',
        "GLOBAL_COMMONDATA_EXTERN": "// CPU passes commondata by reference, no global needed.\n",
        "BHAH_DEVICE_SYNC()": "#define BHAH_DEVICE_SYNC() do {} while(0)\n",
    }

    # Match the combined numerical spacetime payload layout produced by the
    # trusted two-black-hole raytracing pipeline, which stores NGHOSTS = 3.
    BHaH_defines_h.output_BHaH_defines_h(
        project_dir=project_dir,
        enable_rfm_precompute=False,
        fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
        supplemental_defines_dict=cpu_macros,
    )

    compiler = "gcc"
    cflags = [
        "-fopenmp",
        "-O3",
        "-fno-omit-frame-pointer",
        "-DDEBUG",
        "-Wno-stringop-truncation",
    ]
    libs = ["-lm"]
    ext = "c"

    print(" -> Generating Makefile ")

    # Set optimization option string
    opt_option = "fast"

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
    # ##########################################################################
    # PART 2: FINALIZE
    # ##########################################################################

    # Define the directory containing the visualization assets relative to the repository root
    vis_dir = os.path.join("nrpy", "examples", "geodesic_visualizations")

    # Locate the visualization script
    vis_script_src = os.path.join(vis_dir, "visualize_lensed_image.py")
    config_src = os.path.join(vis_dir, "blueprint_config_and_schema.py")
    render_src = os.path.join(vis_dir, "render_lensed_image.py")
    blueprint_analysis_src = os.path.join(vis_dir, "blueprint_analysis.py")
    combined_helper_src = os.path.join(
        "nrpy",
        "examples",
        "combined_raytracing_bin_helper.py",
    )

    for script_src in (
        vis_script_src,
        config_src,
        render_src,
        blueprint_analysis_src,
    ):
        shutil.copy(script_src, project_dir)

    # The inner disk radius ensures the texture mapping aligns with the computed initial conditions.
    c_r_min = float(par.glb_code_params_dict["source_r_min"].defaultvalue)

    # The outer disk radius ensures the mathematical bounds match the physical source geometry.
    c_r_max = float(par.glb_code_params_dict["source_r_max"].defaultvalue)

    # The camera window width scales the horizontal projection bounds.
    c_window_width = float(par.glb_code_params_dict["window_width"].defaultvalue)

    # The camera window height scales the vertical projection bounds.
    c_window_height = float(par.glb_code_params_dict["window_height"].defaultvalue)

    # The terminal execution string now references the local copy of the script.
    # Get tile dimensions to map the partitioned geometry blocks.
    c_tiles_width = int(par.glb_code_params_dict["window_tiles_width"].defaultvalue)
    c_tiles_height = int(par.glb_code_params_dict["window_tiles_height"].defaultvalue)

    # Set the baseline pixel resolution for the output render.
    c_pixel_width = 600
    data_request_file = "numerical_spacetime_data_request.json"
    data_prep_command = f"python3 combined_raytracing_bin_helper.py {data_request_file}"

    vis_command = (
        f"python3 visualize_lensed_image.py "
        f"--source_r_min {c_r_min} "
        f"--source_r_max {c_r_max} "
        f"--window_width {c_window_width} "
        f"--window_height {c_window_height} "
        f"--window_tiles_width {c_tiles_width} "
        f"--window_tiles_height {c_tiles_height} "
        f"--pixel_width {c_pixel_width}"
    )

    print(f"Finished! Now go into {project_dir}.")
    print(
        f"    Parameter file can be found at {os.path.join(project_dir, f'{project_name}.par')}\n"
    )
    print("    To prepare the required numerical spacetime data, run:")
    print(f"    {data_prep_command}\n")
    print("    Then build and run the photon executable:")
    print("    make")
    print(f"    ./{exec_name}\n")
    print(
        "    To generate the lensed image after running the C executable, ensure you have the required Python packages:"
    )
    print("    pip install matplotlib numpy numba Pillow\n")
    print(
        "    Then, execute the visualization script directly from the project directory:"
    )
    print(f"    {vis_command}\n")

    blueprint_command = (
        f"python3 blueprint_analysis.py "
        f"--window_tiles_width {c_tiles_width} "
        f"--window_tiles_height {c_tiles_height} "
        f"--window_width {c_window_width} "
        f"--window_height {c_window_height}"
    )

    print("    To run the blueprint diagnostic:")
    print(f"    {blueprint_command}\n")
