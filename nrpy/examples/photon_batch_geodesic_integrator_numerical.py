r"""
Orchestrator for the Split-Pipeline Photon Geodesic Integrator.

This module constructs a high-performance C project for simulating photon trajectories
in curved spacetimes.

The generated code employs a Structure of Arrays (SoA) memory layout and an
adaptive RKF45 integration scheme managed by a lock-free TimeSlotManager system.

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
from nrpy.infrastructures.BHaH.general_relativity.geodesics import normalization_constraint
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation import (
    azimuthal_symmetry_spatial_lagrange_interpolation,
    numerical_interpolation,
    temporal_lagrange_interpolation,
    time_window_manger_numerical,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
    batch_integrator_numerical,
    calculate_and_fill_blueprint_data_universal,
    calculate_ode_rhs_kernel,
    event_detection_manager_kernel,
    find_event_time_and_state,
    handle_source_plane_intersection,
    handle_window_plane_intersection,
    interpolation_kernel,
    main,
    p0_reverse_kernel,
    rkf45_finalize_and_control_kernel,
    rkf45_stage_update,
    set_initial_conditions_kernel,
    time_slot_manager_helpers,
)

# Import python helper function to check/generate required .bin interpolation data
from nrpy.infrastructures.BHaH.diagnostics.combined_raytracing_bin_helper import (
    ensure_required_combined_bin,
)

# Establish the script directory for reliable relative path resolution
script_dir = os.path.dirname(os.path.abspath(__file__))

# ##############################################################################
# PART 1: MAIN CONFIGURATION
# ##############################################################################

if __name__ == "__main__":
    # Step 1: Configure command-line arguments for the generation pipeline
    parser = argparse.ArgumentParser(
        description="Generate the Split-Pipeline Photon Geodesic Integrator (Numerical)."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="The parent directory where the generated C project will reside.",
    )
    args = parser.parse_args()

    # Step 2: Define strict project constants and simulation targets
    project_name = "photon_batch_geodesic_integrator_numerical"
    exec_name = "photon_batch_geodesic_integrator_numerical"
    project_dir = os.path.abspath(os.path.join(args.outdir, project_name))
    blueprint_path = os.path.join(project_dir, "light_blueprint.bin")

    # Define the physical regime: Target spacetime and particle classification
    SPACETIME = "Numerical"
    integrator_mode = SPACETIME
    PARTICLE = "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"

    # Step 3: Initialize the project directory and select the infrastructure backend
    print(f"Initializing project: {project_name}")
    shutil.rmtree(project_dir, ignore_errors=True)
    os.makedirs(project_dir, exist_ok=True)

    # Instruct NRPy to use the BHaH infrastructure for macro expansions and SoA layouts
    par.set_parval_from_str("Infrastructure", "BHaH")

    # Step 4: Acquire Symbolic Physics Expressions
    print(f" -> Acquiring symbolic data for {GEO_KEY}...")
    geodesic_data = geo.Geodesic_Equations[GEO_KEY]

    # ##########################################################################
    # Step 5: Execute Modules and Register C Functions (Split-Pipeline Order)
    # ##########################################################################
    print(" -> Registering C functions and local CodeParameters...")

    # --- Initialization Phase Kernels ---
    set_initial_conditions_kernel.set_initial_conditions_kernel(SPACETIME)

    if geodesic_data.p0_photon is None:
        raise ValueError(f"p0_photon is None for {GEO_KEY}")
    p0_reverse_kernel.p0_reverse_kernel(geodesic_data.p0_photon)

    # --- Fundamental Tensor Calculations (VRAM Persisted) ---
    normalization_constraint.normalization_constraint(
        geodesic_data.norm_constraint_expr, PARTICLE
    )

    # --- Core Pipeline Kernels (The RKF45 Modular Loop) ---
    interpolation_kernel.interpolation_kernel(SPACETIME)
    azimuthal_symmetry_spatial_lagrange_interpolation.azimuthal_symmetry_spatial_lagrange_interpolation()
    temporal_lagrange_interpolationtemporal_lagrange_interpolation()
    numerical_interpolation.numerical_interpolation()
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(
        geodesic_data.geodesic_rhs, geodesic_data.xx
    )
    rkf45_stage_update.rkf45_stage_update()
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel()

    # --- Event Detection & Boundary Management ---
    find_event_time_and_state.find_event_time_and_state()
    handle_source_plane_intersection.handle_source_plane_intersection()
    handle_window_plane_intersection.handle_window_plane_intersection()
    event_detection_manager_kernel.event_detection_manager_kernel()
    calculate_and_fill_blueprint_data_universal.calculate_and_fill_blueprint_data_universal()

    # --- Infrastructure Helpers ---
    time_slot_manager_helpers.time_slot_manager_helpers()
        time_window_manger_numerical.time_window_manger_numerical()
    batch_integrator_numerical.batch_integrator_numerical(SPACETIME)
    main.main(SPACETIME, integrator_mode)

    # --- Native NRPy Cleanup ---
    # Remove the inline helper functions from the global CFunction dictionary.
    # Because the manager kernel incorporates their logic directly via its prefunc,
    # popping them prevents the infrastructure from generating standalone .cu files
    # and prevents 'ghost' prototypes in BHaH_function_prototypes.h. This eliminates
    # nvlink multiple-definition conflicts during relocatable device code linking.
    for internal_func in [
        "find_event_time_and_state",
        "handle_source_plane_intersection",
        "handle_window_plane_intersection",
    ]:
        cfc.CFunction_dict.pop(internal_func, None)

    # ##########################################################################
    # Step 5.5: OVERRIDE DEFAULT CODE PARAMETERS
    # ##########################################################################

    # Step 5.5.a: Ensure the required numerical spacetime data file exists.
    combined_bin_location = os.path.abspath(
        os.path.join(
            args.outdir,
            "raytracing_data",
            "combined_raytracing_data.bin",
        )
    )

    required_combined_bin_metadata = {
        "combined_file": {
            "format_magic": "NRPYRTSTACK4D",
            "combined_format_version": 1,
            "source_format_version": 1,
            "endianness": "little",
            "serialized_real_bytes": 8,
        },
        "grid": {
            "CoordSystem": "Spherical",
            "Nxx": [72, 12, 2],
            "num_grids": 1,
            "payload_includes_ghost_zones": 0,
            "target_basis": "Cartesian",
            "grid_physical_size": 7.5,
        },
        "payload": {
            "format_name": "Cartesian g4DD+Gamma4UDD",
            "payload_layout": "time_major_stage1_aos",
            "loop_order": "i2maj_i0fast",
            "point_record_real_count": 53,
            "point_record_bytes": 424,
            "record_component_count": 53,
            "metric_component_count": 10,
            "christoffel_component_count": 40,
        },
        "time": {
            "t_start": 0.0,
            "t_final": 7.5,
            "dt": 0.25,
            "absolute_tolerance": 1.0e-12,
        },
        "spatial_lookup": {
            "spatial_lookup_mode": "coordinate_table_only",
            "axisymmetry_enabled": True,
            "axisymmetry_axis": "z",
            "requires_axisymmetry_rotation": True,
        },
        "two_blackholes_run": {
            "floating_point_precision": "double",
            "parallelization": "openmp",
            "raytracing_outputs_enabled": True,
            "BH1_mass": 0.5,
            "BH2_mass": 0.5,
            "BH1_posn_z": 0.5,
            "BH2_posn_z": -0.5,
            "GammaDriving_eta": 1.0,
            "outer_bc_type": "radiation",
            "diagnostics_output_every": 0.25,
        },
        "generation": {
            "project_name": "two_blackholes_collide",
            "python_executable": sys.executable,
            "make_command": ["make"],
            "two_blackholes_example_script": os.path.join(
                script_dir,
                "two_blackholes_collide.py",
            ),
            "combine_raytracing_time_slices_script": os.path.join(
                os.path.dirname(script_dir),
                "infrastructures",
                "BHaH",
                "diagnostics",
                "combine_raytracing_time_slices.py",
            ),
            "stage1_raytracing_output_dir": os.path.abspath(
                os.path.join(
                    args.outdir,
                    "two_blackholes_collide",
                )
            ),
            "stage1_raytracing_pattern": "raytracing_data_t*.bin",
            "executable_name": "two_blackholes_collide",
        },
    }

    state_of_bin, combined_bin_location = ensure_required_combined_bin(
        required_metadata=required_combined_bin_metadata,
        combined_bin_location=combined_bin_location,
    )
    print(f" -> Numerical spacetime data: {state_of_bin}")
    print(f" -> Combined raytracing data path: {combined_bin_location}")

    print(" -> Overriding desired CodeParameters before .par generation...")

    # .bin Path for interpolation
    par.glb_code_params_dict["numerical_spacetime_bin_path"].defaultvalue = combined_bin_location
    # Batch Integrator & Numerical Limits
    par.glb_code_params_dict["p_t_max"].defaultvalue = 1000.0
    par.glb_code_params_dict["perform_conservation_check"].defaultvalue = True
    par.glb_code_params_dict["r_escape"].defaultvalue = 100.0
    par.glb_code_params_dict["slot_manager_delta_t"].defaultvalue = 100.0
    par.glb_code_params_dict["slot_manager_t_min"].defaultvalue = -1000.0

    # Source Plane Geometric Mapping
    par.glb_code_params_dict["source_plane_center_x"].defaultvalue = 0.0
    par.glb_code_params_dict["source_plane_center_y"].defaultvalue = 0.0
    par.glb_code_params_dict["source_plane_center_z"].defaultvalue = 0.0
    par.glb_code_params_dict["source_plane_normal_x"].defaultvalue = 0.0
    par.glb_code_params_dict["source_plane_normal_y"].defaultvalue = 0.0
    par.glb_code_params_dict["source_plane_normal_z"].defaultvalue = 1.0
    par.glb_code_params_dict["source_r_max"].defaultvalue = 20.0
    par.glb_code_params_dict["source_r_min"].defaultvalue = 6.0
    par.glb_code_params_dict["source_up_vec_x"].defaultvalue = 0.0
    par.glb_code_params_dict["source_up_vec_y"].defaultvalue = 1.0
    par.glb_code_params_dict["source_up_vec_z"].defaultvalue = 0.0

    # Camera Window Geometric Mapping
    par.glb_code_params_dict["camera_pos_x"].defaultvalue = 51.0
    par.glb_code_params_dict["camera_pos_y"].defaultvalue = 0.0
    par.glb_code_params_dict["camera_pos_z"].defaultvalue = 10.2
    par.glb_code_params_dict["original_window_center_x"].defaultvalue = 50.0
    par.glb_code_params_dict["original_window_center_y"].defaultvalue = 0.0
    par.glb_code_params_dict["original_window_center_z"].defaultvalue = 10.0
    par.glb_code_params_dict["window_height"].defaultvalue = 1.0
    par.glb_code_params_dict["window_up_vec_x"].defaultvalue = 0.0
    par.glb_code_params_dict["window_up_vec_y"].defaultvalue = 0.0
    par.glb_code_params_dict["window_up_vec_z"].defaultvalue = 1.0
    par.glb_code_params_dict["window_width"].defaultvalue = 1.0
    par.glb_code_params_dict["window_tiles_width"].defaultvalue = 2
    par.glb_code_params_dict["window_tiles_height"].defaultvalue = 2

    # RKF45 Adaptive Control Tolerances
    par.glb_code_params_dict["numerical_initial_h"].defaultvalue = 0.1
    par.glb_code_params_dict["rkf45_absolute_error_tolerance"].defaultvalue = 1e-10
    par.glb_code_params_dict["rkf45_error_tolerance"].defaultvalue = 1e-10
    par.glb_code_params_dict["rkf45_h_max"].defaultvalue = 10.0
    par.glb_code_params_dict["rkf45_h_min"].defaultvalue = 1e-15

    # Execution Initial Conditions
    par.glb_code_params_dict["t_start"].defaultvalue = 1000.0
    par.glb_code_params_dict["scan_density"].defaultvalue = 500

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

    # OpenMP 
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

    BHaH_defines_h.output_BHaH_defines_h(
        project_dir=project_dir,
        enable_rfm_precompute=False,
        supplemental_defines_dict=cpu_macros,
    )

    compiler = "gcc"
    cflags = ["-fopenmp", "-O3", "-DDEBUG", "-Wno-stringop-truncation"]
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

    print(
        f"Finished! Now go into project/{project_name} and type `make` to build, then ./{exec_name} to run."
    )
    print(f"    Parameter file can be found in {project_name}.par\n")
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
