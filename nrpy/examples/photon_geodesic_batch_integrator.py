"""
Generates the complete C project for simulating photon geodesics (ray-tracing) in curved spacetime.

This script serves as the top-level orchestrator and code-generation pipeline for a high-performance,
batch-processing photon geodesic integrator. Utilizing the NRPy+ (Numerical Relativity in Python)
infrastructure, it translates symbolic General Relativity equations (specifically the Kerr-Schild metric)
and geodesic equations into highly optimized, hardware-agnostic C code.

Architectural Overview of the Generated C Code:
    1. Memory Management: Utilizes a Structure of Arrays (SoA) layout mapped for unified
       memory to ensure optimal global memory bandwidth and cache alignment on parallel architectures.
    2. Initialization: Populates the Cartesian initial conditions for an arbitrary number of rays
       originating from a defined camera window.
    3. Staged Integration: Employs an embedded Runge-Kutta-Fehlberg 4(5) [RKF45] ODE solver.
       To handle highly variable trajectory lifespans without OpenMP race conditions or GPU
       thread divergence, it uses a lock-free temporal binning system (TimeSlotManager) alongside
       synchronized stream compaction (dense-to-sparse mapping).
    4. Event Processing: Detects intersections with bounding geometries (e.g., source accretion
       disks, camera windows, celestial escape spheres).
    5. Integrity Validation: Computes relative errors in constants of motion (Energy, Angular
       Momentum, Carter Constant) to track numerical drift.

Author: Dalton J. Moone
"""

# ##############################################################################
# PART 0: IMPORTS AND PATH SETUP
# ##############################################################################

import argparse
import os
import shutil
import subprocess
import sys

# NRPy core and helper modules for C code generation
import nrpy.helpers.generic as gh
import nrpy.params as par

# Physics/Math Generators (Symbolic definitions of spacetimes and geodesics)
from nrpy.equations.general_relativity.geodesics import analytic_spacetimes as anasp
from nrpy.equations.general_relativity.geodesics import geodesics as geo

# NRPy BlackHoles@Home (BHaH) infrastructure modules for C project management
from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.infrastructures.BHaH import BHaH_device_defines_h
from nrpy.infrastructures.BHaH import CodeParameters as CPs
from nrpy.infrastructures.BHaH import Makefile_helpers as Makefile
from nrpy.infrastructures.BHaH import cmdline_input_and_parfiles

# C-Code Builder Functions (Registers specific physical and numerical routines)
from nrpy.infrastructures.BHaH.general_relativity.geodesics import (
    connections,
    conserved_quantities,
    g4DD_metric,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
    batch_integrator_numerical,
    calculate_and_fill_blueprint_data_universal,
    calculate_ode_rhs,
    event_detection_manager,
    find_event_time_and_state,
    handle_source_plane_intersection,
    handle_window_plane_intersection,
    main,
    p0_reverse,
    placeholder_interpolation_engine,
    rkf45_helpers_for_header,
    rkf45_update_and_control_helper,
    set_initial_conditions_cartesian,
    time_slot_manager_helpers,
)

# Establish the script directory for reliable relative path resolution during visualization
script_dir = os.path.dirname(os.path.abspath(__file__))

# ##############################################################################
# PART 1: MAIN CONFIGURATION
# ##############################################################################

if __name__ == "__main__":
    # Step 1: Configure command-line arguments for the generation pipeline
    parser = argparse.ArgumentParser(
        description="Generate the Updated Photon Geodesic Integrator (Kerr-Schild)."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="The parent directory where the generated C project will reside.",
    )
    args = parser.parse_args()

    # Step 2: Define strict project constants and simulation targets
    project_name = "photon_geodesic_integrator"
    exec_name = "photon_geodesic_integrator"
    project_dir = os.path.abspath(os.path.join(args.outdir, project_name))

    # Define the physical regime: Target spacetime and particle classification
    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"

    # Step 3: Initialize the project directory and select the infrastructure backend
    print(f"Initializing project: {project_name}")
    shutil.rmtree(project_dir, ignore_errors=True)
    os.makedirs(project_dir, exist_ok=True)

    # Instruct NRPy+ to use the BHaH infrastructure for macro expansions and SoA layouts
    par.set_parval_from_str("Infrastructure", "BHaH")
    # Explicitly set parallelization to CUDA to enable device-specific header generation
    par.set_parval_from_str("parallelization", "cuda")

    # Step 4: Acquire Symbolic Physics Expressions
    # Translates abstract GR tensor mathematics into evaluable SymPy expressions
    print(f" -> Acquiring symbolic data for {GEO_KEY}...")
    metric_data = anasp.Analytic_Spacetimes[SPACETIME]
    geodesic_data = geo.Geodesic_Equations[GEO_KEY]

    # Step 5: Execute Modules and Register C Functions
    # This phase populates the internal C-function registry with generated ASTs (Abstract Syntax Trees)
    print(" -> Registering C functions and local CodeParameters...")

    # Event detection mechanisms (root-finding for geometric boundaries)
    event_detection_manager.event_detection_manager()
    find_event_time_and_state.find_event_time_and_state()
    handle_source_plane_intersection.handle_source_plane_intersection()
    handle_window_plane_intersection.handle_window_plane_intersection()
    calculate_and_fill_blueprint_data_universal.calculate_and_fill_blueprint_data_universal()


    # Fundamental tensor operations
    g4DD_metric.g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)
    connections.connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)
    calculate_ode_rhs.calculate_ode_rhs(geodesic_data.geodesic_rhs, geodesic_data.xx)
    p0_reverse.p0_reverse(geodesic_data.p0_photon)
    conserved_quantities.conserved_quantities(SPACETIME, PARTICLE)



    # Numerical routines: Adaptive ODE integration and temporal management
    rkf45_helpers_for_header.rkf45_helpers_for_header()
    rkf45_update_and_control_helper.rkf45_update_and_control_helper()
    time_slot_manager_helpers.time_slot_manager_helpers()
    placeholder_interpolation_engine.placeholder_interpolation_engine(SPACETIME)

    # Initial Data
    set_initial_conditions_cartesian.set_initial_conditions_cartesian(SPACETIME)


    # Top-level integration orchestrators
    batch_integrator_numerical.batch_integrator_numerical(SPACETIME)
    main.main(SPACETIME)

    # ##########################################################################
    # Step 5.5: OVERRIDE DEFAULT CODE PARAMETERS
    # ##########################################################################
    print(" -> Overriding desired CodeParameters before .par generation...")

    # --- Analytic Spacetime Parameters ---
    par.glb_code_params_dict["M_scale"].defaultvalue = 1.0  # Mass of the black hole
    par.glb_code_params_dict["a_spin"].defaultvalue = (
        0.9  # Dimensionless spin parameter
    )

    # --- Batch Integrator & Numerical Limits ---
    par.glb_code_params_dict["debug_mode"].defaultvalue = True
    par.glb_code_params_dict["p_t_max"].defaultvalue = (
        1000.0  # Maximum allowable temporal momentum before termination
    )
    par.glb_code_params_dict["perform_conservation_check"].defaultvalue = True
    par.glb_code_params_dict["r_escape"].defaultvalue = (
        150.0  # Radial boundary defining the celestial escape sphere
    )
    par.glb_code_params_dict["slot_manager_delta_t"].defaultvalue = 10.0
    par.glb_code_params_dict["slot_manager_t_min"].defaultvalue = -1000.0
    par.glb_code_params_dict["t_integration_max"].defaultvalue = 10000.0

    # --- Source Plane Geometric Mapping (Defines accretion disk bounds) ---
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

    # --- Camera Window Geometric Mapping (Defines virtual observer's frame) ---
    par.glb_code_params_dict["camera_pos_x"].defaultvalue = 51.0
    par.glb_code_params_dict["camera_pos_y"].defaultvalue = 0.0
    par.glb_code_params_dict["camera_pos_z"].defaultvalue = 10.2
    par.glb_code_params_dict["window_center_x"].defaultvalue = 50.0
    par.glb_code_params_dict["window_center_y"].defaultvalue = 0.0
    par.glb_code_params_dict["window_center_z"].defaultvalue = 10.0
    par.glb_code_params_dict["window_height"].defaultvalue = 1.0
    par.glb_code_params_dict["window_up_vec_x"].defaultvalue = 0.0
    par.glb_code_params_dict["window_up_vec_y"].defaultvalue = 0.0
    par.glb_code_params_dict["window_up_vec_z"].defaultvalue = 1.0
    par.glb_code_params_dict["window_width"].defaultvalue = 1.0

    # --- RKF45 Adaptive Control Tolerances (Governs step sizing and local error) ---
    par.glb_code_params_dict["numerical_initial_h"].defaultvalue = 0.1
    par.glb_code_params_dict["rkf45_absolute_error_tolerance"].defaultvalue = 1e-11
    par.glb_code_params_dict["rkf45_error_tolerance"].defaultvalue = 1e-11
    par.glb_code_params_dict["rkf45_h_max"].defaultvalue = 10.0
    par.glb_code_params_dict["rkf45_h_min"].defaultvalue = 1e-10
    par.glb_code_params_dict["rkf45_max_retries"].defaultvalue = 10
    par.glb_code_params_dict["rkf45_safety_factor"].defaultvalue = 0.9

    # --- Execution Initial Conditions ---
    par.glb_code_params_dict["scan_density"].defaultvalue = (
        700  # Resolution of the ray bundle
    )
    par.glb_code_params_dict["t_start"].defaultvalue = 500.0

    # Step 6: Generate C Code for Parameter Handling
    # Translates Python definitions into C header structs and parser implementations
    print(" -> Generating parameter handling code...")
    CPs.write_CodeParameters_h_files(project_dir=project_dir, set_commondata_only=True)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()
    cmdline_input_and_parfiles.generate_default_parfile(
        project_dir=project_dir, project_name=project_name
    )
    cmdline_inputs_list = list(par.glb_code_params_dict.keys())
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
        project_name=project_name, cmdline_inputs=cmdline_inputs_list
    )
# ##########################################################################
    # Step 7: Assemble Final C Project Infrastructure (CUDA Native Version)
    # ##########################################################################
    print(" -> Assembling C project on disk...")
    
    # Python: Generate the core BHaH defines header
    BHaH_defines_h.output_BHaH_defines_h(
        project_dir=project_dir, enable_rfm_precompute=False
    )
    
    
    # --- Inject Native CUDA Cross-Compilation & Warp-Aggregated Macros ---
    cuda_macros = r"""
    #ifdef __CUDACC__
        #define BHAH_HD_INLINE __device__ __inline__
        #define BHAH_HD_FUNC __device__
        #define BHAH_WARP_ATOMIC_ADD(ptr, val) atomicAdd(ptr, val)
    #else
        #define BHAH_HD_INLINE static inline
        #define BHAH_HD_FUNC static
        #define BHAH_WARP_ATOMIC_ADD(ptr, val) (*(ptr) += (val))
    #endif
    """
    with open(os.path.join(project_dir, "BHaH_defines.h"), "a") as f:
        f.write(cuda_macros)

    # Generate device-specific headers ONLY after all physics modules are registered
    BHaH_device_defines_h.output_device_headers(project_dir=project_dir)

    print(" -> Copying hardware intrinsics to project directory...")
    # Python: Copy cuda_intrinsics.h from nrpy/helpers/ to the project root.
    # Python: This resolves the fatal error: cuda_intrinsics.h: No such file or directory.
    gh.copy_files(
        package="nrpy.helpers",
        filenames_list=["cuda_intrinsics.h"],
        project_dir=project_dir,
        subdirectory="."  # Placing in root so #include "cuda_intrinsics.h" works
    )


    print(" -> Generating Makefile for RTX 3080 (sm_86)...")
    # Python: Orchestrate Makefile generation using the dedicated nvcc compiler profile.
    Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=exec_name,
        compiler_opt_option="nvcc",
        addl_CFLAGS=[
            "-lcudart",
            "-DUSE_GPU",
            "-rdc=true" ,
            "-DDEBUG"
        ],
        addl_libraries=["-lm", "-lcudart"],
        CC="nvcc",
        src_code_file_ext="cu"
    )

    # ##########################################################################
    # PART 2: PIPELINE EXECUTION (COMPILE, RUN, VISUALIZE)
    # ##########################################################################
    print("\n--- PHASE 1: Compiling C Code ---")
    try:
        subprocess.run(["make"], cwd=project_dir, check=True)
        print("Compilation successful.")
    except subprocess.CalledProcessError:
        print("Compilation failed. Exiting pipeline.")
        sys.exit(1)

    print("\n--- PHASE 2: Running Ray-Tracer ---")

    # Direct execution path for Linux (no .exe mapping needed)
    exec_path = os.path.join(project_dir, exec_name)

    try:
        subprocess.run([exec_path], cwd=project_dir, check=True)
        print("Ray-tracing complete. Binary blueprint generated.")
    except subprocess.CalledProcessError:
        print("C executable failed. Exiting pipeline.")
        sys.exit(1)

    print("\n--- PHASE 3: Generating Visualizations ---")

    # Dynamically map the visualization modules relative to the script's execution path
    vis_dir = os.path.abspath(
        os.path.join(script_dir, "..", "helpers", "geodesic_visualizations")
    )

    if not os.path.exists(vis_dir):
        print(f"ERROR: Visualization directory not found at {vis_dir}")
        sys.exit(1)

    # Append the resolved visualization directory to Python's execution path
    sys.path.append(vis_dir)

    try:
        # Import dynamically resolved visualization modules
        import config_and_types as cfg  # type: ignore
        import render_lensed_image as rli  # type: ignore

        # Establish standardized file path mappings for the renderer payload
        blueprint_path = os.path.join(project_dir, "light_blueprint.bin")
        starmap_path = os.path.join(vis_dir, cfg.SPHERE_TEXTURE_FILE)
        output_image_path = os.path.join(project_dir, "lensed_output.png")

        # --- DYNAMICALLY EXTRACT NRPY PARAMETERS ---
        c_m_scale = float(par.glb_code_params_dict["M_scale"].defaultvalue)
        c_r_min = float(par.glb_code_params_dict["source_r_min"].defaultvalue)
        c_r_max = float(par.glb_code_params_dict["source_r_max"].defaultvalue)
        c_window_width = float(par.glb_code_params_dict["window_width"].defaultvalue)
        c_window_height = float(par.glb_code_params_dict["window_height"].defaultvalue)

        # Calculate derived bounds for the rendering engine phase space
        actual_window_width = max(c_window_width, c_window_height)
        source_physical_width = 2.0 * c_r_max

        print(f"  -> Extracted M_scale: {c_m_scale}")
        print(f"  -> Extracted Disk Bounds: [{c_r_min}, {c_r_max}]")
        print(f"  -> Calculated Renderer Window FOV: {actual_window_width}")

        # Instantiate procedural disk texture in memory leveraging bounding physics parameters
        disk_texture = rli.generate_source_disk_array(
            disk_physical_width=source_physical_width,
            disk_inner_radius=c_r_min,
            disk_outer_radius=c_r_max,
            colormap=cfg.COLORMAP,
        )

        # Execute chunked rendering pass
        rli.generate_static_lensed_image(
            output_filename=output_image_path,
            output_pixel_width=cfg.STATIC_IMAGE_PIXEL_WIDTH,
            source_image_width=source_physical_width,
            sphere_image=starmap_path,
            source_image=disk_texture,
            blueprint_filename=blueprint_path,
            window_width=actual_window_width,
        )
        print(f"\nPipeline Complete! Output image saved at: {output_image_path}")

    except ImportError as e:
        print(f"Failed to import visualization scripts: {e}")
