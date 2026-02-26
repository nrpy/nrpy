"""
Generates the complete C project for simulating photon geodesics (ray-tracing).

This script serves as the top-level orchestrator for the updated photon geodesic 
integrator project. It generates a standalone C application that evolves photon 
trajectories in the Kerr-Schild spacetime.

Author: Dalton J. Moone (Updated)
"""

# ##############################################################################
# PART 0: IMPORTS AND PATH SETUP
# ##############################################################################

import os
import sys
import shutil
import argparse
import subprocess

# Step 0.a: Add the nrpy root directory to the Python path.
script_dir = os.path.dirname(os.path.abspath(__file__))

# Import nrpy core modules
import nrpy.params as par
import nrpy.helpers.generic as gh

# Import nrpy BHaH infrastructure modules
from nrpy.infrastructures.BHaH import (
    BHaH_defines_h,
    cmdline_input_and_parfiles,
    Makefile_helpers as Makefile,
    CodeParameters as CPs,
)

# Import Physics/Math Generators (Symbolic)
from nrpy.equations.general_relativity.geodesics import analytic_spacetimes as anasp
from nrpy.equations.general_relativity.geodesics import geodesics as geo

# Import C-Code Builder Functions
from nrpy.infrastructures.BHaH.general_relativity.geodesics import g4DD_metric
from nrpy.infrastructures.BHaH.general_relativity.geodesics import connections
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import calculate_ode_rhs
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import p0_reverse
from nrpy.infrastructures.BHaH.general_relativity.geodesics import conserved_quantities

from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import set_initial_conditions_cartesian
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import handle_source_plane_intersection
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import handle_window_plane_intersection
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import event_detection_manager
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import find_event_time_and_state
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import calculate_and_fill_blueprint_data_universal

from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import rkf45_helpers_for_header
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import rkf45_update_and_control_helper
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import time_slot_manager_helpers
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import placeholder_interpolation_engine

from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import batch_integrator_numerical
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import main

# ##############################################################################
# PART 1: MAIN CONFIGURATION
# ##############################################################################

if __name__ == "__main__":
    # Step 1: Set up arguments
    parser = argparse.ArgumentParser(
        description="Generate the Updated Photon Geodesic Integrator (Kerr-Schild)."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="The parent directory where the C project will be generated.",
    )
    args = parser.parse_args()

    # Step 2: Define Project Constants
    project_name = "photon_geodesic_integrator"
    exec_name = "photon_geodesic_integrator"
    project_dir = os.path.abspath(os.path.join(args.outdir, project_name))
    
    # Configuration
    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"

    # Step 3: Setup Directory and Core Infrastructure
    print(f"Initializing project: {project_name}")
    shutil.rmtree(project_dir, ignore_errors=True)
    os.makedirs(project_dir, exist_ok=True)
    
    par.set_parval_from_str("Infrastructure", "BHaH")

    # Step 4: Acquire Symbolic Physics Expressions
    print(f" -> Acquiring symbolic data for {GEO_KEY}...")
    metric_data = anasp.Analytic_Spacetimes[SPACETIME]
    geodesic_data = geo.Geodesic_Equations[GEO_KEY]

    # Step 5: Execute Modules and Register C Functions
    print(" -> Registering C functions and local CodeParameters...")
    g4DD_metric.g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)
    connections.connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)
    calculate_ode_rhs.calculate_ode_rhs(geodesic_data.geodesic_rhs, geodesic_data.xx)
    p0_reverse.p0_reverse(geodesic_data.p0_photon)
    conserved_quantities.conserved_quantities(SPACETIME, PARTICLE)
    set_initial_conditions_cartesian.set_initial_conditions_cartesian(SPACETIME)
    event_detection_manager.event_detection_manager()
    find_event_time_and_state.find_event_time_and_state()
    handle_source_plane_intersection.handle_source_plane_intersection()
    handle_window_plane_intersection.handle_window_plane_intersection()
    calculate_and_fill_blueprint_data_universal.calculate_and_fill_blueprint_data_universal()
    rkf45_helpers_for_header.rkf45_helpers_for_header(SPACETIME)
    rkf45_update_and_control_helper.rkf45_update_and_control_helper()
    time_slot_manager_helpers.time_slot_manager_helpers()
    placeholder_interpolation_engine.placeholder_interpolation_engine(SPACETIME)
    batch_integrator_numerical.batch_integrator_numerical(SPACETIME)
    main.main(SPACETIME)

    # ##########################################################################
    # Step 5.5: OVERRIDE DEFAULT CODE PARAMETERS 
    # ##########################################################################
    print(" -> Overriding desired CodeParameters before .par generation...")
    # Analytic Spacetimes
    par.glb_code_params_dict["M_scale"].defaultvalue = 1.0
    par.glb_code_params_dict["a_spin"].defaultvalue = 0.9

    # Batch Integrator Numerical
    par.glb_code_params_dict["debug_mode"].defaultvalue = True
    par.glb_code_params_dict["p_t_max"].defaultvalue = 1000.0
    par.glb_code_params_dict["perform_conservation_check"].defaultvalue = True
    par.glb_code_params_dict["r_escape"].defaultvalue = 150.0
    par.glb_code_params_dict["slot_manager_delta_t"].defaultvalue = 10.0
    par.glb_code_params_dict["slot_manager_t_min"].defaultvalue = -1000.0
    par.glb_code_params_dict["t_integration_max"].defaultvalue = 10000.0

    # Source Plane Intersection
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

    # Window Plane Intersection
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

    # RKF45 Update and Control Helper
    par.glb_code_params_dict["numerical_initial_h"].defaultvalue = 0.1
    par.glb_code_params_dict["rkf45_absolute_error_tolerance"].defaultvalue = 1e-9
    par.glb_code_params_dict["rkf45_error_tolerance"].defaultvalue = 1e-9
    par.glb_code_params_dict["rkf45_h_max"].defaultvalue = 10.0
    par.glb_code_params_dict["rkf45_h_min"].defaultvalue = 1e-10
    par.glb_code_params_dict["rkf45_max_retries"].defaultvalue = 10
    par.glb_code_params_dict["rkf45_safety_factor"].defaultvalue = 0.9

    # Set Initial Conditions Cartesian
    par.glb_code_params_dict["scan_density"].defaultvalue = 700
    par.glb_code_params_dict["t_start"].defaultvalue = 500.0

    # Step 6: Generate C Code for Parameters
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

    ##########################################################################
    # Step 7: Assemble Final C Project (Linux Optimized)
    # ########################################################################
    print(" -> Assembling C project on disk...")

    # Output BHaH_defines.h using the standard infrastructure
    BHaH_defines_h.output_BHaH_defines_h(
        project_dir=project_dir, 
        enable_rfm_precompute=False
    )

    # Copy SIMD headers required for high-performance physics
    gh.copy_files(
        package="nrpy.helpers",
        filenames_list=["simd_intrinsics.h"],
        project_dir=project_dir,
        subdirectory="intrinsics",
    )
    
    print(" -> Generating Makefile and Function Prototypes...")
    Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=exec_name,
        addl_CFLAGS=["-Wall -Wextra -g -fopenmp -O3 -march=native -Wno-stringop-truncation -Wno-unknown-pragmas"],
        addl_libraries=["-lm", "-fopenmp"], 
    )

    # ##########################################################################
    # PART 2: PIPELINE EXECUTION (COMPILE, RUN, VISUALIZE)
    # ##########################################################################
    print("\n--- PHASE 1: Compiling C Code ---")
    try:
        # Standard Linux multi-threaded build
        subprocess.run(["make", "-j"], cwd=project_dir, check=True)
        print("Compilation successful.")
    except subprocess.CalledProcessError:
        print("Compilation failed. Ensure 'gcc' and 'make' are installed. Exiting.")
        sys.exit(1)

    print("\n--- PHASE 2: Running Ray-Tracer ---")
    
    # Construct the standard Linux execution path for the local binary
    exec_path = os.path.join(".", exec_name)
    
    try:
        subprocess.run([exec_path], cwd=project_dir, check=True)
        print("Ray-tracing complete. Blueprint generated.")
    except subprocess.CalledProcessError as e:
        print(f"C executable failed with error: {e}. Exiting pipeline.")
        sys.exit(1)

    print("\n--- PHASE 3: Generating Visualizations ---")
    
    # Dynamically find the visualization directory relative to this script
    vis_dir = os.path.abspath(os.path.join(script_dir, "..", "helpers", "geodesic_visualizations"))
    
    if not os.path.exists(vis_dir):
        print(f"ERROR: Visualization directory not found at {vis_dir}")
        sys.exit(1)

    # Add the visualization directory to the Python path
    sys.path.append(vis_dir)
    
    try:
        import render_lensed_image as rli
        import config_and_types as cfg
        
        # Define Linux-style paths for the renderer
        blueprint_path = os.path.join(project_dir, "light_blueprint.bin")
        starmap_path = os.path.join(vis_dir, cfg.SPHERE_TEXTURE_FILE)
        output_image_path = os.path.join(project_dir, "lensed_output.png")

        # --- DYNAMICALLY EXTRACT NRPY PARAMETERS ---
        c_m_scale = float(par.glb_code_params_dict["M_scale"].defaultvalue)
        c_r_min = float(par.glb_code_params_dict["source_r_min"].defaultvalue)
        c_r_max = float(par.glb_code_params_dict["source_r_max"].defaultvalue)
        c_window_width = float(par.glb_code_params_dict["window_width"].defaultvalue)
        c_window_height = float(par.glb_code_params_dict["window_height"].defaultvalue)
        
        # Calculate derived physics for rendering
        actual_window_width = max(c_window_width, c_window_height)
        source_physical_width = 2.0 * c_r_max

        print(f"  -> Extracted M_scale: {c_m_scale}")
        print(f"  -> Extracted Disk Bounds: [{c_r_min}, {c_r_max}]")
        print(f"  -> Calculated Renderer Window FOV: {actual_window_width}")

        # Create procedural disk texture in memory using physics parameters
        disk_texture = rli.generate_source_disk_array(
            disk_physical_width=source_physical_width,
            disk_inner_radius=c_r_min,
            disk_outer_radius=c_r_max,
            colormap=cfg.COLORMAP
        )

        # Execute chunked rendering process
        rli.generate_static_lensed_image(
            output_filename=output_image_path,
            output_pixel_width=cfg.STATIC_IMAGE_PIXEL_WIDTH,
            source_image_width=source_physical_width,
            sphere_image=starmap_path,
            source_image=disk_texture,
            blueprint_filename=blueprint_path,
            window_width=actual_window_width
        )
        print(f"\nPipeline Complete! Output image saved at: {output_image_path}")
        
    except ImportError as e:
        print(f"Failed to import visualization scripts: {e}")