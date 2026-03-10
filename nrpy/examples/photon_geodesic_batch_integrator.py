"""
Orchestrator for the Split-Pipeline Photon Geodesic Integrator targeting the RTX 3080.

This script generates a high-performance C project for simulating photon trajectories 
in curved spacetimes. It utilizes a modular, SIMD-batch pipeline to circumvent the 
255-register hardware limit of the Ampere architecture by persisting intermediate 
tensors (metrics, connections, and derivatives) in Global VRAM scratchpads. 

The generated code employs a Structure of Arrays (SoA) memory layout and an 
adaptive RKF45 integration scheme managed by a lock-free TimeSlotManager system.

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
import nrpy.c_function as cfc
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

# Establish the script directory for reliable relative path resolution
script_dir = os.path.dirname(os.path.abspath(__file__))

# ##############################################################################
# PART 1: MAIN CONFIGURATION
# ##############################################################################

if __name__ == "__main__":
    # Step 1: Configure command-line arguments for the generation pipeline
    parser = argparse.ArgumentParser(
        description="Generate the Split-Pipeline Photon Geodesic Integrator (Kerr-Schild)."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="The parent directory where the generated C project will reside.",
    )
    parser.add_argument("--parallelization", type=str, default="cuda", choices=["cuda", "openmp"])
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


    # Instruct NRPy to use the BHaH infrastructure for macro expansions and SoA layouts
    par.set_parval_from_str("Infrastructure", "BHaH")
    # Set parallelization
    par.set_parval_from_str("parallelization", args.parallelization)

    # Step 4: Acquire Symbolic Physics Expressions
    print(f" -> Acquiring symbolic data for {GEO_KEY}...")
    metric_data = anasp.Analytic_Spacetimes[SPACETIME]
    geodesic_data = geo.Geodesic_Equations[GEO_KEY]

    # ##########################################################################
    # Step 5: Execute Modules and Register C Functions (Split-Pipeline Order)
    # ##########################################################################
    print(" -> Registering C functions and local CodeParameters...")

    # --- Initialization Phase Kernels ---
    set_initial_conditions_kernel.set_initial_conditions_kernel(SPACETIME)
    p0_reverse_kernel.p0_reverse_kernel(geodesic_data.p0_photon)

    # --- Fundamental Tensor Calculations (VRAM Persisted) ---
    g4DD_metric.g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)
    connections.connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)
    conserved_quantities.conserved_quantities(SPACETIME, PARTICLE)

    # --- Core Pipeline Kernels (The RKF45 Modular Loop) ---
    interpolation_kernel.interpolation_kernel(SPACETIME)
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(geodesic_data.geodesic_rhs, geodesic_data.xx)
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
    batch_integrator_numerical.batch_integrator_numerical(SPACETIME)
    main.main(SPACETIME)


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
        f"g4DD_metric_{SPACETIME}", 
        f"connections_{SPACETIME}"     
    ]:
        cfc.CFunction_dict.pop(internal_func, None)

    # ##########################################################################
    # Step 5.5: OVERRIDE DEFAULT CODE PARAMETERS
    # ##########################################################################

    # The RTX 3080 (sm_86 architecture) performs optimally with block sizes of 128 or 256 
    # for register-heavy kernels, rather than the infrastructure's default of 32.
    if "DEVICE_THREAD_MACROS" not in par.glb_extras_dict:
        par.glb_extras_dict["DEVICE_THREAD_MACROS"] = {}
    
    par.glb_extras_dict["DEVICE_THREAD_MACROS"].update({
        "BHAH_THREADS_IN_X_DIR_DEFAULT": 256,
        "BHAH_THREADS_IN_Y_DIR_DEFAULT": 1,
        "BHAH_THREADS_IN_Z_DIR_DEFAULT": 1,
    })


    print(" -> Overriding desired CodeParameters before .par generation...")

    # Analytic Spacetime Parameters
    par.glb_code_params_dict["M_scale"].defaultvalue = 1.0
    par.glb_code_params_dict["a_spin"].defaultvalue = 0.9

    # Batch Integrator & Numerical Limits
    par.glb_code_params_dict["p_t_max"].defaultvalue = 1000.0
    par.glb_code_params_dict["perform_conservation_check"].defaultvalue = True
    par.glb_code_params_dict["r_escape"].defaultvalue = 150.0
    par.glb_code_params_dict["slot_manager_delta_t"].defaultvalue = 300.0
    par.glb_code_params_dict["slot_manager_t_min"].defaultvalue = -1000.0
    par.glb_code_params_dict["t_integration_max"].defaultvalue = 10000.0

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
    par.glb_code_params_dict["window_center_x"].defaultvalue = 50.0
    par.glb_code_params_dict["window_center_y"].defaultvalue = 0.0
    par.glb_code_params_dict["window_center_z"].defaultvalue = 10.0
    par.glb_code_params_dict["window_height"].defaultvalue = 1.0
    par.glb_code_params_dict["window_up_vec_x"].defaultvalue = 0.0
    par.glb_code_params_dict["window_up_vec_y"].defaultvalue = 0.0
    par.glb_code_params_dict["window_up_vec_z"].defaultvalue = 1.0
    par.glb_code_params_dict["window_width"].defaultvalue = 1.0

    # RKF45 Adaptive Control Tolerances
    par.glb_code_params_dict["numerical_initial_h"].defaultvalue = 0.1
    par.glb_code_params_dict["rkf45_absolute_error_tolerance"].defaultvalue = 1e-12
    par.glb_code_params_dict["rkf45_error_tolerance"].defaultvalue = 1e-12
    par.glb_code_params_dict["rkf45_h_max"].defaultvalue = 10.0
    par.glb_code_params_dict["rkf45_h_min"].defaultvalue = 1e-15

    # Execution Initial Conditions
    par.glb_code_params_dict["scan_density"].defaultvalue = 500
    par.glb_code_params_dict["t_start"].defaultvalue = 1000.0

    # Step 6: Generate C Code for Parameter Handling
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
    # Step 7: Assemble Final C Project Infrastructure
    # ##########################################################################
    print(" -> Assembling C project on disk...")
    
    parallelization = args.parallelization

    if parallelization == "cuda":
        cuda_macros = {
            "BHAH_HD_FUNC": "#define BHAH_HD_FUNC __device__\n",
            "BHAH_HD_INLINE": "#define BHAH_HD_INLINE __device__ __inline__\n",
            "BHAH_WARP_ATOMIC_ADD(ptr, val)": "#define BHAH_WARP_ATOMIC_ADD(ptr, val) atomicAdd(ptr, val)\n",
            "GLOBAL_COMMONDATA_EXTERN": "extern __constant__ commondata_struct d_commondata;\n"
        }
        BHaH_defines_h.output_BHaH_defines_h(
            project_dir=project_dir, 
            enable_rfm_precompute=False,
            supplemental_defines_dict=cuda_macros
        )
        BHaH_device_defines_h.output_device_headers(project_dir=project_dir)
        
        print(" -> Copying hardware intrinsics...")
        gh.copy_files(
            package="nrpy.helpers",
            filenames_list=["cuda_intrinsics.h"],
            project_dir=project_dir,
            subdirectory="."
        )
        
        compiler = "nvcc"
        cflags = ["-lcudart", "-DUSE_GPU", "-rdc=true", "-DDEBUG"]
        libs = ["-lm", "-lcudart"]
        ext = "cu"
    else:
# OpenMP / CPU Fallback
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
            "BHAH_WARP_ATOMIC_ADD(ptr, val)": "#define BHAH_WARP_ATOMIC_ADD(ptr, val) _Pragma(\"omp atomic\") *(ptr) += (val)\n",
            "GLOBAL_COMMONDATA_EXTERN": "// CPU passes commondata by reference, no global needed.\n",
            "BHAH_DEVICE_SYNC()": "#define BHAH_DEVICE_SYNC() do {} while(0)\n"
        }

        BHaH_defines_h.output_BHaH_defines_h(
            project_dir=project_dir, 
            enable_rfm_precompute=False,
            supplemental_defines_dict=cpu_macros
        )
        
        compiler = "gcc"
        cflags = ["-fopenmp", "-O3", "-DDEBUG"]
        libs = ["-lm"]
        ext = "c"

    print(" -> Generating Makefile ")
    
    # Determine the correct optimization option string
    opt_option = "nvcc" if parallelization == "cuda" else "fast"
    
    Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=exec_name,
        compiler_opt_option=opt_option, 
        addl_CFLAGS=cflags,
        addl_libraries=libs,
        CC=compiler,
        src_code_file_ext=ext
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
    exec_path = os.path.join(project_dir, exec_name)
    try:
        subprocess.run([exec_path], cwd=project_dir, check=True)
        print("Ray-tracing complete. Binary blueprint generated.")
    except subprocess.CalledProcessError:
        print("C executable failed. Exiting pipeline.")
        sys.exit(1)

    print("\n--- PHASE 3: Generating Visualizations ---")
    vis_dir = os.path.abspath(os.path.join(script_dir, "..", "helpers", "geodesic_visualizations"))
    if not os.path.exists(vis_dir):
        print(f"ERROR: Visualization directory not found at {vis_dir}")
        sys.exit(1)

    sys.path.append(vis_dir)

    try:
        import config_and_types as cfg
        import render_lensed_image as rli

        blueprint_path = os.path.join(project_dir, "light_blueprint.bin")
        starmap_path = os.path.join(vis_dir, cfg.SPHERE_TEXTURE_FILE)
        output_image_path = os.path.join(project_dir, "lensed_output.png")
        STATIC_IMAGE_PIXEL_WIDTH = 500

        c_r_min = float(par.glb_code_params_dict["source_r_min"].defaultvalue)
        c_r_max = float(par.glb_code_params_dict["source_r_max"].defaultvalue)
        c_window_width = float(par.glb_code_params_dict["window_width"].defaultvalue)
        c_window_height = float(par.glb_code_params_dict["window_height"].defaultvalue)

        actual_window_width = max(c_window_width, c_window_height)
        source_physical_width = 2.0 * c_r_max

        disk_texture = rli.generate_source_disk_array(
            disk_physical_width=source_physical_width,
            disk_inner_radius=c_r_min,
            disk_outer_radius=c_r_max,
            colormap=cfg.COLORMAP,
        )

        rli.generate_static_lensed_image(
            output_filename=output_image_path,
            output_pixel_width=STATIC_IMAGE_PIXEL_WIDTH,
            source_image_width=source_physical_width,
            sphere_image=starmap_path,
            source_image=disk_texture,
            blueprint_filename=blueprint_path,
            window_width=actual_window_width,
        )
        print(f"\nPipeline Complete! Output image saved at: {output_image_path}")

    except ImportError as e:
        print(f"Failed to import visualization scripts: {e}")