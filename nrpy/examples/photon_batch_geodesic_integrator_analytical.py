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
    parser.add_argument(
        "--cuda",
        action="store_true",
        help="Enable CUDA support. Defaults to OpenMP if omitted.",
    )
    args = parser.parse_args()

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

    # Step 5.a: Register initialization kernels.
    set_initial_conditions_kernel.set_initial_conditions_kernel(SPACETIME)

    if geodesic_data.p0_photon is None:
        raise ValueError(f"p0_photon is None for {GEO_KEY}")
    p0_reverse_kernel.p0_reverse_kernel(geodesic_data.p0_photon)

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
        geodesic_data.geodesic_rhs, geodesic_data.xx
    )
    rkf45_stage_update.rkf45_stage_update()
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel()

    # Step 5.d: Register event-detection and boundary-intersection kernels.
    find_event_time_and_state.find_event_time_and_state()
    handle_source_plane_intersection.handle_source_plane_intersection()
    handle_window_plane_intersection.handle_window_plane_intersection()
    event_detection_manager_kernel.event_detection_manager_kernel()
    calculate_and_fill_blueprint_data_universal.calculate_and_fill_blueprint_data_universal()

    # Step 5.e: Register project-level orchestration helpers.
    time_slot_manager_helpers.time_slot_manager_helpers()
    batch_integrator_analytical.batch_integrator_analytical(SPACETIME)
    main.main(SPACETIME)

    # Step 5.f: Remove helper registrations emitted only through other kernels.
    # The event manager emits its local event helpers through prefunc, so keeping
    # standalone registrations would add unused source files and prototypes. The
    # metric and connection helpers are likewise inlined where they are consumed
    # in this generated project. Removing these entries prevents redundant emitted
    # wrappers in both OpenMP and CUDA builds, and also avoids duplicate CUDA
    # definitions during device linking.
    for internal_func in [
        "find_event_time_and_state",
        "handle_source_plane_intersection",
        "handle_window_plane_intersection",
        f"g4DD_metric_{SPACETIME}",
        f"connections_{SPACETIME}",
    ]:
        cfc.CFunction_dict.pop(internal_func, None)

    # Step 6: Override CodeParameter defaults before parfile generation.
    print(" -> Overriding desired CodeParameters before .par generation...")

    # Step 6.a: Set analytic spacetime defaults.
    par.glb_code_params_dict["M_scale"].defaultvalue = 1.0
    par.glb_code_params_dict["a_spin"].defaultvalue = 0.9

    # Step 6.b: Set batch-integrator and numerical-limit defaults.
    par.glb_code_params_dict["p_t_max"].defaultvalue = 1000.0
    par.glb_code_params_dict["perform_conservation_check"].defaultvalue = True
    par.glb_code_params_dict["r_escape"].defaultvalue = 100.0
    par.glb_code_params_dict["slot_manager_delta_t"].defaultvalue = 100.0
    par.glb_code_params_dict["slot_manager_t_min"].defaultvalue = -1000.0

    # Step 6.c: Set source-plane geometry defaults.
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

    # Step 6.d: Set camera window geometry defaults.
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

    # Step 6.e: Set window tiling defaults for the selected backend.
    if parallelization_mode == "cuda":
        par.glb_code_params_dict["window_tiles_width"].defaultvalue = 1
        par.glb_code_params_dict["window_tiles_height"].defaultvalue = 1
    else:
        par.glb_code_params_dict["window_tiles_width"].defaultvalue = 2
        par.glb_code_params_dict["window_tiles_height"].defaultvalue = 2

    # Step 6.f: Set RKF45 controller defaults.
    par.glb_code_params_dict["numerical_initial_h"].defaultvalue = 0.1
    par.glb_code_params_dict["rkf45_absolute_error_tolerance"].defaultvalue = 1e-10
    par.glb_code_params_dict["rkf45_error_tolerance"].defaultvalue = 1e-10
    par.glb_code_params_dict["rkf45_h_max"].defaultvalue = 10.0
    par.glb_code_params_dict["rkf45_h_min"].defaultvalue = 1e-15

    # Step 6.g: Set initial-condition defaults.
    par.glb_code_params_dict["t_start"].defaultvalue = 1000.0

    # Step 6.h: Set default scan density for the selected backend.
    if parallelization_mode == "cuda":
        par.glb_code_params_dict["scan_density"].defaultvalue = 1000
    else:
        par.glb_code_params_dict["scan_density"].defaultvalue = 500

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
            "source_r_min",
            "source_r_max",
            "window_width",
            "window_height",
            "window_tiles_width",
            "window_tiles_height",
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
    # Step 8: Copy visualization helpers and print usage instructions.
    vis_dir = os.path.join("nrpy", "examples", "geodesic_visualizations")
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

    # Step 8.a: Build visualization command arguments from CodeParameter defaults.
    c_r_min = float(par.glb_code_params_dict["source_r_min"].defaultvalue)
    c_r_max = float(par.glb_code_params_dict["source_r_max"].defaultvalue)
    c_window_width = float(par.glb_code_params_dict["window_width"].defaultvalue)
    c_window_height = float(par.glb_code_params_dict["window_height"].defaultvalue)
    c_tiles_width = int(par.glb_code_params_dict["window_tiles_width"].defaultvalue)
    c_tiles_height = int(par.glb_code_params_dict["window_tiles_height"].defaultvalue)
    c_pixel_width = 600
    parfile_path = os.path.join(project_dir, f"{project_name}.par")

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

    blueprint_command = (
        f"python3 blueprint_analysis.py "
        f"--window_tiles_width {c_tiles_width} "
        f"--window_tiles_height {c_tiles_height} "
        f"--window_width {c_window_width} "
        f"--window_height {c_window_height}"
    )

    print("    To run the blueprint diagnostic:")
    print(f"    {blueprint_command}\n")
