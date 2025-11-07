"""
Generates the complete C project for simulating photon geodesics (ray-tracing).

This script serves as the top-level orchestrator for the photon geodesic integrator
project. It follows the nrpy "gold standard" for an example script. Its primary
responsibilities are:

1.  **Symbolic Workflow**: Imports and instantiates the necessary symbolic "recipe"
    modules from `nrpy/equations/` to get the mathematical expressions for
    metrics, Christoffel symbols, geodesic equations of motion, and conserved
    quantities for massless particles.

2.  **C Code Generation**: Imports and calls the modular "builder" functions from
    the `.../geodesics/photon/` and `.../geodesics/` subdirectories. Each builder
    takes a symbolic recipe and generates a corresponding C function, registering
    it with nrpy's CFunction dictionary.

3.  **Project Assembly**: Orchestrates the creation of the final C project
    directory, generating parameter files, core BHaH header files, and the
    Makefile.

4.  **Post-processing**: Imports and calls utilities from `nrpy.helpers.kdtree_utils`
    to convert raw binary snapshot files from the mass integrator into
    query-optimized .kdtree.bin files, which are used to model the accretion disk.

This script can be run in several modes via command-line arguments:
- Default: Generates the C project.
- --run: Generates, compiles, runs the C code, and then runs post-processing.
- --post-process: Skips C code generation and runs only the post-processing step.

Author: Dalton J. Moone
"""

# ##############################################################################
# PART 0: IMPORTS AND PATH SETUP
# ##############################################################################

# Step 0.a: Add the nrpy root directory to the Python path.
import os
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
nrpy_root_dir = os.path.abspath(os.path.join(script_dir, "..", ".."))
if nrpy_root_dir not in sys.path:
    sys.path.insert(0, nrpy_root_dir)

# Step 0.b: Import standard Python modules
import argparse
import shutil
import subprocess

# Step 0.c: Import nrpy core modules
import nrpy.params as par
import nrpy.helpers.generic as gh
from nrpy.helpers import kdtree_utils

# Step 0.d: Import nrpy BHaH infrastructure modules
from nrpy.infrastructures.BHaH import (
    BHaH_defines_h,
    cmdline_input_and_parfiles,
    Makefile_helpers as Makefile,
    CodeParameters as CPs,
)

# Step 0.e: Import all necessary nrpy/equations modules (the "cookbooks")
from nrpy.equations.general_relativity import analytic_spacetimes as anasp
from nrpy.equations.general_relativity import geodesics as geo
from nrpy.equations.general_relativity.geodesic_diagnostics import conserved_quantities as geodiag

# Step 0.f: Import all C-generating "builder" functions.
# Note the new directory structure: generic builders are in `geodesics`,
# while photon-specific builders are in `geodesics/photon`.

# Generic/Shared Builders
from nrpy.infrastructures.BHaH.general_relativity.geodesics import g4dd_kerr_schild
from nrpy.infrastructures.BHaH.general_relativity.geodesics import g4dd_schwarzschild_cartesian
from nrpy.infrastructures.BHaH.general_relativity.geodesics import con_kerr_schild
from nrpy.infrastructures.BHaH.general_relativity.geodesics import con_schwarzschild_cartesian
from nrpy.infrastructures.BHaH.general_relativity.geodesics import g4dd_metric
from nrpy.infrastructures.BHaH.general_relativity.geodesics import connections
from nrpy.infrastructures.BHaH.general_relativity.geodesics import compare_filenames
from nrpy.infrastructures.BHaH.general_relativity.geodesics import load_kdtree_snapshot
from nrpy.infrastructures.BHaH.general_relativity.geodesics import unload_kdtree_snapshot
from nrpy.infrastructures.BHaH.general_relativity.geodesics import load_all_kdtree_snapshots
from nrpy.infrastructures.BHaH.general_relativity.geodesics import find_n_nearest_neighbors

# Photon-Specific Builders
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import calculate_p0_reverse
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import calculate_ode_rhs
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import check_conservation
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import set_initial_conditions_cartesian
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import radiative_transfer_engine
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import handle_source_plane_intersection
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import handle_disk_intersection
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import lagrange_interp_engine_generic
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import event_detection_manager
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import calculate_and_fill_blueprint_data_universal
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import rkf45_helpers_for_header
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import rkf45_update_and_control_helper
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import placeholder_interpolation_engine
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import time_slot_manager_helpers
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import batch_integrator_numerical
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import main


# ##############################################################################
# PART 1: PROJECT ASSEMBLY AND EXECUTION
# ##############################################################################

def register_custom_structures() -> None:
    """
    Generate C code for all custom structs and enums, then register them.

    This function defines the C `typedef`s for all custom data structures
    used by the photon geodesic project and registers this block of C code
    to be written into `BHaH_defines.h`.
    """
    metric_components = [f"g{i}{j}" for i in range(4) for j in range(i, 4)]
    metric_struct_str = "typedef struct { double " + ", ".join(metric_components) + "; } metric_struct;"

    connection_components = [f"Gamma4UDD{i}{j}{k}" for i in range(4) for j in range(4) for k in range(j, 4)]
    connections_struct_str = "typedef struct { double " + ", ".join(connection_components) + "; } connection_struct;"

    deriv_components = [f"g{i}{j}d{k}" for i in range(4) for j in range(i, 4) for k in range(4)]
    deriv_struct_str = "typedef struct { double " + ", ".join(deriv_components) + "; } g4DD_deriv_struct;"

    consolidated_structs_c_code = rf"""
// Core Metric & Connection Structs
{metric_struct_str}
{connections_struct_str}
{deriv_struct_str}

typedef enum {{ Schwarzschild, Kerr, Numerical, Schwarzschild_Standard }} Metric_t;
typedef struct {{ Metric_t type; }} metric_params;

// Event Detection and Plane Crossing Helpers
typedef struct {{ double n[3]; double d; }} plane_event_params;
typedef double (*event_function_t)(const double y[9], void *event_params);
static inline double plane_event_func(const double y[9], void *event_params) {{
    plane_event_params *params = (plane_event_params *)event_params;
    return y[1]*params->n[0] + y[2]*params->n[1] + y[3]*params->n[2] - params->d;
}}

// K-d Tree and Particle Data Structures
typedef struct {{
    int id; double pos[3]; double u[4]; double lambda_rest; float j_intrinsic;
}} __attribute__((packed)) MassiveParticle;
typedef struct {{
    int32_t* node_metadata; MassiveParticle* particle_data; uint64_t num_particles;
    uint64_t dimensions; void* original_mmap_ptr; size_t file_size;
}} CustomKDTree;
#define MAX_NEIGHBORS 10
typedef struct {{
    int indices[MAX_NEIGHBORS]; double sq_distances[MAX_NEIGHBORS]; int count; int n_wanted;
}} WinnersCircle;

// Batch Integration and Output Structs
typedef struct {{ int photon_id; double pos[4]; }} photon_request_t;
typedef struct {{ bool found; double lambda_event, t_event; double y_event[9]; }} event_data_struct;

typedef enum {{
    FAILURE_PT_TOO_BIG, FAILURE_RKF45_REJECTION_LIMIT, FAILURE_T_MAX_EXCEEDED,
    FAILURE_SLOT_MANAGER_ERROR, TERMINATION_TYPE_FAILURE,
    TERMINATION_TYPE_DISK, TERMINATION_TYPE_SOURCE_PLANE,
    TERMINATION_TYPE_CELESTIAL_SPHERE, ACTIVE,
}} termination_type_t;

typedef struct {{
    termination_type_t termination_type; double y_w, z_w; double stokes_I, lambda_observed;
    double y_s, z_s; double final_theta, final_phi; double L_w, t_w, L_s, t_s;
}} __attribute__((packed)) blueprint_data_t;
#define CACHE_LINE_SIZE 64
#define BUNDLE_CAPACITY 16384

typedef struct {{
    double y[9], y_p[9], y_p_p[9];
    double affine_param, affine_param_p, affine_param_p_p;
    double h;
    termination_type_t status;
    int rejection_retries;
    bool on_positive_side_of_window_prev, on_positive_side_of_source_prev;
    event_data_struct source_event_data, window_event_data;
    MassiveParticle nearest_neighbor;
    char _padding[CACHE_LINE_SIZE - (
        sizeof(double)*31 + sizeof(termination_type_t) + sizeof(int) +
        sizeof(bool)*2 + sizeof(event_data_struct)*2 + sizeof(MassiveParticle)
    ) % CACHE_LINE_SIZE];
}} __attribute__((aligned(CACHE_LINE_SIZE))) PhotonState;
"""
    BHaH_defines_h.register_BHaH_defines("after_general", consolidated_structs_c_code)


if __name__ == "__main__":
    # Step 1.A: Set up command-line argument parser.
    parser = argparse.ArgumentParser(
        description="Generate, compile, and run the photon geodesic integrator, or run post-processing."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="The parent directory where the C project and processed data will be generated.",
    )
    parser.add_argument(
        "--post-process",
        action="store_true",
        help="Run only the post-processing (k-d tree generation) step on existing snapshot files.",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Automatically compile, run the C code, and then run post-processing.",
    )
    args = parser.parse_args()

    # Step 1.B: Define core project names and paths.
    project_name = "photon_geodesic_integrator"
    exec_name = "photon_geodesic_integrator"
    project_dir = os.path.join(args.outdir, project_name)
    snapshot_output_dir = os.path.join(args.outdir, "processed_snapshots")

    # Step 1.A.1: Smart Dependency Check for --run
    # If --run is specified, check if the required processed snapshots exist.
    # If not, run the mass integrator first to generate them.
    if args.run and not os.path.isdir(snapshot_output_dir):
        print("----------------------------------------------------------", flush=True)
        print(f"Required data directory '{snapshot_output_dir}' not found.", flush=True)
        print("Automatically running the mass_geodesic_integrator to generate it...", flush=True)
        print("This may take several minutes.", flush=True)
        print("----------------------------------------------------------", flush=True)

        # Construct the path to the mass integrator script
        mass_integrator_script = os.path.join(script_dir, "mass_geodesic_integrator.py")
        if not os.path.exists(mass_integrator_script):
            print(f"Error: Could not find mass_geodesic_integrator.py at {mass_integrator_script}", file=sys.stderr)
            sys.exit(1)

        # Construct the command to run the mass integrator script
        # We pass the --outdir to ensure it outputs to the same parent directory
        command = [sys.executable, mass_integrator_script, "--run", "--outdir", args.outdir]

        try:
            # Run the mass integrator script as a subprocess, showing its output in real-time
            with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1) as proc:
                if proc.stdout:
                    for line in proc.stdout:
                        print(line, end='', flush=True)
            if proc.returncode != 0:
                raise subprocess.CalledProcessError(proc.returncode, command)

            print("----------------------------------------------------------", flush=True)
            print("Mass integrator finished. Resuming photon integrator...", flush=True)
            print("----------------------------------------------------------", flush=True)
        except (FileNotFoundError, subprocess.CalledProcessError) as e:
            print(f"Error: Failed to run mass_geodesic_integrator.py. Please run it manually.", file=sys.stderr)
            print(f"Details: {e}", file=sys.stderr)
            sys.exit(1)

    if args.post_process:
        # --- POST-PROCESSING-ONLY MODE ---
        # Assumes the mass integrator C code has already been run and produced binary snapshots.
        snapshot_input_dir = os.path.join(args.outdir, "mass_integrator", "output")
        kdtree_utils.run_kdtree_preprocessor_parallel(snapshot_input_dir, snapshot_output_dir)

    else:
        # --- FULL C-CODE GENERATION MODE ---
        # Step 1.C: Set up project directories and core nrpy parameters.
        shutil.rmtree(project_dir, ignore_errors=True)
        os.makedirs(project_dir, exist_ok=True)
        par.set_parval_from_str("Infrastructure", "BHaH")

        # Step 1.D: Register all CodeParameters for the photon integrator.
        param_defs = {
            "REAL": [
                ("M_scale", 1.0), ("a_spin", 0.0),
                ("camera_pos_x", 0.0), ("camera_pos_y", 0.0), ("camera_pos_z", 51.0),
                ("window_center_x", 0.0), ("window_center_y", 0.0), ("window_center_z", 50.0),
                ("window_up_vec_x", 0.0), ("window_up_vec_y", 1.0), ("window_up_vec_z", 0.0),
                ("window_size", 1.5),
                ("source_plane_normal_x", 0.0), ("source_plane_normal_y", 0.0), ("source_plane_normal_z", 1.0),
                ("source_plane_center_x", 0.0), ("source_plane_center_y", 0.0), ("source_plane_center_z", 0.0),
                ("source_up_vec_x", 0.0), ("source_up_vec_y", 1.0), ("source_up_vec_z", 0.0),
                ("source_r_min", 6.0), ("source_r_max", 25.0),
                ("t_start", 2000.0), ("t_integration_max", 10000.0), ("r_escape", 1500.0), ("p_t_max", 1000.0),
                ("mass_snapshot_every_t", 10.0), ("delta_r_max", 2.0),
                ("disk_bounds_x_min", -26.0), ("disk_bounds_x_max", 26.0),
                ("disk_bounds_y_min", -26.0), ("disk_bounds_y_max", 26.0),
                ("disk_bounds_z_min", -1.0), ("disk_bounds_z_max", 1.0),
                ("slot_manager_t_min", -100.0), ("slot_manager_delta_t", 0.1),
                ("numerical_initial_h", 1.0), ("rkf45_error_tolerance", 1e-8),
                ("rkf45_absolute_error_tolerance", 1e-8), ("rkf45_h_min", 1e-10),
                ("rkf45_h_max", 10.0), ("rkf45_safety_factor", 0.9),
                ("log_polar_r_min", 0.1),
            ],
            "int": [
                ("metric_choice", 0), ("scan_density", 512), ("rkf45_max_retries", 10),
                ("window_grid_type", 0), ("log_polar_num_r", 512), ("log_polar_num_phi", 1024),
            ],
            "bool": [
                ("use_numerical_pipeline", True), ("perform_conservation_check", False), ("debug_mode", False),
            ],
        }
        for c_type, params_list in param_defs.items():
            for name, default in params_list:
                par.register_CodeParameter(
                    c_type, __name__, name, default, commondata=True, add_to_parfile=True
                )

        # Step 1.E: Execute the symbolic workflow by instantiating the equation modules.
        print(" -> Instantiating symbolic equation modules for massless particles...")
        kerr_metric = anasp.Analytic_Spacetimes["KerrSchild"]
        schw_metric = anasp.Analytic_Spacetimes["Schwarzschild_Cartesian"]
        kerr_geo_eqs = geo.Geodesic_Equations["KerrSchild_massless"]
        schw_geo_eqs = geo.Geodesic_Equations["Schwarzschild_Cartesian_massless"]
        kerr_diags = geodiag.Geodesic_Diagnostics["KerrSchild_massless"]
        schw_diags = geodiag.Geodesic_Diagnostics["Schwarzschild_Cartesian_massless"]
        print("    ... All symbolic expressions generated and cached.")

        # Step 1.F: Generate C code for parameter handling.
        print(" -> Generating C code for parameter handling...")
        CPs.write_CodeParameters_h_files(project_dir=project_dir)
        CPs.register_CFunctions_params_commondata_struct_set_to_default()
        cmdline_input_and_parfiles.generate_default_parfile(
            project_dir=project_dir, project_name=project_name
        )
        cmdline_inputs_list = [name for _, params in param_defs.items() for name, _ in params]
        cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
            project_name=project_name, cmdline_inputs=cmdline_inputs_list
        )
        print("    ... Parameter handling code generated.")

        # Step 1.G: Call all C-generating Python "builder" functions to register C code.
        print(" -> Registering C functions for project...")
        # Generic Builders
        # Pass the symbolic expressions computed by the GeodesicEquations class
        # to the appropriate C-code builder functions.
        g4dd_kerr_schild.g4dd_kerr_schild(kerr_geo_eqs.g4DD)
        g4dd_schwarzschild_cartesian.g4dd_schwarzschild_cartesian(schw_geo_eqs.g4DD)
        con_kerr_schild.con_kerr_schild(kerr_geo_eqs.Gamma4UDD)
        con_schwarzschild_cartesian.con_schwarzschild_cartesian(schw_geo_eqs.Gamma4UDD)
        g4dd_metric.g4dd_metric()
        connections.connections()
        compare_filenames.compare_filenames()
        load_kdtree_snapshot.load_kdtree_snapshot()
        unload_kdtree_snapshot.unload_kdtree_snapshot()
        load_all_kdtree_snapshots.load_all_kdtree_snapshots()
        find_n_nearest_neighbors.find_n_nearest_neighbors()
        # Photon-Specific Builders
        calculate_p0_reverse.calculate_p0_reverse(kerr_geo_eqs.p0_expr)
        calculate_ode_rhs.calculate_ode_rhs(kerr_geo_eqs.geodesic_rhs)
        check_conservation.check_conservation(
            kerr_expressions=[kerr_diags.E_expr] + kerr_diags.L_exprs + [kerr_diags.Q_expr],
            schw_expressions=[schw_diags.E_expr] + schw_diags.L_exprs + [schw_diags.Q_expr],
        )
        set_initial_conditions_cartesian.set_initial_conditions_cartesian()
        radiative_transfer_engine.radiative_transfer_engine()
        handle_source_plane_intersection.handle_source_plane_intersection()
        handle_disk_intersection.handle_disk_intersection()
        lagrange_interp_engine_generic.lagrange_interp_engine_generic()
        event_detection_manager.event_detection_manager()
        calculate_and_fill_blueprint_data_universal.calculate_and_fill_blueprint_data_universal()
        rkf45_helpers_for_header.rkf45_helpers_for_header()
        rkf45_update_and_control_helper.rkf45_update_and_control_helper()
        placeholder_interpolation_engine.placeholder_interpolation_engine()
        time_slot_manager_helpers.time_slot_manager_helpers()
        batch_integrator_numerical.batch_integrator_numerical()
        main.main()
        print("    ... All C functions registered.")

        # Step 1.H: Assemble the final C project on disk.
        print(" -> Assembling final C project directory...")
        register_custom_structures()
        BHaH_defines_h.output_BHaH_defines_h(project_dir=project_dir, enable_rfm_precompute=False)
        gh.copy_files(
            package="nrpy.helpers",
            filenames_list=["simd_intrinsics.h"],
            project_dir=project_dir,
            subdirectory="intrinsics",
        )
        Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
            project_dir=project_dir,
            project_name=project_name,
            exec_or_library_name=exec_name,
            addl_CFLAGS=["-Wall -Wextra -g -fopenmp"],
            addl_libraries=["-lm -fopenmp"], # No GSL needed for custom RKF45
        )
        print("    ... Project assembly complete.")

        # Step 1.I: Print a final confirmation message.
        print(f"\nFinished! A C project has been generated in '{project_dir}/'")

        # Step 1.J: Optionally compile, run, and post-process the C code.
        if args.run:
            print("\n--- AUTO-COMPILE, RUN, AND POST-PROCESS ---")
            original_directory = os.getcwd()
            c_code_succeeded = False
            try:
                # 1. Compile the C code
                os.chdir(project_dir)
                print(f"[{project_name}] Compiling C code...")
                subprocess.run(["make", "-j"], check=True, capture_output=True, text=True)
                print(f"[{project_name}] Compilation successful.")

                # 2. Run the C executable
                print(f"[{project_name}] Running executable './{exec_name}'...")
                subprocess.run([f"./{exec_name}"], check=True)
                print(f"[{project_name}] C code execution finished.")
                c_code_succeeded = True

            except FileNotFoundError as e:
                print(f"Error: Required command not found. Is 'make' or a C compiler installed? Details: {e}")
            except subprocess.CalledProcessError as e:
                print(f"Error: An error occurred during compile or run.")
                print(f"STDOUT:\n{e.stdout}")
                print(f"STDERR:\n{e.stderr}")
            finally:
                os.chdir(original_directory)

        else:
            # If --run is not specified, print manual instructions.
            print(f"To build, navigate to '{project_dir}/' and type 'make'.")
            print(f"To run, type './{exec_name}'.")
