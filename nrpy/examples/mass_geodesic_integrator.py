"""
Generates the complete C project for simulating massive particle geodesics.

This script serves as the top-level orchestrator for the mass geodesic integrator
project. It follows the nrpy "gold standard" for an example script. Its primary
responsibilities are:

1.  **Symbolic Workflow**: Imports and instantiates the necessary symbolic "recipe"
    modules from `nrpy/equations/` to get the mathematical expressions for
    metrics, Christoffel symbols, geodesic equations of motion, and conserved
    quantities.

2.  **C Code Generation**: Imports and calls the modular "builder" functions from
    the `.../geodesics/` subdirectory. Each builder takes a symbolic recipe
    and generates a corresponding C function, registering it with nrpy's
    CFunction dictionary.

3.  **Project Assembly**: Orchestrates the creation of the final C project
    directory, generating parameter files, core BHaH header files, and the
    Makefile.

4.  **Post-processing**: Imports and calls utilities from `nrpy.helpers.kdtree_utils`
    to convert raw binary snapshot files into query-optimized .kdtree.bin files.

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
#           This allows for clean imports of all nrpy modules.
import os
import sys

# Get the absolute path of the directory containing this script.
script_dir = os.path.dirname(os.path.abspath(__file__))
# Navigate up two levels to get the absolute path of the nrpy root directory.
# e.g., from /path/to/nrpy/nrpy/examples -> /path/to/nrpy
nrpy_root_dir = os.path.abspath(os.path.join(script_dir, "..", ".."))


# Add the nrpy root directory to the front of the Python path.
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
from nrpy.equations.general_relativity.InitialData_Cartesian import InitialData_Cartesian

# Step 0.f: Import all shared C-generating "builder" functions from the geodesics subdirectory
from nrpy.infrastructures.BHaH.general_relativity.geodesics import g4dd_kerr_schild as g4dd_ks
from nrpy.infrastructures.BHaH.general_relativity.geodesics import g4dd_schwarzschild_cartesian as g4dd_sc
from nrpy.infrastructures.BHaH.general_relativity.geodesics import con_kerr_schild as con_ks
from nrpy.infrastructures.BHaH.general_relativity.geodesics import con_schwarzschild_cartesian as con_sc
from nrpy.infrastructures.BHaH.general_relativity.geodesics import g4dd_metric
from nrpy.infrastructures.BHaH.general_relativity.geodesics import connections

# Step 0.g: Import all massive C-generating "builder" functions from the geodesics.massive subdirectory
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import calculate_ode_rhs_massive as calc_ode_rhs
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import check_conservation_massive as check_cons
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import ode_gsl_wrapper_massive as ode_gsl
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import calculate_ut_uphi_from_r as calc_ut_uphi
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import set_initial_conditions_massive as set_ic
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import integrate_single_particle as int_single
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import integrate_single_particle_DEBUG as int_debug
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import generate_disk_initial_conditions as gen_disk_ic
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import generate_spiral_galaxy_initial_conditions as gen_spiral_ic
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import generate_barred_flocculent_spiral_ic as gen_barred_ic
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import run_mass_integrator_production as run_prod
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive import main_c_function



# ##############################################################################
# PART 1: PROJECT ASSEMBLY AND EXECUTION
# ##############################################################################

def register_custom_structures_and_params() -> None:
    """
    Generate C code for all custom structs and enums, then register them.

    This function defines the C `typedef`s for all custom data structures
    and registers this block of C code to be written into `BHaH_defines.h`.
    """
    metric_components = [f"g{nu}{mu}" for nu in range(4) for mu in range(nu, 4)]
    metric_struct_str = (
        "typedef struct { double " + "; double ".join(metric_components) + "; } metric_struct;"
    )
    connection_components = [
        f"Gamma4UDD{i}{j}{k}" for i in range(4) for j in range(4) for k in range(j, 4)
    ]
    connections_struct_str = (
        "typedef struct { double "
        + "; double ".join(connection_components)
        + "; } connection_struct;"
    )
    other_structs = r"""
typedef enum { Schwarzschild, Kerr, Schwarzschild_Standard } Metric_t;
typedef struct { Metric_t type; } metric_params;
typedef struct {
    const commondata_struct *commondata;
    const params_struct *params;
    const metric_params *metric;
} gsl_params;
typedef struct {
    int id; double pos[4]; double u_spatial[3];
} particle_initial_state_t;
typedef struct {
    int id; double pos[3]; double u[4]; double lambda_rest; float j_intrinsic;
} __attribute__((packed)) mass_particle_state_t;
"""
    BHaH_defines_h.register_BHaH_defines(
        "data_structures",
        f"{metric_struct_str}\n{connections_struct_str}\n{other_structs}",
    )


if __name__ == "__main__":
    # Step 1.A: Set up command-line argument parser.
    parser = argparse.ArgumentParser(
        description="Generate, compile, and run the mass geodesic integrator, or run post-processing."
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
    project_name = "mass_integrator"
    exec_name = "mass_integrator"
    project_dir = os.path.join(args.outdir, project_name)
    snapshot_output_dir = os.path.join(args.outdir, "processed_snapshots")

    if args.post_process:
        # --- POST-PROCESSING-ONLY MODE ---
        snapshot_input_dir = os.path.join(project_dir, "output")
        kdtree_utils.run_kdtree_preprocessor_parallel(snapshot_input_dir, snapshot_output_dir)

    else:
        # --- FULL C-CODE GENERATION MODE ---
        # Step 1.C: Set up project directories and core nrpy parameters.
        shutil.rmtree(project_dir, ignore_errors=True)
        os.makedirs(project_dir, exist_ok=True)
        par.set_parval_from_str("Infrastructure", "BHaH")

        # Step 1.D: Register all CodeParameters for the mass integrator.
        param_groups = {
            "REAL": {
                "M_scale": 1.0, "a_spin": 0.9, "t_max_integration": 2000.0,
                "r_escape": 1500.0, "ut_max": 1e3, "disk_lambda_rest_at_r_min": 656.3,
                "disk_r_min": 6.0, "disk_r_max": 25.0, "snapshot_every_t": 5.0,
                "t_final": 100.0, "bar_length": 5.0, "bar_aspect_ratio": 0.25,
                "bulge_radius": 1.5, "arm_particle_density": 0.3,
                "arm_clumpiness_factor": 8.0, "arm_clump_size": 0.5,
                "bar_density_factor": 2.0, "bulge_density_factor": 3.0,
                "spiral_galaxy_arm_tightness": 0.2,
            },
            "int": {"disk_num_r": 100, "disk_num_phi": 200, "spiral_galaxy_num_arms": 2},
            "bool": {
                "generate_kdtree_files": True, "perform_conservation_check": True,
                "run_in_debug_mode": False,
            },
            "char[100]": {"output_folder": "output", "initial_conditions_type": "KeplerianDisk"},
        }
        for c_type, params_dict in param_groups.items():
            for name, default in params_dict.items():
                par.register_CodeParameter(
                    c_type, __name__, name, default, commondata=True, add_to_parfile=True
                )

        # Step 1.E: Execute the symbolic workflow by instantiating the equation modules.
        print(" -> Instantiating symbolic equation modules...")
        kerr_metric = anasp.Analytic_Spacetimes["KerrSchild"]
        schw_metric = anasp.Analytic_Spacetimes["Schwarzschild_Cartesian"]
        kerr_geo_eqs = geo.Geodesic_Equations["KerrSchild_massive"]
        schw_geo_eqs = geo.Geodesic_Equations["Schwarzschild_Cartesian_massive"]
        kerr_diags = geodiag.Geodesic_Diagnostics["KerrSchild_massive"]
        schw_diags = geodiag.Geodesic_Diagnostics["Schwarzschild_Cartesian_massive"]
        circ_orbit_ID = InitialData_Cartesian(
            "MassiveParticle_StableCircularOrbit",
            M_scale=anasp.M_scale,
            a_spin=anasp.a_spin,
        )
        print("    ... All symbolic expressions generated and cached.")

        # Step 1.F: Generate C code for parameter handling.
        print(" -> Generating C code for parameter handling...")
        CPs.write_CodeParameters_h_files(project_dir=project_dir)
        CPs.register_CFunctions_params_commondata_struct_set_to_default()
        cmdline_input_and_parfiles.generate_default_parfile(
            project_dir=project_dir, project_name=project_name
        )
        cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
            project_name=project_name,
            cmdline_inputs=list(param_groups["REAL"].keys())
            + list(param_groups["int"].keys())
            + list(param_groups["bool"].keys())
            + list(param_groups["char[100]"].keys()),
        )
        print("    ... Parameter handling code generated.")

        # Step 1.G: Call all C-generating Python "builder" functions to register C code.
        print(" -> Registering C functions for project...")
        g4dd_ks.g4dd_kerr_schild(kerr_geo_eqs.g4DD)
        g4dd_sc.g4dd_schwarzschild_cartesian(schw_geo_eqs.g4DD)
        con_ks.con_kerr_schild(kerr_geo_eqs.Gamma4UDD)
        con_sc.con_schwarzschild_cartesian(schw_geo_eqs.Gamma4UDD)
        g4dd_metric.g4dd_metric()
        connections.connections()
        calc_ode_rhs.calculate_ode_rhs_massive(kerr_geo_eqs.geodesic_rhs)
        check_cons.check_conservation_massive(
            kerr_E_expr=kerr_diags.E_expr,
            kerr_L_exprs=kerr_diags.L_exprs,
            kerr_Q_expr=kerr_diags.Q_expr,
            schw_E_expr=schw_diags.E_expr,
            schw_L_exprs=schw_diags.L_exprs,
            schw_Q_expr=schw_diags.Q_expr,
        )
        ode_gsl.ode_gsl_wrapper_massive()
        calc_ut_uphi.calculate_ut_uphi_from_r(
            circ_orbit_ID.ut_stable_expr, circ_orbit_ID.uphi_stable_expr
        )
        set_ic.set_initial_conditions_massive()
        gen_disk_ic.generate_disk_initial_conditions()
        gen_spiral_ic.generate_spiral_galaxy_initial_conditions()
        gen_barred_ic.generate_barred_flocculent_spiral_ic()
        int_single.integrate_single_particle()
        int_debug.integrate_single_particle_DEBUG()
        run_prod.run_mass_integrator_production()
        main_c_function.main_c_function()
        print("    ... All C functions registered.")

        # Step 1.H: Assemble the final C project on disk.
        print(" -> Assembling final C project directory...")
        register_custom_structures_and_params()
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
            addl_CFLAGS=["-Wall -Wextra -g $(shell gsl-config --cflags) -fopenmp"],
            addl_libraries=["$(shell gsl-config --libs) -fopenmp"],
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

            # 3. Run the Python post-processing step, ONLY if the C code succeeded.
            if c_code_succeeded:
                print(f"\n[{project_name}] Starting Python post-processing...")
                snapshot_input_dir = os.path.join(project_dir, "output")
                kdtree_utils.run_kdtree_preprocessor_parallel(snapshot_input_dir, snapshot_output_dir)
                print(f"[{project_name}] Post-processing complete.")

        else:
            # If --run is not specified, print manual instructions.
            print(f"To build, navigate to '{project_dir}/' and type 'make'.")
            print(f"To run, type './{exec_name}'.")
