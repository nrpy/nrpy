"""
Construct a complete C project for integrating photon geodesics in curved spacetime.

Project: NRPy+ Standalone Geodesic Integrator
Description:
    This script acts as the primary driver to generate a standalone C application
    that evolves the trajectory of a massless photon test particle. It coordinates the
    generation of spacetime-specific physics kernels (Metric, Christoffel Symbols, ODE RHS)
    and links them with the GNU Scientific Library (GSL) for high-order time integration.

Physics Context:
    The simulation solves the geodesic equation for a photon:
        d(p^mu)/d(lambda) = -Gamma^mu_{alpha beta} p^alpha p^beta
    subject to the normalization constraint p^mu p_mu = 0.

    Numerical fidelity is rigorously validated by monitoring constants of motion
    associated with the spacetime's symmetries (Killing vectors and tensors):
    1. Energy (E): Associated with time-translation invariance (Killing vector dt).
    2. Axial Angular Momentum (Lz): Associated with rotational invariance (Killing vector dphi).
    3. Carter Constant (Q): Associated with the hidden symmetry of the Kerr metric.

Author: Dalton J. Moone
"""

import os
import shutil
import subprocess
import sys

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.params as par
from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
    Analytic_Spacetimes,
)
from nrpy.equations.general_relativity.geodesics.geodesics import Geodesic_Equations
from nrpy.infrastructures.BHaH import cmdline_input_and_parfiles
from nrpy.infrastructures.BHaH.general_relativity.geodesics.connections import (
    connections,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.conserved_quantities import (
    conserved_quantities,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.g4DD_metric import (
    g4DD_metric,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.normalization_constraint import (
    normalization_constraint,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.calculate_ode_rhs import (
    calculate_ode_rhs,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.ode_gsl_wrapper import (
    ode_gsl_wrapper,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.p0_reverse import (
    p0_reverse,
)

# Set codegen and compile-time parameters for the BHaH infrastructure
par.set_parval_from_str("Infrastructure", "BHaH")

# Identifier for the generated C project directory and executable
project_name = "photon_geodesic_integrator"
project_dir = os.path.join("project", project_name)

# Clean the project directory if it already exists to ensure a fresh build
if os.path.exists(project_dir):
    shutil.rmtree(project_dir)

# Spacetime and particle configuration
SPACETIME = "KerrSchild_Cartesian"
PARTICLE = "photon"
GEO_KEY = f"{SPACETIME}_{PARTICLE}"

#########################################################
# Register Physics C Functions
# This generates the computational kernels.

print("Acquiring symbolic data...")
# Extract symbolic expressions for the chosen spacetime and particle type
metric_data = Analytic_Spacetimes[SPACETIME]
geodesic_data = Geodesic_Equations[GEO_KEY]

print("Registering C functions...")
# 1. Metric: Registers the spacetime metric tensor evaluation
g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)

# 2. Connections: Registers Christoffel Symbol evaluation
connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)

# 3. ODE Right-Hand Side: Registers the massless geodesic equation evaluation
calculate_ode_rhs(geodesic_data.geodesic_rhs, metric_data.xx)

# 4. Hamiltonian Constraint Solver: Resolves initial p^0 based on spatial momentum
if geodesic_data.p0_photon is None:
    raise ValueError(f"p0_photon is None for {GEO_KEY}")
p0_reverse(geodesic_data.p0_photon)

# 5. Conserved Quantities: Registers functions to track constants of motion
conserved_quantities(SPACETIME, PARTICLE)

# 6. Normalization Constraint: Tracks deviation from p^mu p_mu = 0
normalization_constraint(geodesic_data.norm_constraint_expr, PARTICLE)

# 7. GSL Wrapper: Provides the C-interface needed by the GNU Scientific Library
ode_gsl_wrapper(SPACETIME)


#########################################################
# Declare the main C function
# This drives the integration logic and file I/O.


def main_c() -> None:
    """Generate the main() function for the photon geodesic integrator."""
    # 1. Define C-Function metadata
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "gsl/gsl_errno.h",
        "gsl/gsl_odeiv2.h",
        "gsl/gsl_matrix.h",
        "gsl/gsl_math.h",
        "string.h",
    ]

    desc = """@brief Main driver function for the photon geodesic integrator.
    Detailed algorithm: Initializes the BHaH infrastructure and establishes initial conditions 
    for a massless test particle. Solves for the initial time-component of the 4-momentum 
    via the Hamiltonian constraint. Iteratively evolves the trajectory using an adaptive 
    Runge-Kutta-Fehlberg (RKF45) integrator from the GNU Scientific Library (GSL), 
    logging position, momentum, path length, and constants of motion to validate numerical fidelity."""

    cfunc_type = "int"
    name = "main"
    params = "int argc, char *argv[]"

    # 2. Build the C body with internal descriptive comments
    body = f"""
    // --- Step 1: Setup Infrastructure ---
    // commondata: Struct containing parameters common to all grids and physics configurations.
    commondata_struct commondata;
    commondata_struct_set_to_default(&commondata);

    // Hardcode parameter overrides for this specific Kerr-Schild test case.
    commondata.M_scale = 1.0;
    commondata.a_spin = 0.9;

    printf("Starting Photon Geodesic Integrator...\\n");
    printf("Spacetime {SPACETIME}, M=%.2f, a=%.2f\\n", commondata.M_scale, commondata.a_spin);

    // --- Step 2: Initial Conditions ---
    // y: State vector of length 9. Components map to:
    // [0]=tau, [1]=x, [2]=y, [3]=z, [4]=p^t, [5]=p^x, [6]=p^y, [7]=p^z, [8]=L (Path Length).
    double y[9];
    y[0] = 0.0;   // Coordinate time t
    y[1] = 10.0;  // Spatial coordinate x
    y[2] = 1.0;   // Spatial coordinate y
    y[3] = 1.0;   // Spatial coordinate z

    // Set initial spatial momentum (p^i).
    y[5] = -0.1; // p^x
    y[6] = 0.33; // p^y
    y[7] = 0.0;  // p^z
    y[8] = 0.0;  // L (initial path length)

    // g4DD_local: Flat array allocated to hold the 10 independent components of the symmetric 4D metric tensor.
    double g4DD_local[10];

    // Calculate metric at initial position y (batch size 1, batch ID 0).
    g4DD_metric_{SPACETIME}(&commondata, y, g4DD_local, 1, 0);

    // p0_val: Variable to store the computed time-component of the 4-momentum (p^t).
    double p0_val = 0.0;
    
    // Solve for p^0 using the pre-calculated metric to enforce the normalization constraint.
    // Signature: (metric, state_vector, num_rays, batch_size, photon_idx, batch_id, p0_out)
    p0_reverse(g4DD_local, y, 1, 1, 0, 0, &p0_val);
    y[4] = p0_val;
    // ---------------------------------------------------------

    printf("Initial State\\n");
    printf("  Pos (%.4f, %.4f, %.4f)\\n", y[1], y[2], y[3]);
    printf("  Mom (%.4f, %.4f, %.4f, %.4f)\\n", y[4], y[5], y[6], y[7]);
    printf("  Length (%.4f)\\n", y[8]);

    // --- Step 3: Pre-Integration Diagnostics ---
    // Variables to hold the initial values of constants of motion for baseline comparison.
    double E_init, Lx_init, Ly_init, Lz_init, Q_init;
    conserved_quantities_{SPACETIME}_{PARTICLE}(&commondata, y, 1, 0,
                                                &E_init, &Lx_init, &Ly_init, &Lz_init, &Q_init);

    printf("Initial Conserved Quantities\\n");
    printf("  E = %.8f, Lz = %.8f, Q = %.8f\\n", E_init, Lz_init, Q_init);

    // --- Step 4: GSL Integrator Setup ---
    // T: GSL stepper type, utilizing the Runge-Kutta-Fehlberg (RKF45) method.
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
    // s: GSL stepper object dynamically allocated for a 9-dimensional state vector.
    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc(T, 9);
    // c: GSL control object to maintain local truncation error limits (absolute error 1e-9).
    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new(1e-9, 0.0);
    // e: GSL evolution object tracking current integration state and dynamically updating step sizes.
    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc(9);

    // sys: Struct binding our generated ODE RHS function to the GSL framework.
    gsl_odeiv2_system sys = {{ode_gsl_wrapper_{SPACETIME}, NULL, 9, &commondata}};

    // lambda: Tracks the affine parameter accumulated by the photon.
    double lambda = 0.0;
    // lambda_max: The predefined boundary limit for affine parameter integration.
    double lambda_max = 20000.0;
    // h: Initial step size guess provided to the adaptive GSL routines.
    double h = 1e-3;

    // --- Step 5: File Output Setup ---
    // fp: File pointer directed to write trajectory data into a delimited text format.
    FILE *fp = fopen("trajectory.txt", "w");
    if (fp == NULL) {{
        fprintf(stderr, "Error opening trajectory.txt\\n");
        return 1;
    }}
    fprintf(fp, "# lambda t x y z p^t p^x p^y p^z L\\n");

    // --- Step 6: Integration Loop ---
    // steps: Counter incremented per successful step to prevent runaway loops.
    int steps = 0;
    // max_steps: Absolute hard ceiling on integration iterations.
    int max_steps = 2000000;

    while (lambda < lambda_max && steps < max_steps) {{
        fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\\n",
                lambda, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]);

        // status: GSL exit code reporting the success or failure of the internal state advancement.
        int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &lambda, lambda_max, &h, y);

        if (status != GSL_SUCCESS) {{
            printf("GSL Error %d\\n", status);
            break;
        }}

        // r: Radial Cartesian distance from the origin to dynamically assess horizon crossing.
        double r = sqrt(y[1]*y[1] + y[2]*y[2] + y[3]*y[3]);
        if (r < 2.0 * commondata.M_scale) {{
            printf("Termination Particle reached r = %.4f < 2M at lambda = %.4f\\n", r, lambda);
            break;
        }}
        steps++;
    }}

    fclose(fp);
    printf("Integration finished after %d steps. Final lambda = %.4f\\n", steps, lambda);

    // --- Step 7: Post-Integration Diagnostics ---
    // Variables capturing the final states of conserved quantities for global error analysis.
    double E_final, Lx_final, Ly_final, Lz_final, Q_final;
    // norm_final: Metric normalization state post-integration (p^mu p_mu = 0).
    double norm_final;
    
    conserved_quantities_{SPACETIME}_{PARTICLE}(&commondata, y, 1, 0,
                                                &E_final, &Lx_final, &Ly_final, &Lz_final, &Q_final);   

    g4DD_metric_{SPACETIME}(&commondata, y, g4DD_local, 1, 0);
    
    // Evaluate Normalization constraint using flat array and SoA parameters.
    normalization_constraint_photon(g4DD_local, y, 1, 1, 0, 0, &norm_final);

    printf("Final norm \\n");
    printf("  norm = %.4e\\n", norm_final);

    printf("Final Conserved Quantities\\n");
    printf("  E = %.8f, Lz = %.8f, Q = %.8f\\n", E_final, Lz_final, Q_final);

    // Absolute errors computed against initial conditions to verify integration stability.
    double E_err = fabs(E_final - E_init);
    double Lz_err = fabs(Lz_final - Lz_init);
    double Q_err = fabs(Q_final - Q_init);

    printf("Conservation Check (Absolute Error)\\n");
    printf("  Delta E  = %.4e\\n", E_err);
    printf("  Delta Lz = %.4e\\n", Lz_err);
    printf("  Delta Q  = %.4e\\n", Q_err);

    // Clean up dynamically allocated GSL memory constructs.
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);

    return 0;
    """

    # 3. Register the C function
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


# Register the main function
main_c()

#########################################################
# Generate Header Files and Makefile
# This creates the actual build files.

print("Generating header files and Makefile...")

# A. CodeParameters Headers
CPs.write_CodeParameters_h_files(set_commondata_only=True, project_dir=project_dir)
CPs.register_CFunctions_params_commondata_struct_set_to_default()

# B. Parameter File Defaults (required by infrastructure, even if unused)
cmdline_input_and_parfiles.generate_default_parfile(
    project_dir=project_dir, project_name=project_name
)
cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
    project_name=project_name
)

# C. BHaH Defines (Includes GSL headers)
additional_includes = [
    "gsl/gsl_vector.h",
    "gsl/gsl_matrix.h",
    "gsl/gsl_odeiv2.h",
    "gsl/gsl_errno.h",
    "gsl/gsl_math.h",
]

# Ensure hardware-agnostic array indexing for standalone GSL testing
macro_defs = """
#ifndef IDX_LOCAL
#define IDX_LOCAL(component, batch_id, batch_size) ((component) * (batch_size) + (batch_id))
#endif

#ifndef IDX_GLOBAL
#define IDX_GLOBAL(component, ray_id, num_rays) ((component) * (num_rays) + (ray_id))
#endif
"""
Bdefines_h.register_BHaH_defines("gpu_batch_macros", macro_defs)

Bdefines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    additional_includes=additional_includes,
    enable_rfm_precompute=False,
)

# D. Generate the Makefile FIRST
# Flags required to link the external GSL library during compilation
addl_cflags = ["$(shell gsl-config --cflags)"]
addl_libs = ["$(shell gsl-config --libs)"]

print(" -> Generating Makefile...")
Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=project_name,
    addl_CFLAGS=addl_cflags,
    addl_libraries=addl_libs,
)

# E. Patch the Makefile
print(" -> Patching Makefile for Windows compatibility...")
local_tmp_path = "tmp"
os.makedirs(os.path.join(project_dir, local_tmp_path), exist_ok=True)

makefile_path = os.path.join(project_dir, "Makefile")

# Read the generated Makefile
with open(makefile_path, "r", encoding="utf-8") as f:
    content = f.read()

# Overwrite it, prepending the temporary directory exports to bypass Windows pathing issues
with open(makefile_path, "w", encoding="utf-8") as f:
    f.write(f"export TMPDIR = $(CURDIR)/{local_tmp_path}\n")
    f.write(f"export TMP = $(CURDIR)/{local_tmp_path}\n")
    f.write(f"export TEMP = $(CURDIR)/{local_tmp_path}\n")
    f.write("\n")
    f.write(content)

print("-" * 50)
print(f"Project generated successfully in {project_dir}")

##########################################################################
# PART 2: PIPELINE EXECUTION (COMPILE, RUN, VISUALIZE)
##########################################################################

if __name__ == "__main__":
    import logging

    import matplotlib.pyplot as plt
    import numpy as np

    # Mute matplotlib's verbose font debugging output
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)

    print("\n" + "=" * 50)
    print("PIPELINE EXECUTION: COMPILE, RUN, VISUALIZE")
    print("=" * 50)

    print("\n--- PHASE 1: Compiling C Code ---")

    # Use absolute path to prevent "Makefile not found" errors in subprocess
    abs_project_dir = os.path.abspath(project_dir)
    try:
        subprocess.run(["make", "-j"], cwd=abs_project_dir, check=True)
        print("Compilation successful.")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Compilation failed: {e}")
        print(
            "Tip: Ensure 'make' (MinGW/MSYS2) is in your PATH and the Makefile was generated."
        )
        sys.exit(1)

    print("\n--- PHASE 2: Running Ray-Tracer ---")

    # Conditionally set the executable extension for Windows environments
    exe_extension = (
        ".exe" if os.name == "nt" or sys.platform in ["win32", "cygwin", "msys"] else ""
    )
    # Full path to the compiled integrator executable
    exec_path = os.path.join(abs_project_dir, f"{project_name}{exe_extension}")

    try:
        subprocess.run([exec_path], cwd=abs_project_dir, check=True)
        print("Ray-tracing complete. Trajectory file generated.")
    except subprocess.CalledProcessError:
        print("C executable failed. Exiting pipeline.")
        sys.exit(1)

    print("\n--- PHASE 3: Visualizing Trajectory ---")

    # Path to the text file containing the output ray-tracing state data
    traj_file = os.path.join(abs_project_dir, "trajectory.txt")
    if not os.path.exists(traj_file):
        print(f"Error: {traj_file} not found.")
        sys.exit(1)

    try:
        # Parse the trajectory metrics into a 2D NumPy array
        data = np.loadtxt(traj_file, comments="#")
        x_pts = data[:, 2]
        y_pts = data[:, 3]
        z_pts = data[:, 4]

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")

        # Plot the primary geodesic path
        ax.plot(
            x_pts, y_pts, z_pts, label="Photon Trajectory", color="blue", linewidth=1.5
        )

        # Drop markers for the simulation start and endpoints
        ax.scatter(
            x_pts[0], y_pts[0], z_pts[0], color="green", marker="o", s=50, label="Start"
        )
        ax.scatter(
            x_pts[-1],
            y_pts[-1],
            z_pts[-1],
            color="red",
            marker="x",
            s=50,
            label="End (r < 2M)",
        )

        # Black Hole Horizon setup
        M_scale = 1.0
        r_horizon = 2.0 * M_scale

        # Parameterize angles for generating the spherical horizon surface
        u_val = np.linspace(0, 2 * np.pi, 20)
        v_val = np.linspace(0, np.pi, 10)
        # Create 2D meshes mapping the spherical coordinates
        u, v = np.meshgrid(u_val, v_val, indexing="ij")

        # Convert to Cartesian coordinates to plot the event horizon
        xh = r_horizon * np.cos(u) * np.sin(v)
        yh = r_horizon * np.sin(u) * np.sin(v)
        zh = r_horizon * np.cos(v)
        ax.plot_surface(xh, yh, zh, color="black", alpha=0.3, label="Horizon")

        ax.set_xlabel("x (M)")
        ax.set_ylabel("y (M)")
        ax.set_zlabel("z (M)")
        ax.set_title("Photon Geodesic in Kerr-Schild Cartesian Spacetime")
        ax.legend()

        # Location to save the rendered matplotlib figure
        plot_path = os.path.join(abs_project_dir, "photon_trajectory.png")
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        print(f"Visualization successfully saved to: {plot_path}")

        # Open the interactive UI
        print("Opening plot window...")
        plt.show()

    except (RuntimeError, ValueError, OSError) as e:
        print(f"Plotting failed: {e}")
        sys.exit(1)
