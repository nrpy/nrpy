"""
Construct a complete C project for integrating photon geodesics in curved spacetime.

This script acts as the primary driver to generate a standalone C application that
evolves the trajectory of a photon test particle. It coordinates the generation of
spacetime-specific physics kernels (Metric, Christoffel Symbols, ODE RHS) and
links them with the GNU Scientific Library (GSL) for high-order time integration.

Physics Context
    The simulation solves the geodesic equation for a photon
        d(p^mu)/d(lambda) = -Gamma^mu_{alpha beta} p^alpha p^beta
    subject to the normalization constraint p^mu p_mu = 0.

    Numerical fidelity is validated by monitoring constants of motion associated
    with the spacetime's symmetries (Killing vectors and tensors)
    1. Energy (E) Associated with time-translation invariance (Killing vector dt).
    2. Axial Angular Momentum (Lz) Associated with rotational invariance (Killing vector dphi).
    3. Carter Constant (Q) Associated with the hidden symmetry of the Kerr metric.

Author Dalton J. Moone
"""

import os
import shutil
import subprocess
import sys

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile

# Import needed Python modules
import nrpy.params as par

# Import physics generation modules
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

# Set codegen and compile-time parameters
par.set_parval_from_str("Infrastructure", "BHaH")

project_name = "photon_geodesic_integrator"
project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
if os.path.exists(project_dir):
    shutil.rmtree(project_dir)

# Spacetime configuration
SPACETIME = "KerrSchild_Cartesian"
PARTICLE = "photon"
GEO_KEY = f"{SPACETIME}_{PARTICLE}"

#########################################################
# Register Physics C Functions
# This generates the computational kernels.

print("Acquiring symbolic data...")
metric_data = Analytic_Spacetimes[SPACETIME]
geodesic_data = Geodesic_Equations[GEO_KEY]

print("Registering C functions...")

# 1. Metric
g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)

# 2. Connections (Christoffel Symbols)
connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)

# 3. ODE Right-Hand Side
calculate_ode_rhs(geodesic_data.geodesic_rhs, metric_data.xx)

# 4. Hamiltonian Constraint Solver (for initial p^0)
if geodesic_data.p0_photon is None:
    raise ValueError(f"p0_photon is None for {GEO_KEY}")
p0_reverse(geodesic_data.p0_photon)

# 5. Conserved Quantities (Diagnostics)
conserved_quantities(SPACETIME, PARTICLE)

# 6. Normalization Constraint
normalization_constraint(geodesic_data.norm_constraint_expr, PARTICLE)

# 7. GSL Wrapper
ode_gsl_wrapper(SPACETIME)


#########################################################
# Declare the main C function
# This drives the integration logic.


def main_c() -> None:
    """Generate the main() function for the geodesic integrator."""
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

        Initializes the BHaH infrastructure, sets initial particle conditions, 
        performs time integration using GSL (RKF45), and outputs trajectory data 
        and conservation checks to stdout/files."""
    cfunc_type = "int"
    name = "main"
    params = "int argc, char *argv[]"

    body = f"""
    // 1. Setup Infrastructure
    commondata_struct commondata;
    commondata_struct_set_to_default(&commondata);

    // Hardcode parameters overrides for this specific test case
    commondata.M_scale = 1.0;
    commondata.a_spin = 0.9;

    printf("Starting Photon Geodesic Integrator...\\n");
    printf("Spacetime {SPACETIME}, M=%.2f, a=%.2f\\n", commondata.M_scale, commondata.a_spin);

    // 2. Initial Conditions
    double y[9];
    y[0] = 0.0;
    y[1] = 10.0;
    y[2] = 1.0;
    y[3] = 1.0;

    // Initial Spatial Momentum
    y[5] = -0.1; // p^x
    y[6] = 0.33;  // p^y
    y[7] = 0.0;   // p^z
    y[8] = 0.0;  // L 

    // A. Declare flat array to hold metric components
    double g4DD_local[10];

    // B. Calculate metric at initial position y (batch size 1, batch ID 0)
    g4DD_metric_{SPACETIME}(&commondata, y, g4DD_local, 1, 0);

    // C. Solve for p^0 using the pre-calculated metric 
    double p0_val = 0.0;
    // Signature (metric, state_vector, num_rays, batch_size, photon_idx, batch_id, p0_out)
    p0_reverse(g4DD_local, y, 1, 1, 0, 0, &p0_val);
    y[4] = p0_val;
    // ---------------------------------------------------------

    printf("Initial State\\n");
    printf("  Pos (%.4f, %.4f, %.4f)\\n", y[1], y[2], y[3]);
    printf("  Mom (%.4f, %.4f, %.4f, %.4f)\\n", y[4], y[5], y[6], y[7]);
    printf("  Length (%.4f)\\n", y[8]);

    // 3. Pre-Integration Diagnostics (batch size 1, photon ID 0)
    double E_init, Lx_init, Ly_init, Lz_init, Q_init;
    conserved_quantities_{SPACETIME}_{PARTICLE}(&commondata, y, 1, 0,
                                                &E_init, &Lx_init, &Ly_init, &Lz_init, &Q_init);

    printf("Initial Conserved Quantities\\n");
    printf("  E = %.8f, Lz = %.8f, Q = %.8f\\n", E_init, Lz_init, Q_init);

    // 4. GSL Setup
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc(T, 9);
    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new(1e-9, 0.0);
    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc(9);

    gsl_odeiv2_system sys = {{ode_gsl_wrapper_{SPACETIME}, NULL, 9, &commondata}};

    // Integration variables
    double lambda = 0.0;
    double lambda_max = 20000.0;
    double h = 1e-3;

    // 5. File Output Setup
    FILE *fp = fopen("trajectory.txt", "w");
    if (fp == NULL) {{
        fprintf(stderr, "Error opening trajectory.txt\\n");
        return 1;
    }}
    fprintf(fp, "# lambda t x y z p^t p^x p^y p^z L\\n");

    // 6. Integration Loop
    int steps = 0;
    int max_steps = 2000000;

    while (lambda < lambda_max && steps < max_steps) {{
        fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\\n",
                lambda, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]);

        int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &lambda, lambda_max, &h, y);

        if (status != GSL_SUCCESS) {{
            printf("GSL Error %d\\n", status);
            break;
        }}

        double r = sqrt(y[1]*y[1] + y[2]*y[2] + y[3]*y[3]);
        if (r < 2.0 * commondata.M_scale) {{
            printf("Termination Particle reached r = %.4f < 2M at lambda = %.4f\\n", r, lambda);
            break;
        }}
        steps++;
    }}

    fclose(fp);
    printf("Integration finished after %d steps. Final lambda = %.4f\\n", steps, lambda);

    // 7. Post-Integration Diagnostics
    double E_final, Lx_final, Ly_final, Lz_final, Q_final, norm_final;
    conserved_quantities_{SPACETIME}_{PARTICLE}(&commondata, y, 1, 0,
                                                &E_final, &Lx_final, &Ly_final, &Lz_final, &Q_final);   

    g4DD_metric_{SPACETIME}(&commondata, y, g4DD_local, 1, 0);
    
    // Evaluate Normalization constraint using flat array and SoA parameters
    normalization_constraint_photon(g4DD_local, y, 1, 1, 0, 0, &norm_final);

    printf("Final norm \\n");
    printf("  norm = %.4e\\n", norm_final);


    printf("Final Conserved Quantities\\n");
    printf("  E = %.8f, Lz = %.8f, Q = %.8f\\n", E_final, Lz_final, Q_final);

    double E_err = fabs(E_final - E_init);
    double Lz_err = fabs(Lz_final - Lz_init);
    double Q_err = fabs(Q_final - Q_init);

    printf("Conservation Check (Absolute Error)\\n");
    printf("  Delta E  = %.4e\\n", E_err);
    printf("  Delta Lz = %.4e\\n", Lz_err);
    printf("  Delta Q  = %.4e\\n", Q_err);

    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);

    return 0;
    """

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
# We use standard strings, not Path objects, to ensure valid C syntax on all OSs
additional_includes = [
    "gsl/gsl_vector.h",
    "gsl/gsl_matrix.h",
    "gsl/gsl_odeiv2.h",
    "gsl/gsl_errno.h",
    "gsl/gsl_math.h",
]

macro_defs = """
// Ensure hardware-agnostic array indexing for standalone GSL test
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

# Read the Makefile that was just generated
with open(makefile_path, "r", encoding="utf-8") as f:
    content = f.read()

# Overwrite it, prepending the temporary directory exports
with open(makefile_path, "w", encoding="utf-8") as f:
    # Using $(CURDIR) ensures these paths are absolute and valid during the build
    f.write(f"export TMPDIR = $(CURDIR)/{local_tmp_path}\n")
    f.write(f"export TMP = $(CURDIR)/{local_tmp_path}\n")
    f.write(f"export TEMP = $(CURDIR)/{local_tmp_path}\n")
    f.write("\n")  # Add a newline for safety
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
    try:
        subprocess.run(["make", "-j"], cwd=project_dir, check=True)
        print("Compilation successful.")
    except subprocess.CalledProcessError:
        print("Compilation failed. Exiting pipeline.")
        sys.exit(1)

    print("\n--- PHASE 2: Running Ray-Tracer ---")

    # Conditionally add .exe for Windows/Git Bash compatibility
    exe_extension = (
        ".exe" if os.name == "nt" or sys.platform in ["win32", "cygwin", "msys"] else ""
    )
    # Use the absolute path to ensure Windows can find the executable
    exec_path = os.path.join(project_dir, f"{project_name}{exe_extension}")

    try:
        subprocess.run([exec_path], cwd=project_dir, check=True)
        print("Ray-tracing complete. Trajectory file generated.")
    except subprocess.CalledProcessError:
        print("C executable failed. Exiting pipeline.")
        sys.exit(1)

    print("\n--- PHASE 3: Visualizing Trajectory ---")
    traj_file = os.path.join(project_dir, "trajectory.txt")
    if not os.path.exists(traj_file):
        print(f"Error: {traj_file} not found.")
        sys.exit(1)

    try:
        # Load trajectory data: lambda, t, x, y, z, p^t, p^x, p^y, p^z, L
        data = np.loadtxt(traj_file, comments="#")
        x_pts = data[:, 2]
        y_pts = data[:, 3]
        z_pts = data[:, 4]

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")

        # Plot the geodesic
        ax.plot(
            x_pts, y_pts, z_pts, label="Photon Trajectory", color="blue", linewidth=1.5
        )

        # Mark Start and End points
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

        # Add a visual sphere for the event horizon (r = 2M, where M = 1.0)
        M_scale = 1.0
        r_horizon = 2.0 * M_scale
        # Use linspace + meshgrid to avoid mypy issues with complex slicing in mgrid
        u_val = np.linspace(0, 2 * np.pi, 20)
        v_val = np.linspace(0, np.pi, 10)
        u, v = np.meshgrid(u_val, v_val, indexing="ij")
        xh = r_horizon * np.cos(u) * np.sin(v)
        yh = r_horizon * np.sin(u) * np.sin(v)
        zh = r_horizon * np.cos(v)
        ax.plot_surface(xh, yh, zh, color="black", alpha=0.3, label="Horizon")

        # Formatting
        ax.set_xlabel("x (M)")
        ax.set_ylabel("y (M)")
        ax.set_zlabel("z (M)")
        ax.set_title("Photon Geodesic in Kerr-Schild Cartesian Spacetime")
        ax.legend()

        # Save Plot
        plot_path = os.path.join(project_dir, "photon_trajectory.png")
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        print(f"Visualization successfully saved to: {plot_path}")

        # Pop up the interactive window
        plt.show()

    except (RuntimeError, ValueError, OSError) as e:
        print(f"Plotting failed: {e}")
        sys.exit(1)
