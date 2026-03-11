"""
Construct a complete C project for integrating massive geodesics in curved spacetime.

Project: NRPy+ Standalone Geodesic Integrator
Description:
    This script generates a standalone C application that evolves the trajectory 
    of a massive test particle. It coordinates the generation of spacetime-specific 
    physics kernels and links them with the GNU Scientific Library (GSL) for 
    high-order time integration.

Physics Context:
    The simulation solves the geodesic equation for a particle with non-zero mass:
        $d(u^\\mu)/d(\\tau) = -\\Gamma^\\mu_{\\alpha \\beta} u^\\alpha u^\\beta$
    subject to the normalization constraint $u^\\mu u_\\mu = -1$.

    Numerical fidelity is rigorously validated by monitoring constants of motion
    associated with the spacetime's symmetries (Killing vectors and tensors).

Author: Dalton J. Moone.
"""

import os
import shutil
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
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive.calculate_ode_rhs_massive import (
    calculate_ode_rhs_massive,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive.ode_gsl_wrapper_massive import (
    ode_gsl_wrapper_massive,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive.u0_massive import (
    u0_massive,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.normalization_constraint import (
    normalization_constraint,
)


def register_struct_definitions() -> None:
    """
    Register the Structure of Arrays (SoA) and diagnostic structures into the global C headers.
    """
    termination_enum_def = """
    // Defines the specific exit condition for a particle's integration loop.
    typedef enum {
        ACTIVE = 0,                        // 0: Particle is currently undergoing integration.
        REJECTED,                          // 1: Particle RKF45 step was rejected.
        FAILURE_RKF45_REJECTION_LIMIT      // 2: Adaptive step-size rejected too many consecutive times.
    } termination_type_t;
    """
    Bdefines_h.register_BHaH_defines("termination_type_t", termination_enum_def)

    photon_soa_def = """
    // ==========================================
    // Flattened SoA Struct (Master Storage)
    // ==========================================
    // Hardware Justification: This Structure of Arrays (SoA) provides universal memory compatibility for shared physics kernels.
    typedef struct {
        double *f; // Flattened state vector mapping components $t, x, y, z, u_t, u_x, u_y, u_z$.
        double *affine_param; // Current proper time $\\tau$ for the trajectory.
    } PhotonStateSoA;
    """
    Bdefines_h.register_BHaH_defines("PhotonStateSoA", photon_soa_def)


def main_c(SPACETIME: str, PARTICLE: str) -> None:
    """
    Generate the main() function for the massive geodesic integrator using GSL.

    :param SPACETIME: The specific background spacetime descriptor.
    :param PARTICLE: The type of test particle being integrated.
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "gsl/gsl_errno.h",
        "gsl/gsl_odeiv2.h",
        "gsl/gsl_matrix.h",
        "gsl/gsl_math.h",
        "string.h",
        "stdlib.h"
    ]

    desc = """@brief Main driver function for the massive geodesic integrator.
    Detailed algorithm: Initializes memory for a single ray and executes continuous 
    integration via the GNU Scientific Library (GSL) adaptive RKF45 stepper. Evaluates 
    boundary constraints and verifies the numerical fidelity using Hamiltonian constraints 
    and conserved quantities."""

    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"

    body = f"""
    // --- Step 1: Structural Setup & Parameters ---
    commondata_struct commondata; // Global parameters struct governing spacetime and numerics.
    commondata_struct_set_to_default(&commondata);
    
    commondata.M_scale = 1.0; // The mass $M$ of the central black hole.
    commondata.a_spin = 0.9; // The dimensionless spin $a$ of the central black hole.

    printf("Starting Mass Geodesic Integrator (CPU)...\\n");
    printf("Spacetime: {SPACETIME}, M=%.2f, a=%.2f\\n", commondata.M_scale, commondata.a_spin);

    // --- Step 2: Global Memory & Initial Conditions ---
    const double p_t_max = 1000;
    const double r_squared_max = 1e4;
    long int num_rays = 1; // Total number of global particle trajectories.

    double y[8]; // State vector of length 8 mapping components $t, x, y, z, u^t, u^x, u^y, u^z$.
    y[0] = 0.0;  // The initial coordinate time $t$.
    y[1] = 10.0; // The initial spatial coordinate $x$.
    y[2] = 1.0;  // The initial spatial coordinate $y$.
    y[3] = 1.0;  // The initial spatial coordinate $z$.

    y[5] = -0.1; // The initial spatial velocity $u^x$.
    y[6] = 0.33; // The initial spatial velocity $u^y$.
    y[7] = 0.0;  // The initial spatial velocity $u^z$.

    double tau = 0.0; // Tracks the proper time $\\tau$ accumulated by the massive particle.

    // Hardware Justification: Struct mapping establishes the SoA layout required by the downstream conserved quantities kernel.
    PhotonStateSoA all_particles; // Master struct holding array pointers for the global dataset.
    all_particles.f = y; // Pointer to the flattened state vector $y$.
    all_particles.affine_param = &tau; // Pointer to the proper time $\\tau$.

    double g4DD_local[10]; // Flat array holding the 10 independent components of the symmetric metric $g_{{\\mu\\nu}}$.

    // Hardware Justification: Calculate metric at initial position $y$ to solve the Hamiltonian constraint.
    g4DD_metric_{SPACETIME}(&commondata, y, g4DD_local);

    double u0_val = 0.0; // Variable to store the computed time-component of the 4-velocity $u^t$.
    u0_massive(g4DD_local, y, 1, 1, 0, 0, &u0_val);
    y[4] = u0_val; // Assign the computed temporal velocity $u^t$.

    printf("Initial State:\\n");
    printf("  Pos: (%.4f, %.4f, %.4f)\\n", y[1], y[2], y[3]);
    printf("  Vel: (%.4f, %.4f, %.4f, %.4f)\\n", y[4], y[5], y[6], y[7]);

    // --- Step 3: Pre-Integration Diagnostics ---
    conserved_quantities_t cq_init; // Struct instance to store the baseline constants of motion.
    calculate_conserved_quantities_universal_{SPACETIME}_{PARTICLE}(&commondata, &all_particles, num_rays, &cq_init);

    printf("Initial Conserved Quantities:\\n");
    printf("  E = %.8f, Lz = %.8f, Q = %.8f\\n", cq_init.E, cq_init.Lz, cq_init.Q);

    // --- Step 4: GSL Integrator Setup ---
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45; // GSL stepper type utilizing the RKF45 method.
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(T, 8); // GSL stepper object dynamically allocated for an 8-dimensional state vector.
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-9, 0.0); // GSL control object to maintain local truncation error limits.
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(8); // GSL evolution object tracking current integration state.

    gsl_odeiv2_system sys = {{ode_gsl_wrapper_massive_{SPACETIME}, NULL, 8, &commondata}}; // Struct binding our ODE RHS function to the GSL framework.

    double tau_max = 20000.0; // The predefined boundary limit for proper time $\\tau$ integration.
    double h = 1e-3; // Initial step size guess $h$ provided to the adaptive GSL routines.

    // --- Step 5: File Output Setup ---
    FILE *fp = fopen("trajectory.txt", "w"); // File pointer directed to write trajectory data.
    if (fp == NULL) {{
        fprintf(stderr, "Error opening trajectory.txt\\n");
        return 1;
    }}
    fprintf(fp, "# tau t x y z u^t u^x u^y u^z\\n");

    // --- Step 6: Integration Loop ---
    int steps = 0; // Counter incremented per successful step to prevent runaway loops.
    int max_steps = 2000000; // Absolute hard ceiling on integration iterations.

    while (tau < tau_max && steps < max_steps) {{
        fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\\n",
                tau, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);

        int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &tau, tau_max, &h, y); // GSL exit code reporting the success or failure of state advancement.

        if (status != GSL_SUCCESS) {{
            printf("GSL Error: %d\\n", status);
            break;
        }}

        if (fabs(y[4]) > p_t_max) {{
            printf("Temporal momentum $|p_t|$ exceeded numerical limit.\\n");
            break;
        }}

        double r_squared = y[1]*y[1] + y[2]*y[2] + y[3]*y[3]; // Radial Cartesian distance $r_squared$ from the origin.
        if (r_squared > r_squared_max ) {{
            printf("Termination: Particle reached r_squared = %.4f at tau = %.4f\\n", r_squared, tau);
            break;
        }}
        steps++;
    }}

    fclose(fp);
    printf("Integration finished after %d steps. Final tau = %.4f\\n", steps, tau);

    // --- Step 7: Post-Integration Diagnostics ---
    conserved_quantities_t cq_final; // Struct holding constants of motion evaluated at the trajectory boundary.
    calculate_conserved_quantities_universal_{SPACETIME}_{PARTICLE}(&commondata, &all_particles, num_rays, &cq_final);
    
    g4DD_metric_{SPACETIME}(&commondata, y, g4DD_local);
    
    normalization_constraint_t norm_final; // Struct tracking the residual of the Hamiltonian constraint $u^\\mu u_\\mu = -1$.
    // Hardware Justification: Evaluates final metric constraint deviation at the boundary limit.
    normalization_constraint_{PARTICLE}(y, g4DD_local, &norm_final, 1, 0);

    printf("Final Normalization Constraint Evaluation (Massive Particle):\\n");
    printf("  Expected g_mu_nu p^mu p^nu = -1.0\\n");
    printf("  Calculated value (C)       = %.10g\\n", norm_final.C);
    printf("  Absolute Error |C + 1|     = %.10e\\n", fabs(norm_final.C + 1.0));

    printf("Final Conserved Quantities:\\n");
    printf("  E = %.8f, Lz = %.8f, Q = %.8f\\n", cq_final.E, cq_final.Lz, cq_final.Q);

    double E_err = fabs(cq_final.E - cq_init.E); // Absolute error in energy $E$.
    double Lz_err = fabs(cq_final.Lz - cq_init.Lz); // Absolute error in angular momentum $L_z$.
    double Q_err = fabs(cq_final.Q - cq_init.Q); // Absolute error in Carter constant $Q$.

    printf("Conservation Check (Absolute Error):\\n");
    printf("  Delta E  = %.4e\\n", E_err);
    printf("  Delta Lz = %.4e\\n", Lz_err);
    printf("  Delta Q  = %.4e\\n", Q_err);

    // --- Memory Cleanup ---
    gsl_odeiv2_evolve_free(e); // Free GSL evolution object.
    gsl_odeiv2_control_free(c); // Free GSL control object.
    gsl_odeiv2_step_free(s); // Free GSL step object.

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


if __name__ == "__main__":
    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")

    project_name = "mass_geodesic_integrator"
    project_dir = os.path.join("project", project_name)

    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "massive"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"

    if os.path.exists(project_dir):
        shutil.rmtree(project_dir)
    os.makedirs(project_dir, exist_ok=True)

    print(f"Acquiring symbolic data for {GEO_KEY}...")
    metric_data = Analytic_Spacetimes[SPACETIME]
    geodesic_data = Geodesic_Equations[GEO_KEY]

    print("Registering Physics Kernels...")
    register_struct_definitions()

    # 1. Fundamental Tensors
    g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)
    connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)
    calculate_ode_rhs_massive(geodesic_data.geodesic_rhs, metric_data.xx)
    
    if geodesic_data.u0_massive is None:
        raise ValueError(f"u0_massive is None for {GEO_KEY}")
    u0_massive(geodesic_data.u0_massive)

    conserved_quantities(SPACETIME, PARTICLE)
    normalization_constraint(geodesic_data.norm_constraint_expr, PARTICLE)
    ode_gsl_wrapper_massive(SPACETIME)

    main_c(SPACETIME, PARTICLE)

    print("Generating Header Files...")
    CPs.write_CodeParameters_h_files(set_commondata_only=True, project_dir=project_dir)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()
    cmdline_input_and_parfiles.generate_default_parfile(project_dir=project_dir, project_name=project_name)
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(project_name=project_name)


    # Required macros for GSL array indexing
    macro_defs = """
    // Ensure hardware-agnostic array indexing for standalone GSL test
    #ifndef IDX_LOCAL
    #define IDX_LOCAL(component, batch_id, batch_size) ((component) * (batch_size) + (batch_id))
    #endif

    #ifndef IDX_GLOBAL
    #define IDX_GLOBAL(component, ray_id, num_rays) ((component) * (num_rays) + (ray_id))
    #endif
    
    // Add BUNDLE_CAPACITY to define chunk size for memory allocation
    #ifndef BUNDLE_CAPACITY
    #define BUNDLE_CAPACITY 1 
    #endif
    """
    Bdefines_h.register_BHaH_defines("gpu_batch_macros", macro_defs)

    additional_includes = [
        "gsl/gsl_vector.h",
        "gsl/gsl_matrix.h",
        "gsl/gsl_odeiv2.h",
        "gsl/gsl_errno.h",
        "gsl/gsl_math.h",
    ]

    cpu_macros = {
        "BHAH_MALLOC_DEVICE(a, sz)": "#define BHAH_MALLOC_DEVICE(a, sz) BHAH_MALLOC(a, sz)",
        "BHAH_FREE_DEVICE(a)": "#define BHAH_FREE_DEVICE(a) BHAH_FREE(a)",
        "BHAH_HD_INLINE": "#define BHAH_HD_INLINE", 
    }

    Bdefines_h.output_BHaH_defines_h(
        project_dir=project_dir,
        additional_includes=additional_includes,
        enable_rfm_precompute=False,
        supplemental_defines_dict=cpu_macros
    )

    print("Generating Makefile...")
    addl_cflags = ["$(shell gsl-config --cflags)"]
    addl_libs = ["$(shell gsl-config --libs)"]

    Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=project_name,
        compiler_opt_option="fast",
        addl_CFLAGS=addl_cflags,
        addl_libraries=addl_libs,
        CC="gcc",
        src_code_file_ext="c"
    )

    print(" -> Patching Makefile for Windows compatibility...")
    local_tmp_path = "tmp"
    os.makedirs(os.path.join(project_dir, local_tmp_path), exist_ok=True)

    makefile_path = os.path.join(project_dir, "Makefile")
    with open(makefile_path, "r", encoding="utf-8") as f:
        content = f.read()

    with open(makefile_path, "w", encoding="utf-8") as f:
        f.write(f"export TMPDIR = $(CURDIR)/{local_tmp_path}\n")
        f.write(f"export TMP = $(CURDIR)/{local_tmp_path}\n")
        f.write(f"export TEMP = $(CURDIR)/{local_tmp_path}\n")
        f.write("\n")
        f.write(content)

    print("-" * 50)
    print(f"Project generated successfully in {project_dir}")

    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")