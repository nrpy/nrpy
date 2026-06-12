# nrpy/examples/massive_single_geodesic_integrator_analytical.py
r"""
Defines a C project for integrating massive geodesics in curved spacetime.

The script provides a standalone C application that evolves the trajectory
of a massive test particle. It coordinates the registration of spacetime-specific
physics kernels and links them with the GNU Scientific Library for time integration.
Structure of Arrays (SoA) layouts provide universal memory compatibility for shared
physics kernels. Struct mapping establishes the layout required by downstream conserved
quantities kernels. The initial metric is evaluated at the starting position so that
u^t can be computed from the massive-particle normalization condition. The final
normalization residual is reported at the last trajectory point.

The simulation solves the geodesic equation for a particle with non-zero mass:
    $d(u^\mu)/d(\tau) = -\Gamma^\mu_{\alpha \beta} u^\alpha u^\beta$
subject to the normalization constraint $u^\mu u_\mu = -1$.

Numerical fidelity is validated by monitoring constants of motion
associated with the spacetime's symmetries, including Killing vectors and tensors.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import os
import shutil

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
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
    """Register the termination enum and shared particle-state SoA definition."""
    termination_enum_def = r"""
    // Defines the specific exit condition for a particle's integration loop.
    typedef enum {
        ACTIVE = 0,                        // 0: Particle is currently undergoing integration.
        REJECTED,                          // 1: Particle RKF45 step was rejected.
        FAILURE_RKF45_REJECTION_LIMIT      // 2: Adaptive step-size rejected too many consecutive times.
    } termination_type_t;
    """
    Bdefines_h.register_BHaH_defines("termination_type_t", termination_enum_def)

    particle_soa_def = r"""
    // ==========================================
    // Flattened SoA Struct (Master Storage)
    // ==========================================
    // Compatibility name retained because shared geodesic kernels expect PhotonStateSoA.
    typedef struct {
        double *f;            // Flattened state vector: t, x, y, z, u^t, u^x, u^y, u^z.
        double *affine_param; // Current proper time $\tau$ for the trajectory.
    } PhotonStateSoA;
    """
    Bdefines_h.register_BHaH_defines("PhotonStateSoA", particle_soa_def)


def register_CFunction_main_c(spacetime: str, particle: str) -> None:
    """
    Register the generated C main() function for the massive geodesic integrator.

    :param spacetime: The specific background spacetime descriptor.
    :param particle: The type of test particle being integrated.
    """
    includes = [
        "math.h",
        "stdio.h",
        "stdlib.h",
        "gsl/gsl_errno.h",
        "gsl/gsl_math.h",
        "gsl/gsl_matrix.h",
        "gsl/gsl_odeiv2.h",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]

    desc = """Main driver function for the massive geodesic integrator.

Initializes a single massive-particle trajectory, computes the initial
four-velocity normalization, advances the state with the GSL RKF45
integrator, writes trajectory samples, and reports final normalization and
conserved-quantity diagnostics.

@param argc Number of command-line arguments.
@param argv Command-line argument array.
@return EXIT_SUCCESS on successful completion; EXIT_FAILURE if setup fails.
"""

    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"

    body = rf"""
    // ==========================================
    // STRUCTURAL SETUP & PARAMETERS
    // ==========================================
    commondata_struct commondata;
    commondata_struct_set_to_default(&commondata);

    commondata.M_scale = 1.0;
    commondata.a_spin = 0.9;

    printf("Starting Massive Geodesic Integrator...\n");
    printf("spacetime: {spacetime}, M=%.2f, a=%.2f\n", commondata.M_scale, commondata.a_spin);

    // ==========================================
    // GLOBAL MEMORY & INITIAL CONDITIONS
    // ==========================================
    const double u0_abs_max = 1000.0;
    const double r_squared_max = 1e4;
    const long int num_rays = 1;

    double y[8];
    y[0] = 0.0;
    y[1] = 10.0;
    y[2] = 1.0;
    y[3] = 1.0;

    y[5] = -0.1;
    y[6] = 0.33;
    y[7] = 0.0;

    double tau = 0.0;

    PhotonStateSoA all_particles;
    all_particles.f = y;
    all_particles.affine_param = &tau;

    double g4dd_local[10];

    g4DD_metric_{spacetime}(&commondata, y, g4dd_local);

    double u0_val = 0.0;
    u0_massive(g4dd_local, y, 1, 1, 0, 0, &u0_val);
    y[4] = u0_val;

    printf("Initial State:\n");
    printf("  Pos: (%.4f, %.4f, %.4f)\n", y[1], y[2], y[3]);
    printf("  Vel: (%.4f, %.4f, %.4f, %.4f)\n", y[4], y[5], y[6], y[7]);

    // ==========================================
    // PRE-INTEGRATION DIAGNOSTICS
    // ==========================================
    conserved_quantities_t cq_init;
    calculate_conserved_quantities_universal_{spacetime}_{particle}(&commondata, &all_particles, num_rays, &cq_init);

    printf("Initial Conserved Quantities:\n");
    printf("  E = %.8f, Lz = %.8f, Q = %.8f\n", cq_init.E, cq_init.Lz, cq_init.Q);

    // ==========================================
    // GSL INTEGRATOR SETUP
    // ==========================================
    const gsl_odeiv2_step_type *step_type = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(step_type, 8);
    gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(1e-9, 0.0);
    gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(8);

    if (step == NULL || control == NULL || evolve == NULL) {{
      if (evolve != NULL) {{
        gsl_odeiv2_evolve_free(evolve);
      }} // END IF: evolve allocation succeeded
      if (control != NULL) {{
        gsl_odeiv2_control_free(control);
      }} // END IF: control allocation succeeded
      if (step != NULL) {{
        gsl_odeiv2_step_free(step);
      }} // END IF: step allocation succeeded
      fprintf(stderr, "Error: failed to allocate GSL integrator objects.\n");
      return EXIT_FAILURE;
    }} // END IF: GSL allocation failed

    gsl_odeiv2_system sys = {{ode_gsl_wrapper_massive_{spacetime}, NULL, 8, &commondata}};

    const double tau_max = 20000.0;
    double h = 1e-3;

    // ==========================================
    // FILE OUTPUT SETUP
    // ==========================================
    FILE *fp = fopen("trajectory.txt", "w");
    if (fp == NULL) {{
      gsl_odeiv2_evolve_free(evolve);
      gsl_odeiv2_control_free(control);
      gsl_odeiv2_step_free(step);
      fprintf(stderr, "Error opening trajectory.txt\n");
      return EXIT_FAILURE;
    }}
    fprintf(fp, "# tau t x y z u^t u^x u^y u^z\n");

    // ==========================================
    // INTEGRATION LOOP
    // ==========================================
    int steps = 0;
    const int max_steps = 2000000;

    while (tau < tau_max && steps < max_steps) {{
        fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
                tau, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);

        // Advance the state by one adaptive GSL step.
        int status = gsl_odeiv2_evolve_apply(
            evolve, control, step, &sys, &tau, tau_max, &h, y
        );

        if (status != GSL_SUCCESS) {{
            printf("GSL Error: %d\n", status);
            break;
        }}

        if (fabs(y[4]) > u0_abs_max) {{
            printf("Temporal four-velocity component |u^t| exceeded numerical limit.\n");
            break;
        }}

        // Squared Cartesian radius from the origin.
        const double r_squared = y[1] * y[1] + y[2] * y[2] + y[3] * y[3];
        if (r_squared > r_squared_max) {{
            printf("Termination: particle reached r_squared = %.4f at tau = %.4f\n", r_squared, tau);
            break;
        }}
        steps++;
    }} // END WHILE: integrate massive particle geodesic

    fclose(fp);
    printf("Integration finished after %d steps. Final tau = %.4f\n", steps, tau);

    // ==========================================
    // POST-INTEGRATION DIAGNOSTICS
    // ==========================================
    conserved_quantities_t cq_final;
    calculate_conserved_quantities_universal_{spacetime}_{particle}(&commondata, &all_particles, num_rays, &cq_final);

    g4DD_metric_{spacetime}(&commondata, y, g4dd_local);

    normalization_constraint_t norm_final;
    normalization_constraint_{particle}(y, g4dd_local, &norm_final, 1, 0);

    printf("Final Normalization Constraint Evaluation (Massive particle):\n");
    printf("  Expected g_mu_nu u^mu u^nu = -1.0\n");
    printf("  Calculated value (C)       = %.10g\n", norm_final.C);
    printf("  Absolute Error |C + 1|     = %.10e\n", fabs(norm_final.C + 1.0));

    printf("Final Conserved Quantities:\n");
    printf("  E = %.8f, Lz = %.8f, Q = %.8f\n", cq_final.E, cq_final.Lz, cq_final.Q);

    const double energy_err = fabs(cq_final.E - cq_init.E);
    const double lz_err = fabs(cq_final.Lz - cq_init.Lz);
    const double carter_q_err = fabs(cq_final.Q - cq_init.Q);

    printf("Conservation Check (Absolute Error):\n");
    printf("  Delta E  = %.4e\n", energy_err);
    printf("  Delta Lz = %.4e\n", lz_err);
    printf("  Delta Q  = %.4e\n", carter_q_err);

    // ==========================================
    // MEMORY CLEANUP
    // ==========================================
    gsl_odeiv2_evolve_free(evolve);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);

    return EXIT_SUCCESS;
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
    # Step P1: Set code-generation parameters and register the geodesic kernels.
    enable_parallel_codegen = True
    if enable_parallel_codegen:
        pcg.do_parallel_codegen()

    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")

    project_name = "massive_single_geodesic_integrator_analytical"
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

    g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)
    connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)
    calculate_ode_rhs_massive(geodesic_data.geodesic_rhs, metric_data.xx)

    if geodesic_data.u0_massive is None:
        raise ValueError(f"u0_massive is None for {GEO_KEY}")
    u0_massive(geodesic_data.u0_massive)

    conserved_quantities(SPACETIME, PARTICLE)
    normalization_constraint(geodesic_data.norm_constraint_expr, PARTICLE)
    ode_gsl_wrapper_massive(SPACETIME)

    register_CFunction_main_c(SPACETIME, PARTICLE)

    print("Generating header files and Makefile...")

    CPs.write_CodeParameters_h_files(set_commondata_only=True, project_dir=project_dir)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()

    cmdline_input_and_parfiles.generate_default_parfile(
        project_dir=project_dir, project_name=project_name
    )
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
        project_name=project_name
    )

    macro_defs = r"""
    // Provide standalone indexing helpers used by the shared geodesic kernels.
    #ifndef IDX_LOCAL
    #define IDX_LOCAL(component, batch_id, batch_size) ((component) * (batch_size) + (batch_id))
    #endif

    #ifndef IDX_GLOBAL
    #define IDX_GLOBAL(component, ray_id, num_rays) ((component) * (num_rays) + (ray_id))
    #endif

    // Single-particle GSL driver: set BUNDLE_CAPACITY to 1 for shared memory-layout macros.
    #ifndef BUNDLE_CAPACITY
    #define BUNDLE_CAPACITY 1
    #endif
    """
    Bdefines_h.register_BHaH_defines("gpu_batch_macros", macro_defs)

    cpu_macros = {
        "BHAH_MALLOC_DEVICE(a, sz)": "#define BHAH_MALLOC_DEVICE(a, sz) BHAH_MALLOC((a), (sz))",
        "BHAH_FREE_DEVICE(a)": "#define BHAH_FREE_DEVICE(a) BHAH_FREE((a))",
        "BHAH_HD_INLINE": "#define BHAH_HD_INLINE",
    }

    additional_includes = [
        "gsl/gsl_vector.h",
        "gsl/gsl_matrix.h",
        "gsl/gsl_odeiv2.h",
        "gsl/gsl_errno.h",
        "gsl/gsl_math.h",
    ]

    Bdefines_h.output_BHaH_defines_h(
        project_dir=project_dir,
        additional_includes=additional_includes,
        enable_rfm_precompute=False,
        supplemental_defines_dict=cpu_macros,
    )

    addl_cflags = ["$(shell gsl-config --cflags)", "-Wno-stringop-truncation"]
    addl_libs = ["$(shell gsl-config --libs)"]

    Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=project_name,
        compiler_opt_option="fast",
        addl_CFLAGS=addl_cflags,
        addl_libraries=addl_libs,
        CC="gcc",
        src_code_file_ext="c",
    )

    # Step P2: Copy the trajectory visualizer into the generated project.
    vis_dir = os.path.join("nrpy", "examples", "geodesic_visualizations")
    vis_script_src = os.path.join(vis_dir, "visualize_trajectory.py")

    if os.path.exists(vis_script_src):
        shutil.copy(vis_script_src, project_dir)
    else:
        print(
            f"Warning: Visualization script not found at {vis_script_src}; trajectory plotting script was not copied."
        )

    print(
        f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
    )
    print(f"    Parameter file can be found in {project_name}.par\n")
    print(
        "    To generate the trajectory plot after running the C executable, ensure you have the required Python packages:"
    )
    print("    pip install matplotlib numpy\n")
    print(
        "    Then, execute the visualization script directly from the project directory:"
    )
    print("    python3 visualize_trajectory.py --particle_type Massive\n")
