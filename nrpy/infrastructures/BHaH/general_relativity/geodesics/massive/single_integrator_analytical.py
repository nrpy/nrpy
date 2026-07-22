# nrpy/infrastructures/BHaH/general_relativity/geodesics/massive/single_integrator_analytical.py
"""
Define the reusable single-particle analytical geodesic integrator orchestrator.

This module registers the C function that evolves one massive test particle in an
analytic spacetime using the GNU Scientific Library Runge-Kutta-Fehlberg 4(5)
integrator. The registered function writes trajectory samples and reports
normalization and conserved-quantity diagnostics while preserving the Structure of
Arrays layout expected by the shared geodesic kernels.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par


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
        double *proper_time;  // Current proper time $\tau$ for the trajectory.
    } PhotonStateSoA;
    """
    Bdefines_h.register_BHaH_defines("PhotonStateSoA", particle_soa_def)


def single_integrator_analytical(spacetime: str, particle: str) -> None:
    """
    Register the single-particle analytical geodesic integrator C function.

    The generated C function initializes one massive-particle trajectory, computes
    its normalized temporal four-velocity, advances the state with the GSL RKF45
    integrator, writes trajectory samples, and reports final normalization and
    conserved-quantity diagnostics.

    :param spacetime: The background spacetime descriptor.
    :param particle: The test-particle type.
    """
    register_struct_definitions()

    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "initial_proper_time",
            "initial_h",
            "initial_p_x",
            "initial_p_y",
            "initial_p_z",
            "initial_t",
            "initial_x",
            "initial_y",
            "initial_z",
            "evolution_measure_max",
            "r_escape",
        ],
        [
            0.0,
            0.1,
            -0.5641,
            0.0,
            0.0,
            0.0,
            4.0123,
            0.0,
            0.0,
            1000.0,
            150.0,
        ],
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameter(
        "int",
        __name__,
        "max_steps",
        200000,
        commondata=True,
        add_to_parfile=True,
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
    Bdefines_h.register_BHaH_defines("single_particle_macros", macro_defs)

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

    desc = """Integrate one massive-particle geodesic with the GSL RKF45 solver.

Initializes a single massive-particle trajectory, computes the initial
four-velocity normalization, advances the state with the GSL RKF45
integrator, writes trajectory samples, and reports final normalization and
conserved-quantity diagnostics.

@return EXIT_SUCCESS on successful completion; EXIT_FAILURE if setup fails.
"""

    cfunc_type = "int"
    name = "single_integrator_analytical"
    params = "int argc, const char *argv[]"

    body = rf"""
    // ==========================================
    // STRUCTURAL SETUP & PARAMETERS
    // ==========================================
    commondata_struct commondata;
    commondata_struct_set_to_default(&commondata);
    cmdline_input_and_parfile_parser(&commondata, argc, argv);

    printf("Starting Massive Geodesic Integrator...\n");
    printf("spacetime: {spacetime}, M=%.2f, a=%.2f\n", commondata.M_scale, commondata.a_spin);

    // ==========================================
    // GLOBAL MEMORY & INITIAL CONDITIONS
    // ==========================================
    const long int num_rays = 1;

    double y[8];
    y[0] = commondata.initial_t;
    y[1] = commondata.initial_x;
    y[2] = commondata.initial_y;
    y[3] = commondata.initial_z;

    y[5] = commondata.initial_p_x;
    y[6] = commondata.initial_p_y;
    y[7] = commondata.initial_p_z;

    double proper_time = commondata.initial_proper_time;

    PhotonStateSoA all_particles;
    all_particles.f = y;
    all_particles.proper_time = &proper_time;

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
      if (evolve != NULL)
        gsl_odeiv2_evolve_free(evolve);
      if (control != NULL)
        gsl_odeiv2_control_free(control);
      if (step != NULL)
        gsl_odeiv2_step_free(step);
      fprintf(stderr, "Error: failed to allocate GSL integrator objects.\n");
      return EXIT_FAILURE;
    }} // END IF: GSL integrator-object allocation failed

    gsl_odeiv2_system sys = {{ode_gsl_wrapper_massive_{spacetime}, NULL, 8, &commondata}};

    const double proper_time_max = 20000.0;
    double h = commondata.initial_h;

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
    }} // END IF: trajectory output file could not be opened
    fprintf(fp, "# proper_time t x y z u^t u^x u^y u^z\n");

    // ==========================================
    // INTEGRATION LOOP
    // ==========================================
    int steps = 0;

    while (proper_time < proper_time_max && steps < commondata.max_steps) {{
      fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
              proper_time, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);

      // Advance the state by one adaptive GSL step.
      int status = gsl_odeiv2_evolve_apply(
          evolve, control, step, &sys, &proper_time, proper_time_max, &h, y
      );

      if (status != GSL_SUCCESS) {{
        printf("GSL Error: %d\n", status);
        break;
      }} // END IF: GSL integration step failed

      if (fabs(y[4]) > commondata.evolution_measure_max) {{
        printf("Temporal four-velocity component |u^t| exceeded numerical limit.\n");
        break;
      }} // END IF: temporal four-velocity exceeded the configured bound

      // Squared Cartesian radius from the origin.
      const double r_squared = y[1] * y[1] + y[2] * y[2] + y[3] * y[3];
      if (r_squared > commondata.r_escape * commondata.r_escape) {{
        printf("Termination: particle reached r_squared = %.4f at proper_time = %.4f\n", r_squared, proper_time);
        break;
      }} // END IF: particle crossed the escape radius
      steps++;
    }} // END WHILE: integrate massive particle geodesic

    fclose(fp);
    printf("Integration finished after %d steps. Final proper_time = %.4f\n", steps, proper_time);

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
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
