# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_analytical.py
"""
Define the reusable single-photon analytical geodesic integrator orchestrator.

This module registers the C function that evolves one massless test particle in an
analytic spacetime using the split Runge-Kutta-Fehlberg 4(5) photon pipeline. The
registered function writes trajectory samples and reports normalization and
conserved-quantity diagnostics while preserving the Structure of Arrays layout
expected by the shared geodesic kernels.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.BHaH import BHaH_defines_h


def register_struct_definitions() -> None:
    """Register definitions needed only by the single-photon integrator."""
    termination_enum_def = r"""
    typedef enum {
      ACTIVE = 0,
      REJECTED,
      FAILURE_RKF45_REJECTION_LIMIT
    } termination_type_t; // END ENUM: termination_type_t
    """
    BHaH_defines_h.register_BHaH_defines("termination_type_t", termination_enum_def)

    photon_soa_def = r"""
    // Single-ray compatibility layout used by conserved-quantity diagnostics.
    typedef struct {
      double *f;            // State: t, x, y, z, p^t, p^x, p^y, p^z, Eulerian distance.
      double *integration_param; // Direct-mode affine integration parameter.
    } PhotonStateSoA; // END STRUCT: PhotonStateSoA
    """
    BHaH_defines_h.register_BHaH_defines("PhotonStateSoA", photon_soa_def)

    macro_defs = r"""
    #ifndef BUNDLE_CAPACITY
    #define BUNDLE_CAPACITY 1
    #endif
    """
    BHaH_defines_h.register_BHaH_defines("single_photon_macros", macro_defs)


def single_integrator_analytical(
    spacetime: str, particle: str, normalized_eom: bool = False
) -> None:
    """
    Register the single-photon analytical geodesic integrator C function.

    The generated C function initializes one photon state, evolves it with the
    split RKF45 pipeline, writes trajectory samples, and reports final
    normalization and conserved-quantity diagnostics.

    :param spacetime: The background spacetime descriptor.
    :param particle: The test-particle type.
    :param normalized_eom: Whether to use normalized photon evolution.
    :raises ValueError: If normalized photon evolution is requested.
    """
    if normalized_eom:
        raise ValueError(
            "single_integrator_analytical supports direct geodesic evolution only."
        )

    register_struct_definitions()
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "initial_t",
            "initial_x",
            "initial_y",
            "initial_z",
            "initial_p_x",
            "initial_p_y",
            "initial_p_z",
            "initial_integration_param",
            "initial_eulerian_distance",
            "initial_h",
            "r_escape",
            "energy_max",
        ],
        [0.0] * 12,
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameter(
        "int",
        __name__,
        "max_steps",
        1,
        commondata=True,
        add_to_parfile=True,
    )

    includes = [
        "math.h",
        "stdio.h",
        "stdlib.h",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]

    desc = """Integrate one photon geodesic with the split RKF45 pipeline.

Initializes one photon state, evolves it with the split RKF45 pipeline,
writes trajectory samples, and reports final normalization and
conserved-quantity diagnostics.

@return EXIT_SUCCESS on successful completion; EXIT_FAILURE if setup fails.
"""

    cfunc_type = "int"
    name = "single_integrator_analytical"
    params = "int argc, const char *argv[]"

    body = rf"""
    (void)argc;
    (void)argv;

    // ==========================================
    // STRUCTURAL SETUP & PARAMETERS
    // ==========================================
    commondata_struct commondata;
    commondata_struct_set_to_default(&commondata);

    const long int num_rays = 1;
    const long int chunk_size = 1;
    const int stream_idx = 0;
    const char *status_names[] = {{
      "ACTIVE",
      "REJECTED",
      "FAILURE_RKF45_REJECTION_LIMIT"
    }};

    int exit_status = EXIT_SUCCESS;
    FILE *fp = NULL;

    printf("Starting Split-Pipeline Geodesic Integrator...\n");
    printf("spacetime: {spacetime}, M=%.2f, a=%.2f\n", commondata.M_scale, commondata.a_spin);

    // ==========================================
    // GLOBAL MEMORY ALLOCATION
    // ==========================================
    double *f = NULL;
    double *f_base = NULL;
    double *f_temp = NULL;
    double *metric = NULL;
    double *connection = NULL;
    double *k_bundle = NULL;
    double *integration_param = NULL;
    double *h = NULL;
    int *rejection_retries = NULL;
    termination_type_t *status = NULL;

    BHAH_MALLOC(f, 9 * sizeof(double));
    BHAH_MALLOC(f_base, 9 * sizeof(double));
    BHAH_MALLOC(f_temp, 9 * sizeof(double));
    BHAH_MALLOC(metric, 10 * sizeof(double));
    BHAH_MALLOC(connection, 40 * sizeof(double));
    BHAH_MALLOC(k_bundle, 6 * 9 * sizeof(double));
    BHAH_MALLOC(integration_param, sizeof(double));
    BHAH_MALLOC(h, sizeof(double));
    BHAH_MALLOC(rejection_retries, sizeof(int));
    BHAH_MALLOC(status, sizeof(termination_type_t));

    if (f == NULL || f_base == NULL || f_temp == NULL || metric == NULL ||
        connection == NULL || k_bundle == NULL || integration_param == NULL ||
        h == NULL || rejection_retries == NULL || status == NULL) {{
      fprintf(stderr, "Error: failed to allocate photon state buffers.\n");
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: photon buffer allocation failed

    PhotonStateSoA all_photons;
    all_photons.f = f;
    all_photons.integration_param = integration_param;

    // ==========================================
    // INITIAL CONDITIONS
    // ==========================================
    f[0] = commondata.initial_t;
    f[1] = commondata.initial_x;
    f[2] = commondata.initial_y;
    f[3] = commondata.initial_z;
    f[5] = commondata.initial_p_x;
    f[6] = commondata.initial_p_y;
    f[7] = commondata.initial_p_z;
    f[8] = commondata.initial_eulerian_distance;

    *integration_param = commondata.initial_integration_param;
    *h = commondata.initial_h;
    *rejection_retries = 0;
    *status = ACTIVE;

    // ==========================================
    // INITIAL METRIC EVALUATION
    // ==========================================
    interpolation_kernel_{spacetime}(
      &commondata, f, metric, NULL, chunk_size, stream_idx
    );
    p0_reverse_kernel(f, metric, chunk_size, stream_idx);

    printf("Initial State:\n");
    printf("  Pos (%.4f, %.4f, %.4f)\n", f[1], f[2], f[3]);
    printf("  Mom (%.4f, %.4f, %.4f, %.4f)\n", f[4], f[5], f[6], f[7]);

    // ==========================================
    // PRE-INTEGRATION DIAGNOSTICS
    // ==========================================
    conserved_quantities_t cq_init;
    calculate_conserved_quantities_universal_{spacetime}_{particle}(
      &commondata, &all_photons, num_rays, &cq_init
    );

    fp = fopen("trajectory.txt", "w");
    if (fp == NULL) {{
      fprintf(stderr, "Error: could not open trajectory.txt for writing.\n");
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: trajectory output file could not be opened
    fprintf(fp, "# lambda t x y z energy_measure p_x p_y p_z aux\n");

    // ==========================================
    // MODULAR SPLIT-PIPELINE INTEGRATION LOOP
    // ==========================================
    int steps = 0;

    while (steps < commondata.max_steps) {{
      for (int i = 0; i < 9; i++) {{
        f_base[i] = f[i];
        f_temp[i] = f[i];
      }} // END LOOP: copy the current state into RKF45 work buffers

      for (int stage = 1; stage <= 6; stage++) {{
        interpolation_kernel_{spacetime}(
          &commondata, f_temp, metric, connection, chunk_size, stream_idx
        );
        calculate_ode_rhs_kernel(
          f_temp, metric, connection, k_bundle, stage, chunk_size, stream_idx
        );

        if (stage < 6)
          rkf45_stage_update(
            f, k_bundle, h, stage, chunk_size, f_temp, stream_idx
          );
      }} // END LOOP: execute RKF45 stages 1 through 6

      rkf45_finalize_and_control(
        &commondata,
        f,
        f_base,
        k_bundle,
        h,
        status,
        integration_param,
        rejection_retries,
        chunk_size,
        stream_idx
      );

      if (*rejection_retries == 0) {{
        fprintf(
          fp,
          "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
          *integration_param,
          f[0],
          f[1],
          f[2],
          f[3],
          f[4],
          f[5],
          f[6],
          f[7],
          f[8]
        );
        steps++;
      }} // END IF: accepted step written to the trajectory file

      const double r2 = f[1] * f[1] + f[2] * f[2] + f[3] * f[3];
      if (r2 > commondata.r_escape * commondata.r_escape) {{
        printf("Particle escaped to r > %.2f.\n", commondata.r_escape);
        break;
      }} // END IF: particle crossed the escape radius

      if (fabs(f[4]) > commondata.energy_max) {{
        printf("Energy measure exceeded numerical limit.\n");
        break;
      }} // END IF: energy measure exceeded the configured bound

      if (*status == FAILURE_RKF45_REJECTION_LIMIT) {{
        printf(
          "Integration terminated with status: %s (%d)\n",
          status_names[*status],
          *status
        );
        break;
      }} // END IF: RKF45 rejection limit was reached
    }} // END WHILE: integrate the photon geodesic

    printf(
      "Integration finished after %d steps. Final lambda = %.4f\n",
      steps,
      *integration_param
    );

    // ==========================================
    // POST-INTEGRATION DIAGNOSTICS
    // ==========================================
    conserved_quantities_t cq_final;
    calculate_conserved_quantities_universal_{spacetime}_{particle}(
      &commondata, &all_photons, num_rays, &cq_final
    );

    normalization_constraint_t norm_final;
    interpolation_kernel_{spacetime}(
      &commondata, f, metric, NULL, chunk_size, stream_idx
    );
    normalization_constraint_photon(
      f, metric, &norm_final, chunk_size, stream_idx
    );

    printf("\nFinal normalization constraint: |C| = %.4e\n", fabs(norm_final.C));
    printf("Conservation Absolute Errors:\n");
    printf("  Delta E  = %.4e\n", fabs(cq_final.E - cq_init.E));
    printf("  Delta Lz = %.4e\n", fabs(cq_final.Lz - cq_init.Lz));
    printf("  Delta Q  = %.4e\n", fabs(cq_final.Q - cq_init.Q));

    cleanup:
    if (fp != NULL)
      fclose(fp);
    BHAH_FREE(f);
    BHAH_FREE(f_base);
    BHAH_FREE(f_temp);
    BHAH_FREE(metric);
    BHAH_FREE(connection);
    BHAH_FREE(k_bundle);
    BHAH_FREE(integration_param);
    BHAH_FREE(h);
    BHAH_FREE(rejection_retries);
    BHAH_FREE(status);

    return exit_status;
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )

    par.glb_code_params_dict["rkf45_absolute_error_tolerance"].defaultvalue = 1e-17
    par.glb_code_params_dict["rkf45_error_tolerance"].defaultvalue = 1e-17
    par.glb_code_params_dict["rkf45_h_max"].defaultvalue = 10.0
    par.glb_code_params_dict["rkf45_h_min"].defaultvalue = 1e-20
    par.glb_code_params_dict["rkf45_max_retries"].defaultvalue = 15


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
