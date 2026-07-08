# nrpy/examples/photon_single_geodesic_integrator_analytical.py
r"""
Generate a standalone C project for integrating a single photon geodesic.

The generated project evolves one massless test particle in an analytic
spacetime using the split RKF45 photon pipeline. It writes trajectory samples
and reports normalization and conserved-quantity diagnostics while preserving
the Structure of Arrays layout expected by the shared geodesic kernels.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import os
import shutil

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
    Analytic_Spacetimes,
)
from nrpy.equations.general_relativity.geodesics.geodesics import Geodesic_Equations
from nrpy.infrastructures.BHaH import (
    BHaH_defines_h,
)
from nrpy.infrastructures.BHaH import CodeParameters as CPs
from nrpy.infrastructures.BHaH import Makefile_helpers as Makefile
from nrpy.infrastructures.BHaH import (
    cmdline_input_and_parfiles,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics import (
    connections,
    conserved_quantities,
    g4DD_metric,
    normalization_constraint,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
    calculate_ode_rhs_kernel,
    interpolation_kernel,
    p0_reverse_kernel,
    rkf45_finalize_and_control_kernel,
    rkf45_stage_update,
)


def register_struct_definitions() -> None:
    """Register the termination enum and shared photon-state SoA definition."""
    termination_enum_def = r"""
    // Compatibility names retained because shared photon kernels reference
    // termination_type_t together with ACTIVE/REJECTED failure states.
    typedef enum {
      ACTIVE = 0,
      REJECTED,
      FAILURE_RKF45_REJECTION_LIMIT
    } termination_type_t; // END ENUM: termination_type_t
    """
    BHaH_defines_h.register_BHaH_defines("termination_type_t", termination_enum_def)

    photon_soa_def = r"""
    // Compatibility name retained because shared geodesic kernels expect PhotonStateSoA.
    typedef struct {
      double *f;            // Flattened state vector: t, x, y, z, p_t, p_x, p_y, p_z, aux.
      double *affine_param; // Current affine parameter lambda.
    } PhotonStateSoA; // END STRUCT: PhotonStateSoA
    """
    BHaH_defines_h.register_BHaH_defines("PhotonStateSoA", photon_soa_def)


def register_CFunction_photon_single_main(spacetime: str, particle: str) -> None:
    """
    Register the generated C main() function for the photon geodesic integrator.

    :param spacetime: The specific background spacetime descriptor.
    :param particle: The type of test particle being integrated.
    """
    includes = [
        "math.h",
        "stdio.h",
        "stdlib.h",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]

    desc = """Drive the single-photon geodesic integration example.

Initializes one photon state, evolves it with the split RKF45 pipeline,
writes trajectory samples, and reports final normalization and
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
    commondata.a_spin = 0.0;

    const double r_escape = 150.0;
    const double p_t_max = 1000.0;
    const long int num_rays = 1;
    const long int chunk_size = 1;
    const int stream_idx = 0;
    const int max_steps = 200000;
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
    double *affine_param = NULL;
    double *h = NULL;
    int *rejection_retries = NULL;
    termination_type_t *status = NULL;

    BHAH_MALLOC(f, 9 * sizeof(double));
    BHAH_MALLOC(f_base, 9 * sizeof(double));
    BHAH_MALLOC(f_temp, 9 * sizeof(double));
    BHAH_MALLOC(metric, 10 * sizeof(double));
    BHAH_MALLOC(connection, 40 * sizeof(double));
    BHAH_MALLOC(k_bundle, 6 * 9 * sizeof(double));
    BHAH_MALLOC(affine_param, sizeof(double));
    BHAH_MALLOC(h, sizeof(double));
    BHAH_MALLOC(rejection_retries, sizeof(int));
    BHAH_MALLOC(status, sizeof(termination_type_t));

    if (f == NULL || f_base == NULL || f_temp == NULL || metric == NULL ||
        connection == NULL || k_bundle == NULL || affine_param == NULL ||
        h == NULL || rejection_retries == NULL || status == NULL) {{
      fprintf(stderr, "Error: failed to allocate photon state buffers.\n");
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: photon buffer allocation failed

    PhotonStateSoA all_photons;
    all_photons.f = f;
    all_photons.affine_param = affine_param;

    // ==========================================
    // INITIAL CONDITIONS
    // ==========================================
    f[0] = 0.0;
    f[1] = 4.0123;
    f[2] = 0.0;
    f[3] = 0.0;

    f[5] = -0.5641;
    f[6] = 0.0;
    f[7] = 0.0;
    f[8] = 0.0;

    *affine_param = 0.0;
    *h = 0.1;
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
    fprintf(fp, "# lambda t x y z p_t p_x p_y p_z aux\n");

    // ==========================================
    // MODULAR SPLIT-PIPELINE INTEGRATION LOOP
    // ==========================================
    int steps = 0;

    while (steps < max_steps) {{
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
        affine_param,
        rejection_retries,
        chunk_size,
        stream_idx
      );

      if (*rejection_retries == 0) {{
        fprintf(
          fp,
          "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
          *affine_param,
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
      if (r2 > r_escape * r_escape) {{
        printf("Particle escaped to r > %.2f.\n", r_escape);
        break;
      }} // END IF: particle crossed the escape radius

      if (fabs(f[4]) > p_t_max) {{
        printf("Temporal momentum |p_t| exceeded numerical limit.\n");
        break;
      }} // END IF: temporal momentum exceeded the configured bound

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
      *affine_param
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
    BHAH_FREE(affine_param);
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


if __name__ == "__main__":
    # Step 1: Select the BHaH OpenMP backend for this example.
    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")

    # Step 2: Define project-level constants and output paths.
    project_name = "photon_single_geodesic_integrator_analytical"
    project_dir = os.path.abspath(os.path.join("project", project_name))
    parfile_path = os.path.join(project_dir, f"{project_name}.par")

    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"

    # Step 3: Recreate the generated project directory.
    shutil.rmtree(project_dir, ignore_errors=True)
    os.makedirs(project_dir, exist_ok=True)

    # Step 4: Acquire symbolic data for the selected spacetime.
    print(f"Acquiring symbolic data for {GEO_KEY}...")
    metric_data = Analytic_Spacetimes[SPACETIME]
    geodesic_data = Geodesic_Equations[GEO_KEY]

    # Step 5: Register the shared physics and diagnostics kernels.
    print("Registering Split-Pipeline Physics Kernels...")
    register_struct_definitions()
    g4DD_metric.g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)
    connections.connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)
    conserved_quantities.conserved_quantities(SPACETIME, PARTICLE)
    normalization_constraint.normalization_constraint(
        geodesic_data.norm_constraint_expr, PARTICLE
    )

    if geodesic_data.p0_photon is None:
        raise ValueError(f"p0_photon is None for {GEO_KEY}")
    p0_reverse_kernel.p0_reverse_kernel(geodesic_data.p0_photon)
    interpolation_kernel.interpolation_kernel(SPACETIME)
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(
        geodesic_data.geodesic_rhs, geodesic_data.xx
    )
    rkf45_stage_update.rkf45_stage_update()
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel()

    # Step 5.a: Register the single-ray C main function.
    register_CFunction_photon_single_main(SPACETIME, PARTICLE)

    # Step 5.b: Remove helper registrations emitted only through other kernels.
    # The generated main() calls the metric and connection logic through shared
    # kernels, so emitting standalone wrappers here would add redundant source
    # files and prototypes.
    for internal_func in [f"g4DD_metric_{SPACETIME}", f"connections_{SPACETIME}"]:
        cfc.CFunction_dict.pop(internal_func, None)

    # Step 6: Set relevant CodeParameter defaults.
    par.glb_code_params_dict["rkf45_absolute_error_tolerance"].defaultvalue = 1e-17
    par.glb_code_params_dict["rkf45_error_tolerance"].defaultvalue = 1e-17
    par.glb_code_params_dict["rkf45_h_max"].defaultvalue = 10.0
    par.glb_code_params_dict["rkf45_h_min"].defaultvalue = 1e-20
    par.glb_code_params_dict["rkf45_max_retries"].defaultvalue = 15

    # Step 7: Generate headers, default parameters, and the Makefile.
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
    #ifndef BUNDLE_CAPACITY
    #define BUNDLE_CAPACITY 1
    #endif
    """
    BHaH_defines_h.register_BHaH_defines("single_photon_macros", macro_defs)

    # BHAH_MALLOC expects an assignable lvalue, so keep direct forwarding in
    # the compatibility wrappers below.
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
            "// Commondata is passed by reference in this OpenMP example.\n"
        ),
        "BHAH_DEVICE_SYNC()": "#define BHAH_DEVICE_SYNC() do {} while(0)\n",
    }

    BHaH_defines_h.output_BHaH_defines_h(
        project_dir=project_dir,
        enable_rfm_precompute=False,
        supplemental_defines_dict=cpu_macros,
    )

    addl_cflags = [
        "-fopenmp",
        "-O3",
        "-DDEBUG",
        "-Wno-stringop-truncation",
    ]
    addl_libs = ["-lm"]
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

    # Step 8: Copy the trajectory-visualization helper and print usage.
    vis_dir = os.path.join("nrpy", "examples", "geodesic_visualizations")
    vis_script_src = os.path.join(vis_dir, "visualize_trajectory.py")

    if os.path.exists(vis_script_src):
        shutil.copy(vis_script_src, project_dir)
    else:
        print(
            f"Warning: Visualization script not found at {vis_script_src}. "
            "Please ensure it exists."
        )

    print(
        f"Finished! Now go into {project_dir} and type `make` to build, "
        f"then ./{project_name} to run."
    )
    print(f"    Parameter file can be found at {parfile_path}\n")
    print(
        "    To generate the trajectory plot after running the C executable, "
        "ensure you have the required Python packages:"
    )
    print("    pip install matplotlib numpy\n")
    print(
        "    Then, execute the visualization script directly from the "
        "project directory:"
    )
    print("    python3 visualize_trajectory.py --particle_type Photon\n")
