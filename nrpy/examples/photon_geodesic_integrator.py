r"""
Construct a complete C project for integrating photon geodesics in curved spacetime.

Project: NRPy Standalone Geodesic Integrator (Split-Pipeline CPU)
Description:
    This module is a standalone C application that evolves the trajectory
    of a massless photon test particle.

Physics Context:
    The simulation solves the geodesic equation for a photon:
        $d(p^\mu)/d(\lambda) = -\Gamma^\mu_{\alpha \beta} p^\alpha p^\beta$
    subject to the normalization constraint $p^\mu p_\mu = 0$.

    Numerical fidelity is rigorously validated by monitoring constants of motion
    associated with the spacetime's symmetries (Killing vectors and tensors).

Author: Dalton J. Moone.
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

# Split-Pipeline Physics and Math Kernels
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
    """Register the Structure of Arrays (SoA) and diagnostic structures into the global C headers."""
    termination_enum_def = """
    // Defines the specific exit condition for a photon's integration loop.
    typedef enum {
        ACTIVE = 0,                        // 0: Photon is currently undergoing integration.
        REJECTED,                          // 1: Photon RKF45 step was rejected.
        FAILURE_RKF45_REJECTION_LIMIT      // 2: Adaptive step-size rejected too many consecutive times.
    } termination_type_t;
    """
    Bdefines_h.register_BHaH_defines("termination_type_t", termination_enum_def)

    photon_soa_def = r"""
    // ==========================================
    // Flattened SoA Struct (Master Storage)
    // ==========================================
    // Hardware Justification: This Structure of Arrays (SoA) minimizes memory divergence during parallel execution.
    typedef struct {
        double *f; // Flattened state vector mapping 9 components $t, x, y, z, p_t, p_x, p_y, p_z, \text{aux}$.
        double *affine_param; // Current affine parameter $\lambda$ for the trajectory.
    } PhotonStateSoA;
    """
    Bdefines_h.register_BHaH_defines("PhotonStateSoA", photon_soa_def)


def main_c(spacetime: str, particle: str) -> None:
    """
    Generate the main() function orchestrating the single-ray Split-Pipeline RKF45 integrator.

    :param spacetime: The specific background spacetime descriptor (e.g., KerrSchild_Cartesian).
    :param particle: The type of test particle being integrated (e.g., photon).
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "string.h", "stdlib.h"]

    desc = """@brief Main driver function for the Split-Pipeline photon geodesic integrator.
    Detailed algorithm: Initializes memory for a single ray using SoA layouts and executes the discrete
    RKF45 stages via isolated kernels to prevent register spilling. Evaluates boundary constraints
    and verifies the numerical fidelity using Hamiltonian constraints and conserved quantities."""

    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"

    body = f"""
    // --- Step 1: Structural Setup & Parameters ---
    commondata_struct commondata; // Global parameters struct governing spacetime and numerics.
    commondata_struct_set_to_default(&commondata);
    commondata.M_scale = 1.0; // The mass $M$ of the central black hole.
    commondata.a_spin = 0.9; // The dimensionless spin $a$ of the central black hole.

    double r_escape = 150.0; // The threshold radial boundary $r_{{escape}}$ for exiting the simulation.
    double p_t_max = 1000.0; // The upper limit for the temporal momentum $p_t$ to prevent numerical blowup.

    printf("Starting Split-Pipeline Geodesic Integrator (CPU)...\\n");
    printf("spacetime: {spacetime}, M=%.2f, a=%.2f\\n", commondata.M_scale, commondata.a_spin);

    // --- Step 2: Global Memory Allocation (Host) ---
    long int num_rays = 1; // Total number of global photon trajectories.
    long int chunk_size = 1; // Number of active rays in the current batch processing chunk.
    int stream_idx = 0; // Hardware stream index for synchronized execution.

    // Hardware Justification: Memory is explicitly dynamically allocated to mimic the batch VRAM pipeline layout.
    double *f = (double *)malloc(9 * sizeof(double));            // Persistent 9-component state vector $f^\\mu$.
    double *f_base = (double *)malloc(9 * sizeof(double));       // The pristine anchor state $f^\\mu_{{base}}$ preserved across all RKF45 stages.
    double *f_temp = (double *)malloc(9 * sizeof(double));       // Intermediate RKF45 state buffer $f^\\mu_{{temp}}$.
    double *metric = (double *)malloc(10 * sizeof(double)); // The 10 independent components of the symmetric metric $g_{{\\mu\\nu}}$.
    double *connection = (double *)malloc(40 * sizeof(double)); // The 40 independent Christoffel symbols $\\Gamma^\\alpha_{{\\beta\\gamma}}$.
    double *k_bundle = (double *)malloc(6 * 9 * sizeof(double)); // The massive derivative bundle storing all 6 RKF45 $k_n$ stages.

    double *affine_param = (double *)malloc(sizeof(double)); // The affine parameter $\\lambda$ accumulated along the geodesic.
    double *h = (double *)malloc(sizeof(double)); // The adaptive integration step size $h$.
    int *rejection_retries = (int *)malloc(sizeof(int)); // Counter for step rejections due to truncation error limits.
    termination_type_t *status = (termination_type_t *)malloc(sizeof(termination_type_t)); // The integration status flag for the current ray.

    // Hardware Justification: Struct mapping establishes the SoA layout required by the downstream conserved quantities kernel.
    PhotonStateSoA all_photons; // Master struct holding array pointers for the global dataset.
    all_photons.f = f;
    all_photons.affine_param = affine_param;

    // --- Step 3: Initial Conditions ---
    f[0] = 0.0;   // The coordinate time $t$.
    f[1] = 4.0123;  // The spatial coordinate $x$.
    f[2] = 0.0;   // The spatial coordinate $y$.
    f[3] = 0.0;   // The spatial coordinate $z$.

    f[5] = -0.5641;  // The spatial momentum $p_x$.
    f[6] = 0.7171;  // The spatial momentum $p_y$.
    f[7] = 0.0;   // The spatial momentum $p_z$.
    f[8] = 0.0;   // The auxiliary path length integral $L$.

    *affine_param = 0.0; // Initialize affine parameter $\\lambda$ at the source.
    *h = 0.1; // Establish initial step size guess $h$.
    *rejection_retries = 0; // Reset error retries.
    *status = 0; // Set initial condition to ACTIVE.

    // --- Step 4: Metric Interpolation & Temporal Momentum Reversal ---
    // Hardware Justification: Pre-computes the metric to solve the quadratic Hamiltonian constraint for $p_t$.
    interpolation_kernel_{spacetime}(&commondata, f, metric, NULL, chunk_size, stream_idx);
    p0_reverse_kernel(f, metric, chunk_size, stream_idx);

    printf("Initial State:\\n");
    printf("  Pos (%.4f, %.4f, %.4f)\\n", f[1], f[2], f[3]);
    printf("  Mom (%.4f, %.4f, %.4f, %.4f)\\n", f[4], f[5], f[6], f[7]);

    // --- Step 5: Initial Conserved Quantities Extraction ---
    conserved_quantities_t cq_init; // Struct instance to store the baseline constants of motion.
    calculate_conserved_quantities_universal_{spacetime}_{particle}(&commondata, &all_photons, num_rays, &cq_init);

    FILE *fp = fopen("trajectory.txt", "w"); // File pointer for trajectory data offloading.
    if (fp == NULL) {{
        fprintf(stderr, "Error: Could not open trajectory.txt for writing.\\n");
        return 1;
    }}
    fprintf(fp, "# lambda t x y z p_t p_x p_y p_z aux\\n");

    // --- Step 6: The Modular Split-Pipeline Integration Loop ---
    int steps = 0; // Iteration counter tracking successful integration steps.
    int max_steps = 200000; // Hard ceiling preventing runaway numerical loops.

    while (steps < max_steps) {{

        // Hardware Justification: Deep copy preserves the uncorrupted anchor state $f^\\mu_{{base}}$ to prevent read-after-write hazards during the RKF45 tableau evaluation.
        for (int i = 0; i < 9; i++) {{
        f_base[i] = f[i]; // The pristine anchor state component $f^\\mu_{{base}}$.
        f_temp[i] = f[i]; // The intermediate scratchpad state component $f^\\mu_{{temp}}$.
        }}

        // --- Execute 6-Stage RKF45 Tableau ---
        for (int stage = 1; stage <= 6; stage++) {{
            interpolation_kernel_{spacetime}(&commondata, f_temp, metric, connection, chunk_size, stream_idx);
            calculate_ode_rhs_kernel(f_temp, metric, connection, k_bundle, stage, chunk_size, stream_idx);

            if (stage < 6) {{
                rkf45_stage_update(f, k_bundle, h, stage, chunk_size, f_temp, stream_idx);
            }}
        }}

        // --- Finalize Step and Control Adaptive Step Size ---
        rkf45_finalize_and_control(&commondata, f, f_base, k_bundle, h, status, affine_param, rejection_retries, chunk_size, stream_idx);

        // --- Write Output & Check Break Conditions ---
        if (*rejection_retries == 0) {{
            fprintf(fp, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\\n",
                    *affine_param, f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8]);
            steps++;
        }}

        double r2 = f[1]*f[1] + f[2]*f[2] + f[3]*f[3]; // The squared Cartesian radius $r^2$ from the origin.
        if (r2 > (r_escape * r_escape)) {{
            printf("particle escaped to r > %.2f.\\n", r_escape);
            break;
        }}

        if (fabs(f[4]) > p_t_max) {{
            printf("Temporal momentum $|p_t|$ exceeded numerical limit.\\n");
            break;
        }}

        if (*status == FAILURE_RKF45_REJECTION_LIMIT) {{
            const char *status_names[] = {{
                "ACTIVE",
                "REJECTED",
                "FAILURE_RKF45_REJECTION_LIMIT"
            }};
            printf("Integration terminated organically with status: %s (%d)\\n", status_names[*status], *status);
            break;
        }}
    }}
    fclose(fp);
    printf("Integration finished after %d steps. Final lambda = %.4f\\n", steps, *affine_param);

    // --- Step 7: Post-Integration Error Diagnostics ---
    conserved_quantities_t cq_final; // Struct holding constants of motion evaluated at the trajectory boundary.
    calculate_conserved_quantities_universal_{spacetime}_{particle}(&commondata, &all_photons, num_rays, &cq_final);

    normalization_constraint_t norm_final; // Struct tracking the residual of the Hamiltonian null constraint $p^\\mu p_\\mu = 0$.
    interpolation_kernel_{spacetime}(&commondata, f, metric, NULL, chunk_size, stream_idx);
    normalization_constraint_photon(f, metric, &norm_final, chunk_size, stream_idx);

    printf("\\nFinal Normalization Constraint: |C| = %.4e\\n", fabs(norm_final.C));
    printf("Conservation Absolute Errors:\\n");
    printf("  Delta E  = %.4e\\n", fabs(cq_final.E - cq_init.E));
    printf("  Delta Lz = %.4e\\n", fabs(cq_final.Lz - cq_init.Lz));
    printf("  Delta Q  = %.4e\\n", fabs(cq_final.Q - cq_init.Q));

    // --- Memory Cleanup ---
    free(f);
    free(f_base);
    free(f_temp);
    free(metric);
    free(connection);
    free(k_bundle);
    free(affine_param);
    free(h);
    free(rejection_retries);
    free(status);

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

    project_name = "photon_geodesic_integrator"
    project_dir = os.path.join("project", project_name)

    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"
    enable_parallel_codegen = True

    if os.path.exists(project_dir):
        shutil.rmtree(project_dir)
    os.makedirs(project_dir, exist_ok=True)

    print(f"Acquiring symbolic data for {GEO_KEY}...")
    metric_data = Analytic_Spacetimes[SPACETIME]
    geodesic_data = Geodesic_Equations[GEO_KEY]

    print("Registering Split-Pipeline Physics Kernels...")
    register_struct_definitions()

    # 1. Fundamental Tensors
    g4DD_metric.g4DD_metric(metric_data.g4DD, SPACETIME, PARTICLE)
    connections.connections(geodesic_data.Gamma4UDD, SPACETIME, PARTICLE)
    conserved_quantities.conserved_quantities(SPACETIME, PARTICLE)
    normalization_constraint.normalization_constraint(
        geodesic_data.norm_constraint_expr, PARTICLE
    )

    # 2. Split-Pipeline Modular Kernels
    p0_reverse_kernel.p0_reverse_kernel(geodesic_data.p0_photon)
    interpolation_kernel.interpolation_kernel(SPACETIME)
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(
        geodesic_data.geodesic_rhs, geodesic_data.xx
    )
    rkf45_stage_update.rkf45_stage_update()
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel()

    # Execute the user-provided main_c implementation
    main_c(SPACETIME, PARTICLE)

    # Native NRPy Cleanup
    for internal_func in [f"g4DD_metric_{SPACETIME}", f"connections_{SPACETIME}"]:
        cfc.CFunction_dict.pop(internal_func, None)

    # Set Relevent Code paramters
    par.glb_code_params_dict["rkf45_absolute_error_tolerance"].defaultvalue = 1e-17
    par.glb_code_params_dict["rkf45_error_tolerance"].defaultvalue = 1e-17
    par.glb_code_params_dict["rkf45_h_max"].defaultvalue = 10.0
    par.glb_code_params_dict["rkf45_h_min"].defaultvalue = 1e-20
    par.glb_code_params_dict["rkf45_max_retries"].defaultvalue = 15
    par.glb_code_params_dict["rkf45_safety_factor"].defaultvalue = 0.9
    # ---------------------------------------------------------
    # PART 1: PARALLEL CODE GENERATION & PARAMETER SETUP
    # ---------------------------------------------------------
    if enable_parallel_codegen:
        pcg.do_parallel_codegen()

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

    # Required macros for Single-Ray CPU integration
    macro_defs = """
    #ifndef BUNDLE_CAPACITY
    #define BUNDLE_CAPACITY 1
    #endif
    """
    Bdefines_h.register_BHaH_defines("gpu_batch_macros", macro_defs)

    cpu_macros = {
        "ReadCUDA(ptr)": "#define ReadCUDA(ptr) (*(ptr))",
        "WriteCUDA(ptr, val)": "#define WriteCUDA(ptr, val) (*(ptr) = (val))",
        "BHAH_MALLOC_DEVICE(a, sz)": "#define BHAH_MALLOC_DEVICE(a, sz) BHAH_MALLOC(a, sz)",
        "BHAH_FREE_DEVICE(a)": "#define BHAH_FREE_DEVICE(a) BHAH_FREE(a)",
        "MulCUDA(a, b)": "#define MulCUDA(a, b) ((a) * (b))",
        "DivCUDA(a, b)": "#define DivCUDA(a, b) ((a) / (b))",
        "AddCUDA(a, b)": "#define AddCUDA(a, b) ((a) + (b))",
        "FusedMulAddCUDA(a, b, c)": "#define FusedMulAddCUDA(a, b, c) ((a) * (b) + (c))",
        "AbsCUDA(val)": "#define AbsCUDA(val) fabs(val)",
        "SqrtCUDA(val)": "#define SqrtCUDA(val) sqrt(val)",
        "PowCUDA(a, b)": "#define PowCUDA(a, b) pow(a, b)",
        "BHAH_HD_FUNC": "#define BHAH_HD_FUNC",
        "BHAH_HD_INLINE": "#define BHAH_HD_INLINE static inline",
        "BHAH_WARP_ATOMIC_ADD(ptr, val)": '#define BHAH_WARP_ATOMIC_ADD(ptr, val) _Pragma("omp atomic") *(ptr) += (val)\n',
        "GLOBAL_COMMONDATA_EXTERN": "// CPU passes commondata by reference, no global needed.",
        "BHAH_DEVICE_SYNC()": "#define BHAH_DEVICE_SYNC() do {} while(0)",
    }

    # C. BHaH Defines (Includes GSL headers)
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

    # D. Makefile
    addl_cflags = [
        "$(shell gsl-config --cflags)",
        "-fopenmp",
        "-O3",
        "-DDEBUG",
        "-Wno-stringop-truncation",
    ]
    addl_libs = ["$(shell gsl-config --libs)", "-lm"]

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

    # ---------------------------------------------------------
    # PART 2: FINALIZE
    # ---------------------------------------------------------

    # Define the directory containing the visualization assets relative to the repository root
    vis_dir = os.path.join("nrpy", "helpers", "geodesic_visualizations")

    # Locate the visualization script
    vis_script_src = os.path.join(vis_dir, "visualize_trajectory.py")

    # Copy the visualization script into the generated project directory
    if os.path.exists(vis_script_src):
        shutil.copy(vis_script_src, project_dir)
    else:
        print(
            f"Warning: Visualization script not found at {vis_script_src}. Please ensure it exists."
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
    print("    python3 visualize_trajectory.py --particle_type Photon\n")
