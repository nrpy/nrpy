"""
Generates the main() C function for the Project Singularity-Axiom pipeline.

This module acts as the master orchestrator, handling the initialization of
global simulation parameters, memory allocation for trajectory results, and
the final serialization of the light blueprint to disk.

Author: Dalton J. Moone
"""

import logging
import sys

import nrpy.c_function as cfc

# Configure logging to adhere to static analysis and project standards
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main(spacetime_name: str) -> None:
    """
    Register the master orchestrator C function for the geodesic integrator.

    This function generates the main entry point (main.c). It coordinates
    the transition from parameter parsing to numerical integration and
    final data preservation.

    :param spacetime_name: The identifier for the spacetime metric (e.g., 'Kerr').
    """
    # 1. Define C-Function metadata in order of appearance
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "stdio.h", "stdlib.h"]

    desc = f"""@brief Master entry point for the Project Singularity-Axiom integrator.

    Algorithm:
    1. Initializes core data structures and sets default physical constants.
    2. Parses command-line arguments and parameter files to override defaults.
    3. Allocates a contiguous results buffer for the light blueprint.
    4. Dispatches the batched numerical integrator for the {spacetime_name} metric.
    5. Serializes the final photon termination data to 'light_blueprint.bin'."""

    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"

    # 2. Build the C body with internal descriptive comments and the Preamble pattern
    body = f"""
    // --- CORE STRUCTURE INITIALIZATION ---
    // Populates the master structure with the base metric parameters.
    commondata_struct commondata; // Master parameter structure containing physical constants and simulation flags.

    // --- PARAMETER PARSING ---
    // Sets default metric properties and overrides them based on host system input.
    commondata_struct_set_to_default(&commondata); // Populates the data structure with default numerical values.
    cmdline_input_and_parfile_parser(&commondata, argc, argv); // Overrides default structure values using standard input arguments.

    // --- TELEMETRY AND PARAMETER VERIFICATION ---
    // Outputs core physical and architectural variables to the standard output.
    printf("=============================================\\n");
    printf("  Project Singularity-Axiom: Geodesic Engine  \\n");
    printf("=============================================\\n");
    printf("Spacetime Metric: {spacetime_name}\\n");
    
    printf("--- Analytic Spacetimes ---\\n");
    printf("Mass Scale (M): %.3f\\n", commondata.M_scale);
    printf("Spin Parameter (a): %.3f\\n", commondata.a_spin);

    printf("--- Initial Conditions ---\\n");
    printf("Scan Density: %d\\n", commondata.scan_density);
    printf("Start Time (t_start): %.3f\\n", commondata.t_start);

    printf("--- Batch Integrator Numerical ---\\n");
    printf("Initial Step Size (h_initial): %.3f\\n", commondata.numerical_initial_h);
    printf("Max Momentum (p_t,max): %.3f\\n", commondata.p_t_max);
    printf("Conservation Check: %d\\n", commondata.perform_conservation_check);
    printf("Escape Radius (r_escape): %.3f\\n", commondata.r_escape);
    printf("Slot Delta t: %.3f\\n", commondata.slot_manager_delta_t);
    printf("Slot Min t: %.3f\\n", commondata.slot_manager_t_min);
    printf("Integration Max t: %.3f\\n", commondata.t_integration_max);

    printf("--- RKF45 Control ---\\n");
    printf("Absolute Tolerance: %e\\n", commondata.rkf45_absolute_error_tolerance);
    printf("Relative Tolerance: %e\\n", commondata.rkf45_error_tolerance);
    printf("Max Step (h_max): %.3f\\n", commondata.rkf45_h_max);
    printf("Min Step (h_min): %e\\n", commondata.rkf45_h_min);
    printf("Max Retries: %d\\n", commondata.rkf45_max_retries);
    printf("Safety Factor: %.3f\\n", commondata.rkf45_safety_factor);

    printf("--- Source Plane ---\\n");
    printf("Center (x, y, z): %.3f, %.3f, %.3f\\n", commondata.source_plane_center_x, commondata.source_plane_center_y, commondata.source_plane_center_z);
    printf("Normal (x, y, z): %.3f, %.3f, %.3f\\n", commondata.source_plane_normal_x, commondata.source_plane_normal_y, commondata.source_plane_normal_z);
    printf("Max Radius (r_max): %.3f\\n", commondata.source_r_max);
    printf("Min Radius (r_min): %.3f\\n", commondata.source_r_min);
    printf("Up Vector (x, y, z): %.3f, %.3f, %.3f\\n", commondata.source_up_vec_x, commondata.source_up_vec_y, commondata.source_up_vec_z);

    printf("--- Window Plane ---\\n");
    printf("Camera Pos (x, y, z): %.3f, %.3f, %.3f\\n", commondata.camera_pos_x, commondata.camera_pos_y, commondata.camera_pos_z);
    printf("Window Center (x, y, z): %.3f, %.3f, %.3f\\n", commondata.window_center_x, commondata.window_center_y, commondata.window_center_z);
    printf("Window Height: %.3f\\n", commondata.window_height);
    printf("Window Width: %.3f\\n", commondata.window_width);
    printf("Window Up Vec (x, y, z): %.3f, %.3f, %.3f\\n", commondata.window_up_vec_x, commondata.window_up_vec_y, commondata.window_up_vec_z);

    // --- Step 4: Resource Allocation ---
    const long int num_rays = (long int)commondata.scan_density * commondata.scan_density; // Total global photon count.

    // Performance Optimization: Use 'restrict' to qualify the results buffer to aid compiler optimization.
    blueprint_data_t *restrict results_buffer = (blueprint_data_t *)malloc(sizeof(blueprint_data_t) * num_rays); // Contiguous memory for the light blueprint.

    if (results_buffer == NULL) {{
        fprintf(stderr, "FATAL: Failed to allocate %ld bytes for results_buffer. Check system RAM.\\n", (long int)(sizeof(blueprint_data_t) * num_rays));
        exit(1);
    }}

    // --- Step 5: Execute Numerical Integration Pipeline ---
    // Note: batch_integrator_numerical handles the SoA initialization and GPU offloading internal to its call.
    batch_integrator_numerical(&commondata, num_rays, results_buffer);

    // --- Step 6: Data Serialization ---
    printf("Integration complete. Serializing %ld rays to light_blueprint.bin...\\n", num_rays);
    FILE *fp_blueprint = fopen("light_blueprint.bin", "wb");
    if (fp_blueprint == NULL) {{
        perror("FATAL: Error opening light_blueprint.bin for writing");
        free(results_buffer);
        exit(1);
    }}
    fwrite(results_buffer, sizeof(blueprint_data_t), num_rays, fp_blueprint);
    fclose(fp_blueprint);

    // --- Step 7: Final Cleanup & Shutdown ---
    free(results_buffer);
    printf("Simulation successful. Global exit code: 0.\\n");
    return 0;
    """

    # 3. Register the function with the NRPy environment
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

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
