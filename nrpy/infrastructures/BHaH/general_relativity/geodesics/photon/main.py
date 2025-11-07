"""
Generates the main() C function for the photon geodesic integrator project.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def main() -> None:
    """
    Generate and register the main() C function.

    This function generates the C code for the main() function, which serves as
    the entry point for the entire executable. It acts as the master orchestrator,
    managing the simulation's lifecycle from parameter parsing and data loading
    to executing the core integration loop and saving the final results.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "stdio.h", "stdlib.h"]
    desc = r"""@brief Main entry point for the geodesic integrator.

    Acts as the master dispatcher, selecting the appropriate integration
    pipeline (analytic vs. numerical) based on the 'use_numerical_pipeline'
    parameter. It manages the lifecycle of all major data structures.
    """
    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"

    # The body is algorithmic, not symbolic, so it is defined as a raw C string.
    body = r"""
    // --- Step 1: Initialize Core Data Structures ---
    commondata_struct commondata;
    params_struct params = {0}; // Initialize to zero; unused in this project but required by signatures.
    metric_params metric;

    // --- Step 2: Set Default Parameters and Parse User Input ---
    commondata_struct_set_to_default(&commondata);
    cmdline_input_and_parfile_parser(&commondata, argc, argv);

    // --- Step 3: Set Metric Type Enum Based on User Choice ---
    // This is used by the analytic workers, even when called from the numerical pipeline's placeholder.
    if (commondata.metric_choice == 0) {
        metric.type = (commondata.a_spin == 0.0) ? Schwarzschild : Kerr;
    } else if (commondata.metric_choice == 1) {
        metric.type = Schwarzschild_Standard;
    } else {
        fprintf(stderr, "Error: Invalid metric_choice = %d for this validation build.\n", commondata.metric_choice);
        fprintf(stderr, "       Please use 0 (Kerr-Schild) or 1 (Standard Schwarzschild).\n");
        exit(1);
    }

    // --- Step 4: Load k-d tree snapshot files (used by both pipelines) ---
    CustomKDTree *kdtree_snapshots = NULL;
    double *snapshot_times = NULL;
    int num_snapshots = 0;
    num_snapshots = load_all_kdtree_snapshots(&commondata, &kdtree_snapshots, &snapshot_times);

    // --- Step 5: Print Simulation Banner ---
    printf("=============================================\n");
    printf("  Photon Geodesic Integrator (Batch Mode)  \n");
    printf("=============================================\n");
    if (commondata.use_numerical_pipeline) {
        printf("PIPELINE: NUMERICAL (Validation Mode)\n");
    } else {
        printf("PIPELINE: ANALYTIC\n");
    }
    printf("Metric: %s (a=%.3f)\n", (metric.type == Kerr) ? "Kerr-Schild" : "Schwarzschild-Standard", commondata.a_spin);
    printf("Scan Resolution: %d x %d\n", commondata.scan_density, commondata.scan_density);
    if (num_snapshots > 0) {
        printf("Accretion Disk: ENABLED (%d snapshots loaded)\n", num_snapshots);
    } else {
        printf("Accretion Disk: DISABLED (no snapshots found)\n");
    }

    // --- Step 6: Main Logic Dispatcher ---
    long int num_rays = (long int)commondata.scan_density * commondata.scan_density;
    blueprint_data_t *results_buffer = (blueprint_data_t *)malloc(sizeof(blueprint_data_t) * num_rays);
    if (results_buffer == NULL) { exit(1); }

    if (commondata.use_numerical_pipeline) {
        // Call the NUMERICAL pipeline orchestrator.
        batch_integrator_numerical(&commondata, &params, &metric, num_rays,
                                    num_snapshots, kdtree_snapshots, snapshot_times,
                                    results_buffer);
    } else {
        // The analytic pipeline is currently disabled in favor of the more feature-rich numerical pipeline.
        printf("ANALYTIC INTEGRATOR (CURRENTLY DISABLED/DEPRECATED) CALLED. EXITING.\n");
        exit(1);
    }

    printf("Scan finished. Writing %ld ray results to light_blueprint.bin...\n", num_rays);
    FILE *fp_blueprint = fopen("light_blueprint.bin", "wb");
    if (fp_blueprint == NULL) { perror("Error opening blueprint file"); exit(1); }
    fwrite(results_buffer, sizeof(blueprint_data_t), num_rays, fp_blueprint);
    fclose(fp_blueprint);
    free(results_buffer);

    // --- Step 7: Cleanup ---
    printf("Unloading k-d tree snapshots...\n");
    if (kdtree_snapshots != NULL) {
        for (int i = 0; i < num_snapshots; ++i) {
            unload_kdtree_snapshot(&kdtree_snapshots[i]);
        }
        free(kdtree_snapshots);
    }
    if (snapshot_times != NULL) { free(snapshot_times); }

    printf("\nRun complete.\n");
    return 0;
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )