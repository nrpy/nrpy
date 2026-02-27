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
    // --- Step 1: Initialize Core Data Structures ---
    commondata_struct commondata; // Master parameter structure containing physical constants and simulation flags.

    // --- Step 2: Set Default Parameters and Parse User Input ---
    commondata_struct_set_to_default(&commondata);
    cmdline_input_and_parfile_parser(&commondata, argc, argv);

    // --- Step 3: Simulation Preamble & Telemetry ---
    printf("=============================================\\n");
    printf("  Project Singularity-Axiom: Geodesic Engine  \\n");
    printf("=============================================\\n");
    printf("Spacetime Metric: {spacetime_name}\\n");
    printf("Spin Parameter (a): %.3f\\n", commondata.a_spin);
    printf("Grid Resolution: %d x %d\\n", commondata.scan_density, commondata.scan_density);

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
    # Standard testing routine for the generator
    SPACETIME = "KerrSchild_Cartesian"
    logger.info("Generating main.c orchestrator for spacetime: %s", SPACETIME)

    try:
        main(SPACETIME)
        if "main" in cfc.CFunction_dict:
            with open("main.c", "w", encoding="utf-8") as f:
                f.write(cfc.CFunction_dict["main"].full_function)
            logger.info("Successfully generated main.c")
        else:
            logger.error("Function registration failed.")
            sys.exit(1)
    except (RuntimeError, OSError) as e:
        logger.error("Generation failed: %s", e)
        sys.exit(1)
