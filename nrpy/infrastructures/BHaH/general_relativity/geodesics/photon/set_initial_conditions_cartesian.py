"""
Generates the C function to set the initial state vector for a photon.

Author: Zachariah B. Etienne (Initial), Dalton J. Moone (Refactor)
"""

import nrpy.c_function as cfc


def set_initial_conditions_cartesian() -> None:
    """
    Generate and register the C orchestrator for setting photon initial conditions.

    This function generates the C code that sets the complete 9-component initial
    state vector for a single light ray. It orchestrates a sequence of geometric
    calculations and calls to other C engines to accomplish this.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "stdlib.h", "stdio.h"]
    desc = r"""@brief Sets the full initial state for a light ray in Cartesian coordinates.

    This function orchestrates the setup of the initial state vector y[9] for a
    single photon. It sets the initial position to the camera's location, computes
    the initial spatial momentum from the aiming vector, and then calls other
    engines to compute the initial time momentum by solving the null condition.

    @param[in]  commondata  Pointer to commondata struct with runtime parameters.
    @param[in]  params      Pointer to params struct (unused, for signature compatibility).
    @param[in]  metric      Pointer to the metric_params struct specifying the metric type.
    @param[in]  camera_pos  3D Cartesian coordinates of the camera.
    @param[in]  target_pos  3D Cartesian coordinates of the target pixel on the window plane.
    @param[out] y_out       The 9-component output array for the initial state vector.
    """
    name = "set_initial_conditions_cartesian"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
                const metric_params *restrict metric,
                const double camera_pos[3], const double target_pos[3],
                double y_out[9]"""

    # The body is algorithmic, not symbolic, so it is defined as a raw C string.
    body = r"""
    // --- Step 1: Set the initial position to the camera's location ---
    y_out[0] = 0.0; // t
    y_out[1] = camera_pos[0]; // x
    y_out[2] = camera_pos[1]; // y
    y_out[3] = camera_pos[2]; // z
    y_out[8] = 0.0; // L (integrated path length)

    // --- Step 2: Calculate the aiming vector V and set spatial momentum ---
    // The initial reverse-time spatial momentum p^i is parallel to the aiming vector.
    const double V_x = target_pos[0] - camera_pos[0];
    const double V_y = target_pos[1] - camera_pos[1];
    const double V_z = target_pos[2] - camera_pos[2];
    const double mag_V = sqrt(V_x*V_x + V_y*V_y + V_z*V_z);

    if (mag_V > 1e-12) {
        const double inv_mag_V = 1.0 / mag_V;
        y_out[5] = V_x * inv_mag_V; // p^x
        y_out[6] = V_y * inv_mag_V; // p^y
        y_out[7] = V_z * inv_mag_V; // p^z
    } else {
        // Should not happen, but as a fallback, aim along the x-axis.
        y_out[5] = 1.0; y_out[6] = 0.0; y_out[7] = 0.0;
    }

    // --- Step 3: Calculate the time component p^t using the null condition ---
    // Allocate the metric struct on the HEAP to prevent stack overflow in parallel loops.
    metric_struct *g4DD = (metric_struct *)malloc(sizeof(metric_struct));
    if (g4DD == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for metric_struct in set_initial_conditions.\n");
        exit(1);
    }

    // Call the metric dispatcher to get g_μν at the camera's location.
    g4DD_metric(commondata, params, metric, y_out, g4DD);

    // Call the physics engine to solve the null condition for p^t.
    y_out[4] = calculate_p0_reverse(g4DD, y_out);

    // Free the heap-allocated memory before returning.
    free(g4DD);
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )