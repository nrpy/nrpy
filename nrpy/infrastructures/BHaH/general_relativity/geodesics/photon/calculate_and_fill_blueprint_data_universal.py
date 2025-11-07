"""
Generates the C "finalizer" engine to process a completed ray trace.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def calculate_and_fill_blueprint_data_universal() -> None:
    """
    Generate and register the C finalization engine.

    This function generates the C code that is called once for each photon after
    its integration is complete. It acts as a dispatcher based on the photon's
    final termination status, calling the appropriate helper engines to compute
    the final physical and geometric quantities that are saved to the output file.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""@brief Processes a photon's final state to compute all blueprint quantities.

    This finalizer engine is called once for each completed ray. It dispatches
    to the appropriate helper engine based on the photon's termination status
    to calculate the final physical results (e.g., intensity, texture coordinates)
    and populates the `blueprint_data_t` struct for output.

    @param[in]  commondata     Pointer to commondata struct with runtime parameters.
    @param[in]  params         Pointer to params struct (unused, for signature compatibility).
    @param[in]  metric         Pointer to the metric_params struct specifying the metric type.
    @param[in]  photon         Pointer to the final state of the completed photon.
    @param[in]  window_center  3D Cartesian coordinates of the window plane's center.
    @param[in]  n_x, n_y       Orthonormal basis vectors for the window plane.
    @return A fully populated `blueprint_data_t` struct ready for output.
    """
    name = "calculate_and_fill_blueprint_data_universal"
    cfunc_type = "blueprint_data_t"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
                const metric_params *restrict metric,
                const PhotonState *restrict photon,
                const double window_center[3], const double n_x[3], const double n_y[3]"""

    # The body is algorithmic, not symbolic, so it is defined as a raw C string.
    body = r"""
    // Initialize all fields to zero.
    blueprint_data_t result = {0};
    result.termination_type = photon->status;

    // --- Step 1: Always populate window data if a crossing was found ---
    if (photon->window_event_data.found) {
        const double *y_event = photon->window_event_data.y_event;
        const double pos_w_cart[3] = {y_event[1], y_event[2], y_event[3]};
        const double vec_w[3] = {pos_w_cart[0] - window_center[0], pos_w_cart[1] - window_center[1], pos_w_cart[2] - window_center[2]};
        // Project intersection point onto window's orthonormal basis to get texture coordinates.
        result.y_w = vec_w[0]*n_x[0] + vec_w[1]*n_x[1] + vec_w[2]*n_x[2];
        result.z_w = vec_w[0]*n_y[0] + vec_w[1]*n_y[1] + vec_w[2]*n_y[2];
        result.L_w = y_event[8];
        result.t_w = photon->window_event_data.t_event;
    }

    // --- Step 2: Populate remaining fields based on the specific termination type ---
    if (photon->status == TERMINATION_TYPE_SOURCE_PLANE) {
        // Call the source plane helper to validate the hit and calculate geometric properties.
        handle_source_plane_intersection(&photon->source_event_data, commondata, &result);
    } else if (photon->status == TERMINATION_TYPE_DISK) {
        // Call the dedicated physics engine for disk hits to perform radiative transfer.
        handle_disk_intersection(
            photon->y,                      // Photon's final state vector
            &photon->nearest_neighbor,      // The stored particle data from the k-d tree hit
            commondata, params, metric,
            &result                         // The blueprint to be filled
        );
    } else if (photon->status == TERMINATION_TYPE_CELESTIAL_SPHERE) {
        // Convert the final Cartesian position to spherical polar angles.
        const double *final_y = photon->y;
        const double x = final_y[1];
        const double y = final_y[2];
        const double z = final_y[3];
        const double r = sqrt(SQR(x) + SQR(y) + SQR(z));
        if (r > 1e-9) {
            result.final_theta = acos(z / r);
            result.final_phi = atan2(y, x);
        }
    }
    // For FAILURE types, no other fields need to be set.

    return result;
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )