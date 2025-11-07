"""
Generates the C orchestrator for geometric event detection.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def event_detection_manager() -> None:
    """
    Generate and register the C event detection manager.

    This function generates the C code that is called at each step of the
    integration. It checks if a photon has crossed either the camera's window
    plane or the fallback source plane. If a crossing is detected, it calls the
    high-accuracy interpolation engine to find the precise intersection state.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>"]
    desc = r"""@brief Detects crossings of the window and source planes.

    This orchestrator is called at each integration step. It determines if a
    sign change has occurred in the photon's distance to either the window or
    source plane. If so, it calls the `find_event_time_and_state()` engine to
    accurately interpolate the intersection point.

    @param[in]      y_prev, y_curr, y_next  State vectors at three consecutive steps.
    @param[in]      lambda_prev, lambda_curr, lambda_next Affine parameters for the three steps.
    @param[in]      commondata              Pointer to commondata struct with runtime parameters.
    @param[in,out]  on_positive_side_of_window_prev Pointer to the state of the photon relative to the window plane at the previous step.
    @param[in,out]  on_positive_side_of_source_prev Pointer to the state of the photon relative to the source plane at the previous step.
    @param[out]     window_event            Pointer to the event_data_struct for window crossings.
    @param[out]     source_plane_event      Pointer to the event_data_struct for source plane crossings.
    """
    name = "event_detection_manager"
    params = """
        const double y_prev[9], const double y_curr[9], const double y_next[9],
        double lambda_prev, double lambda_curr, double lambda_next,
        const commondata_struct *restrict commondata,
        bool *on_positive_side_of_window_prev,
        bool *on_positive_side_of_source_prev,
        event_data_struct *restrict window_event,
        event_data_struct *restrict source_plane_event
        """

    # The body is algorithmic, not symbolic, so it is defined as a raw C string.
    body = r"""
    // --- Window Plane Detection ---
    // This logic is only performed if the caller has not already found the event.
    if (!window_event->found) {
        // Define the window plane: n_i x^i = d
        double window_plane_normal[3] = {commondata->window_center_x - commondata->camera_pos_x,
                                         commondata->window_center_y - commondata->camera_pos_y,
                                         commondata->window_center_z - commondata->camera_pos_z};
        const double mag_w_norm = sqrt(SQR(window_plane_normal[0]) + SQR(window_plane_normal[1]) + SQR(window_plane_normal[2]));
        if (mag_w_norm > 1e-12) {
            const double inv_mag_w_norm = 1.0 / mag_w_norm;
            for(int i=0;i<3;i++) window_plane_normal[i] *= inv_mag_w_norm;
        }
        const double window_plane_dist = commondata->window_center_x * window_plane_normal[0] +
                                         commondata->window_center_y * window_plane_normal[1] +
                                         commondata->window_center_z * window_plane_normal[2];

        plane_event_params window_params = {{window_plane_normal[0], window_plane_normal[1], window_plane_normal[2]}, window_plane_dist};
        // Check if the photon crossed the plane in the last full step (y_prev -> y_next)
        bool on_positive_side_curr = (plane_event_func(y_next, &window_params) > 0);
        if (on_positive_side_curr != *on_positive_side_of_window_prev) {
            // A crossing occurred. Call the high-accuracy interpolator.
            find_event_time_and_state(y_prev, y_curr, y_next, lambda_prev, lambda_curr, lambda_next,
                                      plane_event_func, &window_params, window_event);
        }
        *on_positive_side_of_window_prev = on_positive_side_curr;
    }

    // --- Source Plane Detection ---
    // This logic is only performed if the caller has not already found the event.
    if (!source_plane_event->found) {
        const double source_plane_normal[3] = {commondata->source_plane_normal_x,
                                               commondata->source_plane_normal_y,
                                               commondata->source_plane_normal_z};
        const double source_plane_dist = commondata->source_plane_center_x * source_plane_normal[0] +
                                         commondata->source_plane_center_y * source_plane_normal[1] +
                                         commondata->source_plane_center_z * source_plane_normal[2];

        plane_event_params source_params = {{source_plane_normal[0], source_plane_normal[1], source_plane_normal[2]}, source_plane_dist};
        bool on_positive_side_curr = (plane_event_func(y_next, &source_params) > 0);
        if (on_positive_side_curr != *on_positive_side_of_source_prev) {
            find_event_time_and_state(y_prev, y_curr, y_next, lambda_prev, lambda_curr, lambda_next,
                                      plane_event_func, &source_params, source_plane_event);
        }
        *on_positive_side_of_source_prev = on_positive_side_curr;
    }
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )