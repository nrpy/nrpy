"""
Generates the C orchestrator for geometric event detection.

This module provides the high-level logic for detecting crossings of the
observer window and the source emission plane. It evaluates the scalar product
between the photon trajectory vector and the target plane normal, triggering
high-precision root-finding upon a sign transition. The logic executes entirely
on thread-local arrays to map optimally to the streaming bundle architecture.

Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc


def event_detection_manager() -> None:
    """
    Register the event_detection_manager C function with the Batch 4 API.

    :raises SystemError: If C function registration fails within the NRPy+ pipeline.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>", "<stdbool.h>"]
    desc = """@brief GPU-optimized detection of plane crossings using consolidated blueprints.

    @param f_local The thread-local 9-component array for the current state $f^\\mu_{n}$.
    @param f_p_local The thread-local state array for step $f^\\mu_{n-1}$.
    @param f_p_p_local The thread-local state array for step $f^\\mu_{n-2}$.
    @param commondata Pointer to the globally constant parameters struct.
    @param status Pointer to the photon's termination status in the bundle buffer.
    @param blueprint Pointer to the physical blueprint data structure for result persistence.
    @param on_pos_window_prev Pointer to the persistent flag tracking the window plane side.
    @param on_pos_source_prev Pointer to the persistent flag tracking the source plane side."""
    
    cfunc_type = "BHAH_HD_FUNC void"
    name = "event_detection_manager"
    params = (
        "const double *restrict f_local, "
        "const double *restrict f_p_local, "
        "const double *restrict f_p_p_local, "
        "const commondata_struct *restrict commondata, "
        "termination_type_t *restrict status, "
        "blueprint_data_t *restrict blueprint, "
        "bool *restrict on_pos_window_prev, "
        "bool *restrict on_pos_source_prev"
    )
    include_CodeParameters_h = False
    
    body = r"""
    const double x = f_local[1]; 
    const double y = f_local[2]; 
    const double z = f_local[3];

    // --- Observer Window Plane Logic ---
    if (blueprint->L_w < 0.0) {
        double w_normal[3];
        w_normal[0] = commondata->window_center_x - commondata->camera_pos_x;
        w_normal[1] = commondata->window_center_y - commondata->camera_pos_y;
        w_normal[2] = commondata->window_center_z - commondata->camera_pos_z;
        const double mag_inv = 1.0 / sqrt(w_normal[0]*w_normal[0] + w_normal[1]*w_normal[1] + w_normal[2]*w_normal[2]);
        for(int i=0; i<3; i++) w_normal[i] *= mag_inv;
        const double w_dist = commondata->window_center_x*w_normal[0] + commondata->window_center_y*w_normal[1] + commondata->window_center_z*w_normal[2];
        const double w_val = x*w_normal[0] + y*w_normal[1] + z*w_normal[2] - w_dist;
        const bool on_pos_curr = (w_val > 1e-10);

        if (on_pos_curr != *on_pos_window_prev) {
             double f_int[9], lam;
             find_event_time_and_state(f_local, f_p_local, f_p_p_local, w_normal, w_dist, &lam, f_int);
             handle_window_plane_intersection(f_int, commondata, blueprint);
        }
        *on_pos_window_prev = on_pos_curr;
    }

    // --- Physical Source Plane Logic ---
    if (*status == ACTIVE) {
        const double s_normal[3] = {commondata->source_plane_normal_x, commondata->source_plane_normal_y, commondata->source_plane_normal_z};
        const double s_dist = commondata->source_plane_center_x*s_normal[0] + commondata->source_plane_center_y*s_normal[1] + commondata->source_plane_center_z*s_normal[2];
        const double s_val = x*s_normal[0] + y*s_normal[1] + z*s_normal[2] - s_dist;
        const bool on_pos_curr = (s_val > 1e-10);

        if (on_pos_curr != *on_pos_source_prev) {
             double f_int[9], lam;
             find_event_time_and_state(f_local, f_p_local, f_p_p_local, s_normal, s_dist, &lam, f_int);
             if (handle_source_plane_intersection(f_int, commondata, blueprint)) {
                 *status = TERMINATION_TYPE_SOURCE_PLANE;
             }
        }
        *on_pos_source_prev = on_pos_curr;
    }
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=False,
    )