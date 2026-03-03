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
    # Python: Extract physical inline helpers from the global CFunction_dict registry.
    # We retrieve the 'full_function' string (prototype + body) to inject directly
    # into this translation unit, preventing linker errors for inline GPU device code.
    find_event_c_code = cfc.CFunction_dict["find_event_time_and_state"].full_function
    window_c_code = cfc.CFunction_dict["handle_window_plane_intersection"].full_function
    source_c_code = cfc.CFunction_dict["handle_source_plane_intersection"].full_function

    # Python: Concatenate helper functions into the prefunc block.
    # This ensures the compiler sees the definitions before the manager calls them.
    prefunc = find_event_c_code + "\n" + window_c_code + "\n" + source_c_code

    # Hardware Note: BHaH_device_defines.h required for d_commondata visibility.
    includes = ["BHaH_defines.h", "BHaH_device_defines.h", "BHaH_function_prototypes.h", "<math.h>", "<stdbool.h>"]
    
    desc = """@brief GPU-optimized detection of plane crossings using consolidated blueprints.

    @param f_local The thread-local 9-component array for the current state $f^\\mu_{n}$.
    @param f_p_local The thread-local state array for step $f^\\mu_{n-1}$.
    @param f_p_p_local The thread-local state array for step $f^\\mu_{n-2}$.
    @param status Pointer to the photon's termination status in the bundle buffer.
    @param blueprint Pointer to the physical blueprint data structure for result persistence.
    @param on_pos_window_prev Pointer to the persistent flag tracking the window plane side.
    @param on_pos_source_prev Pointer to the persistent flag tracking the source plane side."""
    
    cfunc_type = "BHAH_HD_FUNC void"
    name = "event_detection_manager"
    
    # Optimization: commondata pointer removed. The function accesses the global
    # constant symbol 'd_commondata' directly to reduce register pressure.
    params = (
        "const double *restrict f_local, "
        "const double *restrict f_p_local, "
        "const double *restrict f_p_p_local, "
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
    // Checks if the photon has crossed the mathematical plane defining the camera.
    if (blueprint->L_w < 0.0) {
        double w_normal[3];
        w_normal[0] = d_commondata.window_center_x - d_commondata.camera_pos_x;
        w_normal[1] = d_commondata.window_center_y - d_commondata.camera_pos_y;
        w_normal[2] = d_commondata.window_center_z - d_commondata.camera_pos_z;
        
        const double mag_inv = 1.0 / sqrt(w_normal[0]*w_normal[0] + w_normal[1]*w_normal[1] + w_normal[2]*w_normal[2]);
        for(int i=0; i<3; i++) w_normal[i] *= mag_inv;
        
        const double w_dist = d_commondata.window_center_x*w_normal[0] + d_commondata.window_center_y*w_normal[1] + d_commondata.window_center_z*w_normal[2];
        const double w_val = x*w_normal[0] + y*w_normal[1] + z*w_normal[2] - w_dist;
        
        // Logical check: Is the photon on the 'positive' side of the plane?
        const bool on_pos_curr = (w_val > 1e-10);

        if (on_pos_curr != *on_pos_window_prev) {
             double f_int[9], lam;
             // High-precision root finding to locate the exact crossing time.
             find_event_time_and_state(f_local, f_p_local, f_p_p_local, w_normal, w_dist, &lam, f_int);
             // Update the blueprint with the intersection data.
             handle_window_plane_intersection(f_int, &d_commondata, blueprint);
        }
        *on_pos_window_prev = on_pos_curr;
    }

    // --- Physical Source Plane Logic ---
    // Checks if the photon has collided with the accretion disk geometry.
    if (*status == ACTIVE) {
        const double s_normal[3] = {d_commondata.source_plane_normal_x, d_commondata.source_plane_normal_y, d_commondata.source_plane_normal_z};
        const double s_dist = d_commondata.source_plane_center_x*s_normal[0] + d_commondata.source_plane_center_y*s_normal[1] + d_commondata.source_plane_center_z*s_normal[2];
        const double s_val = x*s_normal[0] + y*s_normal[1] + z*s_normal[2] - s_dist;
        
        // Logical check: Is the photon on the 'positive' side of the plane?
        const bool on_pos_curr = (s_val > 1e-10);

        if (on_pos_curr != *on_pos_source_prev) {
             double f_int[9], lam;
             find_event_time_and_state(f_local, f_p_local, f_p_p_local, s_normal, s_dist, &lam, f_int);
             // If the intersection is within the valid disk radius, terminate the ray.
             if (handle_source_plane_intersection(f_int, &d_commondata, blueprint)) {
                 *status = TERMINATION_TYPE_SOURCE_PLANE;
             }
        }
        *on_pos_source_prev = on_pos_curr;
    }
    """
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=include_CodeParameters_h,
    )