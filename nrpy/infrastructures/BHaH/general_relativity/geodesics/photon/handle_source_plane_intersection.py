"""
Generates the C engine to handle a fallback source plane intersection.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def handle_source_plane_intersection() -> None:
    """
    Generate and register the C engine for processing source plane intersections.

    This function generates the C code that is called when a photon hits the
    fallback source plane. It transforms the 3D intersection point into local
    2D texture coordinates, checks if the hit is within the active "glowing"
    region, and populates the final output struct if it is.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdbool.h>",
    ]
    desc = r"""@brief Handles a source plane intersection by checking bounds and populating the blueprint.
    @param[in]  source_plane_event  Pointer to the event data struct for the intersection.
    @param[in]  commondata          Pointer to commondata struct with runtime parameters.
    @param[out] final_blueprint_data Pointer to the final output struct to be populated.
    @return True if the intersection is valid and processed, false otherwise.
    """
    name = "handle_source_plane_intersection"
    cfunc_type = "bool"
    params = """
    const event_data_struct *restrict source_plane_event,
    const commondata_struct *restrict commondata,
    blueprint_data_t *restrict final_blueprint_data
    """

    # The body is algorithmic, not symbolic, so it is defined as a raw C string.
    body = r"""
    // --- Step 1: Calculate the local (y_s, z_s) coordinates on the plane ---
    const double intersection_pos[3] = {source_plane_event->y_event[1], source_plane_event->y_event[2], source_plane_event->y_event[3]};
    const double source_plane_center[3] = {commondata->source_plane_center_x, commondata->source_plane_center_y, commondata->source_plane_center_z};
    const double source_plane_normal[3] = {commondata->source_plane_normal_x, commondata->source_plane_normal_y, commondata->source_plane_normal_z};
    const double source_up_vector[3] = {commondata->source_up_vec_x, commondata->source_up_vec_y, commondata->source_up_vec_z};

    // --- Step 2: Construct an orthonormal basis for the source plane ---
    double s_z[3] = {source_plane_normal[0], source_plane_normal[1], source_plane_normal[2]};
    // s_x = up x s_z (cross product)
    double s_x[3] = {source_up_vector[1]*s_z[2] - source_up_vector[2]*s_z[1],
                     source_up_vector[2]*s_z[0] - source_up_vector[0]*s_z[2],
                     source_up_vector[0]*s_z[1] - source_up_vector[1]*s_z[0]};
    double mag_s_x = sqrt(SQR(s_x[0]) + SQR(s_x[1]) + SQR(s_x[2]));

    // Handle case where up vector is parallel to the normal
    if (mag_s_x < 1e-9) {
        double alternative_up[3] = {1.0, 0.0, 0.0};
        if (fabs(s_z[0]) > 0.999) {
            alternative_up[0] = 0.0; alternative_up[1] = 1.0;
        }
        s_x[0] = alternative_up[1]*s_z[2] - alternative_up[2]*s_z[1];
        s_x[1] = alternative_up[2]*s_z[0] - alternative_up[0]*s_z[2];
        s_x[2] = alternative_up[0]*s_z[1] - alternative_up[1]*s_z[0];
        mag_s_x = sqrt(SQR(s_x[0]) + SQR(s_x[1]) + SQR(s_x[2]));
    }
    for(int i=0; i<3; i++) s_x[i] /= mag_s_x;

    // s_y = s_z x s_x (cross product)
    double s_y[3] = {s_z[1]*s_x[2] - s_z[2]*s_x[1],
                     s_z[2]*s_x[0] - s_z[0]*s_x[2],
                     s_z[0]*s_x[1] - s_z[1]*s_x[0]};

    // --- Step 3: Project the intersection point onto the plane's basis ---
    const double vec_s[3] = {intersection_pos[0] - source_plane_center[0],
                             intersection_pos[1] - source_plane_center[1],
                             intersection_pos[2] - source_plane_center[2]};

    const double y_s = vec_s[0]*s_x[0] + vec_s[1]*s_x[1] + vec_s[2]*s_x[2];
    const double z_s = vec_s[0]*s_y[0] + vec_s[1]*s_y[1] + vec_s[2]*s_y[2];

    // --- Step 4: Check if the hit is within the active annulus ---
    const double r_s_sq = SQR(y_s) + SQR(z_s);
    if (r_s_sq >= SQR(commondata->source_r_min) && r_s_sq <= SQR(commondata->source_r_max)) {
        // This is a valid hit. Populate the blueprint and return true.
        final_blueprint_data->termination_type = TERMINATION_TYPE_SOURCE_PLANE;
        final_blueprint_data->y_s = y_s;
        final_blueprint_data->z_s = z_s;
        final_blueprint_data->t_s = source_plane_event->t_event;
        final_blueprint_data->L_s = source_plane_event->y_event[8];
        return true;
    }

    // The intersection was outside the valid radial bounds. Return false.
    return false;
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        cfunc_type=cfunc_type,
        params=params,
        body=body,
    )