"""
Generates the C engine to handle a window plane intersection.
"""

import os
import sys
import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par

def handle_window_plane_intersection() -> None:
    """
    Generate and register the C engine for processing window plane intersections.
    """
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "window_center_x", "window_center_y", "window_center_z",
            "camera_pos_x", "camera_pos_y", "camera_pos_z",
            "window_up_vec_x", "window_up_vec_y", "window_up_vec_z",
            "window_width", "window_height",
        ],
        [50.0, 0.0, 0.0, 51.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        commondata=True,
        add_to_parfile=True,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>", "<stdbool.h>"]
    desc = r"""@brief Handles a window plane intersection without terminating the ray."""
    name = "handle_window_plane_intersection"
    cfunc_type = "bool"
    params = """
    const PhotonStateSoA *restrict all_photons,
    const long int num_rays,
    const long int photon_idx,
    const commondata_struct *restrict commondata,
    blueprint_data_t *restrict final_blueprint_data
    """

    body = r"""
    const double intersection_pos[3] = {
        all_photons->window_event_f_intersect[IDX_GLOBAL(1, photon_idx, num_rays)], 
        all_photons->window_event_f_intersect[IDX_GLOBAL(2, photon_idx, num_rays)], 
        all_photons->window_event_f_intersect[IDX_GLOBAL(3, photon_idx, num_rays)]
    };
    
    const double window_center[3] = {commondata->window_center_x, commondata->window_center_y, commondata->window_center_z};
    
    double w_normal[3] = {
        commondata->window_center_x - commondata->camera_pos_x,
        commondata->window_center_y - commondata->camera_pos_y,
        commondata->window_center_z - commondata->camera_pos_z
    };
    
    double mag_w_norm = sqrt(w_normal[0]*w_normal[0] + w_normal[1]*w_normal[1] + w_normal[2]*w_normal[2]);
    if (mag_w_norm > 1e-12) {
        w_normal[0] /= mag_w_norm; w_normal[1] /= mag_w_norm; w_normal[2] /= mag_w_norm;
    }
    
    const double window_up_vector[3] = {commondata->window_up_vec_x, commondata->window_up_vec_y, commondata->window_up_vec_z};

    double w_z[3] = { w_normal[0], w_normal[1], w_normal[2] };
    double w_x[3];
    w_x[0] = window_up_vector[1]*w_z[2] - window_up_vector[2]*w_z[1];
    w_x[1] = window_up_vector[2]*w_z[0] - window_up_vector[0]*w_z[2];
    w_x[2] = window_up_vector[0]*w_z[1] - window_up_vector[1]*w_z[0];
    
    double mag_w_x = sqrt(w_x[0]*w_x[0] + w_x[1]*w_x[1] + w_x[2]*w_x[2]);

    if (mag_w_x < 1e-9) {
        double alternative_up[3] = {1.0, 0.0, 0.0};
        if (fabs(w_z[0]) > 0.999) { alternative_up[0] = 0.0; alternative_up[1] = 1.0; }
        w_x[0] = alternative_up[1]*w_z[2] - alternative_up[2]*w_z[1];
        w_x[1] = alternative_up[2]*w_z[0] - alternative_up[0]*w_z[2];
        w_x[2] = alternative_up[0]*w_z[1] - alternative_up[1]*w_z[0];
        mag_w_x = sqrt(w_x[0]*w_x[0] + w_x[1]*w_x[1] + w_x[2]*w_x[2]);
    }

    double inv_mag_w_x = 1.0 / mag_w_x;
    w_x[0] *= inv_mag_w_x; w_x[1] *= inv_mag_w_x; w_x[2] *= inv_mag_w_x;

    double w_y[3];
    w_y[0] = w_z[1]*w_x[2] - w_z[2]*w_x[1];
    w_y[1] = w_z[2]*w_x[0] - w_z[0]*w_x[2];
    w_y[2] = w_z[0]*w_x[1] - w_z[1]*w_x[0];

    const double vec_w[3] = {
        intersection_pos[0] - window_center[0],
        intersection_pos[1] - window_center[1],
        intersection_pos[2] - window_center[2]
    };

    const double y_w = vec_w[0]*w_x[0] + vec_w[1]*w_x[1] + vec_w[2]*w_x[2];
    const double z_w = vec_w[0]*w_y[0] + vec_w[1]*w_y[1] + vec_w[2]*w_y[2];

    // Validate Bounds and Fill Data
    if (fabs(y_w) <= commondata->window_width / 2.0 && fabs(z_w) <= commondata->window_height / 2.0) {
        final_blueprint_data->y_w = y_w;
        final_blueprint_data->z_w = z_w;
        final_blueprint_data->t_w = all_photons->window_event_f_intersect[IDX_GLOBAL(0, photon_idx, num_rays)];
        final_blueprint_data->L_w = all_photons->window_event_f_intersect[IDX_GLOBAL(8, photon_idx, num_rays)]; 
        return true;j
    }

    return false;
    """

    cfc.register_CFunction(includes=includes, desc=desc, name=name, cfunc_type=cfunc_type, params=params, body=body)