"""
Generates the C engine to handle a fallback source plane intersection.
"""
import os
import sys

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par

def handle_source_plane_intersection() -> None:
    par.register_CodeParameters(
        "REAL", __name__,
        ["source_plane_normal_x", "source_plane_normal_y", "source_plane_normal_z",
         "source_plane_center_x", "source_plane_center_y", "source_plane_center_z",
         "source_up_vec_x", "source_up_vec_y", "source_up_vec_z",
         "source_r_min", "source_r_max"],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 6.0, 25.0],
        commondata=True, add_to_parfile=True,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>", "<stdbool.h>"]
    desc = r"""@brief Handles a source plane intersection using SoA architecture."""
    name = "handle_source_plane_intersection"
    cfunc_type = "bool"
    

    params = """
    const PhotonStateSoA *restrict all_photons,
    long int num_rays,
    long int photon_idx,
    const commondata_struct *restrict commondata,
    blueprint_data_t *restrict final_blueprint_data
    """


    body = r"""
    const double intersection_pos[3] = {
        all_photons->source_event_f_intersect[IDX_GLOBAL(1, photon_idx, num_rays)], 
        all_photons->source_event_f_intersect[IDX_GLOBAL(2, photon_idx, num_rays)], 
        all_photons->source_event_f_intersect[IDX_GLOBAL(3, photon_idx, num_rays)]
    };
    
    const double source_plane_center[3] = { commondata->source_plane_center_x, commondata->source_plane_center_y, commondata->source_plane_center_z };
    const double source_plane_normal[3] = { commondata->source_plane_normal_x, commondata->source_plane_normal_y, commondata->source_plane_normal_z };
    const double source_up_vector[3] = { commondata->source_up_vec_x, commondata->source_up_vec_y, commondata->source_up_vec_z };

    double s_z[3] = { source_plane_normal[0], source_plane_normal[1], source_plane_normal[2] };
    double s_x[3];
    s_x[0] = source_up_vector[1]*s_z[2] - source_up_vector[2]*s_z[1];
    s_x[1] = source_up_vector[2]*s_z[0] - source_up_vector[0]*s_z[2];
    s_x[2] = source_up_vector[0]*s_z[1] - source_up_vector[1]*s_z[0];
    
    double mag_s_x = sqrt(s_x[0]*s_x[0] + s_x[1]*s_x[1] + s_x[2]*s_x[2]);

    if (mag_s_x < 1e-9) {
        double alternative_up[3] = {1.0, 0.0, 0.0};
        if (fabs(s_z[0]) > 0.999) { alternative_up[0] = 0.0; alternative_up[1] = 1.0; }
        s_x[0] = alternative_up[1]*s_z[2] - alternative_up[2]*s_z[1];
        s_x[1] = alternative_up[2]*s_z[0] - alternative_up[0]*s_z[2];
        s_x[2] = alternative_up[0]*s_z[1] - alternative_up[1]*s_z[0];
        mag_s_x = sqrt(s_x[0]*s_x[0] + s_x[1]*s_x[1] + s_x[2]*s_x[2]);
    }

    double inv_mag_s_x = 1.0 / mag_s_x;
    s_x[0] *= inv_mag_s_x; s_x[1] *= inv_mag_s_x; s_x[2] *= inv_mag_s_x;

    double s_y[3];
    s_y[0] = s_z[1]*s_x[2] - s_z[2]*s_x[1];
    s_y[1] = s_z[2]*s_x[0] - s_z[0]*s_x[2];
    s_y[2] = s_z[0]*s_x[1] - s_z[1]*s_x[0];

    const double vec_s[3] = {
        intersection_pos[0] - source_plane_center[0],
        intersection_pos[1] - source_plane_center[1],
        intersection_pos[2] - source_plane_center[2]
    };

    const double y_s = vec_s[0]*s_x[0] + vec_s[1]*s_x[1] + vec_s[2]*s_x[2];
    const double z_s = vec_s[0]*s_y[0] + vec_s[1]*s_y[1] + vec_s[2]*s_y[2];

    const double r_s_sq = y_s*y_s + z_s*z_s;
    const double r_min_sq = commondata->source_r_min * commondata->source_r_min;
    const double r_max_sq = commondata->source_r_max * commondata->source_r_max;

    if (r_s_sq >= r_min_sq && r_s_sq <= r_max_sq) {
        final_blueprint_data->termination_type = TERMINATION_TYPE_SOURCE_PLANE;
        final_blueprint_data->y_s = y_s;
        final_blueprint_data->z_s = z_s;
        
        // UPDATED: Pull time and Lambda directly from the SoA arrays
        final_blueprint_data->t_s = all_photons->source_event_f_intersect[IDX_GLOBAL(0, photon_idx, num_rays)];
        final_blueprint_data->L_s = all_photons->source_event_f_intersect[IDX_GLOBAL(8, photon_idx, num_rays)]; 
        return true;
    }
    return false;
    """

    prefunc = """
    #ifdef USE_GPU
    #pragma omp declare target
    #endif
    """
        
    postfunc = """
    #ifdef USE_GPU
    #pragma omp end declare target
    #endif
    """

    cfc.register_CFunction(
        includes=includes, 
        desc=desc, 
        name=name, 
        cfunc_type=cfunc_type, 
        params=params, 
        body=body,
        prefunc=prefunc,  
        postfunc=postfunc  
    )