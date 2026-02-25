"""
Generates the C orchestrator for geometric event detection.
Optimized for GPU: Uses direct indexing and handles sparse branching logic.
"""

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h

# Define the C-code for event types and geometry
event_mgmt_c_code = r"""
    // Function pointer type for event functions (kept for architectural compatibility)
    typedef double (*event_function_t)(const double *f, long int num_rays, long int photon_idx, const void *params);

    typedef struct {
        double normal[3];
        double d;
    } plane_event_params;
"""

Bdefines_h.register_BHaH_defines("photon_03_event_management", event_mgmt_c_code)

def event_detection_manager() -> None:
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>"]
    desc = "@brief GPU-optimized detection of plane crossings using SoA data."
    name = "event_detection_manager"

    params = """
        PhotonStateSoA *restrict all_photons, 
        const long int num_rays, 
        const long int photon_idx, 
        const commondata_struct *restrict commondata
        """

    body = r"""
    // Optimized Macro: Accesses only the 3 spatial coordinates needed for plane check
    // This avoids loading the full 9-component state into registers.
    #define PLANE_VAL(soa_ptr, r_idx, total_rays, n, d) ( \
        (soa_ptr)[IDX_GLOBAL(1, r_idx, total_rays)]*(n)[0] + \
        (soa_ptr)[IDX_GLOBAL(2, r_idx, total_rays)]*(n)[1] + \
        (soa_ptr)[IDX_GLOBAL(3, r_idx, total_rays)]*(n)[2] - (d) )

    // --- Window Plane Detection ---
    if (!all_photons->window_event_found[photon_idx]) {
        // These parameters are ideally loaded from __constant__ memory.
        double w_normal[3] = {
            commondata->window_center_x - commondata->camera_pos_x,
            commondata->window_center_y - commondata->camera_pos_y,
            commondata->window_center_z - commondata->camera_pos_z
        };
        double mag_inv = 1.0 / sqrt(w_normal[0]*w_normal[0] + w_normal[1]*w_normal[1] + w_normal[2]*w_normal[2]);
        for(int i=0; i<3; i++) w_normal[i] *= mag_inv;
        double w_dist = commondata->window_center_x*w_normal[0] + commondata->window_center_y*w_normal[1] + commondata->window_center_z*w_normal[2];

        double w_val = PLANE_VAL(all_photons->f, photon_idx, num_rays, w_normal, w_dist);
        bool on_pos_curr = (w_val > 1e-10); 
        
        // Standard branching: Saves 99.9% of threads from running the root-finder.
        if (on_pos_curr != all_photons->on_positive_side_of_window_prev[photon_idx]) {
             find_event_time_and_state(all_photons, num_rays, photon_idx, w_normal, w_dist, WINDOW_EVENT);
        }
        all_photons->on_positive_side_of_window_prev[photon_idx] = on_pos_curr;
    }

    // --- Source Plane Detection ---
    if (!all_photons->source_event_found[photon_idx]) {
        double s_normal[3] = {commondata->source_plane_normal_x, commondata->source_plane_normal_y, commondata->source_plane_normal_z};
        double s_dist = commondata->source_plane_center_x*s_normal[0] + commondata->source_plane_center_y*s_normal[1] + commondata->source_plane_center_z*s_normal[2];

        double s_val = PLANE_VAL(all_photons->f, photon_idx, num_rays, s_normal, s_dist);
        bool on_pos_curr = (s_val > 1e-10);
        
        if (on_pos_curr != all_photons->on_positive_side_of_source_prev[photon_idx]) {
             find_event_time_and_state(all_photons, num_rays, photon_idx, s_normal, s_dist, SOURCE_EVENT);
        }
        all_photons->on_positive_side_of_source_prev[photon_idx] = on_pos_curr;
    }
    #undef PLANE_VAL
    """

    portable_body = f"""
    #ifdef USE_GPU
    #pragma omp declare target
    #endif
    {body}
    #ifdef USE_GPU
    #pragma omp end declare target
    #endif
    """

    cfc.register_CFunction(includes=includes, desc=desc, name=name, params=params, body=portable_body)