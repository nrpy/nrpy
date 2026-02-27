"""
Generates the C orchestrator for geometric event detection.

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module provides the high-level logic for detecting crossings of the
observer window and the source emission plane using sparse branching
to optimize GPU thread efficiency.
"""

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h

# Define the C-code for event types and geometry parameters
event_mgmt_c_code = r"""
    /**
     * @brief Defines the geometric parameters for a plane in 3D space.
     * * The plane is defined by the equation: normal[0]*x + normal[1]*y + normal[2]*z - d = 0.
     */
    typedef struct {
        double normal[3]; // The unit normal vector describing the orientation of the plane.
        double d;         // The scalar distance from the origin to the plane along the normal.
    } plane_event_params_t;
"""

# Register definitions to the BHaH header system
Bdefines_h.register_BHaH_defines("photon_03_event_management", event_mgmt_c_code)


def event_detection_manager() -> None:
    """
    Register the event_detection_manager C function.

    This function performs the geometric checks for plane crossings. It is designed
    to be called within the main integration loop for every active photon.

    >>> event_detection_manager()
    """
    # 1. Define C-Function metadata
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"]
    desc = """@brief GPU-optimized detection of plane crossings using SoA data.
    Algorithm: Evaluates the scalar product of the photon's position and the
    plane normal. A transition in the sign of this product indicates a crossing,
    triggering the high-precision root-finder."""
    name = "event_detection_manager"

    params = """
        PhotonStateSoA *restrict all_photons,
        const long int num_rays,
        const long int photon_idx,
        const commondata_struct *restrict commondata
        """

    # 2. Build the C body with the Preamble Pattern
    body = r"""
    // --- Preamble: Descriptive Physical Variable Mapping ---
    // Extract spatial coordinates from the flattened Structure of Arrays (SoA).
    const double x = all_photons->f[IDX_GLOBAL(1, photon_idx, num_rays)]; // Cartesian x-coordinate of the photon.
    const double y = all_photons->f[IDX_GLOBAL(2, photon_idx, num_rays)]; // Cartesian y-coordinate of the photon.
    const double z = all_photons->f[IDX_GLOBAL(3, photon_idx, num_rays)]; // Cartesian z-coordinate of the photon.

    // --- Window Plane Detection Logic ---
    if (!all_photons->window_event_found[photon_idx]) {
        // Calculate the window plane normal vector based on camera and window centers.
        double w_normal[3];
        w_normal[0] = commondata->window_center_x - commondata->camera_pos_x; // Normal component along the x-axis.
        w_normal[1] = commondata->window_center_y - commondata->camera_pos_y; // Normal component along the y-axis.
        w_normal[2] = commondata->window_center_z - commondata->camera_pos_z; // Normal component along the z-axis.

        // Normalize the vector to ensure the plane evaluation is purely geometric.
        const double mag_inv = 1.0 / sqrt(w_normal[0]*w_normal[0] + w_normal[1]*w_normal[1] + w_normal[2]*w_normal[2]);
        for(int i=0; i<3; i++) w_normal[i] *= mag_inv;

        // Calculate the scalar distance (d) from the origin to the window plane.
        const double w_dist = commondata->window_center_x*w_normal[0] + commondata->window_center_y*w_normal[1] + commondata->window_center_z*w_normal[2];

        // Evaluate the plane equation at the current photon position.
        const double w_val = x*w_normal[0] + y*w_normal[1] + z*w_normal[2] - w_dist;
        const bool on_pos_curr = (w_val > 1e-10); // True if the photon is on the 'positive' side of the normal.

        // Zero-Crossing Detection: Compare current state with the previous integration step.
        if (on_pos_curr != all_photons->on_positive_side_of_window_prev[photon_idx]) {
             find_event_time_and_state(all_photons, num_rays, photon_idx, w_normal, w_dist, WINDOW_EVENT);
        }
        all_photons->on_positive_side_of_window_prev[photon_idx] = on_pos_curr;
    }

    // --- Source Plane Detection Logic ---
    if (!all_photons->source_event_found[photon_idx]) {
        // Unpack source plane normal from common data.
        const double s_normal[3] = {commondata->source_plane_normal_x, commondata->source_plane_normal_y, commondata->source_plane_normal_z};

        // Calculate the scalar distance (d) to the physical emission source plane.
        const double s_dist = commondata->source_plane_center_x*s_normal[0] + commondata->source_plane_center_y*s_normal[1] + commondata->source_plane_center_z*s_normal[2];

        // Evaluate the plane equation at the current photon position.
        const double s_val = x*s_normal[0] + y*s_normal[1] + z*s_normal[2] - s_dist;
        const bool on_pos_curr = (s_val > 1e-10);

        if (on_pos_curr != all_photons->on_positive_side_of_source_prev[photon_idx]) {
             find_event_time_and_state(all_photons, num_rays, photon_idx, s_normal, s_dist, SOURCE_EVENT);
        }
        all_photons->on_positive_side_of_source_prev[photon_idx] = on_pos_curr;
    }
    """

    prefunc = "#ifdef USE_GPU\n#pragma omp declare target\n#endif"
    postfunc = "#ifdef USE_GPU\n#pragma omp end declare target\n#endif"

    # 3. Register the C function
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        postfunc=postfunc,
    )
