"""
Generates the C engine to handle a window plane intersection.

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module calculates the local 2D coordinates on the observer's camera window
when a photon crosses the window plane.
"""

import nrpy.c_function as cfc
import nrpy.params as par


def handle_window_plane_intersection() -> None:
    """
    Generate and register the C engine for processing window plane intersections.

    This function defines the geometric transformation from global Cartesian
    coordinates to the local camera frame (n_x, n_y, n_z).
    """
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "window_center_x",
            "window_center_y",
            "window_center_z",
            "camera_pos_x",
            "camera_pos_y",
            "camera_pos_z",
            "window_up_vec_x",
            "window_up_vec_y",
            "window_up_vec_z",
            "window_width",
            "window_height",
        ],
        [50.0, 0.0, 0.0, 51.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        commondata=True,
        add_to_parfile=True,
    )

    # 1. Define C-Function metadata
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdbool.h>",
    ]
    desc = r"""@brief Processes a window plane intersection without terminating the trajectory.

    Algorithm:
    1. Reconstructs the orthonormal camera basis (w_x, w_y, w_z).
    2. Projects the 3D intersection point onto the local window axes.
    3. Validates if the intersection falls within the physical window boundaries."""

    name = "handle_window_plane_intersection"
    cfunc_type = "bool"
    params = """const PhotonStateSoA *restrict all_photons, const long int num_rays,
                const long int photon_idx, const commondata_struct *restrict commondata,
                blueprint_data_t *restrict final_blueprint_data"""

    # 2. Build the C body with the Preamble Pattern and descriptive comments
    body = r"""
    // === Preamble: Unpack Intersection State ===
    // Retrieve the exact 4-vector position at the moment of plane crossing.
    const double t_intersect = all_photons->window_event_f_intersect[IDX_GLOBAL(0, photon_idx, num_rays)]; // Coordinate time $t$ at intersection.
    const double x_intersect = all_photons->window_event_f_intersect[IDX_GLOBAL(1, photon_idx, num_rays)]; // Cartesian $x$ at intersection.
    const double y_intersect = all_photons->window_event_f_intersect[IDX_GLOBAL(2, photon_idx, num_rays)]; // Cartesian $y$ at intersection.
    const double z_intersect = all_photons->window_event_f_intersect[IDX_GLOBAL(3, photon_idx, num_rays)]; // Cartesian $z$ at intersection.
    const double L_intersect = all_photons->window_event_f_intersect[IDX_GLOBAL(8, photon_idx, num_rays)]; // Affine parameter $\lambda$ at intersection.

    // === Step 1: Reconstruct Orthonormal Camera Basis ===
    const double window_center[3] = {commondata->window_center_x, commondata->window_center_y, commondata->window_center_z}; // Geometric center of the camera window.

    double w_z[3] = {
        commondata->window_center_x - commondata->camera_pos_x,
        commondata->window_center_y - commondata->camera_pos_y,
        commondata->window_center_z - commondata->camera_pos_z
    }; // Normal vector $w_z$ pointing from the camera toward the window center.

    double mag_w_z = sqrt(w_z[0]*w_z[0] + w_z[1]*w_z[1] + w_z[2]*w_z[2]); // Magnitude of the normal vector for normalization.
    if (mag_w_z > 1e-12) {
        double inv_mag = 1.0 / mag_w_z;
        w_z[0] *= inv_mag; w_z[1] *= inv_mag; w_z[2] *= inv_mag;
    }

    const double window_up[3] = {commondata->window_up_vec_x, commondata->window_up_vec_y, commondata->window_up_vec_z}; // User-defined 'up' vector for the camera frame.

    double w_x[3]; // Horizontal basis vector $w_x$ (orthogonal to $w_z$ and 'up').
    w_x[0] = window_up[1]*w_z[2] - window_up[2]*w_z[1];
    w_x[1] = window_up[2]*w_z[0] - window_up[0]*w_z[2];
    w_x[2] = window_up[0]*w_z[1] - window_up[1]*w_z[0];

    double mag_w_x = sqrt(w_x[0]*w_x[0] + w_x[1]*w_x[1] + w_x[2]*w_x[2]); // Magnitude of the horizontal basis vector.

    if (mag_w_x < 1e-9) {
        double alt_up[3] = {1.0, 0.0, 0.0}; // Fallback vector if the primary 'up' is parallel to $w_z$.
        if (fabs(w_z[0]) > 0.999) { alt_up[0] = 0.0; alt_up[1] = 1.0; }
        w_x[0] = alt_up[1]*w_z[2] - alt_up[2]*w_z[1];
        w_x[1] = alt_up[2]*w_z[0] - alt_up[0]*w_z[2];
        w_x[2] = alt_up[0]*w_z[1] - alt_up[1]*w_z[0];
        mag_w_x = sqrt(w_x[0]*w_x[0] + w_x[1]*w_x[1] + w_x[2]*w_x[2]);
    }

    double inv_mag_w_x = 1.0 / mag_w_x;
    w_x[0] *= inv_mag_w_x; w_x[1] *= inv_mag_w_x; w_x[2] *= inv_mag_w_x;

    double w_y[3]; // Vertical basis vector $w_y$ completing the right-handed frame.
    w_y[0] = w_z[1]*w_x[2] - w_z[2]*w_x[1];
    w_y[1] = w_z[2]*w_x[0] - w_z[0]*w_x[2];
    w_y[2] = w_z[0]*w_x[1] - w_z[1]*w_x[0];

    // === Step 2: Project Intersection into Local Window Frame ===
    const double relative_pos[3] = {
        x_intersect - window_center[0],
        y_intersect - window_center[1],
        z_intersect - window_center[2]
    }; // Vector from window center to intersection point.

    const double local_y_w = relative_pos[0]*w_x[0] + relative_pos[1]*w_x[1] + relative_pos[2]*w_x[2]; // Projected horizontal coordinate.
    const double local_z_w = relative_pos[0]*w_y[0] + relative_pos[1]*w_y[1] + relative_pos[2]*w_y[2]; // Projected vertical coordinate.

    // === Step 3: Bounds Validation and Data Persistence ===
    if (fabs(local_y_w) <= commondata->window_width / 2.0 && fabs(local_z_w) <= commondata->window_height / 2.0) {
        final_blueprint_data->y_w = local_y_w; // Store local horizontal offset.
        final_blueprint_data->z_w = local_z_w; // Store local vertical offset.
        final_blueprint_data->t_w = t_intersect; // Persistent coordinate time of crossing.
        final_blueprint_data->L_w = L_intersect; // Persistent affine parameter of crossing.
        return true;
    }

    return false;
    """

    prefunc = "#ifdef USE_GPU\n#pragma omp declare target\n#endif"
    postfunc = "#ifdef USE_GPU\n#pragma omp end declare target\n#endif"

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        cfunc_type=cfunc_type,
        params=params,
        body=body,
        prefunc=prefunc,
        postfunc=postfunc,
    )
