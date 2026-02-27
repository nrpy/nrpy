"""
Generates the C engine to handle a source plane intersection.

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module computes the impact parameters $(y_s, z_s)$ on the physical
emission plane for use in ray-tracing visualizations.
"""

import nrpy.c_function as cfc
import nrpy.params as par


def handle_source_plane_intersection() -> None:
    """
    Generate and register the C engine for source plane intersection handling.

    This function defines the local coordinate system of the source plane
    and checks if the photon hit the active region (r_min to r_max).
    """
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "source_plane_normal_x",
            "source_plane_normal_y",
            "source_plane_normal_z",
            "source_plane_center_x",
            "source_plane_center_y",
            "source_plane_center_z",
            "source_up_vec_x",
            "source_up_vec_y",
            "source_up_vec_z",
            "source_r_min",
            "source_r_max",
        ],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 6.0, 25.0],
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
    desc = r"""@brief Processes a terminal intersection with the source emission plane.

    Algorithm:
    1. Reconstructs the source plane orthonormal basis (s_x, s_y, s_z).
    2. Projects the intersection into local 2D coordinates.
    3. Filters based on radial bounds [source_r_min, source_r_max]."""

    name = "handle_source_plane_intersection"
    cfunc_type = "bool"
    params = """const PhotonStateSoA *restrict all_photons, long int num_rays,
                long int photon_idx, const commondata_struct *restrict commondata,
                blueprint_data_t *restrict final_blueprint_data"""

    # 2. Build the C body
    body = r"""
    // === Preamble: Unpack Intersection State ===
    const double t_intersect = all_photons->source_event_f_intersect[IDX_GLOBAL(0, photon_idx, num_rays)]; // Coordinate time $t$ at intersection.
    const double x_intersect = all_photons->source_event_f_intersect[IDX_GLOBAL(1, photon_idx, num_rays)]; // Cartesian $x$ at intersection.
    const double y_intersect = all_photons->source_event_f_intersect[IDX_GLOBAL(2, photon_idx, num_rays)]; // Cartesian $y$ at intersection.
    const double z_intersect = all_photons->source_event_f_intersect[IDX_GLOBAL(3, photon_idx, num_rays)]; // Cartesian $z$ at intersection.
    const double L_intersect = all_photons->source_event_f_intersect[IDX_GLOBAL(8, photon_idx, num_rays)]; // Affine parameter $\lambda$ at intersection.

    // === Step 1: Reconstruct Source Plane Basis ===
    const double s_z[3] = { commondata->source_plane_normal_x, commondata->source_plane_normal_y, commondata->source_plane_normal_z }; // Plane normal vector.
    const double source_up[3] = { commondata->source_up_vec_x, commondata->source_up_vec_y, commondata->source_up_vec_z }; // Reference orientation vector.

    double s_x[3]; // Horizontal basis vector $s_x$ for the source plane.
    s_x[0] = source_up[1]*s_z[2] - source_up[2]*s_z[1];
    s_x[1] = source_up[2]*s_z[0] - source_up[0]*s_z[2];
    s_x[2] = source_up[0]*s_z[1] - source_up[1]*s_z[0];

    double mag_s_x = sqrt(s_x[0]*s_x[0] + s_x[1]*s_x[1] + s_x[2]*s_x[2]); // Normalization factor for the $s_x$ axis.

    if (mag_s_x < 1e-9) {
        double alt_up[3] = {1.0, 0.0, 0.0}; // Alternative orientation if primary 'up' is degenerate.
        if (fabs(s_z[0]) > 0.999) { alt_up[0] = 0.0; alt_up[1] = 1.0; }
        s_x[0] = alt_up[1]*s_z[2] - alt_up[2]*s_z[1];
        s_x[1] = alt_up[2]*s_z[0] - alt_up[0]*s_z[2];
        s_x[2] = alt_up[0]*s_z[1] - alt_up[1]*s_z[0];
        mag_s_x = sqrt(s_x[0]*s_x[0] + s_x[1]*s_x[1] + s_x[2]*s_x[2]);
    }

    double inv_mag_s_x = 1.0 / mag_s_x;
    s_x[0] *= inv_mag_s_x; s_x[1] *= inv_mag_s_x; s_x[2] *= inv_mag_s_x;

    double s_y[3]; // Vertical basis vector $s_y$ completing the plane coordinate system.
    s_y[0] = s_z[1]*s_x[2] - s_z[2]*s_x[1];
    s_y[1] = s_z[2]*s_x[0] - s_z[0]*s_x[2];
    s_y[2] = s_z[0]*s_x[1] - s_z[1]*s_x[0];

    // === Step 2: Project into Local Source Coordinates ===
    const double relative_pos[3] = {
        x_intersect - commondata->source_plane_center_x,
        y_intersect - commondata->source_plane_center_y,
        z_intersect - commondata->source_plane_center_z
    }; // Vector from plane center to intersection point.

    const double local_y_s = relative_pos[0]*s_x[0] + relative_pos[1]*s_x[1] + relative_pos[2]*s_x[2]; // Impact parameter along $s_x$.
    const double local_z_s = relative_pos[0]*s_y[0] + relative_pos[1]*s_y[1] + relative_pos[2]*s_y[2]; // Impact parameter along $s_y$.

    // === Step 3: Radial Filtering and Termination ===
    const double r_sq = local_y_s*local_y_s + local_z_s*local_z_s; // Squared radial distance from the source center.
    const double r_min_sq = commondata->source_r_min * commondata->source_r_min; // Lower radial cutoff (e.g., ISCO).
    const double r_max_sq = commondata->source_r_max * commondata->source_r_max; // Upper radial cutoff.

    if (r_sq >= r_min_sq && r_sq <= r_max_sq) {
        final_blueprint_data->termination_type = TERMINATION_TYPE_SOURCE_PLANE; // Mark trajectory as successfully terminated at source.
        final_blueprint_data->y_s = local_y_s;
        final_blueprint_data->z_s = local_z_s;
        final_blueprint_data->t_s = t_intersect;
        final_blueprint_data->L_s = L_intersect;
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
