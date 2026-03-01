"""
Generates the C engine to handle a source plane intersection.

This module computes the impact parameters $(y_s, z_s)$ on the physical
emission plane for use in ray-tracing visualizations. It unpacks coordinates
strictly from thread-local memory to avoid slow global memory accesses,
satisfying the 255-register limitation of the NVIDIA Ampere sm_86 architecture.
Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc
import nrpy.params as par


def handle_source_plane_intersection() -> None:
    """
    Generate and register the C engine for source plane intersection handling.

    :param: None
    :raises SystemError: If C function registration fails within the NRPy+ pipeline.
    """
    # Register physical CodeParameters necessary for emission plane coordinates
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

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdbool.h>",
    ]
    desc = """@brief Processes a terminal intersection with the source emission plane entirely in thread-local memory.

    @param source_event_f_intersect Thread-local state array holding the 9-component intersection state $f^\mu$.
    @param commondata Pointer to the globally constant parameters struct.
    @param final_blueprint_data Pointer to the local blueprint structure $b_i$ to store impact geometries.

    Algorithm:
    1. Reconstructs the source plane orthonormal basis $(s_x, s_y, s_z)$.
    2. Projects the intersection state $x^\mu$ into local 2D coordinates.
    3. Filters based on physical radial bounds $[r_{min}, r_{max}]$.
    
    By consuming `source_event_f_intersect` as a local 9-element array instead of reading from 
    `PhotonStateSoA` via `IDX_GLOBAL`, this kernel bypasses VRAM stalls and guarantees compliance 
    with the 255-register limits on Ampere sm_86 hardware."""
    cfunc_type = "BHAH_HD_FUNC bool"
    name = "handle_source_plane_intersection"
    params = (
        "const double *restrict source_event_f_intersect, "
        "const commondata_struct *restrict commondata, "
        "blueprint_data_t *restrict final_blueprint_data"
    )
    include_CodeParameters_h = False
    body = r"""
    // --- PREAMBLE: UNPACK INTERSECTION STATE ---
    // Reads intersection data from thread-local registers, circumventing IDX_GLOBAL global reads.

    const double t_intersect = source_event_f_intersect[0]; // Coordinate time $t$ at intersection.
    const double x_intersect = source_event_f_intersect[1]; // Cartesian $x^1$ at intersection.
    const double y_intersect = source_event_f_intersect[2]; // Cartesian $x^2$ at intersection.
    const double z_intersect = source_event_f_intersect[3]; // Cartesian $x^3$ at intersection.
    const double L_intersect = source_event_f_intersect[8]; // Affine parameter $\lambda$ at intersection.

    // --- STEP 1: RECONSTRUCT SOURCE PLANE BASIS ---
    const double s_z[3] = { commondata->source_plane_normal_x, commondata->source_plane_normal_y, commondata->source_plane_normal_z }; // Plane normal vector $n_i$.
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

    // Normalize the horizontal basis vector $s_x \cdot s_x = 1$.
    double inv_mag_s_x = 1.0 / mag_s_x;
    s_x[0] *= inv_mag_s_x; s_x[1] *= inv_mag_s_x; s_x[2] *= inv_mag_s_x;

    double s_y[3]; // Vertical basis vector $s_y$ completing the local orthogonal system.
    s_y[0] = s_z[1]*s_x[2] - s_z[2]*s_x[1];
    s_y[1] = s_z[2]*s_x[0] - s_z[0]*s_x[2];
    s_y[2] = s_z[0]*s_x[1] - s_z[1]*s_x[0];

    // --- STEP 2: PROJECT INTO LOCAL SOURCE COORDINATES ---
    const double relative_pos[3] = {
        x_intersect - commondata->source_plane_center_x,
        y_intersect - commondata->source_plane_center_y,
        z_intersect - commondata->source_plane_center_z
    }; // Separation vector $\Delta x^i$ from plane center to intersection point.

    const double local_y_s = relative_pos[0]*s_x[0] + relative_pos[1]*s_x[1] + relative_pos[2]*s_x[2]; // Impact parameter along the $s_x$ axis.
    const double local_z_s = relative_pos[0]*s_y[0] + relative_pos[1]*s_y[1] + relative_pos[2]*s_y[2]; // Impact parameter along the $s_y$ axis.

    // --- STEP 3: RADIAL FILTERING AND TERMINATION ---
    const double r_sq = local_y_s*local_y_s + local_z_s*local_z_s; // Squared radial distance $r^2$ from the source center.
    const double r_min_sq = commondata->source_r_min * commondata->source_r_min; // Lower radial cutoff to simulate physical bounding.
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

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body
    )