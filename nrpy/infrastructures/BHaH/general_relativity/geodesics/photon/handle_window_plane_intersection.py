"""
Generates the C engine to handle a window plane intersection.

This module calculates the local 2D coordinates on the observer's camera window
when a photon crosses the window plane. It relies entirely on constant memory
to preserve hardware register constraints.

Author: Dalton J. Moone.
"""
import nrpy.c_function as cfc
import nrpy.params as par


def handle_window_plane_intersection() -> None:
    """
    Generate and register the C engine for processing window plane intersections.

    :raises SystemError: If C function registration fails during the pipeline compilation.
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

    includes = [
        "BHaH_defines.h",
        "BHaH_device_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdbool.h>",
        "cuda_intrinsics.h"
    ]

    desc = r"""@brief Processes a window plane intersection without terminating the trajectory.

    @param f_local Thread-local array containing the 9-component photon state $f^\mu$.
    @param final_blueprint_data Pointer to the blueprint structure $b_i$ for data persistence.

    Algorithm:
    1. Reconstructs the orthonormal camera basis ($w_x$, $w_y$, $w_z$).
    2. Projects the 3D intersection point $x^i$ onto the local window axes.
    3. Validates if the intersection falls within the physical window boundaries."""

    cfunc_type = "BHAH_HD_FUNC bool"
    name = "handle_window_plane_intersection"

    params = (
        "const double *restrict f_local, "
        "blueprint_data_t *restrict final_blueprint_data"
    )

    include_CodeParameters_h = False

    body = r"""
    // --- STATE UNPACKING ---
    // Maps global states to immediate thread-local variables to minimize VRAM latency.
    const double t_intersect = f_local[0]; // Coordinate time $t$ at intersection.
    const double x_intersect = f_local[1]; // Cartesian $x$ at intersection.
    const double y_intersect = f_local[2]; // Cartesian $y$ at intersection.
    const double z_intersect = f_local[3]; // Cartesian $z$ at intersection.
    const double L_intersect = f_local[8]; // Affine parameter $\lambda$ at intersection.

    // --- CAMERA BASIS RECONSTRUCTION ---
    // Reconstructs the orthonormal frame to map the physical intersection onto the virtual pixel grid.
    const double window_center[3] = {d_commondata.window_center_x, d_commondata.window_center_y, d_commondata.window_center_z}; 

    double w_z[3] = {
        d_commondata.window_center_x - d_commondata.camera_pos_x,
        d_commondata.window_center_y - d_commondata.camera_pos_y,
        d_commondata.window_center_z - d_commondata.camera_pos_z
    }; 

    double mag_w_z = SqrtCUDA(w_z[0]*w_z[0] + w_z[1]*w_z[1] + w_z[2]*w_z[2]); 
    if (mag_w_z > 1e-12) {
        double inv_mag = 1.0 / mag_w_z; // Inverse multiplication is utilized to optimize floating-point arithmetic speed.
        w_z[0] *= inv_mag; w_z[1] *= inv_mag; w_z[2] *= inv_mag;
    }

    const double window_up[3] = {d_commondata.window_up_vec_x, d_commondata.window_up_vec_y, d_commondata.window_up_vec_z}; 

    double w_x[3]; // Horizontal basis vector $w_x$ (orthogonal to $w_z$ and 'up').
    w_x[0] = window_up[1]*w_z[2] - window_up[2]*w_z[1];
    w_x[1] = window_up[2]*w_z[0] - window_up[0]*w_z[2];
    w_x[2] = window_up[0]*w_z[1] - window_up[1]*w_z[0];

    double mag_w_x = SqrtCUDA(w_x[0]*w_x[0] + w_x[1]*w_x[1] + w_x[2]*w_x[2]); 

    if (mag_w_x < 1e-9) {
        double alt_up[3] = {1.0, 0.0, 0.0}; 
        if (AbsCUDA(w_z[0]) > 0.999) { alt_up[0] = 0.0; alt_up[1] = 1.0; }
        w_x[0] = alt_up[1]*w_z[2] - alt_up[2]*w_z[1];
        w_x[1] = alt_up[2]*w_z[0] - alt_up[0]*w_z[2];
        w_x[2] = alt_up[0]*w_z[1] - alt_up[1]*w_z[0];
        mag_w_x = SqrtCUDA(w_x[0]*w_x[0] + w_x[1]*w_x[1] + w_x[2]*w_x[2]);
    }

    double inv_mag_w_x = 1.0 / mag_w_x;
    w_x[0] *= inv_mag_w_x; w_x[1] *= inv_mag_w_x; w_x[2] *= inv_mag_w_x;

    double w_y[3]; // Vertical basis vector $w_y$ completing the right-handed frame.
    w_y[0] = w_z[1]*w_x[2] - w_z[2]*w_x[1];
    w_y[1] = w_z[2]*w_x[0] - w_z[0]*w_x[2];
    w_y[2] = w_z[0]*w_x[1] - w_z[1]*w_x[0];

    // --- PROJECTION & VALIDATION ---
    // Transforms the global spatial intersection into the local 2D window coordinate system.
    const double relative_pos[3] = {
        x_intersect - window_center[0],
        y_intersect - window_center[1],
        z_intersect - window_center[2]
    }; 

    const double local_y_w = relative_pos[0]*w_x[0] + relative_pos[1]*w_x[1] + relative_pos[2]*w_x[2]; 
    const double local_z_w = relative_pos[0]*w_y[0] + relative_pos[1]*w_y[1] + relative_pos[2]*w_y[2]; 

    if (AbsCUDA(local_y_w) <= d_commondata.window_width / 2.0 && AbsCUDA(local_z_w) <= d_commondata.window_height / 2.0) {
        final_blueprint_data->y_w = local_y_w; // Store local horizontal offset.
        final_blueprint_data->z_w = local_z_w; // Store local vertical offset.
        final_blueprint_data->t_w = t_intersect; // Persistent coordinate time of crossing.
        final_blueprint_data->L_w = L_intersect; // Persistent affine parameter of crossing.
        return true;
    }

    return false;
    """

    # Enforcing strict Master Order. Omitted prefunc/postfunc to prevent compilation bloat.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body
    )