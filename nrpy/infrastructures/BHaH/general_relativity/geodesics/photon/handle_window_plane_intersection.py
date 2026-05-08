# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/handle_window_plane_intersection.py
"""
Defines the C engine that handles a window plane intersection.

This module provides the metaprogramming orchestration that defines the C function
responsible for calculating the local 2D coordinates on the observer's camera window
when a photon crosses the window plane. It maps global Cartesian intersections to a
reconstructed orthonormal camera basis by projecting the relative 3D position onto
the horizontal and vertical basis vectors. Upon a valid bounds check, it permanently
records the spatial coordinates, coordinate time, and affine parameter into the
persistent blueprint structure.

Reconstructing the orthonormal frame using the immutable global original window center
prevents perspective warping or rotation across peripheral tiles. The mathematical
implementation optimizes division overhead by utilizing inverse multiplication for
vector normalization. It also mitigates cross-product singularities by failing over to
an alternative up-vector when necessary. Furthermore, the dynamically shifted local
window center acts as the mapping anchor to place intersections within the active tile
boundaries, while global state values are mapped to immediate local variables to
minimize memory latency.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par


def handle_window_plane_intersection() -> None:
    """Generate and register the C engine for processing window plane intersections."""
    parallelization = par.parval_from_str("parallelization")

    # Add the access variable
    cd_access = parallel_utils.get_commondata_access(parallelization)

    # Define inclusion headers for the C function compilation.
    includes = ["BHaH_defines.h"]

    if parallelization == "cuda":
        includes.append("BHaH_device_defines.h")
        includes.append("cuda_intrinsics.h")

    # Define the Doxygen-formatted description for the C function.
    desc = r""" Processes a window plane intersection without terminating the trajectory.

    @param f_local Thread-local array containing the 9-component photon state $f^\mu$.
    @param lam_intersect The explicit affine parameter $\lambda$ of the intersection.
    @param final_blueprint_data Pointer to the blueprint structure $b_i$ for data persistence.

    Algorithm:
    1. Reconstructs the orthonormal camera basis ($w_x$, $w_y$, $w_z$).
    2. Projects the 3D intersection point $x^i$ onto the local window axes.
    3. Validates if the intersection falls within the physical window boundaries."""

    # Specify the signature components of the C function.
    cfunc_type = "BHAH_HD_INLINE bool"
    name = "handle_window_plane_intersection"

    # Define the specific C arguments passed into the kernel.
    params = (
        "const double *restrict f_local, "
        "const double lam_intersect, "
        "blueprint_data_t *restrict final_blueprint_data"
    )
    if parallelization != "cuda":
        params += ", const commondata_struct *restrict commondata"

    # Toggle generation of CodeParameters.h inclusion.
    include_CodeParameters_h = False

    # Construct the C-string for the function body.
    body = r"""
    //==========================================
    // STATE UNPACKING
    //==========================================
    // Maps global states to immediate thread-local variables to minimize VRAM latency.
    const double t_intersect = f_local[0]; // Coordinate time $t$ at intersection.
    const double x_intersect = f_local[1]; // Cartesian $x$ at intersection.
    const double y_intersect = f_local[2]; // Cartesian $y$ at intersection.
    const double z_intersect = f_local[3]; // Cartesian $z$ at intersection.

    //==========================================
    // CAMERA BASIS RECONSTRUCTION
    //==========================================
    // Reconstructs the orthonormal frame using the immutable global original_window_center
    // to prevent perspective warping or rotation across peripheral tiles.
    double w_z[3] = {
        d_commondata.original_window_center_x - d_commondata.camera_pos_x,
        d_commondata.original_window_center_y - d_commondata.camera_pos_y,
        d_commondata.original_window_center_z - d_commondata.camera_pos_z
    }; // Global line-of-sight vector pointing from the camera position to the original window center.

    double mag_w_z = SqrtCUDA(w_z[0]*w_z[0] + w_z[1]*w_z[1] + w_z[2]*w_z[2]); // Magnitude of the $w_z$ vector.
    if (mag_w_z > 1e-12) {
        double inv_mag = 1.0 / mag_w_z; // Inverse multiplication optimizes floating-point arithmetic speed.
        w_z[0] *= inv_mag; w_z[1] *= inv_mag; w_z[2] *= inv_mag;
    } // END IF: mag_w_z > 1e-12

    const double window_up[3] = {d_commondata.window_up_vec_x, d_commondata.window_up_vec_y, d_commondata.window_up_vec_z}; // Camera up-vector defining vertical orientation.

    double w_x[3]; // Horizontal basis vector $w_x$ orthogonal to $w_z$ and up-vector.
    w_x[0] = window_up[1]*w_z[2] - window_up[2]*w_z[1];
    w_x[1] = window_up[2]*w_z[0] - window_up[0]*w_z[2];
    w_x[2] = window_up[0]*w_z[1] - window_up[1]*w_z[0];

    double mag_w_x = SqrtCUDA(w_x[0]*w_x[0] + w_x[1]*w_x[1] + w_x[2]*w_x[2]); // Magnitude of the horizontal basis vector $w_x$.

    if (mag_w_x < 1e-9) {
        double alt_up[3] = {1.0, 0.0, 0.0}; // Alternative up-vector mitigates cross product singularities.
        if (AbsCUDA(w_z[0]) > 0.999) { alt_up[0] = 0.0; alt_up[1] = 1.0; }
        w_x[0] = alt_up[1]*w_z[2] - alt_up[2]*w_z[1];
        w_x[1] = alt_up[2]*w_z[0] - alt_up[0]*w_z[2];
        w_x[2] = alt_up[0]*w_z[1] - alt_up[1]*w_z[0];
        mag_w_x = SqrtCUDA(w_x[0]*w_x[0] + w_x[1]*w_x[1] + w_x[2]*w_x[2]);
    } // END IF: mag_w_x < 1e-9

    double inv_mag_w_x = 1.0 / mag_w_x; // Inverse magnitude of $w_x$ optimizes division overhead.
    w_x[0] *= inv_mag_w_x; w_x[1] *= inv_mag_w_x; w_x[2] *= inv_mag_w_x;

    double w_y[3]; // Vertical basis vector $w_y$ completing the right-handed frame.
    w_y[0] = w_z[1]*w_x[2] - w_z[2]*w_x[1];
    w_y[1] = w_z[2]*w_x[0] - w_z[0]*w_x[2];
    w_y[2] = w_z[0]*w_x[1] - w_z[1]*w_x[0];

    //==========================================
    // PROJECTION & VALIDATION
    //==========================================
    // Transforms the global spatial intersection into the local 2D window coordinate system.
    // The dynamically shifted local_window_center is used here as the anchor to map hits within the active tile bounds.
    const double local_window_center[3] = {d_commondata.window_center_x, d_commondata.window_center_y, d_commondata.window_center_z};

    const double relative_pos[3] = {
        x_intersect - local_window_center[0],
        y_intersect - local_window_center[1],
        z_intersect - local_window_center[2]
    }; // Vector pointing from the local tile window center to the 3D intersection.

    const double local_y_w = relative_pos[0]*w_x[0] + relative_pos[1]*w_x[1] + relative_pos[2]*w_x[2]; // Projected local horizontal coordinate on the window.
    const double local_z_w = relative_pos[0]*w_y[0] + relative_pos[1]*w_y[1] + relative_pos[2]*w_y[2]; // Projected local vertical coordinate on the window.

    if (AbsCUDA(local_y_w) <= d_commondata.window_width / 2.0 && AbsCUDA(local_z_w) <= d_commondata.window_height / 2.0) {
        final_blueprint_data->y_w = local_y_w; // Store local horizontal offset in persistent struct.
        final_blueprint_data->z_w = local_z_w; // Store local vertical offset in persistent struct.
        final_blueprint_data->t_w = t_intersect; // Persistent coordinate time of crossing.
        final_blueprint_data->L_w = lam_intersect; // Explicit $\lambda$ mapped to persistent blueprint.
        return true;
    } // END IF: intersection is within active window bounds

    return false;
    """

    # Inject the string replacement right before registration
    body = body.replace("d_commondata.", cd_access)

    # Register the C function using the defined components.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
