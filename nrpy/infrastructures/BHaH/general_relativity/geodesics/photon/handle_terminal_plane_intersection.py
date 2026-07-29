# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/handle_terminal_plane_intersection.py
"""
Defines the computation for impact parameters on the physical terminal plane.

This module constructs the C engine responsible for evaluating terminal interactions
with the terminal surface. It accepts the fully interpolated nine-component tensor
state and affine parameter at the exact boundary crossing alongside a pointer to the
output blueprint structure. The mathematical mechanism constructs a local orthonormal
basis using the parameterized terminal plane normal and up-vectors, with a fallback
strategy for degenerate cross products. The global Cartesian intersection coordinates
are projected into this local two-dimensional coordinate system and filtered against
the minimum and maximum radial disk bounds. Valid impacts are stored persistently in
the blueprint structures. It relies on constant memory caching to process physical
impact parameters.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par


def register_terminal_plane_parameters() -> None:
    """Register the independent terminal-plane geometry parameters."""
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "terminal_plane_normal_x",
            "terminal_plane_normal_y",
            "terminal_plane_normal_z",
            "terminal_plane_center_x",
            "terminal_plane_center_y",
            "terminal_plane_center_z",
            "terminal_plane_up_x",
            "terminal_plane_up_y",
            "terminal_plane_up_z",
            "terminal_plane_min_coord_radius",
            "terminal_plane_max_coord_radius",
        ],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 6.0, 25.0],
        commondata=True,
        add_to_parfile=True,
    )


def handle_terminal_plane_intersection() -> None:
    """Register the C engine for terminal plane intersection handling."""
    parallelization = par.parval_from_str("parallelization")

    register_terminal_plane_parameters()

    # Add the access variable
    cd_access = parallel_utils.get_commondata_access(parallelization)

    includes = ["BHaH_defines.h"]

    if parallelization == "cuda":
        includes.append("BHaH_device_defines.h")
        includes.append("cuda_intrinsics.h")

    desc = r""" Processes a terminal intersection with the configured terminal plane.

    @param[in] terminal_plane_event_f_intersect Thread-local intersection state; slot zero contains coordinate time.
    @param physical_lambda Physical affine parameter $\lambda$ at the intersection.
    @param[in,out] final_blueprint_data Local blueprint structure $b_i$ updated with impact geometry.
    @param[in] commondata Host-only commondata pointer for non-CUDA builds.

    Algorithm:
    1. Reconstructs the terminal plane orthonormal basis ($s_x$, $s_y$, $s_z$).
    2. Projects the intersection state $x^\mu$ into local 2D coordinates.
    3. Filters based on terminal-plane coordinate-radius bounds
       [$r_{min}$, $r_{max}$]."""

    cfunc_type = "BHAH_HD_INLINE bool"

    name = "handle_terminal_plane_intersection"

    params = (
        "const double *restrict terminal_plane_event_f_intersect, "
        "const double physical_lambda, "
        "blueprint_data_t *restrict final_blueprint_data"
    )
    if parallelization != "cuda":
        params += ", const commondata_struct *restrict commondata"

    include_CodeParameters_h = False

    body = r"""
    //==========================================
    // UNPACK INTERSECTION STATE
    //==========================================
    // Reads intersection data from thread-local registers, circumventing global reads.
    const double t_intersect = terminal_plane_event_f_intersect[0]; // Coordinate time $t$ at intersection.
    const double x_intersect = terminal_plane_event_f_intersect[1]; // Cartesian $x^1$ at intersection.
    const double y_intersect = terminal_plane_event_f_intersect[2]; // Cartesian $x^2$ at intersection.
    const double z_intersect = terminal_plane_event_f_intersect[3]; // Cartesian $x^3$ at intersection.

    //==========================================
    // RECONSTRUCT TERMINAL-PLANE BASIS
    //==========================================
    // Normalize the independent terminal-plane normal before constructing the
    // local basis. This keeps local coordinates and radial bounds invariant
    // under rescaling of the configured normal.
    double terminal_plane_normal[3] = {
        d_commondata.terminal_plane_normal_x,
        d_commondata.terminal_plane_normal_y,
        d_commondata.terminal_plane_normal_z
    };
    const double normal_sq =
        terminal_plane_normal[0] * terminal_plane_normal[0] +
        terminal_plane_normal[1] * terminal_plane_normal[1] +
        terminal_plane_normal[2] * terminal_plane_normal[2];
    if (!isfinite(normal_sq) || normal_sq <= 1.0e-28) {
        return false;
    }
    const double inverse_normal_norm = 1.0 / SqrtCUDA(normal_sq);
    terminal_plane_normal[0] *= inverse_normal_norm;
    terminal_plane_normal[1] *= inverse_normal_norm;
    terminal_plane_normal[2] *= inverse_normal_norm;

    // Orthogonalize the supplied up seed against the plane normal. If it is
    // parallel to that normal, use a coordinate-axis fallback and repeat the
    // same projection.
    double terminal_plane_up[3] = {
        d_commondata.terminal_plane_up_x,
        d_commondata.terminal_plane_up_y,
        d_commondata.terminal_plane_up_z
    };
    double up_dot_normal =
        terminal_plane_up[0] * terminal_plane_normal[0] +
        terminal_plane_up[1] * terminal_plane_normal[1] +
        terminal_plane_up[2] * terminal_plane_normal[2];
    terminal_plane_up[0] -= up_dot_normal * terminal_plane_normal[0];
    terminal_plane_up[1] -= up_dot_normal * terminal_plane_normal[1];
    terminal_plane_up[2] -= up_dot_normal * terminal_plane_normal[2];

    double up_sq =
        terminal_plane_up[0] * terminal_plane_up[0] +
        terminal_plane_up[1] * terminal_plane_up[1] +
        terminal_plane_up[2] * terminal_plane_up[2];
    if (!isfinite(up_sq) || up_sq <= 1.0e-18) {
        const double fallback_axis[3] = {
            AbsCUDA(terminal_plane_normal[0]) < 0.9 ? 1.0 : 0.0,
            AbsCUDA(terminal_plane_normal[1]) < 0.9 ? 1.0 : 0.0,
            AbsCUDA(terminal_plane_normal[2]) < 0.9 ? 1.0 : 0.0
        };
        up_dot_normal =
            fallback_axis[0] * terminal_plane_normal[0] +
            fallback_axis[1] * terminal_plane_normal[1] +
            fallback_axis[2] * terminal_plane_normal[2];
        terminal_plane_up[0] = fallback_axis[0] - up_dot_normal * terminal_plane_normal[0];
        terminal_plane_up[1] = fallback_axis[1] - up_dot_normal * terminal_plane_normal[1];
        terminal_plane_up[2] = fallback_axis[2] - up_dot_normal * terminal_plane_normal[2];
        up_sq =
            terminal_plane_up[0] * terminal_plane_up[0] +
            terminal_plane_up[1] * terminal_plane_up[1] +
            terminal_plane_up[2] * terminal_plane_up[2];
    }
    if (!isfinite(up_sq) || up_sq <= 1.0e-18) {
        return false;
    }
    const double inverse_up_norm = 1.0 / SqrtCUDA(up_sq);
    terminal_plane_up[0] *= inverse_up_norm;
    terminal_plane_up[1] *= inverse_up_norm;
    terminal_plane_up[2] *= inverse_up_norm;

    // right = up x normal. The stored local coordinates are (right, up).
    const double terminal_plane_right[3] = {
        terminal_plane_up[1] * terminal_plane_normal[2] -
            terminal_plane_up[2] * terminal_plane_normal[1],
        terminal_plane_up[2] * terminal_plane_normal[0] -
            terminal_plane_up[0] * terminal_plane_normal[2],
        terminal_plane_up[0] * terminal_plane_normal[1] -
            terminal_plane_up[1] * terminal_plane_normal[0]
    };

    //==========================================
    // PROJECT INTO LOCAL TERMINAL-PLANE COORDINATES
    //==========================================
    // Projects the global Cartesian intersection state into the local 2D terminal plane.
    const double relative_pos[3] = {
        x_intersect - d_commondata.terminal_plane_center_x,
        y_intersect - d_commondata.terminal_plane_center_y,
        z_intersect - d_commondata.terminal_plane_center_z
    }; // Relative displacement vector from the terminal plane center.

    const double y_t =
        relative_pos[0] * terminal_plane_right[0] +
        relative_pos[1] * terminal_plane_right[1] +
        relative_pos[2] * terminal_plane_right[2]; // Local width coordinate.
    const double z_t =
        relative_pos[0] * terminal_plane_up[0] +
        relative_pos[1] * terminal_plane_up[1] +
        relative_pos[2] * terminal_plane_up[2]; // Local height coordinate.

    //==========================================
    // RADIAL FILTERING AND TERMINATION
    //==========================================
    // Filters photons that fall outside the physical bounds of the accretion disk.
    const double r_sq = y_t*y_t + z_t*z_t; // Squared radial distance $r^2$ from the terminal-plane center.
    const double r_min = d_commondata.terminal_plane_min_coord_radius;
    const double r_max = d_commondata.terminal_plane_max_coord_radius;
    if (!isfinite(r_min) || !isfinite(r_max) || r_min < 0.0 || r_max < r_min) {
        return false;
    }
    const double r_min_sq = r_min * r_min; // Squared minimum radius.
    const double r_max_sq = r_max * r_max; // Squared maximum radius.

    if (r_sq >= r_min_sq && r_sq <= r_max_sq) {
        final_blueprint_data->termination_type = STOP_CONDITION_TERMINAL_PLANE; // Flags the ray as stopped on the terminal plane.
        final_blueprint_data->y_t = y_t;
        final_blueprint_data->z_t = z_t;
        // For terminal-plane hits, the intersection is the termination event.
        final_blueprint_data->t_f = t_intersect;
        final_blueprint_data->L_f = physical_lambda; // Explicit $\lambda$ mapped to persistent blueprint.
        return true;
    } // END IF: radial filtering
    return false;
    """

    # Inject the string replacement right before registration
    body = body.replace("d_commondata.", cd_access)

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
