# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/handle_non_terminal_plane_intersection.py
"""
Defines the C engine that handles a nonterminal-plane intersection.

This module provides the metaprogramming orchestration that defines the C function
responsible for calculating local 2D coordinates on the nonterminal plane when a photon
crosses that plane. It uses the plane's own center, normal, and up seed; it does not
derive its basis from observer setup and it never applies a finite rectangle filter.
The crossing is recorded as a diagnostic and integration continues.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par


def register_non_terminal_plane_parameters() -> None:
    """Register the independent nonterminal-plane geometry parameters."""
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "non_terminal_plane_center_x",
            "non_terminal_plane_center_y",
            "non_terminal_plane_center_z",
            "non_terminal_plane_normal_x",
            "non_terminal_plane_normal_y",
            "non_terminal_plane_normal_z",
            "non_terminal_plane_up_x",
            "non_terminal_plane_up_y",
            "non_terminal_plane_up_z",
        ],
        [50.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        commondata=True,
        add_to_parfile=True,
    )


def handle_non_terminal_plane_intersection() -> None:
    """Generate and register the C engine for nonterminal-plane intersections."""
    parallelization = par.parval_from_str("parallelization")

    register_non_terminal_plane_parameters()

    # Add the access variable
    cd_access = parallel_utils.get_commondata_access(parallelization)

    # Define inclusion headers for the C function compilation.
    includes = ["BHaH_defines.h"]

    if parallelization == "cuda":
        includes.append("BHaH_device_defines.h")
        includes.append("cuda_intrinsics.h")

    # Define the Doxygen-formatted description for the C function.
    desc = r""" Processes a nonterminal-plane intersection without terminating the trajectory.

    @param[in] f_local Thread-local intersection state; slot zero contains coordinate time.
    @param physical_lambda Physical affine parameter $\lambda$ at the intersection.
    @param[in,out] final_blueprint_data Blueprint structure $b_i$ updated with the intersection.
    @param[in] commondata Host-only commondata pointer for non-CUDA builds.

    Algorithm:
    1. Reconstructs the independent nonterminal-plane basis ($w_x$, $w_y$, $w_z$).
    2. Projects the 3D intersection point $x^i$ onto local plane axes.
    3. Stores $y_{nt}$, $z_{nt}$, $t_{nt}$, and the affine parameter."""

    # Specify the signature components of the C function.
    cfunc_type = "BHAH_HD_INLINE bool"
    name = "handle_non_terminal_plane_intersection"

    # Define the specific C arguments passed into the kernel.
    params = (
        "const double *restrict f_local, "
        "const double physical_lambda, "
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
    const double t_intersect = f_local[0];
    const double x_intersect = f_local[1];
    const double y_intersect = f_local[2];
    const double z_intersect = f_local[3];

    //==========================================
    // INDEPENDENT NONTERMINAL-PLANE FRAME
    //==========================================
    // The nonterminal plane owns its normal, up seed, and center. It is not
    // derived from observer position or observer look-forward direction.
    double plane_normal[3] = {
        d_commondata.non_terminal_plane_normal_x,
        d_commondata.non_terminal_plane_normal_y,
        d_commondata.non_terminal_plane_normal_z
    };
    const double normal_sq =
        plane_normal[0] * plane_normal[0] +
        plane_normal[1] * plane_normal[1] +
        plane_normal[2] * plane_normal[2];
    if (!isfinite(normal_sq) || normal_sq <= 1.0e-28) {
        return false;
    }
    const double inverse_normal_norm = 1.0 / SqrtCUDA(normal_sq);
    plane_normal[0] *= inverse_normal_norm;
    plane_normal[1] *= inverse_normal_norm;
    plane_normal[2] *= inverse_normal_norm;

    double plane_up[3] = {
        d_commondata.non_terminal_plane_up_x,
        d_commondata.non_terminal_plane_up_y,
        d_commondata.non_terminal_plane_up_z
    };
    const double up_dot_normal =
        plane_up[0] * plane_normal[0] +
        plane_up[1] * plane_normal[1] +
        plane_up[2] * plane_normal[2];
    plane_up[0] -= up_dot_normal * plane_normal[0];
    plane_up[1] -= up_dot_normal * plane_normal[1];
    plane_up[2] -= up_dot_normal * plane_normal[2];

    double up_sq =
        plane_up[0] * plane_up[0] +
        plane_up[1] * plane_up[1] +
        plane_up[2] * plane_up[2];
    if (!isfinite(up_sq) || up_sq <= 1.0e-18) {
        // Select a coordinate-axis fallback least aligned with the supplied
        // normal, then apply the same Euclidean plane projection.
        const double fallback_axis[3] = {
            AbsCUDA(plane_normal[0]) < 0.9 ? 1.0 : 0.0,
            AbsCUDA(plane_normal[1]) < 0.9 ? 1.0 : 0.0,
            AbsCUDA(plane_normal[2]) < 0.9 ? 1.0 : 0.0
        };
        const double fallback_dot =
            fallback_axis[0] * plane_normal[0] +
            fallback_axis[1] * plane_normal[1] +
            fallback_axis[2] * plane_normal[2];
        plane_up[0] = fallback_axis[0] - fallback_dot * plane_normal[0];
        plane_up[1] = fallback_axis[1] - fallback_dot * plane_normal[1];
        plane_up[2] = fallback_axis[2] - fallback_dot * plane_normal[2];
        up_sq =
            plane_up[0] * plane_up[0] +
            plane_up[1] * plane_up[1] +
            plane_up[2] * plane_up[2];
    }
    if (!isfinite(up_sq) || up_sq <= 1.0e-18) {
        return false;
    }
    const double inverse_up_norm = 1.0 / SqrtCUDA(up_sq);
    plane_up[0] *= inverse_up_norm;
    plane_up[1] *= inverse_up_norm;
    plane_up[2] *= inverse_up_norm;

    // Right direction completes the right-handed local plane frame:
    // right = up x normal, and the stored coordinates are (right, up).
    const double plane_right[3] = {
        plane_up[1] * plane_normal[2] - plane_up[2] * plane_normal[1],
        plane_up[2] * plane_normal[0] - plane_up[0] * plane_normal[2],
        plane_up[0] * plane_normal[1] - plane_up[1] * plane_normal[0]
    };

    //==========================================
    // LOCAL COORDINATES AND NONTERMINAL RECORD
    //==========================================
    const double relative_pos[3] = {
        x_intersect - d_commondata.non_terminal_plane_center_x,
        y_intersect - d_commondata.non_terminal_plane_center_y,
        z_intersect - d_commondata.non_terminal_plane_center_z
    };
    final_blueprint_data->y_nt =
        relative_pos[0] * plane_right[0] +
        relative_pos[1] * plane_right[1] +
        relative_pos[2] * plane_right[2];
    final_blueprint_data->z_nt =
        relative_pos[0] * plane_up[0] +
        relative_pos[1] * plane_up[1] +
        relative_pos[2] * plane_up[2];
    final_blueprint_data->non_terminal_plane_t = t_intersect;
    final_blueprint_data->non_terminal_plane_lambda = physical_lambda;

    // A nonterminal crossing is diagnostic only; it never changes status.
    return true;
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
