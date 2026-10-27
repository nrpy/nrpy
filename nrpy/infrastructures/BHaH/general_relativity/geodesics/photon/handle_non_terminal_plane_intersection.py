# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/handle_non_terminal_plane_intersection.py
"""
Defines the C engine that handles a nonterminal-plane intersection.

The generated C function projects a refined photon state onto the nonterminal
plane. It uses that plane's center, normal, and up seed to compute local
coordinates. The crossing result leaves the photon trajectory unchanged.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par


def register_non_terminal_plane_crossing_result() -> None:
    """Register the state and measured quantities at a nonterminal crossing."""
    if "non_terminal_plane_crossing_result" in par.glb_extras_dict.get(
        "BHaH_defines", {}
    ):
        return
    Bdefines_h.register_BHaH_defines(
        "non_terminal_plane_crossing_result",
        r"""
    typedef struct {
        bool found;
        double state[9];
        double integration_param;
        double coordinate_time;
        double physical_lambda;
        double y_local;
        double z_local;
        double plane_distance;
        int brent_iterations;
        bool null_constraint_recorded;
        double null_constraint_abs;
        bool Lz_recorded;
        double Lz;
    } non_terminal_plane_crossing_result_t;
""",
    )


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


def handle_non_terminal_plane_intersection(normalized_eom: bool = False) -> None:
    """
    Generate local nonterminal-plane quantities from a refined root.

    :param normalized_eom: Whether coordinate time is the integration parameter.
    """
    parallelization = par.parval_from_str("parallelization")

    register_non_terminal_plane_crossing_result()
    register_non_terminal_plane_parameters()

    # Add the access variable
    cd_access = parallel_utils.get_commondata_access(parallelization)

    # Define inclusion headers for the C function compilation.
    includes = ["BHaH_defines.h"]

    if parallelization == "cuda":
        includes.append("BHaH_device_defines.h")
        includes.append("cuda_intrinsics.h")

    # Define the Doxygen-formatted description for the C function.
    desc = "Project a refined photon state onto the nonterminal plane."

    # Specify the signature components of the C function.
    cfunc_type = "__device__ bool" if parallelization == "cuda" else "bool"
    name = "handle_non_terminal_plane_intersection"

    # Define the specific C arguments passed into the kernel.
    params = (
        "const photon_plane_crossing_root_t *restrict root, "
        "non_terminal_plane_crossing_result_t *restrict result"
    )
    if parallelization != "cuda":
        params += ", const commondata_struct *restrict commondata"

    # Toggle generation of CodeParameters.h inclusion.
    include_CodeParameters_h = False

    # Construct the C-string for the function body.
    body = r"""
    *result = (non_terminal_plane_crossing_result_t){0};
    if (!root->found) return false;
    const double *f_local = root->state;
    //==========================================
    // STATE UNPACKING
    //==========================================
    const double t_intersect = @COORDINATE_TIME@;
    const double physical_lambda = @PHYSICAL_LAMBDA@;
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
    } // END IF: invalid nonterminal plane normal
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

    const double up_sq =
        plane_up[0] * plane_up[0] +
        plane_up[1] * plane_up[1] +
        plane_up[2] * plane_up[2];
    if (!isfinite(up_sq) || up_sq <= 1.0e-18) {
        return false;
    } // END IF: invalid projected up vector
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
    // LOCAL COORDINATES AND NONTERMINAL RESULT
    //==========================================
    const double relative_pos[3] = {
        x_intersect - d_commondata.non_terminal_plane_center_x,
        y_intersect - d_commondata.non_terminal_plane_center_y,
        z_intersect - d_commondata.non_terminal_plane_center_z
    };
    result->y_local =
        relative_pos[0] * plane_right[0] +
        relative_pos[1] * plane_right[1] +
        relative_pos[2] * plane_right[2];
    result->z_local =
        relative_pos[0] * plane_up[0] +
        relative_pos[1] * plane_up[1] +
        relative_pos[2] * plane_up[2];
    if (!isfinite(result->y_local) || !isfinite(result->z_local) ||
        !isfinite(t_intersect) || !isfinite(physical_lambda)) return false;
    result->found = true;
    result->integration_param = root->parameter;
    result->coordinate_time = t_intersect;
    result->physical_lambda = physical_lambda;
    result->plane_distance = root->plane_distance;
    result->brent_iterations = root->iterations;
    result->null_constraint_abs = NAN;
    result->Lz = NAN;
    for (int component = 0; component < 9; ++component) {
        result->state[component] = f_local[component];
    }
    return true;
    """

    body = body.replace(
        "@COORDINATE_TIME@", "root->parameter" if normalized_eom else "f_local[0]"
    ).replace(
        "@PHYSICAL_LAMBDA@", "f_local[0]" if normalized_eom else "root->parameter"
    )
    host_body = body.replace("d_commondata.", "commondata->").replace(
        "SqrtCUDA(", "sqrt("
    )
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
    if parallelization == "cuda":
        cfc.register_CFunction(
            includes=["BHaH_defines.h", "<math.h>"],
            desc="Project a refined photon state onto the nonterminal plane on the host.",
            cfunc_type="bool",
            name="handle_non_terminal_plane_intersection_host",
            params=(
                "const photon_plane_crossing_root_t *restrict root, "
                "non_terminal_plane_crossing_result_t *restrict result, "
                "const commondata_struct *restrict commondata"
            ),
            include_CodeParameters_h=False,
            body=host_body,
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
