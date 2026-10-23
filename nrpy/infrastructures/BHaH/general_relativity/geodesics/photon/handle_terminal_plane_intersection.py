# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/handle_terminal_plane_intersection.py
"""
Defines the computation for impact parameters on the physical terminal plane.

The generated C function projects a refined photon state onto the terminal
plane. It constructs the plane basis from its normal and up vector, computes
local coordinates, and records whether the crossing lies inside the radial
bounds. The wrapper selects one in-range crossing for each photon.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par


def register_terminal_plane_crossing_result() -> None:
    """Register the selected terminal crossing and its plane quantities."""
    if "terminal_plane_crossing_result" in par.glb_extras_dict.get("BHaH_defines", {}):
        return
    Bdefines_h.register_BHaH_defines(
        "terminal_plane_crossing_result",
        r"""
    typedef struct {
        bool found;
        bool inside_radial_bounds;
        int selected_candidate_number; // Zero means no terminal candidate selected.
        termination_type_t trajectory_stop_reason;
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
    } terminal_plane_crossing_result_t;
""",
    )


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


def handle_terminal_plane_intersection(normalized_eom: bool = False) -> None:
    """Register local terminal-plane coordinates and radial bounds."""
    parallelization = par.parval_from_str("parallelization")

    register_terminal_plane_crossing_result()
    register_terminal_plane_parameters()

    # Add the access variable
    cd_access = parallel_utils.get_commondata_access(parallelization)

    includes = ["BHaH_defines.h"]

    if parallelization == "cuda":
        includes.append("BHaH_device_defines.h")
        includes.append("cuda_intrinsics.h")

    desc = "Project a refined photon state onto the terminal plane and check radius."

    cfunc_type = "__device__ bool" if parallelization == "cuda" else "bool"

    name = "handle_terminal_plane_intersection"

    params = (
        "const photon_plane_crossing_root_t *restrict root, "
        "terminal_plane_crossing_result_t *restrict result"
    )
    if parallelization != "cuda":
        params += ", const commondata_struct *restrict commondata"

    include_CodeParameters_h = False

    body = r"""
    *result = (terminal_plane_crossing_result_t){0};
    if (!root->found) return false;
    const double *terminal_plane_event_f_intersect = root->state;
    //==========================================
    // UNPACK INTERSECTION STATE
    //==========================================
    // Root parameter stores coordinate time for normalized equations and
    // affine parameter for direct equations.
    const double t_intersect = @COORDINATE_TIME@;
    const double physical_lambda = @PHYSICAL_LAMBDA@;
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
    } // END IF: invalid terminal plane normal
    const double inverse_normal_norm = 1.0 / SqrtCUDA(normal_sq);
    terminal_plane_normal[0] *= inverse_normal_norm;
    terminal_plane_normal[1] *= inverse_normal_norm;
    terminal_plane_normal[2] *= inverse_normal_norm;

    // Project the supplied up vector perpendicular to the plane normal.
    double terminal_plane_up[3] = {
        d_commondata.terminal_plane_up_x,
        d_commondata.terminal_plane_up_y,
        d_commondata.terminal_plane_up_z
    };
    const double up_dot_normal =
        terminal_plane_up[0] * terminal_plane_normal[0] +
        terminal_plane_up[1] * terminal_plane_normal[1] +
        terminal_plane_up[2] * terminal_plane_normal[2];
    terminal_plane_up[0] -= up_dot_normal * terminal_plane_normal[0];
    terminal_plane_up[1] -= up_dot_normal * terminal_plane_normal[1];
    terminal_plane_up[2] -= up_dot_normal * terminal_plane_normal[2];

    const double up_sq =
        terminal_plane_up[0] * terminal_plane_up[0] +
        terminal_plane_up[1] * terminal_plane_up[1] +
        terminal_plane_up[2] * terminal_plane_up[2];
    if (!isfinite(up_sq) || up_sq <= 1.0e-18) {
        return false;
    } // END IF: invalid projected up vector
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
    // TERMINAL PLANE RADIAL BOUNDS
    //==========================================
    // Filters photons that fall outside the physical bounds of the accretion disk.
    const double r_sq = y_t*y_t + z_t*z_t; // Squared radial distance $r^2$ from the terminal-plane center.
    if (!isfinite(r_sq) || !isfinite(t_intersect) ||
        !isfinite(physical_lambda)) return false;
    const double r_min = d_commondata.terminal_plane_min_coord_radius;
    const double r_max = d_commondata.terminal_plane_max_coord_radius;
    if (!isfinite(r_min) || !isfinite(r_max) || r_min < 0.0 || r_max < r_min) {
        return false;
    } // END IF: invalid radius bounds
    const double r_min_sq = r_min * r_min; // Squared minimum radius.
    const double r_max_sq = r_max * r_max; // Squared maximum radius.

    result->found = true;
    result->inside_radial_bounds = r_sq >= r_min_sq && r_sq <= r_max_sq;
    result->integration_param = root->parameter;
    result->coordinate_time = t_intersect;
    result->physical_lambda = physical_lambda;
    result->plane_distance = root->plane_distance;
    result->brent_iterations = root->iterations;
    result->y_local = y_t;
    result->z_local = z_t;
    result->null_constraint_abs = NAN;
    result->Lz = NAN;
    for (int component = 0; component < 9; ++component) {
        result->state[component] = terminal_plane_event_f_intersect[component];
    }
    return true;
    """

    body = body.replace(
        "@COORDINATE_TIME@",
        "root->parameter" if normalized_eom else "terminal_plane_event_f_intersect[0]",
    ).replace(
        "@PHYSICAL_LAMBDA@",
        "terminal_plane_event_f_intersect[0]" if normalized_eom else "root->parameter",
    )
    host_body = body.replace("d_commondata.", "commondata->").replace(
        "SqrtCUDA(", "sqrt("
    )
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
    if parallelization == "cuda":
        cfc.register_CFunction(
            includes=["BHaH_defines.h", "<math.h>"],
            desc="Project a refined photon state onto the terminal plane on the host.",
            cfunc_type="bool",
            name="handle_terminal_plane_intersection_host",
            params=(
                "const photon_plane_crossing_root_t *restrict root, "
                "terminal_plane_crossing_result_t *restrict result, "
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
