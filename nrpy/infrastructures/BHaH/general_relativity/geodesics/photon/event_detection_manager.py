"""
Generates the C orchestrator for geometric event detection.

This module provides the high-level logic for detecting crossings of the
observer window and the source emission plane. It evaluates the scalar product
between the photon trajectory vector and the target plane normal, triggering
high-precision root-finding upon a sign transition. The logic executes entirely
on thread-local arrays to map optimally to the streaming bundle architecture.
Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h

# Define the C-code for event types and geometry parameters mapped to pure C structs
event_mgmt_c_code = r"""
    /**
     * @brief Defines the geometric parameters for a plane in 3D space.
     * The plane is defined by the equation: $n_x x + n_y y + n_z z - d = 0$.
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

    :raises SystemError: If C function registration fails within the NRPy+ pipeline.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>", "<stdbool.h>"]
    desc = """@brief GPU-optimized detection of plane crossings using thread-local data.

    @param f_local The thread-local 9-component array for the current state $f^\mu_{n}$.
    @param f_p_local The thread-local state array for step $f^\mu_{n-1}$.
    @param f_p_p_local The thread-local state array for step $f^\mu_{n-2}$.
    @param window_event_found Pointer to the local boolean flag for window intersection.
    @param on_pos_window_prev Pointer to the local boolean flag tracking previous window side.
    @param window_event_f_intersect Thread-local array storing the interpolated window state $f^\mu_{w}$.
    @param window_event_lambda Pointer to the local affine parameter $\lambda_w$ for the window hit.
    @param source_event_found Pointer to the local boolean flag for source intersection.
    @param on_pos_source_prev Pointer to the local boolean flag tracking previous source side.
    @param source_event_f_intersect Thread-local array storing the interpolated source state $f^\mu_{s}$.
    @param source_event_lambda Pointer to the local affine parameter $\lambda_s$ for the source hit.
    @param commondata Pointer to the globally constant parameters struct.

    Detailed algorithm: Evaluates the scalar product of the photon's position $x^i$ and the
    plane normal $n_i$. A transition in the sign of this product indicates a crossing,
    triggering the high-precision root-finder $f(x) = 0$. Mapping logic strictly to `f_local` 
    preserves the register limits for the sm_86 architecture by ensuring all calculations remain 
    within the 255-register per thread constraint, entirely circumventing global memory reads."""
    cfunc_type = "BHAH_HD_FUNC void"
    name = "event_detection_manager"
    params = (
        "const double *restrict f_local, "
        "const double *restrict f_p_local, "
        "const double *restrict f_p_p_local, "
        "bool *restrict window_event_found, "
        "bool *restrict on_pos_window_prev, "
        "double *restrict window_event_f_intersect, "
        "double *restrict window_event_lambda, "
        "bool *restrict source_event_found, "
        "bool *restrict on_pos_source_prev, "
        "double *restrict source_event_f_intersect, "
        "double *restrict source_event_lambda, "
        "const commondata_struct *restrict commondata"
    )
    include_CodeParameters_h = False
    body = r"""
    // --- EVENT DETECTION & TERMINATION CHECKS ---
    // Evaluates terminal bounding geometries to trigger root-finding logic.
    // Extracts spatial coordinates directly from the thread-local state array to guarantee
    // rapid execution on native GPU registers without global VRAM penalties.

    const double x = f_local[1]; // Cartesian $x$-coordinate of the photon $x^1$.
    const double y = f_local[2]; // Cartesian $y$-coordinate of the photon $x^2$.
    const double z = f_local[3]; // Cartesian $z$-coordinate of the photon $x^3$.

    // --- Observer Window Plane Logic ---
    if (!(*window_event_found)) {
        // Calculate the window plane normal vector $n_i$ based on camera and window centers.
        double w_normal[3];
        w_normal[0] = commondata->window_center_x - commondata->camera_pos_x; // Normal component along $x$-axis.
        w_normal[1] = commondata->window_center_y - commondata->camera_pos_y; // Normal component along $y$-axis.
        w_normal[2] = commondata->window_center_z - commondata->camera_pos_z; // Normal component along $z$-axis.

        // Normalize the vector $n_i n^i = 1$ to ensure the plane evaluation is purely geometric.
        const double mag_inv = 1.0 / sqrt(w_normal[0]*w_normal[0] + w_normal[1]*w_normal[1] + w_normal[2]*w_normal[2]);
        for(int i=0; i<3; i++) w_normal[i] *= mag_inv;

        // Calculate the scalar distance $d$ from the origin to the window plane.
        const double w_dist = commondata->window_center_x*w_normal[0] + commondata->window_center_y*w_normal[1] + commondata->window_center_z*w_normal[2];

        // Evaluate the plane equation $n_i x^i - d = 0$ at the current photon position.
        const double w_val = x*w_normal[0] + y*w_normal[1] + z*w_normal[2] - w_dist;
        const bool on_pos_curr = (w_val > 1e-10); // True if the photon is on the positive side of the normal.

        // Zero-Crossing Detection: Compare current state with the previous integration step.
        if (on_pos_curr != *on_pos_window_prev) {
             find_event_time_and_state(f_local, f_p_local, f_p_p_local, w_normal, w_dist, window_event_lambda, window_event_f_intersect);
             *window_event_found = true;
        }
        *on_pos_window_prev = on_pos_curr;
    }

    // --- Physical Source Plane Logic ---
    if (!(*source_event_found)) {
        // Unpack source plane normal $n_i$ from read-only common constant memory.
        const double s_normal[3] = {commondata->source_plane_normal_x, commondata->source_plane_normal_y, commondata->source_plane_normal_z};

        // Calculate the scalar distance $d$ to the physical emission source plane.
        const double s_dist = commondata->source_plane_center_x*s_normal[0] + commondata->source_plane_center_y*s_normal[1] + commondata->source_plane_center_z*s_normal[2];

        // Evaluate the plane equation $n_i x^i - d = 0$ at the current photon position.
        const double s_val = x*s_normal[0] + y*s_normal[1] + z*s_normal[2] - s_dist;
        const bool on_pos_curr = (s_val > 1e-10); // True if the photon is on the positive side.

        // Zero-Crossing Detection: Triggers the root-finder upon crossing the physical emission disk.
        if (on_pos_curr != *on_pos_source_prev) {
             find_event_time_and_state(f_local, f_p_local, f_p_p_local, s_normal, s_dist, source_event_lambda, source_event_f_intersect);
             *source_event_found = true;
        }
        *on_pos_source_prev = on_pos_curr;
    }
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )