"""
Generates the C helper function for the RKF45 adaptive step-size controller.

Author: Dalton J. Moone
"""

from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h


def rkf45_update_and_control_helper() -> None:
    """
    Generate and register the RKF45 controller for injection into BHaH_defines.h.

    This function generates the C code for the `static inline` helper function
    `update_photon_state_and_stepsize()`. This function implements a robust,
    GSL-style error controller that decides whether to accept or reject a step
    and computes the optimal size for the next step.
    """
    # The entire C code block to be injected into the header.
    # This is algorithmic, not symbolic.
    c_code_for_header = r"""
// --- Adaptive Step-Size Controller (Robust, GSL-Style) ---
// Manages the adaptive step-size control for the RKF45 stepper.
static inline bool update_photon_state_and_stepsize(
    PhotonState *restrict photon,         // Pointer to the photon's full state
    const double y_start[9],              // The state at the START of the step (for scaling)
    const double y_out[9],                // The 5th-order result from the kernel
    const double y_err[9],                // The error vector from the kernel
    const commondata_struct *restrict commondata // For accessing control parameters
) {
    // This function implements a robust error control mechanism similar to that
    // used in the GNU Scientific Library (GSL). It computes a scale factor for
    // each component of the state vector and calculates a weighted error norm.

    const double rtol = commondata->rkf45_error_tolerance;
    const double atol = commondata->rkf45_absolute_error_tolerance;
    double err_norm_sq = 0.0;

    // --- Step 1: Calculate the squared error norm, treating components differently ---

    // For physical components (spatial position and 4-momentum, y[1]..y[7]),
    // use a scale that combines absolute and relative tolerances.
    for (int i = 1; i < 8; ++i) {
        const double scale_y = atol + rtol * fabs(y_start[i]);
        const double ratio = y_err[i] / scale_y;
        err_norm_sq += ratio * ratio;
    }

    // For coordinate time (y[0]) and path length (y[8]), which can grow
    // indefinitely, use a purely absolute tolerance to prevent their large
    // magnitudes from dominating the error and making the controller insensitive
    // to physical errors.
    const double scale_t = atol;
    const double ratio_t = y_err[0] / scale_t;
    err_norm_sq += ratio_t * ratio_t;

    const double scale_L = atol;
    const double ratio_L = y_err[8] / scale_L;
    err_norm_sq += ratio_L * ratio_L;

    // Final error norm is the root-mean-square of the weighted errors.
    const double err_norm = sqrt(err_norm_sq / 9.0);

    // --- Step 2: Accept or reject the step ---
    const double tolerance = 1.0; // The target for our normalized error is 1.0
    bool step_accepted = (err_norm <= tolerance);

    // --- Step 3: Calculate the optimal size for the next step ---
    double h_new;
    if (err_norm > 1e-15) {
        // Standard formula for step-size adjustment.
        // The exponent 0.2 is 1/5, appropriate for a 4(5) method.
        h_new = commondata->rkf45_safety_factor * photon->h * pow(tolerance / err_norm, 0.2);
    } else {
        // If error is zero or tiny, increase step size by a fixed, safe factor.
        h_new = 2.0 * photon->h;
    }

    // Enforce minimum and maximum step sizes.
    h_new = fmax(h_new, commondata->rkf45_h_min);
    h_new = fmin(h_new, commondata->rkf45_h_max);

    // --- Step 4: Update the photon's state based on the decision ---
    if (step_accepted) {
        // If step is accepted, update the state and reset retry counter.
        for (int i = 0; i < 9; ++i) {
            photon->y[i] = y_out[i];
        }
        photon->affine_param += photon->h;
        photon->rejection_retries = 0;
    } else {
        // If step is rejected, increment retry counter. The state is NOT updated.
        photon->rejection_retries++;
    }

    // The step size for the *next* attempt is always updated.
    photon->h = h_new;
    return step_accepted;
}
"""

    # Register this C code block to be injected into the BHaH_defines.h header file.
    Bdefines_h.register_BHaH_defines("rkf45_update_control", c_code_for_header)