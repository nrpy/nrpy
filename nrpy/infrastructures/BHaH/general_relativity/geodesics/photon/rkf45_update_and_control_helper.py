"""
Generates the C adaptive step-size controller, updated for BHaH integration.

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
Features a high-performance L_infinity norm and L1 momentum floor for
relativistic ray tracing. Adapted for pure thread-local GPU execution.

Author: Dalton J. Moone.
"""

import nrpy.params as par
import nrpy.c_function as cfc


def rkf45_update_and_control_helper() -> None:
    """
    Register the high-performance adaptive step-size controller for RKF45.

    Generates the C function `update_photon_state_and_stepsize`, handling
    local truncation error calculation, step acceptance logic, and time
    step adaptation based on L-infinity and L1 bounds exclusively within
    thread-local registers.

    :raises TypeError: If incorrect parameters are passed to the code generation functions.
    """
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "rkf45_error_tolerance",
            "rkf45_absolute_error_tolerance",
            "rkf45_h_min",
            "rkf45_h_max",
            "rkf45_safety_factor",
            "numerical_initial_h",
        ],
        [1e-8, 1e-8, 1e-10, 10.0, 0.9, 1.0],
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameter(
        "int", __name__, "rkf45_max_retries", 10, commondata=True, add_to_parfile=True
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>", "<stdbool.h>"]
    desc = r"""@brief Evaluates step validity and computes the adaptive time step $h$.
    
    @param f_local Thread-local array for the photon's primary state vector $f^\mu$.
    @param f_start_local Thread-local array capturing the state at the beginning of the RKF45 step.
    @param f_out_local Thread-local array capturing the proposed 5th-order physical state.
    @param f_err_local Thread-local array capturing the calculated truncation error.
    @param h_local Pointer to the thread-local integration step size $h$.
    @param affine_param_local Pointer to the thread-local affine parameter $\lambda$.
    @param retries_local Pointer to the thread-local variable tracking rejection count.
    @param commondata Global read-only struct containing integration tolerance parameters.

    Algorithm:
    1. Constructs composite scaling bounds utilizing an absolute tolerance floor and mixed norms.
    2. Identifies the maximum truncation error across all physical components.
    3. Updates local register states or increments retry counters without touching global VRAM."""
    cfunc_type = "BHAH_HD_FUNC bool"
    name = "update_photon_state_and_stepsize"
    params = (
        "double *restrict f_local, "
        "const double *restrict f_start_local, "
        "const double *restrict f_out_local, "
        "const double *restrict f_err_local, "
        "double *restrict h_local, "
        "double *restrict affine_param_local, "
        "int *restrict retries_local, "
        "const commondata_struct *restrict commondata"
    )
    include_CodeParameters_h = False
    body = r"""
    // --- TOLERANCE Bounding & ERROR CALCULATION ---
    const double rtol = commondata->rkf45_error_tolerance; // The relative error tolerance threshold.
    const double atol = commondata->rkf45_absolute_error_tolerance; // The absolute error tolerance floor, preventing division by zero.

    double err_norm = 0.0; // The L_infinity (maximum) error norm acting as the strict acceptance threshold.
    double ratio; // The local ratio of truncation error to allowed tolerance for a specific component.

    // 1. Coordinate Time ($f[0]$): Pure absolute tolerance to prevent secular drift
    ratio = fabs(f_err_local[0]) / atol;
    err_norm = fmax(err_norm, ratio);

    // 2. Spatial Coordinates ($f[1]$, $f[2]$, $f[3]$): Standard mixed norm
    for (int i = 1; i <= 3; ++i) {
        const double scale = atol + rtol * fabs(f_start_local[i]); // The composite scaling factor bounding spatial position variation.
        ratio = fabs(f_err_local[i]) / scale;
        err_norm = fmax(err_norm, ratio);
    }

    // 3. Temporal Momentum $p^t$ ($f[4]$): Standard mixed norm
    // Represents energy evolution; does not suffer from zero-crossings
    const double scale_pt = atol + rtol * fabs(f_start_local[4]); // The composite scaling factor bounding relativistic energy variation.
    ratio = fabs(f_err_local[4]) / scale_pt;
    err_norm = fmax(err_norm, ratio);

    // 4. Spatial Momenta $p^x$, $p^y$, $p^z$ ($f[5]$, $f[6]$, $f[7]$): L1 floored mixed norm
    const double p_L1 = fabs(f_start_local[5]) +
                        fabs(f_start_local[6]) +
                        fabs(f_start_local[7]); // The L1 Manhattan distance, serving as a robust physical floor for momentum scaling.

    for (int i = 5; i <= 7; ++i) {
        const double scale = atol + rtol * p_L1; // The composite scaling factor bounded by the L1 momentum floor.
        ratio = fabs(f_err_local[i]) / scale;
        err_norm = fmax(err_norm, ratio);
    }

    // NOTE: Path length $f[8]$ is explicitly excluded from the error controller.
    const bool step_accepted = (err_norm <= 1.0); // Boolean indicating if the worst-case error falls within the mandated threshold.

    // --- STEP ADAPTATION ---
    double h_new; // The newly calculated adaptive step size proposed for the subsequent integration step.
    if (err_norm > 1e-15) {
        h_new = commondata->rkf45_safety_factor * (*h_local) * pow(1.0 / err_norm, 0.2);
    } else {
        h_new = 2.0 * (*h_local);
    }

    h_new = fmax(h_new, commondata->rkf45_h_min);
    h_new = fmin(h_new, commondata->rkf45_h_max);

    // --- LOCAL REGISTER COMMITMENT ---
    // Safely updates the photon trajectory exclusively within thread-local arrays to bypass global VRAM bottlenecks.
    if (step_accepted) {
        for (int i = 0; i < 9; ++i) {
            f_local[i] = f_out_local[i]; // Update the core state vector.
        }
        *affine_param_local += *h_local; // Advance local affine parameter.
        *retries_local = 0; // Reset rejection retries count.
    } else {
        *retries_local += 1; // Increment rejection retries count.
    }

    *h_local = h_new; // Store the new step size $h$ into thread-local scope.
    return step_accepted;
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