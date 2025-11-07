"""
Generates the C helper functions for the RKF45 computational kernels.

Author: Dalton J. Moone
"""

from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h


def rkf45_helpers_for_header() -> None:
    """
    Generate and register the RKF45 kernels for injection into BHaH_defines.h.

    This function generates the C code for two `static inline` helper functions:
    1. `calculate_intermediate_state()`: Computes the state vector at one of the
       six intermediate stages of an RKF45 step.
    2. `rkf45_kernel()`: Performs the final summations to compute the 4th and
       5th-order solutions and the error estimate.
    """
    # The entire C code block to be injected into the header.
    # This is algorithmic, not symbolic.
    c_code_for_header = r"""
// =============================================
// NRPy-Generated RKF45 Stepper Helpers
// =============================================

// --- Intermediate State Calculator ---
// Computes the intermediate state vector for a specific RKF45 stage.
static inline void calculate_intermediate_state(
    const int stage,                // The stage to compute (1-6, for k1-k6)
    const double y_in[9],           // The state at the beginning of the step
    const double k_array[6][9],     // Array of the k vectors (only previous stages are valid)
    const double h,                 // The step size
    double y_temp[9]                // Output: the temporary state vector for the stage
) {
    switch (stage) {
        case 1: // k1 is evaluated at y_in
            for (int i = 0; i < 9; ++i) y_temp[i] = y_in[i];
            break;
        case 2: // k2
            for (int i = 0; i < 9; ++i) y_temp[i] = y_in[i] + h * (1.0/4.0) * k_array[0][i];
            break;
        case 3: // k3
            for (int i = 0; i < 9; ++i) y_temp[i] = y_in[i] + h * ( (3.0/32.0)*k_array[0][i] + (9.0/32.0)*k_array[1][i] );
            break;
        case 4: // k4
            for (int i = 0; i < 9; ++i) y_temp[i] = y_in[i] + h * ( (1932.0/2197.0)*k_array[0][i] - (7200.0/2197.0)*k_array[1][i] + (7296.0/2197.0)*k_array[2][i] );
            break;
        case 5: // k5
            for (int i = 0; i < 9; ++i) y_temp[i] = y_in[i] + h * ( (439.0/216.0)*k_array[0][i] - 8.0*k_array[1][i] + (3680.0/513.0)*k_array[2][i] - (845.0/4104.0)*k_array[3][i] );
            break;
        case 6: // k6
            for (int i = 0; i < 9; ++i) y_temp[i] = y_in[i] + h * ( -(8.0/27.0)*k_array[0][i] + 2.0*k_array[1][i] - (3544.0/2565.0)*k_array[2][i] + (1859.0/4104.0)*k_array[3][i] - (11.0/40.0)*k_array[4][i] );
            break;
    }
}

// --- RKF45 Kernel ---
// Pure computational kernel for the RKF45 method.
static inline void rkf45_kernel(
    const double y_in[9],           // The state at the beginning of the step
    const double k_array[6][9],     // Array of the 6 pre-computed k vectors
    const double h,                 // The step size attempted
    double y_out[9],                // Output: the final 5th-order state
    double y_err[9]                 // Output: the error vector (y_5th - y_4th)
) {
    // Calculate the 4th-Order Accurate Solution for the error estimate.
    double y_4th[9];
    for (int i = 0; i < 9; ++i) {
        y_4th[i] = y_in[i] + h * ( (25.0/216.0) * k_array[0][i] +
                                   (1408.0/2565.0) * k_array[2][i] +
                                   (2197.0/4104.0) * k_array[3][i] -
                                   (1.0/5.0) * k_array[4][i] );
    }

    // Calculate the 5th-Order Accurate Solution for the final state.
    for (int i = 0; i < 9; ++i) {
        y_out[i] = y_in[i] + h * ( (16.0/135.0) * k_array[0][i] +
                                   (6656.0/12825.0) * k_array[2][i] +
                                   (28561.0/56430.0) * k_array[3][i] -
                                   (9.0/50.0) * k_array[4][i] +
                                   (2.0/55.0) * k_array[5][i] );
    }

    // Calculate the Error Vector.
    for (int i = 0; i < 9; ++i) {
        y_err[i] = y_out[i] - y_4th[i];
    }
}
"""

    # Register this C code block to be injected into the BHaH_defines.h header file.
    Bdefines_h.register_BHaH_defines("rkf45_helpers", c_code_for_header)