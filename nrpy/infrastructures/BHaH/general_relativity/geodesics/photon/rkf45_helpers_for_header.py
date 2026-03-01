"""
Generates the C helper functions for the RKF45 computational kernels.

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module provides the core mathematical stages and error estimation logic
for the Runge-Kutta-Fehlberg 4(5) embedded integration scheme.

Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc


def rkf45_helpers_for_header() -> None:
    """
    Register the Runge-Kutta-Fehlberg 4(5) intermediate stage and kernel functions.

    Generates purely thread-local C functions responsible for evaluating the intermediate
    stages of the RKF45 algorithm and computing the final candidate state
    along with its local truncation error.

    :raises TypeError: If incorrect parameters are passed to the code generation functions.
    """
    # --- Function 1: calculate_rkf45_stage_f_temp ---
    includes_stage = ["BHaH_defines.h"]
    desc_stage = r"""@brief Calculates the intermediate state for a specific RKF45 stage.
    
    @param stage The integer RKF45 stage identifier (1 through 6).
    @param f_in_local Thread-local array containing the initial 9-component physical state $f^\mu$.
    @param k_array_local Thread-local array holding intermediate $k$ evaluation derivatives. Maps directly to registers to avoid sm_86 limits.
    @param h The local integration step size $h$.
    @param f_temp_out_local Thread-local array to store the computed intermediate stage state.

    Algorithm:
    1. Evaluates the intermediate Runge-Kutta stages using Fehlberg's optimal coefficients.
    2. Performs strictly thread-local scalar additions to compute candidate vectors."""
    cfunc_type_stage = "BHAH_HD_FUNC void"
    name_stage = "calculate_rkf45_stage_f_temp"
    params_stage = (
        "const int stage, "
        "const double *restrict f_in_local, "
        "const double *restrict k_array_local, "
        "const double h, "
        "double *restrict f_temp_out_local"
    )
    include_CodeParameters_h_stage = False
    body_stage = r"""
    // --- STAGE EVALUATION ROUTINE ---
    // The physical 9-component state vector $f^\mu$ is strictly mapped to thread-local variables.
    // 0: $t$, 1-3: $x,y,z$, 4: $p_t$, 5-7: $p_x,p_y,p_z$, 8: $\lambda$
    switch (stage) {
        case 1:
            for (int i = 0; i < 9; ++i) {
                f_temp_out_local[i] = f_in_local[i]; // Stage 1 state transfer.
            }
            break;
        case 2:
            for (int i = 0; i < 9; ++i) {
                f_temp_out_local[i] = f_in_local[i] +
                    h * (1.0/4.0) * k_array_local[0*9 + i]; // Stage 2 evaluation step.
            }
            break;
        case 3:
            for (int i = 0; i < 9; ++i) {
                f_temp_out_local[i] = f_in_local[i] +
                    h * ( (3.0/32.0)*k_array_local[0*9 + i] +
                          (9.0/32.0)*k_array_local[1*9 + i] ); // Stage 3 evaluation step.
            }
            break;
        case 4:
            for (int i = 0; i < 9; ++i) {
                f_temp_out_local[i] = f_in_local[i] +
                    h * ( (1932.0/2197.0)*k_array_local[0*9 + i] -
                          (7200.0/2197.0)*k_array_local[1*9 + i] +
                          (7296.0/2197.0)*k_array_local[2*9 + i] ); // Stage 4 evaluation step.
            }
            break;
        case 5:
            for (int i = 0; i < 9; ++i) {
                f_temp_out_local[i] = f_in_local[i] +
                    h * ( (439.0/216.0)*k_array_local[0*9 + i] -
                          8.0*k_array_local[1*9 + i] +
                          (3680.0/513.0)*k_array_local[2*9 + i] -
                          (845.0/4104.0)*k_array_local[3*9 + i] ); // Stage 5 evaluation step.
            }
            break;
        case 6:
            for (int i = 0; i < 9; ++i) {
                f_temp_out_local[i] = f_in_local[i] +
                    h * ( -(8.0/27.0)*k_array_local[0*9 + i] +
                          2.0*k_array_local[1*9 + i] -
                          (3544.0/2565.0)*k_array_local[2*9 + i] +
                          (1859.0/4104.0)*k_array_local[3*9 + i] -
                          (11.0/40.0)*k_array_local[4*9 + i] ); // Stage 6 evaluation step.
            }
            break;
    }
    """

    cfc.register_CFunction(
        includes=includes_stage,
        desc=desc_stage,
        cfunc_type=cfunc_type_stage,
        name=name_stage,
        params=params_stage,
        include_CodeParameters_h=include_CodeParameters_h_stage,
        body=body_stage,
    )

    # --- Function 2: rkf45_kernel ---
    includes_kernel = ["BHaH_defines.h"]
    desc_kernel = r"""@brief Finalizes the RKF45 candidate states and computes the local truncation error.
    
    @param f_in_local Thread-local array containing the initial physical state $f^\mu$.
    @param k_array_local Thread-local array holding intermediate $k$ evaluation derivatives.
    @param h The local integration step size $h$.
    @param f_out_local Thread-local array to store the newly proposed 5th-order physical state.
    @param f_err_local Thread-local array computing the calculated local truncation error.

    Algorithm:
    1. Evaluates the 4th-order candidate state strictly for error baseline approximation.
    2. Evaluates the 5th-order candidate state representing the physical step progression.
    3. Computes the truncation error disparity to feed the step-size controller."""
    cfunc_type_kernel = "BHAH_HD_FUNC void"
    name_kernel = "rkf45_kernel"
    params_kernel = (
        "const double *restrict f_in_local, "
        "const double *restrict k_array_local, "
        "const double h, "
        "double *restrict f_out_local, "
        "double *restrict f_err_local"
    )
    include_CodeParameters_h_kernel = False
    body_kernel = r"""
    // --- RKF45 ERROR ESTIMATION & STATE UPDATE ---
    for (int i = 0; i < 9; ++i) {
        // The 4th-order candidate state, used exclusively as the baseline for error estimation.
        const double f_4th = f_in_local[i] + h * (
                                       (25.0/216.0) * k_array_local[0*9 + i] +
                                       (1408.0/2565.0) * k_array_local[2*9 + i] +
                                       (2197.0/4104.0) * k_array_local[3*9 + i] -
                                       (1.0/5.0) * k_array_local[4*9 + i] );

        // The 5th-order candidate state, representing the newly proposed physical configuration.
        f_out_local[i] = f_in_local[i] + h * (
                                   (16.0/135.0) * k_array_local[0*9 + i] +
                                   (6656.0/12825.0) * k_array_local[2*9 + i] +
                                   (28561.0/56430.0) * k_array_local[3*9 + i] -
                                   (9.0/50.0) * k_array_local[4*9 + i] +
                                   (2.0/55.0) * k_array_local[5*9 + i] );

        // The calculated local truncation error derived from the embedded method's order disparity.
        f_err_local[i] = f_out_local[i] - f_4th;
    }
    """

    cfc.register_CFunction(
        includes=includes_kernel,
        desc=desc_kernel,
        cfunc_type=cfunc_type_kernel,
        name=name_kernel,
        params=params_kernel,
        include_CodeParameters_h=include_CodeParameters_h_kernel,
        body=body_kernel,
    )