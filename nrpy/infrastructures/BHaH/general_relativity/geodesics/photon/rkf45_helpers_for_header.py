"""
Generates the C helper functions for the RKF45 computational kernels.

Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
This module provides the core mathematical stages and error estimation logic 
for the Runge-Kutta-Fehlberg 4(5) embedded integration scheme.
"""

from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h


def rkf45_helpers_for_header(spacetime_name: str) -> None:
    """
    Registers the Runge-Kutta-Fehlberg 4(5) intermediate stage and kernel functions.

    Generates inline C functions responsible for evaluating the intermediate 
    stages of the RKF45 algorithm and computing the final candidate state 
    along with its local truncation error.

    Args:
        spacetime_name (str): The specific metric/spacetime identifier.

    >>> rkf45_helpers_for_header("Kerr")
    """
    c_code_for_header = r"""
// --- RKF45 Stage State Calculator ---
static inline void calculate_rkf45_stage_f_temp(
    const int stage, const double *restrict f_in, const double *restrict k_array, 
    const double h, double *restrict f_temp_out, const int batch_size, const int batch_id
) {
    // The physical 9-component state vector f is strictly mapped as:
    // 0: t, 1-3: x,y,z, 4: p_t, 5-7: p_x,p_y,p_z, 8: lambda
    switch (stage) {
        case 1:
            #pragma unroll
            for (int i = 0; i < 9; ++i) 
                f_temp_out[IDX_LOCAL(i, batch_id, batch_size)] = f_in[IDX_LOCAL(i, batch_id, batch_size)];
            break;
        case 2:
            #pragma unroll
            for (int i = 0; i < 9; ++i)
                f_temp_out[IDX_LOCAL(i, batch_id, batch_size)] = f_in[IDX_LOCAL(i, batch_id, batch_size)] + 
                    h * (1.0/4.0) * k_array[IDX_LOCAL(0*9 + i, batch_id, batch_size)];
            break;
        case 3:
            #pragma unroll
            for (int i = 0; i < 9; ++i)
                f_temp_out[IDX_LOCAL(i, batch_id, batch_size)] = f_in[IDX_LOCAL(i, batch_id, batch_size)] + 
                    h * ( (3.0/32.0)*k_array[IDX_LOCAL(0*9 + i, batch_id, batch_size)] + 
                          (9.0/32.0)*k_array[IDX_LOCAL(1*9 + i, batch_id, batch_size)] );
            break;
        case 4:
            #pragma unroll
            for (int i = 0; i < 9; ++i)
                f_temp_out[IDX_LOCAL(i, batch_id, batch_size)] = f_in[IDX_LOCAL(i, batch_id, batch_size)] + 
                    h * ( (1932.0/2197.0)*k_array[IDX_LOCAL(0*9 + i, batch_id, batch_size)] - 
                          (7200.0/2197.0)*k_array[IDX_LOCAL(1*9 + i, batch_id, batch_size)] + 
                          (7296.0/2197.0)*k_array[IDX_LOCAL(2*9 + i, batch_id, batch_size)] );
            break;
        case 5:
            #pragma unroll
            for (int i = 0; i < 9; ++i)
                f_temp_out[IDX_LOCAL(i, batch_id, batch_size)] = f_in[IDX_LOCAL(i, batch_id, batch_size)] + 
                    h * ( (439.0/216.0)*k_array[IDX_LOCAL(0*9 + i, batch_id, batch_size)] - 
                          8.0*k_array[IDX_LOCAL(1*9 + i, batch_id, batch_size)] + 
                          (3680.0/513.0)*k_array[IDX_LOCAL(2*9 + i, batch_id, batch_size)] - 
                          (845.0/4104.0)*k_array[IDX_LOCAL(3*9 + i, batch_id, batch_size)] );
            break;
        case 6:
            #pragma unroll
            for (int i = 0; i < 9; ++i)
                f_temp_out[IDX_LOCAL(i, batch_id, batch_size)] = f_in[IDX_LOCAL(i, batch_id, batch_size)] + 
                    h * ( -(8.0/27.0)*k_array[IDX_LOCAL(0*9 + i, batch_id, batch_size)] + 
                          2.0*k_array[IDX_LOCAL(1*9 + i, batch_id, batch_size)] - 
                          (3544.0/2565.0)*k_array[IDX_LOCAL(2*9 + i, batch_id, batch_size)] + 
                          (1859.0/4104.0)*k_array[IDX_LOCAL(3*9 + i, batch_id, batch_size)] - 
                          (11.0/40.0)*k_array[IDX_LOCAL(4*9 + i, batch_id, batch_size)] );
            break;
    }
}

// --- RKF45 Kernel ---
static inline void rkf45_kernel(
    const double *restrict f_in, const double *restrict k_array, const double h, 
    double *restrict f_out, double *restrict f_err, const int batch_size, const int batch_id
) {
    #pragma unroll
    for (int i = 0; i < 9; ++i) {
        // The 4th-order candidate state, used exclusively as the baseline for error estimation.
        const double f_4th = f_in[IDX_LOCAL(i, batch_id, batch_size)] + h * ( 
                                       (25.0/216.0) * k_array[IDX_LOCAL(0*9 + i, batch_id, batch_size)] +
                                       (1408.0/2565.0) * k_array[IDX_LOCAL(2*9 + i, batch_id, batch_size)] +
                                       (2197.0/4104.0) * k_array[IDX_LOCAL(3*9 + i, batch_id, batch_size)] -
                                       (1.0/5.0) * k_array[IDX_LOCAL(4*9 + i, batch_id, batch_size)] );

        // The 5th-order candidate state, representing the newly proposed physical configuration.
        f_out[IDX_LOCAL(i, batch_id, batch_size)] = f_in[IDX_LOCAL(i, batch_id, batch_size)] + h * ( 
                                   (16.0/135.0) * k_array[IDX_LOCAL(0*9 + i, batch_id, batch_size)] +
                                   (6656.0/12825.0) * k_array[IDX_LOCAL(2*9 + i, batch_id, batch_size)] +
                                   (28561.0/56430.0) * k_array[IDX_LOCAL(3*9 + i, batch_id, batch_size)] -
                                   (9.0/50.0) * k_array[IDX_LOCAL(4*9 + i, batch_id, batch_size)] +
                                   (2.0/55.0) * k_array[IDX_LOCAL(5*9 + i, batch_id, batch_size)] );

        // The calculated local truncation error derived from the embedded method's order disparity.
        f_err[IDX_LOCAL(i, batch_id, batch_size)] = f_out[IDX_LOCAL(i, batch_id, batch_size)] - f_4th;
    }
}
"""
    gpu_c_code_for_header = (
        "#pragma omp declare target\n"
        + c_code_for_header
        + "\n#pragma omp end declare target"
    )
    Bdefines_h.register_BHaH_defines("rkf45_helpers", gpu_c_code_for_header)