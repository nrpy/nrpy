"""
Generates the C adaptive step-size controller, updated for BHaH integration.
Features a high-performance L_infinity norm and L1 momentum floor for relativistic ray tracing.
Adapted for batched SoA global arrays to ensure optimal memory throughput.

Author: Dalton J. Moone
"""

import nrpy.params as par
from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h


def rkf45_update_and_control_helper() -> None:
    """Generate the high-performance adaptive controller."""
    # --- Local Parameter Registration ---
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
    # ------------------------------------
    # Raw string prevents python from trying to evaluate C blocks
    c_code_for_header = r"""
static inline bool update_photon_state_and_stepsize(
    PhotonStateSoA *restrict all_photons, 
    const long int num_rays, 
    const long int photon_idx, 
    const double *restrict f_start, 
    const double *restrict f_out, 
    const double *restrict f_err, 
    const int batch_size, 
    const int batch_id, 
    const commondata_struct *restrict commondata
) {
    const double rtol = commondata->rkf45_error_tolerance;
    const double atol = commondata->rkf45_absolute_error_tolerance;
    
    // Initialize the L_infinity (maximum) error norm
    double err_norm = 0.0;
    double ratio;

    // 1. Coordinate Time (f[0]): Pure absolute tolerance to prevent secular drift
    ratio = fabs(f_err[IDX_LOCAL(0, batch_id, batch_size)]) / atol;
    err_norm = fmax(err_norm, ratio);

    // 2. Spatial Coordinates (f[1], f[2], f[3]): Standard mixed norm
    for (int i = 1; i <= 3; ++i) {
        double scale = atol + rtol * fabs(f_start[IDX_LOCAL(i, batch_id, batch_size)]);
        ratio = fabs(f_err[IDX_LOCAL(i, batch_id, batch_size)]) / scale;
        err_norm = fmax(err_norm, ratio);
    }

    // 3. Temporal Momentum p^t (f[4]): Standard mixed norm
    // Represents energy evolution; does not suffer from zero-crossings
    double scale_pt = atol + rtol * fabs(f_start[IDX_LOCAL(4, batch_id, batch_size)]);
    ratio = fabs(f_err[IDX_LOCAL(4, batch_id, batch_size)]) / scale_pt;
    err_norm = fmax(err_norm, ratio);

    // 4. Spatial Momenta p^x, p^y, p^z (f[5], f[6], f[7]): L1 floored mixed norm
    // Calculate the L1 Manhattan distance to serve as a robust, computationally cheap physical floor
    double p_L1 = fabs(f_start[IDX_LOCAL(5, batch_id, batch_size)]) + 
                  fabs(f_start[IDX_LOCAL(6, batch_id, batch_size)]) + 
                  fabs(f_start[IDX_LOCAL(7, batch_id, batch_size)]);
    
    for (int i = 5; i <= 7; ++i) {
        // p_L1 inherently bounds the individual components, eliminating the need for fmax() here
        double scale = atol + rtol * p_L1; 
        ratio = fabs(f_err[IDX_LOCAL(i, batch_id, batch_size)]) / scale;
        err_norm = fmax(err_norm, ratio);
    }

    // NOTE: Path length f[8] is explicitly excluded from the error controller.

    // Step acceptance is determined strictly by the single worst-performing component
    bool step_accepted = (err_norm <= 1.0);

    // Compute next h using the unclamped 1/5 power law
    double h_new;
    if (err_norm > 1e-15) {
        h_new = commondata->rkf45_safety_factor * all_photons->h[photon_idx] * pow(1.0 / err_norm, 0.2);
    } else {
        h_new = 2.0 * all_photons->h[photon_idx];
    }

    // Apply strict global boundaries
    h_new = fmax(h_new, commondata->rkf45_h_min);
    h_new = fmin(h_new, commondata->rkf45_h_max);

    // Memory commitment to SoA Global Arrays
    if (step_accepted) {
        for (int i = 0; i < 9; ++i) {
            all_photons->f[IDX_GLOBAL(i, photon_idx, num_rays)] = f_out[IDX_LOCAL(i, batch_id, batch_size)];
        }
        all_photons->affine_param[photon_idx] += all_photons->h[photon_idx];
        all_photons->rejection_retries[photon_idx] = 0;
    } else {
        all_photons->rejection_retries[photon_idx]++;
    }

    all_photons->h[photon_idx] = h_new;
    return step_accepted;
}
"""
    gpu_c_code_for_header = "#pragma omp declare target\n" + c_code_for_header + "\n#pragma omp end declare target"

    Bdefines_h.register_BHaH_defines("rkf45_update_control", gpu_c_code_for_header)