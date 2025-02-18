#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0_host.
 * Kernel to precompute metric quantity f0_of_xx0.
 */
static void rfm_precompute_defines__f0_of_xx0_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLRHO = params->AMPLRHO;
  const REAL SINHWRHO = params->SINHWRHO;
  const REAL rho_slope = params->rho_slope;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0[i0] = rho_slope * xx0 + powf(xx0, 2) * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) - expf(-xx0 / SINHWRHO)) /
                                                     (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__D0_host.
 * Kernel to precompute metric quantity f0_of_xx0__D0.
 */
static void rfm_precompute_defines__f0_of_xx0__D0_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                       const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLRHO = params->AMPLRHO;
  const REAL SINHWRHO = params->SINHWRHO;
  const REAL rho_slope = params->rho_slope;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__D0[i0] =
        rho_slope +
        powf(xx0, 2) * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / SINHWRHO + expf(-xx0 / SINHWRHO) / SINHWRHO) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO)) +
        2 * xx0 * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) - expf(-xx0 / SINHWRHO)) / (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__D0_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DD00_host.
 * Kernel to precompute metric quantity f0_of_xx0__DD00.
 */
static void rfm_precompute_defines__f0_of_xx0__DD00_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLRHO = params->AMPLRHO;
  const REAL SINHWRHO = params->SINHWRHO;
  const REAL rho_slope = params->rho_slope;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DD00[i0] =
        powf(xx0, 2) * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / powf(SINHWRHO, 2) - expf(-xx0 / SINHWRHO) / powf(SINHWRHO, 2)) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO)) +
        4 * xx0 * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / SINHWRHO + expf(-xx0 / SINHWRHO) / SINHWRHO) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO)) +
        2 * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) - expf(-xx0 / SINHWRHO)) / (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DD00_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DDD000_host.
 * Kernel to precompute metric quantity f0_of_xx0__DDD000.
 */
static void rfm_precompute_defines__f0_of_xx0__DDD000_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                           const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLRHO = params->AMPLRHO;
  const REAL SINHWRHO = params->SINHWRHO;
  const REAL rho_slope = params->rho_slope;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DDD000[i0] =
        powf(xx0, 2) * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / powf(SINHWRHO, 3) + expf(-xx0 / SINHWRHO) / powf(SINHWRHO, 3)) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO)) +
        6 * xx0 * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / powf(SINHWRHO, 2) - expf(-xx0 / SINHWRHO) / powf(SINHWRHO, 2)) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO)) +
        6 * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / SINHWRHO + expf(-xx0 / SINHWRHO) / SINHWRHO) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DDD000_host
/**
 * Kernel: rfm_precompute_defines__f3_of_xx2_host.
 * Kernel to precompute metric quantity f3_of_xx2.
 */
static void rfm_precompute_defines__f3_of_xx2_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct, const REAL *restrict x2) {
  // Temporary parameters
  const REAL AMPLZ = params->AMPLZ;
  const REAL SINHWZ = params->SINHWZ;
  const REAL z_slope = params->z_slope;
  for (int i2 = 0; i2 < params->Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = x2[i2];
    rfmstruct->f3_of_xx2[i2] =
        powf(xx2, 2) * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / SINHWZ + expf(-xx2 / SINHWZ) / SINHWZ) / (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) +
        2 * xx2 * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) - expf(-xx2 / SINHWZ)) / (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) + z_slope;
  }
} // END FUNCTION rfm_precompute_defines__f3_of_xx2_host
/**
 * Kernel: rfm_precompute_defines__f3_of_xx2__D2_host.
 * Kernel to precompute metric quantity f3_of_xx2__D2.
 */
static void rfm_precompute_defines__f3_of_xx2__D2_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                       const REAL *restrict x2) {
  // Temporary parameters
  const REAL AMPLZ = params->AMPLZ;
  const REAL SINHWZ = params->SINHWZ;
  const REAL z_slope = params->z_slope;
  for (int i2 = 0; i2 < params->Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = x2[i2];
    rfmstruct->f3_of_xx2__D2[i2] =
        powf(xx2, 2) * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / powf(SINHWZ, 2) - expf(-xx2 / SINHWZ) / powf(SINHWZ, 2)) /
            (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) +
        4 * xx2 * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / SINHWZ + expf(-xx2 / SINHWZ) / SINHWZ) / (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) +
        2 * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) - expf(-xx2 / SINHWZ)) / (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ));
  }
} // END FUNCTION rfm_precompute_defines__f3_of_xx2__D2_host
/**
 * Kernel: rfm_precompute_defines__f3_of_xx2__DD22_host.
 * Kernel to precompute metric quantity f3_of_xx2__DD22.
 */
static void rfm_precompute_defines__f3_of_xx2__DD22_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict x2) {
  // Temporary parameters
  const REAL AMPLZ = params->AMPLZ;
  const REAL SINHWZ = params->SINHWZ;
  const REAL z_slope = params->z_slope;
  for (int i2 = 0; i2 < params->Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = x2[i2];
    rfmstruct->f3_of_xx2__DD22[i2] =
        powf(xx2, 2) * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / powf(SINHWZ, 3) + expf(-xx2 / SINHWZ) / powf(SINHWZ, 3)) /
            (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) +
        6 * xx2 * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / powf(SINHWZ, 2) - expf(-xx2 / SINHWZ) / powf(SINHWZ, 2)) /
            (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) +
        6 * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / SINHWZ + expf(-xx2 / SINHWZ) / SINHWZ) / (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ));
  }
} // END FUNCTION rfm_precompute_defines__f3_of_xx2__DD22_host

/**
 * rfm_precompute_defines: reference metric precomputed lookup arrays: defines
 */
void rfm_precompute_defines__rfm__SinhCylindricalv2n2(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                      rfm_struct *restrict rfmstruct, REAL *restrict xx[3]) {
  MAYBE_UNUSED const REAL *restrict x0 = xx[0];
  MAYBE_UNUSED const REAL *restrict x1 = xx[1];
  MAYBE_UNUSED const REAL *restrict x2 = xx[2];
  rfm_precompute_defines__f0_of_xx0_host(params, rfmstruct, x0);
  rfm_precompute_defines__f0_of_xx0__D0_host(params, rfmstruct, x0);
  rfm_precompute_defines__f0_of_xx0__DD00_host(params, rfmstruct, x0);
  rfm_precompute_defines__f0_of_xx0__DDD000_host(params, rfmstruct, x0);
  rfm_precompute_defines__f3_of_xx2_host(params, rfmstruct, x2);
  rfm_precompute_defines__f3_of_xx2__D2_host(params, rfmstruct, x2);
  rfm_precompute_defines__f3_of_xx2__DD22_host(params, rfmstruct, x2);
} // END FUNCTION rfm_precompute_defines__rfm__SinhCylindricalv2n2
