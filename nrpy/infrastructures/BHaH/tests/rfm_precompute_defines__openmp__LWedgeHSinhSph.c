#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0_host.
 * Kernel to precompute metric quantity f0_of_xx0.
 */
static void rfm_precompute_defines__f0_of_xx0_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPL = params->AMPL;
  const REAL SINHW = params->SINHW;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0[i0] = AMPL * (expf(xx0 / SINHW) - expf(-xx0 / SINHW)) / (expf(1.0F / SINHW) - expf(-1 / SINHW));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__D0_host.
 * Kernel to precompute metric quantity f0_of_xx0__D0.
 */
static void rfm_precompute_defines__f0_of_xx0__D0_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                       const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPL = params->AMPL;
  const REAL SINHW = params->SINHW;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__D0[i0] = AMPL * (expf(xx0 / SINHW) / SINHW + expf(-xx0 / SINHW) / SINHW) / (expf(1.0F / SINHW) - expf(-1 / SINHW));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__D0_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DD00_host.
 * Kernel to precompute metric quantity f0_of_xx0__DD00.
 */
static void rfm_precompute_defines__f0_of_xx0__DD00_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPL = params->AMPL;
  const REAL SINHW = params->SINHW;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DD00[i0] =
        AMPL * (expf(xx0 / SINHW) / powf(SINHW, 2) - expf(-xx0 / SINHW) / powf(SINHW, 2)) / (expf(1.0F / SINHW) - expf(-1 / SINHW));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DD00_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DDD000_host.
 * Kernel to precompute metric quantity f0_of_xx0__DDD000.
 */
static void rfm_precompute_defines__f0_of_xx0__DDD000_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                           const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPL = params->AMPL;
  const REAL SINHW = params->SINHW;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DDD000[i0] =
        AMPL * (expf(xx0 / SINHW) / powf(SINHW, 3) + expf(-xx0 / SINHW) / powf(SINHW, 3)) / (expf(1.0F / SINHW) - expf(-1 / SINHW));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DDD000_host
/**
 * Kernel: rfm_precompute_defines__f1_of_xx1_host.
 * Kernel to precompute metric quantity f1_of_xx1.
 */
static void rfm_precompute_defines__f1_of_xx1_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct, const REAL *restrict x1) {
  // Temporary parameters
  for (int i1 = 0; i1 < params->Nxx_plus_2NGHOSTS1; i1++) {
    const REAL xx1 = x1[i1];
    rfmstruct->f1_of_xx1[i1] = sinf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f1_of_xx1_host
/**
 * Kernel: rfm_precompute_defines__f1_of_xx1__D1_host.
 * Kernel to precompute metric quantity f1_of_xx1__D1.
 */
static void rfm_precompute_defines__f1_of_xx1__D1_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                       const REAL *restrict x1) {
  // Temporary parameters
  for (int i1 = 0; i1 < params->Nxx_plus_2NGHOSTS1; i1++) {
    const REAL xx1 = x1[i1];
    rfmstruct->f1_of_xx1__D1[i1] = cosf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f1_of_xx1__D1_host
/**
 * Kernel: rfm_precompute_defines__f1_of_xx1__DD11_host.
 * Kernel to precompute metric quantity f1_of_xx1__DD11.
 */
static void rfm_precompute_defines__f1_of_xx1__DD11_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict x1) {
  // Temporary parameters
  for (int i1 = 0; i1 < params->Nxx_plus_2NGHOSTS1; i1++) {
    const REAL xx1 = x1[i1];
    rfmstruct->f1_of_xx1__DD11[i1] = -sinf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f1_of_xx1__DD11_host

/**
 * rfm_precompute_defines: reference metric precomputed lookup arrays: defines
 */
void rfm_precompute_defines__rfm__LWedgeHSinhSph(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                 rfm_struct *restrict rfmstruct, REAL *restrict xx[3]) {
  MAYBE_UNUSED const REAL *restrict x0 = xx[0];
  MAYBE_UNUSED const REAL *restrict x1 = xx[1];
  MAYBE_UNUSED const REAL *restrict x2 = xx[2];
  rfm_precompute_defines__f0_of_xx0_host(params, rfmstruct, x0);
  rfm_precompute_defines__f0_of_xx0__D0_host(params, rfmstruct, x0);
  rfm_precompute_defines__f0_of_xx0__DD00_host(params, rfmstruct, x0);
  rfm_precompute_defines__f0_of_xx0__DDD000_host(params, rfmstruct, x0);
  rfm_precompute_defines__f1_of_xx1_host(params, rfmstruct, x1);
  rfm_precompute_defines__f1_of_xx1__D1_host(params, rfmstruct, x1);
  rfm_precompute_defines__f1_of_xx1__DD11_host(params, rfmstruct, x1);
} // END FUNCTION rfm_precompute_defines__rfm__LWedgeHSinhSph
