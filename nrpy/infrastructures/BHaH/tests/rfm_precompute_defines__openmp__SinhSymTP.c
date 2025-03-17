#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0_host.
 * Kernel to precompute metric quantity f0_of_xx0.
 */
static void rfm_precompute_defines__f0_of_xx0_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = params->AMAX;
  const REAL SINHWAA = params->SINHWAA;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0[i0] = AMAX * (expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA)) / (expf(1.0F / SINHWAA) - expf(-1 / SINHWAA));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__D0_host.
 * Kernel to precompute metric quantity f0_of_xx0__D0.
 */
static void rfm_precompute_defines__f0_of_xx0__D0_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                       const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = params->AMAX;
  const REAL SINHWAA = params->SINHWAA;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__D0[i0] =
        AMAX * (expf(xx0 / SINHWAA) / SINHWAA + expf(-xx0 / SINHWAA) / SINHWAA) / (expf(1.0F / SINHWAA) - expf(-1 / SINHWAA));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__D0_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DD00_host.
 * Kernel to precompute metric quantity f0_of_xx0__DD00.
 */
static void rfm_precompute_defines__f0_of_xx0__DD00_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = params->AMAX;
  const REAL SINHWAA = params->SINHWAA;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DD00[i0] =
        AMAX * (expf(xx0 / SINHWAA) / powf(SINHWAA, 2) - expf(-xx0 / SINHWAA) / powf(SINHWAA, 2)) / (expf(1.0F / SINHWAA) - expf(-1 / SINHWAA));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DD00_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DDD000_host.
 * Kernel to precompute metric quantity f0_of_xx0__DDD000.
 */
static void rfm_precompute_defines__f0_of_xx0__DDD000_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                           const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = params->AMAX;
  const REAL SINHWAA = params->SINHWAA;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DDD000[i0] =
        AMAX * (expf(xx0 / SINHWAA) / powf(SINHWAA, 3) + expf(-xx0 / SINHWAA) / powf(SINHWAA, 3)) / (expf(1.0F / SINHWAA) - expf(-1 / SINHWAA));
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
 * Kernel: rfm_precompute_defines__f2_of_xx0_host.
 * Kernel to precompute metric quantity f2_of_xx0.
 */
static void rfm_precompute_defines__f2_of_xx0_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = params->AMAX;
  const REAL SINHWAA = params->SINHWAA;
  const REAL bScale = params->bScale;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f2_of_xx0[i0] = sqrtf(
        powf(AMAX, 2) * powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) / powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2) + powf(bScale, 2));
  }
} // END FUNCTION rfm_precompute_defines__f2_of_xx0_host
/**
 * Kernel: rfm_precompute_defines__f2_of_xx0__D0_host.
 * Kernel to precompute metric quantity f2_of_xx0__D0.
 */
static void rfm_precompute_defines__f2_of_xx0__D0_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                       const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = params->AMAX;
  const REAL SINHWAA = params->SINHWAA;
  const REAL bScale = params->bScale;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f2_of_xx0__D0[i0] =
        (1.0F / 2.0F) * powf(AMAX, 2) * (2 * expf(xx0 / SINHWAA) / SINHWAA + 2 * expf(-xx0 / SINHWAA) / SINHWAA) *
        (expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA)) /
        (sqrtf(powf(AMAX, 2) * powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) / powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2) +
               powf(bScale, 2)) *
         powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2));
  }
} // END FUNCTION rfm_precompute_defines__f2_of_xx0__D0_host
/**
 * Kernel: rfm_precompute_defines__f2_of_xx0__DD00_host.
 * Kernel to precompute metric quantity f2_of_xx0__DD00.
 */
static void rfm_precompute_defines__f2_of_xx0__DD00_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = params->AMAX;
  const REAL SINHWAA = params->SINHWAA;
  const REAL bScale = params->bScale;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f2_of_xx0__DD00[i0] =
        -1.0F / 4.0F * powf(AMAX, 4) * powf(2 * expf(xx0 / SINHWAA) / SINHWAA + 2 * expf(-xx0 / SINHWAA) / SINHWAA, 2) *
            powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) /
            (powf(powf(AMAX, 2) * powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) / powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2) +
                      powf(bScale, 2),
                  3.0F / 2.0F) *
             powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 4)) +
        (1.0F / 2.0F) * powf(AMAX, 2) * (2 * expf(xx0 / SINHWAA) / powf(SINHWAA, 2) - 2 * expf(-xx0 / SINHWAA) / powf(SINHWAA, 2)) *
            (expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA)) /
            (sqrtf(powf(AMAX, 2) * powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) / powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2) +
                   powf(bScale, 2)) *
             powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2)) +
        (1.0F / 2.0F) * powf(AMAX, 2) * (expf(xx0 / SINHWAA) / SINHWAA + expf(-xx0 / SINHWAA) / SINHWAA) *
            (2 * expf(xx0 / SINHWAA) / SINHWAA + 2 * expf(-xx0 / SINHWAA) / SINHWAA) /
            (sqrtf(powf(AMAX, 2) * powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) / powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2) +
                   powf(bScale, 2)) *
             powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2));
  }
} // END FUNCTION rfm_precompute_defines__f2_of_xx0__DD00_host
/**
 * Kernel: rfm_precompute_defines__f4_of_xx1_host.
 * Kernel to precompute metric quantity f4_of_xx1.
 */
static void rfm_precompute_defines__f4_of_xx1_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct, const REAL *restrict x1) {
  // Temporary parameters
  const REAL bScale = params->bScale;
  for (int i1 = 0; i1 < params->Nxx_plus_2NGHOSTS1; i1++) {
    const REAL xx1 = x1[i1];
    rfmstruct->f4_of_xx1[i1] = bScale * sinf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f4_of_xx1_host
/**
 * Kernel: rfm_precompute_defines__f4_of_xx1__D1_host.
 * Kernel to precompute metric quantity f4_of_xx1__D1.
 */
static void rfm_precompute_defines__f4_of_xx1__D1_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                       const REAL *restrict x1) {
  // Temporary parameters
  const REAL bScale = params->bScale;
  for (int i1 = 0; i1 < params->Nxx_plus_2NGHOSTS1; i1++) {
    const REAL xx1 = x1[i1];
    rfmstruct->f4_of_xx1__D1[i1] = bScale * cosf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f4_of_xx1__D1_host
/**
 * Kernel: rfm_precompute_defines__f4_of_xx1__DD11_host.
 * Kernel to precompute metric quantity f4_of_xx1__DD11.
 */
static void rfm_precompute_defines__f4_of_xx1__DD11_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict x1) {
  // Temporary parameters
  const REAL bScale = params->bScale;
  for (int i1 = 0; i1 < params->Nxx_plus_2NGHOSTS1; i1++) {
    const REAL xx1 = x1[i1];
    rfmstruct->f4_of_xx1__DD11[i1] = -bScale * sinf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f4_of_xx1__DD11_host

/**
 * rfm_precompute_defines: reference metric precomputed lookup arrays: defines
 */
void rfm_precompute_defines__rfm__SinhSymTP(const commondata_struct *restrict commondata, const params_struct *restrict params,
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
  rfm_precompute_defines__f2_of_xx0_host(params, rfmstruct, x0);
  rfm_precompute_defines__f2_of_xx0__D0_host(params, rfmstruct, x0);
  rfm_precompute_defines__f2_of_xx0__DD00_host(params, rfmstruct, x0);
  rfm_precompute_defines__f4_of_xx1_host(params, rfmstruct, x1);
  rfm_precompute_defines__f4_of_xx1__D1_host(params, rfmstruct, x1);
  rfm_precompute_defines__f4_of_xx1__DD11_host(params, rfmstruct, x1);
} // END FUNCTION rfm_precompute_defines__rfm__SinhSymTP
