#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0_host.
 * Kernel to precompute metric quantity f0_of_xx0.
 */
static void rfm_precompute_defines__f0_of_xx0_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLXYZ = params->AMPLXYZ;
  const REAL SINHWXYZ = params->SINHWXYZ;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0[i0] =
        AMPLXYZ * (expf(xx0 / SINHWXYZ) / SINHWXYZ + expf(-xx0 / SINHWXYZ) / SINHWXYZ) / (expf(1.0F / SINHWXYZ) - expf(-1 / SINHWXYZ));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__D0_host.
 * Kernel to precompute metric quantity f0_of_xx0__D0.
 */
static void rfm_precompute_defines__f0_of_xx0__D0_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                       const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLXYZ = params->AMPLXYZ;
  const REAL SINHWXYZ = params->SINHWXYZ;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__D0[i0] = AMPLXYZ * (expf(xx0 / SINHWXYZ) / powf(SINHWXYZ, 2) - expf(-xx0 / SINHWXYZ) / powf(SINHWXYZ, 2)) /
                                   (expf(1.0F / SINHWXYZ) - expf(-1 / SINHWXYZ));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__D0_host
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DD00_host.
 * Kernel to precompute metric quantity f0_of_xx0__DD00.
 */
static void rfm_precompute_defines__f0_of_xx0__DD00_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLXYZ = params->AMPLXYZ;
  const REAL SINHWXYZ = params->SINHWXYZ;
  for (int i0 = 0; i0 < params->Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DD00[i0] = AMPLXYZ * (expf(xx0 / SINHWXYZ) / powf(SINHWXYZ, 3) + expf(-xx0 / SINHWXYZ) / powf(SINHWXYZ, 3)) /
                                     (expf(1.0F / SINHWXYZ) - expf(-1 / SINHWXYZ));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DD00_host
/**
 * Kernel: rfm_precompute_defines__f1_of_xx1_host.
 * Kernel to precompute metric quantity f1_of_xx1.
 */
static void rfm_precompute_defines__f1_of_xx1_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct, const REAL *restrict x1) {
  // Temporary parameters
  const REAL AMPLXYZ = params->AMPLXYZ;
  const REAL SINHWXYZ = params->SINHWXYZ;
  for (int i1 = 0; i1 < params->Nxx_plus_2NGHOSTS1; i1++) {
    const REAL xx1 = x1[i1];
    rfmstruct->f1_of_xx1[i1] =
        AMPLXYZ * (expf(xx1 / SINHWXYZ) / SINHWXYZ + expf(-xx1 / SINHWXYZ) / SINHWXYZ) / (expf(1.0F / SINHWXYZ) - expf(-1 / SINHWXYZ));
  }
} // END FUNCTION rfm_precompute_defines__f1_of_xx1_host
/**
 * Kernel: rfm_precompute_defines__f1_of_xx1__D1_host.
 * Kernel to precompute metric quantity f1_of_xx1__D1.
 */
static void rfm_precompute_defines__f1_of_xx1__D1_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                       const REAL *restrict x1) {
  // Temporary parameters
  const REAL AMPLXYZ = params->AMPLXYZ;
  const REAL SINHWXYZ = params->SINHWXYZ;
  for (int i1 = 0; i1 < params->Nxx_plus_2NGHOSTS1; i1++) {
    const REAL xx1 = x1[i1];
    rfmstruct->f1_of_xx1__D1[i1] = AMPLXYZ * (expf(xx1 / SINHWXYZ) / powf(SINHWXYZ, 2) - expf(-xx1 / SINHWXYZ) / powf(SINHWXYZ, 2)) /
                                   (expf(1.0F / SINHWXYZ) - expf(-1 / SINHWXYZ));
  }
} // END FUNCTION rfm_precompute_defines__f1_of_xx1__D1_host
/**
 * Kernel: rfm_precompute_defines__f1_of_xx1__DD11_host.
 * Kernel to precompute metric quantity f1_of_xx1__DD11.
 */
static void rfm_precompute_defines__f1_of_xx1__DD11_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict x1) {
  // Temporary parameters
  const REAL AMPLXYZ = params->AMPLXYZ;
  const REAL SINHWXYZ = params->SINHWXYZ;
  for (int i1 = 0; i1 < params->Nxx_plus_2NGHOSTS1; i1++) {
    const REAL xx1 = x1[i1];
    rfmstruct->f1_of_xx1__DD11[i1] = AMPLXYZ * (expf(xx1 / SINHWXYZ) / powf(SINHWXYZ, 3) + expf(-xx1 / SINHWXYZ) / powf(SINHWXYZ, 3)) /
                                     (expf(1.0F / SINHWXYZ) - expf(-1 / SINHWXYZ));
  }
} // END FUNCTION rfm_precompute_defines__f1_of_xx1__DD11_host
/**
 * Kernel: rfm_precompute_defines__f3_of_xx2_host.
 * Kernel to precompute metric quantity f3_of_xx2.
 */
static void rfm_precompute_defines__f3_of_xx2_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct, const REAL *restrict x2) {
  // Temporary parameters
  const REAL AMPLXYZ = params->AMPLXYZ;
  const REAL SINHWXYZ = params->SINHWXYZ;
  for (int i2 = 0; i2 < params->Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = x2[i2];
    rfmstruct->f3_of_xx2[i2] =
        AMPLXYZ * (expf(xx2 / SINHWXYZ) / SINHWXYZ + expf(-xx2 / SINHWXYZ) / SINHWXYZ) / (expf(1.0F / SINHWXYZ) - expf(-1 / SINHWXYZ));
  }
} // END FUNCTION rfm_precompute_defines__f3_of_xx2_host
/**
 * Kernel: rfm_precompute_defines__f3_of_xx2__D2_host.
 * Kernel to precompute metric quantity f3_of_xx2__D2.
 */
static void rfm_precompute_defines__f3_of_xx2__D2_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                       const REAL *restrict x2) {
  // Temporary parameters
  const REAL AMPLXYZ = params->AMPLXYZ;
  const REAL SINHWXYZ = params->SINHWXYZ;
  for (int i2 = 0; i2 < params->Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = x2[i2];
    rfmstruct->f3_of_xx2__D2[i2] = AMPLXYZ * (expf(xx2 / SINHWXYZ) / powf(SINHWXYZ, 2) - expf(-xx2 / SINHWXYZ) / powf(SINHWXYZ, 2)) /
                                   (expf(1.0F / SINHWXYZ) - expf(-1 / SINHWXYZ));
  }
} // END FUNCTION rfm_precompute_defines__f3_of_xx2__D2_host
/**
 * Kernel: rfm_precompute_defines__f3_of_xx2__DD22_host.
 * Kernel to precompute metric quantity f3_of_xx2__DD22.
 */
static void rfm_precompute_defines__f3_of_xx2__DD22_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                                                         const REAL *restrict x2) {
  // Temporary parameters
  const REAL AMPLXYZ = params->AMPLXYZ;
  const REAL SINHWXYZ = params->SINHWXYZ;
  for (int i2 = 0; i2 < params->Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = x2[i2];
    rfmstruct->f3_of_xx2__DD22[i2] = AMPLXYZ * (expf(xx2 / SINHWXYZ) / powf(SINHWXYZ, 3) + expf(-xx2 / SINHWXYZ) / powf(SINHWXYZ, 3)) /
                                     (expf(1.0F / SINHWXYZ) - expf(-1 / SINHWXYZ));
  }
} // END FUNCTION rfm_precompute_defines__f3_of_xx2__DD22_host

/**
 * rfm_precompute_defines: reference metric precomputed lookup arrays: defines
 */
void rfm_precompute_defines__rfm__SinhCartesian(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                rfm_struct *restrict rfmstruct, REAL *restrict xx[3]) {
  MAYBE_UNUSED const REAL *restrict x0 = xx[0];
  MAYBE_UNUSED const REAL *restrict x1 = xx[1];
  MAYBE_UNUSED const REAL *restrict x2 = xx[2];
  rfm_precompute_defines__f0_of_xx0_host(params, rfmstruct, x0);
  rfm_precompute_defines__f0_of_xx0__D0_host(params, rfmstruct, x0);
  rfm_precompute_defines__f0_of_xx0__DD00_host(params, rfmstruct, x0);
  rfm_precompute_defines__f1_of_xx1_host(params, rfmstruct, x1);
  rfm_precompute_defines__f1_of_xx1__D1_host(params, rfmstruct, x1);
  rfm_precompute_defines__f1_of_xx1__DD11_host(params, rfmstruct, x1);
  rfm_precompute_defines__f3_of_xx2_host(params, rfmstruct, x2);
  rfm_precompute_defines__f3_of_xx2__D2_host(params, rfmstruct, x2);
  rfm_precompute_defines__f3_of_xx2__DD22_host(params, rfmstruct, x2);
} // END FUNCTION rfm_precompute_defines__rfm__SinhCartesian
