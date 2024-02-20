#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"

// ADM variables in the Cartesian basis:
typedef struct __ADM_Cart_basis_struct__ {
  REAL alpha, betaU0, betaU1, betaU2, BU0, BU1, BU2;
  REAL gammaDD00, gammaDD01, gammaDD02, gammaDD11, gammaDD12, gammaDD22;
  REAL KDD00, KDD01, KDD02, KDD11, KDD12, KDD22;
} ADM_Cart_basis_struct;

// BSSN variables in the Cartesian basis:
typedef struct __BSSN_Cart_basis_struct__ {
  REAL alpha, betaU0, betaU1, betaU2, BU0, BU1, BU2;
  REAL cf, trK;
  REAL gammabarDD00, gammabarDD01, gammabarDD02, gammabarDD11, gammabarDD12, gammabarDD22;
  REAL AbarDD00, AbarDD01, AbarDD02, AbarDD11, AbarDD12, AbarDD22;
} BSSN_Cart_basis_struct;

// Rescaled BSSN variables in the rfm basis:
typedef struct __rescaled_BSSN_rfm_basis_struct__ {
  REAL alpha, vetU0, vetU1, vetU2, betU0, betU1, betU2;
  REAL cf, trK;
  REAL hDD00, hDD01, hDD02, hDD11, hDD12, hDD22;
  REAL aDD00, aDD01, aDD02, aDD11, aDD12, aDD22;
} rescaled_BSSN_rfm_basis_struct;
/*
 * Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis
 */
static void ADM_SphorCart_to_Cart(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                                  const initial_data_struct *restrict initial_data, ADM_Cart_basis_struct *restrict ADM_Cart_basis) {
#include "../set_CodeParameters.h"
  // Unpack initial_data for ADM vectors/tensors
  const REAL betaSphorCartU0 = initial_data->betaSphorCartU0;
  const REAL betaSphorCartU1 = initial_data->betaSphorCartU1;
  const REAL betaSphorCartU2 = initial_data->betaSphorCartU2;

  const REAL BSphorCartU0 = initial_data->BSphorCartU0;
  const REAL BSphorCartU1 = initial_data->BSphorCartU1;
  const REAL BSphorCartU2 = initial_data->BSphorCartU2;

  const REAL gammaSphorCartDD00 = initial_data->gammaSphorCartDD00;
  const REAL gammaSphorCartDD01 = initial_data->gammaSphorCartDD01;
  const REAL gammaSphorCartDD02 = initial_data->gammaSphorCartDD02;
  const REAL gammaSphorCartDD11 = initial_data->gammaSphorCartDD11;
  const REAL gammaSphorCartDD12 = initial_data->gammaSphorCartDD12;
  const REAL gammaSphorCartDD22 = initial_data->gammaSphorCartDD22;

  const REAL KSphorCartDD00 = initial_data->KSphorCartDD00;
  const REAL KSphorCartDD01 = initial_data->KSphorCartDD01;
  const REAL KSphorCartDD02 = initial_data->KSphorCartDD02;
  const REAL KSphorCartDD11 = initial_data->KSphorCartDD11;
  const REAL KSphorCartDD12 = initial_data->KSphorCartDD12;
  const REAL KSphorCartDD22 = initial_data->KSphorCartDD22;

  ADM_Cart_basis->BU0 = BSphorCartU0;
  ADM_Cart_basis->BU1 = BSphorCartU1;
  ADM_Cart_basis->BU2 = BSphorCartU2;
  ADM_Cart_basis->KDD00 = KSphorCartDD00;
  ADM_Cart_basis->KDD01 = KSphorCartDD01;
  ADM_Cart_basis->KDD02 = KSphorCartDD02;
  ADM_Cart_basis->KDD11 = KSphorCartDD11;
  ADM_Cart_basis->KDD12 = KSphorCartDD12;
  ADM_Cart_basis->KDD22 = KSphorCartDD22;
  ADM_Cart_basis->alpha = initial_data->alpha;
  ADM_Cart_basis->betaU0 = betaSphorCartU0;
  ADM_Cart_basis->betaU1 = betaSphorCartU1;
  ADM_Cart_basis->betaU2 = betaSphorCartU2;
  ADM_Cart_basis->gammaDD00 = gammaSphorCartDD00;
  ADM_Cart_basis->gammaDD01 = gammaSphorCartDD01;
  ADM_Cart_basis->gammaDD02 = gammaSphorCartDD02;
  ADM_Cart_basis->gammaDD11 = gammaSphorCartDD11;
  ADM_Cart_basis->gammaDD12 = gammaSphorCartDD12;
  ADM_Cart_basis->gammaDD22 = gammaSphorCartDD22;
}
/*
 * Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis
 */
static void ADM_Cart_to_BSSN_Cart(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                                  const ADM_Cart_basis_struct *restrict ADM_Cart_basis, BSSN_Cart_basis_struct *restrict BSSN_Cart_basis) {

  // *In the Cartesian basis*, convert ADM quantities gammaDD & KDD
  //   into BSSN gammabarDD, AbarDD, cf, and trK.
  BSSN_Cart_basis->alpha = ADM_Cart_basis->alpha;
  BSSN_Cart_basis->betaU0 = ADM_Cart_basis->betaU0;
  BSSN_Cart_basis->betaU1 = ADM_Cart_basis->betaU1;
  BSSN_Cart_basis->betaU2 = ADM_Cart_basis->betaU2;
  BSSN_Cart_basis->BU0 = ADM_Cart_basis->BU0;
  BSSN_Cart_basis->BU1 = ADM_Cart_basis->BU1;
  BSSN_Cart_basis->BU2 = ADM_Cart_basis->BU2;
  const REAL tmp1 = ADM_Cart_basis->gammaDD00 * ADM_Cart_basis->gammaDD11 * ADM_Cart_basis->gammaDD22;
  const REAL tmp3 = 2 * ADM_Cart_basis->gammaDD01 * ADM_Cart_basis->gammaDD02 * ADM_Cart_basis->gammaDD12;
  const REAL tmp5 = ADM_Cart_basis->gammaDD00 * ((ADM_Cart_basis->gammaDD12) * (ADM_Cart_basis->gammaDD12));
  const REAL tmp7 = ((ADM_Cart_basis->gammaDD01) * (ADM_Cart_basis->gammaDD01)) * ADM_Cart_basis->gammaDD22;
  const REAL tmp9 = ((ADM_Cart_basis->gammaDD02) * (ADM_Cart_basis->gammaDD02)) * ADM_Cart_basis->gammaDD11;
  const REAL tmp10 = tmp1 + tmp3 - tmp5 - tmp7 - tmp9;
  const REAL tmp11 = (1.0 / (tmp10));
  const REAL tmp12 = cbrt(tmp11);
  const REAL tmp13 = 2 * tmp11;
  const REAL tmp14 = ADM_Cart_basis->KDD00 * tmp11 *
                         (ADM_Cart_basis->gammaDD11 * ADM_Cart_basis->gammaDD22 - ((ADM_Cart_basis->gammaDD12) * (ADM_Cart_basis->gammaDD12))) +
                     ADM_Cart_basis->KDD01 * tmp13 *
                         (-ADM_Cart_basis->gammaDD01 * ADM_Cart_basis->gammaDD22 + ADM_Cart_basis->gammaDD02 * ADM_Cart_basis->gammaDD12) +
                     ADM_Cart_basis->KDD02 * tmp13 *
                         (ADM_Cart_basis->gammaDD01 * ADM_Cart_basis->gammaDD12 - ADM_Cart_basis->gammaDD02 * ADM_Cart_basis->gammaDD11) +
                     ADM_Cart_basis->KDD11 * tmp11 *
                         (ADM_Cart_basis->gammaDD00 * ADM_Cart_basis->gammaDD22 - ((ADM_Cart_basis->gammaDD02) * (ADM_Cart_basis->gammaDD02))) +
                     ADM_Cart_basis->KDD12 * tmp13 *
                         (-ADM_Cart_basis->gammaDD00 * ADM_Cart_basis->gammaDD12 + ADM_Cart_basis->gammaDD01 * ADM_Cart_basis->gammaDD02) +
                     ADM_Cart_basis->KDD22 * tmp11 *
                         (ADM_Cart_basis->gammaDD00 * ADM_Cart_basis->gammaDD11 - ((ADM_Cart_basis->gammaDD01) * (ADM_Cart_basis->gammaDD01)));
  const REAL tmp15 = (1.0 / 3.0) * tmp14;
  BSSN_Cart_basis->AbarDD00 = tmp12 * (ADM_Cart_basis->KDD00 - ADM_Cart_basis->gammaDD00 * tmp15);
  BSSN_Cart_basis->AbarDD01 = tmp12 * (ADM_Cart_basis->KDD01 - ADM_Cart_basis->gammaDD01 * tmp15);
  BSSN_Cart_basis->AbarDD02 = tmp12 * (ADM_Cart_basis->KDD02 - ADM_Cart_basis->gammaDD02 * tmp15);
  BSSN_Cart_basis->AbarDD11 = tmp12 * (ADM_Cart_basis->KDD11 - ADM_Cart_basis->gammaDD11 * tmp15);
  BSSN_Cart_basis->AbarDD12 = tmp12 * (ADM_Cart_basis->KDD12 - ADM_Cart_basis->gammaDD12 * tmp15);
  BSSN_Cart_basis->AbarDD22 = tmp12 * (ADM_Cart_basis->KDD22 - ADM_Cart_basis->gammaDD22 * tmp15);
  BSSN_Cart_basis->cf = pow(tmp10 / (tmp1 * tmp11 + tmp11 * tmp3 - tmp11 * tmp5 - tmp11 * tmp7 - tmp11 * tmp9), -1.0 / 6.0);
  BSSN_Cart_basis->gammabarDD00 = ADM_Cart_basis->gammaDD00 * tmp12;
  BSSN_Cart_basis->gammabarDD01 = ADM_Cart_basis->gammaDD01 * tmp12;
  BSSN_Cart_basis->gammabarDD02 = ADM_Cart_basis->gammaDD02 * tmp12;
  BSSN_Cart_basis->gammabarDD11 = ADM_Cart_basis->gammaDD11 * tmp12;
  BSSN_Cart_basis->gammabarDD12 = ADM_Cart_basis->gammaDD12 * tmp12;
  BSSN_Cart_basis->gammabarDD22 = ADM_Cart_basis->gammaDD22 * tmp12;
  BSSN_Cart_basis->trK = tmp14;
}
/*
 * Cartesian -> SinhSpherical basis transformation of BSSN vectors/tensors *except* lambda^i.
 * After the basis transform, all BSSN quantities are rescaled.
 */
static void BSSN_Cart_to_rescaled_BSSN_rfm(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                                           const BSSN_Cart_basis_struct *restrict BSSN_Cart_basis,
                                           rescaled_BSSN_rfm_basis_struct *restrict rescaled_BSSN_rfm_basis) {
#include "../set_CodeParameters.h"

  REAL xx0, xx1, xx2 __attribute__((unused)); // xx2 might be unused in the case of axisymmetric initial data.
  {
    int unused_Cart_to_i0i1i2[3];
    REAL xx[3];
    Cart_to_xx_and_nearest_i0i1i2(commondata, params, xCart, xx, unused_Cart_to_i0i1i2);
    xx0 = xx[0];
    xx1 = xx[1];
    xx2 = xx[2];
  }
  const REAL tmp0 = cos(xx1);
  const REAL tmp2 = (1.0 / (SINHW));
  const REAL tmp14 = sin(xx1);
  const REAL tmp17 = cos(xx2);
  const REAL tmp20 = sin(xx2);
  const REAL tmp1 = ((tmp0) * (tmp0));
  const REAL tmp9 = exp(tmp2) - exp(-tmp2);
  const REAL tmp15 = tmp0 * tmp14;
  const REAL tmp18 = 2 * tmp17;
  const REAL tmp23 = ((tmp14) * (tmp14));
  const REAL tmp25 = ((tmp17) * (tmp17));
  const REAL tmp27 = ((tmp20) * (tmp20));
  const REAL tmp53 = (1.0 / (tmp14));
  const REAL tmp4 = exp(tmp2 * xx0);
  const REAL tmp5 = exp(-tmp2 * xx0);
  const REAL tmp22 = 2 * BSSN_Cart_basis->AbarDD12 * tmp20;
  const REAL tmp26 = BSSN_Cart_basis->AbarDD00 * tmp25;
  const REAL tmp28 = BSSN_Cart_basis->AbarDD11 * tmp27;
  const REAL tmp29 = BSSN_Cart_basis->AbarDD01 * tmp18;
  const REAL tmp78 = AMPL / tmp9;
  const REAL tmp88 = 2 * BSSN_Cart_basis->gammabarDD12 * tmp20;
  const REAL tmp91 = BSSN_Cart_basis->gammabarDD01 * tmp18;
  const REAL tmp6 = tmp2 * tmp4 + tmp2 * tmp5;
  const REAL tmp11 = ((AMPL) * (AMPL)) / ((tmp9) * (tmp9));
  const REAL tmp31 = ((tmp9) * (tmp9)) / ((AMPL) * (AMPL));
  const REAL tmp33 = tmp4 - tmp5;
  const REAL tmp34 = tmp11 * tmp33;
  const REAL tmp38 = tmp17 * tmp6;
  const REAL tmp40 = tmp20 * tmp6;
  const REAL tmp46 = tmp31 / (tmp33 * tmp6);
  const REAL tmp48 = tmp25 * tmp6;
  const REAL tmp51 = tmp27 * tmp6;
  const REAL tmp55 = ((tmp33) * (tmp33));
  const REAL tmp83 = tmp33 * tmp78;
  const REAL tmp12 = tmp11 * ((tmp6) * (tmp6));
  const REAL tmp32 = tmp31 / ((tmp6) * (tmp6));
  const REAL tmp35 = tmp15 * tmp34;
  const REAL tmp42 = tmp23 * tmp34;
  const REAL tmp56 = tmp11 * tmp55;
  const REAL tmp60 = tmp31 / tmp55;
  const REAL tmp70 = ((AMPL) * (AMPL) * (AMPL)) * tmp55 / ((tmp9) * (tmp9) * (tmp9));
  const REAL tmp16 = tmp12 * tmp15;
  const REAL tmp24 = tmp12 * tmp23;
  const REAL tmp36 = tmp35 * tmp6;
  const REAL tmp39 = tmp1 * tmp34 * tmp38;
  const REAL tmp41 = tmp1 * tmp34 * tmp40;
  const REAL tmp43 = tmp38 * tmp42;
  const REAL tmp44 = tmp40 * tmp42;
  const REAL tmp45 = tmp35 * tmp40;
  const REAL tmp47 = tmp35 * tmp38;
  const REAL tmp49 = tmp42 * tmp48;
  const REAL tmp52 = tmp42 * tmp51;
  const REAL tmp57 = tmp23 * tmp56;
  const REAL tmp58 = tmp15 * tmp56;
  const REAL tmp59 = tmp1 * tmp56;
  const REAL tmp69 = tmp60 / tmp23;
  const REAL tmp71 = ((tmp14) * (tmp14) * (tmp14)) * tmp70;
  const REAL tmp72 = tmp1 * tmp14 * tmp70;
  const REAL tmp50 = tmp17 * tmp44;
  const REAL tmp61 = tmp20 * tmp57;
  const REAL tmp62 = tmp17 * tmp57;
  const REAL tmp63 = tmp25 * tmp58;
  const REAL tmp64 = tmp17 * tmp20 * tmp58;
  const REAL tmp65 = tmp27 * tmp58;
  const REAL tmp73 = (1.0 / (tmp48 * tmp71 + tmp48 * tmp72 + tmp51 * tmp71 + tmp51 * tmp72));
  const REAL tmp82 = -tmp49 - tmp52;
  const REAL tmp84 = tmp73 * (tmp39 + tmp43);
  const REAL tmp85 = tmp73 * (-tmp41 - tmp44);
  rescaled_BSSN_rfm_basis->aDD00 = tmp32 * (BSSN_Cart_basis->AbarDD02 * tmp16 * tmp18 + BSSN_Cart_basis->AbarDD22 * tmp1 * tmp12 + tmp16 * tmp22 +
                                            tmp20 * tmp24 * tmp29 + tmp24 * tmp26 + tmp24 * tmp28);
  rescaled_BSSN_rfm_basis->aDD01 =
      tmp46 * (BSSN_Cart_basis->AbarDD02 * tmp39 - BSSN_Cart_basis->AbarDD02 * tmp43 + BSSN_Cart_basis->AbarDD12 * tmp41 -
               BSSN_Cart_basis->AbarDD12 * tmp44 - BSSN_Cart_basis->AbarDD22 * tmp36 + tmp26 * tmp36 + tmp28 * tmp36 + tmp29 * tmp45);
  rescaled_BSSN_rfm_basis->aDD02 = tmp46 * tmp53 *
                                   (-BSSN_Cart_basis->AbarDD00 * tmp50 + BSSN_Cart_basis->AbarDD01 * tmp49 - BSSN_Cart_basis->AbarDD01 * tmp52 -
                                    BSSN_Cart_basis->AbarDD02 * tmp45 + BSSN_Cart_basis->AbarDD11 * tmp50 + BSSN_Cart_basis->AbarDD12 * tmp47);
  rescaled_BSSN_rfm_basis->aDD11 = tmp60 * (-BSSN_Cart_basis->AbarDD02 * tmp18 * tmp58 + BSSN_Cart_basis->AbarDD22 * tmp57 + tmp20 * tmp29 * tmp59 -
                                            tmp22 * tmp58 + tmp26 * tmp59 + tmp28 * tmp59);
  rescaled_BSSN_rfm_basis->aDD12 = tmp53 * tmp60 *
                                   (-BSSN_Cart_basis->AbarDD00 * tmp64 + BSSN_Cart_basis->AbarDD01 * tmp63 - BSSN_Cart_basis->AbarDD01 * tmp65 +
                                    BSSN_Cart_basis->AbarDD02 * tmp61 + BSSN_Cart_basis->AbarDD11 * tmp64 - BSSN_Cart_basis->AbarDD12 * tmp62);
  rescaled_BSSN_rfm_basis->aDD22 = tmp69 * (BSSN_Cart_basis->AbarDD00 * tmp27 * tmp57 + BSSN_Cart_basis->AbarDD11 * tmp25 * tmp57 - tmp29 * tmp61);
  rescaled_BSSN_rfm_basis->alpha = BSSN_Cart_basis->alpha;
  rescaled_BSSN_rfm_basis->betU0 =
      tmp6 * tmp78 * (BSSN_Cart_basis->BU0 * tmp62 * tmp73 + BSSN_Cart_basis->BU1 * tmp61 * tmp73 + BSSN_Cart_basis->BU2 * tmp73 * (tmp63 + tmp65));
  rescaled_BSSN_rfm_basis->betU1 =
      tmp83 * (BSSN_Cart_basis->BU0 * tmp47 * tmp73 + BSSN_Cart_basis->BU1 * tmp45 * tmp73 + BSSN_Cart_basis->BU2 * tmp73 * tmp82);
  rescaled_BSSN_rfm_basis->betU2 = tmp14 * tmp83 * (BSSN_Cart_basis->BU0 * tmp85 + BSSN_Cart_basis->BU1 * tmp84);
  rescaled_BSSN_rfm_basis->cf = BSSN_Cart_basis->cf;
  rescaled_BSSN_rfm_basis->hDD00 = tmp32 * (BSSN_Cart_basis->gammabarDD00 * tmp24 * tmp25 + BSSN_Cart_basis->gammabarDD02 * tmp16 * tmp18 +
                                            BSSN_Cart_basis->gammabarDD11 * tmp24 * tmp27 + BSSN_Cart_basis->gammabarDD22 * tmp1 * tmp12 - tmp12 +
                                            tmp16 * tmp88 + tmp20 * tmp24 * tmp91);
  rescaled_BSSN_rfm_basis->hDD01 =
      tmp46 * (BSSN_Cart_basis->gammabarDD00 * tmp35 * tmp48 + BSSN_Cart_basis->gammabarDD02 * tmp39 - BSSN_Cart_basis->gammabarDD02 * tmp43 +
               BSSN_Cart_basis->gammabarDD11 * tmp35 * tmp51 + BSSN_Cart_basis->gammabarDD12 * tmp41 - BSSN_Cart_basis->gammabarDD12 * tmp44 -
               BSSN_Cart_basis->gammabarDD22 * tmp36 + tmp45 * tmp91);
  rescaled_BSSN_rfm_basis->hDD02 =
      tmp46 * tmp53 *
      (-BSSN_Cart_basis->gammabarDD00 * tmp50 + BSSN_Cart_basis->gammabarDD01 * tmp49 - BSSN_Cart_basis->gammabarDD01 * tmp52 -
       BSSN_Cart_basis->gammabarDD02 * tmp45 + BSSN_Cart_basis->gammabarDD11 * tmp50 + BSSN_Cart_basis->gammabarDD12 * tmp47);
  rescaled_BSSN_rfm_basis->hDD11 =
      tmp60 * (BSSN_Cart_basis->gammabarDD00 * tmp25 * tmp59 - BSSN_Cart_basis->gammabarDD02 * tmp18 * tmp58 +
               BSSN_Cart_basis->gammabarDD11 * tmp27 * tmp59 + BSSN_Cart_basis->gammabarDD22 * tmp57 + tmp20 * tmp59 * tmp91 - tmp56 - tmp58 * tmp88);
  rescaled_BSSN_rfm_basis->hDD12 =
      tmp53 * tmp60 *
      (-BSSN_Cart_basis->gammabarDD00 * tmp64 + BSSN_Cart_basis->gammabarDD01 * tmp63 - BSSN_Cart_basis->gammabarDD01 * tmp65 +
       BSSN_Cart_basis->gammabarDD02 * tmp61 + BSSN_Cart_basis->gammabarDD11 * tmp64 - BSSN_Cart_basis->gammabarDD12 * tmp62);
  rescaled_BSSN_rfm_basis->hDD22 =
      tmp69 * (BSSN_Cart_basis->gammabarDD00 * tmp27 * tmp57 + BSSN_Cart_basis->gammabarDD11 * tmp25 * tmp57 - tmp57 - tmp61 * tmp91);
  rescaled_BSSN_rfm_basis->trK = BSSN_Cart_basis->trK;
  rescaled_BSSN_rfm_basis->vetU0 =
      tmp6 * tmp78 *
      (BSSN_Cart_basis->betaU0 * tmp62 * tmp73 + BSSN_Cart_basis->betaU1 * tmp61 * tmp73 + BSSN_Cart_basis->betaU2 * tmp73 * (tmp63 + tmp65));
  rescaled_BSSN_rfm_basis->vetU1 =
      tmp83 * (BSSN_Cart_basis->betaU0 * tmp47 * tmp73 + BSSN_Cart_basis->betaU1 * tmp45 * tmp73 + BSSN_Cart_basis->betaU2 * tmp73 * tmp82);
  rescaled_BSSN_rfm_basis->vetU2 = tmp14 * tmp83 * (BSSN_Cart_basis->betaU0 * tmp85 + BSSN_Cart_basis->betaU1 * tmp84);
}
/*
 * Compute lambdaU in SinhSpherical coordinates
 */
static void initial_data_lambdaU_grid_interior(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                               REAL *restrict xx[3], REAL *restrict in_gfs) {
#include "../set_CodeParameters.h"
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
      const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
        const REAL xx0 = xx[0][i0]; /*
                                     * NRPy+-Generated GF Access/FD Code, Step 1 of 2:
                                     * Read gridfunction(s) from main memory and compute FD stencils as needed.
                                     */
        const REAL hDD00_i1m4 = in_gfs[IDX4(HDD00GF, i0, i1 - 4, i2)];
        const REAL hDD00_i1m3 = in_gfs[IDX4(HDD00GF, i0, i1 - 3, i2)];
        const REAL hDD00_i1m2 = in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2)];
        const REAL hDD00_i1m1 = in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2)];
        const REAL hDD00_i0m4 = in_gfs[IDX4(HDD00GF, i0 - 4, i1, i2)];
        const REAL hDD00_i0m3 = in_gfs[IDX4(HDD00GF, i0 - 3, i1, i2)];
        const REAL hDD00_i0m2 = in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2)];
        const REAL hDD00_i0m1 = in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2)];
        const REAL hDD00 = in_gfs[IDX4(HDD00GF, i0, i1, i2)];
        const REAL hDD00_i0p1 = in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2)];
        const REAL hDD00_i0p2 = in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2)];
        const REAL hDD00_i0p3 = in_gfs[IDX4(HDD00GF, i0 + 3, i1, i2)];
        const REAL hDD00_i0p4 = in_gfs[IDX4(HDD00GF, i0 + 4, i1, i2)];
        const REAL hDD00_i1p1 = in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2)];
        const REAL hDD00_i1p2 = in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2)];
        const REAL hDD00_i1p3 = in_gfs[IDX4(HDD00GF, i0, i1 + 3, i2)];
        const REAL hDD00_i1p4 = in_gfs[IDX4(HDD00GF, i0, i1 + 4, i2)];
        const REAL hDD01_i1m4 = in_gfs[IDX4(HDD01GF, i0, i1 - 4, i2)];
        const REAL hDD01_i1m3 = in_gfs[IDX4(HDD01GF, i0, i1 - 3, i2)];
        const REAL hDD01_i1m2 = in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2)];
        const REAL hDD01_i1m1 = in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2)];
        const REAL hDD01_i0m4 = in_gfs[IDX4(HDD01GF, i0 - 4, i1, i2)];
        const REAL hDD01_i0m3 = in_gfs[IDX4(HDD01GF, i0 - 3, i1, i2)];
        const REAL hDD01_i0m2 = in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2)];
        const REAL hDD01_i0m1 = in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2)];
        const REAL hDD01 = in_gfs[IDX4(HDD01GF, i0, i1, i2)];
        const REAL hDD01_i0p1 = in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2)];
        const REAL hDD01_i0p2 = in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2)];
        const REAL hDD01_i0p3 = in_gfs[IDX4(HDD01GF, i0 + 3, i1, i2)];
        const REAL hDD01_i0p4 = in_gfs[IDX4(HDD01GF, i0 + 4, i1, i2)];
        const REAL hDD01_i1p1 = in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2)];
        const REAL hDD01_i1p2 = in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2)];
        const REAL hDD01_i1p3 = in_gfs[IDX4(HDD01GF, i0, i1 + 3, i2)];
        const REAL hDD01_i1p4 = in_gfs[IDX4(HDD01GF, i0, i1 + 4, i2)];
        const REAL hDD02_i1m4 = in_gfs[IDX4(HDD02GF, i0, i1 - 4, i2)];
        const REAL hDD02_i1m3 = in_gfs[IDX4(HDD02GF, i0, i1 - 3, i2)];
        const REAL hDD02_i1m2 = in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2)];
        const REAL hDD02_i1m1 = in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2)];
        const REAL hDD02_i0m4 = in_gfs[IDX4(HDD02GF, i0 - 4, i1, i2)];
        const REAL hDD02_i0m3 = in_gfs[IDX4(HDD02GF, i0 - 3, i1, i2)];
        const REAL hDD02_i0m2 = in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2)];
        const REAL hDD02_i0m1 = in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2)];
        const REAL hDD02 = in_gfs[IDX4(HDD02GF, i0, i1, i2)];
        const REAL hDD02_i0p1 = in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2)];
        const REAL hDD02_i0p2 = in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2)];
        const REAL hDD02_i0p3 = in_gfs[IDX4(HDD02GF, i0 + 3, i1, i2)];
        const REAL hDD02_i0p4 = in_gfs[IDX4(HDD02GF, i0 + 4, i1, i2)];
        const REAL hDD02_i1p1 = in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2)];
        const REAL hDD02_i1p2 = in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2)];
        const REAL hDD02_i1p3 = in_gfs[IDX4(HDD02GF, i0, i1 + 3, i2)];
        const REAL hDD02_i1p4 = in_gfs[IDX4(HDD02GF, i0, i1 + 4, i2)];
        const REAL hDD11_i1m4 = in_gfs[IDX4(HDD11GF, i0, i1 - 4, i2)];
        const REAL hDD11_i1m3 = in_gfs[IDX4(HDD11GF, i0, i1 - 3, i2)];
        const REAL hDD11_i1m2 = in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2)];
        const REAL hDD11_i1m1 = in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2)];
        const REAL hDD11_i0m4 = in_gfs[IDX4(HDD11GF, i0 - 4, i1, i2)];
        const REAL hDD11_i0m3 = in_gfs[IDX4(HDD11GF, i0 - 3, i1, i2)];
        const REAL hDD11_i0m2 = in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2)];
        const REAL hDD11_i0m1 = in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2)];
        const REAL hDD11 = in_gfs[IDX4(HDD11GF, i0, i1, i2)];
        const REAL hDD11_i0p1 = in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2)];
        const REAL hDD11_i0p2 = in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2)];
        const REAL hDD11_i0p3 = in_gfs[IDX4(HDD11GF, i0 + 3, i1, i2)];
        const REAL hDD11_i0p4 = in_gfs[IDX4(HDD11GF, i0 + 4, i1, i2)];
        const REAL hDD11_i1p1 = in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2)];
        const REAL hDD11_i1p2 = in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2)];
        const REAL hDD11_i1p3 = in_gfs[IDX4(HDD11GF, i0, i1 + 3, i2)];
        const REAL hDD11_i1p4 = in_gfs[IDX4(HDD11GF, i0, i1 + 4, i2)];
        const REAL hDD12_i1m4 = in_gfs[IDX4(HDD12GF, i0, i1 - 4, i2)];
        const REAL hDD12_i1m3 = in_gfs[IDX4(HDD12GF, i0, i1 - 3, i2)];
        const REAL hDD12_i1m2 = in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2)];
        const REAL hDD12_i1m1 = in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2)];
        const REAL hDD12_i0m4 = in_gfs[IDX4(HDD12GF, i0 - 4, i1, i2)];
        const REAL hDD12_i0m3 = in_gfs[IDX4(HDD12GF, i0 - 3, i1, i2)];
        const REAL hDD12_i0m2 = in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2)];
        const REAL hDD12_i0m1 = in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2)];
        const REAL hDD12 = in_gfs[IDX4(HDD12GF, i0, i1, i2)];
        const REAL hDD12_i0p1 = in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2)];
        const REAL hDD12_i0p2 = in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2)];
        const REAL hDD12_i0p3 = in_gfs[IDX4(HDD12GF, i0 + 3, i1, i2)];
        const REAL hDD12_i0p4 = in_gfs[IDX4(HDD12GF, i0 + 4, i1, i2)];
        const REAL hDD12_i1p1 = in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2)];
        const REAL hDD12_i1p2 = in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2)];
        const REAL hDD12_i1p3 = in_gfs[IDX4(HDD12GF, i0, i1 + 3, i2)];
        const REAL hDD12_i1p4 = in_gfs[IDX4(HDD12GF, i0, i1 + 4, i2)];
        const REAL hDD22_i1m4 = in_gfs[IDX4(HDD22GF, i0, i1 - 4, i2)];
        const REAL hDD22_i1m3 = in_gfs[IDX4(HDD22GF, i0, i1 - 3, i2)];
        const REAL hDD22_i1m2 = in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2)];
        const REAL hDD22_i1m1 = in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2)];
        const REAL hDD22_i0m4 = in_gfs[IDX4(HDD22GF, i0 - 4, i1, i2)];
        const REAL hDD22_i0m3 = in_gfs[IDX4(HDD22GF, i0 - 3, i1, i2)];
        const REAL hDD22_i0m2 = in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2)];
        const REAL hDD22_i0m1 = in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2)];
        const REAL hDD22 = in_gfs[IDX4(HDD22GF, i0, i1, i2)];
        const REAL hDD22_i0p1 = in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2)];
        const REAL hDD22_i0p2 = in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2)];
        const REAL hDD22_i0p3 = in_gfs[IDX4(HDD22GF, i0 + 3, i1, i2)];
        const REAL hDD22_i0p4 = in_gfs[IDX4(HDD22GF, i0 + 4, i1, i2)];
        const REAL hDD22_i1p1 = in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2)];
        const REAL hDD22_i1p2 = in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2)];
        const REAL hDD22_i1p3 = in_gfs[IDX4(HDD22GF, i0, i1 + 3, i2)];
        const REAL hDD22_i1p4 = in_gfs[IDX4(HDD22GF, i0, i1 + 4, i2)];
        const REAL FDPart1_Rational_4_5 = 4.0 / 5.0;
        const REAL FDPart1_Rational_4_105 = 4.0 / 105.0;
        const REAL FDPart1_Rational_1_5 = 1.0 / 5.0;
        const REAL FDPart1_Rational_1_280 = 1.0 / 280.0;
        const REAL hDD_dD000 = invdxx0 * (FDPart1_Rational_1_280 * (hDD00_i0m4 - hDD00_i0p4) + FDPart1_Rational_1_5 * (hDD00_i0m2 - hDD00_i0p2) +
                                          FDPart1_Rational_4_105 * (-hDD00_i0m3 + hDD00_i0p3) + FDPart1_Rational_4_5 * (-hDD00_i0m1 + hDD00_i0p1));
        const REAL hDD_dD001 = invdxx1 * (FDPart1_Rational_1_280 * (hDD00_i1m4 - hDD00_i1p4) + FDPart1_Rational_1_5 * (hDD00_i1m2 - hDD00_i1p2) +
                                          FDPart1_Rational_4_105 * (-hDD00_i1m3 + hDD00_i1p3) + FDPart1_Rational_4_5 * (-hDD00_i1m1 + hDD00_i1p1));
        const REAL hDD_dD010 = invdxx0 * (FDPart1_Rational_1_280 * (hDD01_i0m4 - hDD01_i0p4) + FDPart1_Rational_1_5 * (hDD01_i0m2 - hDD01_i0p2) +
                                          FDPart1_Rational_4_105 * (-hDD01_i0m3 + hDD01_i0p3) + FDPart1_Rational_4_5 * (-hDD01_i0m1 + hDD01_i0p1));
        const REAL hDD_dD011 = invdxx1 * (FDPart1_Rational_1_280 * (hDD01_i1m4 - hDD01_i1p4) + FDPart1_Rational_1_5 * (hDD01_i1m2 - hDD01_i1p2) +
                                          FDPart1_Rational_4_105 * (-hDD01_i1m3 + hDD01_i1p3) + FDPart1_Rational_4_5 * (-hDD01_i1m1 + hDD01_i1p1));
        const REAL hDD_dD020 = invdxx0 * (FDPart1_Rational_1_280 * (hDD02_i0m4 - hDD02_i0p4) + FDPart1_Rational_1_5 * (hDD02_i0m2 - hDD02_i0p2) +
                                          FDPart1_Rational_4_105 * (-hDD02_i0m3 + hDD02_i0p3) + FDPart1_Rational_4_5 * (-hDD02_i0m1 + hDD02_i0p1));
        const REAL hDD_dD021 = invdxx1 * (FDPart1_Rational_1_280 * (hDD02_i1m4 - hDD02_i1p4) + FDPart1_Rational_1_5 * (hDD02_i1m2 - hDD02_i1p2) +
                                          FDPart1_Rational_4_105 * (-hDD02_i1m3 + hDD02_i1p3) + FDPart1_Rational_4_5 * (-hDD02_i1m1 + hDD02_i1p1));
        const REAL hDD_dD110 = invdxx0 * (FDPart1_Rational_1_280 * (hDD11_i0m4 - hDD11_i0p4) + FDPart1_Rational_1_5 * (hDD11_i0m2 - hDD11_i0p2) +
                                          FDPart1_Rational_4_105 * (-hDD11_i0m3 + hDD11_i0p3) + FDPart1_Rational_4_5 * (-hDD11_i0m1 + hDD11_i0p1));
        const REAL hDD_dD111 = invdxx1 * (FDPart1_Rational_1_280 * (hDD11_i1m4 - hDD11_i1p4) + FDPart1_Rational_1_5 * (hDD11_i1m2 - hDD11_i1p2) +
                                          FDPart1_Rational_4_105 * (-hDD11_i1m3 + hDD11_i1p3) + FDPart1_Rational_4_5 * (-hDD11_i1m1 + hDD11_i1p1));
        const REAL hDD_dD120 = invdxx0 * (FDPart1_Rational_1_280 * (hDD12_i0m4 - hDD12_i0p4) + FDPart1_Rational_1_5 * (hDD12_i0m2 - hDD12_i0p2) +
                                          FDPart1_Rational_4_105 * (-hDD12_i0m3 + hDD12_i0p3) + FDPart1_Rational_4_5 * (-hDD12_i0m1 + hDD12_i0p1));
        const REAL hDD_dD121 = invdxx1 * (FDPart1_Rational_1_280 * (hDD12_i1m4 - hDD12_i1p4) + FDPart1_Rational_1_5 * (hDD12_i1m2 - hDD12_i1p2) +
                                          FDPart1_Rational_4_105 * (-hDD12_i1m3 + hDD12_i1p3) + FDPart1_Rational_4_5 * (-hDD12_i1m1 + hDD12_i1p1));
        const REAL hDD_dD220 = invdxx0 * (FDPart1_Rational_1_280 * (hDD22_i0m4 - hDD22_i0p4) + FDPart1_Rational_1_5 * (hDD22_i0m2 - hDD22_i0p2) +
                                          FDPart1_Rational_4_105 * (-hDD22_i0m3 + hDD22_i0p3) + FDPart1_Rational_4_5 * (-hDD22_i0m1 + hDD22_i0p1));
        const REAL hDD_dD221 = invdxx1 * (FDPart1_Rational_1_280 * (hDD22_i1m4 - hDD22_i1p4) + FDPart1_Rational_1_5 * (hDD22_i1m2 - hDD22_i1p2) +
                                          FDPart1_Rational_4_105 * (-hDD22_i1m3 + hDD22_i1p3) + FDPart1_Rational_4_5 * (-hDD22_i1m1 + hDD22_i1p1));

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        const REAL FDPart3tmp0 = (1.0 / (SINHW));
        const REAL FDPart3tmp7 = sin(xx1);
        const REAL FDPart3tmp16 = cos(xx1);
        const REAL FDPart3tmp85 = (1.0 / ((SINHW) * (SINHW)));
        const REAL FDPart3tmp8 = ((FDPart3tmp7) * (FDPart3tmp7));
        const REAL FDPart3tmp11 = exp(FDPart3tmp0) - exp(-FDPart3tmp0);
        const REAL FDPart3tmp2 = exp(FDPart3tmp0 * xx0);
        const REAL FDPart3tmp4 = exp(-FDPart3tmp0 * xx0);
        const REAL FDPart3tmp12 = ((AMPL) * (AMPL)) / ((FDPart3tmp11) * (FDPart3tmp11));
        const REAL FDPart3tmp22 = ((AMPL) * (AMPL) * (AMPL) * (AMPL)) / ((FDPart3tmp11) * (FDPart3tmp11) * (FDPart3tmp11) * (FDPart3tmp11));
        const REAL FDPart3tmp94 = AMPL / FDPart3tmp11;
        const REAL FDPart3tmp9 = FDPart3tmp2 - FDPart3tmp4;
        const REAL FDPart3tmp6 = FDPart3tmp0 * FDPart3tmp2 + FDPart3tmp0 * FDPart3tmp4;
        const REAL FDPart3tmp31 = ((FDPart3tmp9) * (FDPart3tmp9) * (FDPart3tmp9) * (FDPart3tmp9));
        const REAL FDPart3tmp47 = 2 * FDPart3tmp0 * FDPart3tmp2 + 2 * FDPart3tmp0 * FDPart3tmp4;
        const REAL FDPart3tmp88 = 2 * FDPart3tmp2 * FDPart3tmp85 - 2 * FDPart3tmp4 * FDPart3tmp85;
        const REAL FDPart3tmp13 = FDPart3tmp12 * ((FDPart3tmp9) * (FDPart3tmp9));
        const REAL FDPart3tmp24 = FDPart3tmp6 * ((FDPart3tmp9) * (FDPart3tmp9) * (FDPart3tmp9));
        const REAL FDPart3tmp30 = ((FDPart3tmp6) * (FDPart3tmp6));
        const REAL FDPart3tmp36 = FDPart3tmp22 * FDPart3tmp31 * FDPart3tmp8 * ((hDD12) * (hDD12));
        const REAL FDPart3tmp49 = FDPart3tmp12 * FDPart3tmp47 * FDPart3tmp9;
        const REAL FDPart3tmp89 = FDPart3tmp12 * FDPart3tmp9 * (FDPart3tmp2 * FDPart3tmp85 - FDPart3tmp4 * FDPart3tmp85);
        const REAL FDPart3tmp14 = FDPart3tmp13 * FDPart3tmp8;
        const REAL FDPart3tmp25 = FDPart3tmp13 * hDD11 + FDPart3tmp13;
        const REAL FDPart3tmp27 = FDPart3tmp12 * FDPart3tmp6 * FDPart3tmp9;
        const REAL FDPart3tmp33 = FDPart3tmp12 * FDPart3tmp30;
        const REAL FDPart3tmp37 = FDPart3tmp22 * FDPart3tmp30 * ((FDPart3tmp9) * (FDPart3tmp9));
        const REAL FDPart3tmp45 = FDPart3tmp13 * FDPart3tmp7;
        const REAL FDPart3tmp50 = FDPart3tmp49 * FDPart3tmp7 * hDD12;
        const REAL FDPart3tmp61 = FDPart3tmp49 * FDPart3tmp8;
        const REAL FDPart3tmp70 = (1.0 / 2.0) * FDPart3tmp47 * FDPart3tmp9 / FDPart3tmp30;
        const REAL FDPart3tmp80 = FDPart3tmp13 * hDD_dD111;
        const REAL FDPart3tmp91 = FDPart3tmp12 * FDPart3tmp6 * FDPart3tmp88;
        const REAL FDPart3tmp97 = -1.0 / 2.0 * FDPart3tmp47 / FDPart3tmp9;
        const REAL FDPart3tmp19 = 2 * FDPart3tmp13 * FDPart3tmp16 * FDPart3tmp7;
        const REAL FDPart3tmp34 = FDPart3tmp33 * hDD00 + FDPart3tmp33;
        const REAL FDPart3tmp38 = FDPart3tmp37 * FDPart3tmp8 * ((hDD02) * (hDD02));
        const REAL FDPart3tmp39 = FDPart3tmp14 * hDD22 + FDPart3tmp14;
        const REAL FDPart3tmp40 = FDPart3tmp37 * ((hDD01) * (hDD01));
        const REAL FDPart3tmp52 = FDPart3tmp27 * FDPart3tmp7 * hDD_dD021;
        const REAL FDPart3tmp74 = FDPart3tmp33 * hDD_dD001;
        const REAL FDPart3tmp77 = FDPart3tmp13 * hDD_dD110 + FDPart3tmp49 * hDD11 + FDPart3tmp49;
        const REAL FDPart3tmp81 = 2 * FDPart3tmp13 * FDPart3tmp16 * hDD12 + 2 * FDPart3tmp45 * hDD_dD121;
        const REAL FDPart3tmp82 = 2 * FDPart3tmp27;
        const REAL FDPart3tmp92 = FDPart3tmp33 * hDD_dD000 + FDPart3tmp91 * hDD00 + FDPart3tmp91;
        const REAL FDPart3tmp29 = FDPart3tmp22 * FDPart3tmp24 * FDPart3tmp7 * hDD01 * hDD12 - FDPart3tmp25 * FDPart3tmp27 * FDPart3tmp7 * hDD02;
        const REAL FDPart3tmp51 = FDPart3tmp16 * FDPart3tmp27 * hDD02;
        const REAL FDPart3tmp57 = -FDPart3tmp34 * FDPart3tmp45 * hDD12 + FDPart3tmp37 * FDPart3tmp7 * hDD01 * hDD02;
        const REAL FDPart3tmp63 = FDPart3tmp14 * hDD_dD220 + FDPart3tmp61 * hDD22 + FDPart3tmp61;
        const REAL FDPart3tmp64 = FDPart3tmp45 * hDD_dD120 + FDPart3tmp50;
        const REAL FDPart3tmp66 = FDPart3tmp22 * FDPart3tmp24 * FDPart3tmp8 * hDD02 * hDD12 - FDPart3tmp27 * FDPart3tmp39 * hDD01;
        const REAL FDPart3tmp72 = -FDPart3tmp14 * hDD_dD220 - FDPart3tmp61 * hDD22 - FDPart3tmp61;
        const REAL FDPart3tmp83 = -FDPart3tmp13 * hDD_dD110 - FDPart3tmp49 * hDD11 - FDPart3tmp49 + FDPart3tmp82 * hDD_dD011;
        const REAL FDPart3tmp90 = FDPart3tmp7 * FDPart3tmp82 * hDD_dD020 + 2 * hDD02 * (FDPart3tmp33 * FDPart3tmp7 + FDPart3tmp7 * FDPart3tmp89);
        const REAL FDPart3tmp93 = -FDPart3tmp74 + FDPart3tmp82 * hDD_dD010 + 2 * hDD01 * (FDPart3tmp33 + FDPart3tmp89);
        const REAL FDPart3tmp21 = FDPart3tmp14 * hDD_dD221 + FDPart3tmp19 * hDD22 + FDPart3tmp19;
        const REAL FDPart3tmp42 =
            (1.0 /
             (2 * pow(AMPL, 6) * FDPart3tmp30 * FDPart3tmp31 * FDPart3tmp8 * hDD01 * hDD02 * hDD12 / pow(FDPart3tmp11, 6) +
              FDPart3tmp25 * FDPart3tmp34 * FDPart3tmp39 - FDPart3tmp25 * FDPart3tmp38 - FDPart3tmp34 * FDPart3tmp36 - FDPart3tmp39 * FDPart3tmp40));
        const REAL FDPart3tmp65 = -FDPart3tmp51 - FDPart3tmp52 + FDPart3tmp64;
        const REAL FDPart3tmp71 = -FDPart3tmp14 * hDD_dD221 - FDPart3tmp19 * hDD22 - FDPart3tmp19;
        const REAL FDPart3tmp43 = (1.0 / 2.0) * FDPart3tmp42;
        const REAL FDPart3tmp54 = -FDPart3tmp45 * hDD_dD120 - FDPart3tmp50 + FDPart3tmp51 + FDPart3tmp52;
        const REAL FDPart3tmp55 = FDPart3tmp42 * (FDPart3tmp25 * FDPart3tmp39 - FDPart3tmp36);
        const REAL FDPart3tmp58 = 2 * FDPart3tmp42;
        const REAL FDPart3tmp73 = FDPart3tmp42 * (FDPart3tmp25 * FDPart3tmp34 - FDPart3tmp40);
        const REAL FDPart3tmp78 = FDPart3tmp51 + FDPart3tmp52 + FDPart3tmp64;
        const REAL FDPart3tmp84 = FDPart3tmp42 * (FDPart3tmp34 * FDPart3tmp39 - FDPart3tmp38);
        const REAL FDPart3tmp44 = FDPart3tmp29 * FDPart3tmp43;
        const REAL FDPart3tmp56 = (1.0 / 2.0) * FDPart3tmp55;
        const REAL FDPart3tmp59 = FDPart3tmp57 * FDPart3tmp58;
        const REAL FDPart3tmp67 = FDPart3tmp43 * FDPart3tmp66;
        const REAL FDPart3tmp68 = FDPart3tmp29 * FDPart3tmp58;
        const REAL FDPart3tmp79 = FDPart3tmp58 * FDPart3tmp66;
        const REAL FDPart3tmp95 = FDPart3tmp43 * FDPart3tmp57;
        const REAL FDPart3tmp96 = (1.0 / 2.0) * FDPart3tmp84;
        const REAL FDPart3tmp99 = (1.0 / 2.0) * FDPart3tmp73;
        in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2)] =
            FDPart3tmp6 * FDPart3tmp94 *
            (FDPart3tmp55 *
                 (FDPart3tmp44 * FDPart3tmp90 + FDPart3tmp56 * FDPart3tmp92 + FDPart3tmp67 * FDPart3tmp93 - 1.0 / 2.0 * FDPart3tmp88 / FDPart3tmp6) +
             FDPart3tmp59 * (FDPart3tmp21 * FDPart3tmp44 + FDPart3tmp54 * FDPart3tmp56) +
             FDPart3tmp68 * (FDPart3tmp44 * FDPart3tmp63 + FDPart3tmp65 * FDPart3tmp67) +
             FDPart3tmp73 * (FDPart3tmp56 * FDPart3tmp72 + FDPart3tmp67 * FDPart3tmp71 + FDPart3tmp70 * FDPart3tmp8) +
             FDPart3tmp79 * (FDPart3tmp44 * FDPart3tmp78 + FDPart3tmp56 * FDPart3tmp74 + FDPart3tmp67 * FDPart3tmp77) +
             FDPart3tmp84 * (FDPart3tmp44 * FDPart3tmp81 + FDPart3tmp56 * FDPart3tmp83 + FDPart3tmp67 * FDPart3tmp80 + FDPart3tmp70));
        in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2)] =
            FDPart3tmp9 * FDPart3tmp94 *
            (FDPart3tmp55 * (FDPart3tmp67 * FDPart3tmp92 + FDPart3tmp90 * FDPart3tmp95 + FDPart3tmp93 * FDPart3tmp96) +
             FDPart3tmp59 * (FDPart3tmp21 * FDPart3tmp95 + FDPart3tmp54 * FDPart3tmp67) +
             FDPart3tmp68 * (FDPart3tmp63 * FDPart3tmp95 + FDPart3tmp65 * FDPart3tmp96) +
             FDPart3tmp73 * (FDPart3tmp16 * FDPart3tmp7 + FDPart3tmp67 * FDPart3tmp72 + FDPart3tmp71 * FDPart3tmp96) +
             FDPart3tmp79 * (FDPart3tmp67 * FDPart3tmp74 + FDPart3tmp77 * FDPart3tmp96 + FDPart3tmp78 * FDPart3tmp95 + FDPart3tmp97) +
             FDPart3tmp84 * (FDPart3tmp67 * FDPart3tmp83 + FDPart3tmp80 * FDPart3tmp96 + FDPart3tmp81 * FDPart3tmp95));
        in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2)] =
            FDPart3tmp7 * FDPart3tmp9 * FDPart3tmp94 *
            (FDPart3tmp55 * (FDPart3tmp44 * FDPart3tmp92 + FDPart3tmp90 * FDPart3tmp99 + FDPart3tmp93 * FDPart3tmp95) +
             FDPart3tmp59 * (-FDPart3tmp16 / FDPart3tmp7 + FDPart3tmp21 * FDPart3tmp99 + FDPart3tmp44 * FDPart3tmp54) +
             FDPart3tmp68 * (FDPart3tmp63 * FDPart3tmp99 + FDPart3tmp65 * FDPart3tmp95 + FDPart3tmp97) +
             FDPart3tmp73 * (FDPart3tmp44 * FDPart3tmp72 + FDPart3tmp71 * FDPart3tmp95) +
             FDPart3tmp79 * (FDPart3tmp44 * FDPart3tmp74 + FDPart3tmp77 * FDPart3tmp95 + FDPart3tmp78 * FDPart3tmp99) +
             FDPart3tmp84 * (FDPart3tmp44 * FDPart3tmp83 + FDPart3tmp80 * FDPart3tmp95 + FDPart3tmp81 * FDPart3tmp99));

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}

/*
 * Read ADM data in the Cartesian basis, and output rescaled BSSN data in the SinhSpherical basis
 */
void initial_data_reader__convert_ADM_Cartesian_to_BSSN__rfm__SinhSpherical(
    const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct,
    MoL_gridfunctions_struct *restrict gridfuncs, ID_persist_struct *restrict ID_persist,
    void ID_function(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                     const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data)) {

  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  LOOP_OMP("omp parallel for", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
    // xCart is the global Cartesian coordinate, which accounts for any grid offsets from the origin.
    REAL xCart[3];
    xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

    // Read or compute initial data at destination point xCart
    initial_data_struct initial_data;
    ID_function(commondata, params, xCart, ID_persist, &initial_data);

    ADM_Cart_basis_struct ADM_Cart_basis;
    ADM_SphorCart_to_Cart(commondata, params, xCart, &initial_data, &ADM_Cart_basis);

    BSSN_Cart_basis_struct BSSN_Cart_basis;
    ADM_Cart_to_BSSN_Cart(commondata, params, xCart, &ADM_Cart_basis, &BSSN_Cart_basis);

    rescaled_BSSN_rfm_basis_struct rescaled_BSSN_rfm_basis;
    BSSN_Cart_to_rescaled_BSSN_rfm(commondata, params, xCart, &BSSN_Cart_basis, &rescaled_BSSN_rfm_basis);

    const int idx3 = IDX3(i0, i1, i2);
    gridfuncs->y_n_gfs[IDX4pt(ADD00GF, idx3)] = rescaled_BSSN_rfm_basis.aDD00;
    gridfuncs->y_n_gfs[IDX4pt(ADD01GF, idx3)] = rescaled_BSSN_rfm_basis.aDD01;
    gridfuncs->y_n_gfs[IDX4pt(ADD02GF, idx3)] = rescaled_BSSN_rfm_basis.aDD02;
    gridfuncs->y_n_gfs[IDX4pt(ADD11GF, idx3)] = rescaled_BSSN_rfm_basis.aDD11;
    gridfuncs->y_n_gfs[IDX4pt(ADD12GF, idx3)] = rescaled_BSSN_rfm_basis.aDD12;
    gridfuncs->y_n_gfs[IDX4pt(ADD22GF, idx3)] = rescaled_BSSN_rfm_basis.aDD22;
    gridfuncs->y_n_gfs[IDX4pt(ALPHAGF, idx3)] = rescaled_BSSN_rfm_basis.alpha;
    gridfuncs->y_n_gfs[IDX4pt(BETU0GF, idx3)] = rescaled_BSSN_rfm_basis.betU0;
    gridfuncs->y_n_gfs[IDX4pt(BETU1GF, idx3)] = rescaled_BSSN_rfm_basis.betU1;
    gridfuncs->y_n_gfs[IDX4pt(BETU2GF, idx3)] = rescaled_BSSN_rfm_basis.betU2;
    gridfuncs->y_n_gfs[IDX4pt(CFGF, idx3)] = rescaled_BSSN_rfm_basis.cf;
    gridfuncs->y_n_gfs[IDX4pt(HDD00GF, idx3)] = rescaled_BSSN_rfm_basis.hDD00;
    gridfuncs->y_n_gfs[IDX4pt(HDD01GF, idx3)] = rescaled_BSSN_rfm_basis.hDD01;
    gridfuncs->y_n_gfs[IDX4pt(HDD02GF, idx3)] = rescaled_BSSN_rfm_basis.hDD02;
    gridfuncs->y_n_gfs[IDX4pt(HDD11GF, idx3)] = rescaled_BSSN_rfm_basis.hDD11;
    gridfuncs->y_n_gfs[IDX4pt(HDD12GF, idx3)] = rescaled_BSSN_rfm_basis.hDD12;
    gridfuncs->y_n_gfs[IDX4pt(HDD22GF, idx3)] = rescaled_BSSN_rfm_basis.hDD22;
    gridfuncs->y_n_gfs[IDX4pt(TRKGF, idx3)] = rescaled_BSSN_rfm_basis.trK;
    gridfuncs->y_n_gfs[IDX4pt(VETU0GF, idx3)] = rescaled_BSSN_rfm_basis.vetU0;
    gridfuncs->y_n_gfs[IDX4pt(VETU1GF, idx3)] = rescaled_BSSN_rfm_basis.vetU1;
    gridfuncs->y_n_gfs[IDX4pt(VETU2GF, idx3)] = rescaled_BSSN_rfm_basis.vetU2;

    // Initialize lambdaU to zero
    gridfuncs->y_n_gfs[IDX4pt(LAMBDAU0GF, idx3)] = 0.0;
    gridfuncs->y_n_gfs[IDX4pt(LAMBDAU1GF, idx3)] = 0.0;
    gridfuncs->y_n_gfs[IDX4pt(LAMBDAU2GF, idx3)] = 0.0;
  } // END LOOP over all gridpoints on given grid

  // Now we've set all but lambda^i, which will be computed via a finite-difference of hDD.
  //    However, hDD is not correctly set in inner boundary points so we apply inner bcs first.

  // Apply inner bcs to get correct values of all tensor quantities across symmetry boundaries;
  //    BSSN_Cart_to_rescaled_BSSN_rfm() converts each xCart->xx, which guarantees a mapping
  //    to the grid interior. It therefore does not account for parity conditions across
  //    symmetry boundaries being correct.
  apply_bcs_inner_only(commondata, params, bcstruct, gridfuncs->y_n_gfs);

  initial_data_lambdaU_grid_interior(commondata, params, xx, gridfuncs->y_n_gfs);
}
