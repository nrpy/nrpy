#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/*
 * Finite difference function for operator dD0, with FD accuracy order 8.
 */
static REAL fd_function_dD0_fdorder8(const REAL FDPROTO_i0m1, const REAL FDPROTO_i0m2, const REAL FDPROTO_i0m3, const REAL FDPROTO_i0m4,
                                     const REAL FDPROTO_i0p1, const REAL FDPROTO_i0p2, const REAL FDPROTO_i0p3, const REAL FDPROTO_i0p4,
                                     const REAL invdxx0) {
  const REAL FDPart1_Rational_4_5 = 4.0 / 5.0;
  const REAL FDPart1_Rational_4_105 = 4.0 / 105.0;
  const REAL FDPart1_Rational_1_5 = 1.0 / 5.0;
  const REAL FDPart1_Rational_1_280 = 1.0 / 280.0;
  const REAL FD_result = invdxx0 * (FDPart1_Rational_1_280 * (FDPROTO_i0m4 - FDPROTO_i0p4) + FDPart1_Rational_1_5 * (FDPROTO_i0m2 - FDPROTO_i0p2) +
                                    FDPart1_Rational_4_105 * (-FDPROTO_i0m3 + FDPROTO_i0p3) + FDPart1_Rational_4_5 * (-FDPROTO_i0m1 + FDPROTO_i0p1));

  return FD_result;
}
/*
 * Finite difference function for operator dD1, with FD accuracy order 8.
 */
static REAL fd_function_dD1_fdorder8(const REAL FDPROTO_i1m1, const REAL FDPROTO_i1m2, const REAL FDPROTO_i1m3, const REAL FDPROTO_i1m4,
                                     const REAL FDPROTO_i1p1, const REAL FDPROTO_i1p2, const REAL FDPROTO_i1p3, const REAL FDPROTO_i1p4,
                                     const REAL invdxx1) {
  const REAL FDPart1_Rational_4_5 = 4.0 / 5.0;
  const REAL FDPart1_Rational_4_105 = 4.0 / 105.0;
  const REAL FDPart1_Rational_1_5 = 1.0 / 5.0;
  const REAL FDPart1_Rational_1_280 = 1.0 / 280.0;
  const REAL FD_result = invdxx1 * (FDPart1_Rational_1_280 * (FDPROTO_i1m4 - FDPROTO_i1p4) + FDPart1_Rational_1_5 * (FDPROTO_i1m2 - FDPROTO_i1p2) +
                                    FDPart1_Rational_4_105 * (-FDPROTO_i1m3 + FDPROTO_i1p3) + FDPart1_Rational_4_5 * (-FDPROTO_i1m1 + FDPROTO_i1p1));

  return FD_result;
}
/*
 * Finite difference function for operator dDD00, with FD accuracy order 8.
 */
static REAL fd_function_dDD00_fdorder8(const REAL FDPROTO, const REAL FDPROTO_i0m1, const REAL FDPROTO_i0m2, const REAL FDPROTO_i0m3,
                                       const REAL FDPROTO_i0m4, const REAL FDPROTO_i0p1, const REAL FDPROTO_i0p2, const REAL FDPROTO_i0p3,
                                       const REAL FDPROTO_i0p4, const REAL invdxx0) {
  const REAL FDPart1_Rational_205_72 = 205.0 / 72.0;
  const REAL FDPart1_Rational_1_5 = 1.0 / 5.0;
  const REAL FDPart1_Rational_1_560 = 1.0 / 560.0;
  const REAL FDPart1_Rational_8_5 = 8.0 / 5.0;
  const REAL FDPart1_Rational_8_315 = 8.0 / 315.0;
  const REAL FD_result =
      ((invdxx0) * (invdxx0)) * (-FDPROTO * FDPart1_Rational_205_72 + FDPart1_Rational_1_5 * (-FDPROTO_i0m2 - FDPROTO_i0p2) +
                                 FDPart1_Rational_1_560 * (-FDPROTO_i0m4 - FDPROTO_i0p4) + FDPart1_Rational_8_315 * (FDPROTO_i0m3 + FDPROTO_i0p3) +
                                 FDPart1_Rational_8_5 * (FDPROTO_i0m1 + FDPROTO_i0p1));

  return FD_result;
}
/*
 * Finite difference function for operator dDD01, with FD accuracy order 8.
 */
static REAL
fd_function_dDD01_fdorder8(const REAL FDPROTO_i0m1_i1m1, const REAL FDPROTO_i0m1_i1m2, const REAL FDPROTO_i0m1_i1m3, const REAL FDPROTO_i0m1_i1m4,
                           const REAL FDPROTO_i0m1_i1p1, const REAL FDPROTO_i0m1_i1p2, const REAL FDPROTO_i0m1_i1p3, const REAL FDPROTO_i0m1_i1p4,
                           const REAL FDPROTO_i0m2_i1m1, const REAL FDPROTO_i0m2_i1m2, const REAL FDPROTO_i0m2_i1m3, const REAL FDPROTO_i0m2_i1m4,
                           const REAL FDPROTO_i0m2_i1p1, const REAL FDPROTO_i0m2_i1p2, const REAL FDPROTO_i0m2_i1p3, const REAL FDPROTO_i0m2_i1p4,
                           const REAL FDPROTO_i0m3_i1m1, const REAL FDPROTO_i0m3_i1m2, const REAL FDPROTO_i0m3_i1m3, const REAL FDPROTO_i0m3_i1m4,
                           const REAL FDPROTO_i0m3_i1p1, const REAL FDPROTO_i0m3_i1p2, const REAL FDPROTO_i0m3_i1p3, const REAL FDPROTO_i0m3_i1p4,
                           const REAL FDPROTO_i0m4_i1m1, const REAL FDPROTO_i0m4_i1m2, const REAL FDPROTO_i0m4_i1m3, const REAL FDPROTO_i0m4_i1m4,
                           const REAL FDPROTO_i0m4_i1p1, const REAL FDPROTO_i0m4_i1p2, const REAL FDPROTO_i0m4_i1p3, const REAL FDPROTO_i0m4_i1p4,
                           const REAL FDPROTO_i0p1_i1m1, const REAL FDPROTO_i0p1_i1m2, const REAL FDPROTO_i0p1_i1m3, const REAL FDPROTO_i0p1_i1m4,
                           const REAL FDPROTO_i0p1_i1p1, const REAL FDPROTO_i0p1_i1p2, const REAL FDPROTO_i0p1_i1p3, const REAL FDPROTO_i0p1_i1p4,
                           const REAL FDPROTO_i0p2_i1m1, const REAL FDPROTO_i0p2_i1m2, const REAL FDPROTO_i0p2_i1m3, const REAL FDPROTO_i0p2_i1m4,
                           const REAL FDPROTO_i0p2_i1p1, const REAL FDPROTO_i0p2_i1p2, const REAL FDPROTO_i0p2_i1p3, const REAL FDPROTO_i0p2_i1p4,
                           const REAL FDPROTO_i0p3_i1m1, const REAL FDPROTO_i0p3_i1m2, const REAL FDPROTO_i0p3_i1m3, const REAL FDPROTO_i0p3_i1m4,
                           const REAL FDPROTO_i0p3_i1p1, const REAL FDPROTO_i0p3_i1p2, const REAL FDPROTO_i0p3_i1p3, const REAL FDPROTO_i0p3_i1p4,
                           const REAL FDPROTO_i0p4_i1m1, const REAL FDPROTO_i0p4_i1m2, const REAL FDPROTO_i0p4_i1m3, const REAL FDPROTO_i0p4_i1m4,
                           const REAL FDPROTO_i0p4_i1p1, const REAL FDPROTO_i0p4_i1p2, const REAL FDPROTO_i0p4_i1p3, const REAL FDPROTO_i0p4_i1p4,
                           const REAL invdxx0, const REAL invdxx1) {
  const REAL FDPart1_Rational_16_25 = 16.0 / 25.0;
  const REAL FDPart1_Rational_16_525 = 16.0 / 525.0;
  const REAL FDPart1_Rational_16_11025 = 16.0 / 11025.0;
  const REAL FDPart1_Rational_4_25 = 4.0 / 25.0;
  const REAL FDPart1_Rational_4_525 = 4.0 / 525.0;
  const REAL FDPart1_Rational_1_25 = 1.0 / 25.0;
  const REAL FDPart1_Rational_1_350 = 1.0 / 350.0;
  const REAL FDPart1_Rational_1_1400 = 1.0 / 1400.0;
  const REAL FDPart1_Rational_1_7350 = 1.0 / 7350.0;
  const REAL FDPart1_Rational_1_78400 = 1.0 / 78400.0;
  const REAL FD_result = invdxx0 * invdxx1 *
                         (FDPart1_Rational_16_11025 * (FDPROTO_i0m3_i1m3 - FDPROTO_i0m3_i1p3 - FDPROTO_i0p3_i1m3 + FDPROTO_i0p3_i1p3) +
                          FDPart1_Rational_16_25 * (FDPROTO_i0m1_i1m1 - FDPROTO_i0m1_i1p1 - FDPROTO_i0p1_i1m1 + FDPROTO_i0p1_i1p1) +
                          FDPart1_Rational_16_525 * (FDPROTO_i0m1_i1m3 - FDPROTO_i0m1_i1p3 + FDPROTO_i0m3_i1m1 - FDPROTO_i0m3_i1p1 -
                                                     FDPROTO_i0p1_i1m3 + FDPROTO_i0p1_i1p3 - FDPROTO_i0p3_i1m1 + FDPROTO_i0p3_i1p1) +
                          FDPart1_Rational_1_1400 * (FDPROTO_i0m2_i1m4 - FDPROTO_i0m2_i1p4 + FDPROTO_i0m4_i1m2 - FDPROTO_i0m4_i1p2 -
                                                     FDPROTO_i0p2_i1m4 + FDPROTO_i0p2_i1p4 - FDPROTO_i0p4_i1m2 + FDPROTO_i0p4_i1p2) +
                          FDPart1_Rational_1_25 * (FDPROTO_i0m2_i1m2 - FDPROTO_i0m2_i1p2 - FDPROTO_i0p2_i1m2 + FDPROTO_i0p2_i1p2) +
                          FDPart1_Rational_1_350 * (-FDPROTO_i0m1_i1m4 + FDPROTO_i0m1_i1p4 - FDPROTO_i0m4_i1m1 + FDPROTO_i0m4_i1p1 +
                                                    FDPROTO_i0p1_i1m4 - FDPROTO_i0p1_i1p4 + FDPROTO_i0p4_i1m1 - FDPROTO_i0p4_i1p1) +
                          FDPart1_Rational_1_7350 * (-FDPROTO_i0m3_i1m4 + FDPROTO_i0m3_i1p4 - FDPROTO_i0m4_i1m3 + FDPROTO_i0m4_i1p3 +
                                                     FDPROTO_i0p3_i1m4 - FDPROTO_i0p3_i1p4 + FDPROTO_i0p4_i1m3 - FDPROTO_i0p4_i1p3) +
                          FDPart1_Rational_1_78400 * (FDPROTO_i0m4_i1m4 - FDPROTO_i0m4_i1p4 - FDPROTO_i0p4_i1m4 + FDPROTO_i0p4_i1p4) +
                          FDPart1_Rational_4_25 * (-FDPROTO_i0m1_i1m2 + FDPROTO_i0m1_i1p2 - FDPROTO_i0m2_i1m1 + FDPROTO_i0m2_i1p1 +
                                                   FDPROTO_i0p1_i1m2 - FDPROTO_i0p1_i1p2 + FDPROTO_i0p2_i1m1 - FDPROTO_i0p2_i1p1) +
                          FDPart1_Rational_4_525 * (-FDPROTO_i0m2_i1m3 + FDPROTO_i0m2_i1p3 - FDPROTO_i0m3_i1m2 + FDPROTO_i0m3_i1p2 +
                                                    FDPROTO_i0p2_i1m3 - FDPROTO_i0p2_i1p3 + FDPROTO_i0p3_i1m2 - FDPROTO_i0p3_i1p2));

  return FD_result;
}
/*
 * Finite difference function for operator dDD11, with FD accuracy order 8.
 */
static REAL fd_function_dDD11_fdorder8(const REAL FDPROTO, const REAL FDPROTO_i1m1, const REAL FDPROTO_i1m2, const REAL FDPROTO_i1m3,
                                       const REAL FDPROTO_i1m4, const REAL FDPROTO_i1p1, const REAL FDPROTO_i1p2, const REAL FDPROTO_i1p3,
                                       const REAL FDPROTO_i1p4, const REAL invdxx1) {
  const REAL FDPart1_Rational_205_72 = 205.0 / 72.0;
  const REAL FDPart1_Rational_1_5 = 1.0 / 5.0;
  const REAL FDPart1_Rational_1_560 = 1.0 / 560.0;
  const REAL FDPart1_Rational_8_5 = 8.0 / 5.0;
  const REAL FDPart1_Rational_8_315 = 8.0 / 315.0;
  const REAL FD_result =
      ((invdxx1) * (invdxx1)) * (-FDPROTO * FDPart1_Rational_205_72 + FDPart1_Rational_1_5 * (-FDPROTO_i1m2 - FDPROTO_i1p2) +
                                 FDPart1_Rational_1_560 * (-FDPROTO_i1m4 - FDPROTO_i1p4) + FDPart1_Rational_8_315 * (FDPROTO_i1m3 + FDPROTO_i1p3) +
                                 FDPart1_Rational_8_5 * (FDPROTO_i1m1 + FDPROTO_i1p1));

  return FD_result;
}

/*
 * Compute psi4 at all interior gridpoints, part 0
 */
void psi4_part0__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
                                    const REAL *restrict in_gfs, REAL *restrict diagnostic_output_gfs) {
#include "../set_CodeParameters.h"
  if (!(Cart_originx == 0 && Cart_originy == 0 && Cart_originz == 0)) {
    fprintf(stderr, "Error: psi4 output assumes that the grid is centered on the origin.\n");
    fprintf(stderr, "       Good news: check out the C code for proposed modifications.\n");
    exit(1);
  }
#pragma omp parallel for collapse(2)
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
      for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {

        REAL mre4U0, mre4U1, mre4U2, mre4U3, mim4U0, mim4U1, mim4U2, mim4U3, n4U0, n4U1, n4U2, n4U3;
        REAL xx0, xx1, xx2;
        const int idx3 = IDX3(i0, i1, i2);
        psi4_tetrad(commondata, params, in_gfs[IDX4pt(CFGF, idx3)], in_gfs[IDX4pt(HDD00GF, idx3)], in_gfs[IDX4pt(HDD01GF, idx3)],
                    in_gfs[IDX4pt(HDD02GF, idx3)], in_gfs[IDX4pt(HDD11GF, idx3)], in_gfs[IDX4pt(HDD12GF, idx3)], in_gfs[IDX4pt(HDD22GF, idx3)],
                    &mre4U0, &mre4U1, &mre4U2, &mre4U3, &mim4U0, &mim4U1, &mim4U2, &mim4U3, &n4U0, &n4U1, &n4U2, &n4U3, xx, i0, i1, i2);
        {
          xx0 = xx[0][i0];
          xx1 = xx[1][i1];

          /* PROPOSED MODIFICATIONS FOR COMPUTING PSI4 ON GRIDS NOT CENTERED ON THE ORIGIN
              REAL xCart_rel_to_globalgrid_center[3];
              xx_to_Cart(commondata, params, xx, i0, i1, i2,  xCart_rel_to_globalgrid_center);
              int ignore_Cart_to_i0i1i2[3];  REAL xx_rel_to_globalgridorigin[3];
              Cart_to_xx_and_nearest_i0i1i2_global_grid_center(commondata, params,
             xCart_rel_to_globalgrid_center,xx_rel_to_globalgridorigin,ignore_Cart_to_i0i1i2); xx0=xx_rel_to_globalgridorigin[0];
              xx1=xx_rel_to_globalgridorigin[1];
              xx2=xx_rel_to_globalgridorigin[2];
          */
        }
        /*
         * NRPy+-Generated GF Access/FD Code, Step 1 of 2:
         * Read gridfunction(s) from main memory and compute FD stencils as needed.
         */
        const REAL aDD00 = in_gfs[IDX4(ADD00GF, i0, i1, i2)];
        const REAL aDD01 = in_gfs[IDX4(ADD01GF, i0, i1, i2)];
        const REAL aDD02 = in_gfs[IDX4(ADD02GF, i0, i1, i2)];
        const REAL aDD11 = in_gfs[IDX4(ADD11GF, i0, i1, i2)];
        const REAL aDD12 = in_gfs[IDX4(ADD12GF, i0, i1, i2)];
        const REAL aDD22 = in_gfs[IDX4(ADD22GF, i0, i1, i2)];
        const REAL cf_i0m4_i1m4 = in_gfs[IDX4(CFGF, i0 - 4, i1 - 4, i2)];
        const REAL cf_i0m3_i1m4 = in_gfs[IDX4(CFGF, i0 - 3, i1 - 4, i2)];
        const REAL cf_i0m2_i1m4 = in_gfs[IDX4(CFGF, i0 - 2, i1 - 4, i2)];
        const REAL cf_i0m1_i1m4 = in_gfs[IDX4(CFGF, i0 - 1, i1 - 4, i2)];
        const REAL cf_i1m4 = in_gfs[IDX4(CFGF, i0, i1 - 4, i2)];
        const REAL cf_i0p1_i1m4 = in_gfs[IDX4(CFGF, i0 + 1, i1 - 4, i2)];
        const REAL cf_i0p2_i1m4 = in_gfs[IDX4(CFGF, i0 + 2, i1 - 4, i2)];
        const REAL cf_i0p3_i1m4 = in_gfs[IDX4(CFGF, i0 + 3, i1 - 4, i2)];
        const REAL cf_i0p4_i1m4 = in_gfs[IDX4(CFGF, i0 + 4, i1 - 4, i2)];
        const REAL cf_i0m4_i1m3 = in_gfs[IDX4(CFGF, i0 - 4, i1 - 3, i2)];
        const REAL cf_i0m3_i1m3 = in_gfs[IDX4(CFGF, i0 - 3, i1 - 3, i2)];
        const REAL cf_i0m2_i1m3 = in_gfs[IDX4(CFGF, i0 - 2, i1 - 3, i2)];
        const REAL cf_i0m1_i1m3 = in_gfs[IDX4(CFGF, i0 - 1, i1 - 3, i2)];
        const REAL cf_i1m3 = in_gfs[IDX4(CFGF, i0, i1 - 3, i2)];
        const REAL cf_i0p1_i1m3 = in_gfs[IDX4(CFGF, i0 + 1, i1 - 3, i2)];
        const REAL cf_i0p2_i1m3 = in_gfs[IDX4(CFGF, i0 + 2, i1 - 3, i2)];
        const REAL cf_i0p3_i1m3 = in_gfs[IDX4(CFGF, i0 + 3, i1 - 3, i2)];
        const REAL cf_i0p4_i1m3 = in_gfs[IDX4(CFGF, i0 + 4, i1 - 3, i2)];
        const REAL cf_i0m4_i1m2 = in_gfs[IDX4(CFGF, i0 - 4, i1 - 2, i2)];
        const REAL cf_i0m3_i1m2 = in_gfs[IDX4(CFGF, i0 - 3, i1 - 2, i2)];
        const REAL cf_i0m2_i1m2 = in_gfs[IDX4(CFGF, i0 - 2, i1 - 2, i2)];
        const REAL cf_i0m1_i1m2 = in_gfs[IDX4(CFGF, i0 - 1, i1 - 2, i2)];
        const REAL cf_i1m2 = in_gfs[IDX4(CFGF, i0, i1 - 2, i2)];
        const REAL cf_i0p1_i1m2 = in_gfs[IDX4(CFGF, i0 + 1, i1 - 2, i2)];
        const REAL cf_i0p2_i1m2 = in_gfs[IDX4(CFGF, i0 + 2, i1 - 2, i2)];
        const REAL cf_i0p3_i1m2 = in_gfs[IDX4(CFGF, i0 + 3, i1 - 2, i2)];
        const REAL cf_i0p4_i1m2 = in_gfs[IDX4(CFGF, i0 + 4, i1 - 2, i2)];
        const REAL cf_i0m4_i1m1 = in_gfs[IDX4(CFGF, i0 - 4, i1 - 1, i2)];
        const REAL cf_i0m3_i1m1 = in_gfs[IDX4(CFGF, i0 - 3, i1 - 1, i2)];
        const REAL cf_i0m2_i1m1 = in_gfs[IDX4(CFGF, i0 - 2, i1 - 1, i2)];
        const REAL cf_i0m1_i1m1 = in_gfs[IDX4(CFGF, i0 - 1, i1 - 1, i2)];
        const REAL cf_i1m1 = in_gfs[IDX4(CFGF, i0, i1 - 1, i2)];
        const REAL cf_i0p1_i1m1 = in_gfs[IDX4(CFGF, i0 + 1, i1 - 1, i2)];
        const REAL cf_i0p2_i1m1 = in_gfs[IDX4(CFGF, i0 + 2, i1 - 1, i2)];
        const REAL cf_i0p3_i1m1 = in_gfs[IDX4(CFGF, i0 + 3, i1 - 1, i2)];
        const REAL cf_i0p4_i1m1 = in_gfs[IDX4(CFGF, i0 + 4, i1 - 1, i2)];
        const REAL cf_i0m4 = in_gfs[IDX4(CFGF, i0 - 4, i1, i2)];
        const REAL cf_i0m3 = in_gfs[IDX4(CFGF, i0 - 3, i1, i2)];
        const REAL cf_i0m2 = in_gfs[IDX4(CFGF, i0 - 2, i1, i2)];
        const REAL cf_i0m1 = in_gfs[IDX4(CFGF, i0 - 1, i1, i2)];
        const REAL cf = in_gfs[IDX4(CFGF, i0, i1, i2)];
        const REAL cf_i0p1 = in_gfs[IDX4(CFGF, i0 + 1, i1, i2)];
        const REAL cf_i0p2 = in_gfs[IDX4(CFGF, i0 + 2, i1, i2)];
        const REAL cf_i0p3 = in_gfs[IDX4(CFGF, i0 + 3, i1, i2)];
        const REAL cf_i0p4 = in_gfs[IDX4(CFGF, i0 + 4, i1, i2)];
        const REAL cf_i0m4_i1p1 = in_gfs[IDX4(CFGF, i0 - 4, i1 + 1, i2)];
        const REAL cf_i0m3_i1p1 = in_gfs[IDX4(CFGF, i0 - 3, i1 + 1, i2)];
        const REAL cf_i0m2_i1p1 = in_gfs[IDX4(CFGF, i0 - 2, i1 + 1, i2)];
        const REAL cf_i0m1_i1p1 = in_gfs[IDX4(CFGF, i0 - 1, i1 + 1, i2)];
        const REAL cf_i1p1 = in_gfs[IDX4(CFGF, i0, i1 + 1, i2)];
        const REAL cf_i0p1_i1p1 = in_gfs[IDX4(CFGF, i0 + 1, i1 + 1, i2)];
        const REAL cf_i0p2_i1p1 = in_gfs[IDX4(CFGF, i0 + 2, i1 + 1, i2)];
        const REAL cf_i0p3_i1p1 = in_gfs[IDX4(CFGF, i0 + 3, i1 + 1, i2)];
        const REAL cf_i0p4_i1p1 = in_gfs[IDX4(CFGF, i0 + 4, i1 + 1, i2)];
        const REAL cf_i0m4_i1p2 = in_gfs[IDX4(CFGF, i0 - 4, i1 + 2, i2)];
        const REAL cf_i0m3_i1p2 = in_gfs[IDX4(CFGF, i0 - 3, i1 + 2, i2)];
        const REAL cf_i0m2_i1p2 = in_gfs[IDX4(CFGF, i0 - 2, i1 + 2, i2)];
        const REAL cf_i0m1_i1p2 = in_gfs[IDX4(CFGF, i0 - 1, i1 + 2, i2)];
        const REAL cf_i1p2 = in_gfs[IDX4(CFGF, i0, i1 + 2, i2)];
        const REAL cf_i0p1_i1p2 = in_gfs[IDX4(CFGF, i0 + 1, i1 + 2, i2)];
        const REAL cf_i0p2_i1p2 = in_gfs[IDX4(CFGF, i0 + 2, i1 + 2, i2)];
        const REAL cf_i0p3_i1p2 = in_gfs[IDX4(CFGF, i0 + 3, i1 + 2, i2)];
        const REAL cf_i0p4_i1p2 = in_gfs[IDX4(CFGF, i0 + 4, i1 + 2, i2)];
        const REAL cf_i0m4_i1p3 = in_gfs[IDX4(CFGF, i0 - 4, i1 + 3, i2)];
        const REAL cf_i0m3_i1p3 = in_gfs[IDX4(CFGF, i0 - 3, i1 + 3, i2)];
        const REAL cf_i0m2_i1p3 = in_gfs[IDX4(CFGF, i0 - 2, i1 + 3, i2)];
        const REAL cf_i0m1_i1p3 = in_gfs[IDX4(CFGF, i0 - 1, i1 + 3, i2)];
        const REAL cf_i1p3 = in_gfs[IDX4(CFGF, i0, i1 + 3, i2)];
        const REAL cf_i0p1_i1p3 = in_gfs[IDX4(CFGF, i0 + 1, i1 + 3, i2)];
        const REAL cf_i0p2_i1p3 = in_gfs[IDX4(CFGF, i0 + 2, i1 + 3, i2)];
        const REAL cf_i0p3_i1p3 = in_gfs[IDX4(CFGF, i0 + 3, i1 + 3, i2)];
        const REAL cf_i0p4_i1p3 = in_gfs[IDX4(CFGF, i0 + 4, i1 + 3, i2)];
        const REAL cf_i0m4_i1p4 = in_gfs[IDX4(CFGF, i0 - 4, i1 + 4, i2)];
        const REAL cf_i0m3_i1p4 = in_gfs[IDX4(CFGF, i0 - 3, i1 + 4, i2)];
        const REAL cf_i0m2_i1p4 = in_gfs[IDX4(CFGF, i0 - 2, i1 + 4, i2)];
        const REAL cf_i0m1_i1p4 = in_gfs[IDX4(CFGF, i0 - 1, i1 + 4, i2)];
        const REAL cf_i1p4 = in_gfs[IDX4(CFGF, i0, i1 + 4, i2)];
        const REAL cf_i0p1_i1p4 = in_gfs[IDX4(CFGF, i0 + 1, i1 + 4, i2)];
        const REAL cf_i0p2_i1p4 = in_gfs[IDX4(CFGF, i0 + 2, i1 + 4, i2)];
        const REAL cf_i0p3_i1p4 = in_gfs[IDX4(CFGF, i0 + 3, i1 + 4, i2)];
        const REAL cf_i0p4_i1p4 = in_gfs[IDX4(CFGF, i0 + 4, i1 + 4, i2)];
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
        const REAL hDD01_i0m4_i1m4 = in_gfs[IDX4(HDD01GF, i0 - 4, i1 - 4, i2)];
        const REAL hDD01_i0m3_i1m4 = in_gfs[IDX4(HDD01GF, i0 - 3, i1 - 4, i2)];
        const REAL hDD01_i0m2_i1m4 = in_gfs[IDX4(HDD01GF, i0 - 2, i1 - 4, i2)];
        const REAL hDD01_i0m1_i1m4 = in_gfs[IDX4(HDD01GF, i0 - 1, i1 - 4, i2)];
        const REAL hDD01_i1m4 = in_gfs[IDX4(HDD01GF, i0, i1 - 4, i2)];
        const REAL hDD01_i0p1_i1m4 = in_gfs[IDX4(HDD01GF, i0 + 1, i1 - 4, i2)];
        const REAL hDD01_i0p2_i1m4 = in_gfs[IDX4(HDD01GF, i0 + 2, i1 - 4, i2)];
        const REAL hDD01_i0p3_i1m4 = in_gfs[IDX4(HDD01GF, i0 + 3, i1 - 4, i2)];
        const REAL hDD01_i0p4_i1m4 = in_gfs[IDX4(HDD01GF, i0 + 4, i1 - 4, i2)];
        const REAL hDD01_i0m4_i1m3 = in_gfs[IDX4(HDD01GF, i0 - 4, i1 - 3, i2)];
        const REAL hDD01_i0m3_i1m3 = in_gfs[IDX4(HDD01GF, i0 - 3, i1 - 3, i2)];
        const REAL hDD01_i0m2_i1m3 = in_gfs[IDX4(HDD01GF, i0 - 2, i1 - 3, i2)];
        const REAL hDD01_i0m1_i1m3 = in_gfs[IDX4(HDD01GF, i0 - 1, i1 - 3, i2)];
        const REAL hDD01_i1m3 = in_gfs[IDX4(HDD01GF, i0, i1 - 3, i2)];
        const REAL hDD01_i0p1_i1m3 = in_gfs[IDX4(HDD01GF, i0 + 1, i1 - 3, i2)];
        const REAL hDD01_i0p2_i1m3 = in_gfs[IDX4(HDD01GF, i0 + 2, i1 - 3, i2)];
        const REAL hDD01_i0p3_i1m3 = in_gfs[IDX4(HDD01GF, i0 + 3, i1 - 3, i2)];
        const REAL hDD01_i0p4_i1m3 = in_gfs[IDX4(HDD01GF, i0 + 4, i1 - 3, i2)];
        const REAL hDD01_i0m4_i1m2 = in_gfs[IDX4(HDD01GF, i0 - 4, i1 - 2, i2)];
        const REAL hDD01_i0m3_i1m2 = in_gfs[IDX4(HDD01GF, i0 - 3, i1 - 2, i2)];
        const REAL hDD01_i0m2_i1m2 = in_gfs[IDX4(HDD01GF, i0 - 2, i1 - 2, i2)];
        const REAL hDD01_i0m1_i1m2 = in_gfs[IDX4(HDD01GF, i0 - 1, i1 - 2, i2)];
        const REAL hDD01_i1m2 = in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2)];
        const REAL hDD01_i0p1_i1m2 = in_gfs[IDX4(HDD01GF, i0 + 1, i1 - 2, i2)];
        const REAL hDD01_i0p2_i1m2 = in_gfs[IDX4(HDD01GF, i0 + 2, i1 - 2, i2)];
        const REAL hDD01_i0p3_i1m2 = in_gfs[IDX4(HDD01GF, i0 + 3, i1 - 2, i2)];
        const REAL hDD01_i0p4_i1m2 = in_gfs[IDX4(HDD01GF, i0 + 4, i1 - 2, i2)];
        const REAL hDD01_i0m4_i1m1 = in_gfs[IDX4(HDD01GF, i0 - 4, i1 - 1, i2)];
        const REAL hDD01_i0m3_i1m1 = in_gfs[IDX4(HDD01GF, i0 - 3, i1 - 1, i2)];
        const REAL hDD01_i0m2_i1m1 = in_gfs[IDX4(HDD01GF, i0 - 2, i1 - 1, i2)];
        const REAL hDD01_i0m1_i1m1 = in_gfs[IDX4(HDD01GF, i0 - 1, i1 - 1, i2)];
        const REAL hDD01_i1m1 = in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2)];
        const REAL hDD01_i0p1_i1m1 = in_gfs[IDX4(HDD01GF, i0 + 1, i1 - 1, i2)];
        const REAL hDD01_i0p2_i1m1 = in_gfs[IDX4(HDD01GF, i0 + 2, i1 - 1, i2)];
        const REAL hDD01_i0p3_i1m1 = in_gfs[IDX4(HDD01GF, i0 + 3, i1 - 1, i2)];
        const REAL hDD01_i0p4_i1m1 = in_gfs[IDX4(HDD01GF, i0 + 4, i1 - 1, i2)];
        const REAL hDD01_i0m4 = in_gfs[IDX4(HDD01GF, i0 - 4, i1, i2)];
        const REAL hDD01_i0m3 = in_gfs[IDX4(HDD01GF, i0 - 3, i1, i2)];
        const REAL hDD01_i0m2 = in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2)];
        const REAL hDD01_i0m1 = in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2)];
        const REAL hDD01 = in_gfs[IDX4(HDD01GF, i0, i1, i2)];
        const REAL hDD01_i0p1 = in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2)];
        const REAL hDD01_i0p2 = in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2)];
        const REAL hDD01_i0p3 = in_gfs[IDX4(HDD01GF, i0 + 3, i1, i2)];
        const REAL hDD01_i0p4 = in_gfs[IDX4(HDD01GF, i0 + 4, i1, i2)];
        const REAL hDD01_i0m4_i1p1 = in_gfs[IDX4(HDD01GF, i0 - 4, i1 + 1, i2)];
        const REAL hDD01_i0m3_i1p1 = in_gfs[IDX4(HDD01GF, i0 - 3, i1 + 1, i2)];
        const REAL hDD01_i0m2_i1p1 = in_gfs[IDX4(HDD01GF, i0 - 2, i1 + 1, i2)];
        const REAL hDD01_i0m1_i1p1 = in_gfs[IDX4(HDD01GF, i0 - 1, i1 + 1, i2)];
        const REAL hDD01_i1p1 = in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2)];
        const REAL hDD01_i0p1_i1p1 = in_gfs[IDX4(HDD01GF, i0 + 1, i1 + 1, i2)];
        const REAL hDD01_i0p2_i1p1 = in_gfs[IDX4(HDD01GF, i0 + 2, i1 + 1, i2)];
        const REAL hDD01_i0p3_i1p1 = in_gfs[IDX4(HDD01GF, i0 + 3, i1 + 1, i2)];
        const REAL hDD01_i0p4_i1p1 = in_gfs[IDX4(HDD01GF, i0 + 4, i1 + 1, i2)];
        const REAL hDD01_i0m4_i1p2 = in_gfs[IDX4(HDD01GF, i0 - 4, i1 + 2, i2)];
        const REAL hDD01_i0m3_i1p2 = in_gfs[IDX4(HDD01GF, i0 - 3, i1 + 2, i2)];
        const REAL hDD01_i0m2_i1p2 = in_gfs[IDX4(HDD01GF, i0 - 2, i1 + 2, i2)];
        const REAL hDD01_i0m1_i1p2 = in_gfs[IDX4(HDD01GF, i0 - 1, i1 + 2, i2)];
        const REAL hDD01_i1p2 = in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2)];
        const REAL hDD01_i0p1_i1p2 = in_gfs[IDX4(HDD01GF, i0 + 1, i1 + 2, i2)];
        const REAL hDD01_i0p2_i1p2 = in_gfs[IDX4(HDD01GF, i0 + 2, i1 + 2, i2)];
        const REAL hDD01_i0p3_i1p2 = in_gfs[IDX4(HDD01GF, i0 + 3, i1 + 2, i2)];
        const REAL hDD01_i0p4_i1p2 = in_gfs[IDX4(HDD01GF, i0 + 4, i1 + 2, i2)];
        const REAL hDD01_i0m4_i1p3 = in_gfs[IDX4(HDD01GF, i0 - 4, i1 + 3, i2)];
        const REAL hDD01_i0m3_i1p3 = in_gfs[IDX4(HDD01GF, i0 - 3, i1 + 3, i2)];
        const REAL hDD01_i0m2_i1p3 = in_gfs[IDX4(HDD01GF, i0 - 2, i1 + 3, i2)];
        const REAL hDD01_i0m1_i1p3 = in_gfs[IDX4(HDD01GF, i0 - 1, i1 + 3, i2)];
        const REAL hDD01_i1p3 = in_gfs[IDX4(HDD01GF, i0, i1 + 3, i2)];
        const REAL hDD01_i0p1_i1p3 = in_gfs[IDX4(HDD01GF, i0 + 1, i1 + 3, i2)];
        const REAL hDD01_i0p2_i1p3 = in_gfs[IDX4(HDD01GF, i0 + 2, i1 + 3, i2)];
        const REAL hDD01_i0p3_i1p3 = in_gfs[IDX4(HDD01GF, i0 + 3, i1 + 3, i2)];
        const REAL hDD01_i0p4_i1p3 = in_gfs[IDX4(HDD01GF, i0 + 4, i1 + 3, i2)];
        const REAL hDD01_i0m4_i1p4 = in_gfs[IDX4(HDD01GF, i0 - 4, i1 + 4, i2)];
        const REAL hDD01_i0m3_i1p4 = in_gfs[IDX4(HDD01GF, i0 - 3, i1 + 4, i2)];
        const REAL hDD01_i0m2_i1p4 = in_gfs[IDX4(HDD01GF, i0 - 2, i1 + 4, i2)];
        const REAL hDD01_i0m1_i1p4 = in_gfs[IDX4(HDD01GF, i0 - 1, i1 + 4, i2)];
        const REAL hDD01_i1p4 = in_gfs[IDX4(HDD01GF, i0, i1 + 4, i2)];
        const REAL hDD01_i0p1_i1p4 = in_gfs[IDX4(HDD01GF, i0 + 1, i1 + 4, i2)];
        const REAL hDD01_i0p2_i1p4 = in_gfs[IDX4(HDD01GF, i0 + 2, i1 + 4, i2)];
        const REAL hDD01_i0p3_i1p4 = in_gfs[IDX4(HDD01GF, i0 + 3, i1 + 4, i2)];
        const REAL hDD01_i0p4_i1p4 = in_gfs[IDX4(HDD01GF, i0 + 4, i1 + 4, i2)];
        const REAL hDD02_i0m4_i1m4 = in_gfs[IDX4(HDD02GF, i0 - 4, i1 - 4, i2)];
        const REAL hDD02_i0m3_i1m4 = in_gfs[IDX4(HDD02GF, i0 - 3, i1 - 4, i2)];
        const REAL hDD02_i0m2_i1m4 = in_gfs[IDX4(HDD02GF, i0 - 2, i1 - 4, i2)];
        const REAL hDD02_i0m1_i1m4 = in_gfs[IDX4(HDD02GF, i0 - 1, i1 - 4, i2)];
        const REAL hDD02_i1m4 = in_gfs[IDX4(HDD02GF, i0, i1 - 4, i2)];
        const REAL hDD02_i0p1_i1m4 = in_gfs[IDX4(HDD02GF, i0 + 1, i1 - 4, i2)];
        const REAL hDD02_i0p2_i1m4 = in_gfs[IDX4(HDD02GF, i0 + 2, i1 - 4, i2)];
        const REAL hDD02_i0p3_i1m4 = in_gfs[IDX4(HDD02GF, i0 + 3, i1 - 4, i2)];
        const REAL hDD02_i0p4_i1m4 = in_gfs[IDX4(HDD02GF, i0 + 4, i1 - 4, i2)];
        const REAL hDD02_i0m4_i1m3 = in_gfs[IDX4(HDD02GF, i0 - 4, i1 - 3, i2)];
        const REAL hDD02_i0m3_i1m3 = in_gfs[IDX4(HDD02GF, i0 - 3, i1 - 3, i2)];
        const REAL hDD02_i0m2_i1m3 = in_gfs[IDX4(HDD02GF, i0 - 2, i1 - 3, i2)];
        const REAL hDD02_i0m1_i1m3 = in_gfs[IDX4(HDD02GF, i0 - 1, i1 - 3, i2)];
        const REAL hDD02_i1m3 = in_gfs[IDX4(HDD02GF, i0, i1 - 3, i2)];
        const REAL hDD02_i0p1_i1m3 = in_gfs[IDX4(HDD02GF, i0 + 1, i1 - 3, i2)];
        const REAL hDD02_i0p2_i1m3 = in_gfs[IDX4(HDD02GF, i0 + 2, i1 - 3, i2)];
        const REAL hDD02_i0p3_i1m3 = in_gfs[IDX4(HDD02GF, i0 + 3, i1 - 3, i2)];
        const REAL hDD02_i0p4_i1m3 = in_gfs[IDX4(HDD02GF, i0 + 4, i1 - 3, i2)];
        const REAL hDD02_i0m4_i1m2 = in_gfs[IDX4(HDD02GF, i0 - 4, i1 - 2, i2)];
        const REAL hDD02_i0m3_i1m2 = in_gfs[IDX4(HDD02GF, i0 - 3, i1 - 2, i2)];
        const REAL hDD02_i0m2_i1m2 = in_gfs[IDX4(HDD02GF, i0 - 2, i1 - 2, i2)];
        const REAL hDD02_i0m1_i1m2 = in_gfs[IDX4(HDD02GF, i0 - 1, i1 - 2, i2)];
        const REAL hDD02_i1m2 = in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2)];
        const REAL hDD02_i0p1_i1m2 = in_gfs[IDX4(HDD02GF, i0 + 1, i1 - 2, i2)];
        const REAL hDD02_i0p2_i1m2 = in_gfs[IDX4(HDD02GF, i0 + 2, i1 - 2, i2)];
        const REAL hDD02_i0p3_i1m2 = in_gfs[IDX4(HDD02GF, i0 + 3, i1 - 2, i2)];
        const REAL hDD02_i0p4_i1m2 = in_gfs[IDX4(HDD02GF, i0 + 4, i1 - 2, i2)];
        const REAL hDD02_i0m4_i1m1 = in_gfs[IDX4(HDD02GF, i0 - 4, i1 - 1, i2)];
        const REAL hDD02_i0m3_i1m1 = in_gfs[IDX4(HDD02GF, i0 - 3, i1 - 1, i2)];
        const REAL hDD02_i0m2_i1m1 = in_gfs[IDX4(HDD02GF, i0 - 2, i1 - 1, i2)];
        const REAL hDD02_i0m1_i1m1 = in_gfs[IDX4(HDD02GF, i0 - 1, i1 - 1, i2)];
        const REAL hDD02_i1m1 = in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2)];
        const REAL hDD02_i0p1_i1m1 = in_gfs[IDX4(HDD02GF, i0 + 1, i1 - 1, i2)];
        const REAL hDD02_i0p2_i1m1 = in_gfs[IDX4(HDD02GF, i0 + 2, i1 - 1, i2)];
        const REAL hDD02_i0p3_i1m1 = in_gfs[IDX4(HDD02GF, i0 + 3, i1 - 1, i2)];
        const REAL hDD02_i0p4_i1m1 = in_gfs[IDX4(HDD02GF, i0 + 4, i1 - 1, i2)];
        const REAL hDD02_i0m4 = in_gfs[IDX4(HDD02GF, i0 - 4, i1, i2)];
        const REAL hDD02_i0m3 = in_gfs[IDX4(HDD02GF, i0 - 3, i1, i2)];
        const REAL hDD02_i0m2 = in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2)];
        const REAL hDD02_i0m1 = in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2)];
        const REAL hDD02 = in_gfs[IDX4(HDD02GF, i0, i1, i2)];
        const REAL hDD02_i0p1 = in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2)];
        const REAL hDD02_i0p2 = in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2)];
        const REAL hDD02_i0p3 = in_gfs[IDX4(HDD02GF, i0 + 3, i1, i2)];
        const REAL hDD02_i0p4 = in_gfs[IDX4(HDD02GF, i0 + 4, i1, i2)];
        const REAL hDD02_i0m4_i1p1 = in_gfs[IDX4(HDD02GF, i0 - 4, i1 + 1, i2)];
        const REAL hDD02_i0m3_i1p1 = in_gfs[IDX4(HDD02GF, i0 - 3, i1 + 1, i2)];
        const REAL hDD02_i0m2_i1p1 = in_gfs[IDX4(HDD02GF, i0 - 2, i1 + 1, i2)];
        const REAL hDD02_i0m1_i1p1 = in_gfs[IDX4(HDD02GF, i0 - 1, i1 + 1, i2)];
        const REAL hDD02_i1p1 = in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2)];
        const REAL hDD02_i0p1_i1p1 = in_gfs[IDX4(HDD02GF, i0 + 1, i1 + 1, i2)];
        const REAL hDD02_i0p2_i1p1 = in_gfs[IDX4(HDD02GF, i0 + 2, i1 + 1, i2)];
        const REAL hDD02_i0p3_i1p1 = in_gfs[IDX4(HDD02GF, i0 + 3, i1 + 1, i2)];
        const REAL hDD02_i0p4_i1p1 = in_gfs[IDX4(HDD02GF, i0 + 4, i1 + 1, i2)];
        const REAL hDD02_i0m4_i1p2 = in_gfs[IDX4(HDD02GF, i0 - 4, i1 + 2, i2)];
        const REAL hDD02_i0m3_i1p2 = in_gfs[IDX4(HDD02GF, i0 - 3, i1 + 2, i2)];
        const REAL hDD02_i0m2_i1p2 = in_gfs[IDX4(HDD02GF, i0 - 2, i1 + 2, i2)];
        const REAL hDD02_i0m1_i1p2 = in_gfs[IDX4(HDD02GF, i0 - 1, i1 + 2, i2)];
        const REAL hDD02_i1p2 = in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2)];
        const REAL hDD02_i0p1_i1p2 = in_gfs[IDX4(HDD02GF, i0 + 1, i1 + 2, i2)];
        const REAL hDD02_i0p2_i1p2 = in_gfs[IDX4(HDD02GF, i0 + 2, i1 + 2, i2)];
        const REAL hDD02_i0p3_i1p2 = in_gfs[IDX4(HDD02GF, i0 + 3, i1 + 2, i2)];
        const REAL hDD02_i0p4_i1p2 = in_gfs[IDX4(HDD02GF, i0 + 4, i1 + 2, i2)];
        const REAL hDD02_i0m4_i1p3 = in_gfs[IDX4(HDD02GF, i0 - 4, i1 + 3, i2)];
        const REAL hDD02_i0m3_i1p3 = in_gfs[IDX4(HDD02GF, i0 - 3, i1 + 3, i2)];
        const REAL hDD02_i0m2_i1p3 = in_gfs[IDX4(HDD02GF, i0 - 2, i1 + 3, i2)];
        const REAL hDD02_i0m1_i1p3 = in_gfs[IDX4(HDD02GF, i0 - 1, i1 + 3, i2)];
        const REAL hDD02_i1p3 = in_gfs[IDX4(HDD02GF, i0, i1 + 3, i2)];
        const REAL hDD02_i0p1_i1p3 = in_gfs[IDX4(HDD02GF, i0 + 1, i1 + 3, i2)];
        const REAL hDD02_i0p2_i1p3 = in_gfs[IDX4(HDD02GF, i0 + 2, i1 + 3, i2)];
        const REAL hDD02_i0p3_i1p3 = in_gfs[IDX4(HDD02GF, i0 + 3, i1 + 3, i2)];
        const REAL hDD02_i0p4_i1p3 = in_gfs[IDX4(HDD02GF, i0 + 4, i1 + 3, i2)];
        const REAL hDD02_i0m4_i1p4 = in_gfs[IDX4(HDD02GF, i0 - 4, i1 + 4, i2)];
        const REAL hDD02_i0m3_i1p4 = in_gfs[IDX4(HDD02GF, i0 - 3, i1 + 4, i2)];
        const REAL hDD02_i0m2_i1p4 = in_gfs[IDX4(HDD02GF, i0 - 2, i1 + 4, i2)];
        const REAL hDD02_i0m1_i1p4 = in_gfs[IDX4(HDD02GF, i0 - 1, i1 + 4, i2)];
        const REAL hDD02_i1p4 = in_gfs[IDX4(HDD02GF, i0, i1 + 4, i2)];
        const REAL hDD02_i0p1_i1p4 = in_gfs[IDX4(HDD02GF, i0 + 1, i1 + 4, i2)];
        const REAL hDD02_i0p2_i1p4 = in_gfs[IDX4(HDD02GF, i0 + 2, i1 + 4, i2)];
        const REAL hDD02_i0p3_i1p4 = in_gfs[IDX4(HDD02GF, i0 + 3, i1 + 4, i2)];
        const REAL hDD02_i0p4_i1p4 = in_gfs[IDX4(HDD02GF, i0 + 4, i1 + 4, i2)];
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
        const REAL hDD12_i0m4_i1m4 = in_gfs[IDX4(HDD12GF, i0 - 4, i1 - 4, i2)];
        const REAL hDD12_i0m3_i1m4 = in_gfs[IDX4(HDD12GF, i0 - 3, i1 - 4, i2)];
        const REAL hDD12_i0m2_i1m4 = in_gfs[IDX4(HDD12GF, i0 - 2, i1 - 4, i2)];
        const REAL hDD12_i0m1_i1m4 = in_gfs[IDX4(HDD12GF, i0 - 1, i1 - 4, i2)];
        const REAL hDD12_i1m4 = in_gfs[IDX4(HDD12GF, i0, i1 - 4, i2)];
        const REAL hDD12_i0p1_i1m4 = in_gfs[IDX4(HDD12GF, i0 + 1, i1 - 4, i2)];
        const REAL hDD12_i0p2_i1m4 = in_gfs[IDX4(HDD12GF, i0 + 2, i1 - 4, i2)];
        const REAL hDD12_i0p3_i1m4 = in_gfs[IDX4(HDD12GF, i0 + 3, i1 - 4, i2)];
        const REAL hDD12_i0p4_i1m4 = in_gfs[IDX4(HDD12GF, i0 + 4, i1 - 4, i2)];
        const REAL hDD12_i0m4_i1m3 = in_gfs[IDX4(HDD12GF, i0 - 4, i1 - 3, i2)];
        const REAL hDD12_i0m3_i1m3 = in_gfs[IDX4(HDD12GF, i0 - 3, i1 - 3, i2)];
        const REAL hDD12_i0m2_i1m3 = in_gfs[IDX4(HDD12GF, i0 - 2, i1 - 3, i2)];
        const REAL hDD12_i0m1_i1m3 = in_gfs[IDX4(HDD12GF, i0 - 1, i1 - 3, i2)];
        const REAL hDD12_i1m3 = in_gfs[IDX4(HDD12GF, i0, i1 - 3, i2)];
        const REAL hDD12_i0p1_i1m3 = in_gfs[IDX4(HDD12GF, i0 + 1, i1 - 3, i2)];
        const REAL hDD12_i0p2_i1m3 = in_gfs[IDX4(HDD12GF, i0 + 2, i1 - 3, i2)];
        const REAL hDD12_i0p3_i1m3 = in_gfs[IDX4(HDD12GF, i0 + 3, i1 - 3, i2)];
        const REAL hDD12_i0p4_i1m3 = in_gfs[IDX4(HDD12GF, i0 + 4, i1 - 3, i2)];
        const REAL hDD12_i0m4_i1m2 = in_gfs[IDX4(HDD12GF, i0 - 4, i1 - 2, i2)];
        const REAL hDD12_i0m3_i1m2 = in_gfs[IDX4(HDD12GF, i0 - 3, i1 - 2, i2)];
        const REAL hDD12_i0m2_i1m2 = in_gfs[IDX4(HDD12GF, i0 - 2, i1 - 2, i2)];
        const REAL hDD12_i0m1_i1m2 = in_gfs[IDX4(HDD12GF, i0 - 1, i1 - 2, i2)];
        const REAL hDD12_i1m2 = in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2)];
        const REAL hDD12_i0p1_i1m2 = in_gfs[IDX4(HDD12GF, i0 + 1, i1 - 2, i2)];
        const REAL hDD12_i0p2_i1m2 = in_gfs[IDX4(HDD12GF, i0 + 2, i1 - 2, i2)];
        const REAL hDD12_i0p3_i1m2 = in_gfs[IDX4(HDD12GF, i0 + 3, i1 - 2, i2)];
        const REAL hDD12_i0p4_i1m2 = in_gfs[IDX4(HDD12GF, i0 + 4, i1 - 2, i2)];
        const REAL hDD12_i0m4_i1m1 = in_gfs[IDX4(HDD12GF, i0 - 4, i1 - 1, i2)];
        const REAL hDD12_i0m3_i1m1 = in_gfs[IDX4(HDD12GF, i0 - 3, i1 - 1, i2)];
        const REAL hDD12_i0m2_i1m1 = in_gfs[IDX4(HDD12GF, i0 - 2, i1 - 1, i2)];
        const REAL hDD12_i0m1_i1m1 = in_gfs[IDX4(HDD12GF, i0 - 1, i1 - 1, i2)];
        const REAL hDD12_i1m1 = in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2)];
        const REAL hDD12_i0p1_i1m1 = in_gfs[IDX4(HDD12GF, i0 + 1, i1 - 1, i2)];
        const REAL hDD12_i0p2_i1m1 = in_gfs[IDX4(HDD12GF, i0 + 2, i1 - 1, i2)];
        const REAL hDD12_i0p3_i1m1 = in_gfs[IDX4(HDD12GF, i0 + 3, i1 - 1, i2)];
        const REAL hDD12_i0p4_i1m1 = in_gfs[IDX4(HDD12GF, i0 + 4, i1 - 1, i2)];
        const REAL hDD12_i0m4 = in_gfs[IDX4(HDD12GF, i0 - 4, i1, i2)];
        const REAL hDD12_i0m3 = in_gfs[IDX4(HDD12GF, i0 - 3, i1, i2)];
        const REAL hDD12_i0m2 = in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2)];
        const REAL hDD12_i0m1 = in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2)];
        const REAL hDD12 = in_gfs[IDX4(HDD12GF, i0, i1, i2)];
        const REAL hDD12_i0p1 = in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2)];
        const REAL hDD12_i0p2 = in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2)];
        const REAL hDD12_i0p3 = in_gfs[IDX4(HDD12GF, i0 + 3, i1, i2)];
        const REAL hDD12_i0p4 = in_gfs[IDX4(HDD12GF, i0 + 4, i1, i2)];
        const REAL hDD12_i0m4_i1p1 = in_gfs[IDX4(HDD12GF, i0 - 4, i1 + 1, i2)];
        const REAL hDD12_i0m3_i1p1 = in_gfs[IDX4(HDD12GF, i0 - 3, i1 + 1, i2)];
        const REAL hDD12_i0m2_i1p1 = in_gfs[IDX4(HDD12GF, i0 - 2, i1 + 1, i2)];
        const REAL hDD12_i0m1_i1p1 = in_gfs[IDX4(HDD12GF, i0 - 1, i1 + 1, i2)];
        const REAL hDD12_i1p1 = in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2)];
        const REAL hDD12_i0p1_i1p1 = in_gfs[IDX4(HDD12GF, i0 + 1, i1 + 1, i2)];
        const REAL hDD12_i0p2_i1p1 = in_gfs[IDX4(HDD12GF, i0 + 2, i1 + 1, i2)];
        const REAL hDD12_i0p3_i1p1 = in_gfs[IDX4(HDD12GF, i0 + 3, i1 + 1, i2)];
        const REAL hDD12_i0p4_i1p1 = in_gfs[IDX4(HDD12GF, i0 + 4, i1 + 1, i2)];
        const REAL hDD12_i0m4_i1p2 = in_gfs[IDX4(HDD12GF, i0 - 4, i1 + 2, i2)];
        const REAL hDD12_i0m3_i1p2 = in_gfs[IDX4(HDD12GF, i0 - 3, i1 + 2, i2)];
        const REAL hDD12_i0m2_i1p2 = in_gfs[IDX4(HDD12GF, i0 - 2, i1 + 2, i2)];
        const REAL hDD12_i0m1_i1p2 = in_gfs[IDX4(HDD12GF, i0 - 1, i1 + 2, i2)];
        const REAL hDD12_i1p2 = in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2)];
        const REAL hDD12_i0p1_i1p2 = in_gfs[IDX4(HDD12GF, i0 + 1, i1 + 2, i2)];
        const REAL hDD12_i0p2_i1p2 = in_gfs[IDX4(HDD12GF, i0 + 2, i1 + 2, i2)];
        const REAL hDD12_i0p3_i1p2 = in_gfs[IDX4(HDD12GF, i0 + 3, i1 + 2, i2)];
        const REAL hDD12_i0p4_i1p2 = in_gfs[IDX4(HDD12GF, i0 + 4, i1 + 2, i2)];
        const REAL hDD12_i0m4_i1p3 = in_gfs[IDX4(HDD12GF, i0 - 4, i1 + 3, i2)];
        const REAL hDD12_i0m3_i1p3 = in_gfs[IDX4(HDD12GF, i0 - 3, i1 + 3, i2)];
        const REAL hDD12_i0m2_i1p3 = in_gfs[IDX4(HDD12GF, i0 - 2, i1 + 3, i2)];
        const REAL hDD12_i0m1_i1p3 = in_gfs[IDX4(HDD12GF, i0 - 1, i1 + 3, i2)];
        const REAL hDD12_i1p3 = in_gfs[IDX4(HDD12GF, i0, i1 + 3, i2)];
        const REAL hDD12_i0p1_i1p3 = in_gfs[IDX4(HDD12GF, i0 + 1, i1 + 3, i2)];
        const REAL hDD12_i0p2_i1p3 = in_gfs[IDX4(HDD12GF, i0 + 2, i1 + 3, i2)];
        const REAL hDD12_i0p3_i1p3 = in_gfs[IDX4(HDD12GF, i0 + 3, i1 + 3, i2)];
        const REAL hDD12_i0p4_i1p3 = in_gfs[IDX4(HDD12GF, i0 + 4, i1 + 3, i2)];
        const REAL hDD12_i0m4_i1p4 = in_gfs[IDX4(HDD12GF, i0 - 4, i1 + 4, i2)];
        const REAL hDD12_i0m3_i1p4 = in_gfs[IDX4(HDD12GF, i0 - 3, i1 + 4, i2)];
        const REAL hDD12_i0m2_i1p4 = in_gfs[IDX4(HDD12GF, i0 - 2, i1 + 4, i2)];
        const REAL hDD12_i0m1_i1p4 = in_gfs[IDX4(HDD12GF, i0 - 1, i1 + 4, i2)];
        const REAL hDD12_i1p4 = in_gfs[IDX4(HDD12GF, i0, i1 + 4, i2)];
        const REAL hDD12_i0p1_i1p4 = in_gfs[IDX4(HDD12GF, i0 + 1, i1 + 4, i2)];
        const REAL hDD12_i0p2_i1p4 = in_gfs[IDX4(HDD12GF, i0 + 2, i1 + 4, i2)];
        const REAL hDD12_i0p3_i1p4 = in_gfs[IDX4(HDD12GF, i0 + 3, i1 + 4, i2)];
        const REAL hDD12_i0p4_i1p4 = in_gfs[IDX4(HDD12GF, i0 + 4, i1 + 4, i2)];
        const REAL hDD22_i0m4_i1m4 = in_gfs[IDX4(HDD22GF, i0 - 4, i1 - 4, i2)];
        const REAL hDD22_i0m3_i1m4 = in_gfs[IDX4(HDD22GF, i0 - 3, i1 - 4, i2)];
        const REAL hDD22_i0m2_i1m4 = in_gfs[IDX4(HDD22GF, i0 - 2, i1 - 4, i2)];
        const REAL hDD22_i0m1_i1m4 = in_gfs[IDX4(HDD22GF, i0 - 1, i1 - 4, i2)];
        const REAL hDD22_i1m4 = in_gfs[IDX4(HDD22GF, i0, i1 - 4, i2)];
        const REAL hDD22_i0p1_i1m4 = in_gfs[IDX4(HDD22GF, i0 + 1, i1 - 4, i2)];
        const REAL hDD22_i0p2_i1m4 = in_gfs[IDX4(HDD22GF, i0 + 2, i1 - 4, i2)];
        const REAL hDD22_i0p3_i1m4 = in_gfs[IDX4(HDD22GF, i0 + 3, i1 - 4, i2)];
        const REAL hDD22_i0p4_i1m4 = in_gfs[IDX4(HDD22GF, i0 + 4, i1 - 4, i2)];
        const REAL hDD22_i0m4_i1m3 = in_gfs[IDX4(HDD22GF, i0 - 4, i1 - 3, i2)];
        const REAL hDD22_i0m3_i1m3 = in_gfs[IDX4(HDD22GF, i0 - 3, i1 - 3, i2)];
        const REAL hDD22_i0m2_i1m3 = in_gfs[IDX4(HDD22GF, i0 - 2, i1 - 3, i2)];
        const REAL hDD22_i0m1_i1m3 = in_gfs[IDX4(HDD22GF, i0 - 1, i1 - 3, i2)];
        const REAL hDD22_i1m3 = in_gfs[IDX4(HDD22GF, i0, i1 - 3, i2)];
        const REAL hDD22_i0p1_i1m3 = in_gfs[IDX4(HDD22GF, i0 + 1, i1 - 3, i2)];
        const REAL hDD22_i0p2_i1m3 = in_gfs[IDX4(HDD22GF, i0 + 2, i1 - 3, i2)];
        const REAL hDD22_i0p3_i1m3 = in_gfs[IDX4(HDD22GF, i0 + 3, i1 - 3, i2)];
        const REAL hDD22_i0p4_i1m3 = in_gfs[IDX4(HDD22GF, i0 + 4, i1 - 3, i2)];
        const REAL hDD22_i0m4_i1m2 = in_gfs[IDX4(HDD22GF, i0 - 4, i1 - 2, i2)];
        const REAL hDD22_i0m3_i1m2 = in_gfs[IDX4(HDD22GF, i0 - 3, i1 - 2, i2)];
        const REAL hDD22_i0m2_i1m2 = in_gfs[IDX4(HDD22GF, i0 - 2, i1 - 2, i2)];
        const REAL hDD22_i0m1_i1m2 = in_gfs[IDX4(HDD22GF, i0 - 1, i1 - 2, i2)];
        const REAL hDD22_i1m2 = in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2)];
        const REAL hDD22_i0p1_i1m2 = in_gfs[IDX4(HDD22GF, i0 + 1, i1 - 2, i2)];
        const REAL hDD22_i0p2_i1m2 = in_gfs[IDX4(HDD22GF, i0 + 2, i1 - 2, i2)];
        const REAL hDD22_i0p3_i1m2 = in_gfs[IDX4(HDD22GF, i0 + 3, i1 - 2, i2)];
        const REAL hDD22_i0p4_i1m2 = in_gfs[IDX4(HDD22GF, i0 + 4, i1 - 2, i2)];
        const REAL hDD22_i0m4_i1m1 = in_gfs[IDX4(HDD22GF, i0 - 4, i1 - 1, i2)];
        const REAL hDD22_i0m3_i1m1 = in_gfs[IDX4(HDD22GF, i0 - 3, i1 - 1, i2)];
        const REAL hDD22_i0m2_i1m1 = in_gfs[IDX4(HDD22GF, i0 - 2, i1 - 1, i2)];
        const REAL hDD22_i0m1_i1m1 = in_gfs[IDX4(HDD22GF, i0 - 1, i1 - 1, i2)];
        const REAL hDD22_i1m1 = in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2)];
        const REAL hDD22_i0p1_i1m1 = in_gfs[IDX4(HDD22GF, i0 + 1, i1 - 1, i2)];
        const REAL hDD22_i0p2_i1m1 = in_gfs[IDX4(HDD22GF, i0 + 2, i1 - 1, i2)];
        const REAL hDD22_i0p3_i1m1 = in_gfs[IDX4(HDD22GF, i0 + 3, i1 - 1, i2)];
        const REAL hDD22_i0p4_i1m1 = in_gfs[IDX4(HDD22GF, i0 + 4, i1 - 1, i2)];
        const REAL hDD22_i0m4 = in_gfs[IDX4(HDD22GF, i0 - 4, i1, i2)];
        const REAL hDD22_i0m3 = in_gfs[IDX4(HDD22GF, i0 - 3, i1, i2)];
        const REAL hDD22_i0m2 = in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2)];
        const REAL hDD22_i0m1 = in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2)];
        const REAL hDD22 = in_gfs[IDX4(HDD22GF, i0, i1, i2)];
        const REAL hDD22_i0p1 = in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2)];
        const REAL hDD22_i0p2 = in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2)];
        const REAL hDD22_i0p3 = in_gfs[IDX4(HDD22GF, i0 + 3, i1, i2)];
        const REAL hDD22_i0p4 = in_gfs[IDX4(HDD22GF, i0 + 4, i1, i2)];
        const REAL hDD22_i0m4_i1p1 = in_gfs[IDX4(HDD22GF, i0 - 4, i1 + 1, i2)];
        const REAL hDD22_i0m3_i1p1 = in_gfs[IDX4(HDD22GF, i0 - 3, i1 + 1, i2)];
        const REAL hDD22_i0m2_i1p1 = in_gfs[IDX4(HDD22GF, i0 - 2, i1 + 1, i2)];
        const REAL hDD22_i0m1_i1p1 = in_gfs[IDX4(HDD22GF, i0 - 1, i1 + 1, i2)];
        const REAL hDD22_i1p1 = in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2)];
        const REAL hDD22_i0p1_i1p1 = in_gfs[IDX4(HDD22GF, i0 + 1, i1 + 1, i2)];
        const REAL hDD22_i0p2_i1p1 = in_gfs[IDX4(HDD22GF, i0 + 2, i1 + 1, i2)];
        const REAL hDD22_i0p3_i1p1 = in_gfs[IDX4(HDD22GF, i0 + 3, i1 + 1, i2)];
        const REAL hDD22_i0p4_i1p1 = in_gfs[IDX4(HDD22GF, i0 + 4, i1 + 1, i2)];
        const REAL hDD22_i0m4_i1p2 = in_gfs[IDX4(HDD22GF, i0 - 4, i1 + 2, i2)];
        const REAL hDD22_i0m3_i1p2 = in_gfs[IDX4(HDD22GF, i0 - 3, i1 + 2, i2)];
        const REAL hDD22_i0m2_i1p2 = in_gfs[IDX4(HDD22GF, i0 - 2, i1 + 2, i2)];
        const REAL hDD22_i0m1_i1p2 = in_gfs[IDX4(HDD22GF, i0 - 1, i1 + 2, i2)];
        const REAL hDD22_i1p2 = in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2)];
        const REAL hDD22_i0p1_i1p2 = in_gfs[IDX4(HDD22GF, i0 + 1, i1 + 2, i2)];
        const REAL hDD22_i0p2_i1p2 = in_gfs[IDX4(HDD22GF, i0 + 2, i1 + 2, i2)];
        const REAL hDD22_i0p3_i1p2 = in_gfs[IDX4(HDD22GF, i0 + 3, i1 + 2, i2)];
        const REAL hDD22_i0p4_i1p2 = in_gfs[IDX4(HDD22GF, i0 + 4, i1 + 2, i2)];
        const REAL hDD22_i0m4_i1p3 = in_gfs[IDX4(HDD22GF, i0 - 4, i1 + 3, i2)];
        const REAL hDD22_i0m3_i1p3 = in_gfs[IDX4(HDD22GF, i0 - 3, i1 + 3, i2)];
        const REAL hDD22_i0m2_i1p3 = in_gfs[IDX4(HDD22GF, i0 - 2, i1 + 3, i2)];
        const REAL hDD22_i0m1_i1p3 = in_gfs[IDX4(HDD22GF, i0 - 1, i1 + 3, i2)];
        const REAL hDD22_i1p3 = in_gfs[IDX4(HDD22GF, i0, i1 + 3, i2)];
        const REAL hDD22_i0p1_i1p3 = in_gfs[IDX4(HDD22GF, i0 + 1, i1 + 3, i2)];
        const REAL hDD22_i0p2_i1p3 = in_gfs[IDX4(HDD22GF, i0 + 2, i1 + 3, i2)];
        const REAL hDD22_i0p3_i1p3 = in_gfs[IDX4(HDD22GF, i0 + 3, i1 + 3, i2)];
        const REAL hDD22_i0p4_i1p3 = in_gfs[IDX4(HDD22GF, i0 + 4, i1 + 3, i2)];
        const REAL hDD22_i0m4_i1p4 = in_gfs[IDX4(HDD22GF, i0 - 4, i1 + 4, i2)];
        const REAL hDD22_i0m3_i1p4 = in_gfs[IDX4(HDD22GF, i0 - 3, i1 + 4, i2)];
        const REAL hDD22_i0m2_i1p4 = in_gfs[IDX4(HDD22GF, i0 - 2, i1 + 4, i2)];
        const REAL hDD22_i0m1_i1p4 = in_gfs[IDX4(HDD22GF, i0 - 1, i1 + 4, i2)];
        const REAL hDD22_i1p4 = in_gfs[IDX4(HDD22GF, i0, i1 + 4, i2)];
        const REAL hDD22_i0p1_i1p4 = in_gfs[IDX4(HDD22GF, i0 + 1, i1 + 4, i2)];
        const REAL hDD22_i0p2_i1p4 = in_gfs[IDX4(HDD22GF, i0 + 2, i1 + 4, i2)];
        const REAL hDD22_i0p3_i1p4 = in_gfs[IDX4(HDD22GF, i0 + 3, i1 + 4, i2)];
        const REAL hDD22_i0p4_i1p4 = in_gfs[IDX4(HDD22GF, i0 + 4, i1 + 4, i2)];
        const REAL trK = in_gfs[IDX4(TRKGF, i0, i1, i2)];
        const REAL cf_dD0 = fd_function_dD0_fdorder8(cf_i0m1, cf_i0m2, cf_i0m3, cf_i0m4, cf_i0p1, cf_i0p2, cf_i0p3, cf_i0p4, invdxx0);
        const REAL cf_dD1 = fd_function_dD1_fdorder8(cf_i1m1, cf_i1m2, cf_i1m3, cf_i1m4, cf_i1p1, cf_i1p2, cf_i1p3, cf_i1p4, invdxx1);
        const REAL cf_dDD00 = fd_function_dDD00_fdorder8(cf, cf_i0m1, cf_i0m2, cf_i0m3, cf_i0m4, cf_i0p1, cf_i0p2, cf_i0p3, cf_i0p4, invdxx0);
        const REAL cf_dDD01 = fd_function_dDD01_fdorder8(
            cf_i0m1_i1m1, cf_i0m1_i1m2, cf_i0m1_i1m3, cf_i0m1_i1m4, cf_i0m1_i1p1, cf_i0m1_i1p2, cf_i0m1_i1p3, cf_i0m1_i1p4, cf_i0m2_i1m1,
            cf_i0m2_i1m2, cf_i0m2_i1m3, cf_i0m2_i1m4, cf_i0m2_i1p1, cf_i0m2_i1p2, cf_i0m2_i1p3, cf_i0m2_i1p4, cf_i0m3_i1m1, cf_i0m3_i1m2,
            cf_i0m3_i1m3, cf_i0m3_i1m4, cf_i0m3_i1p1, cf_i0m3_i1p2, cf_i0m3_i1p3, cf_i0m3_i1p4, cf_i0m4_i1m1, cf_i0m4_i1m2, cf_i0m4_i1m3,
            cf_i0m4_i1m4, cf_i0m4_i1p1, cf_i0m4_i1p2, cf_i0m4_i1p3, cf_i0m4_i1p4, cf_i0p1_i1m1, cf_i0p1_i1m2, cf_i0p1_i1m3, cf_i0p1_i1m4,
            cf_i0p1_i1p1, cf_i0p1_i1p2, cf_i0p1_i1p3, cf_i0p1_i1p4, cf_i0p2_i1m1, cf_i0p2_i1m2, cf_i0p2_i1m3, cf_i0p2_i1m4, cf_i0p2_i1p1,
            cf_i0p2_i1p2, cf_i0p2_i1p3, cf_i0p2_i1p4, cf_i0p3_i1m1, cf_i0p3_i1m2, cf_i0p3_i1m3, cf_i0p3_i1m4, cf_i0p3_i1p1, cf_i0p3_i1p2,
            cf_i0p3_i1p3, cf_i0p3_i1p4, cf_i0p4_i1m1, cf_i0p4_i1m2, cf_i0p4_i1m3, cf_i0p4_i1m4, cf_i0p4_i1p1, cf_i0p4_i1p2, cf_i0p4_i1p3,
            cf_i0p4_i1p4, invdxx0, invdxx1);
        const REAL cf_dDD11 = fd_function_dDD11_fdorder8(cf, cf_i1m1, cf_i1m2, cf_i1m3, cf_i1m4, cf_i1p1, cf_i1p2, cf_i1p3, cf_i1p4, invdxx1);
        const REAL hDD_dD000 =
            fd_function_dD0_fdorder8(hDD00_i0m1, hDD00_i0m2, hDD00_i0m3, hDD00_i0m4, hDD00_i0p1, hDD00_i0p2, hDD00_i0p3, hDD00_i0p4, invdxx0);
        const REAL hDD_dD001 =
            fd_function_dD1_fdorder8(hDD00_i1m1, hDD00_i1m2, hDD00_i1m3, hDD00_i1m4, hDD00_i1p1, hDD00_i1p2, hDD00_i1p3, hDD00_i1p4, invdxx1);
        const REAL hDD_dD010 =
            fd_function_dD0_fdorder8(hDD01_i0m1, hDD01_i0m2, hDD01_i0m3, hDD01_i0m4, hDD01_i0p1, hDD01_i0p2, hDD01_i0p3, hDD01_i0p4, invdxx0);
        const REAL hDD_dD011 =
            fd_function_dD1_fdorder8(hDD01_i1m1, hDD01_i1m2, hDD01_i1m3, hDD01_i1m4, hDD01_i1p1, hDD01_i1p2, hDD01_i1p3, hDD01_i1p4, invdxx1);
        const REAL hDD_dD020 =
            fd_function_dD0_fdorder8(hDD02_i0m1, hDD02_i0m2, hDD02_i0m3, hDD02_i0m4, hDD02_i0p1, hDD02_i0p2, hDD02_i0p3, hDD02_i0p4, invdxx0);
        const REAL hDD_dD021 =
            fd_function_dD1_fdorder8(hDD02_i1m1, hDD02_i1m2, hDD02_i1m3, hDD02_i1m4, hDD02_i1p1, hDD02_i1p2, hDD02_i1p3, hDD02_i1p4, invdxx1);
        const REAL hDD_dD110 =
            fd_function_dD0_fdorder8(hDD11_i0m1, hDD11_i0m2, hDD11_i0m3, hDD11_i0m4, hDD11_i0p1, hDD11_i0p2, hDD11_i0p3, hDD11_i0p4, invdxx0);
        const REAL hDD_dD111 =
            fd_function_dD1_fdorder8(hDD11_i1m1, hDD11_i1m2, hDD11_i1m3, hDD11_i1m4, hDD11_i1p1, hDD11_i1p2, hDD11_i1p3, hDD11_i1p4, invdxx1);
        const REAL hDD_dD120 =
            fd_function_dD0_fdorder8(hDD12_i0m1, hDD12_i0m2, hDD12_i0m3, hDD12_i0m4, hDD12_i0p1, hDD12_i0p2, hDD12_i0p3, hDD12_i0p4, invdxx0);
        const REAL hDD_dD121 =
            fd_function_dD1_fdorder8(hDD12_i1m1, hDD12_i1m2, hDD12_i1m3, hDD12_i1m4, hDD12_i1p1, hDD12_i1p2, hDD12_i1p3, hDD12_i1p4, invdxx1);
        const REAL hDD_dD220 =
            fd_function_dD0_fdorder8(hDD22_i0m1, hDD22_i0m2, hDD22_i0m3, hDD22_i0m4, hDD22_i0p1, hDD22_i0p2, hDD22_i0p3, hDD22_i0p4, invdxx0);
        const REAL hDD_dD221 =
            fd_function_dD1_fdorder8(hDD22_i1m1, hDD22_i1m2, hDD22_i1m3, hDD22_i1m4, hDD22_i1p1, hDD22_i1p2, hDD22_i1p3, hDD22_i1p4, invdxx1);
        const REAL hDD_dDD0011 = fd_function_dDD11_fdorder8(hDD00, hDD00_i1m1, hDD00_i1m2, hDD00_i1m3, hDD00_i1m4, hDD00_i1p1, hDD00_i1p2, hDD00_i1p3,
                                                            hDD00_i1p4, invdxx1);
        const REAL hDD_dDD0101 = fd_function_dDD01_fdorder8(
            hDD01_i0m1_i1m1, hDD01_i0m1_i1m2, hDD01_i0m1_i1m3, hDD01_i0m1_i1m4, hDD01_i0m1_i1p1, hDD01_i0m1_i1p2, hDD01_i0m1_i1p3, hDD01_i0m1_i1p4,
            hDD01_i0m2_i1m1, hDD01_i0m2_i1m2, hDD01_i0m2_i1m3, hDD01_i0m2_i1m4, hDD01_i0m2_i1p1, hDD01_i0m2_i1p2, hDD01_i0m2_i1p3, hDD01_i0m2_i1p4,
            hDD01_i0m3_i1m1, hDD01_i0m3_i1m2, hDD01_i0m3_i1m3, hDD01_i0m3_i1m4, hDD01_i0m3_i1p1, hDD01_i0m3_i1p2, hDD01_i0m3_i1p3, hDD01_i0m3_i1p4,
            hDD01_i0m4_i1m1, hDD01_i0m4_i1m2, hDD01_i0m4_i1m3, hDD01_i0m4_i1m4, hDD01_i0m4_i1p1, hDD01_i0m4_i1p2, hDD01_i0m4_i1p3, hDD01_i0m4_i1p4,
            hDD01_i0p1_i1m1, hDD01_i0p1_i1m2, hDD01_i0p1_i1m3, hDD01_i0p1_i1m4, hDD01_i0p1_i1p1, hDD01_i0p1_i1p2, hDD01_i0p1_i1p3, hDD01_i0p1_i1p4,
            hDD01_i0p2_i1m1, hDD01_i0p2_i1m2, hDD01_i0p2_i1m3, hDD01_i0p2_i1m4, hDD01_i0p2_i1p1, hDD01_i0p2_i1p2, hDD01_i0p2_i1p3, hDD01_i0p2_i1p4,
            hDD01_i0p3_i1m1, hDD01_i0p3_i1m2, hDD01_i0p3_i1m3, hDD01_i0p3_i1m4, hDD01_i0p3_i1p1, hDD01_i0p3_i1p2, hDD01_i0p3_i1p3, hDD01_i0p3_i1p4,
            hDD01_i0p4_i1m1, hDD01_i0p4_i1m2, hDD01_i0p4_i1m3, hDD01_i0p4_i1m4, hDD01_i0p4_i1p1, hDD01_i0p4_i1p2, hDD01_i0p4_i1p3, hDD01_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL hDD_dDD0201 = fd_function_dDD01_fdorder8(
            hDD02_i0m1_i1m1, hDD02_i0m1_i1m2, hDD02_i0m1_i1m3, hDD02_i0m1_i1m4, hDD02_i0m1_i1p1, hDD02_i0m1_i1p2, hDD02_i0m1_i1p3, hDD02_i0m1_i1p4,
            hDD02_i0m2_i1m1, hDD02_i0m2_i1m2, hDD02_i0m2_i1m3, hDD02_i0m2_i1m4, hDD02_i0m2_i1p1, hDD02_i0m2_i1p2, hDD02_i0m2_i1p3, hDD02_i0m2_i1p4,
            hDD02_i0m3_i1m1, hDD02_i0m3_i1m2, hDD02_i0m3_i1m3, hDD02_i0m3_i1m4, hDD02_i0m3_i1p1, hDD02_i0m3_i1p2, hDD02_i0m3_i1p3, hDD02_i0m3_i1p4,
            hDD02_i0m4_i1m1, hDD02_i0m4_i1m2, hDD02_i0m4_i1m3, hDD02_i0m4_i1m4, hDD02_i0m4_i1p1, hDD02_i0m4_i1p2, hDD02_i0m4_i1p3, hDD02_i0m4_i1p4,
            hDD02_i0p1_i1m1, hDD02_i0p1_i1m2, hDD02_i0p1_i1m3, hDD02_i0p1_i1m4, hDD02_i0p1_i1p1, hDD02_i0p1_i1p2, hDD02_i0p1_i1p3, hDD02_i0p1_i1p4,
            hDD02_i0p2_i1m1, hDD02_i0p2_i1m2, hDD02_i0p2_i1m3, hDD02_i0p2_i1m4, hDD02_i0p2_i1p1, hDD02_i0p2_i1p2, hDD02_i0p2_i1p3, hDD02_i0p2_i1p4,
            hDD02_i0p3_i1m1, hDD02_i0p3_i1m2, hDD02_i0p3_i1m3, hDD02_i0p3_i1m4, hDD02_i0p3_i1p1, hDD02_i0p3_i1p2, hDD02_i0p3_i1p3, hDD02_i0p3_i1p4,
            hDD02_i0p4_i1m1, hDD02_i0p4_i1m2, hDD02_i0p4_i1m3, hDD02_i0p4_i1m4, hDD02_i0p4_i1p1, hDD02_i0p4_i1p2, hDD02_i0p4_i1p3, hDD02_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL hDD_dDD0211 = fd_function_dDD11_fdorder8(hDD02, hDD02_i1m1, hDD02_i1m2, hDD02_i1m3, hDD02_i1m4, hDD02_i1p1, hDD02_i1p2, hDD02_i1p3,
                                                            hDD02_i1p4, invdxx1);
        const REAL hDD_dDD1100 = fd_function_dDD00_fdorder8(hDD11, hDD11_i0m1, hDD11_i0m2, hDD11_i0m3, hDD11_i0m4, hDD11_i0p1, hDD11_i0p2, hDD11_i0p3,
                                                            hDD11_i0p4, invdxx0);
        const REAL hDD_dDD1200 = fd_function_dDD00_fdorder8(hDD12, hDD12_i0m1, hDD12_i0m2, hDD12_i0m3, hDD12_i0m4, hDD12_i0p1, hDD12_i0p2, hDD12_i0p3,
                                                            hDD12_i0p4, invdxx0);
        const REAL hDD_dDD1201 = fd_function_dDD01_fdorder8(
            hDD12_i0m1_i1m1, hDD12_i0m1_i1m2, hDD12_i0m1_i1m3, hDD12_i0m1_i1m4, hDD12_i0m1_i1p1, hDD12_i0m1_i1p2, hDD12_i0m1_i1p3, hDD12_i0m1_i1p4,
            hDD12_i0m2_i1m1, hDD12_i0m2_i1m2, hDD12_i0m2_i1m3, hDD12_i0m2_i1m4, hDD12_i0m2_i1p1, hDD12_i0m2_i1p2, hDD12_i0m2_i1p3, hDD12_i0m2_i1p4,
            hDD12_i0m3_i1m1, hDD12_i0m3_i1m2, hDD12_i0m3_i1m3, hDD12_i0m3_i1m4, hDD12_i0m3_i1p1, hDD12_i0m3_i1p2, hDD12_i0m3_i1p3, hDD12_i0m3_i1p4,
            hDD12_i0m4_i1m1, hDD12_i0m4_i1m2, hDD12_i0m4_i1m3, hDD12_i0m4_i1m4, hDD12_i0m4_i1p1, hDD12_i0m4_i1p2, hDD12_i0m4_i1p3, hDD12_i0m4_i1p4,
            hDD12_i0p1_i1m1, hDD12_i0p1_i1m2, hDD12_i0p1_i1m3, hDD12_i0p1_i1m4, hDD12_i0p1_i1p1, hDD12_i0p1_i1p2, hDD12_i0p1_i1p3, hDD12_i0p1_i1p4,
            hDD12_i0p2_i1m1, hDD12_i0p2_i1m2, hDD12_i0p2_i1m3, hDD12_i0p2_i1m4, hDD12_i0p2_i1p1, hDD12_i0p2_i1p2, hDD12_i0p2_i1p3, hDD12_i0p2_i1p4,
            hDD12_i0p3_i1m1, hDD12_i0p3_i1m2, hDD12_i0p3_i1m3, hDD12_i0p3_i1m4, hDD12_i0p3_i1p1, hDD12_i0p3_i1p2, hDD12_i0p3_i1p3, hDD12_i0p3_i1p4,
            hDD12_i0p4_i1m1, hDD12_i0p4_i1m2, hDD12_i0p4_i1m3, hDD12_i0p4_i1m4, hDD12_i0p4_i1p1, hDD12_i0p4_i1p2, hDD12_i0p4_i1p3, hDD12_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL hDD_dDD2200 = fd_function_dDD00_fdorder8(hDD22, hDD22_i0m1, hDD22_i0m2, hDD22_i0m3, hDD22_i0m4, hDD22_i0p1, hDD22_i0p2, hDD22_i0p3,
                                                            hDD22_i0p4, invdxx0);
        const REAL hDD_dDD2201 = fd_function_dDD01_fdorder8(
            hDD22_i0m1_i1m1, hDD22_i0m1_i1m2, hDD22_i0m1_i1m3, hDD22_i0m1_i1m4, hDD22_i0m1_i1p1, hDD22_i0m1_i1p2, hDD22_i0m1_i1p3, hDD22_i0m1_i1p4,
            hDD22_i0m2_i1m1, hDD22_i0m2_i1m2, hDD22_i0m2_i1m3, hDD22_i0m2_i1m4, hDD22_i0m2_i1p1, hDD22_i0m2_i1p2, hDD22_i0m2_i1p3, hDD22_i0m2_i1p4,
            hDD22_i0m3_i1m1, hDD22_i0m3_i1m2, hDD22_i0m3_i1m3, hDD22_i0m3_i1m4, hDD22_i0m3_i1p1, hDD22_i0m3_i1p2, hDD22_i0m3_i1p3, hDD22_i0m3_i1p4,
            hDD22_i0m4_i1m1, hDD22_i0m4_i1m2, hDD22_i0m4_i1m3, hDD22_i0m4_i1m4, hDD22_i0m4_i1p1, hDD22_i0m4_i1p2, hDD22_i0m4_i1p3, hDD22_i0m4_i1p4,
            hDD22_i0p1_i1m1, hDD22_i0p1_i1m2, hDD22_i0p1_i1m3, hDD22_i0p1_i1m4, hDD22_i0p1_i1p1, hDD22_i0p1_i1p2, hDD22_i0p1_i1p3, hDD22_i0p1_i1p4,
            hDD22_i0p2_i1m1, hDD22_i0p2_i1m2, hDD22_i0p2_i1m3, hDD22_i0p2_i1m4, hDD22_i0p2_i1p1, hDD22_i0p2_i1p2, hDD22_i0p2_i1p3, hDD22_i0p2_i1p4,
            hDD22_i0p3_i1m1, hDD22_i0p3_i1m2, hDD22_i0p3_i1m3, hDD22_i0p3_i1m4, hDD22_i0p3_i1p1, hDD22_i0p3_i1p2, hDD22_i0p3_i1p3, hDD22_i0p3_i1p4,
            hDD22_i0p4_i1m1, hDD22_i0p4_i1m2, hDD22_i0p4_i1m3, hDD22_i0p4_i1m4, hDD22_i0p4_i1p1, hDD22_i0p4_i1p2, hDD22_i0p4_i1p3, hDD22_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL hDD_dDD2211 = fd_function_dDD11_fdorder8(hDD22, hDD22_i1m1, hDD22_i1m2, hDD22_i1m3, hDD22_i1m4, hDD22_i1p1, hDD22_i1p2, hDD22_i1p3,
                                                            hDD22_i1p4, invdxx1);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        const REAL FDPart3tmp0 = ((n4U3) * (n4U3));
        const REAL FDPart3tmp1 = mre4U2 * mre4U3;
        const REAL FDPart3tmp2 = mim4U2 * mim4U3;
        const REAL FDPart3tmp3 = sin(xx1);
        const REAL FDPart3tmp5 = (1.0 / (SINHW));
        const REAL FDPart3tmp16 = (1.0 / ((cf) * (cf) * (cf)));
        const REAL FDPart3tmp20 = (1.0 / ((cf) * (cf)));
        const REAL FDPart3tmp21 = cos(xx1);
        const REAL FDPart3tmp33 = (1.0 / ((cf) * (cf) * (cf) * (cf)));
        const REAL FDPart3tmp40 = pow(cf, -6);
        const REAL FDPart3tmp105 = ((mre4U3) * (mre4U3));
        const REAL FDPart3tmp106 = n4U2 * n4U3;
        const REAL FDPart3tmp107 = ((mim4U3) * (mim4U3));
        const REAL FDPart3tmp108 = mre4U1 * mre4U3;
        const REAL FDPart3tmp109 = mim4U1 * mim4U3;
        const REAL FDPart3tmp123 = n4U1 * n4U3;
        const REAL FDPart3tmp139 = ((n4U2) * (n4U2));
        const REAL FDPart3tmp161 = ((mre4U2) * (mre4U2));
        const REAL FDPart3tmp162 = ((mim4U2) * (mim4U2));
        const REAL FDPart3tmp163 = mre4U1 * mre4U2;
        const REAL FDPart3tmp164 = mim4U1 * mim4U2;
        const REAL FDPart3tmp188 = n4U1 * n4U2;
        const REAL FDPart3tmp205 = ((n4U1) * (n4U1));
        const REAL FDPart3tmp206 = (1.0 / ((SINHW) * (SINHW)));
        const REAL FDPart3tmp230 = ((mre4U1) * (mre4U1));
        const REAL FDPart3tmp231 = ((mim4U1) * (mim4U1));
        const REAL FDPart3tmp246 = (1.0 / 3.0) * trK;
        const REAL FDPart3tmp251 = (1.0 / 2.0) / cf;
        const REAL FDPart3tmp460 = 2 * mim4U3 * mre4U3;
        const REAL FDPart3tmp461 = mim4U2 * mre4U3;
        const REAL FDPart3tmp462 = mim4U3 * mre4U2;
        const REAL FDPart3tmp463 = mim4U1 * mre4U3;
        const REAL FDPart3tmp464 = mim4U3 * mre4U1;
        const REAL FDPart3tmp467 = mim4U2 * mre4U2;
        const REAL FDPart3tmp469 = mim4U1 * mre4U2;
        const REAL FDPart3tmp470 = mim4U2 * mre4U1;
        const REAL FDPart3tmp475 = mim4U1 * mre4U1;
        const REAL FDPart3tmp4 = ((FDPart3tmp3) * (FDPart3tmp3));
        const REAL FDPart3tmp11 = exp(FDPart3tmp5) - exp(-FDPart3tmp5);
        const REAL FDPart3tmp17 = FDPart3tmp16 * cf_dD1;
        const REAL FDPart3tmp52 = FDPart3tmp16 * cf_dD0;
        const REAL FDPart3tmp124 = FDPart3tmp106 * FDPart3tmp108 - FDPart3tmp106 * FDPart3tmp109;
        const REAL FDPart3tmp138 = FDPart3tmp1 * FDPart3tmp123 - FDPart3tmp123 * FDPart3tmp2;
        const REAL FDPart3tmp143 = 2 * FDPart3tmp20;
        const REAL FDPart3tmp165 = FDPart3tmp106 * FDPart3tmp163 - FDPart3tmp106 * FDPart3tmp164;
        const REAL FDPart3tmp189 = FDPart3tmp1 * FDPart3tmp188 - FDPart3tmp188 * FDPart3tmp2;
        const REAL FDPart3tmp190 = FDPart3tmp123 * FDPart3tmp163 - FDPart3tmp123 * FDPart3tmp164;
        const REAL FDPart3tmp204 = FDPart3tmp108 * FDPart3tmp188 - FDPart3tmp109 * FDPart3tmp188;
        const REAL FDPart3tmp250 = 4 * FDPart3tmp33;
        const REAL FDPart3tmp252 = (1.0 / 2.0) * FDPart3tmp20;
        const REAL FDPart3tmp253 = 4 * FDPart3tmp20;
        const REAL FDPart3tmp317 = FDPart3tmp0 * FDPart3tmp163 - FDPart3tmp0 * FDPart3tmp164;
        const REAL FDPart3tmp352 = FDPart3tmp105 * FDPart3tmp188 - FDPart3tmp107 * FDPart3tmp188;
        const REAL FDPart3tmp381 = FDPart3tmp108 * FDPart3tmp139 - FDPart3tmp109 * FDPart3tmp139;
        const REAL FDPart3tmp391 = FDPart3tmp123 * FDPart3tmp161 - FDPart3tmp123 * FDPart3tmp162;
        const REAL FDPart3tmp416 = FDPart3tmp106 * FDPart3tmp230 - FDPart3tmp106 * FDPart3tmp231;
        const REAL FDPart3tmp426 = FDPart3tmp1 * FDPart3tmp205 - FDPart3tmp2 * FDPart3tmp205;
        const REAL FDPart3tmp465 = -FDPart3tmp106 * FDPart3tmp463 - FDPart3tmp106 * FDPart3tmp464;
        const REAL FDPart3tmp466 = -FDPart3tmp123 * FDPart3tmp461 - FDPart3tmp123 * FDPart3tmp462;
        const REAL FDPart3tmp468 = 2 * FDPart3tmp467;
        const REAL FDPart3tmp471 = -FDPart3tmp106 * FDPart3tmp469 - FDPart3tmp106 * FDPart3tmp470;
        const REAL FDPart3tmp472 = -FDPart3tmp188 * FDPart3tmp461 - FDPart3tmp188 * FDPart3tmp462;
        const REAL FDPart3tmp473 = -FDPart3tmp123 * FDPart3tmp469 - FDPart3tmp123 * FDPart3tmp470;
        const REAL FDPart3tmp474 = -FDPart3tmp188 * FDPart3tmp463 - FDPart3tmp188 * FDPart3tmp464;
        const REAL FDPart3tmp476 = 2 * FDPart3tmp475;
        const REAL FDPart3tmp479 = -FDPart3tmp0 * FDPart3tmp469 - FDPart3tmp0 * FDPart3tmp470;
        const REAL FDPart3tmp481 = -FDPart3tmp139 * FDPart3tmp463 - FDPart3tmp139 * FDPart3tmp464;
        const REAL FDPart3tmp483 = -FDPart3tmp205 * FDPart3tmp461 - FDPart3tmp205 * FDPart3tmp462;
        const REAL FDPart3tmp7 = exp(FDPart3tmp5 * xx0);
        const REAL FDPart3tmp8 = exp(-FDPart3tmp5 * xx0);
        const REAL FDPart3tmp12 = ((AMPL) * (AMPL)) / ((FDPart3tmp11) * (FDPart3tmp11));
        const REAL FDPart3tmp18 = 2 * FDPart3tmp17;
        const REAL FDPart3tmp34 = ((AMPL) * (AMPL) * (AMPL) * (AMPL)) / ((FDPart3tmp11) * (FDPart3tmp11) * (FDPart3tmp11) * (FDPart3tmp11));
        const REAL FDPart3tmp53 = 2 * FDPart3tmp52;
        const REAL FDPart3tmp254 = (1.0 / 2.0) * FDPart3tmp250 * ((cf_dD1) * (cf_dD1)) +
                                   (1.0 / 2.0) * FDPart3tmp253 * (-FDPart3tmp251 * cf_dDD11 + FDPart3tmp252 * ((cf_dD1) * (cf_dD1)));
        const REAL FDPart3tmp286 = (1.0 / 2.0) * FDPart3tmp250 * ((cf_dD0) * (cf_dD0)) +
                                   (1.0 / 2.0) * FDPart3tmp253 * (-FDPart3tmp251 * cf_dDD00 + FDPart3tmp252 * ((cf_dD0) * (cf_dD0)));
        const REAL FDPart3tmp338 = FDPart3tmp250 * cf_dD0 * cf_dD1 + FDPart3tmp253 * (-FDPart3tmp251 * cf_dDD01 + FDPart3tmp252 * cf_dD0 * cf_dD1);
        const REAL FDPart3tmp9 = FDPart3tmp7 - FDPart3tmp8;
        const REAL FDPart3tmp35 = FDPart3tmp33 * FDPart3tmp34;
        const REAL FDPart3tmp339 = (1.0 / 2.0) * FDPart3tmp338;
        const REAL FDPart3tmp30 = FDPart3tmp5 * FDPart3tmp7 + FDPart3tmp5 * FDPart3tmp8;
        const REAL FDPart3tmp36 = FDPart3tmp3 * FDPart3tmp35 * hDD01;
        const REAL FDPart3tmp41 = FDPart3tmp4 * ((FDPart3tmp9) * (FDPart3tmp9) * (FDPart3tmp9) * (FDPart3tmp9));
        const REAL FDPart3tmp55 = 2 * FDPart3tmp5 * FDPart3tmp7 + 2 * FDPart3tmp5 * FDPart3tmp8;
        const REAL FDPart3tmp56 = FDPart3tmp12 * FDPart3tmp9;
        const REAL FDPart3tmp213 = 2 * FDPart3tmp206 * FDPart3tmp7 - 2 * FDPart3tmp206 * FDPart3tmp8;
        const REAL FDPart3tmp13 = FDPart3tmp12 * ((FDPart3tmp9) * (FDPart3tmp9));
        const REAL FDPart3tmp31 = ((FDPart3tmp30) * (FDPart3tmp30));
        const REAL FDPart3tmp43 = FDPart3tmp41 * ((hDD12) * (hDD12));
        const REAL FDPart3tmp57 = FDPart3tmp55 * FDPart3tmp56;
        const REAL FDPart3tmp62 = FDPart3tmp30 * ((FDPart3tmp9) * (FDPart3tmp9) * (FDPart3tmp9));
        const REAL FDPart3tmp63 = FDPart3tmp30 * FDPart3tmp56;
        const REAL FDPart3tmp209 = FDPart3tmp56 * (FDPart3tmp206 * FDPart3tmp7 - FDPart3tmp206 * FDPart3tmp8);
        const REAL FDPart3tmp293 = FDPart3tmp213 * FDPart3tmp56;
        const REAL FDPart3tmp14 = FDPart3tmp13 * FDPart3tmp4;
        const REAL FDPart3tmp22 = FDPart3tmp13 * FDPart3tmp3;
        const REAL FDPart3tmp32 = FDPart3tmp31 * ((FDPart3tmp9) * (FDPart3tmp9));
        const REAL FDPart3tmp37 = FDPart3tmp12 * FDPart3tmp31;
        const REAL FDPart3tmp45 = FDPart3tmp13 * hDD11 + FDPart3tmp13;
        const REAL FDPart3tmp58 = FDPart3tmp4 * FDPart3tmp57;
        const REAL FDPart3tmp73 = FDPart3tmp3 * FDPart3tmp63;
        const REAL FDPart3tmp84 = FDPart3tmp63 * hDD01;
        const REAL FDPart3tmp90 = FDPart3tmp20 * FDPart3tmp63;
        const REAL FDPart3tmp98 = FDPart3tmp13 * FDPart3tmp20;
        const REAL FDPart3tmp148 = FDPart3tmp13 * hDD_dD110 + FDPart3tmp57 * hDD11 + FDPart3tmp57;
        const REAL FDPart3tmp215 = FDPart3tmp12 * FDPart3tmp213 * FDPart3tmp30;
        const REAL FDPart3tmp294 = FDPart3tmp12 * FDPart3tmp30 * FDPart3tmp55;
        const REAL FDPart3tmp15 = FDPart3tmp14 * hDD22 + FDPart3tmp14;
        const REAL FDPart3tmp38 = FDPart3tmp37 * hDD00 + FDPart3tmp37;
        const REAL FDPart3tmp39 = FDPart3tmp22 * hDD12;
        const REAL FDPart3tmp46 = FDPart3tmp32 * ((hDD02) * (hDD02));
        const REAL FDPart3tmp48 = FDPart3tmp32 * ((hDD01) * (hDD01));
        const REAL FDPart3tmp59 = FDPart3tmp14 * hDD_dD220 + FDPart3tmp58 * hDD22 + FDPart3tmp58;
        const REAL FDPart3tmp65 = FDPart3tmp3 * FDPart3tmp63 * hDD02;
        const REAL FDPart3tmp70 = FDPart3tmp22 * hDD_dD120 + FDPart3tmp3 * FDPart3tmp57 * hDD12;
        const REAL FDPart3tmp74 = FDPart3tmp21 * FDPart3tmp63 * hDD02 + FDPart3tmp73 * hDD_dD021;
        const REAL FDPart3tmp103 = FDPart3tmp90 * hDD01;
        const REAL FDPart3tmp140 = -FDPart3tmp18 * FDPart3tmp45 + FDPart3tmp98 * hDD_dD111;
        const REAL FDPart3tmp142 = FDPart3tmp13 * FDPart3tmp21 * hDD12 + FDPart3tmp22 * hDD_dD121;
        const REAL FDPart3tmp210 = FDPart3tmp209 * FDPart3tmp3 + FDPart3tmp3 * FDPart3tmp37;
        const REAL FDPart3tmp257 = 2 * FDPart3tmp13 * ((FDPart3tmp21) * (FDPart3tmp21)) - 2 * FDPart3tmp14;
        const REAL FDPart3tmp259 = FDPart3tmp20 * FDPart3tmp45;
        const REAL FDPart3tmp295 = FDPart3tmp293 * FDPart3tmp4 + FDPart3tmp294 * FDPart3tmp4;
        const REAL FDPart3tmp334 = 2 * FDPart3tmp3 * FDPart3tmp57;
        const REAL FDPart3tmp366 = FDPart3tmp252 * (FDPart3tmp13 * FDPart3tmp21 * hDD_dD120 + FDPart3tmp21 * FDPart3tmp57 * hDD12 +
                                                    FDPart3tmp22 * hDD_dDD1201 + FDPart3tmp3 * FDPart3tmp57 * hDD_dD121);
        const REAL FDPart3tmp427 = FDPart3tmp252 * FDPart3tmp37 * hDD_dDD0011;
        const REAL FDPart3tmp428 = FDPart3tmp18 * FDPart3tmp37 * hDD_dD001;
        const REAL FDPart3tmp429 = FDPart3tmp53 * FDPart3tmp63 * hDD_dD011;
        const REAL FDPart3tmp24 = 2 * FDPart3tmp21 * FDPart3tmp22;
        const REAL FDPart3tmp79 = FDPart3tmp15 * FDPart3tmp33;
        const REAL FDPart3tmp92 = FDPart3tmp3 * FDPart3tmp90 * hDD02;
        const REAL FDPart3tmp100 = FDPart3tmp3 * FDPart3tmp98 * hDD12;
        const REAL FDPart3tmp144 = FDPart3tmp142 * FDPart3tmp143 - 4 * FDPart3tmp17 * FDPart3tmp39;
        const REAL FDPart3tmp150 =
            FDPart3tmp143 * FDPart3tmp63 * hDD_dD011 - FDPart3tmp148 * FDPart3tmp20 - 4 * FDPart3tmp17 * FDPart3tmp84 + FDPart3tmp45 * FDPart3tmp53;
        const REAL FDPart3tmp167 = FDPart3tmp20 * FDPart3tmp37 * hDD_dD001;
        const REAL FDPart3tmp170 = FDPart3tmp148 * FDPart3tmp20 - FDPart3tmp45 * FDPart3tmp53;
        const REAL FDPart3tmp211 = FDPart3tmp210 * hDD02 + FDPart3tmp73 * hDD_dD020;
        const REAL FDPart3tmp216 = FDPart3tmp20 * (FDPart3tmp215 * hDD00 + FDPart3tmp215 + FDPart3tmp37 * hDD_dD000) - FDPart3tmp38 * FDPart3tmp53;
        const REAL FDPart3tmp218 = FDPart3tmp63 * hDD_dD010 + hDD01 * (FDPart3tmp209 + FDPart3tmp37);
        const REAL FDPart3tmp258 =
            FDPart3tmp252 * (FDPart3tmp14 * hDD_dDD2211 + 4 * FDPart3tmp21 * FDPart3tmp22 * hDD_dD221 + FDPart3tmp257 * hDD22 + FDPart3tmp257);
        const REAL FDPart3tmp260 = FDPart3tmp246 * FDPart3tmp259 + FDPart3tmp98 * aDD11;
        const REAL FDPart3tmp261 = FDPart3tmp15 * FDPart3tmp20;
        const REAL FDPart3tmp266 = FDPart3tmp20 * FDPart3tmp38;
        const REAL FDPart3tmp296 =
            FDPart3tmp252 * (FDPart3tmp14 * hDD_dDD2200 + FDPart3tmp295 * hDD22 + FDPart3tmp295 + 2 * FDPart3tmp58 * hDD_dD220);
        const REAL FDPart3tmp330 = FDPart3tmp103 * FDPart3tmp246 + FDPart3tmp90 * aDD01;
        const REAL FDPart3tmp364 = FDPart3tmp252 * (2 * FDPart3tmp21 * FDPart3tmp63 * hDD_dD021 - FDPart3tmp65 + FDPart3tmp73 * hDD_dDD0211);
        const REAL FDPart3tmp399 = FDPart3tmp252 * (FDPart3tmp22 * hDD_dDD1200 + FDPart3tmp334 * hDD_dD120 +
                                                    hDD12 * (FDPart3tmp293 * FDPart3tmp3 + FDPart3tmp294 * FDPart3tmp3));
        const REAL FDPart3tmp401 = FDPart3tmp252 * (FDPart3tmp21 * FDPart3tmp63 * hDD_dD020 + FDPart3tmp210 * hDD_dD021 + FDPart3tmp73 * hDD_dDD0201 +
                                                    hDD02 * (FDPart3tmp209 * FDPart3tmp21 + FDPart3tmp21 * FDPart3tmp37));
        const REAL FDPart3tmp435 = FDPart3tmp20 * (FDPart3tmp63 * hDD_dDD0101 + hDD_dD011 * (FDPart3tmp209 + FDPart3tmp37));
        const REAL FDPart3tmp439 = FDPart3tmp252 * (FDPart3tmp13 * hDD_dDD1100 + FDPart3tmp293 + FDPart3tmp294 + 2 * FDPart3tmp57 * hDD_dD110 +
                                                    hDD11 * (FDPart3tmp293 + FDPart3tmp294));
        const REAL FDPart3tmp25 = FDPart3tmp14 * hDD_dD221 + FDPart3tmp24 * hDD22 + FDPart3tmp24;
        const REAL FDPart3tmp50 =
            (1.0 / 2.0) /
            (2 * pow(AMPL, 6) * FDPart3tmp31 * FDPart3tmp40 * FDPart3tmp41 * hDD01 * hDD02 * hDD12 / pow(FDPart3tmp11, 6) -
             FDPart3tmp15 * FDPart3tmp34 * FDPart3tmp40 * FDPart3tmp48 + FDPart3tmp15 * FDPart3tmp38 * FDPart3tmp40 * FDPart3tmp45 -
             FDPart3tmp34 * FDPart3tmp38 * FDPart3tmp40 * FDPart3tmp43 - FDPart3tmp34 * FDPart3tmp4 * FDPart3tmp40 * FDPart3tmp45 * FDPart3tmp46);
        const REAL FDPart3tmp61 = FDPart3tmp15 * FDPart3tmp53 - FDPart3tmp20 * FDPart3tmp59;
        const REAL FDPart3tmp77 = -FDPart3tmp18 * FDPart3tmp65 + FDPart3tmp20 * FDPart3tmp74;
        const REAL FDPart3tmp110 = -FDPart3tmp15 * FDPart3tmp53 + FDPart3tmp20 * FDPart3tmp59;
        const REAL FDPart3tmp111 = FDPart3tmp20 * FDPart3tmp70 - FDPart3tmp39 * FDPart3tmp53;
        const REAL FDPart3tmp169 = FDPart3tmp167 - FDPart3tmp18 * FDPart3tmp38;
        const REAL FDPart3tmp212 = FDPart3tmp143 * FDPart3tmp211 - 4 * FDPart3tmp52 * FDPart3tmp65;
        const REAL FDPart3tmp219 = FDPart3tmp143 * FDPart3tmp218 - FDPart3tmp167 + FDPart3tmp18 * FDPart3tmp38 - 4 * FDPart3tmp52 * FDPart3tmp84;
        const REAL FDPart3tmp247 = FDPart3tmp100 * FDPart3tmp246 + FDPart3tmp3 * FDPart3tmp98 * aDD12;
        const REAL FDPart3tmp262 = FDPart3tmp14 * FDPart3tmp20 * aDD22 + FDPart3tmp246 * FDPart3tmp261;
        const REAL FDPart3tmp288 = FDPart3tmp246 * FDPart3tmp92 + FDPart3tmp3 * FDPart3tmp90 * aDD02;
        const REAL FDPart3tmp291 = FDPart3tmp20 * FDPart3tmp37 * aDD00 + FDPart3tmp246 * FDPart3tmp266;
        const REAL FDPart3tmp336 = FDPart3tmp252 * (FDPart3tmp14 * hDD_dDD2201 + FDPart3tmp21 * FDPart3tmp334 * hDD22 + FDPart3tmp21 * FDPart3tmp334 +
                                                    FDPart3tmp24 * hDD_dD220 + FDPart3tmp58 * hDD_dD221);
        const REAL FDPart3tmp51 = FDPart3tmp50 * (FDPart3tmp32 * FDPart3tmp36 * hDD02 - FDPart3tmp33 * FDPart3tmp38 * FDPart3tmp39);
        const REAL FDPart3tmp66 = FDPart3tmp50 * (-FDPart3tmp33 * FDPart3tmp45 * FDPart3tmp65 + FDPart3tmp36 * FDPart3tmp62 * hDD12);
        const REAL FDPart3tmp78 = -FDPart3tmp20 * FDPart3tmp70 + FDPart3tmp39 * FDPart3tmp53 + FDPart3tmp77;
        const REAL FDPart3tmp80 = FDPart3tmp50 * (-FDPart3tmp35 * FDPart3tmp43 + FDPart3tmp45 * FDPart3tmp79);
        const REAL FDPart3tmp85 = FDPart3tmp50 * (FDPart3tmp35 * FDPart3tmp4 * FDPart3tmp62 * hDD02 * hDD12 - FDPart3tmp79 * FDPart3tmp84);
        const REAL FDPart3tmp87 = FDPart3tmp50 * (FDPart3tmp33 * FDPart3tmp38 * FDPart3tmp45 - FDPart3tmp35 * FDPart3tmp48);
        const REAL FDPart3tmp95 = FDPart3tmp50 * (-FDPart3tmp35 * FDPart3tmp4 * FDPart3tmp46 + FDPart3tmp38 * FDPart3tmp79);
        const REAL FDPart3tmp112 = FDPart3tmp111 + FDPart3tmp18 * FDPart3tmp65 - FDPart3tmp20 * FDPart3tmp74;
        const REAL FDPart3tmp171 = FDPart3tmp111 + FDPart3tmp77;
        const REAL FDPart3tmp27 = FDPart3tmp15 * FDPart3tmp18 - FDPart3tmp20 * FDPart3tmp25;
        const REAL FDPart3tmp68 = -FDPart3tmp15 * FDPart3tmp18 + FDPart3tmp20 * FDPart3tmp25;
        const REAL FDPart3tmp113 = FDPart3tmp110 * FDPart3tmp66 + FDPart3tmp112 * FDPart3tmp85;
        const REAL FDPart3tmp115 = FDPart3tmp110 * FDPart3tmp87 + FDPart3tmp112 * FDPart3tmp51;
        const REAL FDPart3tmp117 = FDPart3tmp110 * FDPart3tmp51 + FDPart3tmp112 * FDPart3tmp95;
        const REAL FDPart3tmp151 = FDPart3tmp140 * FDPart3tmp85 + FDPart3tmp144 * FDPart3tmp66 + FDPart3tmp150 * FDPart3tmp80;
        const REAL FDPart3tmp153 = FDPart3tmp140 * FDPart3tmp95 + FDPart3tmp144 * FDPart3tmp51 + FDPart3tmp150 * FDPart3tmp85;
        const REAL FDPart3tmp155 = FDPart3tmp140 * FDPart3tmp51 + FDPart3tmp144 * FDPart3tmp87 + FDPart3tmp150 * FDPart3tmp66;
        const REAL FDPart3tmp172 = FDPart3tmp169 * FDPart3tmp80 + FDPart3tmp170 * FDPart3tmp85 + FDPart3tmp171 * FDPart3tmp66;
        const REAL FDPart3tmp175 = FDPart3tmp169 * FDPart3tmp85 + FDPart3tmp170 * FDPart3tmp95 + FDPart3tmp171 * FDPart3tmp51;
        const REAL FDPart3tmp178 = FDPart3tmp169 * FDPart3tmp66 + FDPart3tmp170 * FDPart3tmp51 + FDPart3tmp171 * FDPart3tmp87;
        const REAL FDPart3tmp220 = FDPart3tmp212 * FDPart3tmp87 + FDPart3tmp216 * FDPart3tmp66 + FDPart3tmp219 * FDPart3tmp51;
        const REAL FDPart3tmp222 = FDPart3tmp212 * FDPart3tmp66 + FDPart3tmp216 * FDPart3tmp80 + FDPart3tmp219 * FDPart3tmp85;
        const REAL FDPart3tmp225 = FDPart3tmp212 * FDPart3tmp51 + FDPart3tmp216 * FDPart3tmp85 + FDPart3tmp219 * FDPart3tmp95;
        const REAL FDPart3tmp67 = FDPart3tmp27 * FDPart3tmp51 + FDPart3tmp61 * FDPart3tmp66;
        const REAL FDPart3tmp81 = FDPart3tmp66 * FDPart3tmp68 + FDPart3tmp78 * FDPart3tmp80;
        const REAL FDPart3tmp86 = FDPart3tmp27 * FDPart3tmp85 + FDPart3tmp61 * FDPart3tmp80;
        const REAL FDPart3tmp88 = FDPart3tmp66 * FDPart3tmp78 + FDPart3tmp68 * FDPart3tmp87;
        const REAL FDPart3tmp93 = FDPart3tmp51 * FDPart3tmp68 + FDPart3tmp78 * FDPart3tmp85;
        const REAL FDPart3tmp96 = FDPart3tmp27 * FDPart3tmp95 + FDPart3tmp61 * FDPart3tmp85;
        const REAL FDPart3tmp191 = FDPart3tmp113 * FDPart3tmp178;
        const REAL FDPart3tmp192 = FDPart3tmp115 * FDPart3tmp172;
        const REAL FDPart3tmp195 = FDPart3tmp117 * FDPart3tmp172;
        const REAL FDPart3tmp197 = FDPart3tmp113 * FDPart3tmp175;
        const REAL FDPart3tmp199 = FDPart3tmp117 * FDPart3tmp178;
        const REAL FDPart3tmp200 = FDPart3tmp115 * FDPart3tmp175;
        const REAL FDPart3tmp303 = FDPart3tmp113 * FDPart3tmp115;
        const REAL FDPart3tmp307 = FDPart3tmp115 * FDPart3tmp117;
        const REAL FDPart3tmp310 = FDPart3tmp113 * FDPart3tmp117;
        const REAL FDPart3tmp356 = FDPart3tmp117 * FDPart3tmp151;
        const REAL FDPart3tmp357 = FDPart3tmp113 * FDPart3tmp153;
        const REAL FDPart3tmp358 = FDPart3tmp115 * FDPart3tmp153;
        const REAL FDPart3tmp359 = FDPart3tmp117 * FDPart3tmp155;
        const REAL FDPart3tmp360 = FDPart3tmp115 * FDPart3tmp151;
        const REAL FDPart3tmp361 = FDPart3tmp113 * FDPart3tmp155;
        const REAL FDPart3tmp446 = FDPart3tmp172 * FDPart3tmp178;
        const REAL FDPart3tmp449 = FDPart3tmp172 * FDPart3tmp175;
        const REAL FDPart3tmp452 = FDPart3tmp175 * FDPart3tmp178;
        const REAL FDPart3tmp125 = FDPart3tmp113 * FDPart3tmp93;
        const REAL FDPart3tmp126 = FDPart3tmp117 * FDPart3tmp81;
        const REAL FDPart3tmp129 = FDPart3tmp115 * FDPart3tmp81;
        const REAL FDPart3tmp131 = FDPart3tmp113 * FDPart3tmp88;
        const REAL FDPart3tmp133 = FDPart3tmp115 * FDPart3tmp93;
        const REAL FDPart3tmp135 = FDPart3tmp117 * FDPart3tmp88;
        const REAL FDPart3tmp173 = FDPart3tmp172 * FDPart3tmp93;
        const REAL FDPart3tmp176 = FDPart3tmp175 * FDPart3tmp81;
        const REAL FDPart3tmp179 = FDPart3tmp178 * FDPart3tmp81;
        const REAL FDPart3tmp180 = FDPart3tmp172 * FDPart3tmp88;
        const REAL FDPart3tmp183 = FDPart3tmp178 * FDPart3tmp93;
        const REAL FDPart3tmp184 = FDPart3tmp175 * FDPart3tmp88;
        const REAL FDPart3tmp229 = FDPart3tmp100 * (-FDPart3tmp115 * FDPart3tmp225 + FDPart3tmp117 * FDPart3tmp220) +
                                   FDPart3tmp100 * (FDPart3tmp115 * FDPart3tmp225 - FDPart3tmp117 * FDPart3tmp220) +
                                   FDPart3tmp103 * (-FDPart3tmp113 * FDPart3tmp225 + FDPart3tmp117 * FDPart3tmp222) +
                                   FDPart3tmp103 * (FDPart3tmp113 * FDPart3tmp225 - FDPart3tmp117 * FDPart3tmp222) +
                                   FDPart3tmp92 * (-FDPart3tmp113 * FDPart3tmp220 + FDPart3tmp115 * FDPart3tmp222) +
                                   FDPart3tmp92 * (FDPart3tmp113 * FDPart3tmp220 - FDPart3tmp115 * FDPart3tmp222);
        const REAL FDPart3tmp238 = FDPart3tmp100 * (-FDPart3tmp153 * FDPart3tmp178 + FDPart3tmp155 * FDPart3tmp175) +
                                   FDPart3tmp100 * (FDPart3tmp153 * FDPart3tmp178 - FDPart3tmp155 * FDPart3tmp175) +
                                   FDPart3tmp103 * (-FDPart3tmp151 * FDPart3tmp175 + FDPart3tmp153 * FDPart3tmp172) +
                                   FDPart3tmp103 * (FDPart3tmp151 * FDPart3tmp175 - FDPart3tmp153 * FDPart3tmp172) +
                                   FDPart3tmp92 * (-FDPart3tmp151 * FDPart3tmp178 + FDPart3tmp155 * FDPart3tmp172) +
                                   FDPart3tmp92 * (FDPart3tmp151 * FDPart3tmp178 - FDPart3tmp155 * FDPart3tmp172);
        const REAL FDPart3tmp245 = FDPart3tmp100 * (-FDPart3tmp175 * FDPart3tmp220 + FDPart3tmp178 * FDPart3tmp225) +
                                   FDPart3tmp100 * (FDPart3tmp175 * FDPart3tmp220 - FDPart3tmp178 * FDPart3tmp225) +
                                   FDPart3tmp103 * (-FDPart3tmp172 * FDPart3tmp225 + FDPart3tmp175 * FDPart3tmp222) +
                                   FDPart3tmp103 * (FDPart3tmp172 * FDPart3tmp225 - FDPart3tmp175 * FDPart3tmp222) +
                                   FDPart3tmp92 * (-FDPart3tmp172 * FDPart3tmp220 + FDPart3tmp178 * FDPart3tmp222) +
                                   FDPart3tmp92 * (FDPart3tmp172 * FDPart3tmp220 - FDPart3tmp178 * FDPart3tmp222);
        const REAL FDPart3tmp271 = FDPart3tmp81 * FDPart3tmp88;
        const REAL FDPart3tmp274 = FDPart3tmp88 * FDPart3tmp93;
        const REAL FDPart3tmp277 = FDPart3tmp81 * FDPart3tmp93;
        const REAL FDPart3tmp318 = FDPart3tmp175 * FDPart3tmp86;
        const REAL FDPart3tmp320 = FDPart3tmp172 * FDPart3tmp96;
        const REAL FDPart3tmp322 = FDPart3tmp178 * FDPart3tmp96;
        const REAL FDPart3tmp324 = FDPart3tmp175 * FDPart3tmp67;
        const REAL FDPart3tmp326 = FDPart3tmp178 * FDPart3tmp86;
        const REAL FDPart3tmp328 = FDPart3tmp172 * FDPart3tmp67;
        const REAL FDPart3tmp392 = FDPart3tmp222 * FDPart3tmp93;
        const REAL FDPart3tmp393 = FDPart3tmp225 * FDPart3tmp81;
        const REAL FDPart3tmp394 = FDPart3tmp225 * FDPart3tmp88;
        const REAL FDPart3tmp395 = FDPart3tmp220 * FDPart3tmp93;
        const REAL FDPart3tmp396 = FDPart3tmp222 * FDPart3tmp88;
        const REAL FDPart3tmp397 = FDPart3tmp220 * FDPart3tmp81;
        const REAL FDPart3tmp455 =
            FDPart3tmp100 * (-FDPart3tmp153 * FDPart3tmp220 + FDPart3tmp452) + FDPart3tmp100 * (-FDPart3tmp155 * FDPart3tmp225 + FDPart3tmp452) +
            FDPart3tmp103 * (-FDPart3tmp151 * FDPart3tmp225 + FDPart3tmp449) + FDPart3tmp103 * (-FDPart3tmp153 * FDPart3tmp222 + FDPart3tmp449) +
            FDPart3tmp148 * FDPart3tmp53 - FDPart3tmp18 * FDPart3tmp218 - FDPart3tmp254 * FDPart3tmp38 +
            FDPart3tmp259 * (-FDPart3tmp153 * FDPart3tmp225 + ((FDPart3tmp175) * (FDPart3tmp175))) + FDPart3tmp260 * FDPart3tmp291 +
            FDPart3tmp261 * (-FDPart3tmp155 * FDPart3tmp220 + ((FDPart3tmp178) * (FDPart3tmp178))) +
            FDPart3tmp266 * (-FDPart3tmp151 * FDPart3tmp222 + ((FDPart3tmp172) * (FDPart3tmp172))) - FDPart3tmp286 * FDPart3tmp45 -
            ((FDPart3tmp330) * (FDPart3tmp330)) + FDPart3tmp338 * FDPart3tmp84 - FDPart3tmp427 + FDPart3tmp428 - FDPart3tmp429 + FDPart3tmp435 -
            FDPart3tmp439 + FDPart3tmp92 * (-FDPart3tmp151 * FDPart3tmp220 + FDPart3tmp446) +
            FDPart3tmp92 * (-FDPart3tmp155 * FDPart3tmp222 + FDPart3tmp446);
        const REAL FDPart3tmp104 = FDPart3tmp100 * (-FDPart3tmp67 * FDPart3tmp93 + FDPart3tmp88 * FDPart3tmp96) +
                                   FDPart3tmp100 * (FDPart3tmp67 * FDPart3tmp93 - FDPart3tmp88 * FDPart3tmp96) +
                                   FDPart3tmp103 * (-FDPart3tmp81 * FDPart3tmp96 + FDPart3tmp86 * FDPart3tmp93) +
                                   FDPart3tmp103 * (FDPart3tmp81 * FDPart3tmp96 - FDPart3tmp86 * FDPart3tmp93) +
                                   FDPart3tmp92 * (-FDPart3tmp67 * FDPart3tmp81 + FDPart3tmp86 * FDPart3tmp88) +
                                   FDPart3tmp92 * (FDPart3tmp67 * FDPart3tmp81 - FDPart3tmp86 * FDPart3tmp88);
        const REAL FDPart3tmp122 = FDPart3tmp100 * (-FDPart3tmp115 * FDPart3tmp96 + FDPart3tmp117 * FDPart3tmp67) +
                                   FDPart3tmp100 * (FDPart3tmp115 * FDPart3tmp96 - FDPart3tmp117 * FDPart3tmp67) +
                                   FDPart3tmp103 * (-FDPart3tmp113 * FDPart3tmp96 + FDPart3tmp117 * FDPart3tmp86) +
                                   FDPart3tmp103 * (FDPart3tmp113 * FDPart3tmp96 - FDPart3tmp117 * FDPart3tmp86) +
                                   FDPart3tmp92 * (-FDPart3tmp113 * FDPart3tmp67 + FDPart3tmp115 * FDPart3tmp86) +
                                   FDPart3tmp92 * (FDPart3tmp113 * FDPart3tmp67 - FDPart3tmp115 * FDPart3tmp86);
        const REAL FDPart3tmp160 = FDPart3tmp100 * (-FDPart3tmp153 * FDPart3tmp88 + FDPart3tmp155 * FDPart3tmp93) +
                                   FDPart3tmp100 * (FDPart3tmp153 * FDPart3tmp88 - FDPart3tmp155 * FDPart3tmp93) +
                                   FDPart3tmp103 * (-FDPart3tmp151 * FDPart3tmp93 + FDPart3tmp153 * FDPart3tmp81) +
                                   FDPart3tmp103 * (FDPart3tmp151 * FDPart3tmp93 - FDPart3tmp153 * FDPart3tmp81) +
                                   FDPart3tmp92 * (-FDPart3tmp151 * FDPart3tmp88 + FDPart3tmp155 * FDPart3tmp81) +
                                   FDPart3tmp92 * (FDPart3tmp151 * FDPart3tmp88 - FDPart3tmp155 * FDPart3tmp81);
        const REAL FDPart3tmp203 = FDPart3tmp100 * (-FDPart3tmp199 + FDPart3tmp200) + FDPart3tmp100 * (FDPart3tmp199 - FDPart3tmp200) +
                                   FDPart3tmp103 * (-FDPart3tmp195 + FDPart3tmp197) + FDPart3tmp103 * (FDPart3tmp195 - FDPart3tmp197) +
                                   FDPart3tmp92 * (-FDPart3tmp191 + FDPart3tmp192) + FDPart3tmp92 * (FDPart3tmp191 - FDPart3tmp192);
        const REAL FDPart3tmp284 =
            FDPart3tmp100 * (-FDPart3tmp153 * FDPart3tmp67 + FDPart3tmp274) + FDPart3tmp100 * (-FDPart3tmp155 * FDPart3tmp96 + FDPart3tmp274) +
            FDPart3tmp103 * (-FDPart3tmp151 * FDPart3tmp96 + FDPart3tmp277) + FDPart3tmp103 * (-FDPart3tmp153 * FDPart3tmp86 + FDPart3tmp277) -
            FDPart3tmp15 * FDPart3tmp254 + FDPart3tmp18 * FDPart3tmp25 - ((FDPart3tmp247) * (FDPart3tmp247)) - FDPart3tmp258 +
            FDPart3tmp259 * (-FDPart3tmp153 * FDPart3tmp96 + ((FDPart3tmp93) * (FDPart3tmp93))) + FDPart3tmp260 * FDPart3tmp262 +
            FDPart3tmp261 * (-FDPart3tmp155 * FDPart3tmp67 + ((FDPart3tmp88) * (FDPart3tmp88))) +
            FDPart3tmp266 * (-FDPart3tmp151 * FDPart3tmp86 + ((FDPart3tmp81) * (FDPart3tmp81))) +
            FDPart3tmp92 * (-FDPart3tmp151 * FDPart3tmp67 + FDPart3tmp271) + FDPart3tmp92 * (-FDPart3tmp155 * FDPart3tmp86 + FDPart3tmp271);
        const REAL FDPart3tmp315 =
            2 * FDPart3tmp100 * (FDPart3tmp220 * FDPart3tmp96 - FDPart3tmp307) + 2 * FDPart3tmp100 * (FDPart3tmp225 * FDPart3tmp67 - FDPart3tmp307) +
            2 * FDPart3tmp103 * (FDPart3tmp222 * FDPart3tmp96 - FDPart3tmp310) + 2 * FDPart3tmp103 * (FDPart3tmp225 * FDPart3tmp86 - FDPart3tmp310) +
            2 * FDPart3tmp15 * FDPart3tmp286 + 2 * FDPart3tmp259 * (-((FDPart3tmp117) * (FDPart3tmp117)) + FDPart3tmp225 * FDPart3tmp96) +
            2 * FDPart3tmp261 * (-((FDPart3tmp115) * (FDPart3tmp115)) + FDPart3tmp220 * FDPart3tmp67) - 2 * FDPart3tmp262 * FDPart3tmp291 +
            2 * FDPart3tmp266 * (-((FDPart3tmp113) * (FDPart3tmp113)) + FDPart3tmp222 * FDPart3tmp86) + 2 * ((FDPart3tmp288) * (FDPart3tmp288)) +
            2 * FDPart3tmp296 - 2 * FDPart3tmp53 * FDPart3tmp59 + 2 * FDPart3tmp92 * (FDPart3tmp220 * FDPart3tmp86 - FDPart3tmp303) +
            2 * FDPart3tmp92 * (FDPart3tmp222 * FDPart3tmp67 - FDPart3tmp303);
        const REAL FDPart3tmp316 =
            FDPart3tmp100 * (-FDPart3tmp220 * FDPart3tmp96 + FDPart3tmp307) + FDPart3tmp100 * (-FDPart3tmp225 * FDPart3tmp67 + FDPart3tmp307) +
            FDPart3tmp103 * (-FDPart3tmp222 * FDPart3tmp96 + FDPart3tmp310) + FDPart3tmp103 * (-FDPart3tmp225 * FDPart3tmp86 + FDPart3tmp310) -
            FDPart3tmp15 * FDPart3tmp286 + FDPart3tmp259 * (((FDPart3tmp117) * (FDPart3tmp117)) - FDPart3tmp225 * FDPart3tmp96) +
            FDPart3tmp261 * (((FDPart3tmp115) * (FDPart3tmp115)) - FDPart3tmp220 * FDPart3tmp67) + FDPart3tmp262 * FDPart3tmp291 +
            FDPart3tmp266 * (((FDPart3tmp113) * (FDPart3tmp113)) - FDPart3tmp222 * FDPart3tmp86) - ((FDPart3tmp288) * (FDPart3tmp288)) -
            FDPart3tmp296 + FDPart3tmp53 * FDPart3tmp59 + FDPart3tmp92 * (-FDPart3tmp220 * FDPart3tmp86 + FDPart3tmp303) +
            FDPart3tmp92 * (-FDPart3tmp222 * FDPart3tmp67 + FDPart3tmp303);
        const REAL FDPart3tmp349 = -FDPart3tmp15 * FDPart3tmp339 + FDPart3tmp17 * FDPart3tmp59 - FDPart3tmp247 * FDPart3tmp288 +
                                   FDPart3tmp25 * FDPart3tmp52 + FDPart3tmp259 * (FDPart3tmp117 * FDPart3tmp93 - FDPart3tmp175 * FDPart3tmp96) +
                                   FDPart3tmp261 * (FDPart3tmp115 * FDPart3tmp88 - FDPart3tmp178 * FDPart3tmp67) + FDPart3tmp262 * FDPart3tmp330 +
                                   FDPart3tmp266 * (FDPart3tmp113 * FDPart3tmp81 - FDPart3tmp172 * FDPart3tmp86) - FDPart3tmp336;
        const REAL FDPart3tmp353 = FDPart3tmp15 * FDPart3tmp339 - FDPart3tmp17 * FDPart3tmp59 + FDPart3tmp247 * FDPart3tmp288 -
                                   FDPart3tmp25 * FDPart3tmp52 + FDPart3tmp259 * (-FDPart3tmp117 * FDPart3tmp93 + FDPart3tmp175 * FDPart3tmp96) +
                                   FDPart3tmp261 * (-FDPart3tmp115 * FDPart3tmp88 + FDPart3tmp178 * FDPart3tmp67) - FDPart3tmp262 * FDPart3tmp330 +
                                   FDPart3tmp266 * (-FDPart3tmp113 * FDPart3tmp81 + FDPart3tmp172 * FDPart3tmp86) + FDPart3tmp336;
        const REAL FDPart3tmp378 = FDPart3tmp142 * FDPart3tmp52 + FDPart3tmp17 * FDPart3tmp70 - FDPart3tmp18 * FDPart3tmp74 +
                                   FDPart3tmp247 * FDPart3tmp330 + FDPart3tmp254 * FDPart3tmp65 +
                                   FDPart3tmp259 * (FDPart3tmp117 * FDPart3tmp153 - FDPart3tmp175 * FDPart3tmp93) - FDPart3tmp260 * FDPart3tmp288 +
                                   FDPart3tmp261 * (FDPart3tmp115 * FDPart3tmp155 - FDPart3tmp178 * FDPart3tmp88) +
                                   FDPart3tmp266 * (FDPart3tmp113 * FDPart3tmp151 - FDPart3tmp172 * FDPart3tmp81) - FDPart3tmp339 * FDPart3tmp39 +
                                   FDPart3tmp364 - FDPart3tmp366;
        const REAL FDPart3tmp388 = -FDPart3tmp142 * FDPart3tmp52 - FDPart3tmp17 * FDPart3tmp70 + FDPart3tmp18 * FDPart3tmp74 -
                                   FDPart3tmp247 * FDPart3tmp330 - FDPart3tmp254 * FDPart3tmp65 +
                                   FDPart3tmp259 * (-FDPart3tmp117 * FDPart3tmp153 + FDPart3tmp175 * FDPart3tmp93) + FDPart3tmp260 * FDPart3tmp288 +
                                   FDPart3tmp261 * (-FDPart3tmp115 * FDPart3tmp155 + FDPart3tmp178 * FDPart3tmp88) +
                                   FDPart3tmp266 * (-FDPart3tmp113 * FDPart3tmp151 + FDPart3tmp172 * FDPart3tmp81) + FDPart3tmp339 * FDPart3tmp39 -
                                   FDPart3tmp364 + FDPart3tmp366;
        const REAL FDPart3tmp413 = FDPart3tmp17 * FDPart3tmp211 - FDPart3tmp247 * FDPart3tmp291 +
                                   FDPart3tmp259 * (-FDPart3tmp117 * FDPart3tmp175 + FDPart3tmp225 * FDPart3tmp93) +
                                   FDPart3tmp261 * (-FDPart3tmp115 * FDPart3tmp178 + FDPart3tmp220 * FDPart3tmp88) +
                                   FDPart3tmp266 * (-FDPart3tmp113 * FDPart3tmp172 + FDPart3tmp222 * FDPart3tmp81) + FDPart3tmp286 * FDPart3tmp39 +
                                   FDPart3tmp288 * FDPart3tmp330 - FDPart3tmp339 * FDPart3tmp65 + FDPart3tmp399 - FDPart3tmp401 +
                                   FDPart3tmp52 * FDPart3tmp74 - FDPart3tmp53 * FDPart3tmp70;
        const REAL FDPart3tmp423 = -FDPart3tmp17 * FDPart3tmp211 + FDPart3tmp247 * FDPart3tmp291 +
                                   FDPart3tmp259 * (FDPart3tmp117 * FDPart3tmp175 - FDPart3tmp225 * FDPart3tmp93) +
                                   FDPart3tmp261 * (FDPart3tmp115 * FDPart3tmp178 - FDPart3tmp220 * FDPart3tmp88) +
                                   FDPart3tmp266 * (FDPart3tmp113 * FDPart3tmp172 - FDPart3tmp222 * FDPart3tmp81) - FDPart3tmp286 * FDPart3tmp39 -
                                   FDPart3tmp288 * FDPart3tmp330 + FDPart3tmp339 * FDPart3tmp65 - FDPart3tmp399 + FDPart3tmp401 -
                                   FDPart3tmp52 * FDPart3tmp74 + FDPart3tmp53 * FDPart3tmp70;
        const REAL FDPart3tmp459 =
            2 * FDPart3tmp100 * (FDPart3tmp153 * FDPart3tmp220 - FDPart3tmp452) +
            2 * FDPart3tmp100 * (FDPart3tmp155 * FDPart3tmp225 - FDPart3tmp452) +
            2 * FDPart3tmp103 * (FDPart3tmp151 * FDPart3tmp225 - FDPart3tmp449) +
            2 * FDPart3tmp103 * (FDPart3tmp153 * FDPart3tmp222 - FDPart3tmp449) - 2 * FDPart3tmp148 * FDPart3tmp53 +
            2 * FDPart3tmp18 * FDPart3tmp218 + 2 * FDPart3tmp254 * FDPart3tmp38 +
            2 * FDPart3tmp259 * (FDPart3tmp153 * FDPart3tmp225 - ((FDPart3tmp175) * (FDPart3tmp175))) - 2 * FDPart3tmp260 * FDPart3tmp291 +
            2 * FDPart3tmp261 * (FDPart3tmp155 * FDPart3tmp220 - ((FDPart3tmp178) * (FDPart3tmp178))) +
            2 * FDPart3tmp266 * (FDPart3tmp151 * FDPart3tmp222 - ((FDPart3tmp172) * (FDPart3tmp172))) + 2 * FDPart3tmp286 * FDPart3tmp45 +
            2 * ((FDPart3tmp330) * (FDPart3tmp330)) - 2 * FDPart3tmp338 * FDPart3tmp84 + 2 * FDPart3tmp427 - 2 * FDPart3tmp428 + 2 * FDPart3tmp429 -
            2 * FDPart3tmp435 + 2 * FDPart3tmp439 + 2 * FDPart3tmp92 * (FDPart3tmp151 * FDPart3tmp220 - FDPart3tmp446) +
            2 * FDPart3tmp92 * (FDPart3tmp155 * FDPart3tmp222 - FDPart3tmp446);
        const REAL FDPart3tmp137 = FDPart3tmp100 * (-FDPart3tmp133 + FDPart3tmp135) + FDPart3tmp100 * (FDPart3tmp133 - FDPart3tmp135) +
                                   FDPart3tmp103 * (-FDPart3tmp125 + FDPart3tmp126) + FDPart3tmp103 * (FDPart3tmp125 - FDPart3tmp126) +
                                   FDPart3tmp92 * (-FDPart3tmp129 + FDPart3tmp131) + FDPart3tmp92 * (FDPart3tmp129 - FDPart3tmp131);
        const REAL FDPart3tmp187 = FDPart3tmp100 * (-FDPart3tmp183 + FDPart3tmp184) + FDPart3tmp100 * (FDPart3tmp183 - FDPart3tmp184) +
                                   FDPart3tmp103 * (-FDPart3tmp173 + FDPart3tmp176) + FDPart3tmp103 * (FDPart3tmp173 - FDPart3tmp176) +
                                   FDPart3tmp92 * (-FDPart3tmp179 + FDPart3tmp180) + FDPart3tmp92 * (FDPart3tmp179 - FDPart3tmp180);
        const REAL FDPart3tmp283 =
            2 * FDPart3tmp100 * (FDPart3tmp153 * FDPart3tmp67 - FDPart3tmp274) + 2 * FDPart3tmp100 * (FDPart3tmp155 * FDPart3tmp96 - FDPart3tmp274) +
            2 * FDPart3tmp103 * (FDPart3tmp151 * FDPart3tmp96 - FDPart3tmp277) + 2 * FDPart3tmp103 * (FDPart3tmp153 * FDPart3tmp86 - FDPart3tmp277) +
            2 * FDPart3tmp15 * FDPart3tmp254 - 2 * FDPart3tmp18 * FDPart3tmp25 + 2 * ((FDPart3tmp247) * (FDPart3tmp247)) + 2 * FDPart3tmp258 +
            2 * FDPart3tmp259 * (FDPart3tmp153 * FDPart3tmp96 - ((FDPart3tmp93) * (FDPart3tmp93))) - 2 * FDPart3tmp260 * FDPart3tmp262 +
            2 * FDPart3tmp261 * (FDPart3tmp155 * FDPart3tmp67 - ((FDPart3tmp88) * (FDPart3tmp88))) +
            2 * FDPart3tmp266 * (FDPart3tmp151 * FDPart3tmp86 - ((FDPart3tmp81) * (FDPart3tmp81))) +
            2 * FDPart3tmp92 * (FDPart3tmp151 * FDPart3tmp67 - FDPart3tmp271) + 2 * FDPart3tmp92 * (FDPart3tmp155 * FDPart3tmp86 - FDPart3tmp271);
        const REAL FDPart3tmp350 = FDPart3tmp100 * (FDPart3tmp133 - FDPart3tmp322) + FDPart3tmp100 * (FDPart3tmp135 - FDPart3tmp324) +
                                   FDPart3tmp103 * (FDPart3tmp125 - FDPart3tmp320) + FDPart3tmp103 * (FDPart3tmp126 - FDPart3tmp318) + FDPart3tmp349 +
                                   FDPart3tmp92 * (FDPart3tmp129 - FDPart3tmp326) + FDPart3tmp92 * (FDPart3tmp131 - FDPart3tmp328);
        const REAL FDPart3tmp351 = FDPart3tmp100 * (FDPart3tmp133 - FDPart3tmp324) + FDPart3tmp100 * (FDPart3tmp135 - FDPart3tmp322) +
                                   FDPart3tmp103 * (FDPart3tmp125 - FDPart3tmp318) + FDPart3tmp103 * (FDPart3tmp126 - FDPart3tmp320) + FDPart3tmp349 +
                                   FDPart3tmp92 * (FDPart3tmp129 - FDPart3tmp328) + FDPart3tmp92 * (FDPart3tmp131 - FDPart3tmp326);
        const REAL FDPart3tmp354 = FDPart3tmp100 * (-FDPart3tmp133 + FDPart3tmp322) + FDPart3tmp100 * (-FDPart3tmp135 + FDPart3tmp324) +
                                   FDPart3tmp103 * (-FDPart3tmp125 + FDPart3tmp320) + FDPart3tmp103 * (-FDPart3tmp126 + FDPart3tmp318) +
                                   FDPart3tmp353 + FDPart3tmp92 * (-FDPart3tmp129 + FDPart3tmp326) + FDPart3tmp92 * (-FDPart3tmp131 + FDPart3tmp328);
        const REAL FDPart3tmp355 = FDPart3tmp100 * (-FDPart3tmp133 + FDPart3tmp324) + FDPart3tmp100 * (-FDPart3tmp135 + FDPart3tmp322) +
                                   FDPart3tmp103 * (-FDPart3tmp125 + FDPart3tmp318) + FDPart3tmp103 * (-FDPart3tmp126 + FDPart3tmp320) +
                                   FDPart3tmp353 + FDPart3tmp92 * (-FDPart3tmp129 + FDPart3tmp328) + FDPart3tmp92 * (-FDPart3tmp131 + FDPart3tmp326);
        const REAL FDPart3tmp379 = FDPart3tmp100 * (-FDPart3tmp183 + FDPart3tmp359) + FDPart3tmp100 * (-FDPart3tmp184 + FDPart3tmp358) +
                                   FDPart3tmp103 * (-FDPart3tmp173 + FDPart3tmp356) + FDPart3tmp103 * (-FDPart3tmp176 + FDPart3tmp357) +
                                   FDPart3tmp378 + FDPart3tmp92 * (-FDPart3tmp179 + FDPart3tmp361) + FDPart3tmp92 * (-FDPart3tmp180 + FDPart3tmp360);
        const REAL FDPart3tmp380 = FDPart3tmp100 * (-FDPart3tmp183 + FDPart3tmp358) + FDPart3tmp100 * (-FDPart3tmp184 + FDPart3tmp359) +
                                   FDPart3tmp103 * (-FDPart3tmp173 + FDPart3tmp357) + FDPart3tmp103 * (-FDPart3tmp176 + FDPart3tmp356) +
                                   FDPart3tmp378 + FDPart3tmp92 * (-FDPart3tmp179 + FDPart3tmp360) + FDPart3tmp92 * (-FDPart3tmp180 + FDPart3tmp361);
        const REAL FDPart3tmp389 = FDPart3tmp100 * (FDPart3tmp183 - FDPart3tmp359) + FDPart3tmp100 * (FDPart3tmp184 - FDPart3tmp358) +
                                   FDPart3tmp103 * (FDPart3tmp173 - FDPart3tmp356) + FDPart3tmp103 * (FDPart3tmp176 - FDPart3tmp357) + FDPart3tmp388 +
                                   FDPart3tmp92 * (FDPart3tmp179 - FDPart3tmp361) + FDPart3tmp92 * (FDPart3tmp180 - FDPart3tmp360);
        const REAL FDPart3tmp390 = FDPart3tmp100 * (FDPart3tmp183 - FDPart3tmp358) + FDPart3tmp100 * (FDPart3tmp184 - FDPart3tmp359) +
                                   FDPart3tmp103 * (FDPart3tmp173 - FDPart3tmp357) + FDPart3tmp103 * (FDPart3tmp176 - FDPart3tmp356) + FDPart3tmp388 +
                                   FDPart3tmp92 * (FDPart3tmp179 - FDPart3tmp360) + FDPart3tmp92 * (FDPart3tmp180 - FDPart3tmp361);
        const REAL FDPart3tmp414 = FDPart3tmp100 * (-FDPart3tmp199 + FDPart3tmp395) + FDPart3tmp100 * (-FDPart3tmp200 + FDPart3tmp394) +
                                   FDPart3tmp103 * (-FDPart3tmp195 + FDPart3tmp392) + FDPart3tmp103 * (-FDPart3tmp197 + FDPart3tmp393) +
                                   FDPart3tmp413 + FDPart3tmp92 * (-FDPart3tmp191 + FDPart3tmp397) + FDPart3tmp92 * (-FDPart3tmp192 + FDPart3tmp396);
        const REAL FDPart3tmp415 = FDPart3tmp100 * (-FDPart3tmp199 + FDPart3tmp394) + FDPart3tmp100 * (-FDPart3tmp200 + FDPart3tmp395) +
                                   FDPart3tmp103 * (-FDPart3tmp195 + FDPart3tmp393) + FDPart3tmp103 * (-FDPart3tmp197 + FDPart3tmp392) +
                                   FDPart3tmp413 + FDPart3tmp92 * (-FDPart3tmp191 + FDPart3tmp396) + FDPart3tmp92 * (-FDPart3tmp192 + FDPart3tmp397);
        const REAL FDPart3tmp424 = FDPart3tmp100 * (FDPart3tmp199 - FDPart3tmp395) + FDPart3tmp100 * (FDPart3tmp200 - FDPart3tmp394) +
                                   FDPart3tmp103 * (FDPart3tmp195 - FDPart3tmp392) + FDPart3tmp103 * (FDPart3tmp197 - FDPart3tmp393) + FDPart3tmp423 +
                                   FDPart3tmp92 * (FDPart3tmp191 - FDPart3tmp397) + FDPart3tmp92 * (FDPart3tmp192 - FDPart3tmp396);
        const REAL FDPart3tmp425 = FDPart3tmp100 * (FDPart3tmp199 - FDPart3tmp394) + FDPart3tmp100 * (FDPart3tmp200 - FDPart3tmp395) +
                                   FDPart3tmp103 * (FDPart3tmp195 - FDPart3tmp393) + FDPart3tmp103 * (FDPart3tmp197 - FDPart3tmp392) + FDPart3tmp423 +
                                   FDPart3tmp92 * (FDPart3tmp191 - FDPart3tmp396) + FDPart3tmp92 * (FDPart3tmp192 - FDPart3tmp397);
        diagnostic_output_gfs[IDX4(PSI4_PART0REGF, i0, i1, i2)] =
            FDPart3tmp104 * (FDPart3tmp0 * FDPart3tmp1 - FDPart3tmp0 * FDPart3tmp2) +
            FDPart3tmp104 * (FDPart3tmp105 * FDPart3tmp106 - FDPart3tmp106 * FDPart3tmp107) +
            FDPart3tmp122 * (FDPart3tmp0 * FDPart3tmp108 - FDPart3tmp0 * FDPart3tmp109) +
            FDPart3tmp122 * (FDPart3tmp105 * FDPart3tmp123 - FDPart3tmp107 * FDPart3tmp123) + FDPart3tmp124 * FDPart3tmp137 +
            FDPart3tmp124 * FDPart3tmp354 + FDPart3tmp124 * FDPart3tmp355 + FDPart3tmp137 * FDPart3tmp138 + FDPart3tmp138 * FDPart3tmp354 +
            FDPart3tmp138 * FDPart3tmp355 + FDPart3tmp160 * (FDPart3tmp1 * FDPart3tmp139 - FDPart3tmp139 * FDPart3tmp2) +
            FDPart3tmp160 * (FDPart3tmp106 * FDPart3tmp161 - FDPart3tmp106 * FDPart3tmp162) + FDPart3tmp165 * FDPart3tmp187 +
            FDPart3tmp165 * FDPart3tmp379 + FDPart3tmp165 * FDPart3tmp380 + FDPart3tmp187 * FDPart3tmp189 + FDPart3tmp189 * FDPart3tmp379 +
            FDPart3tmp189 * FDPart3tmp380 + FDPart3tmp190 * FDPart3tmp203 + FDPart3tmp190 * FDPart3tmp414 + FDPart3tmp190 * FDPart3tmp415 +
            FDPart3tmp203 * FDPart3tmp204 + FDPart3tmp204 * FDPart3tmp414 + FDPart3tmp204 * FDPart3tmp415 +
            FDPart3tmp229 * (FDPart3tmp108 * FDPart3tmp205 - FDPart3tmp109 * FDPart3tmp205) +
            FDPart3tmp229 * (FDPart3tmp123 * FDPart3tmp230 - FDPart3tmp123 * FDPart3tmp231) +
            FDPart3tmp238 * (FDPart3tmp139 * FDPart3tmp163 - FDPart3tmp139 * FDPart3tmp164) +
            FDPart3tmp238 * (FDPart3tmp161 * FDPart3tmp188 - FDPart3tmp162 * FDPart3tmp188) +
            FDPart3tmp245 * (FDPart3tmp163 * FDPart3tmp205 - FDPart3tmp164 * FDPart3tmp205) +
            FDPart3tmp245 * (FDPart3tmp188 * FDPart3tmp230 - FDPart3tmp188 * FDPart3tmp231) +
            FDPart3tmp283 * (FDPart3tmp1 * FDPart3tmp106 - FDPart3tmp106 * FDPart3tmp2) +
            FDPart3tmp284 * (FDPart3tmp0 * FDPart3tmp161 - FDPart3tmp0 * FDPart3tmp162) +
            FDPart3tmp284 * (FDPart3tmp105 * FDPart3tmp139 - FDPart3tmp107 * FDPart3tmp139) +
            FDPart3tmp315 * (FDPart3tmp108 * FDPart3tmp123 - FDPart3tmp109 * FDPart3tmp123) +
            FDPart3tmp316 * (FDPart3tmp0 * FDPart3tmp230 - FDPart3tmp0 * FDPart3tmp231) +
            FDPart3tmp316 * (FDPart3tmp105 * FDPart3tmp205 - FDPart3tmp107 * FDPart3tmp205) + FDPart3tmp317 * FDPart3tmp350 +
            FDPart3tmp317 * FDPart3tmp351 + FDPart3tmp350 * FDPart3tmp352 + FDPart3tmp351 * FDPart3tmp352 + FDPart3tmp381 * FDPart3tmp389 +
            FDPart3tmp381 * FDPart3tmp390 + FDPart3tmp389 * FDPart3tmp391 + FDPart3tmp390 * FDPart3tmp391 + FDPart3tmp416 * FDPart3tmp424 +
            FDPart3tmp416 * FDPart3tmp425 + FDPart3tmp424 * FDPart3tmp426 + FDPart3tmp425 * FDPart3tmp426 +
            FDPart3tmp455 * (FDPart3tmp139 * FDPart3tmp230 - FDPart3tmp139 * FDPart3tmp231) +
            FDPart3tmp455 * (FDPart3tmp161 * FDPart3tmp205 - FDPart3tmp162 * FDPart3tmp205) +
            FDPart3tmp459 * (FDPart3tmp163 * FDPart3tmp188 - FDPart3tmp164 * FDPart3tmp188);
        diagnostic_output_gfs[IDX4(PSI4_PART0IMGF, i0, i1, i2)] =
            -2 * FDPart3tmp0 * FDPart3tmp284 * FDPart3tmp467 - 2 * FDPart3tmp0 * FDPart3tmp316 * FDPart3tmp475 -
            FDPart3tmp104 * FDPart3tmp106 * FDPart3tmp460 + FDPart3tmp104 * (-FDPart3tmp0 * FDPart3tmp461 - FDPart3tmp0 * FDPart3tmp462) -
            FDPart3tmp106 * FDPart3tmp160 * FDPart3tmp468 - FDPart3tmp106 * FDPart3tmp424 * FDPart3tmp476 -
            FDPart3tmp106 * FDPart3tmp425 * FDPart3tmp476 - FDPart3tmp122 * FDPart3tmp123 * FDPart3tmp460 +
            FDPart3tmp122 * (-FDPart3tmp0 * FDPart3tmp463 - FDPart3tmp0 * FDPart3tmp464) - FDPart3tmp123 * FDPart3tmp229 * FDPart3tmp476 -
            FDPart3tmp123 * FDPart3tmp389 * FDPart3tmp468 - FDPart3tmp123 * FDPart3tmp390 * FDPart3tmp468 + FDPart3tmp137 * FDPart3tmp465 +
            FDPart3tmp137 * FDPart3tmp466 - FDPart3tmp139 * FDPart3tmp284 * FDPart3tmp460 - 2 * FDPart3tmp139 * FDPart3tmp455 * FDPart3tmp475 +
            FDPart3tmp160 * (-FDPart3tmp139 * FDPart3tmp461 - FDPart3tmp139 * FDPart3tmp462) + FDPart3tmp187 * FDPart3tmp471 +
            FDPart3tmp187 * FDPart3tmp472 - FDPart3tmp188 * FDPart3tmp238 * FDPart3tmp468 - FDPart3tmp188 * FDPart3tmp245 * FDPart3tmp476 -
            FDPart3tmp188 * FDPart3tmp350 * FDPart3tmp460 - FDPart3tmp188 * FDPart3tmp351 * FDPart3tmp460 + FDPart3tmp203 * FDPart3tmp473 +
            FDPart3tmp203 * FDPart3tmp474 - FDPart3tmp205 * FDPart3tmp316 * FDPart3tmp460 - 2 * FDPart3tmp205 * FDPart3tmp455 * FDPart3tmp467 +
            FDPart3tmp229 * (-FDPart3tmp205 * FDPart3tmp463 - FDPart3tmp205 * FDPart3tmp464) +
            FDPart3tmp238 * (-FDPart3tmp139 * FDPart3tmp469 - FDPart3tmp139 * FDPart3tmp470) +
            FDPart3tmp245 * (-FDPart3tmp205 * FDPart3tmp469 - FDPart3tmp205 * FDPart3tmp470) +
            FDPart3tmp283 * (-FDPart3tmp106 * FDPart3tmp461 - FDPart3tmp106 * FDPart3tmp462) +
            FDPart3tmp315 * (-FDPart3tmp123 * FDPart3tmp463 - FDPart3tmp123 * FDPart3tmp464) + FDPart3tmp350 * FDPart3tmp479 +
            FDPart3tmp351 * FDPart3tmp479 + FDPart3tmp354 * FDPart3tmp465 + FDPart3tmp354 * FDPart3tmp466 + FDPart3tmp355 * FDPart3tmp465 +
            FDPart3tmp355 * FDPart3tmp466 + FDPart3tmp379 * FDPart3tmp471 + FDPart3tmp379 * FDPart3tmp472 + FDPart3tmp380 * FDPart3tmp471 +
            FDPart3tmp380 * FDPart3tmp472 + FDPart3tmp389 * FDPart3tmp481 + FDPart3tmp390 * FDPart3tmp481 + FDPart3tmp414 * FDPart3tmp473 +
            FDPart3tmp414 * FDPart3tmp474 + FDPart3tmp415 * FDPart3tmp473 + FDPart3tmp415 * FDPart3tmp474 + FDPart3tmp424 * FDPart3tmp483 +
            FDPart3tmp425 * FDPart3tmp483 + FDPart3tmp459 * (-FDPart3tmp188 * FDPart3tmp469 - FDPart3tmp188 * FDPart3tmp470);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
