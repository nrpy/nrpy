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
 * Compute psi4 at all interior gridpoints, part 2
 */
void psi4_part2__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
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
        const REAL FDPart3tmp0 = (1.0 / 2.0) * mim4U0;
        const REAL FDPart3tmp1 = n4U0 * n4U3;
        const REAL FDPart3tmp3 = (1.0 / 2.0) * mre4U0;
        const REAL FDPart3tmp5 = ((mim4U0) * (mim4U0));
        const REAL FDPart3tmp10 = ((mre4U0) * (mre4U0));
        const REAL FDPart3tmp11 = (1.0 / (SINHW));
        const REAL FDPart3tmp22 = sin(xx1);
        const REAL FDPart3tmp23 = (1.0 / ((cf) * (cf) * (cf) * (cf)));
        const REAL FDPart3tmp35 = pow(cf, -6);
        const REAL FDPart3tmp49 = (1.0 / ((cf) * (cf)));
        const REAL FDPart3tmp52 = (1.0 / 3.0) * trK;
        const REAL FDPart3tmp114 = (1.0 / ((cf) * (cf) * (cf)));
        const REAL FDPart3tmp118 = cos(xx1);
        const REAL FDPart3tmp183 = (1.0 / 2.0) / cf;
        const REAL FDPart3tmp226 = (1.0 / ((SINHW) * (SINHW)));
        const REAL FDPart3tmp313 = (1.0 / 4.0) * n4U2;
        const REAL FDPart3tmp407 = (1.0 / 4.0) * n4U1;
        const REAL FDPart3tmp7 = (1.0 / 4.0) * ((n4U3) * (n4U3));
        const REAL FDPart3tmp9 = (1.0 / 4.0) * ((n4U0) * (n4U0));
        const REAL FDPart3tmp24 = exp(FDPart3tmp11) - exp(-FDPart3tmp11);
        const REAL FDPart3tmp36 = ((FDPart3tmp22) * (FDPart3tmp22));
        const REAL FDPart3tmp115 = FDPart3tmp114 * cf_dD1;
        const REAL FDPart3tmp126 = FDPart3tmp114 * cf_dD0;
        const REAL FDPart3tmp182 = 4 * FDPart3tmp23;
        const REAL FDPart3tmp184 = (1.0 / 2.0) * FDPart3tmp49;
        const REAL FDPart3tmp185 = 4 * FDPart3tmp49;
        const REAL FDPart3tmp195 = 2 * FDPart3tmp49;
        const REAL FDPart3tmp315 = FDPart3tmp313 * mim4U0 * n4U0;
        const REAL FDPart3tmp318 = FDPart3tmp313 * mre4U0 * n4U0;
        const REAL FDPart3tmp408 = FDPart3tmp407 * mim4U0 * n4U0;
        const REAL FDPart3tmp409 = FDPart3tmp407 * mre4U0 * n4U0;
        const REAL FDPart3tmp467 = FDPart3tmp0 * n4U0 * n4U2;
        const REAL FDPart3tmp468 = FDPart3tmp3 * n4U0 * n4U2;
        const REAL FDPart3tmp470 = (1.0 / 4.0) * ((n4U2) * (n4U2));
        const REAL FDPart3tmp514 = FDPart3tmp0 * n4U0 * n4U1;
        const REAL FDPart3tmp515 = FDPart3tmp3 * n4U0 * n4U1;
        const REAL FDPart3tmp517 = (1.0 / 4.0) * ((n4U1) * (n4U1));
        const REAL FDPart3tmp527 = FDPart3tmp0 * mre4U0;
        const REAL FDPart3tmp528 = (1.0 / 2.0) * ((n4U0) * (n4U0));
        const REAL FDPart3tmp13 = exp(FDPart3tmp11 * xx0);
        const REAL FDPart3tmp15 = exp(-FDPart3tmp11 * xx0);
        const REAL FDPart3tmp25 = ((AMPL) * (AMPL) * (AMPL) * (AMPL)) / ((FDPart3tmp24) * (FDPart3tmp24) * (FDPart3tmp24) * (FDPart3tmp24));
        const REAL FDPart3tmp28 = ((AMPL) * (AMPL)) / ((FDPart3tmp24) * (FDPart3tmp24));
        const REAL FDPart3tmp140 = 2 * FDPart3tmp126;
        const REAL FDPart3tmp145 = 2 * FDPart3tmp115;
        const REAL FDPart3tmp186 = (1.0 / 2.0) * FDPart3tmp182 * ((cf_dD1) * (cf_dD1)) +
                                   (1.0 / 2.0) * FDPart3tmp185 * (-FDPart3tmp183 * cf_dDD11 + FDPart3tmp184 * ((cf_dD1) * (cf_dD1)));
        const REAL FDPart3tmp223 = (1.0 / 2.0) * FDPart3tmp182 * ((cf_dD0) * (cf_dD0)) +
                                   (1.0 / 2.0) * FDPart3tmp185 * (-FDPart3tmp183 * cf_dDD00 + FDPart3tmp184 * ((cf_dD0) * (cf_dD0)));
        const REAL FDPart3tmp295 = FDPart3tmp182 * cf_dD0 * cf_dD1 + FDPart3tmp185 * (-FDPart3tmp183 * cf_dDD01 + FDPart3tmp184 * cf_dD0 * cf_dD1);
        const REAL FDPart3tmp312 = (1.0 / 4.0) * FDPart3tmp1 * mim4U0;
        const REAL FDPart3tmp316 = (1.0 / 4.0) * FDPart3tmp1 * mre4U0;
        const REAL FDPart3tmp411 = FDPart3tmp9 * mim4U1;
        const REAL FDPart3tmp412 = FDPart3tmp9 * mre4U1;
        const REAL FDPart3tmp19 = FDPart3tmp13 - FDPart3tmp15;
        const REAL FDPart3tmp26 = FDPart3tmp23 * FDPart3tmp25;
        const REAL FDPart3tmp296 = (1.0 / 2.0) * FDPart3tmp295;
        const REAL FDPart3tmp322 = FDPart3tmp10 * FDPart3tmp313 * n4U3 + FDPart3tmp312 * mim4U2 - FDPart3tmp313 * FDPart3tmp5 * n4U3 +
                                   FDPart3tmp315 * mim4U3 - FDPart3tmp316 * mre4U2 - FDPart3tmp318 * mre4U3 - FDPart3tmp9 * mim4U2 * mim4U3 +
                                   FDPart3tmp9 * mre4U2 * mre4U3;
        const REAL FDPart3tmp413 = FDPart3tmp10 * FDPart3tmp407 * n4U3 + FDPart3tmp312 * mim4U1 - FDPart3tmp316 * mre4U1 -
                                   FDPart3tmp407 * FDPart3tmp5 * n4U3 + FDPart3tmp408 * mim4U3 - FDPart3tmp409 * mre4U3 - FDPart3tmp411 * mim4U3 +
                                   FDPart3tmp412 * mre4U3;
        const REAL FDPart3tmp507 = FDPart3tmp10 * FDPart3tmp407 * n4U2 + FDPart3tmp315 * mim4U1 - FDPart3tmp318 * mre4U1 -
                                   FDPart3tmp407 * FDPart3tmp5 * n4U2 + FDPart3tmp408 * mim4U2 - FDPart3tmp409 * mre4U2 - FDPart3tmp411 * mim4U2 +
                                   FDPart3tmp412 * mre4U2;
        const REAL FDPart3tmp529 = FDPart3tmp312 * mre4U2 + FDPart3tmp315 * mre4U3 + FDPart3tmp316 * mim4U2 + FDPart3tmp318 * mim4U3 -
                                   FDPart3tmp527 * n4U2 * n4U3 - FDPart3tmp9 * mim4U2 * mre4U3 - FDPart3tmp9 * mim4U3 * mre4U2;
        const REAL FDPart3tmp531 = FDPart3tmp312 * mre4U1 + FDPart3tmp316 * mim4U1 + FDPart3tmp408 * mre4U3 + FDPart3tmp409 * mim4U3 -
                                   FDPart3tmp411 * mre4U3 - FDPart3tmp412 * mim4U3 - FDPart3tmp527 * n4U1 * n4U3;
        const REAL FDPart3tmp532 = FDPart3tmp315 * mre4U1 + FDPart3tmp318 * mim4U1 + FDPart3tmp408 * mre4U2 + FDPart3tmp409 * mim4U2 -
                                   FDPart3tmp411 * mre4U2 - FDPart3tmp412 * mim4U2 - FDPart3tmp527 * n4U1 * n4U2;
        const REAL FDPart3tmp17 = FDPart3tmp11 * FDPart3tmp13 + FDPart3tmp11 * FDPart3tmp15;
        const REAL FDPart3tmp27 = FDPart3tmp22 * FDPart3tmp26 * hDD01;
        const REAL FDPart3tmp37 = ((FDPart3tmp19) * (FDPart3tmp19) * (FDPart3tmp19) * (FDPart3tmp19)) * FDPart3tmp36;
        const REAL FDPart3tmp66 = FDPart3tmp19 * FDPart3tmp28;
        const REAL FDPart3tmp128 = 2 * FDPart3tmp11 * FDPart3tmp13 + 2 * FDPart3tmp11 * FDPart3tmp15;
        const REAL FDPart3tmp229 = FDPart3tmp28 * (2 * FDPart3tmp13 * FDPart3tmp226 - 2 * FDPart3tmp15 * FDPart3tmp226);
        const REAL FDPart3tmp18 = ((FDPart3tmp17) * (FDPart3tmp17));
        const REAL FDPart3tmp31 = ((FDPart3tmp19) * (FDPart3tmp19)) * FDPart3tmp28;
        const REAL FDPart3tmp39 = FDPart3tmp37 * ((hDD12) * (hDD12));
        const REAL FDPart3tmp65 = FDPart3tmp17 * ((FDPart3tmp19) * (FDPart3tmp19) * (FDPart3tmp19));
        const REAL FDPart3tmp67 = FDPart3tmp17 * FDPart3tmp66;
        const REAL FDPart3tmp129 = FDPart3tmp128 * FDPart3tmp66;
        const REAL FDPart3tmp230 = FDPart3tmp19 * FDPart3tmp229;
        const REAL FDPart3tmp231 = FDPart3tmp128 * FDPart3tmp17 * FDPart3tmp28;
        const REAL FDPart3tmp235 = FDPart3tmp66 * (FDPart3tmp13 * FDPart3tmp226 - FDPart3tmp15 * FDPart3tmp226);
        const REAL FDPart3tmp21 = FDPart3tmp18 * ((FDPart3tmp19) * (FDPart3tmp19));
        const REAL FDPart3tmp29 = FDPart3tmp18 * FDPart3tmp28;
        const REAL FDPart3tmp32 = FDPart3tmp22 * FDPart3tmp31;
        const REAL FDPart3tmp41 = FDPart3tmp31 * hDD11 + FDPart3tmp31;
        const REAL FDPart3tmp44 = FDPart3tmp31 * FDPart3tmp36;
        const REAL FDPart3tmp58 = FDPart3tmp31 * FDPart3tmp49;
        const REAL FDPart3tmp68 = FDPart3tmp67 * hDD01;
        const REAL FDPart3tmp70 = FDPart3tmp49 * FDPart3tmp67;
        const REAL FDPart3tmp81 = FDPart3tmp22 * FDPart3tmp67;
        const REAL FDPart3tmp130 = FDPart3tmp129 * FDPart3tmp36;
        const REAL FDPart3tmp142 = FDPart3tmp118 * FDPart3tmp67;
        const REAL FDPart3tmp200 = FDPart3tmp129 * hDD11 + FDPart3tmp129 + FDPart3tmp31 * hDD_dD110;
        const REAL FDPart3tmp232 = FDPart3tmp230 * FDPart3tmp36 + FDPart3tmp231 * FDPart3tmp36;
        const REAL FDPart3tmp30 = FDPart3tmp29 * hDD00 + FDPart3tmp29;
        const REAL FDPart3tmp33 = FDPart3tmp32 * hDD12;
        const REAL FDPart3tmp43 = FDPart3tmp21 * ((hDD02) * (hDD02));
        const REAL FDPart3tmp45 = FDPart3tmp44 * hDD22 + FDPart3tmp44;
        const REAL FDPart3tmp46 = FDPart3tmp21 * ((hDD01) * (hDD01));
        const REAL FDPart3tmp50 = FDPart3tmp44 * FDPart3tmp49 * aDD22;
        const REAL FDPart3tmp79 = FDPart3tmp23 * FDPart3tmp41;
        const REAL FDPart3tmp82 = FDPart3tmp81 * hDD02;
        const REAL FDPart3tmp95 = FDPart3tmp70 * hDD01;
        const REAL FDPart3tmp108 = FDPart3tmp41 * FDPart3tmp49;
        const REAL FDPart3tmp131 = FDPart3tmp130 * hDD22 + FDPart3tmp130 + FDPart3tmp44 * hDD_dD220;
        const REAL FDPart3tmp138 = FDPart3tmp129 * FDPart3tmp22 * hDD12 + FDPart3tmp32 * hDD_dD120;
        const REAL FDPart3tmp143 = FDPart3tmp142 * hDD02 + FDPart3tmp81 * hDD_dD021;
        const REAL FDPart3tmp189 = 2 * ((FDPart3tmp118) * (FDPart3tmp118)) * FDPart3tmp31 - 2 * FDPart3tmp44;
        const REAL FDPart3tmp192 = -FDPart3tmp145 * FDPart3tmp41 + FDPart3tmp58 * hDD_dD111;
        const REAL FDPart3tmp194 = FDPart3tmp118 * FDPart3tmp31 * hDD12 + FDPart3tmp32 * hDD_dD121;
        const REAL FDPart3tmp233 =
            FDPart3tmp184 * (2 * FDPart3tmp130 * hDD_dD220 + FDPart3tmp232 * hDD22 + FDPart3tmp232 + FDPart3tmp44 * hDD_dDD2200);
        const REAL FDPart3tmp236 = FDPart3tmp22 * FDPart3tmp235 + FDPart3tmp22 * FDPart3tmp29;
        const REAL FDPart3tmp291 = 2 * FDPart3tmp129 * FDPart3tmp22;
        const REAL FDPart3tmp379 = FDPart3tmp184 * (FDPart3tmp118 * FDPart3tmp129 * hDD12 + FDPart3tmp118 * FDPart3tmp31 * hDD_dD120 +
                                                    FDPart3tmp129 * FDPart3tmp22 * hDD_dD121 + FDPart3tmp32 * hDD_dDD1201);
        const REAL FDPart3tmp478 = FDPart3tmp184 * FDPart3tmp29 * hDD_dDD0011;
        const REAL FDPart3tmp479 = FDPart3tmp145 * FDPart3tmp29 * hDD_dD001;
        const REAL FDPart3tmp480 = FDPart3tmp140 * FDPart3tmp67 * hDD_dD011;
        const REAL FDPart3tmp488 = FDPart3tmp184 * (2 * FDPart3tmp129 * hDD_dD110 + FDPart3tmp230 + FDPart3tmp231 + FDPart3tmp31 * hDD_dDD1100 +
                                                    hDD11 * (FDPart3tmp230 + FDPart3tmp231));
        const REAL FDPart3tmp34 = FDPart3tmp21 * FDPart3tmp27 * hDD02 - FDPart3tmp23 * FDPart3tmp30 * FDPart3tmp33;
        const REAL FDPart3tmp51 = FDPart3tmp45 * FDPart3tmp49;
        const REAL FDPart3tmp60 = FDPart3tmp22 * FDPart3tmp58 * aDD12;
        const REAL FDPart3tmp61 = FDPart3tmp22 * FDPart3tmp58 * hDD12;
        const REAL FDPart3tmp69 = -FDPart3tmp23 * FDPart3tmp45 * FDPart3tmp68 + FDPart3tmp26 * FDPart3tmp36 * FDPart3tmp65 * hDD02 * hDD12;
        const REAL FDPart3tmp72 = FDPart3tmp22 * FDPart3tmp70 * aDD02;
        const REAL FDPart3tmp73 = FDPart3tmp22 * FDPart3tmp70 * hDD02;
        const REAL FDPart3tmp80 = -FDPart3tmp26 * FDPart3tmp39 + FDPart3tmp45 * FDPart3tmp79;
        const REAL FDPart3tmp83 = FDPart3tmp27 * FDPart3tmp65 * hDD12 - FDPart3tmp79 * FDPart3tmp82;
        const REAL FDPart3tmp89 = -FDPart3tmp26 * FDPart3tmp46 + FDPart3tmp30 * FDPart3tmp79;
        const REAL FDPart3tmp101 = FDPart3tmp29 * FDPart3tmp49 * aDD00;
        const REAL FDPart3tmp102 = FDPart3tmp30 * FDPart3tmp49;
        const REAL FDPart3tmp120 = 2 * FDPart3tmp118 * FDPart3tmp32;
        const REAL FDPart3tmp190 =
            FDPart3tmp184 * (4 * FDPart3tmp118 * FDPart3tmp32 * hDD_dD221 + FDPart3tmp189 * hDD22 + FDPart3tmp189 + FDPart3tmp44 * hDD_dDD2211);
        const REAL FDPart3tmp196 = -4 * FDPart3tmp115 * FDPart3tmp33 + FDPart3tmp194 * FDPart3tmp195;
        const REAL FDPart3tmp202 = -4 * FDPart3tmp115 * FDPart3tmp68 + FDPart3tmp140 * FDPart3tmp41 + FDPart3tmp195 * FDPart3tmp67 * hDD_dD011 -
                                   FDPart3tmp200 * FDPart3tmp49;
        const REAL FDPart3tmp237 = FDPart3tmp236 * hDD02 + FDPart3tmp81 * hDD_dD020;
        const REAL FDPart3tmp242 = FDPart3tmp29 * FDPart3tmp49 * hDD_dD001;
        const REAL FDPart3tmp245 = FDPart3tmp67 * hDD_dD010 + hDD01 * (FDPart3tmp235 + FDPart3tmp29);
        const REAL FDPart3tmp268 = -FDPart3tmp140 * FDPart3tmp41 + FDPart3tmp200 * FDPart3tmp49;
        const REAL FDPart3tmp345 = FDPart3tmp184 * (FDPart3tmp142 * hDD_dD020 + FDPart3tmp236 * hDD_dD021 + FDPart3tmp81 * hDD_dDD0201 +
                                                    hDD02 * (FDPart3tmp118 * FDPart3tmp235 + FDPart3tmp118 * FDPart3tmp29));
        const REAL FDPart3tmp346 = FDPart3tmp184 * (FDPart3tmp291 * hDD_dD120 + FDPart3tmp32 * hDD_dDD1200 +
                                                    hDD12 * (FDPart3tmp22 * FDPart3tmp230 + FDPart3tmp22 * FDPart3tmp231));
        const REAL FDPart3tmp378 = FDPart3tmp184 * (2 * FDPart3tmp142 * hDD_dD021 + FDPart3tmp81 * hDD_dDD0211 - FDPart3tmp82);
        const REAL FDPart3tmp485 = FDPart3tmp49 * (FDPart3tmp67 * hDD_dDD0101 + hDD_dD011 * (FDPart3tmp235 + FDPart3tmp29));
        const REAL FDPart3tmp48 =
            (1.0 /
             (2 * pow(AMPL, 6) * FDPart3tmp18 * FDPart3tmp35 * FDPart3tmp37 * hDD01 * hDD02 * hDD12 / pow(FDPart3tmp24, 6) -
              FDPart3tmp25 * FDPart3tmp30 * FDPart3tmp35 * FDPart3tmp39 - FDPart3tmp25 * FDPart3tmp35 * FDPart3tmp36 * FDPart3tmp41 * FDPart3tmp43 -
              FDPart3tmp25 * FDPart3tmp35 * FDPart3tmp45 * FDPart3tmp46 + FDPart3tmp30 * FDPart3tmp35 * FDPart3tmp41 * FDPart3tmp45));
        const REAL FDPart3tmp57 = FDPart3tmp23 * FDPart3tmp30 * FDPart3tmp45 - FDPart3tmp26 * FDPart3tmp36 * FDPart3tmp43;
        const REAL FDPart3tmp97 = FDPart3tmp52 * FDPart3tmp95 + FDPart3tmp70 * aDD01;
        const REAL FDPart3tmp110 = FDPart3tmp108 * FDPart3tmp52 + FDPart3tmp58 * aDD11;
        const REAL FDPart3tmp117 = 2 * FDPart3tmp115 * FDPart3tmp45;
        const REAL FDPart3tmp121 = FDPart3tmp120 * hDD22 + FDPart3tmp120 + FDPart3tmp44 * hDD_dD221;
        const REAL FDPart3tmp127 = 2 * FDPart3tmp126 * FDPart3tmp45;
        const REAL FDPart3tmp147 = FDPart3tmp143 * FDPart3tmp49 - FDPart3tmp145 * FDPart3tmp82;
        const REAL FDPart3tmp168 = FDPart3tmp138 * FDPart3tmp49 - FDPart3tmp140 * FDPart3tmp33;
        const REAL FDPart3tmp238 = -4 * FDPart3tmp126 * FDPart3tmp82 + FDPart3tmp195 * FDPart3tmp237;
        const REAL FDPart3tmp241 = -2 * FDPart3tmp126 * FDPart3tmp30 +
                                   FDPart3tmp49 * (FDPart3tmp17 * FDPart3tmp229 * hDD00 + FDPart3tmp17 * FDPart3tmp229 + FDPart3tmp29 * hDD_dD000);
        const REAL FDPart3tmp243 = 2 * FDPart3tmp115 * FDPart3tmp30;
        const REAL FDPart3tmp293 = FDPart3tmp184 * (FDPart3tmp118 * FDPart3tmp291 * hDD22 + FDPart3tmp118 * FDPart3tmp291 +
                                                    FDPart3tmp120 * hDD_dD220 + FDPart3tmp130 * hDD_dD221 + FDPart3tmp44 * hDD_dDD2201);
        const REAL FDPart3tmp439 = -4 * FDPart3tmp108 * FDPart3tmp52 - 4 * FDPart3tmp58 * aDD11;
        const REAL FDPart3tmp440 = -FDPart3tmp52 * FDPart3tmp95 - FDPart3tmp70 * aDD01;
        const REAL FDPart3tmp54 = FDPart3tmp50 + FDPart3tmp51 * FDPart3tmp52;
        const REAL FDPart3tmp63 = FDPart3tmp52 * FDPart3tmp61 + FDPart3tmp60;
        const REAL FDPart3tmp75 = FDPart3tmp52 * FDPart3tmp73 + FDPart3tmp72;
        const REAL FDPart3tmp78 = -4 * FDPart3tmp52 * FDPart3tmp61 - 4 * FDPart3tmp60;
        const REAL FDPart3tmp85 = -FDPart3tmp52 * FDPart3tmp73 - FDPart3tmp72;
        const REAL FDPart3tmp90 = FDPart3tmp48 * FDPart3tmp89;
        const REAL FDPart3tmp93 = -4 * FDPart3tmp50 - 4 * FDPart3tmp51 * FDPart3tmp52;
        const REAL FDPart3tmp98 = FDPart3tmp48 * FDPart3tmp97;
        const REAL FDPart3tmp104 = FDPart3tmp101 + FDPart3tmp102 * FDPart3tmp52;
        const REAL FDPart3tmp111 = FDPart3tmp110 * FDPart3tmp48;
        const REAL FDPart3tmp124 = (1.0 / 2.0) * FDPart3tmp48;
        const REAL FDPart3tmp133 = FDPart3tmp127 - FDPart3tmp131 * FDPart3tmp49;
        const REAL FDPart3tmp148 = -FDPart3tmp138 * FDPart3tmp49 + FDPart3tmp140 * FDPart3tmp33 + FDPart3tmp147;
        const REAL FDPart3tmp165 = 4 * FDPart3tmp48;
        const REAL FDPart3tmp167 = -FDPart3tmp127 + FDPart3tmp131 * FDPart3tmp49;
        const REAL FDPart3tmp169 = -FDPart3tmp143 * FDPart3tmp49 + FDPart3tmp145 * FDPart3tmp82 + FDPart3tmp168;
        const REAL FDPart3tmp246 = -4 * FDPart3tmp126 * FDPart3tmp68 + FDPart3tmp195 * FDPart3tmp245 - FDPart3tmp242 + FDPart3tmp243;
        const REAL FDPart3tmp267 = FDPart3tmp242 - FDPart3tmp243;
        const REAL FDPart3tmp269 = FDPart3tmp147 + FDPart3tmp168;
        const REAL FDPart3tmp441 = 4 * FDPart3tmp440;
        const REAL FDPart3tmp453 = -FDPart3tmp101 - FDPart3tmp102 * FDPart3tmp52;
        const REAL FDPart3tmp64 = FDPart3tmp48 * FDPart3tmp63;
        const REAL FDPart3tmp76 = FDPart3tmp48 * FDPart3tmp75;
        const REAL FDPart3tmp86 = 4 * FDPart3tmp85;
        const REAL FDPart3tmp99 = FDPart3tmp69 * FDPart3tmp98;
        const REAL FDPart3tmp105 = FDPart3tmp104 * FDPart3tmp48;
        const REAL FDPart3tmp123 = FDPart3tmp117 - FDPart3tmp121 * FDPart3tmp49;
        const REAL FDPart3tmp125 = FDPart3tmp124 * FDPart3tmp34;
        const REAL FDPart3tmp134 = FDPart3tmp124 * FDPart3tmp83;
        const REAL FDPart3tmp136 = -FDPart3tmp117 + FDPart3tmp121 * FDPart3tmp49;
        const REAL FDPart3tmp149 = FDPart3tmp124 * FDPart3tmp80;
        const REAL FDPart3tmp152 = FDPart3tmp124 * FDPart3tmp69;
        const REAL FDPart3tmp154 = (1.0 / 2.0) * FDPart3tmp90;
        const REAL FDPart3tmp159 = FDPart3tmp124 * FDPart3tmp57;
        const REAL FDPart3tmp166 = FDPart3tmp165 * FDPart3tmp34;
        const REAL FDPart3tmp180 = FDPart3tmp165 * FDPart3tmp83;
        const REAL FDPart3tmp221 = FDPart3tmp165 * FDPart3tmp57;
        const REAL FDPart3tmp265 = FDPart3tmp165 * FDPart3tmp80;
        const REAL FDPart3tmp308 = FDPart3tmp165 * FDPart3tmp69;
        const REAL FDPart3tmp323 = FDPart3tmp111 * FDPart3tmp34 + FDPart3tmp63 * FDPart3tmp90 + FDPart3tmp83 * FDPart3tmp98;
        const REAL FDPart3tmp326 = 4 * FDPart3tmp90;
        const REAL FDPart3tmp77 = FDPart3tmp34 * FDPart3tmp48 * FDPart3tmp54 + FDPart3tmp57 * FDPart3tmp64 + FDPart3tmp69 * FDPart3tmp76;
        const REAL FDPart3tmp84 = FDPart3tmp48 * FDPart3tmp54 * FDPart3tmp83 + FDPart3tmp64 * FDPart3tmp69 + FDPart3tmp76 * FDPart3tmp80;
        const REAL FDPart3tmp87 = FDPart3tmp34 * FDPart3tmp64;
        const REAL FDPart3tmp88 = FDPart3tmp76 * FDPart3tmp83;
        const REAL FDPart3tmp135 = FDPart3tmp123 * FDPart3tmp125 + FDPart3tmp133 * FDPart3tmp134;
        const REAL FDPart3tmp150 = FDPart3tmp134 * FDPart3tmp136 + FDPart3tmp148 * FDPart3tmp149;
        const REAL FDPart3tmp153 = FDPart3tmp123 * FDPart3tmp152 + FDPart3tmp133 * FDPart3tmp149;
        const REAL FDPart3tmp155 = FDPart3tmp134 * FDPart3tmp148 + FDPart3tmp136 * FDPart3tmp154;
        const REAL FDPart3tmp157 = FDPart3tmp125 * FDPart3tmp136 + FDPart3tmp148 * FDPart3tmp152;
        const REAL FDPart3tmp160 = FDPart3tmp123 * FDPart3tmp159 + FDPart3tmp133 * FDPart3tmp152;
        const REAL FDPart3tmp170 = FDPart3tmp134 * FDPart3tmp167 + FDPart3tmp152 * FDPart3tmp169;
        const REAL FDPart3tmp172 = FDPart3tmp125 * FDPart3tmp169 + FDPart3tmp154 * FDPart3tmp167;
        const REAL FDPart3tmp174 = FDPart3tmp125 * FDPart3tmp167 + FDPart3tmp159 * FDPart3tmp169;
        const REAL FDPart3tmp203 = FDPart3tmp134 * FDPart3tmp196 + FDPart3tmp149 * FDPart3tmp202 + FDPart3tmp152 * FDPart3tmp192;
        const REAL FDPart3tmp206 = FDPart3tmp125 * FDPart3tmp192 + FDPart3tmp134 * FDPart3tmp202 + FDPart3tmp154 * FDPart3tmp196;
        const REAL FDPart3tmp209 = FDPart3tmp125 * FDPart3tmp196 + FDPart3tmp152 * FDPart3tmp202 + FDPart3tmp159 * FDPart3tmp192;
        const REAL FDPart3tmp247 = FDPart3tmp134 * FDPart3tmp238 + FDPart3tmp149 * FDPart3tmp241 + FDPart3tmp152 * FDPart3tmp246;
        const REAL FDPart3tmp250 = FDPart3tmp125 * FDPart3tmp246 + FDPart3tmp134 * FDPart3tmp241 + FDPart3tmp154 * FDPart3tmp238;
        const REAL FDPart3tmp253 = FDPart3tmp125 * FDPart3tmp238 + FDPart3tmp152 * FDPart3tmp241 + FDPart3tmp159 * FDPart3tmp246;
        const REAL FDPart3tmp270 = FDPart3tmp125 * FDPart3tmp269 + FDPart3tmp152 * FDPart3tmp267 + FDPart3tmp159 * FDPart3tmp268;
        const REAL FDPart3tmp274 = FDPart3tmp134 * FDPart3tmp269 + FDPart3tmp149 * FDPart3tmp267 + FDPart3tmp152 * FDPart3tmp268;
        const REAL FDPart3tmp278 = FDPart3tmp125 * FDPart3tmp268 + FDPart3tmp134 * FDPart3tmp267 + FDPart3tmp154 * FDPart3tmp269;
        const REAL FDPart3tmp325 = FDPart3tmp111 * FDPart3tmp69 + FDPart3tmp64 * FDPart3tmp83 + FDPart3tmp80 * FDPart3tmp98;
        const REAL FDPart3tmp414 = FDPart3tmp105 * FDPart3tmp83 + FDPart3tmp34 * FDPart3tmp98 + FDPart3tmp76 * FDPart3tmp89;
        const REAL FDPart3tmp415 = FDPart3tmp105 * FDPart3tmp69 + FDPart3tmp34 * FDPart3tmp76 + FDPart3tmp57 * FDPart3tmp98;
        const REAL FDPart3tmp92 = FDPart3tmp54 * FDPart3tmp90 + FDPart3tmp87 + FDPart3tmp88;
        const REAL FDPart3tmp113 = 4 * FDPart3tmp105 * FDPart3tmp80 + 4 * FDPart3tmp111 * FDPart3tmp57 + 4 * FDPart3tmp54 * FDPart3tmp90 +
                                   8 * FDPart3tmp87 + 8 * FDPart3tmp88 + 8 * FDPart3tmp99;
        const REAL FDPart3tmp211 = FDPart3tmp150 * FDPart3tmp155;
        const REAL FDPart3tmp213 = FDPart3tmp155 * FDPart3tmp157;
        const REAL FDPart3tmp215 = FDPart3tmp150 * FDPart3tmp157;
        const REAL FDPart3tmp255 = FDPart3tmp170 * FDPart3tmp172;
        const REAL FDPart3tmp258 = FDPart3tmp172 * FDPart3tmp174;
        const REAL FDPart3tmp260 = FDPart3tmp170 * FDPart3tmp174;
        const REAL FDPart3tmp266 = FDPart3tmp150 * FDPart3tmp174;
        const REAL FDPart3tmp271 = FDPart3tmp153 * FDPart3tmp270;
        const REAL FDPart3tmp273 = FDPart3tmp157 * FDPart3tmp170;
        const REAL FDPart3tmp275 = FDPart3tmp160 * FDPart3tmp274;
        const REAL FDPart3tmp277 = FDPart3tmp157 * FDPart3tmp172;
        const REAL FDPart3tmp279 = FDPart3tmp160 * FDPart3tmp278;
        const REAL FDPart3tmp281 = FDPart3tmp155 * FDPart3tmp174;
        const REAL FDPart3tmp282 = FDPart3tmp135 * FDPart3tmp270;
        const REAL FDPart3tmp284 = FDPart3tmp150 * FDPart3tmp172;
        const REAL FDPart3tmp285 = FDPart3tmp153 * FDPart3tmp278;
        const REAL FDPart3tmp287 = FDPart3tmp155 * FDPart3tmp170;
        const REAL FDPart3tmp288 = FDPart3tmp135 * FDPart3tmp274;
        const REAL FDPart3tmp324 = FDPart3tmp111 * FDPart3tmp57 + FDPart3tmp87 + FDPart3tmp99;
        const REAL FDPart3tmp327 = FDPart3tmp174 * FDPart3tmp274;
        const REAL FDPart3tmp328 = FDPart3tmp157 * FDPart3tmp247;
        const REAL FDPart3tmp330 = FDPart3tmp170 * FDPart3tmp270;
        const REAL FDPart3tmp331 = FDPart3tmp150 * FDPart3tmp253;
        const REAL FDPart3tmp333 = FDPart3tmp172 * FDPart3tmp270;
        const REAL FDPart3tmp334 = FDPart3tmp155 * FDPart3tmp253;
        const REAL FDPart3tmp336 = FDPart3tmp174 * FDPart3tmp278;
        const REAL FDPart3tmp337 = FDPart3tmp157 * FDPart3tmp250;
        const REAL FDPart3tmp339 = FDPart3tmp172 * FDPart3tmp274;
        const REAL FDPart3tmp340 = FDPart3tmp155 * FDPart3tmp247;
        const REAL FDPart3tmp342 = FDPart3tmp170 * FDPart3tmp278;
        const REAL FDPart3tmp343 = FDPart3tmp150 * FDPart3tmp250;
        const REAL FDPart3tmp360 = FDPart3tmp174 * FDPart3tmp203;
        const REAL FDPart3tmp361 = FDPart3tmp157 * FDPart3tmp274;
        const REAL FDPart3tmp363 = FDPart3tmp170 * FDPart3tmp209;
        const REAL FDPart3tmp364 = FDPart3tmp150 * FDPart3tmp270;
        const REAL FDPart3tmp366 = FDPart3tmp172 * FDPart3tmp209;
        const REAL FDPart3tmp367 = FDPart3tmp155 * FDPart3tmp270;
        const REAL FDPart3tmp369 = FDPart3tmp174 * FDPart3tmp206;
        const REAL FDPart3tmp370 = FDPart3tmp157 * FDPart3tmp278;
        const REAL FDPart3tmp372 = FDPart3tmp172 * FDPart3tmp203;
        const REAL FDPart3tmp373 = FDPart3tmp155 * FDPart3tmp274;
        const REAL FDPart3tmp375 = FDPart3tmp170 * FDPart3tmp206;
        const REAL FDPart3tmp376 = FDPart3tmp150 * FDPart3tmp278;
        const REAL FDPart3tmp416 = 4 * FDPart3tmp105 * FDPart3tmp80 + 4 * FDPart3tmp88 + 4 * FDPart3tmp99;
        const REAL FDPart3tmp495 = FDPart3tmp274 * FDPart3tmp278;
        const REAL FDPart3tmp498 = FDPart3tmp270 * FDPart3tmp274;
        const REAL FDPart3tmp501 = FDPart3tmp270 * FDPart3tmp278;
        const REAL FDPart3tmp164 = FDPart3tmp61 * (-FDPart3tmp135 * FDPart3tmp157 + FDPart3tmp155 * FDPart3tmp160) +
                                   FDPart3tmp61 * (FDPart3tmp135 * FDPart3tmp157 - FDPart3tmp155 * FDPart3tmp160) +
                                   FDPart3tmp73 * (-FDPart3tmp135 * FDPart3tmp150 + FDPart3tmp153 * FDPart3tmp155) +
                                   FDPart3tmp73 * (FDPart3tmp135 * FDPart3tmp150 - FDPart3tmp153 * FDPart3tmp155) +
                                   FDPart3tmp95 * (-FDPart3tmp150 * FDPart3tmp160 + FDPart3tmp153 * FDPart3tmp157) +
                                   FDPart3tmp95 * (FDPart3tmp150 * FDPart3tmp160 - FDPart3tmp153 * FDPart3tmp157);
        const REAL FDPart3tmp179 = FDPart3tmp61 * (-FDPart3tmp135 * FDPart3tmp174 + FDPart3tmp160 * FDPart3tmp172) +
                                   FDPart3tmp61 * (FDPart3tmp135 * FDPart3tmp174 - FDPart3tmp160 * FDPart3tmp172) +
                                   FDPart3tmp73 * (-FDPart3tmp135 * FDPart3tmp170 + FDPart3tmp153 * FDPart3tmp172) +
                                   FDPart3tmp73 * (FDPart3tmp135 * FDPart3tmp170 - FDPart3tmp153 * FDPart3tmp172) +
                                   FDPart3tmp95 * (-FDPart3tmp153 * FDPart3tmp174 + FDPart3tmp160 * FDPart3tmp170) +
                                   FDPart3tmp95 * (FDPart3tmp153 * FDPart3tmp174 - FDPart3tmp160 * FDPart3tmp170);
        const REAL FDPart3tmp220 =
            FDPart3tmp102 * (((FDPart3tmp150) * (FDPart3tmp150)) - FDPart3tmp153 * FDPart3tmp203) +
            FDPart3tmp108 * (((FDPart3tmp157) * (FDPart3tmp157)) - FDPart3tmp160 * FDPart3tmp209) + FDPart3tmp121 * FDPart3tmp145 -
            FDPart3tmp186 * FDPart3tmp45 - FDPart3tmp190 + FDPart3tmp51 * (-FDPart3tmp135 * FDPart3tmp206 + ((FDPart3tmp155) * (FDPart3tmp155))) +
            FDPart3tmp61 * (-FDPart3tmp135 * FDPart3tmp209 + FDPart3tmp213) + FDPart3tmp61 * (-FDPart3tmp160 * FDPart3tmp206 + FDPart3tmp213) +
            FDPart3tmp73 * (-FDPart3tmp135 * FDPart3tmp203 + FDPart3tmp211) + FDPart3tmp73 * (-FDPart3tmp153 * FDPart3tmp206 + FDPart3tmp211) +
            FDPart3tmp95 * (-FDPart3tmp153 * FDPart3tmp209 + FDPart3tmp215) + FDPart3tmp95 * (-FDPart3tmp160 * FDPart3tmp203 + FDPart3tmp215);
        const REAL FDPart3tmp264 =
            FDPart3tmp102 * (-FDPart3tmp153 * FDPart3tmp247 + ((FDPart3tmp170) * (FDPart3tmp170))) +
            FDPart3tmp108 * (-FDPart3tmp160 * FDPart3tmp253 + ((FDPart3tmp174) * (FDPart3tmp174))) + FDPart3tmp131 * FDPart3tmp140 -
            FDPart3tmp223 * FDPart3tmp45 - FDPart3tmp233 + FDPart3tmp51 * (-FDPart3tmp135 * FDPart3tmp250 + ((FDPart3tmp172) * (FDPart3tmp172))) +
            FDPart3tmp61 * (-FDPart3tmp135 * FDPart3tmp253 + FDPart3tmp258) + FDPart3tmp61 * (-FDPart3tmp160 * FDPart3tmp250 + FDPart3tmp258) +
            FDPart3tmp73 * (-FDPart3tmp135 * FDPart3tmp247 + FDPart3tmp255) + FDPart3tmp73 * (-FDPart3tmp153 * FDPart3tmp250 + FDPart3tmp255) +
            FDPart3tmp95 * (-FDPart3tmp153 * FDPart3tmp253 + FDPart3tmp260) + FDPart3tmp95 * (-FDPart3tmp160 * FDPart3tmp247 + FDPart3tmp260);
        const REAL FDPart3tmp306 = FDPart3tmp102 * (FDPart3tmp150 * FDPart3tmp170 - FDPart3tmp153 * FDPart3tmp274) +
                                   FDPart3tmp108 * (FDPart3tmp157 * FDPart3tmp174 - FDPart3tmp160 * FDPart3tmp270) + FDPart3tmp115 * FDPart3tmp131 +
                                   FDPart3tmp121 * FDPart3tmp126 - FDPart3tmp293 - FDPart3tmp296 * FDPart3tmp45 +
                                   FDPart3tmp51 * (-FDPart3tmp135 * FDPart3tmp278 + FDPart3tmp155 * FDPart3tmp172);
        const REAL FDPart3tmp358 = FDPart3tmp102 * (-FDPart3tmp150 * FDPart3tmp247 + FDPart3tmp170 * FDPart3tmp274) +
                                   FDPart3tmp108 * (-FDPart3tmp157 * FDPart3tmp253 + FDPart3tmp174 * FDPart3tmp270) - FDPart3tmp115 * FDPart3tmp237 -
                                   FDPart3tmp126 * FDPart3tmp143 + FDPart3tmp138 * FDPart3tmp140 - FDPart3tmp223 * FDPart3tmp33 +
                                   FDPart3tmp296 * FDPart3tmp82 + FDPart3tmp345 - FDPart3tmp346 +
                                   FDPart3tmp51 * (-FDPart3tmp155 * FDPart3tmp250 + FDPart3tmp172 * FDPart3tmp278);
        const REAL FDPart3tmp391 = FDPart3tmp102 * (-FDPart3tmp150 * FDPart3tmp274 + FDPart3tmp170 * FDPart3tmp203) +
                                   FDPart3tmp108 * (-FDPart3tmp157 * FDPart3tmp270 + FDPart3tmp174 * FDPart3tmp209) + FDPart3tmp115 * FDPart3tmp138 +
                                   FDPart3tmp126 * FDPart3tmp194 - FDPart3tmp143 * FDPart3tmp145 + FDPart3tmp186 * FDPart3tmp82 -
                                   FDPart3tmp296 * FDPart3tmp33 + FDPart3tmp378 - FDPart3tmp379 +
                                   FDPart3tmp51 * (-FDPart3tmp155 * FDPart3tmp278 + FDPart3tmp172 * FDPart3tmp206);
        const REAL FDPart3tmp400 = FDPart3tmp102 * (-FDPart3tmp150 * FDPart3tmp170 + FDPart3tmp153 * FDPart3tmp274) +
                                   FDPart3tmp108 * (-FDPart3tmp157 * FDPart3tmp174 + FDPart3tmp160 * FDPart3tmp270) - FDPart3tmp115 * FDPart3tmp131 -
                                   FDPart3tmp121 * FDPart3tmp126 + FDPart3tmp293 + FDPart3tmp296 * FDPart3tmp45 +
                                   FDPart3tmp51 * (FDPart3tmp135 * FDPart3tmp278 - FDPart3tmp155 * FDPart3tmp172);
        const REAL FDPart3tmp423 = FDPart3tmp102 * (FDPart3tmp150 * FDPart3tmp274 - FDPart3tmp170 * FDPart3tmp203) +
                                   FDPart3tmp108 * (FDPart3tmp157 * FDPart3tmp270 - FDPart3tmp174 * FDPart3tmp209) - FDPart3tmp115 * FDPart3tmp138 -
                                   FDPart3tmp126 * FDPart3tmp194 + FDPart3tmp143 * FDPart3tmp145 - FDPart3tmp186 * FDPart3tmp82 +
                                   FDPart3tmp296 * FDPart3tmp33 - FDPart3tmp378 + FDPart3tmp379 +
                                   FDPart3tmp51 * (FDPart3tmp155 * FDPart3tmp278 - FDPart3tmp172 * FDPart3tmp206);
        const REAL FDPart3tmp431 = FDPart3tmp102 * (FDPart3tmp150 * FDPart3tmp247 - FDPart3tmp170 * FDPart3tmp274) +
                                   FDPart3tmp108 * (FDPart3tmp157 * FDPart3tmp253 - FDPart3tmp174 * FDPart3tmp270) + FDPart3tmp115 * FDPart3tmp237 +
                                   FDPart3tmp126 * FDPart3tmp143 - FDPart3tmp138 * FDPart3tmp140 + FDPart3tmp223 * FDPart3tmp33 -
                                   FDPart3tmp296 * FDPart3tmp82 - FDPart3tmp345 + FDPart3tmp346 +
                                   FDPart3tmp51 * (FDPart3tmp155 * FDPart3tmp250 - FDPart3tmp172 * FDPart3tmp278);
        const REAL FDPart3tmp449 = FDPart3tmp61 * (-FDPart3tmp155 * FDPart3tmp209 + FDPart3tmp157 * FDPart3tmp206) +
                                   FDPart3tmp61 * (FDPart3tmp155 * FDPart3tmp209 - FDPart3tmp157 * FDPart3tmp206) +
                                   FDPart3tmp73 * (-FDPart3tmp150 * FDPart3tmp206 + FDPart3tmp155 * FDPart3tmp203) +
                                   FDPart3tmp73 * (FDPart3tmp150 * FDPart3tmp206 - FDPart3tmp155 * FDPart3tmp203) +
                                   FDPart3tmp95 * (-FDPart3tmp150 * FDPart3tmp209 + FDPart3tmp157 * FDPart3tmp203) +
                                   FDPart3tmp95 * (FDPart3tmp150 * FDPart3tmp209 - FDPart3tmp157 * FDPart3tmp203);
        const REAL FDPart3tmp461 = FDPart3tmp61 * (-FDPart3tmp172 * FDPart3tmp253 + FDPart3tmp174 * FDPart3tmp250) +
                                   FDPart3tmp61 * (FDPart3tmp172 * FDPart3tmp253 - FDPart3tmp174 * FDPart3tmp250) +
                                   FDPart3tmp73 * (-FDPart3tmp170 * FDPart3tmp250 + FDPart3tmp172 * FDPart3tmp247) +
                                   FDPart3tmp73 * (FDPart3tmp170 * FDPart3tmp250 - FDPart3tmp172 * FDPart3tmp247) +
                                   FDPart3tmp95 * (-FDPart3tmp170 * FDPart3tmp253 + FDPart3tmp174 * FDPart3tmp247) +
                                   FDPart3tmp95 * (FDPart3tmp170 * FDPart3tmp253 - FDPart3tmp174 * FDPart3tmp247);
        const REAL FDPart3tmp477 = FDPart3tmp61 * (-FDPart3tmp206 * FDPart3tmp270 + FDPart3tmp209 * FDPart3tmp278) +
                                   FDPart3tmp61 * (FDPart3tmp206 * FDPart3tmp270 - FDPart3tmp209 * FDPart3tmp278) +
                                   FDPart3tmp73 * (-FDPart3tmp203 * FDPart3tmp278 + FDPart3tmp206 * FDPart3tmp274) +
                                   FDPart3tmp73 * (FDPart3tmp203 * FDPart3tmp278 - FDPart3tmp206 * FDPart3tmp274) +
                                   FDPart3tmp95 * (-FDPart3tmp203 * FDPart3tmp270 + FDPart3tmp209 * FDPart3tmp274) +
                                   FDPart3tmp95 * (FDPart3tmp203 * FDPart3tmp270 - FDPart3tmp209 * FDPart3tmp274);
        const REAL FDPart3tmp504 =
            FDPart3tmp165 *
            (FDPart3tmp102 * (-FDPart3tmp203 * FDPart3tmp247 + ((FDPart3tmp274) * (FDPart3tmp274))) +
             FDPart3tmp108 * (-FDPart3tmp209 * FDPart3tmp253 + ((FDPart3tmp270) * (FDPart3tmp270))) + FDPart3tmp140 * FDPart3tmp200 -
             FDPart3tmp145 * FDPart3tmp245 - FDPart3tmp186 * FDPart3tmp30 - FDPart3tmp223 * FDPart3tmp41 + FDPart3tmp295 * FDPart3tmp68 -
             FDPart3tmp478 + FDPart3tmp479 - FDPart3tmp480 + FDPart3tmp485 - FDPart3tmp488 +
             FDPart3tmp51 * (-FDPart3tmp206 * FDPart3tmp250 + ((FDPart3tmp278) * (FDPart3tmp278))) +
             FDPart3tmp61 * (-FDPart3tmp206 * FDPart3tmp253 + FDPart3tmp501) + FDPart3tmp61 * (-FDPart3tmp209 * FDPart3tmp250 + FDPart3tmp501) +
             FDPart3tmp73 * (-FDPart3tmp203 * FDPart3tmp250 + FDPart3tmp495) + FDPart3tmp73 * (-FDPart3tmp206 * FDPart3tmp247 + FDPart3tmp495) +
             FDPart3tmp95 * (-FDPart3tmp203 * FDPart3tmp253 + FDPart3tmp498) + FDPart3tmp95 * (-FDPart3tmp209 * FDPart3tmp247 + FDPart3tmp498));
        const REAL FDPart3tmp524 = FDPart3tmp61 * (-FDPart3tmp250 * FDPart3tmp270 + FDPart3tmp253 * FDPart3tmp278) +
                                   FDPart3tmp61 * (FDPart3tmp250 * FDPart3tmp270 - FDPart3tmp253 * FDPart3tmp278) +
                                   FDPart3tmp73 * (-FDPart3tmp247 * FDPart3tmp278 + FDPart3tmp250 * FDPart3tmp274) +
                                   FDPart3tmp73 * (FDPart3tmp247 * FDPart3tmp278 - FDPart3tmp250 * FDPart3tmp274) +
                                   FDPart3tmp95 * (-FDPart3tmp247 * FDPart3tmp270 + FDPart3tmp253 * FDPart3tmp274) +
                                   FDPart3tmp95 * (FDPart3tmp247 * FDPart3tmp270 - FDPart3tmp253 * FDPart3tmp274);
        const REAL FDPart3tmp307 = FDPart3tmp306 + FDPart3tmp61 * (FDPart3tmp277 - FDPart3tmp279) + FDPart3tmp61 * (FDPart3tmp281 - FDPart3tmp282) +
                                   FDPart3tmp73 * (FDPart3tmp284 - FDPart3tmp285) + FDPart3tmp73 * (FDPart3tmp287 - FDPart3tmp288) +
                                   FDPart3tmp95 * (FDPart3tmp266 - FDPart3tmp271) + FDPart3tmp95 * (FDPart3tmp273 - FDPart3tmp275);
        const REAL FDPart3tmp309 = FDPart3tmp306 + FDPart3tmp61 * (FDPart3tmp277 - FDPart3tmp282) + FDPart3tmp61 * (-FDPart3tmp279 + FDPart3tmp281) +
                                   FDPart3tmp73 * (FDPart3tmp284 - FDPart3tmp288) + FDPart3tmp73 * (-FDPart3tmp285 + FDPart3tmp287) +
                                   FDPart3tmp95 * (FDPart3tmp266 - FDPart3tmp275) + FDPart3tmp95 * (-FDPart3tmp271 + FDPart3tmp273);
        const REAL FDPart3tmp359 = FDPart3tmp358 + FDPart3tmp61 * (FDPart3tmp333 - FDPart3tmp334) + FDPart3tmp61 * (FDPart3tmp336 - FDPart3tmp337) +
                                   FDPart3tmp73 * (FDPart3tmp339 - FDPart3tmp340) + FDPart3tmp73 * (FDPart3tmp342 - FDPart3tmp343) +
                                   FDPart3tmp95 * (FDPart3tmp327 - FDPart3tmp328) + FDPart3tmp95 * (FDPart3tmp330 - FDPart3tmp331);
        const REAL FDPart3tmp392 = FDPart3tmp391 + FDPart3tmp61 * (FDPart3tmp366 - FDPart3tmp367) + FDPart3tmp61 * (FDPart3tmp369 - FDPart3tmp370) +
                                   FDPart3tmp73 * (FDPart3tmp372 - FDPart3tmp373) + FDPart3tmp73 * (FDPart3tmp375 - FDPart3tmp376) +
                                   FDPart3tmp95 * (FDPart3tmp360 - FDPart3tmp361) + FDPart3tmp95 * (FDPart3tmp363 - FDPart3tmp364);
        const REAL FDPart3tmp399 = FDPart3tmp61 * (-FDPart3tmp277 + FDPart3tmp281) + FDPart3tmp61 * (FDPart3tmp277 - FDPart3tmp281) +
                                   FDPart3tmp73 * (-FDPart3tmp284 + FDPart3tmp287) + FDPart3tmp73 * (FDPart3tmp284 - FDPart3tmp287) +
                                   FDPart3tmp95 * (-FDPart3tmp266 + FDPart3tmp273) + FDPart3tmp95 * (FDPart3tmp266 - FDPart3tmp273);
        const REAL FDPart3tmp401 = FDPart3tmp400 + FDPart3tmp61 * (-FDPart3tmp277 + FDPart3tmp279) + FDPart3tmp61 * (-FDPart3tmp281 + FDPart3tmp282) +
                                   FDPart3tmp73 * (-FDPart3tmp284 + FDPart3tmp285) + FDPart3tmp73 * (-FDPart3tmp287 + FDPart3tmp288) +
                                   FDPart3tmp95 * (-FDPart3tmp266 + FDPart3tmp271) + FDPart3tmp95 * (-FDPart3tmp273 + FDPart3tmp275);
        const REAL FDPart3tmp405 =
            FDPart3tmp113 * FDPart3tmp63 +
            FDPart3tmp166 *
                (FDPart3tmp102 * (-((FDPart3tmp150) * (FDPart3tmp150)) + FDPart3tmp153 * FDPart3tmp203) +
                 FDPart3tmp108 * (-((FDPart3tmp157) * (FDPart3tmp157)) + FDPart3tmp160 * FDPart3tmp209) - FDPart3tmp121 * FDPart3tmp145 +
                 FDPart3tmp186 * FDPart3tmp45 + FDPart3tmp190 + FDPart3tmp51 * (FDPart3tmp135 * FDPart3tmp206 - ((FDPart3tmp155) * (FDPart3tmp155))) +
                 FDPart3tmp61 * (FDPart3tmp135 * FDPart3tmp209 - FDPart3tmp213) + FDPart3tmp61 * (FDPart3tmp160 * FDPart3tmp206 - FDPart3tmp213) +
                 FDPart3tmp73 * (FDPart3tmp135 * FDPart3tmp203 - FDPart3tmp211) + FDPart3tmp73 * (FDPart3tmp153 * FDPart3tmp206 - FDPart3tmp211) +
                 FDPart3tmp95 * (FDPart3tmp153 * FDPart3tmp209 - FDPart3tmp215) + FDPart3tmp95 * (FDPart3tmp160 * FDPart3tmp203 - FDPart3tmp215));
        const REAL FDPart3tmp424 = FDPart3tmp423 + FDPart3tmp61 * (-FDPart3tmp366 + FDPart3tmp367) + FDPart3tmp61 * (-FDPart3tmp369 + FDPart3tmp370) +
                                   FDPart3tmp73 * (-FDPart3tmp372 + FDPart3tmp373) + FDPart3tmp73 * (-FDPart3tmp375 + FDPart3tmp376) +
                                   FDPart3tmp95 * (-FDPart3tmp360 + FDPart3tmp361) + FDPart3tmp95 * (-FDPart3tmp363 + FDPart3tmp364);
        const REAL FDPart3tmp432 = FDPart3tmp431 + FDPart3tmp61 * (-FDPart3tmp333 + FDPart3tmp334) + FDPart3tmp61 * (-FDPart3tmp336 + FDPart3tmp337) +
                                   FDPart3tmp73 * (-FDPart3tmp339 + FDPart3tmp340) + FDPart3tmp73 * (-FDPart3tmp342 + FDPart3tmp343) +
                                   FDPart3tmp95 * (-FDPart3tmp327 + FDPart3tmp328) + FDPart3tmp95 * (-FDPart3tmp330 + FDPart3tmp331);
        const REAL FDPart3tmp433 = FDPart3tmp400 + FDPart3tmp61 * (-FDPart3tmp277 + FDPart3tmp282) + FDPart3tmp61 * (FDPart3tmp279 - FDPart3tmp281) +
                                   FDPart3tmp73 * (-FDPart3tmp284 + FDPart3tmp288) + FDPart3tmp73 * (FDPart3tmp285 - FDPart3tmp287) +
                                   FDPart3tmp95 * (-FDPart3tmp266 + FDPart3tmp275) + FDPart3tmp95 * (FDPart3tmp271 - FDPart3tmp273);
        const REAL FDPart3tmp437 =
            FDPart3tmp113 * FDPart3tmp75 +
            FDPart3tmp180 *
                (FDPart3tmp102 * (FDPart3tmp153 * FDPart3tmp247 - ((FDPart3tmp170) * (FDPart3tmp170))) +
                 FDPart3tmp108 * (FDPart3tmp160 * FDPart3tmp253 - ((FDPart3tmp174) * (FDPart3tmp174))) - FDPart3tmp131 * FDPart3tmp140 +
                 FDPart3tmp223 * FDPart3tmp45 + FDPart3tmp233 + FDPart3tmp51 * (FDPart3tmp135 * FDPart3tmp250 - ((FDPart3tmp172) * (FDPart3tmp172))) +
                 FDPart3tmp61 * (FDPart3tmp135 * FDPart3tmp253 - FDPart3tmp258) + FDPart3tmp61 * (FDPart3tmp160 * FDPart3tmp250 - FDPart3tmp258) +
                 FDPart3tmp73 * (FDPart3tmp135 * FDPart3tmp247 - FDPart3tmp255) + FDPart3tmp73 * (FDPart3tmp153 * FDPart3tmp250 - FDPart3tmp255) +
                 FDPart3tmp95 * (FDPart3tmp153 * FDPart3tmp253 - FDPart3tmp260) + FDPart3tmp95 * (FDPart3tmp160 * FDPart3tmp247 - FDPart3tmp260));
        const REAL FDPart3tmp442 = FDPart3tmp358 + FDPart3tmp61 * (FDPart3tmp333 - FDPart3tmp337) + FDPart3tmp61 * (-FDPart3tmp334 + FDPart3tmp336) +
                                   FDPart3tmp73 * (FDPart3tmp339 - FDPart3tmp343) + FDPart3tmp73 * (-FDPart3tmp340 + FDPart3tmp342) +
                                   FDPart3tmp95 * (FDPart3tmp327 - FDPart3tmp331) + FDPart3tmp95 * (-FDPart3tmp328 + FDPart3tmp330);
        const REAL FDPart3tmp450 = FDPart3tmp61 * (-FDPart3tmp367 + FDPart3tmp370) + FDPart3tmp61 * (FDPart3tmp367 - FDPart3tmp370) +
                                   FDPart3tmp73 * (-FDPart3tmp373 + FDPart3tmp376) + FDPart3tmp73 * (FDPart3tmp373 - FDPart3tmp376) +
                                   FDPart3tmp95 * (-FDPart3tmp361 + FDPart3tmp364) + FDPart3tmp95 * (FDPart3tmp361 - FDPart3tmp364);
        const REAL FDPart3tmp451 = FDPart3tmp391 + FDPart3tmp61 * (FDPart3tmp366 - FDPart3tmp370) + FDPart3tmp61 * (-FDPart3tmp367 + FDPart3tmp369) +
                                   FDPart3tmp73 * (FDPart3tmp372 - FDPart3tmp376) + FDPart3tmp73 * (-FDPart3tmp373 + FDPart3tmp375) +
                                   FDPart3tmp95 * (FDPart3tmp360 - FDPart3tmp364) + FDPart3tmp95 * (-FDPart3tmp361 + FDPart3tmp363);
        const REAL FDPart3tmp462 = FDPart3tmp423 + FDPart3tmp61 * (-FDPart3tmp366 + FDPart3tmp370) + FDPart3tmp61 * (FDPart3tmp367 - FDPart3tmp369) +
                                   FDPart3tmp73 * (-FDPart3tmp372 + FDPart3tmp376) + FDPart3tmp73 * (FDPart3tmp373 - FDPart3tmp375) +
                                   FDPart3tmp95 * (-FDPart3tmp360 + FDPart3tmp364) + FDPart3tmp95 * (FDPart3tmp361 - FDPart3tmp363);
        const REAL FDPart3tmp463 = FDPart3tmp61 * (-FDPart3tmp333 + FDPart3tmp336) + FDPart3tmp61 * (FDPart3tmp333 - FDPart3tmp336) +
                                   FDPart3tmp73 * (-FDPart3tmp339 + FDPart3tmp342) + FDPart3tmp73 * (FDPart3tmp339 - FDPart3tmp342) +
                                   FDPart3tmp95 * (-FDPart3tmp327 + FDPart3tmp330) + FDPart3tmp95 * (FDPart3tmp327 - FDPart3tmp330);
        const REAL FDPart3tmp464 = FDPart3tmp431 + FDPart3tmp61 * (-FDPart3tmp333 + FDPart3tmp337) + FDPart3tmp61 * (FDPart3tmp334 - FDPart3tmp336) +
                                   FDPart3tmp73 * (-FDPart3tmp339 + FDPart3tmp343) + FDPart3tmp73 * (FDPart3tmp340 - FDPart3tmp342) +
                                   FDPart3tmp95 * (-FDPart3tmp327 + FDPart3tmp331) + FDPart3tmp95 * (FDPart3tmp328 - FDPart3tmp330);
        const REAL FDPart3tmp511 =
            FDPart3tmp113 * FDPart3tmp97 +
            FDPart3tmp308 *
                (FDPart3tmp102 * (FDPart3tmp203 * FDPart3tmp247 - ((FDPart3tmp274) * (FDPart3tmp274))) +
                 FDPart3tmp108 * (FDPart3tmp209 * FDPart3tmp253 - ((FDPart3tmp270) * (FDPart3tmp270))) - FDPart3tmp140 * FDPart3tmp200 +
                 FDPart3tmp145 * FDPart3tmp245 + FDPart3tmp186 * FDPart3tmp30 + FDPart3tmp223 * FDPart3tmp41 - FDPart3tmp295 * FDPart3tmp68 +
                 FDPart3tmp478 - FDPart3tmp479 + FDPart3tmp480 - FDPart3tmp485 + FDPart3tmp488 +
                 FDPart3tmp51 * (FDPart3tmp206 * FDPart3tmp250 - ((FDPart3tmp278) * (FDPart3tmp278))) +
                 FDPart3tmp61 * (FDPart3tmp206 * FDPart3tmp253 - FDPart3tmp501) + FDPart3tmp61 * (FDPart3tmp209 * FDPart3tmp250 - FDPart3tmp501) +
                 FDPart3tmp73 * (FDPart3tmp203 * FDPart3tmp250 - FDPart3tmp495) + FDPart3tmp73 * (FDPart3tmp206 * FDPart3tmp247 - FDPart3tmp495) +
                 FDPart3tmp95 * (FDPart3tmp203 * FDPart3tmp253 - FDPart3tmp498) + FDPart3tmp95 * (FDPart3tmp209 * FDPart3tmp247 - FDPart3tmp498));
        const REAL FDPart3tmp310 = FDPart3tmp113 * FDPart3tmp54 + FDPart3tmp164 * FDPart3tmp166 + FDPart3tmp179 * FDPart3tmp180 +
                                   FDPart3tmp220 * FDPart3tmp221 + FDPart3tmp264 * FDPart3tmp265 + FDPart3tmp307 * FDPart3tmp308 +
                                   FDPart3tmp308 * FDPart3tmp309 + FDPart3tmp77 * FDPart3tmp78 + FDPart3tmp84 * FDPart3tmp86 +
                                   FDPart3tmp92 * FDPart3tmp93;
        const REAL FDPart3tmp406 = FDPart3tmp164 * FDPart3tmp326 + FDPart3tmp180 * FDPart3tmp399 + FDPart3tmp180 * FDPart3tmp401 +
                                   FDPart3tmp265 * FDPart3tmp359 + FDPart3tmp308 * FDPart3tmp392 + FDPart3tmp323 * FDPart3tmp93 +
                                   FDPart3tmp324 * FDPart3tmp78 + FDPart3tmp325 * FDPart3tmp86 + FDPart3tmp405;
        const REAL FDPart3tmp438 = FDPart3tmp166 * FDPart3tmp399 + FDPart3tmp166 * FDPart3tmp433 + FDPart3tmp179 * FDPart3tmp326 +
                                   FDPart3tmp221 * FDPart3tmp424 + FDPart3tmp308 * FDPart3tmp432 + FDPart3tmp414 * FDPart3tmp93 +
                                   FDPart3tmp415 * FDPart3tmp78 + FDPart3tmp416 * FDPart3tmp85 + FDPart3tmp437;
        const REAL FDPart3tmp452 = FDPart3tmp180 * FDPart3tmp433 + FDPart3tmp221 * FDPart3tmp449 + FDPart3tmp265 * FDPart3tmp442 +
                                   FDPart3tmp308 * FDPart3tmp450 + FDPart3tmp308 * FDPart3tmp451 + FDPart3tmp405 + FDPart3tmp439 * FDPart3tmp77 +
                                   FDPart3tmp441 * FDPart3tmp84 + FDPart3tmp78 * FDPart3tmp92;
        const REAL FDPart3tmp465 = FDPart3tmp166 * FDPart3tmp401 + FDPart3tmp221 * FDPart3tmp462 + FDPart3tmp265 * FDPart3tmp461 +
                                   FDPart3tmp308 * FDPart3tmp463 + FDPart3tmp308 * FDPart3tmp464 + FDPart3tmp437 + FDPart3tmp441 * FDPart3tmp77 +
                                   4 * FDPart3tmp453 * FDPart3tmp84 + FDPart3tmp86 * FDPart3tmp92;
        const REAL FDPart3tmp505 = FDPart3tmp110 * FDPart3tmp113 + FDPart3tmp166 * FDPart3tmp449 + FDPart3tmp180 * FDPart3tmp424 +
                                   FDPart3tmp180 * FDPart3tmp462 + FDPart3tmp220 * FDPart3tmp326 + FDPart3tmp308 * FDPart3tmp477 +
                                   FDPart3tmp323 * FDPart3tmp78 + FDPart3tmp324 * FDPart3tmp439 + FDPart3tmp325 * FDPart3tmp441 +
                                   FDPart3tmp504 * FDPart3tmp80;
        const REAL FDPart3tmp512 = FDPart3tmp166 * FDPart3tmp392 + FDPart3tmp166 * FDPart3tmp450 + FDPart3tmp180 * FDPart3tmp464 +
                                   FDPart3tmp221 * FDPart3tmp477 + FDPart3tmp309 * FDPart3tmp326 + FDPart3tmp414 * FDPart3tmp78 +
                                   FDPart3tmp415 * FDPart3tmp439 + FDPart3tmp416 * FDPart3tmp440 + FDPart3tmp511;
        const REAL FDPart3tmp525 = FDPart3tmp104 * FDPart3tmp113 + FDPart3tmp166 * FDPart3tmp359 + FDPart3tmp166 * FDPart3tmp442 +
                                   FDPart3tmp180 * FDPart3tmp461 + FDPart3tmp264 * FDPart3tmp326 + FDPart3tmp308 * FDPart3tmp524 +
                                   FDPart3tmp414 * FDPart3tmp86 + FDPart3tmp415 * FDPart3tmp441 + FDPart3tmp416 * FDPart3tmp453 +
                                   FDPart3tmp504 * FDPart3tmp57;
        const REAL FDPart3tmp526 = FDPart3tmp166 * FDPart3tmp451 + FDPart3tmp180 * FDPart3tmp432 + FDPart3tmp180 * FDPart3tmp463 +
                                   FDPart3tmp265 * FDPart3tmp524 + FDPart3tmp307 * FDPart3tmp326 + FDPart3tmp323 * FDPart3tmp86 +
                                   FDPart3tmp324 * FDPart3tmp441 + 4 * FDPart3tmp325 * FDPart3tmp453 + FDPart3tmp511;
        diagnostic_output_gfs[IDX4(PSI4_PART2REGF, i0, i1, i2)] =
            FDPart3tmp310 * (FDPart3tmp0 * FDPart3tmp1 * mim4U3 - FDPart3tmp1 * FDPart3tmp3 * mre4U3 + FDPart3tmp10 * FDPart3tmp7 -
                             FDPart3tmp5 * FDPart3tmp7 - FDPart3tmp9 * ((mim4U3) * (mim4U3)) + FDPart3tmp9 * ((mre4U3) * (mre4U3))) +
            FDPart3tmp322 * FDPart3tmp406 + FDPart3tmp322 * FDPart3tmp452 + FDPart3tmp413 * FDPart3tmp438 + FDPart3tmp413 * FDPart3tmp465 +
            FDPart3tmp505 * (FDPart3tmp10 * FDPart3tmp470 + FDPart3tmp467 * mim4U2 - FDPart3tmp468 * mre4U2 - FDPart3tmp470 * FDPart3tmp5 -
                             FDPart3tmp9 * ((mim4U2) * (mim4U2)) + FDPart3tmp9 * ((mre4U2) * (mre4U2))) +
            FDPart3tmp507 * FDPart3tmp512 + FDPart3tmp507 * FDPart3tmp526 +
            FDPart3tmp525 * (FDPart3tmp10 * FDPart3tmp517 - FDPart3tmp5 * FDPart3tmp517 + FDPart3tmp514 * mim4U1 - FDPart3tmp515 * mre4U1 -
                             FDPart3tmp9 * ((mim4U1) * (mim4U1)) + FDPart3tmp9 * ((mre4U1) * (mre4U1)));
        diagnostic_output_gfs[IDX4(PSI4_PART2IMGF, i0, i1, i2)] =
            FDPart3tmp310 * (FDPart3tmp0 * FDPart3tmp1 * mre4U3 + FDPart3tmp1 * FDPart3tmp3 * mim4U3 - FDPart3tmp527 * ((n4U3) * (n4U3)) -
                             FDPart3tmp528 * mim4U3 * mre4U3) +
            FDPart3tmp406 * FDPart3tmp529 + FDPart3tmp438 * FDPart3tmp531 + FDPart3tmp452 * FDPart3tmp529 + FDPart3tmp465 * FDPart3tmp531 +
            FDPart3tmp505 * (FDPart3tmp467 * mre4U2 + FDPart3tmp468 * mim4U2 - FDPart3tmp527 * ((n4U2) * (n4U2)) - FDPart3tmp528 * mim4U2 * mre4U2) +
            FDPart3tmp512 * FDPart3tmp532 +
            FDPart3tmp525 * (FDPart3tmp514 * mre4U1 + FDPart3tmp515 * mim4U1 - FDPart3tmp527 * ((n4U1) * (n4U1)) - FDPart3tmp528 * mim4U1 * mre4U1) +
            FDPart3tmp526 * FDPart3tmp532;

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
