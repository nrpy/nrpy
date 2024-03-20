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
 * Compute psi4 at all interior gridpoints, part 1
 */
void psi4_part1__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
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
        const REAL aDD00_i1m4 = in_gfs[IDX4(ADD00GF, i0, i1 - 4, i2)];
        const REAL aDD00_i1m3 = in_gfs[IDX4(ADD00GF, i0, i1 - 3, i2)];
        const REAL aDD00_i1m2 = in_gfs[IDX4(ADD00GF, i0, i1 - 2, i2)];
        const REAL aDD00_i1m1 = in_gfs[IDX4(ADD00GF, i0, i1 - 1, i2)];
        const REAL aDD00 = in_gfs[IDX4(ADD00GF, i0, i1, i2)];
        const REAL aDD00_i1p1 = in_gfs[IDX4(ADD00GF, i0, i1 + 1, i2)];
        const REAL aDD00_i1p2 = in_gfs[IDX4(ADD00GF, i0, i1 + 2, i2)];
        const REAL aDD00_i1p3 = in_gfs[IDX4(ADD00GF, i0, i1 + 3, i2)];
        const REAL aDD00_i1p4 = in_gfs[IDX4(ADD00GF, i0, i1 + 4, i2)];
        const REAL aDD01_i1m4 = in_gfs[IDX4(ADD01GF, i0, i1 - 4, i2)];
        const REAL aDD01_i1m3 = in_gfs[IDX4(ADD01GF, i0, i1 - 3, i2)];
        const REAL aDD01_i1m2 = in_gfs[IDX4(ADD01GF, i0, i1 - 2, i2)];
        const REAL aDD01_i1m1 = in_gfs[IDX4(ADD01GF, i0, i1 - 1, i2)];
        const REAL aDD01_i0m4 = in_gfs[IDX4(ADD01GF, i0 - 4, i1, i2)];
        const REAL aDD01_i0m3 = in_gfs[IDX4(ADD01GF, i0 - 3, i1, i2)];
        const REAL aDD01_i0m2 = in_gfs[IDX4(ADD01GF, i0 - 2, i1, i2)];
        const REAL aDD01_i0m1 = in_gfs[IDX4(ADD01GF, i0 - 1, i1, i2)];
        const REAL aDD01 = in_gfs[IDX4(ADD01GF, i0, i1, i2)];
        const REAL aDD01_i0p1 = in_gfs[IDX4(ADD01GF, i0 + 1, i1, i2)];
        const REAL aDD01_i0p2 = in_gfs[IDX4(ADD01GF, i0 + 2, i1, i2)];
        const REAL aDD01_i0p3 = in_gfs[IDX4(ADD01GF, i0 + 3, i1, i2)];
        const REAL aDD01_i0p4 = in_gfs[IDX4(ADD01GF, i0 + 4, i1, i2)];
        const REAL aDD01_i1p1 = in_gfs[IDX4(ADD01GF, i0, i1 + 1, i2)];
        const REAL aDD01_i1p2 = in_gfs[IDX4(ADD01GF, i0, i1 + 2, i2)];
        const REAL aDD01_i1p3 = in_gfs[IDX4(ADD01GF, i0, i1 + 3, i2)];
        const REAL aDD01_i1p4 = in_gfs[IDX4(ADD01GF, i0, i1 + 4, i2)];
        const REAL aDD02_i1m4 = in_gfs[IDX4(ADD02GF, i0, i1 - 4, i2)];
        const REAL aDD02_i1m3 = in_gfs[IDX4(ADD02GF, i0, i1 - 3, i2)];
        const REAL aDD02_i1m2 = in_gfs[IDX4(ADD02GF, i0, i1 - 2, i2)];
        const REAL aDD02_i1m1 = in_gfs[IDX4(ADD02GF, i0, i1 - 1, i2)];
        const REAL aDD02_i0m4 = in_gfs[IDX4(ADD02GF, i0 - 4, i1, i2)];
        const REAL aDD02_i0m3 = in_gfs[IDX4(ADD02GF, i0 - 3, i1, i2)];
        const REAL aDD02_i0m2 = in_gfs[IDX4(ADD02GF, i0 - 2, i1, i2)];
        const REAL aDD02_i0m1 = in_gfs[IDX4(ADD02GF, i0 - 1, i1, i2)];
        const REAL aDD02 = in_gfs[IDX4(ADD02GF, i0, i1, i2)];
        const REAL aDD02_i0p1 = in_gfs[IDX4(ADD02GF, i0 + 1, i1, i2)];
        const REAL aDD02_i0p2 = in_gfs[IDX4(ADD02GF, i0 + 2, i1, i2)];
        const REAL aDD02_i0p3 = in_gfs[IDX4(ADD02GF, i0 + 3, i1, i2)];
        const REAL aDD02_i0p4 = in_gfs[IDX4(ADD02GF, i0 + 4, i1, i2)];
        const REAL aDD02_i1p1 = in_gfs[IDX4(ADD02GF, i0, i1 + 1, i2)];
        const REAL aDD02_i1p2 = in_gfs[IDX4(ADD02GF, i0, i1 + 2, i2)];
        const REAL aDD02_i1p3 = in_gfs[IDX4(ADD02GF, i0, i1 + 3, i2)];
        const REAL aDD02_i1p4 = in_gfs[IDX4(ADD02GF, i0, i1 + 4, i2)];
        const REAL aDD11_i0m4 = in_gfs[IDX4(ADD11GF, i0 - 4, i1, i2)];
        const REAL aDD11_i0m3 = in_gfs[IDX4(ADD11GF, i0 - 3, i1, i2)];
        const REAL aDD11_i0m2 = in_gfs[IDX4(ADD11GF, i0 - 2, i1, i2)];
        const REAL aDD11_i0m1 = in_gfs[IDX4(ADD11GF, i0 - 1, i1, i2)];
        const REAL aDD11 = in_gfs[IDX4(ADD11GF, i0, i1, i2)];
        const REAL aDD11_i0p1 = in_gfs[IDX4(ADD11GF, i0 + 1, i1, i2)];
        const REAL aDD11_i0p2 = in_gfs[IDX4(ADD11GF, i0 + 2, i1, i2)];
        const REAL aDD11_i0p3 = in_gfs[IDX4(ADD11GF, i0 + 3, i1, i2)];
        const REAL aDD11_i0p4 = in_gfs[IDX4(ADD11GF, i0 + 4, i1, i2)];
        const REAL aDD12_i1m4 = in_gfs[IDX4(ADD12GF, i0, i1 - 4, i2)];
        const REAL aDD12_i1m3 = in_gfs[IDX4(ADD12GF, i0, i1 - 3, i2)];
        const REAL aDD12_i1m2 = in_gfs[IDX4(ADD12GF, i0, i1 - 2, i2)];
        const REAL aDD12_i1m1 = in_gfs[IDX4(ADD12GF, i0, i1 - 1, i2)];
        const REAL aDD12_i0m4 = in_gfs[IDX4(ADD12GF, i0 - 4, i1, i2)];
        const REAL aDD12_i0m3 = in_gfs[IDX4(ADD12GF, i0 - 3, i1, i2)];
        const REAL aDD12_i0m2 = in_gfs[IDX4(ADD12GF, i0 - 2, i1, i2)];
        const REAL aDD12_i0m1 = in_gfs[IDX4(ADD12GF, i0 - 1, i1, i2)];
        const REAL aDD12 = in_gfs[IDX4(ADD12GF, i0, i1, i2)];
        const REAL aDD12_i0p1 = in_gfs[IDX4(ADD12GF, i0 + 1, i1, i2)];
        const REAL aDD12_i0p2 = in_gfs[IDX4(ADD12GF, i0 + 2, i1, i2)];
        const REAL aDD12_i0p3 = in_gfs[IDX4(ADD12GF, i0 + 3, i1, i2)];
        const REAL aDD12_i0p4 = in_gfs[IDX4(ADD12GF, i0 + 4, i1, i2)];
        const REAL aDD12_i1p1 = in_gfs[IDX4(ADD12GF, i0, i1 + 1, i2)];
        const REAL aDD12_i1p2 = in_gfs[IDX4(ADD12GF, i0, i1 + 2, i2)];
        const REAL aDD12_i1p3 = in_gfs[IDX4(ADD12GF, i0, i1 + 3, i2)];
        const REAL aDD12_i1p4 = in_gfs[IDX4(ADD12GF, i0, i1 + 4, i2)];
        const REAL aDD22_i1m4 = in_gfs[IDX4(ADD22GF, i0, i1 - 4, i2)];
        const REAL aDD22_i1m3 = in_gfs[IDX4(ADD22GF, i0, i1 - 3, i2)];
        const REAL aDD22_i1m2 = in_gfs[IDX4(ADD22GF, i0, i1 - 2, i2)];
        const REAL aDD22_i1m1 = in_gfs[IDX4(ADD22GF, i0, i1 - 1, i2)];
        const REAL aDD22_i0m4 = in_gfs[IDX4(ADD22GF, i0 - 4, i1, i2)];
        const REAL aDD22_i0m3 = in_gfs[IDX4(ADD22GF, i0 - 3, i1, i2)];
        const REAL aDD22_i0m2 = in_gfs[IDX4(ADD22GF, i0 - 2, i1, i2)];
        const REAL aDD22_i0m1 = in_gfs[IDX4(ADD22GF, i0 - 1, i1, i2)];
        const REAL aDD22 = in_gfs[IDX4(ADD22GF, i0, i1, i2)];
        const REAL aDD22_i0p1 = in_gfs[IDX4(ADD22GF, i0 + 1, i1, i2)];
        const REAL aDD22_i0p2 = in_gfs[IDX4(ADD22GF, i0 + 2, i1, i2)];
        const REAL aDD22_i0p3 = in_gfs[IDX4(ADD22GF, i0 + 3, i1, i2)];
        const REAL aDD22_i0p4 = in_gfs[IDX4(ADD22GF, i0 + 4, i1, i2)];
        const REAL aDD22_i1p1 = in_gfs[IDX4(ADD22GF, i0, i1 + 1, i2)];
        const REAL aDD22_i1p2 = in_gfs[IDX4(ADD22GF, i0, i1 + 2, i2)];
        const REAL aDD22_i1p3 = in_gfs[IDX4(ADD22GF, i0, i1 + 3, i2)];
        const REAL aDD22_i1p4 = in_gfs[IDX4(ADD22GF, i0, i1 + 4, i2)];
        const REAL cf_i1m4 = in_gfs[IDX4(CFGF, i0, i1 - 4, i2)];
        const REAL cf_i1m3 = in_gfs[IDX4(CFGF, i0, i1 - 3, i2)];
        const REAL cf_i1m2 = in_gfs[IDX4(CFGF, i0, i1 - 2, i2)];
        const REAL cf_i1m1 = in_gfs[IDX4(CFGF, i0, i1 - 1, i2)];
        const REAL cf_i0m4 = in_gfs[IDX4(CFGF, i0 - 4, i1, i2)];
        const REAL cf_i0m3 = in_gfs[IDX4(CFGF, i0 - 3, i1, i2)];
        const REAL cf_i0m2 = in_gfs[IDX4(CFGF, i0 - 2, i1, i2)];
        const REAL cf_i0m1 = in_gfs[IDX4(CFGF, i0 - 1, i1, i2)];
        const REAL cf = in_gfs[IDX4(CFGF, i0, i1, i2)];
        const REAL cf_i0p1 = in_gfs[IDX4(CFGF, i0 + 1, i1, i2)];
        const REAL cf_i0p2 = in_gfs[IDX4(CFGF, i0 + 2, i1, i2)];
        const REAL cf_i0p3 = in_gfs[IDX4(CFGF, i0 + 3, i1, i2)];
        const REAL cf_i0p4 = in_gfs[IDX4(CFGF, i0 + 4, i1, i2)];
        const REAL cf_i1p1 = in_gfs[IDX4(CFGF, i0, i1 + 1, i2)];
        const REAL cf_i1p2 = in_gfs[IDX4(CFGF, i0, i1 + 2, i2)];
        const REAL cf_i1p3 = in_gfs[IDX4(CFGF, i0, i1 + 3, i2)];
        const REAL cf_i1p4 = in_gfs[IDX4(CFGF, i0, i1 + 4, i2)];
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
        const REAL trK_i1m4 = in_gfs[IDX4(TRKGF, i0, i1 - 4, i2)];
        const REAL trK_i1m3 = in_gfs[IDX4(TRKGF, i0, i1 - 3, i2)];
        const REAL trK_i1m2 = in_gfs[IDX4(TRKGF, i0, i1 - 2, i2)];
        const REAL trK_i1m1 = in_gfs[IDX4(TRKGF, i0, i1 - 1, i2)];
        const REAL trK_i0m4 = in_gfs[IDX4(TRKGF, i0 - 4, i1, i2)];
        const REAL trK_i0m3 = in_gfs[IDX4(TRKGF, i0 - 3, i1, i2)];
        const REAL trK_i0m2 = in_gfs[IDX4(TRKGF, i0 - 2, i1, i2)];
        const REAL trK_i0m1 = in_gfs[IDX4(TRKGF, i0 - 1, i1, i2)];
        const REAL trK = in_gfs[IDX4(TRKGF, i0, i1, i2)];
        const REAL trK_i0p1 = in_gfs[IDX4(TRKGF, i0 + 1, i1, i2)];
        const REAL trK_i0p2 = in_gfs[IDX4(TRKGF, i0 + 2, i1, i2)];
        const REAL trK_i0p3 = in_gfs[IDX4(TRKGF, i0 + 3, i1, i2)];
        const REAL trK_i0p4 = in_gfs[IDX4(TRKGF, i0 + 4, i1, i2)];
        const REAL trK_i1p1 = in_gfs[IDX4(TRKGF, i0, i1 + 1, i2)];
        const REAL trK_i1p2 = in_gfs[IDX4(TRKGF, i0, i1 + 2, i2)];
        const REAL trK_i1p3 = in_gfs[IDX4(TRKGF, i0, i1 + 3, i2)];
        const REAL trK_i1p4 = in_gfs[IDX4(TRKGF, i0, i1 + 4, i2)];
        const REAL aDD_dD001 =
            fd_function_dD1_fdorder8(aDD00_i1m1, aDD00_i1m2, aDD00_i1m3, aDD00_i1m4, aDD00_i1p1, aDD00_i1p2, aDD00_i1p3, aDD00_i1p4, invdxx1);
        const REAL aDD_dD010 =
            fd_function_dD0_fdorder8(aDD01_i0m1, aDD01_i0m2, aDD01_i0m3, aDD01_i0m4, aDD01_i0p1, aDD01_i0p2, aDD01_i0p3, aDD01_i0p4, invdxx0);
        const REAL aDD_dD011 =
            fd_function_dD1_fdorder8(aDD01_i1m1, aDD01_i1m2, aDD01_i1m3, aDD01_i1m4, aDD01_i1p1, aDD01_i1p2, aDD01_i1p3, aDD01_i1p4, invdxx1);
        const REAL aDD_dD020 =
            fd_function_dD0_fdorder8(aDD02_i0m1, aDD02_i0m2, aDD02_i0m3, aDD02_i0m4, aDD02_i0p1, aDD02_i0p2, aDD02_i0p3, aDD02_i0p4, invdxx0);
        const REAL aDD_dD021 =
            fd_function_dD1_fdorder8(aDD02_i1m1, aDD02_i1m2, aDD02_i1m3, aDD02_i1m4, aDD02_i1p1, aDD02_i1p2, aDD02_i1p3, aDD02_i1p4, invdxx1);
        const REAL aDD_dD110 =
            fd_function_dD0_fdorder8(aDD11_i0m1, aDD11_i0m2, aDD11_i0m3, aDD11_i0m4, aDD11_i0p1, aDD11_i0p2, aDD11_i0p3, aDD11_i0p4, invdxx0);
        const REAL aDD_dD120 =
            fd_function_dD0_fdorder8(aDD12_i0m1, aDD12_i0m2, aDD12_i0m3, aDD12_i0m4, aDD12_i0p1, aDD12_i0p2, aDD12_i0p3, aDD12_i0p4, invdxx0);
        const REAL aDD_dD121 =
            fd_function_dD1_fdorder8(aDD12_i1m1, aDD12_i1m2, aDD12_i1m3, aDD12_i1m4, aDD12_i1p1, aDD12_i1p2, aDD12_i1p3, aDD12_i1p4, invdxx1);
        const REAL aDD_dD220 =
            fd_function_dD0_fdorder8(aDD22_i0m1, aDD22_i0m2, aDD22_i0m3, aDD22_i0m4, aDD22_i0p1, aDD22_i0p2, aDD22_i0p3, aDD22_i0p4, invdxx0);
        const REAL aDD_dD221 =
            fd_function_dD1_fdorder8(aDD22_i1m1, aDD22_i1m2, aDD22_i1m3, aDD22_i1m4, aDD22_i1p1, aDD22_i1p2, aDD22_i1p3, aDD22_i1p4, invdxx1);
        const REAL cf_dD0 = fd_function_dD0_fdorder8(cf_i0m1, cf_i0m2, cf_i0m3, cf_i0m4, cf_i0p1, cf_i0p2, cf_i0p3, cf_i0p4, invdxx0);
        const REAL cf_dD1 = fd_function_dD1_fdorder8(cf_i1m1, cf_i1m2, cf_i1m3, cf_i1m4, cf_i1p1, cf_i1p2, cf_i1p3, cf_i1p4, invdxx1);
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
        const REAL trK_dD0 = fd_function_dD0_fdorder8(trK_i0m1, trK_i0m2, trK_i0m3, trK_i0m4, trK_i0p1, trK_i0p2, trK_i0p3, trK_i0p4, invdxx0);
        const REAL trK_dD1 = fd_function_dD1_fdorder8(trK_i1m1, trK_i1m2, trK_i1m3, trK_i1m4, trK_i1p1, trK_i1p2, trK_i1p3, trK_i1p4, invdxx1);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        const REAL FDPart3tmp0 = n4U0 * n4U3;
        const REAL FDPart3tmp2 = ((n4U3) * (n4U3));
        const REAL FDPart3tmp3 = mim4U0 * mim4U2;
        const REAL FDPart3tmp5 = mre4U0 * mre4U2;
        const REAL FDPart3tmp6 = sin(xx1);
        const REAL FDPart3tmp8 = (1.0 / (SINHW));
        const REAL FDPart3tmp19 = (1.0 / ((cf) * (cf) * (cf)));
        const REAL FDPart3tmp23 = (1.0 / ((cf) * (cf)));
        const REAL FDPart3tmp26 = (2.0 / 3.0) * trK_dD1;
        const REAL FDPart3tmp28 = cos(xx1);
        const REAL FDPart3tmp38 = (2.0 / 3.0) * trK;
        const REAL FDPart3tmp45 = (1.0 / ((cf) * (cf) * (cf) * (cf)));
        const REAL FDPart3tmp54 = pow(cf, -6);
        const REAL FDPart3tmp80 = (1.0 / 3.0) * trK;
        const REAL FDPart3tmp121 = mim4U0 * mim4U3;
        const REAL FDPart3tmp122 = n4U2 * n4U3;
        const REAL FDPart3tmp124 = n4U0 * n4U2;
        const REAL FDPart3tmp125 = mre4U0 * mre4U3;
        const REAL FDPart3tmp128 = mim4U0 * mim4U1;
        const REAL FDPart3tmp129 = mre4U0 * mre4U1;
        const REAL FDPart3tmp132 = (2.0 / 3.0) * trK_dD0;
        const REAL FDPart3tmp153 = n4U1 * n4U3;
        const REAL FDPart3tmp154 = n4U0 * n4U1;
        const REAL FDPart3tmp186 = ((n4U2) * (n4U2));
        const REAL FDPart3tmp212 = mre4U1 * mre4U2;
        const REAL FDPart3tmp213 = mim4U1 * mim4U2;
        const REAL FDPart3tmp230 = n4U1 * n4U2;
        const REAL FDPart3tmp236 = ((n4U1) * (n4U1));
        const REAL FDPart3tmp239 = (1.0 / ((SINHW) * (SINHW)));
        const REAL FDPart3tmp299 = mim4U0 * mre4U3;
        const REAL FDPart3tmp300 = mim4U3 * mre4U0;
        const REAL FDPart3tmp303 = mim4U0 * mre4U2;
        const REAL FDPart3tmp304 = mim4U2 * mre4U0;
        const REAL FDPart3tmp306 = mim4U0 * mre4U1;
        const REAL FDPart3tmp307 = mim4U1 * mre4U0;
        const REAL FDPart3tmp318 = mim4U1 * mre4U2;
        const REAL FDPart3tmp319 = mim4U2 * mre4U1;
        const REAL FDPart3tmp1 = FDPart3tmp0 * mre4U3;
        const REAL FDPart3tmp4 = FDPart3tmp0 * mim4U3;
        const REAL FDPart3tmp7 = ((FDPart3tmp6) * (FDPart3tmp6));
        const REAL FDPart3tmp14 = exp(FDPart3tmp8) - exp(-FDPart3tmp8);
        const REAL FDPart3tmp31 = 2 * FDPart3tmp23;
        const REAL FDPart3tmp33 = 2 * FDPart3tmp19;
        const REAL FDPart3tmp156 = FDPart3tmp122 * FDPart3tmp128 - FDPart3tmp122 * FDPart3tmp129;
        const REAL FDPart3tmp159 = -FDPart3tmp124 * mim4U1 * mim4U3 + FDPart3tmp124 * mre4U1 * mre4U3;
        const REAL FDPart3tmp179 = FDPart3tmp153 * FDPart3tmp3 - FDPart3tmp153 * FDPart3tmp5;
        const REAL FDPart3tmp182 = -FDPart3tmp154 * mim4U2 * mim4U3 + FDPart3tmp154 * mre4U2 * mre4U3;
        const REAL FDPart3tmp214 = FDPart3tmp0 * FDPart3tmp212 - FDPart3tmp0 * FDPart3tmp213;
        const REAL FDPart3tmp231 = FDPart3tmp121 * FDPart3tmp230 - FDPart3tmp125 * FDPart3tmp230;
        const REAL FDPart3tmp308 = FDPart3tmp122 * FDPart3tmp306 + FDPart3tmp122 * FDPart3tmp307;
        const REAL FDPart3tmp311 = -FDPart3tmp124 * mim4U1 * mre4U3 - FDPart3tmp124 * mim4U3 * mre4U1;
        const REAL FDPart3tmp312 = FDPart3tmp153 * FDPart3tmp303 + FDPart3tmp153 * FDPart3tmp304;
        const REAL FDPart3tmp315 = -FDPart3tmp154 * mim4U2 * mre4U3 - FDPart3tmp154 * mim4U3 * mre4U2;
        const REAL FDPart3tmp320 = -FDPart3tmp0 * FDPart3tmp318 - FDPart3tmp0 * FDPart3tmp319;
        const REAL FDPart3tmp321 = FDPart3tmp230 * FDPart3tmp299 + FDPart3tmp230 * FDPart3tmp300;
        const REAL FDPart3tmp10 = exp(FDPart3tmp8 * xx0);
        const REAL FDPart3tmp11 = exp(-FDPart3tmp8 * xx0);
        const REAL FDPart3tmp15 = ((AMPL) * (AMPL)) / ((FDPart3tmp14) * (FDPart3tmp14));
        const REAL FDPart3tmp21 = 4 * FDPart3tmp19 * cf_dD1;
        const REAL FDPart3tmp34 = FDPart3tmp33 * cf_dD1;
        const REAL FDPart3tmp46 = ((AMPL) * (AMPL) * (AMPL) * (AMPL)) / ((FDPart3tmp14) * (FDPart3tmp14) * (FDPart3tmp14) * (FDPart3tmp14));
        const REAL FDPart3tmp72 = FDPart3tmp33 * cf_dD0;
        const REAL FDPart3tmp130 = 4 * FDPart3tmp19 * cf_dD0;
        const REAL FDPart3tmp12 = FDPart3tmp10 - FDPart3tmp11;
        const REAL FDPart3tmp47 = FDPart3tmp45 * FDPart3tmp46;
        const REAL FDPart3tmp42 = FDPart3tmp10 * FDPart3tmp8 + FDPart3tmp11 * FDPart3tmp8;
        const REAL FDPart3tmp55 = ((FDPart3tmp12) * (FDPart3tmp12) * (FDPart3tmp12) * (FDPart3tmp12)) * FDPart3tmp7;
        const REAL FDPart3tmp67 = FDPart3tmp12 * FDPart3tmp15;
        const REAL FDPart3tmp16 = ((FDPart3tmp12) * (FDPart3tmp12)) * FDPart3tmp15;
        const REAL FDPart3tmp43 = ((FDPart3tmp42) * (FDPart3tmp42));
        const REAL FDPart3tmp57 = FDPart3tmp55 * ((hDD12) * (hDD12));
        const REAL FDPart3tmp66 = ((FDPart3tmp12) * (FDPart3tmp12) * (FDPart3tmp12)) * FDPart3tmp42;
        const REAL FDPart3tmp68 = FDPart3tmp42 * FDPart3tmp67;
        const REAL FDPart3tmp74 = FDPart3tmp67 * (2 * FDPart3tmp10 * FDPart3tmp8 + 2 * FDPart3tmp11 * FDPart3tmp8);
        const REAL FDPart3tmp242 = FDPart3tmp67 * (FDPart3tmp10 * FDPart3tmp239 - FDPart3tmp11 * FDPart3tmp239);
        const REAL FDPart3tmp251 = FDPart3tmp15 * FDPart3tmp42 * (2 * FDPart3tmp10 * FDPart3tmp239 - 2 * FDPart3tmp11 * FDPart3tmp239);
        const REAL FDPart3tmp17 = FDPart3tmp16 * FDPart3tmp7;
        const REAL FDPart3tmp29 = FDPart3tmp16 * FDPart3tmp28;
        const REAL FDPart3tmp44 = ((FDPart3tmp12) * (FDPart3tmp12)) * FDPart3tmp43;
        const REAL FDPart3tmp49 = FDPart3tmp15 * FDPart3tmp43;
        const REAL FDPart3tmp51 = FDPart3tmp16 * FDPart3tmp6;
        const REAL FDPart3tmp59 = FDPart3tmp16 * hDD11 + FDPart3tmp16;
        const REAL FDPart3tmp69 = FDPart3tmp6 * FDPart3tmp68;
        const REAL FDPart3tmp75 = FDPart3tmp7 * FDPart3tmp74;
        const REAL FDPart3tmp86 = FDPart3tmp68 * hDD01;
        const REAL FDPart3tmp90 = FDPart3tmp23 * FDPart3tmp68;
        const REAL FDPart3tmp200 = FDPart3tmp23 * (FDPart3tmp16 * hDD_dD110 + FDPart3tmp74 * hDD11 + FDPart3tmp74);
        const REAL FDPart3tmp266 = FDPart3tmp130 * FDPart3tmp16 * aDD11;
        const REAL FDPart3tmp272 = FDPart3tmp31 * (FDPart3tmp16 * aDD_dD110 + FDPart3tmp74 * aDD11);
        const REAL FDPart3tmp18 = FDPart3tmp17 * aDD22;
        const REAL FDPart3tmp24 = FDPart3tmp17 * hDD22 + FDPart3tmp17;
        const REAL FDPart3tmp30 = 2 * FDPart3tmp29 * FDPart3tmp6;
        const REAL FDPart3tmp50 = FDPart3tmp49 * hDD00 + FDPart3tmp49;
        const REAL FDPart3tmp52 = FDPart3tmp51 * hDD12;
        const REAL FDPart3tmp60 = FDPart3tmp44 * ((hDD02) * (hDD02));
        const REAL FDPart3tmp62 = FDPart3tmp44 * ((hDD01) * (hDD01));
        const REAL FDPart3tmp70 = FDPart3tmp69 * hDD02;
        const REAL FDPart3tmp76 = FDPart3tmp23 * (FDPart3tmp17 * hDD_dD220 + FDPart3tmp75 * hDD22 + FDPart3tmp75);
        const REAL FDPart3tmp91 = FDPart3tmp90 * hDD01;
        const REAL FDPart3tmp96 = FDPart3tmp23 * (FDPart3tmp51 * hDD_dD120 + FDPart3tmp6 * FDPart3tmp74 * hDD12);
        const REAL FDPart3tmp100 = FDPart3tmp23 * (FDPart3tmp28 * FDPart3tmp68 * hDD02 + FDPart3tmp69 * hDD_dD021);
        const REAL FDPart3tmp134 = FDPart3tmp31 * (FDPart3tmp17 * aDD_dD220 + FDPart3tmp75 * aDD22);
        const REAL FDPart3tmp160 = FDPart3tmp31 * (FDPart3tmp51 * aDD_dD120 + FDPart3tmp6 * FDPart3tmp74 * aDD12);
        const REAL FDPart3tmp169 = FDPart3tmp31 * (FDPart3tmp28 * FDPart3tmp68 * aDD02 + FDPart3tmp69 * aDD_dD021);
        const REAL FDPart3tmp189 = FDPart3tmp31 * (FDPart3tmp29 * aDD12 + FDPart3tmp51 * aDD_dD121);
        const REAL FDPart3tmp190 = FDPart3tmp23 * (FDPart3tmp29 * hDD12 + FDPart3tmp51 * hDD_dD121);
        const REAL FDPart3tmp243 = FDPart3tmp242 * FDPart3tmp6 + FDPart3tmp49 * FDPart3tmp6;
        const REAL FDPart3tmp267 = FDPart3tmp31 * FDPart3tmp68 * aDD_dD011;
        const REAL FDPart3tmp269 = FDPart3tmp21 * FDPart3tmp68 * aDD01;
        const REAL FDPart3tmp273 = FDPart3tmp38 * (-FDPart3tmp34 * FDPart3tmp86 + FDPart3tmp90 * hDD_dD011);
        const REAL FDPart3tmp283 = FDPart3tmp31 * FDPart3tmp49 * aDD_dD001;
        const REAL FDPart3tmp284 = FDPart3tmp21 * FDPart3tmp49 * aDD00;
        const REAL FDPart3tmp285 = FDPart3tmp130 * FDPart3tmp68 * aDD01;
        const REAL FDPart3tmp25 = FDPart3tmp23 * FDPart3tmp24;
        const REAL FDPart3tmp32 = FDPart3tmp31 * (FDPart3tmp17 * aDD_dD221 + FDPart3tmp30 * aDD22);
        const REAL FDPart3tmp36 = FDPart3tmp23 * (FDPart3tmp17 * hDD_dD221 + FDPart3tmp30 * hDD22 + FDPart3tmp30);
        const REAL FDPart3tmp53 = FDPart3tmp44 * FDPart3tmp47 * FDPart3tmp6 * hDD01 * hDD02 - FDPart3tmp45 * FDPart3tmp50 * FDPart3tmp52;
        const REAL FDPart3tmp71 = -FDPart3tmp45 * FDPart3tmp59 * FDPart3tmp70 + FDPart3tmp47 * FDPart3tmp6 * FDPart3tmp66 * hDD01 * hDD12;
        const REAL FDPart3tmp81 = FDPart3tmp23 * FDPart3tmp51 * hDD12;
        const REAL FDPart3tmp85 = FDPart3tmp24 * FDPart3tmp45;
        const REAL FDPart3tmp92 = 2 * FDPart3tmp80 * FDPart3tmp91 + 2 * FDPart3tmp90 * aDD01;
        const REAL FDPart3tmp109 = 2 * FDPart3tmp16 * FDPart3tmp23 * aDD11 + 2 * FDPart3tmp23 * FDPart3tmp59 * FDPart3tmp80;
        const REAL FDPart3tmp113 = FDPart3tmp6 * FDPart3tmp90 * hDD02;
        const REAL FDPart3tmp166 = FDPart3tmp130 * FDPart3tmp51 * aDD12;
        const REAL FDPart3tmp175 = FDPart3tmp21 * FDPart3tmp69 * aDD02;
        const REAL FDPart3tmp187 = FDPart3tmp21 * FDPart3tmp51 * aDD12;
        const REAL FDPart3tmp191 = FDPart3tmp38 * (FDPart3tmp190 - FDPart3tmp34 * FDPart3tmp52);
        const REAL FDPart3tmp199 = FDPart3tmp33 * FDPart3tmp59 * cf_dD0;
        const REAL FDPart3tmp215 = FDPart3tmp23 * FDPart3tmp49 * hDD_dD001;
        const REAL FDPart3tmp237 = FDPart3tmp130 * FDPart3tmp69 * aDD02;
        const REAL FDPart3tmp244 = FDPart3tmp31 * (FDPart3tmp243 * aDD02 + FDPart3tmp69 * aDD_dD020);
        const REAL FDPart3tmp245 = FDPart3tmp23 * (FDPart3tmp243 * hDD02 + FDPart3tmp69 * hDD_dD020);
        const REAL FDPart3tmp254 = FDPart3tmp23 * (FDPart3tmp68 * hDD_dD010 + hDD01 * (FDPart3tmp242 + FDPart3tmp49));
        const REAL FDPart3tmp271 = FDPart3tmp132 * FDPart3tmp23 * FDPart3tmp59;
        const REAL FDPart3tmp289 = FDPart3tmp31 * (FDPart3tmp68 * aDD_dD010 + aDD01 * (FDPart3tmp242 + FDPart3tmp49));
        const REAL FDPart3tmp37 = -FDPart3tmp24 * FDPart3tmp34 + FDPart3tmp36;
        const REAL FDPart3tmp64 =
            (1.0 / 2.0) /
            (2 * pow(AMPL, 6) * FDPart3tmp43 * FDPart3tmp54 * FDPart3tmp55 * hDD01 * hDD02 * hDD12 / pow(FDPart3tmp14, 6) -
             FDPart3tmp24 * FDPart3tmp46 * FDPart3tmp54 * FDPart3tmp62 + FDPart3tmp24 * FDPart3tmp50 * FDPart3tmp54 * FDPart3tmp59 -
             FDPart3tmp46 * FDPart3tmp50 * FDPart3tmp54 * FDPart3tmp57 - FDPart3tmp46 * FDPart3tmp54 * FDPart3tmp59 * FDPart3tmp60 * FDPart3tmp7);
        const REAL FDPart3tmp82 = 2 * FDPart3tmp23 * FDPart3tmp51 * aDD12 + 2 * FDPart3tmp80 * FDPart3tmp81;
        const REAL FDPart3tmp87 = FDPart3tmp47 * FDPart3tmp66 * FDPart3tmp7 * hDD02 * hDD12 - FDPart3tmp85 * FDPart3tmp86;
        const REAL FDPart3tmp88 = -FDPart3tmp47 * FDPart3tmp57 + FDPart3tmp59 * FDPart3tmp85;
        const REAL FDPart3tmp101 = FDPart3tmp100 - FDPart3tmp34 * FDPart3tmp70;
        const REAL FDPart3tmp105 = -FDPart3tmp47 * FDPart3tmp60 * FDPart3tmp7 + FDPart3tmp50 * FDPart3tmp85;
        const REAL FDPart3tmp114 = 2 * FDPart3tmp113 * FDPart3tmp80 + 2 * FDPart3tmp6 * FDPart3tmp90 * aDD02;
        const REAL FDPart3tmp116 = FDPart3tmp45 * FDPart3tmp50 * FDPart3tmp59 - FDPart3tmp47 * FDPart3tmp62;
        const REAL FDPart3tmp118 = 2 * FDPart3tmp18 * FDPart3tmp23 + 2 * FDPart3tmp25 * FDPart3tmp80;
        const REAL FDPart3tmp135 = -FDPart3tmp24 * FDPart3tmp72 + FDPart3tmp76;
        const REAL FDPart3tmp140 = 2 * FDPart3tmp23 * FDPart3tmp49 * aDD00 + 2 * FDPart3tmp23 * FDPart3tmp50 * FDPart3tmp80;
        const REAL FDPart3tmp144 = -FDPart3tmp52 * FDPart3tmp72 + FDPart3tmp96;
        const REAL FDPart3tmp217 = FDPart3tmp215 - FDPart3tmp34 * FDPart3tmp50;
        const REAL FDPart3tmp219 = -FDPart3tmp199 + FDPart3tmp200;
        const REAL FDPart3tmp246 = FDPart3tmp38 * (FDPart3tmp245 - FDPart3tmp70 * FDPart3tmp72);
        const REAL FDPart3tmp287 = FDPart3tmp23 * FDPart3tmp26 * FDPart3tmp50;
        const REAL FDPart3tmp290 = FDPart3tmp38 * (FDPart3tmp254 - FDPart3tmp72 * FDPart3tmp86);
        const REAL FDPart3tmp65 = FDPart3tmp64 * (FDPart3tmp24 * FDPart3tmp34 - FDPart3tmp36);
        const REAL FDPart3tmp77 = FDPart3tmp64 * (FDPart3tmp24 * FDPart3tmp72 - FDPart3tmp76);
        const REAL FDPart3tmp94 = FDPart3tmp37 * FDPart3tmp64;
        const REAL FDPart3tmp102 = FDPart3tmp64 * (FDPart3tmp101 + FDPart3tmp52 * FDPart3tmp72 - FDPart3tmp96);
        const REAL FDPart3tmp143 = FDPart3tmp135 * FDPart3tmp64;
        const REAL FDPart3tmp145 = FDPart3tmp64 * (-FDPart3tmp100 + FDPart3tmp144 + FDPart3tmp34 * FDPart3tmp70);
        const REAL FDPart3tmp195 = FDPart3tmp64 * (2 * FDPart3tmp190 - FDPart3tmp21 * FDPart3tmp52);
        const REAL FDPart3tmp197 = FDPart3tmp64 * (FDPart3tmp16 * FDPart3tmp23 * hDD_dD111 - FDPart3tmp33 * FDPart3tmp59 * cf_dD1);
        const REAL FDPart3tmp201 =
            FDPart3tmp64 * (FDPart3tmp199 - FDPart3tmp200 - FDPart3tmp21 * FDPart3tmp86 + FDPart3tmp31 * FDPart3tmp68 * hDD_dD011);
        const REAL FDPart3tmp218 = FDPart3tmp217 * FDPart3tmp64;
        const REAL FDPart3tmp220 = FDPart3tmp219 * FDPart3tmp64;
        const REAL FDPart3tmp221 = FDPart3tmp64 * (FDPart3tmp101 + FDPart3tmp144);
        const REAL FDPart3tmp250 = FDPart3tmp64 * (-FDPart3tmp130 * FDPart3tmp70 + 2 * FDPart3tmp245);
        const REAL FDPart3tmp252 =
            FDPart3tmp64 * (FDPart3tmp23 * (FDPart3tmp251 * hDD00 + FDPart3tmp251 + FDPart3tmp49 * hDD_dD000) - FDPart3tmp50 * FDPart3tmp72);
        const REAL FDPart3tmp255 = FDPart3tmp64 * (-FDPart3tmp130 * FDPart3tmp86 - FDPart3tmp215 + 2 * FDPart3tmp254 + FDPart3tmp34 * FDPart3tmp50);
        const REAL FDPart3tmp78 = FDPart3tmp53 * FDPart3tmp65 + FDPart3tmp71 * FDPart3tmp77;
        const REAL FDPart3tmp89 = FDPart3tmp65 * FDPart3tmp87 + FDPart3tmp77 * FDPart3tmp88;
        const REAL FDPart3tmp103 = FDPart3tmp102 * FDPart3tmp87 + FDPart3tmp53 * FDPart3tmp94;
        const REAL FDPart3tmp106 = FDPart3tmp105 * FDPart3tmp65 + FDPart3tmp77 * FDPart3tmp87;
        const REAL FDPart3tmp111 = FDPart3tmp102 * FDPart3tmp88 + FDPart3tmp71 * FDPart3tmp94;
        const REAL FDPart3tmp117 = FDPart3tmp102 * FDPart3tmp71 + FDPart3tmp116 * FDPart3tmp94;
        const REAL FDPart3tmp146 = FDPart3tmp105 * FDPart3tmp145 + FDPart3tmp143 * FDPart3tmp53;
        const REAL FDPart3tmp148 = FDPart3tmp143 * FDPart3tmp71 + FDPart3tmp145 * FDPart3tmp87;
        const REAL FDPart3tmp150 = FDPart3tmp116 * FDPart3tmp143 + FDPart3tmp145 * FDPart3tmp53;
        const REAL FDPart3tmp202 = FDPart3tmp105 * FDPart3tmp197 + FDPart3tmp195 * FDPart3tmp53 + FDPart3tmp201 * FDPart3tmp87;
        const REAL FDPart3tmp204 = FDPart3tmp195 * FDPart3tmp71 + FDPart3tmp197 * FDPart3tmp87 + FDPart3tmp201 * FDPart3tmp88;
        const REAL FDPart3tmp206 = FDPart3tmp116 * FDPart3tmp195 + FDPart3tmp197 * FDPart3tmp53 + FDPart3tmp201 * FDPart3tmp71;
        const REAL FDPart3tmp222 = FDPart3tmp116 * FDPart3tmp221 + FDPart3tmp218 * FDPart3tmp71 + FDPart3tmp220 * FDPart3tmp53;
        const REAL FDPart3tmp224 = FDPart3tmp105 * FDPart3tmp220 + FDPart3tmp218 * FDPart3tmp87 + FDPart3tmp221 * FDPart3tmp53;
        const REAL FDPart3tmp226 = FDPart3tmp218 * FDPart3tmp88 + FDPart3tmp220 * FDPart3tmp87 + FDPart3tmp221 * FDPart3tmp71;
        const REAL FDPart3tmp256 = FDPart3tmp105 * FDPart3tmp255 + FDPart3tmp250 * FDPart3tmp53 + FDPart3tmp252 * FDPart3tmp87;
        const REAL FDPart3tmp258 = FDPart3tmp250 * FDPart3tmp71 + FDPart3tmp252 * FDPart3tmp88 + FDPart3tmp255 * FDPart3tmp87;
        const REAL FDPart3tmp260 = FDPart3tmp116 * FDPart3tmp250 + FDPart3tmp252 * FDPart3tmp71 + FDPart3tmp255 * FDPart3tmp53;
        const REAL FDPart3tmp120 = FDPart3tmp103 * FDPart3tmp82 - FDPart3tmp106 * FDPart3tmp109 + FDPart3tmp111 * FDPart3tmp114 +
                                   FDPart3tmp117 * FDPart3tmp118 + FDPart3tmp18 * FDPart3tmp21 - FDPart3tmp25 * FDPart3tmp26 - FDPart3tmp32 -
                                   FDPart3tmp37 * FDPart3tmp38 - FDPart3tmp78 * FDPart3tmp82 - FDPart3tmp89 * FDPart3tmp92;
        const REAL FDPart3tmp127 = -FDPart3tmp103 * FDPart3tmp82 + FDPart3tmp106 * FDPart3tmp109 - FDPart3tmp111 * FDPart3tmp114 -
                                   FDPart3tmp117 * FDPart3tmp118 - FDPart3tmp18 * FDPart3tmp21 + FDPart3tmp25 * FDPart3tmp26 + FDPart3tmp32 +
                                   FDPart3tmp37 * FDPart3tmp38 + FDPart3tmp78 * FDPart3tmp82 + FDPart3tmp89 * FDPart3tmp92;
        const REAL FDPart3tmp152 = -FDPart3tmp106 * FDPart3tmp92 + FDPart3tmp114 * FDPart3tmp148 - FDPart3tmp114 * FDPart3tmp78 +
                                   FDPart3tmp118 * FDPart3tmp150 + FDPart3tmp130 * FDPart3tmp18 - FDPart3tmp132 * FDPart3tmp25 - FDPart3tmp134 -
                                   FDPart3tmp135 * FDPart3tmp38 - FDPart3tmp140 * FDPart3tmp89 + FDPart3tmp146 * FDPart3tmp82;
        const REAL FDPart3tmp155 = FDPart3tmp106 * FDPart3tmp92 - FDPart3tmp114 * FDPart3tmp148 + FDPart3tmp114 * FDPart3tmp78 -
                                   FDPart3tmp118 * FDPart3tmp150 - FDPart3tmp130 * FDPart3tmp18 + FDPart3tmp132 * FDPart3tmp25 + FDPart3tmp134 +
                                   FDPart3tmp135 * FDPart3tmp38 + FDPart3tmp140 * FDPart3tmp89 - FDPart3tmp146 * FDPart3tmp82;
        const REAL FDPart3tmp168 = -FDPart3tmp103 * FDPart3tmp92 - FDPart3tmp111 * FDPart3tmp140 - FDPart3tmp114 * FDPart3tmp117 -
                                   FDPart3tmp132 * FDPart3tmp81 - FDPart3tmp144 * FDPart3tmp38 - FDPart3tmp160 + FDPart3tmp166;
        const REAL FDPart3tmp177 = FDPart3tmp101 * FDPart3tmp38 + FDPart3tmp109 * FDPart3tmp146 + FDPart3tmp113 * FDPart3tmp26 +
                                   FDPart3tmp148 * FDPart3tmp92 + FDPart3tmp150 * FDPart3tmp82 + FDPart3tmp169 - FDPart3tmp175;
        const REAL FDPart3tmp183 = -FDPart3tmp101 * FDPart3tmp38 - FDPart3tmp109 * FDPart3tmp146 - FDPart3tmp113 * FDPart3tmp26 -
                                   FDPart3tmp148 * FDPart3tmp92 - FDPart3tmp150 * FDPart3tmp82 - FDPart3tmp169 + FDPart3tmp175;
        const REAL FDPart3tmp184 = FDPart3tmp103 * FDPart3tmp92 + FDPart3tmp111 * FDPart3tmp140 + FDPart3tmp114 * FDPart3tmp117 +
                                   FDPart3tmp132 * FDPart3tmp81 + FDPart3tmp144 * FDPart3tmp38 + FDPart3tmp160 - FDPart3tmp166;
        const REAL FDPart3tmp208 = FDPart3tmp103 * FDPart3tmp109 + FDPart3tmp111 * FDPart3tmp92 - FDPart3tmp114 * FDPart3tmp204 +
                                   FDPart3tmp117 * FDPart3tmp82 - FDPart3tmp118 * FDPart3tmp206 - FDPart3tmp187 + FDPart3tmp189 + FDPart3tmp191 -
                                   FDPart3tmp202 * FDPart3tmp82 + FDPart3tmp26 * FDPart3tmp81;
        const REAL FDPart3tmp211 = -FDPart3tmp103 * FDPart3tmp109 - FDPart3tmp111 * FDPart3tmp92 + FDPart3tmp114 * FDPart3tmp204 -
                                   FDPart3tmp117 * FDPart3tmp82 + FDPart3tmp118 * FDPart3tmp206 + FDPart3tmp187 - FDPart3tmp189 - FDPart3tmp191 +
                                   FDPart3tmp202 * FDPart3tmp82 - FDPart3tmp26 * FDPart3tmp81;
        const REAL FDPart3tmp228 = FDPart3tmp114 * FDPart3tmp226 + FDPart3tmp118 * FDPart3tmp222 + FDPart3tmp224 * FDPart3tmp82;
        const REAL FDPart3tmp232 = -FDPart3tmp114 * FDPart3tmp226 - FDPart3tmp118 * FDPart3tmp222 - FDPart3tmp224 * FDPart3tmp82;
        const REAL FDPart3tmp262 = FDPart3tmp113 * FDPart3tmp132 + FDPart3tmp114 * FDPart3tmp150 - FDPart3tmp114 * FDPart3tmp258 -
                                   FDPart3tmp118 * FDPart3tmp260 + FDPart3tmp140 * FDPart3tmp148 + FDPart3tmp146 * FDPart3tmp92 - FDPart3tmp237 +
                                   FDPart3tmp244 + FDPart3tmp246 - FDPart3tmp256 * FDPart3tmp82;
        const REAL FDPart3tmp265 = -FDPart3tmp113 * FDPart3tmp132 - FDPart3tmp114 * FDPart3tmp150 + FDPart3tmp114 * FDPart3tmp258 +
                                   FDPart3tmp118 * FDPart3tmp260 - FDPart3tmp140 * FDPart3tmp148 - FDPart3tmp146 * FDPart3tmp92 + FDPart3tmp237 -
                                   FDPart3tmp244 - FDPart3tmp246 + FDPart3tmp256 * FDPart3tmp82;
        const REAL FDPart3tmp281 = FDPart3tmp109 * FDPart3tmp224 - FDPart3tmp114 * FDPart3tmp206 - FDPart3tmp140 * FDPart3tmp204 -
                                   FDPart3tmp202 * FDPart3tmp92 - FDPart3tmp219 * FDPart3tmp38 + FDPart3tmp222 * FDPart3tmp82 +
                                   FDPart3tmp226 * FDPart3tmp92 + FDPart3tmp26 * FDPart3tmp91 + FDPart3tmp266 + FDPart3tmp267 - FDPart3tmp269 -
                                   FDPart3tmp271 - FDPart3tmp272 + FDPart3tmp273;
        const REAL FDPart3tmp282 = -FDPart3tmp109 * FDPart3tmp224 + FDPart3tmp114 * FDPart3tmp206 + FDPart3tmp140 * FDPart3tmp204 +
                                   FDPart3tmp202 * FDPart3tmp92 + FDPart3tmp219 * FDPart3tmp38 - FDPart3tmp222 * FDPart3tmp82 -
                                   FDPart3tmp226 * FDPart3tmp92 - FDPart3tmp26 * FDPart3tmp91 - FDPart3tmp266 - FDPart3tmp267 + FDPart3tmp269 +
                                   FDPart3tmp271 + FDPart3tmp272 - FDPart3tmp273;
        const REAL FDPart3tmp297 = -FDPart3tmp109 * FDPart3tmp256 + FDPart3tmp114 * FDPart3tmp222 + FDPart3tmp132 * FDPart3tmp91 +
                                   FDPart3tmp140 * FDPart3tmp226 - FDPart3tmp217 * FDPart3tmp38 + FDPart3tmp224 * FDPart3tmp92 -
                                   FDPart3tmp258 * FDPart3tmp92 - FDPart3tmp260 * FDPart3tmp82 - FDPart3tmp283 + FDPart3tmp284 - FDPart3tmp285 -
                                   FDPart3tmp287 + FDPart3tmp289 + FDPart3tmp290;
        const REAL FDPart3tmp298 = FDPart3tmp109 * FDPart3tmp256 - FDPart3tmp114 * FDPart3tmp222 - FDPart3tmp132 * FDPart3tmp91 -
                                   FDPart3tmp140 * FDPart3tmp226 + FDPart3tmp217 * FDPart3tmp38 - FDPart3tmp224 * FDPart3tmp92 +
                                   FDPart3tmp258 * FDPart3tmp92 + FDPart3tmp260 * FDPart3tmp82 + FDPart3tmp283 - FDPart3tmp284 + FDPart3tmp285 +
                                   FDPart3tmp287 - FDPart3tmp289 - FDPart3tmp290;
        diagnostic_output_gfs[IDX4(PSI4_PART1REGF, i0, i1, i2)] =
            FDPart3tmp120 * (FDPart3tmp1 * mre4U2 + FDPart3tmp2 * FDPart3tmp3 - FDPart3tmp2 * FDPart3tmp5 - FDPart3tmp4 * mim4U2) +
            FDPart3tmp127 * (FDPart3tmp121 * FDPart3tmp122 - FDPart3tmp122 * FDPart3tmp125 - FDPart3tmp124 * ((mim4U3) * (mim4U3)) +
                             FDPart3tmp124 * ((mre4U3) * (mre4U3))) +
            FDPart3tmp152 * (FDPart3tmp1 * mre4U1 + FDPart3tmp128 * FDPart3tmp2 - FDPart3tmp129 * FDPart3tmp2 - FDPart3tmp4 * mim4U1) +
            FDPart3tmp155 * (FDPart3tmp121 * FDPart3tmp153 - FDPart3tmp125 * FDPart3tmp153 - FDPart3tmp154 * ((mim4U3) * (mim4U3)) +
                             FDPart3tmp154 * ((mre4U3) * (mre4U3))) +
            FDPart3tmp208 *
                (FDPart3tmp121 * FDPart3tmp186 - FDPart3tmp124 * mim4U2 * mim4U3 + FDPart3tmp124 * mre4U2 * mre4U3 - FDPart3tmp125 * FDPart3tmp186) +
            FDPart3tmp211 * (-FDPart3tmp0 * ((mim4U2) * (mim4U2)) + FDPart3tmp0 * ((mre4U2) * (mre4U2)) + FDPart3tmp122 * FDPart3tmp3 -
                             FDPart3tmp122 * FDPart3tmp5) +
            FDPart3tmp262 *
                (FDPart3tmp121 * FDPart3tmp236 - FDPart3tmp125 * FDPart3tmp236 - FDPart3tmp154 * mim4U1 * mim4U3 + FDPart3tmp154 * mre4U1 * mre4U3) +
            FDPart3tmp265 * (-FDPart3tmp0 * ((mim4U1) * (mim4U1)) + FDPart3tmp0 * ((mre4U1) * (mre4U1)) + FDPart3tmp128 * FDPart3tmp153 -
                             FDPart3tmp129 * FDPart3tmp153) +
            FDPart3tmp281 *
                (FDPart3tmp124 * FDPart3tmp212 - FDPart3tmp124 * FDPart3tmp213 + FDPart3tmp128 * FDPart3tmp186 - FDPart3tmp129 * FDPart3tmp186) +
            FDPart3tmp282 * (-FDPart3tmp154 * ((mim4U2) * (mim4U2)) + FDPart3tmp154 * ((mre4U2) * (mre4U2)) + FDPart3tmp230 * FDPart3tmp3 -
                             FDPart3tmp230 * FDPart3tmp5) +
            FDPart3tmp297 *
                (FDPart3tmp154 * FDPart3tmp212 - FDPart3tmp154 * FDPart3tmp213 + FDPart3tmp236 * FDPart3tmp3 - FDPart3tmp236 * FDPart3tmp5) +
            FDPart3tmp298 * (-FDPart3tmp124 * ((mim4U1) * (mim4U1)) + FDPart3tmp124 * ((mre4U1) * (mre4U1)) + FDPart3tmp128 * FDPart3tmp230 -
                             FDPart3tmp129 * FDPart3tmp230) +
            (FDPart3tmp156 + FDPart3tmp159) * (FDPart3tmp168 + FDPart3tmp177) + (FDPart3tmp156 + FDPart3tmp214) * (FDPart3tmp168 + FDPart3tmp228) +
            (FDPart3tmp159 + FDPart3tmp231) * (FDPart3tmp177 + FDPart3tmp232) + (FDPart3tmp179 + FDPart3tmp182) * (FDPart3tmp183 + FDPart3tmp184) +
            (FDPart3tmp179 + FDPart3tmp214) * (FDPart3tmp183 + FDPart3tmp228) + (FDPart3tmp182 + FDPart3tmp231) * (FDPart3tmp184 + FDPart3tmp232);
        diagnostic_output_gfs[IDX4(PSI4_PART1IMGF, i0, i1, i2)] =
            FDPart3tmp120 * (-FDPart3tmp1 * mim4U2 + FDPart3tmp2 * FDPart3tmp303 + FDPart3tmp2 * FDPart3tmp304 - FDPart3tmp4 * mre4U2) +
            FDPart3tmp127 * (FDPart3tmp122 * FDPart3tmp299 + FDPart3tmp122 * FDPart3tmp300 - 2 * FDPart3tmp124 * mim4U3 * mre4U3) +
            FDPart3tmp152 * (-FDPart3tmp1 * mim4U1 + FDPart3tmp2 * FDPart3tmp306 + FDPart3tmp2 * FDPart3tmp307 - FDPart3tmp4 * mre4U1) +
            FDPart3tmp155 * (FDPart3tmp153 * FDPart3tmp299 + FDPart3tmp153 * FDPart3tmp300 - 2 * FDPart3tmp154 * mim4U3 * mre4U3) +
            FDPart3tmp208 *
                (-FDPart3tmp124 * mim4U2 * mre4U3 - FDPart3tmp124 * mim4U3 * mre4U2 + FDPart3tmp186 * FDPart3tmp299 + FDPart3tmp186 * FDPart3tmp300) +
            FDPart3tmp211 * (-2 * FDPart3tmp0 * mim4U2 * mre4U2 + FDPart3tmp122 * FDPart3tmp303 + FDPart3tmp122 * FDPart3tmp304) +
            FDPart3tmp262 *
                (-FDPart3tmp154 * mim4U1 * mre4U3 - FDPart3tmp154 * mim4U3 * mre4U1 + FDPart3tmp236 * FDPart3tmp299 + FDPart3tmp236 * FDPart3tmp300) +
            FDPart3tmp265 * (-2 * FDPart3tmp0 * mim4U1 * mre4U1 + FDPart3tmp153 * FDPart3tmp306 + FDPart3tmp153 * FDPart3tmp307) +
            FDPart3tmp281 *
                (-FDPart3tmp124 * FDPart3tmp318 - FDPart3tmp124 * FDPart3tmp319 + FDPart3tmp186 * FDPart3tmp306 + FDPart3tmp186 * FDPart3tmp307) +
            FDPart3tmp282 * (-2 * FDPart3tmp154 * mim4U2 * mre4U2 + FDPart3tmp230 * FDPart3tmp303 + FDPart3tmp230 * FDPart3tmp304) +
            FDPart3tmp297 *
                (-FDPart3tmp154 * FDPart3tmp318 - FDPart3tmp154 * FDPart3tmp319 + FDPart3tmp236 * FDPart3tmp303 + FDPart3tmp236 * FDPart3tmp304) +
            FDPart3tmp298 * (-2 * FDPart3tmp124 * mim4U1 * mre4U1 + FDPart3tmp230 * FDPart3tmp306 + FDPart3tmp230 * FDPart3tmp307) +
            (FDPart3tmp168 + FDPart3tmp177) * (FDPart3tmp308 + FDPart3tmp311) + (FDPart3tmp168 + FDPart3tmp228) * (FDPart3tmp308 + FDPart3tmp320) +
            (FDPart3tmp177 + FDPart3tmp232) * (FDPart3tmp311 + FDPart3tmp321) + (FDPart3tmp183 + FDPart3tmp184) * (FDPart3tmp312 + FDPart3tmp315) +
            (FDPart3tmp183 + FDPart3tmp228) * (FDPart3tmp312 + FDPart3tmp320) + (FDPart3tmp184 + FDPart3tmp232) * (FDPart3tmp315 + FDPart3tmp321);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
