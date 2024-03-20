#include "../BHaH_defines.h"
#include "../simd/simd_intrinsics.h"
/*
 * Finite difference function for operator dD0, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_dD0_fdorder8(const REAL_SIMD_ARRAY FDPROTO_i0m1, const REAL_SIMD_ARRAY FDPROTO_i0m2,
                                                     const REAL_SIMD_ARRAY FDPROTO_i0m3, const REAL_SIMD_ARRAY FDPROTO_i0m4,
                                                     const REAL_SIMD_ARRAY FDPROTO_i0p1, const REAL_SIMD_ARRAY FDPROTO_i0p2,
                                                     const REAL_SIMD_ARRAY FDPROTO_i0p3, const REAL_SIMD_ARRAY FDPROTO_i0p4,
                                                     const REAL_SIMD_ARRAY invdxx0) {
  const double dblFDPart1_Rational_1_280 = 1.0 / 280.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_280 = ConstSIMD(dblFDPart1_Rational_1_280);

  const double dblFDPart1_Rational_1_5 = 1.0 / 5.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_5 = ConstSIMD(dblFDPart1_Rational_1_5);

  const double dblFDPart1_Rational_4_105 = 4.0 / 105.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_105 = ConstSIMD(dblFDPart1_Rational_4_105);

  const double dblFDPart1_Rational_4_5 = 4.0 / 5.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_5 = ConstSIMD(dblFDPart1_Rational_4_5);

  const REAL_SIMD_ARRAY FD_result =
      MulSIMD(invdxx0, FusedMulAddSIMD(FDPart1_Rational_4_105, SubSIMD(FDPROTO_i0p3, FDPROTO_i0m3),
                                       FusedMulAddSIMD(FDPart1_Rational_4_5, SubSIMD(FDPROTO_i0p1, FDPROTO_i0m1),
                                                       FusedMulAddSIMD(FDPart1_Rational_1_280, SubSIMD(FDPROTO_i0m4, FDPROTO_i0p4),
                                                                       MulSIMD(FDPart1_Rational_1_5, SubSIMD(FDPROTO_i0m2, FDPROTO_i0p2))))));

  return FD_result;
}
/*
 * Finite difference function for operator dD1, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_dD1_fdorder8(const REAL_SIMD_ARRAY FDPROTO_i1m1, const REAL_SIMD_ARRAY FDPROTO_i1m2,
                                                     const REAL_SIMD_ARRAY FDPROTO_i1m3, const REAL_SIMD_ARRAY FDPROTO_i1m4,
                                                     const REAL_SIMD_ARRAY FDPROTO_i1p1, const REAL_SIMD_ARRAY FDPROTO_i1p2,
                                                     const REAL_SIMD_ARRAY FDPROTO_i1p3, const REAL_SIMD_ARRAY FDPROTO_i1p4,
                                                     const REAL_SIMD_ARRAY invdxx1) {
  const double dblFDPart1_Rational_1_280 = 1.0 / 280.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_280 = ConstSIMD(dblFDPart1_Rational_1_280);

  const double dblFDPart1_Rational_1_5 = 1.0 / 5.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_5 = ConstSIMD(dblFDPart1_Rational_1_5);

  const double dblFDPart1_Rational_4_105 = 4.0 / 105.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_105 = ConstSIMD(dblFDPart1_Rational_4_105);

  const double dblFDPart1_Rational_4_5 = 4.0 / 5.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_5 = ConstSIMD(dblFDPart1_Rational_4_5);

  const REAL_SIMD_ARRAY FD_result =
      MulSIMD(invdxx1, FusedMulAddSIMD(FDPart1_Rational_4_105, SubSIMD(FDPROTO_i1p3, FDPROTO_i1m3),
                                       FusedMulAddSIMD(FDPart1_Rational_4_5, SubSIMD(FDPROTO_i1p1, FDPROTO_i1m1),
                                                       FusedMulAddSIMD(FDPart1_Rational_1_280, SubSIMD(FDPROTO_i1m4, FDPROTO_i1p4),
                                                                       MulSIMD(FDPart1_Rational_1_5, SubSIMD(FDPROTO_i1m2, FDPROTO_i1p2))))));

  return FD_result;
}
/*
 * Finite difference function for operator dDD00, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_dDD00_fdorder8(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i0m1,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0m2, const REAL_SIMD_ARRAY FDPROTO_i0m3,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0m4, const REAL_SIMD_ARRAY FDPROTO_i0p1,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0p2, const REAL_SIMD_ARRAY FDPROTO_i0p3,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0p4, const REAL_SIMD_ARRAY invdxx0) {
  const double dblFDPart1_Rational_1_5 = 1.0 / 5.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_5 = ConstSIMD(dblFDPart1_Rational_1_5);

  const double dblFDPart1_Rational_1_560 = 1.0 / 560.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_560 = ConstSIMD(dblFDPart1_Rational_1_560);

  const double dblFDPart1_Rational_205_72 = 205.0 / 72.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_205_72 = ConstSIMD(dblFDPart1_Rational_205_72);

  const double dblFDPart1_Rational_8_315 = 8.0 / 315.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_8_315 = ConstSIMD(dblFDPart1_Rational_8_315);

  const double dblFDPart1_Rational_8_5 = 8.0 / 5.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_8_5 = ConstSIMD(dblFDPart1_Rational_8_5);

  const REAL_SIMD_ARRAY FD_result =
      MulSIMD(MulSIMD(invdxx0, invdxx0),
              FusedMulAddSIMD(FDPart1_Rational_8_5, AddSIMD(FDPROTO_i0m1, FDPROTO_i0p1),
                              FusedMulSubSIMD(FDPart1_Rational_8_315, AddSIMD(FDPROTO_i0m3, FDPROTO_i0p3),
                                              FusedMulAddSIMD(FDPart1_Rational_1_5, AddSIMD(FDPROTO_i0m2, FDPROTO_i0p2),
                                                              FusedMulAddSIMD(FDPart1_Rational_1_560, AddSIMD(FDPROTO_i0m4, FDPROTO_i0p4),
                                                                              MulSIMD(FDPROTO, FDPart1_Rational_205_72))))));

  return FD_result;
}
/*
 * Finite difference function for operator dDD01, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_dDD01_fdorder8(
    const REAL_SIMD_ARRAY FDPROTO_i0m1_i1m1, const REAL_SIMD_ARRAY FDPROTO_i0m1_i1m2, const REAL_SIMD_ARRAY FDPROTO_i0m1_i1m3,
    const REAL_SIMD_ARRAY FDPROTO_i0m1_i1m4, const REAL_SIMD_ARRAY FDPROTO_i0m1_i1p1, const REAL_SIMD_ARRAY FDPROTO_i0m1_i1p2,
    const REAL_SIMD_ARRAY FDPROTO_i0m1_i1p3, const REAL_SIMD_ARRAY FDPROTO_i0m1_i1p4, const REAL_SIMD_ARRAY FDPROTO_i0m2_i1m1,
    const REAL_SIMD_ARRAY FDPROTO_i0m2_i1m2, const REAL_SIMD_ARRAY FDPROTO_i0m2_i1m3, const REAL_SIMD_ARRAY FDPROTO_i0m2_i1m4,
    const REAL_SIMD_ARRAY FDPROTO_i0m2_i1p1, const REAL_SIMD_ARRAY FDPROTO_i0m2_i1p2, const REAL_SIMD_ARRAY FDPROTO_i0m2_i1p3,
    const REAL_SIMD_ARRAY FDPROTO_i0m2_i1p4, const REAL_SIMD_ARRAY FDPROTO_i0m3_i1m1, const REAL_SIMD_ARRAY FDPROTO_i0m3_i1m2,
    const REAL_SIMD_ARRAY FDPROTO_i0m3_i1m3, const REAL_SIMD_ARRAY FDPROTO_i0m3_i1m4, const REAL_SIMD_ARRAY FDPROTO_i0m3_i1p1,
    const REAL_SIMD_ARRAY FDPROTO_i0m3_i1p2, const REAL_SIMD_ARRAY FDPROTO_i0m3_i1p3, const REAL_SIMD_ARRAY FDPROTO_i0m3_i1p4,
    const REAL_SIMD_ARRAY FDPROTO_i0m4_i1m1, const REAL_SIMD_ARRAY FDPROTO_i0m4_i1m2, const REAL_SIMD_ARRAY FDPROTO_i0m4_i1m3,
    const REAL_SIMD_ARRAY FDPROTO_i0m4_i1m4, const REAL_SIMD_ARRAY FDPROTO_i0m4_i1p1, const REAL_SIMD_ARRAY FDPROTO_i0m4_i1p2,
    const REAL_SIMD_ARRAY FDPROTO_i0m4_i1p3, const REAL_SIMD_ARRAY FDPROTO_i0m4_i1p4, const REAL_SIMD_ARRAY FDPROTO_i0p1_i1m1,
    const REAL_SIMD_ARRAY FDPROTO_i0p1_i1m2, const REAL_SIMD_ARRAY FDPROTO_i0p1_i1m3, const REAL_SIMD_ARRAY FDPROTO_i0p1_i1m4,
    const REAL_SIMD_ARRAY FDPROTO_i0p1_i1p1, const REAL_SIMD_ARRAY FDPROTO_i0p1_i1p2, const REAL_SIMD_ARRAY FDPROTO_i0p1_i1p3,
    const REAL_SIMD_ARRAY FDPROTO_i0p1_i1p4, const REAL_SIMD_ARRAY FDPROTO_i0p2_i1m1, const REAL_SIMD_ARRAY FDPROTO_i0p2_i1m2,
    const REAL_SIMD_ARRAY FDPROTO_i0p2_i1m3, const REAL_SIMD_ARRAY FDPROTO_i0p2_i1m4, const REAL_SIMD_ARRAY FDPROTO_i0p2_i1p1,
    const REAL_SIMD_ARRAY FDPROTO_i0p2_i1p2, const REAL_SIMD_ARRAY FDPROTO_i0p2_i1p3, const REAL_SIMD_ARRAY FDPROTO_i0p2_i1p4,
    const REAL_SIMD_ARRAY FDPROTO_i0p3_i1m1, const REAL_SIMD_ARRAY FDPROTO_i0p3_i1m2, const REAL_SIMD_ARRAY FDPROTO_i0p3_i1m3,
    const REAL_SIMD_ARRAY FDPROTO_i0p3_i1m4, const REAL_SIMD_ARRAY FDPROTO_i0p3_i1p1, const REAL_SIMD_ARRAY FDPROTO_i0p3_i1p2,
    const REAL_SIMD_ARRAY FDPROTO_i0p3_i1p3, const REAL_SIMD_ARRAY FDPROTO_i0p3_i1p4, const REAL_SIMD_ARRAY FDPROTO_i0p4_i1m1,
    const REAL_SIMD_ARRAY FDPROTO_i0p4_i1m2, const REAL_SIMD_ARRAY FDPROTO_i0p4_i1m3, const REAL_SIMD_ARRAY FDPROTO_i0p4_i1m4,
    const REAL_SIMD_ARRAY FDPROTO_i0p4_i1p1, const REAL_SIMD_ARRAY FDPROTO_i0p4_i1p2, const REAL_SIMD_ARRAY FDPROTO_i0p4_i1p3,
    const REAL_SIMD_ARRAY FDPROTO_i0p4_i1p4, const REAL_SIMD_ARRAY invdxx0, const REAL_SIMD_ARRAY invdxx1) {
  const double dblFDPart1_Rational_16_11025 = 16.0 / 11025.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_16_11025 = ConstSIMD(dblFDPart1_Rational_16_11025);

  const double dblFDPart1_Rational_16_25 = 16.0 / 25.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_16_25 = ConstSIMD(dblFDPart1_Rational_16_25);

  const double dblFDPart1_Rational_16_525 = 16.0 / 525.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_16_525 = ConstSIMD(dblFDPart1_Rational_16_525);

  const double dblFDPart1_Rational_1_1400 = 1.0 / 1400.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_1400 = ConstSIMD(dblFDPart1_Rational_1_1400);

  const double dblFDPart1_Rational_1_25 = 1.0 / 25.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_25 = ConstSIMD(dblFDPart1_Rational_1_25);

  const double dblFDPart1_Rational_1_350 = 1.0 / 350.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_350 = ConstSIMD(dblFDPart1_Rational_1_350);

  const double dblFDPart1_Rational_1_7350 = 1.0 / 7350.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_7350 = ConstSIMD(dblFDPart1_Rational_1_7350);

  const double dblFDPart1_Rational_1_78400 = 1.0 / 78400.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_78400 = ConstSIMD(dblFDPart1_Rational_1_78400);

  const double dblFDPart1_Rational_4_25 = 4.0 / 25.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_25 = ConstSIMD(dblFDPart1_Rational_4_25);

  const double dblFDPart1_Rational_4_525 = 4.0 / 525.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_525 = ConstSIMD(dblFDPart1_Rational_4_525);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx0,
      MulSIMD(
          invdxx1,
          FusedMulAddSIMD(
              FDPart1_Rational_1_7350,
              AddSIMD(AddSIMD(FDPROTO_i0m4_i1p3, FDPROTO_i0p3_i1m4),
                      AddSIMD(FDPROTO_i0p4_i1m3, SubSIMD(FDPROTO_i0m3_i1p4, AddSIMD(AddSIMD(FDPROTO_i0m3_i1m4, FDPROTO_i0m4_i1m3),
                                                                                    AddSIMD(FDPROTO_i0p3_i1p4, FDPROTO_i0p4_i1p3))))),
              FusedMulAddSIMD(
                  FDPart1_Rational_1_78400, AddSIMD(FDPROTO_i0p4_i1p4, SubSIMD(FDPROTO_i0m4_i1m4, AddSIMD(FDPROTO_i0m4_i1p4, FDPROTO_i0p4_i1m4))),
                  FusedMulAddSIMD(
                      FDPart1_Rational_1_25, AddSIMD(FDPROTO_i0p2_i1p2, SubSIMD(FDPROTO_i0m2_i1m2, AddSIMD(FDPROTO_i0m2_i1p2, FDPROTO_i0p2_i1m2))),
                      FusedMulAddSIMD(
                          FDPart1_Rational_1_350,
                          AddSIMD(AddSIMD(FDPROTO_i0m4_i1p1, FDPROTO_i0p1_i1m4),
                                  AddSIMD(FDPROTO_i0p4_i1m1, SubSIMD(FDPROTO_i0m1_i1p4, AddSIMD(AddSIMD(FDPROTO_i0m1_i1m4, FDPROTO_i0m4_i1m1),
                                                                                                AddSIMD(FDPROTO_i0p1_i1p4, FDPROTO_i0p4_i1p1))))),
                          FusedMulAddSIMD(
                              FDPart1_Rational_16_525,
                              AddSIMD(AddSIMD(FDPROTO_i0m3_i1m1, FDPROTO_i0p1_i1p3),
                                      AddSIMD(FDPROTO_i0p3_i1p1, SubSIMD(FDPROTO_i0m1_i1m3, AddSIMD(AddSIMD(FDPROTO_i0m1_i1p3, FDPROTO_i0m3_i1p1),
                                                                                                    AddSIMD(FDPROTO_i0p1_i1m3, FDPROTO_i0p3_i1m1))))),
                              FusedMulAddSIMD(
                                  FDPart1_Rational_1_1400,
                                  AddSIMD(
                                      AddSIMD(FDPROTO_i0m4_i1m2, FDPROTO_i0p2_i1p4),
                                      AddSIMD(FDPROTO_i0p4_i1p2, SubSIMD(FDPROTO_i0m2_i1m4, AddSIMD(AddSIMD(FDPROTO_i0m2_i1p4, FDPROTO_i0m4_i1p2),
                                                                                                    AddSIMD(FDPROTO_i0p2_i1m4, FDPROTO_i0p4_i1m2))))),
                                  FusedMulAddSIMD(
                                      FDPart1_Rational_4_25,
                                      AddSIMD(AddSIMD(FDPROTO_i0m2_i1p1, FDPROTO_i0p1_i1m2),
                                              AddSIMD(FDPROTO_i0p2_i1m1,
                                                      SubSIMD(FDPROTO_i0m1_i1p2, AddSIMD(AddSIMD(FDPROTO_i0m1_i1m2, FDPROTO_i0m2_i1m1),
                                                                                         AddSIMD(FDPROTO_i0p1_i1p2, FDPROTO_i0p2_i1p1))))),
                                      FusedMulAddSIMD(
                                          FDPart1_Rational_4_525,
                                          AddSIMD(AddSIMD(FDPROTO_i0m3_i1p2, FDPROTO_i0p2_i1m3),
                                                  AddSIMD(FDPROTO_i0p3_i1m2,
                                                          SubSIMD(FDPROTO_i0m2_i1p3, AddSIMD(AddSIMD(FDPROTO_i0m2_i1m3, FDPROTO_i0m3_i1m2),
                                                                                             AddSIMD(FDPROTO_i0p2_i1p3, FDPROTO_i0p3_i1p2))))),
                                          FusedMulAddSIMD(
                                              FDPart1_Rational_16_11025,
                                              AddSIMD(FDPROTO_i0p3_i1p3, SubSIMD(FDPROTO_i0m3_i1m3, AddSIMD(FDPROTO_i0m3_i1p3, FDPROTO_i0p3_i1m3))),
                                              MulSIMD(FDPart1_Rational_16_25,
                                                      AddSIMD(FDPROTO_i0p1_i1p1,
                                                              SubSIMD(FDPROTO_i0m1_i1m1, AddSIMD(FDPROTO_i0m1_i1p1, FDPROTO_i0p1_i1m1)))))))))))))));

  return FD_result;
}
/*
 * Finite difference function for operator dDD11, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_dDD11_fdorder8(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i1m1,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1m2, const REAL_SIMD_ARRAY FDPROTO_i1m3,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1m4, const REAL_SIMD_ARRAY FDPROTO_i1p1,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1p2, const REAL_SIMD_ARRAY FDPROTO_i1p3,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1p4, const REAL_SIMD_ARRAY invdxx1) {
  const double dblFDPart1_Rational_1_5 = 1.0 / 5.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_5 = ConstSIMD(dblFDPart1_Rational_1_5);

  const double dblFDPart1_Rational_1_560 = 1.0 / 560.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_560 = ConstSIMD(dblFDPart1_Rational_1_560);

  const double dblFDPart1_Rational_205_72 = 205.0 / 72.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_205_72 = ConstSIMD(dblFDPart1_Rational_205_72);

  const double dblFDPart1_Rational_8_315 = 8.0 / 315.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_8_315 = ConstSIMD(dblFDPart1_Rational_8_315);

  const double dblFDPart1_Rational_8_5 = 8.0 / 5.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_8_5 = ConstSIMD(dblFDPart1_Rational_8_5);

  const REAL_SIMD_ARRAY FD_result =
      MulSIMD(MulSIMD(invdxx1, invdxx1),
              FusedMulAddSIMD(FDPart1_Rational_8_5, AddSIMD(FDPROTO_i1m1, FDPROTO_i1p1),
                              FusedMulSubSIMD(FDPart1_Rational_8_315, AddSIMD(FDPROTO_i1m3, FDPROTO_i1p3),
                                              FusedMulAddSIMD(FDPart1_Rational_1_5, AddSIMD(FDPROTO_i1m2, FDPROTO_i1p2),
                                                              FusedMulAddSIMD(FDPart1_Rational_1_560, AddSIMD(FDPROTO_i1m4, FDPROTO_i1p4),
                                                                              MulSIMD(FDPROTO, FDPart1_Rational_205_72))))));

  return FD_result;
}

/*
 * Evaluate BSSN constraints.
 */
void constraints_eval__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                          const rfm_struct *restrict rfmstruct, const REAL *restrict in_gfs, const REAL *restrict auxevol_gfs,
                                          REAL *restrict diagnostic_output_gfs) {
#include "../set_CodeParameters-simd.h"
#pragma omp parallel for collapse(2)
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
      const double NOSIMDf1_of_xx1 = rfmstruct->f1_of_xx1[i1];
      const REAL_SIMD_ARRAY f1_of_xx1 = ConstSIMD(NOSIMDf1_of_xx1);
      const double NOSIMDf1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
      const REAL_SIMD_ARRAY f1_of_xx1__D1 = ConstSIMD(NOSIMDf1_of_xx1__D1);
      const double NOSIMDf1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];
      const REAL_SIMD_ARRAY f1_of_xx1__DD11 = ConstSIMD(NOSIMDf1_of_xx1__DD11);

      for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0 += simd_width) {
        const REAL_SIMD_ARRAY f0_of_xx0 = ReadSIMD(&rfmstruct->f0_of_xx0[i0]);
        const REAL_SIMD_ARRAY f0_of_xx0__D0 = ReadSIMD(&rfmstruct->f0_of_xx0__D0[i0]);
        const REAL_SIMD_ARRAY f0_of_xx0__DD00 = ReadSIMD(&rfmstruct->f0_of_xx0__DD00[i0]);
        const REAL_SIMD_ARRAY f0_of_xx0__DDD000 = ReadSIMD(&rfmstruct->f0_of_xx0__DDD000[i0]);
        /*
         * NRPy+-Generated GF Access/FD Code, Step 1 of 2:
         * Read gridfunction(s) from main memory and compute FD stencils as needed.
         */
        const REAL_SIMD_ARRAY RbarDD00 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD01 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD02 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD11 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD12 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD22 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY cf_i1m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY cf_i1m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY cf = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY cf_i1p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY cf_i1p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY cf_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY trK_i1m4 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY trK_i1m3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY trK_i1m2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY trK_i1m1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m4 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY trK = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p4 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i1p1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY trK_i1p2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY trK_i1p3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY trK_i1p4 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD_dD000 =
            SIMD_fd_function_dD0_fdorder8(aDD00_i0m1, aDD00_i0m2, aDD00_i0m3, aDD00_i0m4, aDD00_i0p1, aDD00_i0p2, aDD00_i0p3, aDD00_i0p4, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD001 =
            SIMD_fd_function_dD1_fdorder8(aDD00_i1m1, aDD00_i1m2, aDD00_i1m3, aDD00_i1m4, aDD00_i1p1, aDD00_i1p2, aDD00_i1p3, aDD00_i1p4, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD010 =
            SIMD_fd_function_dD0_fdorder8(aDD01_i0m1, aDD01_i0m2, aDD01_i0m3, aDD01_i0m4, aDD01_i0p1, aDD01_i0p2, aDD01_i0p3, aDD01_i0p4, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD011 =
            SIMD_fd_function_dD1_fdorder8(aDD01_i1m1, aDD01_i1m2, aDD01_i1m3, aDD01_i1m4, aDD01_i1p1, aDD01_i1p2, aDD01_i1p3, aDD01_i1p4, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD020 =
            SIMD_fd_function_dD0_fdorder8(aDD02_i0m1, aDD02_i0m2, aDD02_i0m3, aDD02_i0m4, aDD02_i0p1, aDD02_i0p2, aDD02_i0p3, aDD02_i0p4, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD021 =
            SIMD_fd_function_dD1_fdorder8(aDD02_i1m1, aDD02_i1m2, aDD02_i1m3, aDD02_i1m4, aDD02_i1p1, aDD02_i1p2, aDD02_i1p3, aDD02_i1p4, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD110 =
            SIMD_fd_function_dD0_fdorder8(aDD11_i0m1, aDD11_i0m2, aDD11_i0m3, aDD11_i0m4, aDD11_i0p1, aDD11_i0p2, aDD11_i0p3, aDD11_i0p4, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD111 =
            SIMD_fd_function_dD1_fdorder8(aDD11_i1m1, aDD11_i1m2, aDD11_i1m3, aDD11_i1m4, aDD11_i1p1, aDD11_i1p2, aDD11_i1p3, aDD11_i1p4, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD120 =
            SIMD_fd_function_dD0_fdorder8(aDD12_i0m1, aDD12_i0m2, aDD12_i0m3, aDD12_i0m4, aDD12_i0p1, aDD12_i0p2, aDD12_i0p3, aDD12_i0p4, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD121 =
            SIMD_fd_function_dD1_fdorder8(aDD12_i1m1, aDD12_i1m2, aDD12_i1m3, aDD12_i1m4, aDD12_i1p1, aDD12_i1p2, aDD12_i1p3, aDD12_i1p4, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD220 =
            SIMD_fd_function_dD0_fdorder8(aDD22_i0m1, aDD22_i0m2, aDD22_i0m3, aDD22_i0m4, aDD22_i0p1, aDD22_i0p2, aDD22_i0p3, aDD22_i0p4, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD221 =
            SIMD_fd_function_dD1_fdorder8(aDD22_i1m1, aDD22_i1m2, aDD22_i1m3, aDD22_i1m4, aDD22_i1p1, aDD22_i1p2, aDD22_i1p3, aDD22_i1p4, invdxx1);
        const REAL_SIMD_ARRAY cf_dD0 = SIMD_fd_function_dD0_fdorder8(cf_i0m1, cf_i0m2, cf_i0m3, cf_i0m4, cf_i0p1, cf_i0p2, cf_i0p3, cf_i0p4, invdxx0);
        const REAL_SIMD_ARRAY cf_dD1 = SIMD_fd_function_dD1_fdorder8(cf_i1m1, cf_i1m2, cf_i1m3, cf_i1m4, cf_i1p1, cf_i1p2, cf_i1p3, cf_i1p4, invdxx1);
        const REAL_SIMD_ARRAY cf_dDD00 =
            SIMD_fd_function_dDD00_fdorder8(cf, cf_i0m1, cf_i0m2, cf_i0m3, cf_i0m4, cf_i0p1, cf_i0p2, cf_i0p3, cf_i0p4, invdxx0);
        const REAL_SIMD_ARRAY cf_dDD01 = SIMD_fd_function_dDD01_fdorder8(
            cf_i0m1_i1m1, cf_i0m1_i1m2, cf_i0m1_i1m3, cf_i0m1_i1m4, cf_i0m1_i1p1, cf_i0m1_i1p2, cf_i0m1_i1p3, cf_i0m1_i1p4, cf_i0m2_i1m1,
            cf_i0m2_i1m2, cf_i0m2_i1m3, cf_i0m2_i1m4, cf_i0m2_i1p1, cf_i0m2_i1p2, cf_i0m2_i1p3, cf_i0m2_i1p4, cf_i0m3_i1m1, cf_i0m3_i1m2,
            cf_i0m3_i1m3, cf_i0m3_i1m4, cf_i0m3_i1p1, cf_i0m3_i1p2, cf_i0m3_i1p3, cf_i0m3_i1p4, cf_i0m4_i1m1, cf_i0m4_i1m2, cf_i0m4_i1m3,
            cf_i0m4_i1m4, cf_i0m4_i1p1, cf_i0m4_i1p2, cf_i0m4_i1p3, cf_i0m4_i1p4, cf_i0p1_i1m1, cf_i0p1_i1m2, cf_i0p1_i1m3, cf_i0p1_i1m4,
            cf_i0p1_i1p1, cf_i0p1_i1p2, cf_i0p1_i1p3, cf_i0p1_i1p4, cf_i0p2_i1m1, cf_i0p2_i1m2, cf_i0p2_i1m3, cf_i0p2_i1m4, cf_i0p2_i1p1,
            cf_i0p2_i1p2, cf_i0p2_i1p3, cf_i0p2_i1p4, cf_i0p3_i1m1, cf_i0p3_i1m2, cf_i0p3_i1m3, cf_i0p3_i1m4, cf_i0p3_i1p1, cf_i0p3_i1p2,
            cf_i0p3_i1p3, cf_i0p3_i1p4, cf_i0p4_i1m1, cf_i0p4_i1m2, cf_i0p4_i1m3, cf_i0p4_i1m4, cf_i0p4_i1p1, cf_i0p4_i1p2, cf_i0p4_i1p3,
            cf_i0p4_i1p4, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY cf_dDD11 =
            SIMD_fd_function_dDD11_fdorder8(cf, cf_i1m1, cf_i1m2, cf_i1m3, cf_i1m4, cf_i1p1, cf_i1p2, cf_i1p3, cf_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD000 =
            SIMD_fd_function_dD0_fdorder8(hDD00_i0m1, hDD00_i0m2, hDD00_i0m3, hDD00_i0m4, hDD00_i0p1, hDD00_i0p2, hDD00_i0p3, hDD00_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD001 =
            SIMD_fd_function_dD1_fdorder8(hDD00_i1m1, hDD00_i1m2, hDD00_i1m3, hDD00_i1m4, hDD00_i1p1, hDD00_i1p2, hDD00_i1p3, hDD00_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD010 =
            SIMD_fd_function_dD0_fdorder8(hDD01_i0m1, hDD01_i0m2, hDD01_i0m3, hDD01_i0m4, hDD01_i0p1, hDD01_i0p2, hDD01_i0p3, hDD01_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD011 =
            SIMD_fd_function_dD1_fdorder8(hDD01_i1m1, hDD01_i1m2, hDD01_i1m3, hDD01_i1m4, hDD01_i1p1, hDD01_i1p2, hDD01_i1p3, hDD01_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD020 =
            SIMD_fd_function_dD0_fdorder8(hDD02_i0m1, hDD02_i0m2, hDD02_i0m3, hDD02_i0m4, hDD02_i0p1, hDD02_i0p2, hDD02_i0p3, hDD02_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD021 =
            SIMD_fd_function_dD1_fdorder8(hDD02_i1m1, hDD02_i1m2, hDD02_i1m3, hDD02_i1m4, hDD02_i1p1, hDD02_i1p2, hDD02_i1p3, hDD02_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD110 =
            SIMD_fd_function_dD0_fdorder8(hDD11_i0m1, hDD11_i0m2, hDD11_i0m3, hDD11_i0m4, hDD11_i0p1, hDD11_i0p2, hDD11_i0p3, hDD11_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD111 =
            SIMD_fd_function_dD1_fdorder8(hDD11_i1m1, hDD11_i1m2, hDD11_i1m3, hDD11_i1m4, hDD11_i1p1, hDD11_i1p2, hDD11_i1p3, hDD11_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD120 =
            SIMD_fd_function_dD0_fdorder8(hDD12_i0m1, hDD12_i0m2, hDD12_i0m3, hDD12_i0m4, hDD12_i0p1, hDD12_i0p2, hDD12_i0p3, hDD12_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD121 =
            SIMD_fd_function_dD1_fdorder8(hDD12_i1m1, hDD12_i1m2, hDD12_i1m3, hDD12_i1m4, hDD12_i1p1, hDD12_i1p2, hDD12_i1p3, hDD12_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD220 =
            SIMD_fd_function_dD0_fdorder8(hDD22_i0m1, hDD22_i0m2, hDD22_i0m3, hDD22_i0m4, hDD22_i0p1, hDD22_i0p2, hDD22_i0p3, hDD22_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD221 =
            SIMD_fd_function_dD1_fdorder8(hDD22_i1m1, hDD22_i1m2, hDD22_i1m3, hDD22_i1m4, hDD22_i1p1, hDD22_i1p2, hDD22_i1p3, hDD22_i1p4, invdxx1);
        const REAL_SIMD_ARRAY trK_dD0 =
            SIMD_fd_function_dD0_fdorder8(trK_i0m1, trK_i0m2, trK_i0m3, trK_i0m4, trK_i0p1, trK_i0p2, trK_i0p3, trK_i0p4, invdxx0);
        const REAL_SIMD_ARRAY trK_dD1 =
            SIMD_fd_function_dD1_fdorder8(trK_i1m1, trK_i1m2, trK_i1m3, trK_i1m4, trK_i1p1, trK_i1p2, trK_i1p3, trK_i1p4, invdxx1);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        const double dblFDPart3_Integer_1 = 1.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(dblFDPart3_Integer_1);

        const double dblFDPart3_Integer_12 = 12.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_12 = ConstSIMD(dblFDPart3_Integer_12);

        const double dblFDPart3_Integer_16 = 16.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_16 = ConstSIMD(dblFDPart3_Integer_16);

        const double dblFDPart3_Integer_2 = 2.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(dblFDPart3_Integer_2);

        const double dblFDPart3_Integer_4 = 4.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_4 = ConstSIMD(dblFDPart3_Integer_4);

        const double dblFDPart3_Integer_6 = 6.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_6 = ConstSIMD(dblFDPart3_Integer_6);

        const double dblFDPart3_Integer_8 = 8.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_8 = ConstSIMD(dblFDPart3_Integer_8);

        const double dblFDPart3_NegativeOne_ = -1.0;
        const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(dblFDPart3_NegativeOne_);

        const double dblFDPart3_Rational_1_2 = 1.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_2 = ConstSIMD(dblFDPart3_Rational_1_2);

        const double dblFDPart3_Rational_2_3 = 2.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_2_3 = ConstSIMD(dblFDPart3_Rational_2_3);

        const REAL_SIMD_ARRAY FDPart3tmp0 = MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp2 = MulSIMD(f1_of_xx1, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp6 = MulSIMD(f0_of_xx0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp17 = MulSIMD(f0_of_xx0__D0, MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0));
        const REAL_SIMD_ARRAY FDPart3tmp18 = MulSIMD(f0_of_xx0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp107 = DivSIMD(FDPart3_Integer_1, cf);
        const REAL_SIMD_ARRAY FDPart3tmp152 = MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf);
        const REAL_SIMD_ARRAY FDPart3tmp3 = MulSIMD(FDPart3tmp2, MulSIMD(MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0), f0_of_xx0));
        const REAL_SIMD_ARRAY FDPart3tmp4 = FusedMulAddSIMD(FDPart3tmp0, hDD00, FDPart3tmp0);
        const REAL_SIMD_ARRAY FDPart3tmp7 = FusedMulAddSIMD(FDPart3tmp6, hDD11, FDPart3tmp6);
        const REAL_SIMD_ARRAY FDPart3tmp8 = MulSIMD(FDPart3tmp2, FDPart3tmp6);
        const REAL_SIMD_ARRAY FDPart3tmp19 = MulSIMD(FDPart3tmp18, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp22 = MulSIMD(FDPart3tmp0, aDD00);
        const REAL_SIMD_ARRAY FDPart3tmp24 = MulSIMD(FDPart3tmp6, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp27 = MulSIMD(FDPart3tmp6, aDD11);
        const REAL_SIMD_ARRAY FDPart3tmp34 = MulSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, f0_of_xx0__D0));
        const REAL_SIMD_ARRAY FDPart3tmp69 = MulSIMD(FDPart3tmp18, aDD01);
        const REAL_SIMD_ARRAY FDPart3tmp103 = DivSIMD(FDPart3_Integer_1, MulSIMD(cf, cf));
        const REAL_SIMD_ARRAY FDPart3tmp132 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp107);
        const REAL_SIMD_ARRAY FDPart3tmp146 = MulSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0__D0, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp148 = FusedMulAddSIMD(f0_of_xx0, f0_of_xx0__DD00, FDPart3tmp0);
        const REAL_SIMD_ARRAY FDPart3tmp5 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp3, MulSIMD(hDD12, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp10 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp0), MulSIMD(FDPart3tmp8, MulSIMD(hDD02, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp11 = FusedMulAddSIMD(FDPart3tmp8, hDD22, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp12 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp0), MulSIMD(FDPart3tmp6, MulSIMD(hDD01, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp20 =
            FusedMulSubSIMD(MulSIMD(FDPart3tmp17, f1_of_xx1), MulSIMD(hDD01, hDD12), MulSIMD(FDPart3tmp19, MulSIMD(FDPart3tmp7, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp25 = FusedMulSubSIMD(f1_of_xx1, MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp6), MulSIMD(hDD01, hDD02)),
                                                             MulSIMD(FDPart3tmp24, MulSIMD(FDPart3tmp4, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp30 = MulSIMD(FDPart3tmp8, aDD22);
        const REAL_SIMD_ARRAY FDPart3tmp35 = MulSIMD(FDPart3tmp34, aDD01);
        const REAL_SIMD_ARRAY FDPart3tmp38 = MulSIMD(FDPart3tmp34, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp40 = MulSIMD(FDPart3_Integer_2, FDPart3tmp24);
        const REAL_SIMD_ARRAY FDPart3tmp65 = MulSIMD(FDPart3tmp19, aDD02);
        const REAL_SIMD_ARRAY FDPart3tmp67 = MulSIMD(FDPart3tmp24, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp109 = MulSIMD(FDPart3tmp2, FDPart3tmp34);
        const REAL_SIMD_ARRAY FDPart3tmp114 = MulSIMD(FDPart3tmp18, MulSIMD(f1_of_xx1__D1, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp134 = FusedMulAddSIMD(FDPart3tmp34, hDD11, FusedMulAddSIMD(FDPart3tmp6, hDD_dD110, FDPart3tmp34));
        const REAL_SIMD_ARRAY FDPart3tmp140 =
            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, f0_of_xx0__D0),
                               FusedMulAddSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, hDD11)),
                                               FusedMulSubSIMD(FDPart3tmp34, hDD_dD011, MulSIMD(FDPart3tmp6, hDD_dD110))));
        const REAL_SIMD_ARRAY FDPart3tmp144 = FusedMulAddSIMD(FDPart3tmp0, f1_of_xx1, MulSIMD(f0_of_xx0, MulSIMD(f0_of_xx0__DD00, f1_of_xx1)));
        const REAL_SIMD_ARRAY FDPart3tmp147 = FusedMulAddSIMD(FDPart3tmp0, hDD_dD000, FusedMulAddSIMD(FDPart3tmp146, hDD00, FDPart3tmp146));
        const REAL_SIMD_ARRAY FDPart3tmp149 = FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp148, hDD01),
                                                              FusedMulSubSIMD(FDPart3tmp34, hDD_dD010, MulSIMD(FDPart3tmp0, hDD_dD001)));
        const REAL_SIMD_ARRAY FDPart3tmp28 = FusedMulAddSIMD(FDPart3tmp4, FDPart3tmp7, FDPart3tmp12);
        const REAL_SIMD_ARRAY FDPart3tmp39 = MulSIMD(FDPart3tmp38, aDD02);
        const REAL_SIMD_ARRAY FDPart3tmp41 = MulSIMD(FDPart3tmp40, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp44 =
            FusedMulSubSIMD(MulSIMD(FDPart3tmp17, FDPart3tmp2), MulSIMD(hDD02, hDD12), MulSIMD(FDPart3tmp11, MulSIMD(FDPart3tmp18, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp57 = FusedMulAddSIMD(FDPart3tmp11, FDPart3tmp4, FDPart3tmp10);
        const REAL_SIMD_ARRAY FDPart3tmp111 = FusedMulAddSIMD(FDPart3tmp109, hDD22, FusedMulAddSIMD(FDPart3tmp8, hDD_dD220, FDPart3tmp109));
        const REAL_SIMD_ARRAY FDPart3tmp116 = FusedMulAddSIMD(FDPart3tmp24, hDD_dD120, MulSIMD(FDPart3tmp34, MulSIMD(f1_of_xx1, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp120 = MulSIMD(FDPart3tmp40, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp123 = FusedMulAddSIMD(FDPart3tmp19, hDD_dD021, FDPart3tmp114);
        const REAL_SIMD_ARRAY FDPart3tmp128 =
            FusedMulAddSIMD(MulSIMD(FDPart3tmp6, f1_of_xx1), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1__D1, hDD22)),
                            FusedMulSubSIMD(FDPart3tmp6, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, f1_of_xx1__D1)),
                                            MulSIMD(FDPart3tmp8, hDD_dD221)));
        const REAL_SIMD_ARRAY FDPart3tmp129 =
            FusedMulAddSIMD(MulSIMD(FDPart3tmp2, f0_of_xx0), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, hDD22)),
                            FusedMulSubSIMD(FDPart3tmp2, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0)),
                                            MulSIMD(FDPart3tmp8, hDD_dD220)));
        const REAL_SIMD_ARRAY FDPart3tmp139 =
            FusedMulAddSIMD(FDPart3tmp40, hDD_dD121, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp6), MulSIMD(f1_of_xx1__D1, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp145 = FusedMulAddSIMD(FDPart3tmp38, hDD_dD020, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp144, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp14 = FusedMulAddSIMD(
            FDPart3tmp4, FDPart3tmp5,
            FusedMulAddSIMD(FDPart3tmp11, MulSIMD(FDPart3tmp4, FDPart3tmp7),
                            FusedMulAddSIMD(MulSIMD(FDPart3tmp3, hDD01), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp0), MulSIMD(hDD02, hDD12)),
                                            FusedMulAddSIMD(FDPart3tmp10, FDPart3tmp7, MulSIMD(FDPart3tmp11, FDPart3tmp12)))));
        const REAL_SIMD_ARRAY FDPart3tmp47 = FusedMulAddSIMD(FDPart3tmp11, FDPart3tmp7, FDPart3tmp5);
        const REAL_SIMD_ARRAY FDPart3tmp117 = SubSIMD(SubSIMD(FDPart3tmp116, FDPart3tmp114), MulSIMD(FDPart3tmp19, hDD_dD021));
        const REAL_SIMD_ARRAY FDPart3tmp122 = FusedMulAddSIMD(FDPart3tmp120, hDD22, FusedMulAddSIMD(FDPart3tmp8, hDD_dD221, FDPart3tmp120));
        const REAL_SIMD_ARRAY FDPart3tmp124 =
            FusedMulAddSIMD(MulSIMD(f0_of_xx0, f0_of_xx0__D0), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, hDD12)),
                            NegFusedMulAddSIMD(FDPart3tmp24, hDD_dD120, FDPart3tmp123));
        const REAL_SIMD_ARRAY FDPart3tmp135 = AddSIMD(FDPart3tmp116, FDPart3tmp123);
        const REAL_SIMD_ARRAY FDPart3tmp15 = DivSIMD(FDPart3_Integer_1, MulSIMD(FDPart3tmp14, FDPart3tmp14));
        const REAL_SIMD_ARRAY FDPart3tmp101 = DivSIMD(FDPart3_Integer_1, FDPart3tmp14);
        const REAL_SIMD_ARRAY FDPart3tmp21 = MulSIMD(FDPart3tmp15, MulSIMD(FDPart3tmp20, FDPart3tmp20));
        const REAL_SIMD_ARRAY FDPart3tmp26 = MulSIMD(FDPart3tmp15, MulSIMD(FDPart3tmp25, FDPart3tmp25));
        const REAL_SIMD_ARRAY FDPart3tmp31 = MulSIMD(FDPart3tmp15, FDPart3tmp20);
        const REAL_SIMD_ARRAY FDPart3tmp36 = MulSIMD(FDPart3tmp15, FDPart3tmp28);
        const REAL_SIMD_ARRAY FDPart3tmp45 = MulSIMD(FDPart3tmp15, MulSIMD(FDPart3tmp44, FDPart3tmp44));
        const REAL_SIMD_ARRAY FDPart3tmp48 = MulSIMD(FDPart3tmp15, MulSIMD(FDPart3tmp47, FDPart3tmp47));
        const REAL_SIMD_ARRAY FDPart3tmp52 = MulSIMD(FDPart3tmp15, FDPart3tmp47);
        const REAL_SIMD_ARRAY FDPart3tmp58 = MulSIMD(FDPart3tmp15, MulSIMD(FDPart3tmp57, FDPart3tmp57));
        const REAL_SIMD_ARRAY FDPart3tmp61 = MulSIMD(FDPart3tmp15, FDPart3tmp57);
        const REAL_SIMD_ARRAY FDPart3tmp102 = MulSIMD(FDPart3_Integer_2, FDPart3tmp101);
        const REAL_SIMD_ARRAY FDPart3tmp118 =
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp117, FDPart3tmp44)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp111, FDPart3tmp20)));
        const REAL_SIMD_ARRAY FDPart3tmp119 =
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp117, FDPart3tmp57)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp111, FDPart3tmp25)));
        const REAL_SIMD_ARRAY FDPart3tmp125 =
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp124, FDPart3tmp44)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp122, FDPart3tmp25)));
        const REAL_SIMD_ARRAY FDPart3tmp126 =
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp124, FDPart3tmp47)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp122, FDPart3tmp20)));
        const REAL_SIMD_ARRAY FDPart3tmp130 =
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp129, FDPart3tmp47)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp128, FDPart3tmp44)));
        const REAL_SIMD_ARRAY FDPart3tmp131 =
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp129, FDPart3tmp44)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp128, FDPart3tmp57)));
        const REAL_SIMD_ARRAY FDPart3tmp136 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp0, FDPart3tmp101), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp47, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp135, FDPart3tmp20)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp134, FDPart3tmp44))));
        const REAL_SIMD_ARRAY FDPart3tmp137 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp0, FDPart3tmp101), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp44, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp135, FDPart3tmp25)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp134, FDPart3tmp57))));
        const REAL_SIMD_ARRAY FDPart3tmp141 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp101, FDPart3tmp44), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp6, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp140, FDPart3tmp47)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp139, FDPart3tmp20))));
        const REAL_SIMD_ARRAY FDPart3tmp142 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp101, FDPart3tmp57), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp6, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp140, FDPart3tmp44)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp139, FDPart3tmp25))));
        const REAL_SIMD_ARRAY FDPart3tmp150 = FusedMulAddSIMD(
            FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp149, FDPart3tmp44)),
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp147, FDPart3tmp47)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp145, FDPart3tmp20))));
        const REAL_SIMD_ARRAY FDPart3tmp151 = FusedMulAddSIMD(
            FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp149, FDPart3tmp57)),
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp147, FDPart3tmp44)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp145, FDPart3tmp25))));
        const REAL_SIMD_ARRAY FDPart3tmp157 =
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp117, FDPart3tmp25)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp111, FDPart3tmp28)));
        const REAL_SIMD_ARRAY FDPart3tmp160 =
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp124, FDPart3tmp20)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp122, FDPart3tmp28)));
        const REAL_SIMD_ARRAY FDPart3tmp168 =
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp129, FDPart3tmp20)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp128, FDPart3tmp25)));
        const REAL_SIMD_ARRAY FDPart3tmp172 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp0, FDPart3tmp101), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp20, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp135, FDPart3tmp28)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp134, FDPart3tmp25))));
        const REAL_SIMD_ARRAY FDPart3tmp177 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp101, FDPart3tmp25), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp6, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp140, FDPart3tmp20)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp139, FDPart3tmp28))));
        const REAL_SIMD_ARRAY FDPart3tmp180 = FusedMulAddSIMD(
            FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp149, FDPart3tmp25)),
            FusedMulSubSIMD(FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp147, FDPart3tmp20)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp101), MulSIMD(FDPart3tmp145, FDPart3tmp28))));
        const REAL_SIMD_ARRAY FDPart3tmp32 = MulSIMD(FDPart3tmp25, FDPart3tmp31);
        const REAL_SIMD_ARRAY FDPart3tmp37 = MulSIMD(FDPart3tmp20, FDPart3tmp36);
        const REAL_SIMD_ARRAY FDPart3tmp51 = MulSIMD(FDPart3tmp15, MulSIMD(FDPart3tmp20, FDPart3tmp44));
        const REAL_SIMD_ARRAY FDPart3tmp53 = MulSIMD(FDPart3tmp20, FDPart3tmp52);
        const REAL_SIMD_ARRAY FDPart3tmp54 = MulSIMD(FDPart3tmp44, FDPart3tmp52);
        const REAL_SIMD_ARRAY FDPart3tmp60 = MulSIMD(FDPart3tmp15, MulSIMD(FDPart3tmp25, FDPart3tmp44));
        const REAL_SIMD_ARRAY FDPart3tmp62 = MulSIMD(FDPart3tmp25, FDPart3tmp61);
        const REAL_SIMD_ARRAY FDPart3tmp63 = MulSIMD(FDPart3tmp44, FDPart3tmp61);
        const REAL_SIMD_ARRAY FDPart3tmp74 = MulSIMD(FDPart3tmp25, FDPart3tmp52);
        const REAL_SIMD_ARRAY FDPart3tmp76 = MulSIMD(FDPart3tmp36, FDPart3tmp44);
        const REAL_SIMD_ARRAY FDPart3tmp78 = MulSIMD(FDPart3tmp36, FDPart3tmp47);
        const REAL_SIMD_ARRAY FDPart3tmp86 = MulSIMD(FDPart3tmp20, FDPart3tmp61);
        const REAL_SIMD_ARRAY FDPart3tmp89 = MulSIMD(FDPart3tmp36, MulSIMD(FDPart3tmp57, FDPart3tmp67));
        const REAL_SIMD_ARRAY FDPart3tmp98 = MulSIMD(FDPart3tmp52, FDPart3tmp57);
        const REAL_SIMD_ARRAY FDPart3tmp169 =
            MulSIMD(FDPart3tmp15,
                    FusedMulAddSIMD(FDPart3tmp131, FDPart3tmp41,
                                    FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp168, FDPart3tmp30), MulSIMD(FDPart3tmp130, FDPart3tmp39))));
        const REAL_SIMD_ARRAY FDPart3tmp178 =
            FusedMulAddSIMD(FDPart3tmp6, aDD_dD111,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp142, FDPart3tmp27),
                                            FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp35, MulSIMD(FDPart3tmp177, FDPart3tmp41))));
        const REAL_SIMD_ARRAY FDPart3tmp181 =
            FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp35,
                            FusedMulAddSIMD(FDPart3tmp180, FDPart3tmp39,
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp150, FDPart3tmp22),
                                                            FusedMulAddSIMD(FDPart3tmp0, aDD_dD000, MulSIMD(FDPart3tmp146, aDD00)))));
        const REAL_SIMD_ARRAY FDPart3tmp182 =
            FusedMulAddSIMD(FDPart3tmp126, FDPart3tmp22, FusedMulAddSIMD(FDPart3tmp160, FDPart3tmp65, MulSIMD(FDPart3tmp125, FDPart3tmp69)));
        const REAL_SIMD_ARRAY FDPart3tmp183 =
            FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp27, FusedMulAddSIMD(FDPart3tmp157, FDPart3tmp67, MulSIMD(FDPart3tmp118, FDPart3tmp69)));
        const REAL_SIMD_ARRAY FDPart3tmp187 =
            FusedMulAddSIMD(FDPart3tmp137, FDPart3tmp67, FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp30, MulSIMD(FDPart3tmp136, FDPart3tmp65)));
        const REAL_SIMD_ARRAY FDPart3tmp84 = MulSIMD(FDPart3tmp25, MulSIMD(FDPart3tmp30, FDPart3tmp36));
        const REAL_SIMD_ARRAY FDPart3tmp159 =
            FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp41,
                            FusedMulAddSIMD(FDPart3tmp8, aDD_dD220,
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp157, FDPart3tmp30),
                                                            FusedMulAddSIMD(FDPart3tmp109, aDD22, MulSIMD(FDPart3tmp118, FDPart3tmp39)))));
        const REAL_SIMD_ARRAY FDPart3tmp162 =
            FusedMulAddSIMD(FDPart3tmp160, FDPart3tmp41,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp125, FDPart3tmp27), MulSIMD(FDPart3tmp126, FDPart3tmp35)));
        const REAL_SIMD_ARRAY FDPart3tmp164 =
            FusedMulAddSIMD(FDPart3tmp157, FDPart3tmp39,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp118, FDPart3tmp22), MulSIMD(FDPart3tmp119, FDPart3tmp35)));
        const REAL_SIMD_ARRAY FDPart3tmp166 =
            FusedMulAddSIMD(FDPart3tmp126, FDPart3tmp39,
                            FusedMulAddSIMD(FDPart3tmp8, aDD_dD221,
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp160, FDPart3tmp30),
                                                            FusedMulAddSIMD(FDPart3tmp120, aDD22, MulSIMD(FDPart3tmp125, FDPart3tmp41)))));
        const REAL_SIMD_ARRAY FDPart3tmp174 =
            FusedMulAddSIMD(FDPart3tmp34, aDD11,
                            FusedMulAddSIMD(FDPart3tmp6, aDD_dD110,
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp137, FDPart3tmp27),
                                                            FusedMulAddSIMD(FDPart3tmp136, FDPart3tmp35, MulSIMD(FDPart3tmp172, FDPart3tmp41)))));
        const REAL_SIMD_ARRAY FDPart3tmp176 =
            FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp39,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp136, FDPart3tmp22),
                                            FusedMulAddSIMD(FDPart3tmp0, aDD_dD001, MulSIMD(FDPart3tmp137, FDPart3tmp35))));
        const REAL_SIMD_ARRAY FDPart3tmp184 = AddSIMD(FDPart3tmp182, FDPart3tmp183);
        const REAL_SIMD_ARRAY FDPart3tmp185 = FusedMulAddSIMD(
            FDPart3tmp130, FDPart3tmp22,
            FusedMulAddSIMD(FDPart3tmp131, FDPart3tmp69,
                            FusedMulAddSIMD(FDPart3tmp157, FDPart3tmp30,
                                            FusedMulAddSIMD(FDPart3tmp168, FDPart3tmp65,
                                                            FusedMulAddSIMD(FDPart3tmp118, FDPart3tmp65, MulSIMD(FDPart3tmp119, FDPart3tmp67))))));
        const REAL_SIMD_ARRAY FDPart3tmp186 = FusedMulAddSIMD(
            FDPart3tmp130, FDPart3tmp69,
            FusedMulAddSIMD(FDPart3tmp131, FDPart3tmp27,
                            FusedMulAddSIMD(FDPart3tmp160, FDPart3tmp30,
                                            FusedMulAddSIMD(FDPart3tmp168, FDPart3tmp67,
                                                            FusedMulAddSIMD(FDPart3tmp125, FDPart3tmp67, MulSIMD(FDPart3tmp126, FDPart3tmp65))))));
        const REAL_SIMD_ARRAY FDPart3tmp188 = FusedMulAddSIMD(
            FDPart3tmp19, aDD_dD021, FusedMulAddSIMD(FDPart3tmp18, MulSIMD(aDD02, f1_of_xx1__D1), AddSIMD(FDPart3tmp182, FDPart3tmp187)));
        const REAL_SIMD_ARRAY FDPart3tmp189 =
            FusedMulAddSIMD(FDPart3tmp24, aDD_dD120, FusedMulAddSIMD(FDPart3tmp38, aDD12, AddSIMD(FDPart3tmp183, FDPart3tmp187)));
        const REAL_SIMD_ARRAY FDPart3tmp190 = FusedMulAddSIMD(
            FDPart3tmp160, FDPart3tmp67,
            FusedMulAddSIMD(
                FDPart3tmp177, FDPart3tmp30,
                FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp65,
                                FusedMulAddSIMD(FDPart3tmp142, FDPart3tmp67,
                                                FusedMulAddSIMD(FDPart3tmp24, aDD_dD121,
                                                                FusedMulAddSIMD(FDPart3tmp6, MulSIMD(aDD12, f1_of_xx1__D1),
                                                                                FusedMulAddSIMD(FDPart3tmp125, FDPart3tmp27,
                                                                                                MulSIMD(FDPart3tmp126, FDPart3tmp69))))))));
        const REAL_SIMD_ARRAY FDPart3tmp192 = FusedMulAddSIMD(
            FDPart3tmp151, FDPart3tmp67,
            FusedMulAddSIMD(
                FDPart3tmp157, FDPart3tmp65,
                FusedMulAddSIMD(FDPart3tmp144, aDD02,
                                FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp65,
                                                FusedMulAddSIMD(FDPart3tmp180, FDPart3tmp30,
                                                                FusedMulAddSIMD(FDPart3tmp19, aDD_dD020,
                                                                                FusedMulAddSIMD(FDPart3tmp118, FDPart3tmp22,
                                                                                                MulSIMD(FDPart3tmp119, FDPart3tmp69))))))));
        const REAL_SIMD_ARRAY FDPart3tmp193 = FusedMulAddSIMD(
            FDPart3tmp172, FDPart3tmp67,
            FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp22,
                            FusedMulAddSIMD(FDPart3tmp142, FDPart3tmp69,
                                            FusedMulAddSIMD(FDPart3tmp177, FDPart3tmp65,
                                                            FusedMulAddSIMD(FDPart3tmp18, aDD_dD011,
                                                                            FusedMulAddSIMD(FDPart3tmp136, FDPart3tmp69,
                                                                                            MulSIMD(FDPart3tmp137, FDPart3tmp27)))))));
        const REAL_SIMD_ARRAY FDPart3tmp194 = FusedMulAddSIMD(
            FDPart3tmp151, FDPart3tmp27,
            FusedMulAddSIMD(
                FDPart3tmp172, FDPart3tmp65,
                FusedMulAddSIMD(FDPart3tmp148, aDD01,
                                FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp69,
                                                FusedMulAddSIMD(FDPart3tmp18, aDD_dD010,
                                                                FusedMulAddSIMD(FDPart3tmp180, FDPart3tmp67,
                                                                                FusedMulAddSIMD(FDPart3tmp136, FDPart3tmp22,
                                                                                                MulSIMD(FDPart3tmp137, FDPart3tmp69))))))));
        const REAL_SIMD_ARRAY FDPart3tmp171 = FusedMulAddSIMD(
            FDPart3_Integer_6, MulSIMD(FDPart3tmp65, FDPart3tmp74),
            FusedMulAddSIMD(
                FDPart3_Integer_6, MulSIMD(FDPart3tmp51, FDPart3tmp65),
                FusedMulAddSIMD(
                    FDPart3_Integer_6, MulSIMD(FDPart3tmp60, FDPart3tmp67),
                    FusedMulAddSIMD(
                        FDPart3_Integer_6, MulSIMD(FDPart3tmp30, FDPart3tmp32),
                        FusedMulAddSIMD(
                            FDPart3_Integer_6, MulSIMD(FDPart3tmp45, FDPart3tmp69),
                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp67, FDPart3tmp86),
                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp69, FDPart3tmp98),
                                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp22, FDPart3tmp54),
                                                                            MulSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp27, FDPart3tmp63))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp191 = MulSIMD(FDPart3tmp190, FDPart3tmp57);
        const REAL_SIMD_ARRAY FDPart3tmp195 = FusedMulAddSIMD(
            FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp44, trK_dD1)),
            FusedMulAddSIMD(
                FDPart3tmp178, MulSIMD(FDPart3tmp44, FDPart3tmp61),
                FusedMulAddSIMD(
                    FDPart3tmp107,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(cf_dD0,
                                    FusedMulAddSIMD(
                                        FDPart3_Integer_12, MulSIMD(FDPart3tmp54, FDPart3tmp69),
                                        FusedMulAddSIMD(
                                            FDPart3_Integer_6, MulSIMD(FDPart3tmp21, FDPart3tmp30),
                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp22, FDPart3tmp48),
                                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp27, FDPart3tmp45),
                                                                            FusedMulAddSIMD(FDPart3_Integer_12, MulSIMD(FDPart3tmp51, FDPart3tmp67),
                                                                                            MulSIMD(FDPart3_Integer_12,
                                                                                                    MulSIMD(FDPart3tmp53, FDPart3tmp65))))))))),
                    FusedMulAddSIMD(
                        FDPart3tmp166, MulSIMD(FDPart3tmp25, FDPart3tmp31),
                        FusedMulAddSIMD(
                            FDPart3tmp169, MulSIMD(FDPart3tmp20, FDPart3tmp28),
                            FusedMulAddSIMD(
                                FDPart3_Integer_2, MulSIMD(FDPart3tmp192, FDPart3tmp53),
                                FusedMulAddSIMD(
                                    FDPart3_Integer_2, MulSIMD(FDPart3tmp194, FDPart3tmp54),
                                    FusedMulAddSIMD(
                                        FDPart3tmp193, FDPart3tmp98,
                                        FusedMulAddSIMD(
                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp189, FDPart3tmp51),
                                            FusedMulAddSIMD(
                                                FDPart3tmp191, FDPart3tmp31,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp193, FDPart3tmp45,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp188, FDPart3tmp74,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp190, FDPart3tmp60,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp186, FDPart3tmp76,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp188, FDPart3tmp51,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp185, FDPart3tmp78,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp186, FDPart3tmp32,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp184, FDPart3tmp74,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp185, FDPart3tmp21,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp181, FDPart3tmp48,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp184, FDPart3tmp51,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp174, FDPart3tmp45,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp176, FDPart3tmp54,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp162, FDPart3tmp60,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp164, FDPart3tmp53,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp101,
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3_NegativeOne_,
                                                                                                                            FDPart3_Rational_2_3),
                                                                                                                    MulSIMD(FDPart3tmp47, trK_dD0)),
                                                                                                                FusedMulSubSIMD(
                                                                                                                    FDPart3tmp159, FDPart3tmp21,
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3_Rational_1_2,
                                                                                                                                FDPart3tmp107),
                                                                                                                        MulSIMD(
                                                                                                                            FDPart3tmp171,
                                                                                                                            cf_dD1)))))))))))))))))))))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp199 = FusedMulAddSIMD(
            FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp44, trK_dD0)),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp15, FDPart3tmp159), MulSIMD(FDPart3tmp20, FDPart3tmp25),
                FusedMulAddSIMD(
                    FDPart3tmp107,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(cf_dD1,
                                    FusedMulAddSIMD(
                                        FDPart3_Integer_12, MulSIMD(FDPart3tmp63, FDPart3tmp69),
                                        FusedMulAddSIMD(
                                            FDPart3_Integer_6, MulSIMD(FDPart3tmp22, FDPart3tmp45),
                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp26, FDPart3tmp30),
                                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp27, FDPart3tmp58),
                                                                            FusedMulAddSIMD(FDPart3_Integer_12, MulSIMD(FDPart3tmp60, FDPart3tmp65),
                                                                                            MulSIMD(FDPart3_Integer_12,
                                                                                                    MulSIMD(FDPart3tmp62, FDPart3tmp67))))))))),
                    FusedMulAddSIMD(
                        FDPart3tmp186, MulSIMD(FDPart3tmp36, FDPart3tmp57),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3_Integer_2, FDPart3tmp15), MulSIMD(FDPart3tmp191, FDPart3tmp25),
                            FusedMulAddSIMD(
                                FDPart3_Integer_2, MulSIMD(FDPart3tmp193, FDPart3tmp63),
                                FusedMulAddSIMD(
                                    FDPart3tmp169, MulSIMD(FDPart3tmp25, FDPart3tmp28),
                                    FusedMulAddSIMD(
                                        FDPart3tmp194, FDPart3tmp98,
                                        FusedMulAddSIMD(
                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp188, FDPart3tmp60),
                                            FusedMulAddSIMD(
                                                FDPart3tmp192, FDPart3tmp74,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp194, FDPart3tmp45,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp189, FDPart3tmp86,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp192, FDPart3tmp51,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp186, FDPart3tmp26,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp189, FDPart3tmp60,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp185, FDPart3tmp32,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp185, FDPart3tmp76,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp184, FDPart3tmp60,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp184, FDPart3tmp86,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp178, FDPart3tmp58,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp181, FDPart3tmp54,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp174, FDPart3tmp63,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp176, FDPart3tmp45,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp164, FDPart3tmp51,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp166, FDPart3tmp26,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp101,
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3_NegativeOne_,
                                                                                                                            FDPart3_Rational_2_3),
                                                                                                                    MulSIMD(FDPart3tmp57, trK_dD1)),
                                                                                                                FusedMulSubSIMD(
                                                                                                                    FDPart3tmp162, FDPart3tmp62,
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3_Rational_1_2,
                                                                                                                                FDPart3tmp107),
                                                                                                                        MulSIMD(
                                                                                                                            FDPart3tmp171,
                                                                                                                            cf_dD0)))))))))))))))))))))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp200 = FusedMulAddSIMD(
            FDPart3tmp101, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp20, trK_dD0)),
            FusedMulAddSIMD(
                MulSIMD(FDPart3_Integer_2, FDPart3tmp186), MulSIMD(FDPart3tmp25, FDPart3tmp36),
                FusedMulAddSIMD(
                    FDPart3tmp107,
                    MulSIMD(
                        MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                        MulSIMD(cf_dD1, FusedMulAddSIMD(
                                            FDPart3_Integer_6, MulSIMD(FDPart3tmp60, FDPart3tmp69),
                                            FusedMulAddSIMD(
                                                FDPart3_Integer_6, MulSIMD(FDPart3tmp27, FDPart3tmp62),
                                                FusedMulAddSIMD(
                                                    FDPart3_Integer_6, MulSIMD(FDPart3tmp32, FDPart3tmp65),
                                                    FusedMulAddSIMD(
                                                        FDPart3_Integer_6, MulSIMD(FDPart3tmp22, FDPart3tmp51),
                                                        FusedMulAddSIMD(
                                                            FDPart3_Integer_6, MulSIMD(FDPart3tmp26, FDPart3tmp67),
                                                            FusedMulAddSIMD(
                                                                FDPart3_Integer_6, MulSIMD(FDPart3tmp65, FDPart3tmp76),
                                                                FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp69, FDPart3tmp86),
                                                                                FusedMulAddSIMD(FDPart3_Integer_6, FDPart3tmp84,
                                                                                                MulSIMD(FDPart3_Integer_6, FDPart3tmp89))))))))))),
                    FusedMulAddSIMD(
                        FDPart3tmp166, MulSIMD(FDPart3tmp25, FDPart3tmp36),
                        FusedMulAddSIMD(
                            FDPart3tmp178, MulSIMD(FDPart3tmp25, FDPart3tmp61),
                            FusedMulAddSIMD(
                                FDPart3_Integer_2, MulSIMD(FDPart3tmp185, FDPart3tmp37),
                                FusedMulAddSIMD(
                                    FDPart3tmp159, MulSIMD(FDPart3tmp20, FDPart3tmp36),
                                    FusedMulAddSIMD(
                                        FDPart3tmp194, FDPart3tmp74,
                                        FusedMulAddSIMD(
                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp184, FDPart3tmp32),
                                            FusedMulAddSIMD(
                                                FDPart3tmp193, FDPart3tmp86,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp194, FDPart3tmp51,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp192, FDPart3tmp78,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp193, FDPart3tmp60,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp191, FDPart3tmp36,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp192, FDPart3tmp21,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp189, FDPart3tmp76,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp190, FDPart3tmp26,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp188, FDPart3tmp76,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp189, FDPart3tmp32,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp181, FDPart3tmp53,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp188, FDPart3tmp32,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp174, FDPart3tmp60,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp176, FDPart3tmp51,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp164, FDPart3tmp21,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp169,
                                                                                                            MulSIMD(FDPart3tmp28, FDPart3tmp28),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp101,
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3_NegativeOne_,
                                                                                                                            FDPart3_Rational_2_3),
                                                                                                                    MulSIMD(FDPart3tmp25, trK_dD1)),
                                                                                                                FusedMulSubSIMD(
                                                                                                                    FDPart3tmp162, FDPart3tmp26,
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3_Rational_1_2,
                                                                                                                                FDPart3tmp107),
                                                                                                                        MulSIMD(
                                                                                                                            cf_dD0,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3_Integer_6,
                                                                                                                                MulSIMD(FDPart3tmp65,
                                                                                                                                        FDPart3tmp78),
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3_Integer_6,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp32,
                                                                                                                                        FDPart3tmp67),
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3_Integer_6,
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp51,
                                                                                                                                            FDPart3tmp69),
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3_Integer_6,
                                                                                                                                            MulSIMD(
                                                                                                                                                FDPart3tmp27,
                                                                                                                                                FDPart3tmp60),
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3_Integer_6,
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3tmp30,
                                                                                                                                                    FDPart3tmp37),
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3_Integer_6,
                                                                                                                                                    MulSIMD(
                                                                                                                                                        FDPart3tmp67,
                                                                                                                                                        FDPart3tmp76),
                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                        FDPart3_Integer_6,
                                                                                                                                                        MulSIMD(
                                                                                                                                                            FDPart3tmp69,
                                                                                                                                                            FDPart3tmp74),
                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                            FDPart3_Integer_6,
                                                                                                                                                            MulSIMD(
                                                                                                                                                                FDPart3tmp21,
                                                                                                                                                                FDPart3tmp65),
                                                                                                                                                            MulSIMD(
                                                                                                                                                                FDPart3_Integer_6,
                                                                                                                                                                MulSIMD(
                                                                                                                                                                    FDPart3tmp22,
                                                                                                                                                                    FDPart3tmp53)))))))))))))))))))))))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(
            MulSIMD(aDD01, f0_of_xx0),
            MulSIMD(
                MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                MulSIMD(f0_of_xx0__D0,
                        FusedMulAddSIMD(
                            FDPart3tmp65, FDPart3tmp74,
                            FusedMulAddSIMD(
                                FDPart3tmp51, FDPart3tmp65,
                                FusedMulAddSIMD(
                                    FDPart3tmp60, FDPart3tmp67,
                                    FusedMulAddSIMD(
                                        FDPart3tmp30, FDPart3tmp32,
                                        FusedMulAddSIMD(FDPart3tmp45, FDPart3tmp69,
                                                        FusedMulAddSIMD(FDPart3tmp67, FDPart3tmp86,
                                                                        FusedMulAddSIMD(FDPart3tmp69, FDPart3tmp98,
                                                                                        FusedMulAddSIMD(FDPart3tmp22, FDPart3tmp54,
                                                                                                        MulSIMD(FDPart3tmp27, FDPart3tmp63))))))))))),
            NegFusedMulAddSIMD(
                FDPart3tmp0,
                MulSIMD(aDD00,
                        FusedMulAddSIMD(FDPart3tmp27, FDPart3tmp45,
                                        FusedMulAddSIMD(FDPart3tmp35, FDPart3tmp54,
                                                        FusedMulAddSIMD(FDPart3tmp39, FDPart3tmp53,
                                                                        FusedMulAddSIMD(FDPart3tmp41, FDPart3tmp51,
                                                                                        FusedMulAddSIMD(FDPart3tmp21, FDPart3tmp30,
                                                                                                        MulSIMD(FDPart3tmp22, FDPart3tmp48))))))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp6, aDD12),
                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                            MulSIMD(f1_of_xx1,
                                    FusedMulAddSIMD(
                                        FDPart3tmp60, FDPart3tmp69,
                                        FusedMulAddSIMD(
                                            FDPart3tmp27, FDPart3tmp62,
                                            FusedMulAddSIMD(
                                                FDPart3tmp32, FDPart3tmp65,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp22, FDPart3tmp51,
                                                    FusedMulAddSIMD(FDPart3tmp26, FDPart3tmp67,
                                                                    FusedMulAddSIMD(FDPart3tmp65, FDPart3tmp76,
                                                                                    FusedMulAddSIMD(FDPart3tmp69, FDPart3tmp86,
                                                                                                    AddSIMD(FDPart3tmp84, FDPart3tmp89)))))))))),
                    FusedMulAddSIMD(
                        MulSIMD(cf, cf),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3tmp101, FDPart3tmp103),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp47, MulSIMD(cf_dD0, cf_dD0))),
                            FusedMulAddSIMD(
                                MulSIMD(FDPart3tmp101, FDPart3tmp103),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp57, MulSIMD(cf_dD1, cf_dD1))),
                                FusedMulAddSIMD(
                                    FDPart3tmp101,
                                    MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                                            MulSIMD(FDPart3tmp47,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp107,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp151, cf_dD1)),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp132, FusedMulSubSIMD(FDPart3tmp107, MulSIMD(cf_dD0, cf_dD0), cf_dDD00),
                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp107), MulSIMD(FDPart3tmp150, cf_dD0)))))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp101,
                                        MulSIMD(
                                            MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                                            MulSIMD(FDPart3tmp57,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp107,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp142, cf_dD1)),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp132, FusedMulSubSIMD(FDPart3tmp107, MulSIMD(cf_dD1, cf_dD1), cf_dDD11),
                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp107), MulSIMD(FDPart3tmp141, cf_dD0)))))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp101,
                                            MulSIMD(
                                                MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_),
                                                MulSIMD(
                                                    FDPart3tmp44,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp107,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp137, cf_dD1)),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp132, FusedMulSubSIMD(FDPart3tmp107, MulSIMD(cf_dD0, cf_dD1), cf_dDD01),
                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp107), MulSIMD(FDPart3tmp136, cf_dD0)))))),
                                            FusedMulAddSIMD(
                                                FDPart3tmp101,
                                                MulSIMD(
                                                    MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                                                    MulSIMD(FDPart3tmp28, FusedMulSubSIMD(FDPart3tmp107,
                                                                                          MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                                  MulSIMD(FDPart3tmp131, cf_dD1)),
                                                                                          MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp107),
                                                                                                  MulSIMD(FDPart3tmp130, cf_dD0))))),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp102, MulSIMD(FDPart3tmp44, RbarDD01),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp101,
                                                        MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_),
                                                                MulSIMD(FDPart3tmp25,
                                                                        FusedMulSubSIMD(FDPart3tmp107,
                                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                                MulSIMD(FDPart3tmp126, cf_dD0)),
                                                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp107),
                                                                                                MulSIMD(FDPart3tmp125, cf_dD1))))),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp102, MulSIMD(FDPart3tmp20, RbarDD02),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp102, MulSIMD(FDPart3tmp25, RbarDD12),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp101, MulSIMD(FDPart3tmp47, RbarDD00),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp101, MulSIMD(FDPart3tmp57, RbarDD11),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp44,
                                                                            MulSIMD(MulSIMD(FDPart3tmp101, FDPart3tmp103),
                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_),
                                                                                            MulSIMD(cf_dD0, cf_dD1))),
                                                                            FusedMulSubSIMD(
                                                                                FDPart3tmp101, MulSIMD(FDPart3tmp28, RbarDD22),
                                                                                MulSIMD(
                                                                                    MulSIMD(FDPart3_Integer_16, FDPart3tmp101),
                                                                                    MulSIMD(
                                                                                        FDPart3tmp20,
                                                                                        FusedMulSubSIMD(
                                                                                            FDPart3tmp107,
                                                                                            MulSIMD(
                                                                                                MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                                MulSIMD(FDPart3tmp119, cf_dD1)),
                                                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp107),
                                                                                                    MulSIMD(FDPart3tmp118, cf_dD0))))))))))))))))))),
                        SubSIMD(FusedMulAddSIMD(
                                    f0_of_xx0__D0,
                                    MulSIMD(MulSIMD(aDD02, f0_of_xx0),
                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                    MulSIMD(f1_of_xx1,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp65, FDPart3tmp78,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp32, FDPart3tmp67,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp51, FDPart3tmp69,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp27, FDPart3tmp60,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp30, FDPart3tmp37,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp67, FDPart3tmp76,
                                                                                    FusedMulAddSIMD(FDPart3tmp69, FDPart3tmp74,
                                                                                                    FusedMulAddSIMD(FDPart3tmp21, FDPart3tmp65,
                                                                                                                    MulSIMD(FDPart3tmp22,
                                                                                                                            FDPart3tmp53)))))))))))),
                                    FusedMulSubSIMD(
                                        FDPart3_Rational_2_3, MulSIMD(trK, trK),
                                        MulSIMD(FDPart3tmp27,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp27, FDPart3tmp58,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp35, FDPart3tmp63,
                                                        FusedMulAddSIMD(FDPart3tmp39, FDPart3tmp60,
                                                                        FusedMulAddSIMD(FDPart3tmp41, FDPart3tmp62,
                                                                                        FusedMulAddSIMD(FDPart3tmp22, FDPart3tmp45,
                                                                                                        MulSIMD(FDPart3tmp26, FDPart3tmp30))))))))),
                                MulSIMD(FDPart3tmp30,
                                        FusedMulAddSIMD(
                                            FDPart3tmp32, FDPart3tmp35,
                                            FusedMulAddSIMD(
                                                FDPart3tmp37, FDPart3tmp39,
                                                FusedMulAddSIMD(FDPart3tmp15, MulSIMD(FDPart3tmp30, MulSIMD(FDPart3tmp28, FDPart3tmp28)),
                                                                FusedMulAddSIMD(FDPart3tmp25, MulSIMD(FDPart3tmp36, FDPart3tmp41),
                                                                                FusedMulAddSIMD(FDPart3tmp21, FDPart3tmp22,
                                                                                                MulSIMD(FDPart3tmp26, FDPart3tmp27))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(
            FDPart3tmp152, MulSIMD(FDPart3tmp7, MulSIMD(FDPart3tmp199, FDPart3tmp199)),
            FusedMulAddSIMD(
                FDPart3tmp199, MulSIMD(MulSIMD(FDPart3tmp152, FDPart3tmp195), MulSIMD(FDPart3tmp34, hDD01)),
                FusedMulAddSIMD(
                    FDPart3tmp200, MulSIMD(MulSIMD(FDPart3tmp152, FDPart3tmp195), MulSIMD(FDPart3tmp38, hDD02)),
                    FusedMulAddSIMD(FDPart3tmp200, MulSIMD(MulSIMD(FDPart3tmp152, FDPart3tmp199), MulSIMD(FDPart3tmp40, hDD12)),
                                    FusedMulAddSIMD(FDPart3tmp11, MulSIMD(FDPart3tmp152, MulSIMD(FDPart3tmp200, FDPart3tmp200)),
                                                    MulSIMD(FDPart3tmp152, MulSIMD(FDPart3tmp4, MulSIMD(FDPart3tmp195, FDPart3tmp195))))))));

        WriteSIMD(&diagnostic_output_gfs[IDX4(HGF, i0, i1, i2)], __RHS_exp_0);
        WriteSIMD(&diagnostic_output_gfs[IDX4(MSQUAREDGF, i0, i1, i2)], __RHS_exp_1);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += simd_width)
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
