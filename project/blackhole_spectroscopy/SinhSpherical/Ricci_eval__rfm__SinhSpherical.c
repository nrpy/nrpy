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
 * Set Ricci tensor.
 */
void Ricci_eval__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                    const rfm_struct *restrict rfmstruct, const REAL *restrict in_gfs, REAL *restrict auxevol_gfs) {
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
        const REAL_SIMD_ARRAY hDD00_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 4, i2)]);
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
        const REAL_SIMD_ARRAY hDD_dDD0000 = SIMD_fd_function_dDD00_fdorder8(hDD00, hDD00_i0m1, hDD00_i0m2, hDD00_i0m3, hDD00_i0m4, hDD00_i0p1,
                                                                            hDD00_i0p2, hDD00_i0p3, hDD00_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD0001 = SIMD_fd_function_dDD01_fdorder8(
            hDD00_i0m1_i1m1, hDD00_i0m1_i1m2, hDD00_i0m1_i1m3, hDD00_i0m1_i1m4, hDD00_i0m1_i1p1, hDD00_i0m1_i1p2, hDD00_i0m1_i1p3, hDD00_i0m1_i1p4,
            hDD00_i0m2_i1m1, hDD00_i0m2_i1m2, hDD00_i0m2_i1m3, hDD00_i0m2_i1m4, hDD00_i0m2_i1p1, hDD00_i0m2_i1p2, hDD00_i0m2_i1p3, hDD00_i0m2_i1p4,
            hDD00_i0m3_i1m1, hDD00_i0m3_i1m2, hDD00_i0m3_i1m3, hDD00_i0m3_i1m4, hDD00_i0m3_i1p1, hDD00_i0m3_i1p2, hDD00_i0m3_i1p3, hDD00_i0m3_i1p4,
            hDD00_i0m4_i1m1, hDD00_i0m4_i1m2, hDD00_i0m4_i1m3, hDD00_i0m4_i1m4, hDD00_i0m4_i1p1, hDD00_i0m4_i1p2, hDD00_i0m4_i1p3, hDD00_i0m4_i1p4,
            hDD00_i0p1_i1m1, hDD00_i0p1_i1m2, hDD00_i0p1_i1m3, hDD00_i0p1_i1m4, hDD00_i0p1_i1p1, hDD00_i0p1_i1p2, hDD00_i0p1_i1p3, hDD00_i0p1_i1p4,
            hDD00_i0p2_i1m1, hDD00_i0p2_i1m2, hDD00_i0p2_i1m3, hDD00_i0p2_i1m4, hDD00_i0p2_i1p1, hDD00_i0p2_i1p2, hDD00_i0p2_i1p3, hDD00_i0p2_i1p4,
            hDD00_i0p3_i1m1, hDD00_i0p3_i1m2, hDD00_i0p3_i1m3, hDD00_i0p3_i1m4, hDD00_i0p3_i1p1, hDD00_i0p3_i1p2, hDD00_i0p3_i1p3, hDD00_i0p3_i1p4,
            hDD00_i0p4_i1m1, hDD00_i0p4_i1m2, hDD00_i0p4_i1m3, hDD00_i0p4_i1m4, hDD00_i0p4_i1p1, hDD00_i0p4_i1p2, hDD00_i0p4_i1p3, hDD00_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0011 = SIMD_fd_function_dDD11_fdorder8(hDD00, hDD00_i1m1, hDD00_i1m2, hDD00_i1m3, hDD00_i1m4, hDD00_i1p1,
                                                                            hDD00_i1p2, hDD00_i1p3, hDD00_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0100 = SIMD_fd_function_dDD00_fdorder8(hDD01, hDD01_i0m1, hDD01_i0m2, hDD01_i0m3, hDD01_i0m4, hDD01_i0p1,
                                                                            hDD01_i0p2, hDD01_i0p3, hDD01_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD0101 = SIMD_fd_function_dDD01_fdorder8(
            hDD01_i0m1_i1m1, hDD01_i0m1_i1m2, hDD01_i0m1_i1m3, hDD01_i0m1_i1m4, hDD01_i0m1_i1p1, hDD01_i0m1_i1p2, hDD01_i0m1_i1p3, hDD01_i0m1_i1p4,
            hDD01_i0m2_i1m1, hDD01_i0m2_i1m2, hDD01_i0m2_i1m3, hDD01_i0m2_i1m4, hDD01_i0m2_i1p1, hDD01_i0m2_i1p2, hDD01_i0m2_i1p3, hDD01_i0m2_i1p4,
            hDD01_i0m3_i1m1, hDD01_i0m3_i1m2, hDD01_i0m3_i1m3, hDD01_i0m3_i1m4, hDD01_i0m3_i1p1, hDD01_i0m3_i1p2, hDD01_i0m3_i1p3, hDD01_i0m3_i1p4,
            hDD01_i0m4_i1m1, hDD01_i0m4_i1m2, hDD01_i0m4_i1m3, hDD01_i0m4_i1m4, hDD01_i0m4_i1p1, hDD01_i0m4_i1p2, hDD01_i0m4_i1p3, hDD01_i0m4_i1p4,
            hDD01_i0p1_i1m1, hDD01_i0p1_i1m2, hDD01_i0p1_i1m3, hDD01_i0p1_i1m4, hDD01_i0p1_i1p1, hDD01_i0p1_i1p2, hDD01_i0p1_i1p3, hDD01_i0p1_i1p4,
            hDD01_i0p2_i1m1, hDD01_i0p2_i1m2, hDD01_i0p2_i1m3, hDD01_i0p2_i1m4, hDD01_i0p2_i1p1, hDD01_i0p2_i1p2, hDD01_i0p2_i1p3, hDD01_i0p2_i1p4,
            hDD01_i0p3_i1m1, hDD01_i0p3_i1m2, hDD01_i0p3_i1m3, hDD01_i0p3_i1m4, hDD01_i0p3_i1p1, hDD01_i0p3_i1p2, hDD01_i0p3_i1p3, hDD01_i0p3_i1p4,
            hDD01_i0p4_i1m1, hDD01_i0p4_i1m2, hDD01_i0p4_i1m3, hDD01_i0p4_i1m4, hDD01_i0p4_i1p1, hDD01_i0p4_i1p2, hDD01_i0p4_i1p3, hDD01_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0111 = SIMD_fd_function_dDD11_fdorder8(hDD01, hDD01_i1m1, hDD01_i1m2, hDD01_i1m3, hDD01_i1m4, hDD01_i1p1,
                                                                            hDD01_i1p2, hDD01_i1p3, hDD01_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0200 = SIMD_fd_function_dDD00_fdorder8(hDD02, hDD02_i0m1, hDD02_i0m2, hDD02_i0m3, hDD02_i0m4, hDD02_i0p1,
                                                                            hDD02_i0p2, hDD02_i0p3, hDD02_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD0201 = SIMD_fd_function_dDD01_fdorder8(
            hDD02_i0m1_i1m1, hDD02_i0m1_i1m2, hDD02_i0m1_i1m3, hDD02_i0m1_i1m4, hDD02_i0m1_i1p1, hDD02_i0m1_i1p2, hDD02_i0m1_i1p3, hDD02_i0m1_i1p4,
            hDD02_i0m2_i1m1, hDD02_i0m2_i1m2, hDD02_i0m2_i1m3, hDD02_i0m2_i1m4, hDD02_i0m2_i1p1, hDD02_i0m2_i1p2, hDD02_i0m2_i1p3, hDD02_i0m2_i1p4,
            hDD02_i0m3_i1m1, hDD02_i0m3_i1m2, hDD02_i0m3_i1m3, hDD02_i0m3_i1m4, hDD02_i0m3_i1p1, hDD02_i0m3_i1p2, hDD02_i0m3_i1p3, hDD02_i0m3_i1p4,
            hDD02_i0m4_i1m1, hDD02_i0m4_i1m2, hDD02_i0m4_i1m3, hDD02_i0m4_i1m4, hDD02_i0m4_i1p1, hDD02_i0m4_i1p2, hDD02_i0m4_i1p3, hDD02_i0m4_i1p4,
            hDD02_i0p1_i1m1, hDD02_i0p1_i1m2, hDD02_i0p1_i1m3, hDD02_i0p1_i1m4, hDD02_i0p1_i1p1, hDD02_i0p1_i1p2, hDD02_i0p1_i1p3, hDD02_i0p1_i1p4,
            hDD02_i0p2_i1m1, hDD02_i0p2_i1m2, hDD02_i0p2_i1m3, hDD02_i0p2_i1m4, hDD02_i0p2_i1p1, hDD02_i0p2_i1p2, hDD02_i0p2_i1p3, hDD02_i0p2_i1p4,
            hDD02_i0p3_i1m1, hDD02_i0p3_i1m2, hDD02_i0p3_i1m3, hDD02_i0p3_i1m4, hDD02_i0p3_i1p1, hDD02_i0p3_i1p2, hDD02_i0p3_i1p3, hDD02_i0p3_i1p4,
            hDD02_i0p4_i1m1, hDD02_i0p4_i1m2, hDD02_i0p4_i1m3, hDD02_i0p4_i1m4, hDD02_i0p4_i1p1, hDD02_i0p4_i1p2, hDD02_i0p4_i1p3, hDD02_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0211 = SIMD_fd_function_dDD11_fdorder8(hDD02, hDD02_i1m1, hDD02_i1m2, hDD02_i1m3, hDD02_i1m4, hDD02_i1p1,
                                                                            hDD02_i1p2, hDD02_i1p3, hDD02_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD1100 = SIMD_fd_function_dDD00_fdorder8(hDD11, hDD11_i0m1, hDD11_i0m2, hDD11_i0m3, hDD11_i0m4, hDD11_i0p1,
                                                                            hDD11_i0p2, hDD11_i0p3, hDD11_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD1101 = SIMD_fd_function_dDD01_fdorder8(
            hDD11_i0m1_i1m1, hDD11_i0m1_i1m2, hDD11_i0m1_i1m3, hDD11_i0m1_i1m4, hDD11_i0m1_i1p1, hDD11_i0m1_i1p2, hDD11_i0m1_i1p3, hDD11_i0m1_i1p4,
            hDD11_i0m2_i1m1, hDD11_i0m2_i1m2, hDD11_i0m2_i1m3, hDD11_i0m2_i1m4, hDD11_i0m2_i1p1, hDD11_i0m2_i1p2, hDD11_i0m2_i1p3, hDD11_i0m2_i1p4,
            hDD11_i0m3_i1m1, hDD11_i0m3_i1m2, hDD11_i0m3_i1m3, hDD11_i0m3_i1m4, hDD11_i0m3_i1p1, hDD11_i0m3_i1p2, hDD11_i0m3_i1p3, hDD11_i0m3_i1p4,
            hDD11_i0m4_i1m1, hDD11_i0m4_i1m2, hDD11_i0m4_i1m3, hDD11_i0m4_i1m4, hDD11_i0m4_i1p1, hDD11_i0m4_i1p2, hDD11_i0m4_i1p3, hDD11_i0m4_i1p4,
            hDD11_i0p1_i1m1, hDD11_i0p1_i1m2, hDD11_i0p1_i1m3, hDD11_i0p1_i1m4, hDD11_i0p1_i1p1, hDD11_i0p1_i1p2, hDD11_i0p1_i1p3, hDD11_i0p1_i1p4,
            hDD11_i0p2_i1m1, hDD11_i0p2_i1m2, hDD11_i0p2_i1m3, hDD11_i0p2_i1m4, hDD11_i0p2_i1p1, hDD11_i0p2_i1p2, hDD11_i0p2_i1p3, hDD11_i0p2_i1p4,
            hDD11_i0p3_i1m1, hDD11_i0p3_i1m2, hDD11_i0p3_i1m3, hDD11_i0p3_i1m4, hDD11_i0p3_i1p1, hDD11_i0p3_i1p2, hDD11_i0p3_i1p3, hDD11_i0p3_i1p4,
            hDD11_i0p4_i1m1, hDD11_i0p4_i1m2, hDD11_i0p4_i1m3, hDD11_i0p4_i1m4, hDD11_i0p4_i1p1, hDD11_i0p4_i1p2, hDD11_i0p4_i1p3, hDD11_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD1111 = SIMD_fd_function_dDD11_fdorder8(hDD11, hDD11_i1m1, hDD11_i1m2, hDD11_i1m3, hDD11_i1m4, hDD11_i1p1,
                                                                            hDD11_i1p2, hDD11_i1p3, hDD11_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD1200 = SIMD_fd_function_dDD00_fdorder8(hDD12, hDD12_i0m1, hDD12_i0m2, hDD12_i0m3, hDD12_i0m4, hDD12_i0p1,
                                                                            hDD12_i0p2, hDD12_i0p3, hDD12_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD1201 = SIMD_fd_function_dDD01_fdorder8(
            hDD12_i0m1_i1m1, hDD12_i0m1_i1m2, hDD12_i0m1_i1m3, hDD12_i0m1_i1m4, hDD12_i0m1_i1p1, hDD12_i0m1_i1p2, hDD12_i0m1_i1p3, hDD12_i0m1_i1p4,
            hDD12_i0m2_i1m1, hDD12_i0m2_i1m2, hDD12_i0m2_i1m3, hDD12_i0m2_i1m4, hDD12_i0m2_i1p1, hDD12_i0m2_i1p2, hDD12_i0m2_i1p3, hDD12_i0m2_i1p4,
            hDD12_i0m3_i1m1, hDD12_i0m3_i1m2, hDD12_i0m3_i1m3, hDD12_i0m3_i1m4, hDD12_i0m3_i1p1, hDD12_i0m3_i1p2, hDD12_i0m3_i1p3, hDD12_i0m3_i1p4,
            hDD12_i0m4_i1m1, hDD12_i0m4_i1m2, hDD12_i0m4_i1m3, hDD12_i0m4_i1m4, hDD12_i0m4_i1p1, hDD12_i0m4_i1p2, hDD12_i0m4_i1p3, hDD12_i0m4_i1p4,
            hDD12_i0p1_i1m1, hDD12_i0p1_i1m2, hDD12_i0p1_i1m3, hDD12_i0p1_i1m4, hDD12_i0p1_i1p1, hDD12_i0p1_i1p2, hDD12_i0p1_i1p3, hDD12_i0p1_i1p4,
            hDD12_i0p2_i1m1, hDD12_i0p2_i1m2, hDD12_i0p2_i1m3, hDD12_i0p2_i1m4, hDD12_i0p2_i1p1, hDD12_i0p2_i1p2, hDD12_i0p2_i1p3, hDD12_i0p2_i1p4,
            hDD12_i0p3_i1m1, hDD12_i0p3_i1m2, hDD12_i0p3_i1m3, hDD12_i0p3_i1m4, hDD12_i0p3_i1p1, hDD12_i0p3_i1p2, hDD12_i0p3_i1p3, hDD12_i0p3_i1p4,
            hDD12_i0p4_i1m1, hDD12_i0p4_i1m2, hDD12_i0p4_i1m3, hDD12_i0p4_i1m4, hDD12_i0p4_i1p1, hDD12_i0p4_i1p2, hDD12_i0p4_i1p3, hDD12_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD1211 = SIMD_fd_function_dDD11_fdorder8(hDD12, hDD12_i1m1, hDD12_i1m2, hDD12_i1m3, hDD12_i1m4, hDD12_i1p1,
                                                                            hDD12_i1p2, hDD12_i1p3, hDD12_i1p4, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD2200 = SIMD_fd_function_dDD00_fdorder8(hDD22, hDD22_i0m1, hDD22_i0m2, hDD22_i0m3, hDD22_i0m4, hDD22_i0p1,
                                                                            hDD22_i0p2, hDD22_i0p3, hDD22_i0p4, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD2201 = SIMD_fd_function_dDD01_fdorder8(
            hDD22_i0m1_i1m1, hDD22_i0m1_i1m2, hDD22_i0m1_i1m3, hDD22_i0m1_i1m4, hDD22_i0m1_i1p1, hDD22_i0m1_i1p2, hDD22_i0m1_i1p3, hDD22_i0m1_i1p4,
            hDD22_i0m2_i1m1, hDD22_i0m2_i1m2, hDD22_i0m2_i1m3, hDD22_i0m2_i1m4, hDD22_i0m2_i1p1, hDD22_i0m2_i1p2, hDD22_i0m2_i1p3, hDD22_i0m2_i1p4,
            hDD22_i0m3_i1m1, hDD22_i0m3_i1m2, hDD22_i0m3_i1m3, hDD22_i0m3_i1m4, hDD22_i0m3_i1p1, hDD22_i0m3_i1p2, hDD22_i0m3_i1p3, hDD22_i0m3_i1p4,
            hDD22_i0m4_i1m1, hDD22_i0m4_i1m2, hDD22_i0m4_i1m3, hDD22_i0m4_i1m4, hDD22_i0m4_i1p1, hDD22_i0m4_i1p2, hDD22_i0m4_i1p3, hDD22_i0m4_i1p4,
            hDD22_i0p1_i1m1, hDD22_i0p1_i1m2, hDD22_i0p1_i1m3, hDD22_i0p1_i1m4, hDD22_i0p1_i1p1, hDD22_i0p1_i1p2, hDD22_i0p1_i1p3, hDD22_i0p1_i1p4,
            hDD22_i0p2_i1m1, hDD22_i0p2_i1m2, hDD22_i0p2_i1m3, hDD22_i0p2_i1m4, hDD22_i0p2_i1p1, hDD22_i0p2_i1p2, hDD22_i0p2_i1p3, hDD22_i0p2_i1p4,
            hDD22_i0p3_i1m1, hDD22_i0p3_i1m2, hDD22_i0p3_i1m3, hDD22_i0p3_i1m4, hDD22_i0p3_i1p1, hDD22_i0p3_i1p2, hDD22_i0p3_i1p3, hDD22_i0p3_i1p4,
            hDD22_i0p4_i1m1, hDD22_i0p4_i1m2, hDD22_i0p4_i1m3, hDD22_i0p4_i1m4, hDD22_i0p4_i1p1, hDD22_i0p4_i1p2, hDD22_i0p4_i1p3, hDD22_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD2211 = SIMD_fd_function_dDD11_fdorder8(hDD22, hDD22_i1m1, hDD22_i1m2, hDD22_i1m3, hDD22_i1m4, hDD22_i1p1,
                                                                            hDD22_i1p2, hDD22_i1p3, hDD22_i1p4, invdxx1);
        const REAL_SIMD_ARRAY lambdaU_dD00 = SIMD_fd_function_dD0_fdorder8(lambdaU0_i0m1, lambdaU0_i0m2, lambdaU0_i0m3, lambdaU0_i0m4, lambdaU0_i0p1,
                                                                           lambdaU0_i0p2, lambdaU0_i0p3, lambdaU0_i0p4, invdxx0);
        const REAL_SIMD_ARRAY lambdaU_dD01 = SIMD_fd_function_dD1_fdorder8(lambdaU0_i1m1, lambdaU0_i1m2, lambdaU0_i1m3, lambdaU0_i1m4, lambdaU0_i1p1,
                                                                           lambdaU0_i1p2, lambdaU0_i1p3, lambdaU0_i1p4, invdxx1);
        const REAL_SIMD_ARRAY lambdaU_dD10 = SIMD_fd_function_dD0_fdorder8(lambdaU1_i0m1, lambdaU1_i0m2, lambdaU1_i0m3, lambdaU1_i0m4, lambdaU1_i0p1,
                                                                           lambdaU1_i0p2, lambdaU1_i0p3, lambdaU1_i0p4, invdxx0);
        const REAL_SIMD_ARRAY lambdaU_dD11 = SIMD_fd_function_dD1_fdorder8(lambdaU1_i1m1, lambdaU1_i1m2, lambdaU1_i1m3, lambdaU1_i1m4, lambdaU1_i1p1,
                                                                           lambdaU1_i1p2, lambdaU1_i1p3, lambdaU1_i1p4, invdxx1);
        const REAL_SIMD_ARRAY lambdaU_dD20 = SIMD_fd_function_dD0_fdorder8(lambdaU2_i0m1, lambdaU2_i0m2, lambdaU2_i0m3, lambdaU2_i0m4, lambdaU2_i0p1,
                                                                           lambdaU2_i0p2, lambdaU2_i0p3, lambdaU2_i0p4, invdxx0);
        const REAL_SIMD_ARRAY lambdaU_dD21 = SIMD_fd_function_dD1_fdorder8(lambdaU2_i1m1, lambdaU2_i1m2, lambdaU2_i1m3, lambdaU2_i1m4, lambdaU2_i1p1,
                                                                           lambdaU2_i1p2, lambdaU2_i1p3, lambdaU2_i1p4, invdxx1);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        const double dblFDPart3_Integer_1 = 1.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(dblFDPart3_Integer_1);

        const double dblFDPart3_Integer_2 = 2.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(dblFDPart3_Integer_2);

        const double dblFDPart3_Integer_3 = 3.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_3 = ConstSIMD(dblFDPart3_Integer_3);

        const double dblFDPart3_Integer_4 = 4.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_4 = ConstSIMD(dblFDPart3_Integer_4);

        const double dblFDPart3_Integer_6 = 6.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_6 = ConstSIMD(dblFDPart3_Integer_6);

        const double dblFDPart3_NegativeOne_ = -1.0;
        const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(dblFDPart3_NegativeOne_);

        const double dblFDPart3_Rational_1_2 = 1.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_2 = ConstSIMD(dblFDPart3_Rational_1_2);

        const REAL_SIMD_ARRAY FDPart3tmp0 = MulSIMD(f0_of_xx0__D0, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp1 = MulSIMD(f0_of_xx0__D0, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp2 = MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp5 = DivSIMD(FDPart3_Integer_1, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp7 = MulSIMD(f0_of_xx0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp17 = DivSIMD(FDPart3_Integer_1, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp19 = MulSIMD(f0_of_xx0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp23 = MulSIMD(f0_of_xx0, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp27 = MulSIMD(MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0), f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp28 = MulSIMD(f1_of_xx1, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp43 = MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp50 = MulSIMD(f0_of_xx0, f0_of_xx0__DD00);
        const REAL_SIMD_ARRAY FDPart3tmp73 = MulSIMD(f1_of_xx1, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp94 = MulSIMD(FDPart3_Integer_2, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp170 = DivSIMD(FDPart3_Integer_1, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp259 = MulSIMD(f1_of_xx1__D1, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp279 = MulSIMD(FDPart3_Integer_3, MulSIMD(f0_of_xx0__D0, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp331 = MulSIMD(MulSIMD(f1_of_xx1, f1_of_xx1), f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp336 = MulSIMD(f1_of_xx1, f1_of_xx1__DD11);
        const REAL_SIMD_ARRAY FDPart3tmp3 = MulSIMD(FDPart3tmp2, hDD00);
        const REAL_SIMD_ARRAY FDPart3tmp8 = MulSIMD(FDPart3tmp2, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp10 = MulSIMD(FDPart3tmp7, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp20 = MulSIMD(FDPart3tmp19, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp24 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp23), MulSIMD(f0_of_xx0__D0, f1_of_xx1));
        const REAL_SIMD_ARRAY FDPart3tmp33 = MulSIMD(FDPart3tmp7, hDD11);
        const REAL_SIMD_ARRAY FDPart3tmp35 = MulSIMD(FDPart3tmp28, FDPart3tmp7);
        const REAL_SIMD_ARRAY FDPart3tmp45 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(f0_of_xx0, f0_of_xx0__D0));
        const REAL_SIMD_ARRAY FDPart3tmp49 =
            MulSIMD(f1_of_xx1,
                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp17), MulSIMD(hDD02, MulSIMD(MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0))));
        const REAL_SIMD_ARRAY FDPart3tmp51 = MulSIMD(FDPart3tmp50, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp56 = MulSIMD(FDPart3tmp1, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp63 = MulSIMD(f0_of_xx0, MulSIMD(f0_of_xx0__D0, hDD_dD000));
        const REAL_SIMD_ARRAY FDPart3tmp71 = MulSIMD(FDPart3tmp2, hDD_dD001);
        const REAL_SIMD_ARRAY FDPart3tmp74 = MulSIMD(FDPart3tmp0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp81 = AddSIMD(FDPart3tmp2, FDPart3tmp50);
        const REAL_SIMD_ARRAY FDPart3tmp85 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp17, f0_of_xx0__D0));
        const REAL_SIMD_ARRAY FDPart3tmp88 = DivSIMD(FDPart3tmp2, FDPart3tmp7);
        const REAL_SIMD_ARRAY FDPart3tmp97 = DivSIMD(FDPart3_Integer_1, FDPart3tmp2);
        const REAL_SIMD_ARRAY FDPart3tmp100 = MulSIMD(FDPart3tmp94, f0_of_xx0__DD00);
        const REAL_SIMD_ARRAY FDPart3tmp105 = MulSIMD(FDPart3tmp94, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp127 = MulSIMD(FDPart3tmp7, hDD_dD110);
        const REAL_SIMD_ARRAY FDPart3tmp163 = MulSIMD(FDPart3_Integer_2, FDPart3tmp7);
        const REAL_SIMD_ARRAY FDPart3tmp187 = MulSIMD(FDPart3tmp5, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp213 = MulSIMD(FDPart3tmp7, hDD_dD111);
        const REAL_SIMD_ARRAY FDPart3tmp258 = DivSIMD(FDPart3_Integer_1, FDPart3tmp28);
        const REAL_SIMD_ARRAY FDPart3tmp271 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp7, FDPart3tmp73));
        const REAL_SIMD_ARRAY FDPart3tmp366 = MulSIMD(FDPart3tmp43, FDPart3tmp5);
        const REAL_SIMD_ARRAY FDPart3tmp414 = MulSIMD(FDPart3_Integer_2, FDPart3tmp73);
        const REAL_SIMD_ARRAY FDPart3tmp4 = AddSIMD(FDPart3tmp2, FDPart3tmp3);
        const REAL_SIMD_ARRAY FDPart3tmp11 = MulSIMD(FDPart3tmp10, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp16 = MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp2), MulSIMD(f1_of_xx1__D1, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp30 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp2, FDPart3tmp28));
        const REAL_SIMD_ARRAY FDPart3tmp32 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp27), MulSIMD(FDPart3tmp28, MulSIMD(hDD12, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp34 = AddSIMD(FDPart3tmp33, FDPart3tmp7);
        const REAL_SIMD_ARRAY FDPart3tmp37 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp2), MulSIMD(FDPart3tmp35, MulSIMD(hDD02, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp38 = MulSIMD(FDPart3tmp35, hDD22);
        const REAL_SIMD_ARRAY FDPart3tmp40 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp2), MulSIMD(FDPart3tmp7, MulSIMD(hDD01, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp46 = MulSIMD(FDPart3tmp45, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp52 = AddSIMD(FDPart3tmp51, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp59 = MulSIMD(FDPart3tmp45, MulSIMD(f1_of_xx1__D1, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp66 =
            FusedMulAddSIMD(FDPart3tmp19, hDD_dD011, FusedMulAddSIMD(FDPart3tmp45, hDD11, MulSIMD(FDPart3tmp19, hDD00)));
        const REAL_SIMD_ARRAY FDPart3tmp72 = NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp2, hDD01), FDPart3tmp71);
        const REAL_SIMD_ARRAY FDPart3tmp75 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp19, FDPart3tmp28));
        const REAL_SIMD_ARRAY FDPart3tmp78 = MulSIMD(FDPart3tmp45, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp89 = NegFusedMulAddSIMD(FDPart3tmp17, f0_of_xx0__DD00, FDPart3tmp88);
        const REAL_SIMD_ARRAY FDPart3tmp90 = MulSIMD(FDPart3tmp56, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp91 = FusedMulSubSIMD(FDPart3tmp17, f0_of_xx0__DD00, FDPart3tmp88);
        const REAL_SIMD_ARRAY FDPart3tmp98 = MulSIMD(FDPart3tmp97, MulSIMD(f0_of_xx0__DD00, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp101 = FusedMulAddSIMD(FDPart3tmp100, hDD00, MulSIMD(FDPart3tmp2, hDD_dD000));
        const REAL_SIMD_ARRAY FDPart3tmp106 = MulSIMD(FDPart3tmp105, FDPart3tmp28);
        const REAL_SIMD_ARRAY FDPart3tmp107 = MulSIMD(FDPart3tmp35, hDD_dD220);
        const REAL_SIMD_ARRAY FDPart3tmp111 = MulSIMD(FDPart3tmp10, hDD_dD120);
        const REAL_SIMD_ARRAY FDPart3tmp112 = MulSIMD(FDPart3tmp105, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp128 = FusedMulAddSIMD(FDPart3tmp105, hDD11, FDPart3tmp127);
        const REAL_SIMD_ARRAY FDPart3tmp137 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp5, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp164 = MulSIMD(FDPart3tmp163, FDPart3tmp73);
        const REAL_SIMD_ARRAY FDPart3tmp172 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp170, f1_of_xx1__D1));
        const REAL_SIMD_ARRAY FDPart3tmp188 = MulSIMD(FDPart3tmp187, FDPart3tmp28);
        const REAL_SIMD_ARRAY FDPart3tmp214 = MulSIMD(FDPart3tmp10, hDD_dD121);
        const REAL_SIMD_ARRAY FDPart3tmp231 = MulSIMD(FDPart3tmp105, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp243 = FusedMulAddSIMD(FDPart3tmp17, lambdaU0, MulSIMD(FDPart3tmp17, lambdaU_dD11));
        const REAL_SIMD_ARRAY FDPart3tmp245 = FusedMulSubSIMD(FDPart3tmp5, lambdaU_dD01, MulSIMD(FDPart3tmp5, lambdaU1));
        const REAL_SIMD_ARRAY FDPart3tmp251 = FusedMulAddSIMD(FDPart3tmp19, hDD_dDD0101, MulSIMD(FDPart3tmp81, hDD_dD011));
        const REAL_SIMD_ARRAY FDPart3tmp254 = MulSIMD(FDPart3tmp105, MulSIMD(f1_of_xx1__D1, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp261 = FusedMulSubSIMD(FDPart3tmp170, f1_of_xx1__DD11, MulSIMD(FDPart3tmp258, FDPart3tmp259));
        const REAL_SIMD_ARRAY FDPart3tmp272 =
            FusedMulAddSIMD(FDPart3tmp33, FDPart3tmp73, FusedMulAddSIMD(FDPart3tmp35, hDD01, MulSIMD(FDPart3tmp271, hDD22)));
        const REAL_SIMD_ARRAY FDPart3tmp276 = MulSIMD(FDPart3tmp50, FDPart3tmp97);
        const REAL_SIMD_ARRAY FDPart3tmp326 =
            FusedMulAddSIMD(FDPart3tmp17, lambdaU0, MulSIMD(MulSIMD(FDPart3tmp17, FDPart3tmp170), MulSIMD(f1_of_xx1__D1, lambdaU1)));
        const REAL_SIMD_ARRAY FDPart3tmp330 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp19), MulSIMD(FDPart3tmp28, hDD_dD221));
        const REAL_SIMD_ARRAY FDPart3tmp337 = MulSIMD(FDPart3_NegativeOne_, AddSIMD(FDPart3tmp336, FDPart3tmp259));
        const REAL_SIMD_ARRAY FDPart3tmp368 = MulSIMD(FDPart3_Integer_2, FDPart3tmp187);
        const REAL_SIMD_ARRAY FDPart3tmp374 = FusedMulSubSIMD(FDPart3tmp258, FDPart3tmp259, MulSIMD(FDPart3tmp170, f1_of_xx1__DD11));
        const REAL_SIMD_ARRAY FDPart3tmp411 =
            FusedMulAddSIMD(FDPart3tmp35, hDD_dDD2201, MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3tmp19), MulSIMD(FDPart3tmp73, hDD22)));
        const REAL_SIMD_ARRAY FDPart3tmp25 = FusedMulAddSIMD(FDPart3tmp20, hDD_dD021, FDPart3tmp24);
        const REAL_SIMD_ARRAY FDPart3tmp39 = AddSIMD(FDPart3tmp35, FDPart3tmp38);
        const REAL_SIMD_ARRAY FDPart3tmp47 = MulSIMD(FDPart3tmp46, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp57 = FusedMulAddSIMD(FDPart3tmp20, hDD_dD021, MulSIMD(FDPart3tmp56, f1_of_xx1__D1));
        const REAL_SIMD_ARRAY FDPart3tmp60 = AddSIMD(FDPart3tmp24, FDPart3tmp59);
        const REAL_SIMD_ARRAY FDPart3tmp67 =
            MulSIMD(FDPart3tmp17, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp66, f0_of_xx0__D0)));
        const REAL_SIMD_ARRAY FDPart3tmp69 = FusedMulAddSIMD(FDPart3tmp34, FDPart3tmp4, FDPart3tmp40);
        const REAL_SIMD_ARRAY FDPart3tmp76 =
            FusedMulAddSIMD(FDPart3tmp75, hDD22, FusedMulAddSIMD(FDPart3tmp19, MulSIMD(FDPart3tmp28, hDD00), MulSIMD(FDPart3tmp73, FDPart3tmp74)));
        const REAL_SIMD_ARRAY FDPart3tmp83 = FusedMulAddSIMD(FDPart3tmp19, hDD_dD010, MulSIMD(FDPart3tmp81, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp86 = FusedMulAddSIMD(FDPart3tmp2, hDD_dDD0001, MulSIMD(FDPart3tmp72, FDPart3tmp85));
        const REAL_SIMD_ARRAY FDPart3tmp99 = NegFusedMulAddSIMD(FDPart3tmp5, f0_of_xx0__DDD000, FDPart3tmp98);
        const REAL_SIMD_ARRAY FDPart3tmp108 = FusedMulAddSIMD(FDPart3tmp106, hDD22, FDPart3tmp107);
        const REAL_SIMD_ARRAY FDPart3tmp113 = FusedMulAddSIMD(FDPart3tmp112, hDD12, FDPart3tmp111);
        const REAL_SIMD_ARRAY FDPart3tmp129 = AddSIMD(FDPart3tmp105, FDPart3tmp128);
        const REAL_SIMD_ARRAY FDPart3tmp139 = AddSIMD(FDPart3tmp100, FDPart3tmp101);
        const REAL_SIMD_ARRAY FDPart3tmp140 =
            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp81, hDD01), FusedMulSubSIMD(FDPart3tmp105, hDD_dD010, FDPart3tmp71));
        const REAL_SIMD_ARRAY FDPart3tmp166 = FusedMulAddSIMD(FDPart3tmp164, hDD22, MulSIMD(FDPart3tmp35, hDD_dD221));
        const REAL_SIMD_ARRAY FDPart3tmp183 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp28, f0_of_xx0), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, hDD22)),
            FusedMulSubSIMD(FDPart3tmp28, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0)),
                            FDPart3tmp107));
        const REAL_SIMD_ARRAY FDPart3tmp184 =
            FusedMulAddSIMD(MulSIMD(FDPart3tmp7, f1_of_xx1), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1__D1, hDD22)),
                            FusedMulSubSIMD(FDPart3tmp7, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, f1_of_xx1__D1)),
                                            MulSIMD(FDPart3tmp35, hDD_dD221)));
        const REAL_SIMD_ARRAY FDPart3tmp216 =
            FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp214, MulSIMD(FDPart3tmp163, MulSIMD(f1_of_xx1__D1, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp218 =
            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, f0_of_xx0__D0),
                               FusedMulAddSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, hDD11)),
                                               FusedMulSubSIMD(FDPart3tmp105, hDD_dD011, FDPart3tmp127)));
        const REAL_SIMD_ARRAY FDPart3tmp230 = MulSIMD(FDPart3tmp112, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp248 = FusedMulAddSIMD(FDPart3tmp10, hDD02, FDPart3tmp214);
        const REAL_SIMD_ARRAY FDPart3tmp256 = FusedMulAddSIMD(FDPart3tmp7, MulSIMD(f1_of_xx1__D1, hDD12), FDPart3tmp214);
        const REAL_SIMD_ARRAY FDPart3tmp265 = FusedMulAddSIMD(FDPart3tmp163, hDD01, FDPart3tmp213);
        const REAL_SIMD_ARRAY FDPart3tmp333 = FusedMulAddSIMD(FDPart3tmp163, MulSIMD(FDPart3tmp331, hDD02),
                                                              MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp35), MulSIMD(f1_of_xx1__D1, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp339 =
            FusedMulAddSIMD(hDD02, FusedMulAddSIMD(FDPart3tmp2, f1_of_xx1__D1, MulSIMD(FDPart3tmp50, f1_of_xx1__D1)),
                            FusedMulAddSIMD(FDPart3tmp19, MulSIMD(f1_of_xx1__D1, hDD_dD020),
                                            FusedMulAddSIMD(FDPart3tmp20, hDD_dDD0201, MulSIMD(FDPart3tmp52, hDD_dD021))));
        const REAL_SIMD_ARRAY FDPart3tmp340 = FusedMulSubSIMD(FDPart3tmp276, FDPart3tmp28, FDPart3tmp28);
        const REAL_SIMD_ARRAY FDPart3tmp365 = MulSIMD(MulSIMD(FDPart3tmp163, FDPart3tmp170), MulSIMD(FDPart3tmp259, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp371 =
            MulSIMD(FDPart3tmp170, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp272, f1_of_xx1__D1)));
        const REAL_SIMD_ARRAY FDPart3tmp372 = SubSIMD(FDPart3_Integer_1, FDPart3tmp276);
        const REAL_SIMD_ARRAY FDPart3tmp397 = FusedMulAddSIMD(FDPart3tmp10, hDD_dDD1201, FusedMulAddSIMD(FDPart3tmp112, hDD_dD121, FDPart3tmp254));
        const REAL_SIMD_ARRAY FDPart3tmp412 = MulSIMD(FDPart3_Integer_2, FDPart3tmp188);
        const REAL_SIMD_ARRAY FDPart3tmp13 =
            FusedMulSubSIMD(MulSIMD(FDPart3tmp7, FDPart3tmp8), MulSIMD(hDD01, hDD02), MulSIMD(FDPart3tmp11, FDPart3tmp4));
        const REAL_SIMD_ARRAY FDPart3tmp26 =
            MulSIMD(FDPart3tmp17, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp25, f0_of_xx0__D0)));
        const REAL_SIMD_ARRAY FDPart3tmp48 =
            FusedMulAddSIMD(FDPart3tmp34, FDPart3tmp47, MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp43), MulSIMD(f1_of_xx1, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp54 = FusedMulAddSIMD(FDPart3tmp20, hDD_dD020, MulSIMD(FDPart3tmp52, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp61 =
            MulSIMD(FDPart3tmp17, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp60, f0_of_xx0__D0)));
        const REAL_SIMD_ARRAY FDPart3tmp77 =
            MulSIMD(FDPart3tmp17, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp76, f0_of_xx0__D0)));
        const REAL_SIMD_ARRAY FDPart3tmp79 =
            FusedMulAddSIMD(FDPart3tmp39, FDPart3tmp78, MulSIMD(MulSIMD(FDPart3tmp1, FDPart3tmp28), MulSIMD(FDPart3tmp43, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp84 = SubSIMD(NegFusedMulAddSIMD(FDPart3tmp2, hDD01, FDPart3tmp83), MulSIMD(FDPart3tmp50, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp93 = FusedMulAddSIMD(FDPart3tmp34, FDPart3tmp39, FDPart3tmp32);
        const REAL_SIMD_ARRAY FDPart3tmp109 = MulSIMD(FDPart3_Rational_1_2, AddSIMD(FDPart3tmp106, FDPart3tmp108));
        const REAL_SIMD_ARRAY FDPart3tmp114 =
            MulSIMD(FDPart3_Rational_1_2, AddSIMD(FDPart3tmp113, FusedMulAddSIMD(FDPart3tmp46, hDD_dD021, FDPart3tmp59)));
        const REAL_SIMD_ARRAY FDPart3tmp131 = AddSIMD(FDPart3tmp113, FDPart3tmp57);
        const REAL_SIMD_ARRAY FDPart3tmp138 = FusedMulAddSIMD(FDPart3tmp112, hDD_dD020, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp52, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp167 = AddSIMD(FDPart3tmp164, FDPart3tmp166);
        const REAL_SIMD_ARRAY FDPart3tmp168 =
            FusedMulAddSIMD(MulSIMD(f0_of_xx0, f0_of_xx0__D0), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, hDD12)),
                            SubSIMD(FDPart3tmp57, FDPart3tmp111));
        const REAL_SIMD_ARRAY FDPart3tmp262 =
            FusedMulAddSIMD(FDPart3tmp256, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp261, FDPart3tmp47, MulSIMD(FDPart3tmp172, FDPart3tmp57)));
        const REAL_SIMD_ARRAY FDPart3tmp270 = FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp76, MulSIMD(FDPart3tmp66, FDPart3tmp73));
        const REAL_SIMD_ARRAY FDPart3tmp375 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp11, FDPart3tmp261));
        const REAL_SIMD_ARRAY FDPart3tmp42 = DivSIMD(
            FDPart3_Integer_1,
            FusedMulAddSIMD(FDPart3tmp39, FDPart3tmp40,
                            FusedMulAddSIMD(FDPart3tmp34, MulSIMD(FDPart3tmp39, FDPart3tmp4),
                                            FusedMulAddSIMD(hDD01, MulSIMD(MulSIMD(FDPart3tmp27, FDPart3tmp30), MulSIMD(hDD02, hDD12)),
                                                            FusedMulAddSIMD(FDPart3tmp32, FDPart3tmp4, MulSIMD(FDPart3tmp34, FDPart3tmp37))))));
        const REAL_SIMD_ARRAY FDPart3tmp55 = SubSIMD(NegFusedMulAddSIMD(FDPart3tmp51, hDD02, FDPart3tmp54), MulSIMD(FDPart3tmp8, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp68 = FusedMulAddSIMD(FDPart3tmp39, FDPart3tmp4, FDPart3tmp37);
        const REAL_SIMD_ARRAY FDPart3tmp250 = FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp25, MulSIMD(FDPart3tmp248, FDPart3tmp85));
        const REAL_SIMD_ARRAY FDPart3tmp269 =
            FusedMulAddSIMD(FDPart3tmp113, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp54, MulSIMD(FDPart3tmp11, FDPart3tmp89)));
        const REAL_SIMD_ARRAY FDPart3tmp274 = FusedMulAddSIMD(FDPart3tmp188, FDPart3tmp84, MulSIMD(FDPart3tmp272, FDPart3tmp85));
        const REAL_SIMD_ARRAY FDPart3tmp104 = MulSIMD(FDPart3tmp42, FDPart3tmp48);
        const REAL_SIMD_ARRAY FDPart3tmp110 = MulSIMD(FDPart3tmp42, FDPart3tmp79);
        const REAL_SIMD_ARRAY FDPart3tmp117 = MulSIMD(FDPart3tmp13, FDPart3tmp42);
        const REAL_SIMD_ARRAY FDPart3tmp118 = MulSIMD(FDPart3tmp42, FDPart3tmp69);
        const REAL_SIMD_ARRAY FDPart3tmp120 = MulSIMD(FDPart3tmp42, FDPart3tmp68);
        const REAL_SIMD_ARRAY FDPart3tmp125 = MulSIMD(FDPart3tmp42, FDPart3tmp93);
        const REAL_SIMD_ARRAY FDPart3tmp252 = FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp55, MulSIMD(FDPart3tmp46, hDD_dD120));
        const REAL_SIMD_ARRAY FDPart3tmp115 = FusedMulAddSIMD(FDPart3tmp104, FDPart3tmp109, MulSIMD(FDPart3tmp110, FDPart3tmp114));
        const REAL_SIMD_ARRAY FDPart3tmp119 =
            FusedMulAddSIMD(FDPart3tmp109, FDPart3tmp118, FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp117, FDPart3tmp85));
        const REAL_SIMD_ARRAY FDPart3tmp121 = FusedMulAddSIMD(FDPart3tmp109, FDPart3tmp117, MulSIMD(FDPart3tmp114, FDPart3tmp120));
        const REAL_SIMD_ARRAY FDPart3tmp124 = MulSIMD(FDPart3_Integer_3, FDPart3tmp118);
        const REAL_SIMD_ARRAY FDPart3tmp126 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp125);
        const REAL_SIMD_ARRAY FDPart3tmp130 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp110);
        const REAL_SIMD_ARRAY FDPart3tmp132 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp104);
        const REAL_SIMD_ARRAY FDPart3tmp135 = MulSIMD(FDPart3_Integer_3, FDPart3tmp117);
        const REAL_SIMD_ARRAY FDPart3tmp143 = MulSIMD(FDPart3_Integer_3, FDPart3tmp104);
        const REAL_SIMD_ARRAY FDPart3tmp144 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp117);
        const REAL_SIMD_ARRAY FDPart3tmp145 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp118);
        const REAL_SIMD_ARRAY FDPart3tmp148 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp120);
        const REAL_SIMD_ARRAY FDPart3tmp158 = MulSIMD(FDPart3_Integer_3, FDPart3tmp120);
        const REAL_SIMD_ARRAY FDPart3tmp160 = MulSIMD(FDPart3_Integer_3, FDPart3tmp110);
        const REAL_SIMD_ARRAY FDPart3tmp162 = MulSIMD(FDPart3_Integer_3, FDPart3tmp125);
        const REAL_SIMD_ARRAY FDPart3tmp133 =
            FusedMulAddSIMD(FDPart3tmp129, FDPart3tmp130, FusedMulAddSIMD(FDPart3tmp131, FDPart3tmp132, MulSIMD(FDPart3tmp126, FDPart3tmp71)));
        const REAL_SIMD_ARRAY FDPart3tmp141 =
            FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp140,
                            FusedMulAddSIMD(FDPart3tmp132, FDPart3tmp138, FusedMulAddSIMD(FDPart3tmp126, FDPart3tmp139, FDPart3tmp137)));
        const REAL_SIMD_ARRAY FDPart3tmp146 =
            FusedMulAddSIMD(FDPart3tmp131, FDPart3tmp145, FusedMulAddSIMD(FDPart3tmp132, FDPart3tmp71, MulSIMD(FDPart3tmp129, FDPart3tmp144)));
        const REAL_SIMD_ARRAY FDPart3tmp149 = FusedMulAddSIMD(
            FDPart3tmp130, FDPart3tmp71, FusedMulAddSIMD(FDPart3tmp131, FDPart3tmp144, FusedMulAddSIMD(FDPart3tmp129, FDPart3tmp148, FDPart3tmp85)));
        const REAL_SIMD_ARRAY FDPart3tmp152 =
            FusedMulAddSIMD(FDPart3tmp138, FDPart3tmp145, FusedMulAddSIMD(FDPart3tmp140, FDPart3tmp144, MulSIMD(FDPart3tmp132, FDPart3tmp139)));
        const REAL_SIMD_ARRAY FDPart3tmp153 =
            FusedMulAddSIMD(FDPart3tmp138, FDPart3tmp144, FusedMulAddSIMD(FDPart3tmp140, FDPart3tmp148, MulSIMD(FDPart3tmp130, FDPart3tmp139)));
        const REAL_SIMD_ARRAY FDPart3tmp169 = FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp168, MulSIMD(FDPart3tmp144, FDPart3tmp167));
        const REAL_SIMD_ARRAY FDPart3tmp173 =
            FusedMulAddSIMD(FDPart3tmp132, FDPart3tmp168, FusedMulAddSIMD(FDPart3tmp145, FDPart3tmp167, FDPart3tmp172));
        const REAL_SIMD_ARRAY FDPart3tmp174 = FusedMulAddSIMD(FDPart3tmp126, FDPart3tmp168, MulSIMD(FDPart3tmp132, FDPart3tmp167));
        const REAL_SIMD_ARRAY FDPart3tmp177 =
            FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp74, FusedMulAddSIMD(FDPart3tmp121, FDPart3tmp34, MulSIMD(FDPart3tmp11, FDPart3tmp119)));
        const REAL_SIMD_ARRAY FDPart3tmp185 = FusedMulAddSIMD(FDPart3tmp132, FDPart3tmp183, MulSIMD(FDPart3tmp144, FDPart3tmp184));
        const REAL_SIMD_ARRAY FDPart3tmp186 =
            FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp183, FusedMulAddSIMD(FDPart3tmp148, FDPart3tmp184, FDPart3tmp73));
        const REAL_SIMD_ARRAY FDPart3tmp189 =
            FusedMulAddSIMD(FDPart3tmp126, FDPart3tmp183, FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp184, FDPart3tmp188));
        const REAL_SIMD_ARRAY FDPart3tmp219 =
            FusedMulAddSIMD(FDPart3tmp144, FDPart3tmp213, FusedMulAddSIMD(FDPart3tmp145, FDPart3tmp216, MulSIMD(FDPart3tmp132, FDPart3tmp218)));
        const REAL_SIMD_ARRAY FDPart3tmp220 =
            FusedMulAddSIMD(FDPart3tmp144, FDPart3tmp216, FusedMulAddSIMD(FDPart3tmp148, FDPart3tmp213, MulSIMD(FDPart3tmp130, FDPart3tmp218)));
        const REAL_SIMD_ARRAY FDPart3tmp221 =
            FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp213,
                            FusedMulAddSIMD(FDPart3tmp132, FDPart3tmp216, FusedMulAddSIMD(FDPart3tmp126, FDPart3tmp218, FDPart3tmp187)));
        const REAL_SIMD_ARRAY FDPart3tmp122 =
            FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp121, FDPart3tmp74, MulSIMD(FDPart3tmp115, FDPart3tmp4)));
        const REAL_SIMD_ARRAY FDPart3tmp175 =
            FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp4, MulSIMD(FDPart3tmp169, FDPart3tmp74)));
        const REAL_SIMD_ARRAY FDPart3tmp178 = MulSIMD(FDPart3tmp121, FDPart3tmp177);
        const REAL_SIMD_ARRAY FDPart3tmp181 =
            FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp39, MulSIMD(FDPart3tmp11, FDPart3tmp121)));
        const REAL_SIMD_ARRAY FDPart3tmp190 =
            FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp74, FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp4, MulSIMD(FDPart3tmp185, FDPart3tmp90)));
        const REAL_SIMD_ARRAY FDPart3tmp196 = MulSIMD(FDPart3_Integer_2, FDPart3tmp153);
        const REAL_SIMD_ARRAY FDPart3tmp199 = MulSIMD(FDPart3_Integer_2, FDPart3tmp152);
        const REAL_SIMD_ARRAY FDPart3tmp201 =
            FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp39, MulSIMD(FDPart3tmp11, FDPart3tmp149)));
        const REAL_SIMD_ARRAY FDPart3tmp204 =
            FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp39, MulSIMD(FDPart3tmp11, FDPart3tmp153)));
        const REAL_SIMD_ARRAY FDPart3tmp222 =
            FusedMulAddSIMD(FDPart3tmp220, FDPart3tmp74, FusedMulAddSIMD(FDPart3tmp221, FDPart3tmp4, MulSIMD(FDPart3tmp219, FDPart3tmp90)));
        const REAL_SIMD_ARRAY FDPart3tmp225 =
            FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp74, FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp34, MulSIMD(FDPart3tmp11, FDPart3tmp152)));
        const REAL_SIMD_ARRAY FDPart3tmp237 = FusedMulAddSIMD(
            FDPart3tmp117, FDPart3tmp174,
            FusedMulAddSIMD(FDPart3tmp126, FDPart3tmp141,
                            FusedMulAddSIMD(FDPart3tmp145, FDPart3tmp189,
                                            FusedMulAddSIMD(FDPart3tmp148, FDPart3tmp221,
                                                            FusedMulAddSIMD(FDPart3tmp104, FDPart3tmp115, MulSIMD(FDPart3tmp110, FDPart3tmp133))))));
        const REAL_SIMD_ARRAY FDPart3tmp300 =
            FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp34, FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp74, MulSIMD(FDPart3tmp11, FDPart3tmp185)));
        const REAL_SIMD_ARRAY FDPart3tmp309 =
            FusedMulAddSIMD(FDPart3tmp219, FDPart3tmp39, FusedMulAddSIMD(FDPart3tmp221, FDPart3tmp90, MulSIMD(FDPart3tmp11, FDPart3tmp220)));
        const REAL_SIMD_ARRAY FDPart3tmp387 = MulSIMD(FDPart3_Integer_2, FDPart3tmp221);
        const REAL_SIMD_ARRAY FDPart3tmp389 = MulSIMD(FDPart3_Integer_2, FDPart3tmp219);
        const REAL_SIMD_ARRAY FDPart3tmp417 = MulSIMD(FDPart3_Integer_2, FDPart3tmp186);
        const REAL_SIMD_ARRAY FDPart3tmp123 = MulSIMD(FDPart3tmp115, FDPart3tmp122);
        const REAL_SIMD_ARRAY FDPart3tmp134 = MulSIMD(FDPart3tmp122, FDPart3tmp133);
        const REAL_SIMD_ARRAY FDPart3tmp150 =
            FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp74, MulSIMD(FDPart3tmp133, FDPart3tmp4)));
        const REAL_SIMD_ARRAY FDPart3tmp155 =
            FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp74, MulSIMD(FDPart3tmp141, FDPart3tmp4)));
        const REAL_SIMD_ARRAY FDPart3tmp182 = MulSIMD(FDPart3tmp119, FDPart3tmp181);
        const REAL_SIMD_ARRAY FDPart3tmp206 = MulSIMD(FDPart3tmp146, FDPart3tmp201);
        const REAL_SIMD_ARRAY FDPart3tmp211 =
            FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp74, FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp34, MulSIMD(FDPart3tmp11, FDPart3tmp146)));
        const REAL_SIMD_ARRAY FDPart3tmp234 = FusedMulAddSIMD(
            FDPart3tmp117, FDPart3tmp173,
            FusedMulAddSIMD(FDPart3tmp126, FDPart3tmp152,
                            FusedMulAddSIMD(FDPart3tmp145, FDPart3tmp185,
                                            FusedMulAddSIMD(FDPart3tmp148, FDPart3tmp219,
                                                            FusedMulAddSIMD(FDPart3tmp104, FDPart3tmp119, MulSIMD(FDPart3tmp110, FDPart3tmp146))))));
        const REAL_SIMD_ARRAY FDPart3tmp236 = FusedMulAddSIMD(
            FDPart3tmp117, FDPart3tmp169,
            FusedMulAddSIMD(FDPart3tmp126, FDPart3tmp153,
                            FusedMulAddSIMD(FDPart3tmp145, FDPart3tmp186,
                                            FusedMulAddSIMD(FDPart3tmp148, FDPart3tmp220,
                                                            FusedMulAddSIMD(FDPart3tmp104, FDPart3tmp121, MulSIMD(FDPart3tmp110, FDPart3tmp149))))));
        const REAL_SIMD_ARRAY FDPart3tmp283 =
            FusedMulAddSIMD(FDPart3tmp169, FDPart3tmp34, FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp74, MulSIMD(FDPart3tmp11, FDPart3tmp173)));
        const REAL_SIMD_ARRAY FDPart3tmp290 =
            FusedMulAddSIMD(FDPart3tmp220, FDPart3tmp34, FusedMulAddSIMD(FDPart3tmp221, FDPart3tmp74, MulSIMD(FDPart3tmp11, FDPart3tmp219)));
        const REAL_SIMD_ARRAY FDPart3tmp292 = MulSIMD(FDPart3tmp115, FDPart3tmp175);
        const REAL_SIMD_ARRAY FDPart3tmp294 = MulSIMD(FDPart3tmp122, FDPart3tmp174);
        const REAL_SIMD_ARRAY FDPart3tmp297 =
            FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp39, FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp90, MulSIMD(FDPart3tmp11, FDPart3tmp169)));
        const REAL_SIMD_ARRAY FDPart3tmp302 = MulSIMD(FDPart3tmp133, FDPart3tmp175);
        const REAL_SIMD_ARRAY FDPart3tmp344 =
            FusedMulAddSIMD(FDPart3tmp185, FDPart3tmp39, FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp90, MulSIMD(FDPart3tmp11, FDPart3tmp186)));
        const REAL_SIMD_ARRAY FDPart3tmp354 = FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp181, MulSIMD(FDPart3tmp133, FDPart3tmp190));
        const REAL_SIMD_ARRAY FDPart3tmp361 = FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp177, MulSIMD(FDPart3tmp149, FDPart3tmp201));
        const REAL_SIMD_ARRAY FDPart3tmp382 = MulSIMD(FDPart3tmp174, FDPart3tmp175);
        const REAL_SIMD_ARRAY FDPart3tmp151 = MulSIMD(FDPart3tmp115, FDPart3tmp150);
        const REAL_SIMD_ARRAY FDPart3tmp157 = MulSIMD(FDPart3tmp133, FDPart3tmp150);
        const REAL_SIMD_ARRAY FDPart3tmp212 = MulSIMD(FDPart3tmp121, FDPart3tmp211);
        const REAL_SIMD_ARRAY FDPart3tmp226 = MulSIMD(FDPart3tmp149, FDPart3tmp211);
        const REAL_SIMD_ARRAY FDPart3tmp284 = MulSIMD(FDPart3tmp121, FDPart3tmp283);
        const REAL_SIMD_ARRAY FDPart3tmp285 = MulSIMD(FDPart3tmp149, FDPart3tmp283);
        const REAL_SIMD_ARRAY FDPart3tmp298 = MulSIMD(FDPart3tmp119, FDPart3tmp297);
        const REAL_SIMD_ARRAY FDPart3tmp311 = MulSIMD(FDPart3tmp146, FDPart3tmp297);
        const REAL_SIMD_ARRAY FDPart3tmp318 = MulSIMD(FDPart3tmp150, FDPart3tmp174);
        const REAL_SIMD_ARRAY FDPart3tmp345 = MulSIMD(FDPart3_Integer_2, FDPart3tmp344);
        const REAL_SIMD_ARRAY FDPart3tmp346 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp146, FDPart3tmp181));
        const REAL_SIMD_ARRAY FDPart3tmp355 = FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp297, MulSIMD(FDPart3tmp149, FDPart3tmp300));
        const REAL_SIMD_ARRAY FDPart3tmp360 = FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp201, FDPart3tmp302);
        const REAL_SIMD_ARRAY FDPart3tmp363 = FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp204, FDPart3tmp134);
        const REAL_SIMD_ARRAY FDPart3tmp378 = MulSIMD(FDPart3tmp169, FDPart3tmp283);
        const REAL_SIMD_ARRAY FDPart3tmp384 = MulSIMD(FDPart3tmp173, FDPart3tmp297);
        const REAL_SIMD_ARRAY FDPart3tmp362 = FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp309, FDPart3tmp285);
        const REAL_SIMD_ARRAY __RHS_exp_0 =
            FusedMulAddSIMD(
                FDPart3tmp42,
                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                        MulSIMD(FDPart3tmp79,
                                FusedMulAddSIMD(
                                    FDPart3tmp78, FDPart3tmp91,
                                    FusedMulAddSIMD(
                                        FDPart3tmp100, hDD_dD001,
                                        FusedMulAddSIMD(FDPart3tmp74, FDPart3tmp89,
                                                        FusedMulAddSIMD(FDPart3tmp5,
                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                MulSIMD(FDPart3tmp72, f0_of_xx0__DD00)),
                                                                        NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp17),
                                                                                           MulSIMD(FDPart3tmp83, f0_of_xx0__D0), FDPart3tmp86))))))),
                FusedMulAddSIMD(
                    FDPart3tmp42,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp69,
                                    FusedMulAddSIMD(FDPart3tmp28, FDPart3tmp63, FusedMulAddSIMD(FDPart3tmp72, FDPart3tmp73, FDPart3tmp77)))),
                    FusedMulAddSIMD(
                        FDPart3tmp42,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp79, NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp17),
                                                                         MulSIMD(FDPart3tmp84, f0_of_xx0__D0), FDPart3tmp86))),
                        FusedMulAddSIMD(
                            FDPart3tmp42,
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                    MulSIMD(FDPart3tmp48,
                                            FusedMulAddSIMD(
                                                FDPart3tmp47, FDPart3tmp91,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp89, FDPart3tmp90,
                                                    FusedMulAddSIMD(MulSIMD(FDPart3_Integer_4, FDPart3tmp1), MulSIMD(f0_of_xx0__DD00, f1_of_xx1),
                                                                    NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp17),
                                                                                       MulSIMD(FDPart3tmp54, f0_of_xx0__D0), FDPart3tmp49)))))),
                            FusedMulAddSIMD(
                                FDPart3tmp42,
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                        MulSIMD(FDPart3tmp68,
                                                AddSIMD(FDPart3tmp67,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp2, hDD_dDD0011,
                                                            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp2, hDD_dD011), FDPart3tmp63))))),
                                FusedMulAddSIMD(FDPart3tmp13,
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                        MulSIMD(FDPart3tmp42,
                                                                AddSIMD(FDPart3tmp61,
                                                                        NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp17),
                                                                                           MulSIMD(FDPart3tmp57, f0_of_xx0__D0), FDPart3tmp16)))),
                                                FusedMulAddSIMD(FDPart3tmp42,
                                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                        MulSIMD(FDPart3tmp48,
                                                                                NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp17),
                                                                                                   MulSIMD(FDPart3tmp55, f0_of_xx0__D0),
                                                                                                   FDPart3tmp49))),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp141, MulSIMD(FDPart3tmp155, FDPart3tmp162),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp4, MulSIMD(FDPart3tmp5, lambdaU_dD00),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp133, MulSIMD(FDPart3tmp155, FDPart3tmp160),
                                                                            FusedMulAddSIMD(FDPart3tmp141, MulSIMD(FDPart3tmp150, FDPart3tmp160),
                                                                                            FusedMulAddSIMD(FDPart3tmp115,
                                                                                                            MulSIMD(FDPart3tmp143, FDPart3tmp155),
                                                                                                            FusedMulAddSIMD(FDPart3tmp122,
                                                                                                                            MulSIMD(FDPart3tmp141,
                                                                                                                                    FDPart3tmp143),
                                                                                                                            FusedMulAddSIMD(FDPart3tmp236,
                                                                                                                                            FusedMulAddSIMD(FDPart3tmp149,
                                                                                                                                                            FDPart3tmp231,
                                                                                                                                                            FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                            MulSIMD(FDPart3tmp133,
                                                                                                                                                                                    FDPart3tmp4),
                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                FDPart3tmp146,
                                                                                                                                                                                FDPart3tmp230))),
                                                                                                                                            FusedMulAddSIMD(FDPart3tmp237,
                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp153,
                                                                                                                                                                            FDPart3tmp231,
                                                                                                                                                                            FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                FDPart3tmp141,
                                                                                                                                                                                                FDPart3tmp4),
                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                FDPart3tmp152,
                                                                                                                                                                                                FDPart3tmp230))),
                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp157,
                                                                                                                                                                            FDPart3tmp158,
                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp234,
                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp121, FDPart3tmp231, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp115, FDPart3tmp4), MulSIMD(FDPart3tmp119, FDPart3tmp230))),
                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp134,
                                                                                                                                                                                                            FDPart3tmp135,
                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp135,
                                                                                                                                                                                                                            FDPart3tmp151,
                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp125, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp199, MulSIMD(FDPart3tmp152, FDPart3tmp204)),
                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp125,
                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp150,
                                                                                                                                                                                                                                                                            FDPart3tmp196,
                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                FDPart3tmp153,
                                                                                                                                                                                                                                                                                FDPart3tmp225)),
                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp120, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp149, FDPart3tmp222), FDPart3tmp226),
                                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp123,
                                                                                                                                                                                                                                                                                            FDPart3tmp124,
                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp118, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp119, FDPart3tmp190), FDPart3tmp182),
                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                FDPart3tmp120,
                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp146, FDPart3tmp175),
                                                                                                                                                                                                                                                                                                                                FDPart3tmp206),
                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                    FDPart3tmp117,
                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                        FDPart3tmp149, FDPart3tmp177, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp149, FDPart3tmp175))),
                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                        FDPart3tmp118,
                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp121, FDPart3tmp175), FDPart3tmp178),
                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp117, FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp201, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp119, FDPart3tmp175))),
                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp117,
                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp146,
                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp181,
                                                                                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp146,
                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp190))),
                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp110,
                                                                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp153,
                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp211,
                                                                                                                                                                                                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp196,
                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp222)),
                                                                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp117,
                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp121,
                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp222),
                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp212),
                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp149,
                                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp225,
                                                                                                                                                                                                                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp149,
                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp150))),
                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp201, MulSIMD(FDPart3tmp175, FDPart3tmp199)),
                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp177, MulSIMD(FDPart3tmp175, FDPart3tmp196)),
                                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp204, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp122, FDPart3tmp146))),
                                                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp121, FDPart3tmp225, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp121, FDPart3tmp150))),
                                                                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp181, MulSIMD(FDPart3tmp190, FDPart3tmp199)),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp1,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    lambdaU_dD20,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp119, FDPart3tmp204, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp119, FDPart3tmp122))),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp42,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    FDPart3_Rational_1_2),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            MulSIMD(FDPart3tmp93, FusedMulAddSIMD(f0_of_xx0__D0,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  MulSIMD(f0_of_xx0__DD00, hDD_dD000),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      FDPart3tmp3,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      FDPart3tmp99,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      FusedMulAddSIMD(hDD00,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      MulSIMD(f0_of_xx0__DD00, f0_of_xx0__DD00), MulSIMD(FDPart3tmp94, f0_of_xx0__DDD000)),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          FDPart3tmp101,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              MulSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      FDPart3_NegativeOne_),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              MulSIMD(FDPart3tmp5,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      f0_of_xx0__DD00)),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          FusedMulSubSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              FDPart3tmp2,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              hDD_dDD0000,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  FDPart3tmp3,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  FusedMulSubSIMD(FDPart3tmp5, f0_of_xx0__DDD000, FDPart3tmp98))))))))),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulSubSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp0,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        lambdaU_dD10,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                FDPart3_Rational_1_2,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp13),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp42,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                AddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp16, FDPart3tmp26))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(FDPart3tmp42, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp66, FDPart3tmp85, AddSIMD(FusedMulAddSIMD(FDPart3tmp45, hDD_dD110, FDPart3tmp63), SubSIMD(NegFusedMulAddSIMD(FDPart3tmp2, hDD_dD011, FDPart3tmp251), MulSIMD(FDPart3tmp50, hDD_dD011)))))), FusedMulAddSIMD(FDPart3tmp42, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp128, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp137, FDPart3tmp66, AddSIMD(FusedMulAddSIMD(FDPart3tmp101, FDPart3tmp187, FDPart3tmp67), FusedMulAddSIMD(FDPart3tmp33, FDPart3tmp89, NegFusedMulAddSIMD(FDPart3tmp3, AddSIMD(FDPart3_NegativeOne_, FDPart3tmp276), FDPart3tmp251))))))), FusedMulAddSIMD(FDPart3tmp42, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp68, FusedMulAddSIMD(FDPart3tmp19, hDD_dD001, FusedMulAddSIMD(FDPart3tmp19, hDD_dDD0111, FusedMulAddSIMD(FDPart3tmp265, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp45, hDD_dD111, FusedMulAddSIMD(FDPart3tmp187, FDPart3tmp72, MulSIMD(FDPart3tmp187, FDPart3tmp84)))))))), FusedMulAddSIMD(FDPart3tmp42, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp69, AddSIMD(FDPart3tmp270, FDPart3tmp274))), FusedMulAddSIMD(FDPart3tmp42, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp48, FusedMulAddSIMD(FDPart3tmp60, FDPart3tmp85, FDPart3tmp252))), FusedMulAddSIMD(FDPart3tmp42, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp48, AddSIMD(FDPart3tmp269, FusedMulAddSIMD(FDPart3tmp137, FDPart3tmp60, FDPart3tmp61)))), FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, f0_of_xx0), MulSIMD(hDD01, lambdaU_dD00), FusedMulAddSIMD(FDPart3tmp13, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp42, AddSIMD(FDPart3tmp262, FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp60, NegFusedMulAddSIMD(f0_of_xx0__D0, MulSIMD(MulSIMD(FDPart3_Integer_2, f0_of_xx0), MulSIMD(f1_of_xx1, hDD02)), FDPart3tmp254))))), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp245, FDPart3tmp4), FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp17), MulSIMD(FDPart3tmp34, lambdaU_dD10), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp23, lambdaU_dD20), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp243, FDPart3tmp74), FusedMulAddSIMD(FDPart3tmp237, AddSIMD(FDPart3tmp150, FDPart3tmp225), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp1, lambdaU_dD21), FusedMulAddSIMD(FDPart3tmp234, AddSIMD(FDPart3tmp175, FDPart3tmp177), FusedMulAddSIMD(FDPart3tmp236, AddSIMD(FDPart3tmp211, FDPart3tmp222), FusedMulAddSIMD(FDPart3tmp125, FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp177, FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp201, MulSIMD(FDPart3tmp122, FDPart3tmp146))), FusedMulAddSIMD(FDPart3tmp125, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp150, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp225, MulSIMD(FDPart3tmp133, FDPart3tmp155))), FusedMulAddSIMD(FDPart3tmp120, FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp309, FusedMulAddSIMD(FDPart3tmp175, FDPart3tmp219, MulSIMD(FDPart3tmp146, FDPart3tmp283))), FusedMulAddSIMD(FDPart3tmp125, FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp150, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp153, FDPart3tmp211))), FusedMulAddSIMD(FDPart3tmp120, FusedMulAddSIMD(FDPart3tmp220, FDPart3tmp222, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp149, FDPart3tmp290))), FusedMulAddSIMD(FDPart3tmp120, FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp222, FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp221, MulSIMD(FDPart3tmp133, FDPart3tmp211))), FusedMulAddSIMD(FDPart3tmp118, AddSIMD(FDPart3tmp292, FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp177, FDPart3tmp294)), FusedMulAddSIMD(FDPart3tmp118, FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp300, FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp190, FDPart3tmp298)), FusedMulAddSIMD(FDPart3tmp117, FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp309, FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp175, MulSIMD(FDPart3tmp119, FDPart3tmp283))), FusedMulAddSIMD(FDPart3tmp118, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp284, MulSIMD(FDPart3tmp169, FDPart3tmp175)), FusedMulAddSIMD(FDPart3tmp117, FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp300, FusedMulAddSIMD(FDPart3tmp190, FDPart3tmp219, FDPart3tmp311)), FusedMulAddSIMD(FDPart3tmp117, FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp211, FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp222, FDPart3tmp318)), FusedMulAddSIMD(FDPart3tmp117, FusedMulAddSIMD(FDPart3tmp169, FDPart3tmp222, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp121, FDPart3tmp290))), FusedMulAddSIMD(FDPart3tmp117, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp221, FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp177, FDPart3tmp302)), FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp283, FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp309, MulSIMD(FDPart3tmp146, FDPart3tmp175))), FusedMulAddSIMD(FDPart3tmp117, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp285, MulSIMD(FDPart3tmp175, FDPart3tmp220)), FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp211, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp222, FDPart3tmp157)), FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp219, FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp177, FDPart3tmp206)), FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp222, MulSIMD(FDPart3tmp196, FDPart3tmp290)), FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp225, FusedMulAddSIMD(FDPart3tmp155, FDPart3tmp221, FDPart3tmp157)), FusedMulAddSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp297, FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp300, MulSIMD(FDPart3tmp146, FDPart3tmp190))), FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp226, MulSIMD(FDPart3tmp150, FDPart3tmp220)), FusedMulAddSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp225, FusedMulAddSIMD(FDPart3tmp155, FDPart3tmp174, FDPart3tmp151)), FusedMulAddSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp201, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp173, MulSIMD(FDPart3tmp119, FDPart3tmp177))), FusedMulAddSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp175, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp153, FDPart3tmp283))), FusedMulAddSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp175, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp177, FDPart3tmp134)), FusedMulAddSIMD(FDPart3tmp42, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp93, FusedMulAddSIMD(FDPart3tmp84, FDPart3tmp85, FusedMulAddSIMD(hDD01, FusedMulAddSIMD(f0_of_xx0, f0_of_xx0__DDD000, FDPart3tmp279), FusedMulAddSIMD(FDPart3tmp78, FDPart3tmp91, FusedMulAddSIMD(FDPart3tmp83, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp19, hDD_dDD0100, FusedMulAddSIMD(FDPart3tmp74, FDPart3tmp99, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp81, hDD_dD010), FusedMulSubSIMD(FDPart3tmp137, FDPart3tmp83, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp5), MulSIMD(FDPart3tmp84, f0_of_xx0__DD00)))))))))))), FusedMulSubSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp212, MulSIMD(FDPart3tmp150, FDPart3tmp169)), MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp13), MulSIMD(FDPart3tmp42, FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp60, FDPart3tmp250)))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_2 = FusedMulAddSIMD(
    f1_of_xx1, MulSIMD(MulSIMD(FDPart3_Rational_1_2, f0_of_xx0), MulSIMD(hDD02, lambdaU_dD00)),
    FusedMulAddSIMD(
        MulSIMD(FDPart3tmp4, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(f1_of_xx1, lambdaU2)),
        FusedMulAddSIMD(
            FDPart3tmp42,
            MulSIMD(
                MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                MulSIMD(FDPart3tmp93,
                        FusedMulAddSIMD(
                            FDPart3tmp90, FDPart3tmp99,
                            FusedMulAddSIMD(
                                hDD02, FusedMulAddSIMD(FDPart3tmp279, f1_of_xx1, MulSIMD(f0_of_xx0, MulSIMD(f0_of_xx0__DDD000, f1_of_xx1))),
                                FusedMulAddSIMD(
                                    FDPart3tmp54, FDPart3tmp85,
                                    FusedMulAddSIMD(
                                        FDPart3tmp55, FDPart3tmp85,
                                        FusedMulAddSIMD(
                                            FDPart3tmp20, hDD_dDD0200,
                                            FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp91,
                                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp52, hDD_dD020),
                                                                            FusedMulSubSIMD(FDPart3tmp137, FDPart3tmp54,
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp5),
                                                                                                    MulSIMD(FDPart3tmp55, f0_of_xx0__DD00)))))))))))),
            FusedMulAddSIMD(
                FDPart3tmp170, MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp17), MulSIMD(FDPart3tmp39, lambdaU_dD20)),
                FusedMulAddSIMD(
                    FDPart3tmp42,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp79,
                                    AddSIMD(AddSIMD(FDPart3tmp26, FDPart3tmp269), FusedMulAddSIMD(FDPart3tmp137, FDPart3tmp25, FDPart3tmp339)))),
                    FusedMulAddSIMD(
                        FDPart3tmp42,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp137, FDPart3tmp57,
                                                                      FusedMulAddSIMD(FDPart3tmp25, FDPart3tmp85,
                                                                                      FusedMulAddSIMD(FDPart3tmp57, FDPart3tmp85,
                                                                                                      AddSIMD(FDPart3tmp252, FDPart3tmp339)))))),
                        FusedMulAddSIMD(
                            FDPart3tmp42,
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                    MulSIMD(FDPart3tmp68,
                                            FusedMulAddSIMD(
                                                FDPart3tmp187, FDPart3tmp55,
                                                FusedMulAddSIMD(FDPart3tmp20, hDD_dDD0211,
                                                                FusedMulAddSIMD(FDPart3tmp56, f1_of_xx1__DD11,
                                                                                FusedMulAddSIMD(FDPart3tmp105, MulSIMD(f1_of_xx1__D1, hDD_dD021),
                                                                                                AddSIMD(FDPart3tmp250, FDPart3tmp262))))))),
                            FusedMulAddSIMD(
                                FDPart3tmp42,
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                        MulSIMD(FDPart3tmp69,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp25, FDPart3tmp73,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp333, FDPart3tmp85,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp60, FDPart3tmp73,
                                                            FusedMulSubSIMD(FDPart3tmp188, FDPart3tmp55,
                                                                            MulSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp331),
                                                                                                       MulSIMD(f0_of_xx0__D0, hDD02))))))))),
                                FusedMulAddSIMD(
                                    FDPart3tmp42,
                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                            MulSIMD(FDPart3tmp48,
                                                    FusedMulAddSIMD(FDPart3tmp75, hDD_dD220,
                                                                    FusedMulAddSIMD(FDPart3tmp76, FDPart3tmp85,
                                                                                    FusedMulAddSIMD(FDPart3tmp28, FDPart3tmp63,
                                                                                                    MulSIMD(FDPart3tmp73, FDPart3tmp84)))))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp42,
                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                MulSIMD(FDPart3tmp48,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp137, FDPart3tmp76,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp38, FDPart3tmp89,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp101, FDPart3tmp188,
                                                                    FusedMulAddSIMD(FDPart3tmp108, FDPart3tmp85,
                                                                                    FusedMulAddSIMD(FDPart3tmp73, FDPart3tmp83,
                                                                                                    NegFusedMulAddSIMD(FDPart3tmp3, FDPart3tmp340,
                                                                                                                       FDPart3tmp77)))))))),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3_Rational_1_2, FDPart3tmp23), MulSIMD(f1_of_xx1, lambdaU_dD10),
                                            FusedMulAddSIMD(
                                                FDPart3tmp13,
                                                MulSIMD(
                                                    MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                    MulSIMD(FDPart3tmp42,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp105, MulSIMD(FDPart3tmp73, hDD00),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp19, MulSIMD(FDPart3tmp28, hDD_dD001),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp272, FDPart3tmp85,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp337, FDPart3tmp78,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp19, MulSIMD(FDPart3tmp73, hDD_dD011),
                                                                                FusedMulSubSIMD(FDPart3tmp166, FDPart3tmp85,
                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp170),
                                                                                                        MulSIMD(FDPart3tmp76, f1_of_xx1__D1)))))))))),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp237, AddSIMD(FDPart3tmp122, FDPart3tmp204),
                                                    FusedMulAddSIMD(
                                                        FDPart3_Rational_1_2, MulSIMD(FDPart3tmp326, FDPart3tmp90),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp234, AddSIMD(FDPart3tmp181, FDPart3tmp190),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp236, AddSIMD(FDPart3tmp175, FDPart3tmp201),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp125,
                                                                    FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp141,
                                                                                    FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp204,
                                                                                                    MulSIMD(FDPart3tmp115, FDPart3tmp155))),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp125,
                                                                        FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp177,
                                                                                        FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp201,
                                                                                                        MulSIMD(FDPart3tmp121, FDPart3tmp150))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp120,
                                                                            FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp311,
                                                                                            MulSIMD(FDPart3tmp173, FDPart3tmp175)),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp125,
                                                                                FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp122,
                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp152, FDPart3tmp181))),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp120, AddSIMD(FDPart3tmp318, FDPart3tmp360),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp120,
                                                                                        FusedMulAddSIMD(FDPart3tmp169, FDPart3tmp222, FDPart3tmp362),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp118,
                                                                                            FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp190,
                                                                                                            FusedMulAddSIMD(FDPart3tmp122,
                                                                                                                            FDPart3tmp189,
                                                                                                                            MulSIMD(FDPart3tmp115,
                                                                                                                                    FDPart3tmp181))),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp118,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp121, FDPart3tmp300,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp175, FDPart3tmp186,
                                                                                                        MulSIMD(FDPart3tmp121, FDPart3tmp297))),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp117,
                                                                                                    FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp201,
                                                                                                                    FusedMulAddSIMD(FDPart3tmp150,
                                                                                                                                    FDPart3tmp189,
                                                                                                                                    FDPart3tmp292)),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp118,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp119, FDPart3tmp345,
                                                                                                            MulSIMD(FDPart3tmp185, FDPart3tmp190)),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp117,
                                                                                                            FusedMulAddSIMD(FDPart3tmp146,
                                                                                                                            FDPart3tmp345,
                                                                                                                            MulSIMD(FDPart3tmp173,
                                                                                                                                    FDPart3tmp190)),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp117,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp121, FDPart3tmp309,
                                                                                                                    FusedMulAddSIMD(FDPart3tmp186,
                                                                                                                                    FDPart3tmp222,
                                                                                                                                    FDPart3tmp284)),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp117,
                                                                                                                    FusedMulAddSIMD(FDPart3tmp169,
                                                                                                                                    FDPart3tmp175,
                                                                                                                                    FDPart3tmp355),
                                                                                                                    FusedMulAddSIMD(FDPart3tmp117,
                                                                                                                                    FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                    FDPart3tmp298,
                                                                                                                                                    MulSIMD(
                                                                                                                                                        FDPart3tmp175,
                                                                                                                                                        FDPart3tmp185)),
                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp153,
                                                                                                                                                                    FDPart3tmp283,
                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp153,
                                                                                                                                                                                    FDPart3tmp309,
                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                        FDPart3tmp121,
                                                                                                                                                                                        FDPart3tmp222))),
                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp117,
                                                                                                                                                                    AddSIMD(
                                                                                                                                                                        FDPart3tmp294,
                                                                                                                                                                        FDPart3tmp354),
                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp119,
                                                                                                                                                                                                    FDPart3tmp175,
                                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                            FDPart3tmp152,
                                                                                                                                                                                                            FDPart3tmp297))),
                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp175, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp201, FDPart3tmp151)),
                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp150,
                                                                                                                                                                                                                                    FDPart3tmp169,
                                                                                                                                                                                                                                    FDPart3tmp361),
                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp155,
                                                                                                                                                                                                                                                    FDPart3tmp174,
                                                                                                                                                                                                                                                    FDPart3tmp363),
                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                        FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp297, FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp300, MulSIMD(FDPart3tmp121, FDPart3tmp175))),
                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp173, FDPart3tmp346),
                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                            FDPart3tmp104,
                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp141,
                                                                                                                                                                                                                                                                            FDPart3tmp181,
                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                FDPart3tmp141,
                                                                                                                                                                                                                                                                                FDPart3tmp190,
                                                                                                                                                                                                                                                                                FDPart3tmp123)),
                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp121, FDPart3tmp201, FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp186, FDPart3tmp178)),
                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp190, MulSIMD(FDPart3tmp199, FDPart3tmp344)),
                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                    FDPart3tmp115, FDPart3tmp204, FusedMulAddSIMD(FDPart3tmp155, FDPart3tmp189, FDPart3tmp123)),
                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(MulSIMD(
                                                                                                                                                                                                                                                                                                                    f0_of_xx0__D0,
                                                                                                                                                                                                                                                                                                                    f1_of_xx1__D1),
                                                                                                                                                                                                                                                                                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                 hDD01,
                                                                                                                                                                                                                                                                                                                                                                                 lambdaU2)),
                                                                                                                                                                                                                                                                                                                FusedMulSubSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp182, MulSIMD(FDPart3tmp122, FDPart3tmp185)),
                                                                                                                                                                                                                                                                                                                                MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp13), MulSIMD(FDPart3tmp42, AddSIMD(FDPart3tmp270, FusedMulAddSIMD(FDPart3tmp188,
                                                                                                                                                                                                                                                                                                                                                                                                                                                  FDPart3tmp72, FDPart3tmp330))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_3 = FusedMulAddSIMD(
    FDPart3tmp42,
    MulSIMD(
        MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
        MulSIMD(FDPart3tmp79,
                FusedMulAddSIMD(FDPart3tmp7, hDD_dDD1101,
                                FusedMulAddSIMD(FDPart3tmp368, FDPart3tmp83,
                                                FusedMulAddSIMD(FDPart3tmp372, FDPart3tmp74,
                                                                FusedMulAddSIMD(FDPart3tmp78, AddSIMD(FDPart3_NegativeOne_, FDPart3tmp276),
                                                                                FusedMulSubSIMD(FDPart3tmp105, hDD_dD111,
                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_3, FDPart3tmp17),
                                                                                                        MulSIMD(FDPart3tmp265, f0_of_xx0__D0))))))))),
    FusedMulAddSIMD(
        FDPart3tmp42,
        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                MulSIMD(FDPart3tmp69, FusedMulAddSIMD(FDPart3tmp265, FDPart3tmp73,
                                                      FusedMulAddSIMD(FDPart3tmp28, MulSIMD(FDPart3tmp366, hDD_dD110), FDPart3tmp371)))),
        FusedMulAddSIMD(FDPart3tmp42,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp7, hDD_dDD1101,
                                                                      FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp187, FDPart3tmp84),
                                                                                      MulSIMD(FDPart3tmp265, FDPart3tmp85))))),
                        FusedMulAddSIMD(FDPart3tmp42,
                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                MulSIMD(FDPart3tmp48,
                                                        FusedMulSubSIMD(MulSIMD(FDPart3_Integer_6, FDPart3tmp19), MulSIMD(f1_of_xx1__D1, hDD12),
                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp113),
                                                                                MulSIMD(FDPart3tmp170, f1_of_xx1__D1))))),
                                        FusedMulAddSIMD(FDPart3tmp42,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                MulSIMD(FDPart3tmp68,
                                                                        FusedMulAddSIMD(FDPart3tmp368, FDPart3tmp66,
                                                                                        FusedMulAddSIMD(FDPart3tmp7, hDD_dDD1111,
                                                                                                        FusedMulAddSIMD(FDPart3tmp163, hDD_dD011,
                                                                                                                        MulSIMD(FDPart3tmp366,
                                                                                                                                hDD_dD110)))))),
                                                        FusedMulAddSIMD(FDPart3tmp13,
                                                                        MulSIMD(
                                                                            MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                            MulSIMD(FDPart3tmp42, AddSIMD(FusedMulAddSIMD(FDPart3tmp11, FDPart3tmp374,
                                                                                                                          FDPart3tmp375),
                                                                                                          FusedMulAddSIMD(FDPart3tmp368, FDPart3tmp60,
                                                                                                                          NegFusedMulAddSIMD(
                                                                                                                              MulSIMD(FDPart3_Integer_2,
                                                                                                                                      FDPart3tmp170),
                                                                                                                              MulSIMD(FDPart3tmp256, f1_of_xx1__D1), FDPart3tmp365))))),
                                                                        FusedMulAddSIMD(FDPart3tmp42,
                                                                                        MulSIMD(
                                                                                            MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                            MulSIMD(
                                                                                                FDPart3tmp48,
                                                                                                NegFusedMulAddSIMD(
                                                                                                    MulSIMD(FDPart3_Integer_2, FDPart3tmp7),
                                                                                                    MulSIMD(f1_of_xx1__D1, hDD_dD120),
                                                                                                    FDPart3tmp254))),
                                                                                        FusedMulAddSIMD(FDPart3tmp160,
                                                                                                        MulSIMD(FDPart3tmp211, FDPart3tmp220),
                                                                                                        FusedMulAddSIMD(MulSIMD(FDPart3_Integer_3,
                                                                                                                                FDPart3tmp110),
                                                                                                                        MulSIMD(FDPart3tmp149,
                                                                                                                                FDPart3tmp290),
                                                                                                                        FusedMulAddSIMD(FDPart3tmp143,
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp169,
                                                                                                                                            FDPart3tmp211),
                                                                                                                                        FusedMulAddSIMD(FDPart3tmp158,
                                                                                                                                                        MulSIMD(
                                                                                                                                                            FDPart3tmp220,
                                                                                                                                                            FDPart3tmp290),
                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp135,
                                                                                                                                                                        MulSIMD(
                                                                                                                                                                            FDPart3tmp169,
                                                                                                                                                                            FDPart3tmp290),
                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp135,
                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                            FDPart3tmp220,
                                                                                                                                                                                            FDPart3tmp283),
                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp243,
                                                                                                                                                                                                        FDPart3tmp34,
                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp245,
                                                                                                                                                                                                                        FDPart3tmp74,
                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp236,
                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp221,
                                                                                                                                                                                                                                                        FDPart3tmp231,
                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp220, FDPart3tmp34),
                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                            FDPart3tmp11,
                                                                                                                                                                                                                                                                            FDPart3tmp389))),
                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp237,
                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                            FDPart3tmp11,
                                                                                                                                                                                                                                                                            FDPart3tmp146),
                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp149, FDPart3tmp34),
                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                FDPart3tmp133,
                                                                                                                                                                                                                                                                                FDPart3tmp231))),
                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp23,
                                                                                                                                                                                                                                                                        lambdaU_dD21,
                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp234,
                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp11, FDPart3tmp173),
                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp169, FDPart3tmp34), MulSIMD(FDPart3tmp174, FDPart3tmp231))),
                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                            FDPart3tmp143,
                                                                                                                                                                                                                                                                                            FDPart3tmp285,
                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp162,
                                                                                                                                                                                                                                                                                                            FDPart3tmp226,
                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp125,
                                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                FDPart3_Integer_2, MulSIMD(FDPart3tmp133, FDPart3tmp225),
                                                                                                                                                                                                                                                                                                                                FDPart3tmp157),
                                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp125,
                                                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                FDPart3_Integer_2, MulSIMD(FDPart3tmp146, FDPart3tmp177),
                                                                                                                                                                                                                                                                                                                                                FDPart3tmp206),
                                                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                FDPart3tmp120,
                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp219, FDPart3tmp309, MulSIMD(FDPart3tmp283, FDPart3tmp389)),
                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp124,
                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp378,
                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp118, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp173, FDPart3tmp300), FDPart3tmp384),
                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp120,
                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp211,
                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp387,
                                                                                                                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp221,
                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp222)),
                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp117,
                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp219,
                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp297,
                                                                                                                                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp300,
                                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp389)),
                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp118,
                                                                                                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp174,
                                                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp177),
                                                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp382),
                                                                                                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp117, FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp222, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp174, FDPart3tmp211))),
                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp117,
                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp175,
                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp221, MulSIMD(FDPart3tmp177, FDPart3tmp387)),
                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp110,
                                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp177,
                                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp389,
                                                                                                                                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp201,
                                                                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp219)),
                                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp117, FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp309, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp173, FDPart3tmp283))),
                                                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp110, FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp309, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp146, FDPart3tmp283))),
                                                                                                                                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp110,
                                                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp221, MulSIMD(FDPart3tmp225, FDPart3tmp387)),
                                                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp104,
                                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp173,
                                                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp201,
                                                                                                                                                                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp173,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                       FDPart3tmp177))),
                                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp110,
                                                                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp222, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp133, FDPart3tmp211))),
                                                                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp146, FDPart3tmp300),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp311),
                                                                                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp174, FDPart3tmp225),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp318),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp42,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp93, FusedMulAddSIMD(hDD11,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               FDPart3_Integer_2, FDPart3tmp2, MulSIMD(FDPart3_Integer_2, FDPart3tmp50)),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           NegFusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               FDPart3tmp127,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               MulSIMD(FDPart3tmp5, f0_of_xx0__DD00),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               FusedMulAddSIMD(FDPart3tmp33,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               FDPart3tmp89,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               FusedMulAddSIMD(FDPart3tmp7,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               hDD_dDD1100,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               FusedMulAddSIMD(FDPart3tmp128,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       FDPart3_NegativeOne_),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       FDPart3tmp17,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       f0_of_xx0__D0)),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               FusedMulSubSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   FDPart3tmp105,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   hDD_dD110,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       FDPart3tmp33,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       FDPart3tmp91))))))))),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            FusedMulSubSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp104,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp133, FDPart3tmp177),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                FDPart3tmp302),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3_Rational_1_2,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp13),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    MulSIMD(FDPart3tmp42,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp170),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               MulSIMD(FDPart3tmp248, f1_of_xx1__D1), FDPart3tmp365))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_4 = FusedMulAddSIMD(
    FDPart3tmp170, MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp17), MulSIMD(FDPart3tmp39, lambdaU_dD21)),
    FusedMulAddSIMD(
        MulSIMD(FDPart3tmp17, FDPart3tmp34), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(f1_of_xx1__D1, lambdaU2)),
        FusedMulAddSIMD(
            FDPart3tmp42,
            MulSIMD(
                MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                MulSIMD(FDPart3tmp79,
                        FusedMulAddSIMD(FDPart3tmp372, FDPart3tmp90,
                                        FusedMulAddSIMD(FDPart3tmp113, FDPart3tmp172,
                                                        FusedMulAddSIMD(FDPart3tmp187, FDPart3tmp54,
                                                                        FusedMulAddSIMD(FDPart3tmp7, MulSIMD(f1_of_xx1__D1, hDD_dD120),
                                                                                        NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_3, FDPart3tmp17),
                                                                                                           MulSIMD(FDPart3tmp248, f0_of_xx0__D0),
                                                                                                           FDPart3tmp397))))))),
            FusedMulAddSIMD(
                FDPart3tmp42,
                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                        MulSIMD(FDPart3tmp93,
                                FusedMulAddSIMD(
                                    FDPart3tmp112, hDD_dD120,
                                    FusedMulAddSIMD(
                                        hDD12, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp51, MulSIMD(FDPart3_Integer_2, FDPart3tmp8)),
                                        FusedMulAddSIMD(FDPart3tmp11, FDPart3tmp89,
                                                        FusedMulAddSIMD(FDPart3tmp111, FDPart3tmp137,
                                                                        FusedMulAddSIMD(FDPart3tmp113,
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(FDPart3tmp17, f0_of_xx0__D0)),
                                                                                        FusedMulSubSIMD(FDPart3tmp10, hDD_dDD1200,
                                                                                                        MulSIMD(FDPart3tmp11, FDPart3tmp91))))))))),
                FusedMulAddSIMD(
                    FDPart3tmp42,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp69,
                                    FusedMulAddSIMD(
                                        FDPart3tmp188, FDPart3tmp60,
                                        FusedMulAddSIMD(
                                            FDPart3tmp248, FDPart3tmp73,
                                            FusedMulAddSIMD(FDPart3tmp331, MulSIMD(FDPart3tmp366, hDD_dD120),
                                                            FusedMulSubSIMD(FDPart3tmp172, FDPart3tmp333,
                                                                            MulSIMD(FDPart3tmp7, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp259),
                                                                                                         MulSIMD(f1_of_xx1, hDD12))))))))),
                    FusedMulAddSIMD(
                        FDPart3tmp42,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp79,
                                        FusedMulAddSIMD(FDPart3tmp187, FDPart3tmp55,
                                                        FusedMulAddSIMD(FDPart3tmp248, FDPart3tmp85,
                                                                        NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp17),
                                                                                           MulSIMD(FDPart3tmp256, f0_of_xx0__D0), FDPart3tmp397))))),
                        FusedMulAddSIMD(
                            FDPart3tmp42,
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                    MulSIMD(FDPart3tmp48,
                                            FusedMulAddSIMD(
                                                FDPart3tmp128, FDPart3tmp73,
                                                FusedMulAddSIMD(FDPart3tmp188, FDPart3tmp83,
                                                                FusedMulAddSIMD(FDPart3tmp340, FDPart3tmp78,
                                                                                FusedMulSubSIMD(FDPart3tmp108, FDPart3tmp172,
                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_3, FDPart3tmp17),
                                                                                                        MulSIMD(FDPart3tmp272, f0_of_xx0__D0)))))))),
                            FusedMulAddSIMD(
                                FDPart3tmp42,
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                        MulSIMD(FDPart3tmp68,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp163, MulSIMD(f1_of_xx1__D1, hDD_dD121),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp187, FDPart3tmp25,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp187, FDPart3tmp57,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp172, FDPart3tmp248,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp172, FDPart3tmp256,
                                                                    FusedMulAddSIMD(FDPart3tmp366, MulSIMD(f1_of_xx1, hDD_dD120),
                                                                                    FusedMulAddSIMD(FDPart3tmp7, MulSIMD(f1_of_xx1__DD11, hDD12),
                                                                                                    FusedMulAddSIMD(FDPart3tmp10, hDD_dDD1211,
                                                                                                                    FDPart3tmp375)))))))))),
                                FusedMulAddSIMD(
                                    FDPart3tmp13,
                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                            MulSIMD(FDPart3tmp42,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp374, FDPart3tmp38,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp213, FDPart3tmp73,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp35, hDD_dD011,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp166, FDPart3tmp172,
                                                                    FusedMulAddSIMD(FDPart3tmp187, FDPart3tmp76,
                                                                                    FusedMulAddSIMD(FDPart3tmp163, MulSIMD(FDPart3tmp73, hDD01),
                                                                                                    NegFusedMulAddSIMD(FDPart3tmp33, FDPart3tmp337,
                                                                                                                       FDPart3tmp371))))))))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp42,
                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                MulSIMD(FDPart3tmp48, FusedMulAddSIMD(FDPart3tmp127, FDPart3tmp73,
                                                                                      FusedMulAddSIMD(FDPart3tmp271, hDD_dD220, FDPart3tmp274)))),
                                        FusedMulAddSIMD(
                                            FDPart3_Rational_1_2, MulSIMD(FDPart3tmp11, FDPart3tmp326),
                                            FusedMulAddSIMD(
                                                FDPart3_Rational_1_2, MulSIMD(FDPart3tmp245, FDPart3tmp90),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp237, AddSIMD(FDPart3tmp177, FDPart3tmp201),
                                                    FusedMulAddSIMD(
                                                        FDPart3_Rational_1_2, MulSIMD(FDPart3tmp11, FDPart3tmp243),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp234, AddSIMD(FDPart3tmp297, FDPart3tmp300),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp236, AddSIMD(FDPart3tmp283, FDPart3tmp309),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp125, FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp177, FDPart3tmp346),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp125, FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp225, FDPart3tmp363),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp120,
                                                                            FusedMulAddSIMD(FDPart3tmp175, FDPart3tmp221,
                                                                                            FusedMulAddSIMD(FDPart3tmp201, FDPart3tmp221,
                                                                                                            MulSIMD(FDPart3tmp174, FDPart3tmp211))),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp125, AddSIMD(FDPart3tmp212, FDPart3tmp361),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp120,
                                                                                    FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp283,
                                                                                                    MulSIMD(FDPart3_Integer_2,
                                                                                                            MulSIMD(FDPart3tmp219, FDPart3tmp297))),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp120,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp220, FDPart3tmp283,
                                                                                            FusedMulAddSIMD(FDPart3tmp220, FDPart3tmp309,
                                                                                                            MulSIMD(FDPart3tmp169, FDPart3tmp290))),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp118,
                                                                                            FusedMulAddSIMD(FDPart3tmp169, FDPart3tmp300,
                                                                                                            FusedMulAddSIMD(FDPart3tmp186,
                                                                                                                            FDPart3tmp283,
                                                                                                                            MulSIMD(FDPart3tmp169,
                                                                                                                                    FDPart3tmp297))),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp118,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp174, FDPart3tmp190,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp177, FDPart3tmp189,
                                                                                                        MulSIMD(FDPart3tmp174, FDPart3tmp181))),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp117,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp181, FDPart3tmp221,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp190, FDPart3tmp221,
                                                                                                            MulSIMD(FDPart3tmp174, FDPart3tmp177))),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp118,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp173, FDPart3tmp345,
                                                                                                            MulSIMD(FDPart3tmp185, FDPart3tmp300)),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp117,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp220, FDPart3tmp297,
                                                                                                                FusedMulAddSIMD(FDPart3tmp220,
                                                                                                                                FDPart3tmp300,
                                                                                                                                FDPart3tmp378)),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp117,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp174, FDPart3tmp201,
                                                                                                                    FusedMulAddSIMD(FDPart3tmp189,
                                                                                                                                    FDPart3tmp211,
                                                                                                                                    FDPart3tmp382)),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp117,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp173, FDPart3tmp300,
                                                                                                                        MulSIMD(FDPart3tmp219,
                                                                                                                                FDPart3tmp345)),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp117,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp169,
                                                                                                                            FDPart3tmp309,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp186,
                                                                                                                                FDPart3tmp290,
                                                                                                                                FDPart3tmp378)),
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp110,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp177,
                                                                                                                                FDPart3tmp220,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp201,
                                                                                                                                    FDPart3tmp220,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp169,
                                                                                                                                        FDPart3tmp211))),
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp117,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                    FDPart3tmp384,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp185,
                                                                                                                                        FDPart3tmp283)),
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp110,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp173,
                                                                                                                                        FDPart3tmp177,
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp181,
                                                                                                                                            FDPart3tmp389)),
                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp174,
                                                                                                                                                                    FDPart3tmp225,
                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp204,
                                                                                                                                                                                    FDPart3tmp221,
                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                        FDPart3tmp122,
                                                                                                                                                                                        FDPart3tmp221))),
                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp121,
                                                                                                                                                                                    FDPart3tmp290,
                                                                                                                                                                                    FDPart3tmp362),
                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                                    FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                    FDPart3tmp311,
                                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                                        FDPart3tmp119,
                                                                                                                                                                                                        FDPart3tmp283)),
                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp169, FDPart3tmp201, FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp211, MulSIMD(FDPart3tmp169, FDPart3tmp177))),
                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp115,
                                                                                                                                                                                                                                    FDPart3tmp211,
                                                                                                                                                                                                                                    FDPart3tmp360),
                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp177,
                                                                                                                                                                                                                                                    FDPart3tmp185, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp173, FDPart3tmp181))),
                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                        FDPart3tmp104,
                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                            FDPart3tmp174,
                                                                                                                                                                                                                                            FDPart3tmp204,
                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                FDPart3tmp189,
                                                                                                                                                                                                                                                FDPart3tmp225,
                                                                                                                                                                                                                                                FDPart3tmp294)),
                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp177, FDPart3tmp354),
                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp119,
                                                                                                                                                                                                                                                                                        FDPart3tmp300,
                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                            FDPart3tmp146,
                                                                                                                                                                                                                                                                                            FDPart3tmp345)),
                                                                                                                                                                                                                                                                        FusedMulAddSIMD(MulSIMD(f0_of_xx0, f1_of_xx1),
                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                FDPart3_NegativeOne_,
                                                                                                                                                                                                                                                                                                FDPart3_Rational_1_2),
                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                hDD01,
                                                                                                                                                                                                                                                                                                lambdaU2)),
                                                                                                                                                                                                                                                                                        FusedMulSubSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                                        AddSIMD(
                                                                                                                                                                                                                                                                                                            FDPart3tmp284, FDPart3tmp355),
                                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                                FDPart3_Rational_1_2,
                                                                                                                                                                                                                                                                                                                FDPart3tmp13),
                                                                                                                                                                                                                                                                                                            MulSIMD(FDPart3tmp42,
                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp265,
                                                                                                                                                                                                                                                                                                                                    FDPart3tmp73,
                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp271, hDD_dD221,
                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp172,
                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp272, MulSIMD(FDPart3tmp188, FDPart3tmp66))))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_5 = FusedMulAddSIMD(
    FDPart3tmp42,
    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
            MulSIMD(FDPart3tmp93,
                    FusedMulAddSIMD(
                        FDPart3tmp38, FDPart3tmp89,
                        FusedMulAddSIMD(hDD22, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp28, FDPart3tmp50), FDPart3tmp30),
                                        FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp137,
                                                        FusedMulAddSIMD(FDPart3tmp35, hDD_dDD2200,
                                                                        FusedMulAddSIMD(FDPart3tmp108,
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(FDPart3tmp17, f0_of_xx0__D0)),
                                                                                        FusedMulSubSIMD(FDPart3tmp106, hDD_dD220,
                                                                                                        MulSIMD(FDPart3tmp38, FDPart3tmp91))))))))),
    FusedMulAddSIMD(
        FDPart3tmp42,
        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                MulSIMD(FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp19, MulSIMD(FDPart3tmp28, hDD_dD221),
                                                      NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp166),
                                                                         MulSIMD(FDPart3tmp17, f0_of_xx0__D0), FDPart3tmp411)))),
        FusedMulAddSIMD(
            FDPart3tmp42,
            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                    MulSIMD(FDPart3tmp79,
                            AddSIMD(FDPart3tmp411, FusedMulAddSIMD(FDPart3tmp164, hDD_dD220,
                                                                   NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp108),
                                                                                      MulSIMD(FDPart3tmp170, f1_of_xx1__D1), FDPart3tmp330))))),
            FusedMulAddSIMD(
                FDPart3tmp42,
                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                        MulSIMD(FDPart3tmp68,
                                FusedMulAddSIMD(
                                    hDD22, FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp259, MulSIMD(FDPart3tmp163, FDPart3tmp336)),
                                    FusedMulAddSIMD(
                                        FDPart3tmp28, MulSIMD(FDPart3tmp366, hDD_dD220),
                                        FusedMulAddSIMD(FDPart3tmp35, hDD_dDD2211,
                                                        FusedMulAddSIMD(FDPart3tmp374, FDPart3tmp38,
                                                                        FusedMulAddSIMD(FDPart3tmp166,
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(FDPart3tmp170, f1_of_xx1__D1)),
                                                                                        FusedMulSubSIMD(FDPart3tmp164, hDD_dD221,
                                                                                                        MulSIMD(FDPart3tmp261, FDPart3tmp38))))))))),
                FusedMulAddSIMD(
                    FDPart3tmp42,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp69,
                                    FusedMulAddSIMD(
                                        FDPart3tmp366, MulSIMD(hDD_dD220, MulSIMD(MulSIMD(MulSIMD(f1_of_xx1, f1_of_xx1), f1_of_xx1), f1_of_xx1)),
                                        FusedMulAddSIMD(MulSIMD(FDPart3tmp331, FDPart3tmp7), MulSIMD(f1_of_xx1__D1, hDD_dD221),
                                                        FusedMulAddSIMD(FDPart3tmp272, FDPart3tmp414, MulSIMD(FDPart3tmp412, FDPart3tmp76)))))),
                    FusedMulAddSIMD(
                        FDPart3tmp42,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp48,
                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp188, FDPart3tmp55),
                                                        FusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp35), MulSIMD(f1_of_xx1__D1, hDD_dD120),
                                                                        MulSIMD(FDPart3tmp333, FDPart3tmp85))))),
                        FusedMulAddSIMD(
                            FDPart3tmp42,
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                    MulSIMD(FDPart3tmp48,
                                            FusedMulAddSIMD(
                                                FDPart3tmp340, FDPart3tmp47,
                                                FusedMulAddSIMD(FDPart3tmp412, FDPart3tmp54,
                                                                FusedMulAddSIMD(FDPart3tmp90,
                                                                                NegFusedMulAddSIMD(FDPart3tmp276, FDPart3tmp28, FDPart3tmp28),
                                                                                FusedMulSubSIMD(FDPart3tmp113, FDPart3tmp414,
                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_3, FDPart3tmp17),
                                                                                                        MulSIMD(FDPart3tmp333, f0_of_xx0__D0)))))))),
                            FusedMulAddSIMD(
                                FDPart3tmp13,
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                        MulSIMD(FDPart3tmp42, FusedMulAddSIMD(FDPart3tmp25, FDPart3tmp412,
                                                                              FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp248, FDPart3tmp73),
                                                                                              MulSIMD(FDPart3tmp172, FDPart3tmp333))))),
                                FusedMulAddSIMD(
                                    FDPart3tmp13,
                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                            MulSIMD(FDPart3tmp42,
                                                    FusedMulAddSIMD(
                                                        MulSIMD(FDPart3_Integer_4, FDPart3tmp35), MulSIMD(f1_of_xx1__D1, hDD02),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp256, FDPart3tmp414,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp412, FDPart3tmp57,
                                                                FusedMulAddSIMD(FDPart3tmp170,
                                                                                MulSIMD(MulSIMD(FDPart3_Integer_3, FDPart3_NegativeOne_),
                                                                                        MulSIMD(FDPart3tmp333, f1_of_xx1__D1)),
                                                                                FusedMulSubSIMD(FDPart3tmp11, AddSIMD(FDPart3tmp259, FDPart3tmp336),
                                                                                                MulSIMD(FDPart3tmp11, FDPart3tmp337)))))))),
                                    FusedMulAddSIMD(
                                        MulSIMD(FDPart3_Integer_3, FDPart3tmp104), MulSIMD(FDPart3tmp119, FDPart3tmp344),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3_Integer_3, FDPart3tmp117), MulSIMD(FDPart3tmp173, FDPart3tmp344),
                                            FusedMulAddSIMD(
                                                FDPart3tmp143, MulSIMD(FDPart3tmp181, FDPart3tmp185),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp160, MulSIMD(FDPart3tmp173, FDPart3tmp181),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp124, MulSIMD(FDPart3tmp185, FDPart3tmp344),
                                                        FusedMulAddSIMD(FDPart3tmp135, MulSIMD(FDPart3tmp185, FDPart3tmp297),
                                                                        FusedMulAddSIMD(FDPart3tmp237,
                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp11, FDPart3tmp121),
                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                        MulSIMD(FDPart3tmp119,
                                                                                                                                FDPart3tmp39),
                                                                                                                        MulSIMD(FDPart3tmp115,
                                                                                                                                FDPart3tmp230))),
                                                                                        FusedMulAddSIMD(FDPart3tmp326, FDPart3tmp39,
                                                                                                        FusedMulAddSIMD(FDPart3tmp234,
                                                                                                                        FusedMulAddSIMD(FDPart3tmp189,
                                                                                                                                        FDPart3tmp230,
                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                        MulSIMD(
                                                                                                                                                            FDPart3tmp185,
                                                                                                                                                            FDPart3tmp39),
                                                                                                                                                        MulSIMD(
                                                                                                                                                            FDPart3tmp11,
                                                                                                                                                            FDPart3tmp417))),
                                                                                                                        FusedMulAddSIMD(FDPart3tmp236,
                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                        MulSIMD(
                                                                                                                                                            FDPart3tmp11,
                                                                                                                                                            FDPart3tmp169),
                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                        MulSIMD(
                                                                                                                                                                            FDPart3tmp173,
                                                                                                                                                                            FDPart3tmp39),
                                                                                                                                                                        MulSIMD(
                                                                                                                                                                            FDPart3tmp174,
                                                                                                                                                                            FDPart3tmp230))),
                                                                                                                                        FusedMulAddSIMD(FDPart3tmp160,
                                                                                                                                                        FDPart3tmp298,
                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp162,
                                                                                                                                                                        FDPart3tmp182,
                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp125,
                                                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                            FDPart3tmp121,
                                                                                                                                                                                                            FDPart3tmp201),
                                                                                                                                                                                                        FDPart3tmp178),
                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp158,
                                                                                                                                                                                                        FDPart3tmp384,
                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp120,
                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp174, FDPart3tmp201), FDPart3tmp382),
                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp125,
                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                            FDPart3tmp115,
                                                                                                                                                                                                                                                            FDPart3tmp204),
                                                                                                                                                                                                                                                        FDPart3tmp123),
                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp118,
                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp189,
                                                                                                                                                                                                                                                                        FDPart3tmp190,
                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                            FDPart3_Integer_2, MulSIMD(
                                                                                                                                                                                                                                                                                                   FDPart3tmp181,
                                                                                                                                                                                                                                                                                                   FDPart3tmp189))),
                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp120,
                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                            FDPart3tmp169,
                                                                                                                                                                                                                                                                                            FDPart3tmp309),
                                                                                                                                                                                                                                                                                        FDPart3tmp378),
                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp117, FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp283, MulSIMD(FDPart3tmp309, FDPart3tmp417)),
                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                            FDPart3tmp118, FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp300, MulSIMD(FDPart3tmp297, FDPart3tmp417)),
                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                FDPart3tmp117,
                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                    FDPart3tmp174,
                                                                                                                                                                                                                                                                                                    FDPart3tmp190,
                                                                                                                                                                                                                                                                                                    MulSIMD(FDPart3_Integer_2, MulSIMD(
                                                                                                                                                                                                                                                                                                                                   FDPart3tmp174,
                                                                                                                                                                                                                                                                                                                                   FDPart3tmp181))),
                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                    FDPart3tmp117,
                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                        FDPart3tmp175,
                                                                                                                                                                                                                                                                                                        FDPart3tmp189,
                                                                                                                                                                                                                                                                                                        MulSIMD(FDPart3_Integer_2, MulSIMD(
                                                                                                                                                                                                                                                                                                                                       FDPart3tmp189, FDPart3tmp201))),
                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp169,
                                                                                                                                                                                                                                                                                                                                    FDPart3tmp177, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp169, FDPart3tmp201))),
                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp117, FusedMulAddSIMD(FDPart3tmp169, FDPart3tmp300, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp169, FDPart3tmp297))),
                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp115,
                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp201),
                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp292),
                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp174, FDPart3tmp204), FDPart3tmp294),
                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp177, FDPart3tmp186, MulSIMD(FDPart3tmp201, FDPart3tmp417)),
                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp110, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp121, FDPart3tmp309), FDPart3tmp284),
                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp121, FDPart3tmp300, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp121, FDPart3tmp297))),
                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp104,
                                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp189, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp189, FDPart3tmp204))),
                                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp28), MulSIMD(hDD02, lambdaU2)),
                                                                                                                                                                                                                                                                                                                                                                                                                                                    FusedMulSubSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp190, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp115, FDPart3tmp181))),
                                                                                                                                                                                                                                                                                                                                                                                                                                                        MulSIMD(FDPart3tmp23, MulSIMD(FDPart3tmp73,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      lambdaU2)))))))))))))))))))))))))))))))))))))))))));

WriteSIMD(&auxevol_gfs[IDX4(RBARDD00GF, i0, i1, i2)], __RHS_exp_0);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD01GF, i0, i1, i2)], __RHS_exp_1);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD02GF, i0, i1, i2)], __RHS_exp_2);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD11GF, i0, i1, i2)], __RHS_exp_3);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD12GF, i0, i1, i2)], __RHS_exp_4);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD22GF, i0, i1, i2)], __RHS_exp_5);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += simd_width)
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
