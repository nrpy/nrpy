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
 * Finite difference function for operator dKOD0, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_dKOD0_fdorder8(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i0m1,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0m2, const REAL_SIMD_ARRAY FDPROTO_i0m3,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0m4, const REAL_SIMD_ARRAY FDPROTO_i0m5,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0p1, const REAL_SIMD_ARRAY FDPROTO_i0p2,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0p3, const REAL_SIMD_ARRAY FDPROTO_i0p4,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0p5, const REAL_SIMD_ARRAY invdxx0) {
  const double dblFDPart1_Rational_105_512 = 105.0 / 512.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_105_512 = ConstSIMD(dblFDPart1_Rational_105_512);

  const double dblFDPart1_Rational_15_128 = 15.0 / 128.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_15_128 = ConstSIMD(dblFDPart1_Rational_15_128);

  const double dblFDPart1_Rational_1_1024 = 1.0 / 1024.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_1024 = ConstSIMD(dblFDPart1_Rational_1_1024);

  const double dblFDPart1_Rational_45_1024 = 45.0 / 1024.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_45_1024 = ConstSIMD(dblFDPart1_Rational_45_1024);

  const double dblFDPart1_Rational_5_512 = 5.0 / 512.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_512 = ConstSIMD(dblFDPart1_Rational_5_512);

  const double dblFDPart1_Rational_63_256 = 63.0 / 256.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_63_256 = ConstSIMD(dblFDPart1_Rational_63_256);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx0,
      FusedMulAddSIMD(FDPart1_Rational_1_1024, AddSIMD(FDPROTO_i0m5, FDPROTO_i0p5),
                      FusedMulAddSIMD(FDPart1_Rational_45_1024, AddSIMD(FDPROTO_i0m3, FDPROTO_i0p3),
                                      FusedMulSubSIMD(FDPart1_Rational_105_512, AddSIMD(FDPROTO_i0m1, FDPROTO_i0p1),
                                                      FusedMulAddSIMD(FDPart1_Rational_15_128, AddSIMD(FDPROTO_i0m2, FDPROTO_i0p2),
                                                                      FusedMulAddSIMD(FDPart1_Rational_5_512, AddSIMD(FDPROTO_i0m4, FDPROTO_i0p4),
                                                                                      MulSIMD(FDPROTO, FDPart1_Rational_63_256)))))));

  return FD_result;
}
/*
 * Finite difference function for operator dKOD1, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_dKOD1_fdorder8(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i1m1,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1m2, const REAL_SIMD_ARRAY FDPROTO_i1m3,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1m4, const REAL_SIMD_ARRAY FDPROTO_i1m5,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1p1, const REAL_SIMD_ARRAY FDPROTO_i1p2,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1p3, const REAL_SIMD_ARRAY FDPROTO_i1p4,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1p5, const REAL_SIMD_ARRAY invdxx1) {
  const double dblFDPart1_Rational_105_512 = 105.0 / 512.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_105_512 = ConstSIMD(dblFDPart1_Rational_105_512);

  const double dblFDPart1_Rational_15_128 = 15.0 / 128.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_15_128 = ConstSIMD(dblFDPart1_Rational_15_128);

  const double dblFDPart1_Rational_1_1024 = 1.0 / 1024.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_1024 = ConstSIMD(dblFDPart1_Rational_1_1024);

  const double dblFDPart1_Rational_45_1024 = 45.0 / 1024.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_45_1024 = ConstSIMD(dblFDPart1_Rational_45_1024);

  const double dblFDPart1_Rational_5_512 = 5.0 / 512.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_512 = ConstSIMD(dblFDPart1_Rational_5_512);

  const double dblFDPart1_Rational_63_256 = 63.0 / 256.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_63_256 = ConstSIMD(dblFDPart1_Rational_63_256);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx1,
      FusedMulAddSIMD(FDPart1_Rational_1_1024, AddSIMD(FDPROTO_i1m5, FDPROTO_i1p5),
                      FusedMulAddSIMD(FDPart1_Rational_45_1024, AddSIMD(FDPROTO_i1m3, FDPROTO_i1p3),
                                      FusedMulSubSIMD(FDPart1_Rational_105_512, AddSIMD(FDPROTO_i1m1, FDPROTO_i1p1),
                                                      FusedMulAddSIMD(FDPart1_Rational_15_128, AddSIMD(FDPROTO_i1m2, FDPROTO_i1p2),
                                                                      FusedMulAddSIMD(FDPart1_Rational_5_512, AddSIMD(FDPROTO_i1m4, FDPROTO_i1p4),
                                                                                      MulSIMD(FDPROTO, FDPart1_Rational_63_256)))))));

  return FD_result;
}
/*
 * Finite difference function for operator ddnD0, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_ddnD0_fdorder8(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i0m1,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0m2, const REAL_SIMD_ARRAY FDPROTO_i0m3,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0m4, const REAL_SIMD_ARRAY FDPROTO_i0m5,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0p1, const REAL_SIMD_ARRAY FDPROTO_i0p2,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0p3, const REAL_SIMD_ARRAY invdxx0) {
  const double dblFDPart1_Rational_1_14 = 1.0 / 14.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_14 = ConstSIMD(dblFDPart1_Rational_1_14);

  const double dblFDPart1_Rational_1_168 = 1.0 / 168.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_168 = ConstSIMD(dblFDPart1_Rational_1_168);

  const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  const double dblFDPart1_Rational_1_28 = 1.0 / 28.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_28 = ConstSIMD(dblFDPart1_Rational_1_28);

  const double dblFDPart1_Rational_1_280 = 1.0 / 280.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_280 = ConstSIMD(dblFDPart1_Rational_1_280);

  const double dblFDPart1_Rational_1_6 = 1.0 / 6.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_6 = ConstSIMD(dblFDPart1_Rational_1_6);

  const double dblFDPart1_Rational_5_4 = 5.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_4 = ConstSIMD(dblFDPart1_Rational_5_4);

  const double dblFDPart1_Rational_9_20 = 9.0 / 20.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_9_20 = ConstSIMD(dblFDPart1_Rational_9_20);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx0,
      FusedMulAddSIMD(
          FDPROTO_i0m4, FDPart1_Rational_1_28,
          FusedMulAddSIMD(
              FDPROTO_i0p3, FDPart1_Rational_1_168,
              FusedMulAddSIMD(FDPart1_Rational_1_2, AddSIMD(FDPROTO_i0m2, FDPROTO_i0p1),
                              FusedMulSubSIMD(FDPROTO, FDPart1_Rational_9_20,
                                              FusedMulAddSIMD(FDPROTO_i0m5, FDPart1_Rational_1_280,
                                                              FusedMulAddSIMD(FDPROTO_i0p2, FDPart1_Rational_1_14,
                                                                              FusedMulAddSIMD(FDPROTO_i0m1, FDPart1_Rational_5_4,
                                                                                              MulSIMD(FDPROTO_i0m3, FDPart1_Rational_1_6)))))))));

  return FD_result;
}
/*
 * Finite difference function for operator ddnD1, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_ddnD1_fdorder8(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i1m1,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1m2, const REAL_SIMD_ARRAY FDPROTO_i1m3,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1m4, const REAL_SIMD_ARRAY FDPROTO_i1m5,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1p1, const REAL_SIMD_ARRAY FDPROTO_i1p2,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1p3, const REAL_SIMD_ARRAY invdxx1) {
  const double dblFDPart1_Rational_1_14 = 1.0 / 14.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_14 = ConstSIMD(dblFDPart1_Rational_1_14);

  const double dblFDPart1_Rational_1_168 = 1.0 / 168.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_168 = ConstSIMD(dblFDPart1_Rational_1_168);

  const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  const double dblFDPart1_Rational_1_28 = 1.0 / 28.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_28 = ConstSIMD(dblFDPart1_Rational_1_28);

  const double dblFDPart1_Rational_1_280 = 1.0 / 280.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_280 = ConstSIMD(dblFDPart1_Rational_1_280);

  const double dblFDPart1_Rational_1_6 = 1.0 / 6.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_6 = ConstSIMD(dblFDPart1_Rational_1_6);

  const double dblFDPart1_Rational_5_4 = 5.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_4 = ConstSIMD(dblFDPart1_Rational_5_4);

  const double dblFDPart1_Rational_9_20 = 9.0 / 20.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_9_20 = ConstSIMD(dblFDPart1_Rational_9_20);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx1,
      FusedMulAddSIMD(
          FDPROTO_i1m4, FDPart1_Rational_1_28,
          FusedMulAddSIMD(
              FDPROTO_i1p3, FDPart1_Rational_1_168,
              FusedMulAddSIMD(FDPart1_Rational_1_2, AddSIMD(FDPROTO_i1m2, FDPROTO_i1p1),
                              FusedMulSubSIMD(FDPROTO, FDPart1_Rational_9_20,
                                              FusedMulAddSIMD(FDPROTO_i1m5, FDPart1_Rational_1_280,
                                                              FusedMulAddSIMD(FDPROTO_i1p2, FDPart1_Rational_1_14,
                                                                              FusedMulAddSIMD(FDPROTO_i1m1, FDPart1_Rational_5_4,
                                                                                              MulSIMD(FDPROTO_i1m3, FDPart1_Rational_1_6)))))))));

  return FD_result;
}
/*
 * Finite difference function for operator dupD0, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_dupD0_fdorder8(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i0m1,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0m2, const REAL_SIMD_ARRAY FDPROTO_i0m3,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0p1, const REAL_SIMD_ARRAY FDPROTO_i0p2,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0p3, const REAL_SIMD_ARRAY FDPROTO_i0p4,
                                                       const REAL_SIMD_ARRAY FDPROTO_i0p5, const REAL_SIMD_ARRAY invdxx0) {
  const double dblFDPart1_Rational_1_14 = 1.0 / 14.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_14 = ConstSIMD(dblFDPart1_Rational_1_14);

  const double dblFDPart1_Rational_1_168 = 1.0 / 168.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_168 = ConstSIMD(dblFDPart1_Rational_1_168);

  const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  const double dblFDPart1_Rational_1_28 = 1.0 / 28.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_28 = ConstSIMD(dblFDPart1_Rational_1_28);

  const double dblFDPart1_Rational_1_280 = 1.0 / 280.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_280 = ConstSIMD(dblFDPart1_Rational_1_280);

  const double dblFDPart1_Rational_1_6 = 1.0 / 6.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_6 = ConstSIMD(dblFDPart1_Rational_1_6);

  const double dblFDPart1_Rational_5_4 = 5.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_4 = ConstSIMD(dblFDPart1_Rational_5_4);

  const double dblFDPart1_Rational_9_20 = 9.0 / 20.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_9_20 = ConstSIMD(dblFDPart1_Rational_9_20);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx0,
      FusedMulAddSIMD(
          FDPROTO_i0p1, FDPart1_Rational_5_4,
          FusedMulAddSIMD(
              FDPROTO_i0p3, FDPart1_Rational_1_6,
              FusedMulAddSIMD(FDPROTO_i0p5, FDPart1_Rational_1_280,
                              FusedMulSubSIMD(FDPROTO_i0m2, FDPart1_Rational_1_14,
                                              FusedMulAddSIMD(FDPROTO_i0p4, FDPart1_Rational_1_28,
                                                              FusedMulAddSIMD(FDPart1_Rational_1_2, AddSIMD(FDPROTO_i0m1, FDPROTO_i0p2),
                                                                              FusedMulAddSIMD(FDPROTO, FDPart1_Rational_9_20,
                                                                                              MulSIMD(FDPROTO_i0m3, FDPart1_Rational_1_168)))))))));

  return FD_result;
}
/*
 * Finite difference function for operator dupD1, with FD accuracy order 8.
 */
static REAL_SIMD_ARRAY SIMD_fd_function_dupD1_fdorder8(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i1m1,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1m2, const REAL_SIMD_ARRAY FDPROTO_i1m3,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1p1, const REAL_SIMD_ARRAY FDPROTO_i1p2,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1p3, const REAL_SIMD_ARRAY FDPROTO_i1p4,
                                                       const REAL_SIMD_ARRAY FDPROTO_i1p5, const REAL_SIMD_ARRAY invdxx1) {
  const double dblFDPart1_Rational_1_14 = 1.0 / 14.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_14 = ConstSIMD(dblFDPart1_Rational_1_14);

  const double dblFDPart1_Rational_1_168 = 1.0 / 168.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_168 = ConstSIMD(dblFDPart1_Rational_1_168);

  const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  const double dblFDPart1_Rational_1_28 = 1.0 / 28.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_28 = ConstSIMD(dblFDPart1_Rational_1_28);

  const double dblFDPart1_Rational_1_280 = 1.0 / 280.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_280 = ConstSIMD(dblFDPart1_Rational_1_280);

  const double dblFDPart1_Rational_1_6 = 1.0 / 6.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_6 = ConstSIMD(dblFDPart1_Rational_1_6);

  const double dblFDPart1_Rational_5_4 = 5.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_4 = ConstSIMD(dblFDPart1_Rational_5_4);

  const double dblFDPart1_Rational_9_20 = 9.0 / 20.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_9_20 = ConstSIMD(dblFDPart1_Rational_9_20);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx1,
      FusedMulAddSIMD(
          FDPROTO_i1p1, FDPart1_Rational_5_4,
          FusedMulAddSIMD(
              FDPROTO_i1p3, FDPart1_Rational_1_6,
              FusedMulAddSIMD(FDPROTO_i1p5, FDPart1_Rational_1_280,
                              FusedMulSubSIMD(FDPROTO_i1m2, FDPart1_Rational_1_14,
                                              FusedMulAddSIMD(FDPROTO_i1p4, FDPart1_Rational_1_28,
                                                              FusedMulAddSIMD(FDPart1_Rational_1_2, AddSIMD(FDPROTO_i1m1, FDPROTO_i1p2),
                                                                              FusedMulAddSIMD(FDPROTO, FDPart1_Rational_9_20,
                                                                                              MulSIMD(FDPROTO_i1m3, FDPart1_Rational_1_168)))))))));

  return FD_result;
}

/*
 * Set RHSs for the BSSN evolution equations.
 */
void rhs_eval__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                  const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                                  REAL *restrict rhs_gfs) {
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
         * NRPy+-Generated GF Access/FD Code, Step 1 of 3:
         * Read gridfunction(s) from main memory and compute FD stencils as needed.
         */
        const REAL_SIMD_ARRAY RbarDD00 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD01 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD02 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD11 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD12 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD22 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m5 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m5 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p5 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p5 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m5 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m5 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p5 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p5 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m5 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m5 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p5 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p5 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m5 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m5 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p5 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p5 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m5 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m5 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p5 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p5 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m5 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m4 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m5 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m4 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p4 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p5 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p4 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p5 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY alpha_i1m5 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i1m4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i1m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m5 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY alpha = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p5 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i1p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i1p4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY alpha_i1p5 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY betU0_i1m5 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY betU0_i1m4 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY betU0_i1m3 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY betU0_i1m2 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY betU0_i1m1 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0m5 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0m4 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0m3 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0m2 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0m1 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU0 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0p1 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0p2 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0p3 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0p4 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0p5 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i1p1 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY betU0_i1p2 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY betU0_i1p3 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY betU0_i1p4 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY betU0_i1p5 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY betU1_i1m5 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY betU1_i1m4 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY betU1_i1m3 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY betU1_i1m2 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY betU1_i1m1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0m5 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0m4 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0m3 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0m2 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0m1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0p1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0p2 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0p3 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0p4 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0p5 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i1p1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY betU1_i1p2 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY betU1_i1p3 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY betU1_i1p4 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY betU1_i1p5 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY betU2_i1m5 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY betU2_i1m4 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY betU2_i1m3 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY betU2_i1m2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY betU2_i1m1 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0m5 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0m4 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0m3 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0m2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0m1 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0p1 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0p2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0p3 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0p4 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0p5 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i1p1 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY betU2_i1p2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY betU2_i1p3 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY betU2_i1p4 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY betU2_i1p5 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY cf_i1m5 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 5, i2)]);
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
        const REAL_SIMD_ARRAY cf_i0m5 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY cf = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p4 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p5 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 5, i1, i2)]);
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
        const REAL_SIMD_ARRAY cf_i1p5 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m5 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m5 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p5 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p5 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m5 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m5 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p5 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p5 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m5 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m5 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p5 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p5 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m5 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m5 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p5 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p5 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m5 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m5 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p5 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p5 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m5 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m5 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p5 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p4 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p5 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p4 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p5 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY trK_i1m5 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY trK_i1m4 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY trK_i1m3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY trK_i1m2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY trK_i1m1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m5 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m4 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY trK = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p4 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p5 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i1p1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY trK_i1p2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY trK_i1p3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY trK_i1p4 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY trK_i1p5 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1m5 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m5 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p5 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1p5 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1m5 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m5 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p5 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1p5 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1m5 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 5, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m4_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m3_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p3_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 3, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p4_i1m4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 4, i1 - 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m4_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m3_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p3_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 3, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p4_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 4, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m4_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m3_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p3_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 3, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p4_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 4, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m4_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m3_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p3_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 3, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p4_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 4, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m5 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 5, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 4, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 4, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p5 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 5, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m4_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m3_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p3_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 3, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p4_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 4, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m4_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m3_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p3_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 3, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p4_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 4, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m4_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m3_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p3_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 3, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p4_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 4, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m4_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m3_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p3_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 3, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p4_i1p4 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 4, i1 + 4, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1p5 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 5, i2)]);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD000 = SIMD_fd_function_ddnD0_fdorder8(aDD00, aDD00_i0m1, aDD00_i0m2, aDD00_i0m3, aDD00_i0m4,
                                                                                          aDD00_i0m5, aDD00_i0p1, aDD00_i0p2, aDD00_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD001 = SIMD_fd_function_ddnD1_fdorder8(aDD00, aDD00_i1m1, aDD00_i1m2, aDD00_i1m3, aDD00_i1m4,
                                                                                          aDD00_i1m5, aDD00_i1p1, aDD00_i1p2, aDD00_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD010 = SIMD_fd_function_ddnD0_fdorder8(aDD01, aDD01_i0m1, aDD01_i0m2, aDD01_i0m3, aDD01_i0m4,
                                                                                          aDD01_i0m5, aDD01_i0p1, aDD01_i0p2, aDD01_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD011 = SIMD_fd_function_ddnD1_fdorder8(aDD01, aDD01_i1m1, aDD01_i1m2, aDD01_i1m3, aDD01_i1m4,
                                                                                          aDD01_i1m5, aDD01_i1p1, aDD01_i1p2, aDD01_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD020 = SIMD_fd_function_ddnD0_fdorder8(aDD02, aDD02_i0m1, aDD02_i0m2, aDD02_i0m3, aDD02_i0m4,
                                                                                          aDD02_i0m5, aDD02_i0p1, aDD02_i0p2, aDD02_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD021 = SIMD_fd_function_ddnD1_fdorder8(aDD02, aDD02_i1m1, aDD02_i1m2, aDD02_i1m3, aDD02_i1m4,
                                                                                          aDD02_i1m5, aDD02_i1p1, aDD02_i1p2, aDD02_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD110 = SIMD_fd_function_ddnD0_fdorder8(aDD11, aDD11_i0m1, aDD11_i0m2, aDD11_i0m3, aDD11_i0m4,
                                                                                          aDD11_i0m5, aDD11_i0p1, aDD11_i0p2, aDD11_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD111 = SIMD_fd_function_ddnD1_fdorder8(aDD11, aDD11_i1m1, aDD11_i1m2, aDD11_i1m3, aDD11_i1m4,
                                                                                          aDD11_i1m5, aDD11_i1p1, aDD11_i1p2, aDD11_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD120 = SIMD_fd_function_ddnD0_fdorder8(aDD12, aDD12_i0m1, aDD12_i0m2, aDD12_i0m3, aDD12_i0m4,
                                                                                          aDD12_i0m5, aDD12_i0p1, aDD12_i0p2, aDD12_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD121 = SIMD_fd_function_ddnD1_fdorder8(aDD12, aDD12_i1m1, aDD12_i1m2, aDD12_i1m3, aDD12_i1m4,
                                                                                          aDD12_i1m5, aDD12_i1p1, aDD12_i1p2, aDD12_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD220 = SIMD_fd_function_ddnD0_fdorder8(aDD22, aDD22_i0m1, aDD22_i0m2, aDD22_i0m3, aDD22_i0m4,
                                                                                          aDD22_i0m5, aDD22_i0p1, aDD22_i0p2, aDD22_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD221 = SIMD_fd_function_ddnD1_fdorder8(aDD22, aDD22_i1m1, aDD22_i1m2, aDD22_i1m3, aDD22_i1m4,
                                                                                          aDD22_i1m5, aDD22_i1p1, aDD22_i1p2, aDD22_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD000 = SIMD_fd_function_dupD0_fdorder8(aDD00, aDD00_i0m1, aDD00_i0m2, aDD00_i0m3, aDD00_i0p1,
                                                                                          aDD00_i0p2, aDD00_i0p3, aDD00_i0p4, aDD00_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD001 = SIMD_fd_function_dupD1_fdorder8(aDD00, aDD00_i1m1, aDD00_i1m2, aDD00_i1m3, aDD00_i1p1,
                                                                                          aDD00_i1p2, aDD00_i1p3, aDD00_i1p4, aDD00_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD010 = SIMD_fd_function_dupD0_fdorder8(aDD01, aDD01_i0m1, aDD01_i0m2, aDD01_i0m3, aDD01_i0p1,
                                                                                          aDD01_i0p2, aDD01_i0p3, aDD01_i0p4, aDD01_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD011 = SIMD_fd_function_dupD1_fdorder8(aDD01, aDD01_i1m1, aDD01_i1m2, aDD01_i1m3, aDD01_i1p1,
                                                                                          aDD01_i1p2, aDD01_i1p3, aDD01_i1p4, aDD01_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD020 = SIMD_fd_function_dupD0_fdorder8(aDD02, aDD02_i0m1, aDD02_i0m2, aDD02_i0m3, aDD02_i0p1,
                                                                                          aDD02_i0p2, aDD02_i0p3, aDD02_i0p4, aDD02_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD021 = SIMD_fd_function_dupD1_fdorder8(aDD02, aDD02_i1m1, aDD02_i1m2, aDD02_i1m3, aDD02_i1p1,
                                                                                          aDD02_i1p2, aDD02_i1p3, aDD02_i1p4, aDD02_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD110 = SIMD_fd_function_dupD0_fdorder8(aDD11, aDD11_i0m1, aDD11_i0m2, aDD11_i0m3, aDD11_i0p1,
                                                                                          aDD11_i0p2, aDD11_i0p3, aDD11_i0p4, aDD11_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD111 = SIMD_fd_function_dupD1_fdorder8(aDD11, aDD11_i1m1, aDD11_i1m2, aDD11_i1m3, aDD11_i1p1,
                                                                                          aDD11_i1p2, aDD11_i1p3, aDD11_i1p4, aDD11_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD120 = SIMD_fd_function_dupD0_fdorder8(aDD12, aDD12_i0m1, aDD12_i0m2, aDD12_i0m3, aDD12_i0p1,
                                                                                          aDD12_i0p2, aDD12_i0p3, aDD12_i0p4, aDD12_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD121 = SIMD_fd_function_dupD1_fdorder8(aDD12, aDD12_i1m1, aDD12_i1m2, aDD12_i1m3, aDD12_i1p1,
                                                                                          aDD12_i1p2, aDD12_i1p3, aDD12_i1p4, aDD12_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD220 = SIMD_fd_function_dupD0_fdorder8(aDD22, aDD22_i0m1, aDD22_i0m2, aDD22_i0m3, aDD22_i0p1,
                                                                                          aDD22_i0p2, aDD22_i0p3, aDD22_i0p4, aDD22_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD221 = SIMD_fd_function_dupD1_fdorder8(aDD22, aDD22_i1m1, aDD22_i1m2, aDD22_i1m3, aDD22_i1p1,
                                                                                          aDD22_i1p2, aDD22_i1p3, aDD22_i1p4, aDD22_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputalpha_ddnD0 = SIMD_fd_function_ddnD0_fdorder8(alpha, alpha_i0m1, alpha_i0m2, alpha_i0m3, alpha_i0m4,
                                                                                          alpha_i0m5, alpha_i0p1, alpha_i0p2, alpha_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputalpha_ddnD1 = SIMD_fd_function_ddnD1_fdorder8(alpha, alpha_i1m1, alpha_i1m2, alpha_i1m3, alpha_i1m4,
                                                                                          alpha_i1m5, alpha_i1p1, alpha_i1p2, alpha_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputalpha_dupD0 = SIMD_fd_function_dupD0_fdorder8(alpha, alpha_i0m1, alpha_i0m2, alpha_i0m3, alpha_i0p1,
                                                                                          alpha_i0p2, alpha_i0p3, alpha_i0p4, alpha_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputalpha_dupD1 = SIMD_fd_function_dupD1_fdorder8(alpha, alpha_i1m1, alpha_i1m2, alpha_i1m3, alpha_i1p1,
                                                                                          alpha_i1p2, alpha_i1p3, alpha_i1p4, alpha_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD00 = SIMD_fd_function_ddnD0_fdorder8(betU0, betU0_i0m1, betU0_i0m2, betU0_i0m3, betU0_i0m4,
                                                                                          betU0_i0m5, betU0_i0p1, betU0_i0p2, betU0_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD01 = SIMD_fd_function_ddnD1_fdorder8(betU0, betU0_i1m1, betU0_i1m2, betU0_i1m3, betU0_i1m4,
                                                                                          betU0_i1m5, betU0_i1p1, betU0_i1p2, betU0_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD10 = SIMD_fd_function_ddnD0_fdorder8(betU1, betU1_i0m1, betU1_i0m2, betU1_i0m3, betU1_i0m4,
                                                                                          betU1_i0m5, betU1_i0p1, betU1_i0p2, betU1_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD11 = SIMD_fd_function_ddnD1_fdorder8(betU1, betU1_i1m1, betU1_i1m2, betU1_i1m3, betU1_i1m4,
                                                                                          betU1_i1m5, betU1_i1p1, betU1_i1p2, betU1_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD20 = SIMD_fd_function_ddnD0_fdorder8(betU2, betU2_i0m1, betU2_i0m2, betU2_i0m3, betU2_i0m4,
                                                                                          betU2_i0m5, betU2_i0p1, betU2_i0p2, betU2_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD21 = SIMD_fd_function_ddnD1_fdorder8(betU2, betU2_i1m1, betU2_i1m2, betU2_i1m3, betU2_i1m4,
                                                                                          betU2_i1m5, betU2_i1p1, betU2_i1p2, betU2_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD00 = SIMD_fd_function_dupD0_fdorder8(betU0, betU0_i0m1, betU0_i0m2, betU0_i0m3, betU0_i0p1,
                                                                                          betU0_i0p2, betU0_i0p3, betU0_i0p4, betU0_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD01 = SIMD_fd_function_dupD1_fdorder8(betU0, betU0_i1m1, betU0_i1m2, betU0_i1m3, betU0_i1p1,
                                                                                          betU0_i1p2, betU0_i1p3, betU0_i1p4, betU0_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD10 = SIMD_fd_function_dupD0_fdorder8(betU1, betU1_i0m1, betU1_i0m2, betU1_i0m3, betU1_i0p1,
                                                                                          betU1_i0p2, betU1_i0p3, betU1_i0p4, betU1_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD11 = SIMD_fd_function_dupD1_fdorder8(betU1, betU1_i1m1, betU1_i1m2, betU1_i1m3, betU1_i1p1,
                                                                                          betU1_i1p2, betU1_i1p3, betU1_i1p4, betU1_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD20 = SIMD_fd_function_dupD0_fdorder8(betU2, betU2_i0m1, betU2_i0m2, betU2_i0m3, betU2_i0p1,
                                                                                          betU2_i0p2, betU2_i0p3, betU2_i0p4, betU2_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD21 = SIMD_fd_function_dupD1_fdorder8(betU2, betU2_i1m1, betU2_i1m2, betU2_i1m3, betU2_i1p1,
                                                                                          betU2_i1p2, betU2_i1p3, betU2_i1p4, betU2_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputcf_ddnD0 =
            SIMD_fd_function_ddnD0_fdorder8(cf, cf_i0m1, cf_i0m2, cf_i0m3, cf_i0m4, cf_i0m5, cf_i0p1, cf_i0p2, cf_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputcf_ddnD1 =
            SIMD_fd_function_ddnD1_fdorder8(cf, cf_i1m1, cf_i1m2, cf_i1m3, cf_i1m4, cf_i1m5, cf_i1p1, cf_i1p2, cf_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputcf_dupD0 =
            SIMD_fd_function_dupD0_fdorder8(cf, cf_i0m1, cf_i0m2, cf_i0m3, cf_i0p1, cf_i0p2, cf_i0p3, cf_i0p4, cf_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputcf_dupD1 =
            SIMD_fd_function_dupD1_fdorder8(cf, cf_i1m1, cf_i1m2, cf_i1m3, cf_i1p1, cf_i1p2, cf_i1p3, cf_i1p4, cf_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD000 = SIMD_fd_function_ddnD0_fdorder8(hDD00, hDD00_i0m1, hDD00_i0m2, hDD00_i0m3, hDD00_i0m4,
                                                                                          hDD00_i0m5, hDD00_i0p1, hDD00_i0p2, hDD00_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD001 = SIMD_fd_function_ddnD1_fdorder8(hDD00, hDD00_i1m1, hDD00_i1m2, hDD00_i1m3, hDD00_i1m4,
                                                                                          hDD00_i1m5, hDD00_i1p1, hDD00_i1p2, hDD00_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD010 = SIMD_fd_function_ddnD0_fdorder8(hDD01, hDD01_i0m1, hDD01_i0m2, hDD01_i0m3, hDD01_i0m4,
                                                                                          hDD01_i0m5, hDD01_i0p1, hDD01_i0p2, hDD01_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD011 = SIMD_fd_function_ddnD1_fdorder8(hDD01, hDD01_i1m1, hDD01_i1m2, hDD01_i1m3, hDD01_i1m4,
                                                                                          hDD01_i1m5, hDD01_i1p1, hDD01_i1p2, hDD01_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD020 = SIMD_fd_function_ddnD0_fdorder8(hDD02, hDD02_i0m1, hDD02_i0m2, hDD02_i0m3, hDD02_i0m4,
                                                                                          hDD02_i0m5, hDD02_i0p1, hDD02_i0p2, hDD02_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD021 = SIMD_fd_function_ddnD1_fdorder8(hDD02, hDD02_i1m1, hDD02_i1m2, hDD02_i1m3, hDD02_i1m4,
                                                                                          hDD02_i1m5, hDD02_i1p1, hDD02_i1p2, hDD02_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD110 = SIMD_fd_function_ddnD0_fdorder8(hDD11, hDD11_i0m1, hDD11_i0m2, hDD11_i0m3, hDD11_i0m4,
                                                                                          hDD11_i0m5, hDD11_i0p1, hDD11_i0p2, hDD11_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD111 = SIMD_fd_function_ddnD1_fdorder8(hDD11, hDD11_i1m1, hDD11_i1m2, hDD11_i1m3, hDD11_i1m4,
                                                                                          hDD11_i1m5, hDD11_i1p1, hDD11_i1p2, hDD11_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD120 = SIMD_fd_function_ddnD0_fdorder8(hDD12, hDD12_i0m1, hDD12_i0m2, hDD12_i0m3, hDD12_i0m4,
                                                                                          hDD12_i0m5, hDD12_i0p1, hDD12_i0p2, hDD12_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD121 = SIMD_fd_function_ddnD1_fdorder8(hDD12, hDD12_i1m1, hDD12_i1m2, hDD12_i1m3, hDD12_i1m4,
                                                                                          hDD12_i1m5, hDD12_i1p1, hDD12_i1p2, hDD12_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD220 = SIMD_fd_function_ddnD0_fdorder8(hDD22, hDD22_i0m1, hDD22_i0m2, hDD22_i0m3, hDD22_i0m4,
                                                                                          hDD22_i0m5, hDD22_i0p1, hDD22_i0p2, hDD22_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD221 = SIMD_fd_function_ddnD1_fdorder8(hDD22, hDD22_i1m1, hDD22_i1m2, hDD22_i1m3, hDD22_i1m4,
                                                                                          hDD22_i1m5, hDD22_i1p1, hDD22_i1p2, hDD22_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD000 = SIMD_fd_function_dupD0_fdorder8(hDD00, hDD00_i0m1, hDD00_i0m2, hDD00_i0m3, hDD00_i0p1,
                                                                                          hDD00_i0p2, hDD00_i0p3, hDD00_i0p4, hDD00_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD001 = SIMD_fd_function_dupD1_fdorder8(hDD00, hDD00_i1m1, hDD00_i1m2, hDD00_i1m3, hDD00_i1p1,
                                                                                          hDD00_i1p2, hDD00_i1p3, hDD00_i1p4, hDD00_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD010 = SIMD_fd_function_dupD0_fdorder8(hDD01, hDD01_i0m1, hDD01_i0m2, hDD01_i0m3, hDD01_i0p1,
                                                                                          hDD01_i0p2, hDD01_i0p3, hDD01_i0p4, hDD01_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD011 = SIMD_fd_function_dupD1_fdorder8(hDD01, hDD01_i1m1, hDD01_i1m2, hDD01_i1m3, hDD01_i1p1,
                                                                                          hDD01_i1p2, hDD01_i1p3, hDD01_i1p4, hDD01_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD020 = SIMD_fd_function_dupD0_fdorder8(hDD02, hDD02_i0m1, hDD02_i0m2, hDD02_i0m3, hDD02_i0p1,
                                                                                          hDD02_i0p2, hDD02_i0p3, hDD02_i0p4, hDD02_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD021 = SIMD_fd_function_dupD1_fdorder8(hDD02, hDD02_i1m1, hDD02_i1m2, hDD02_i1m3, hDD02_i1p1,
                                                                                          hDD02_i1p2, hDD02_i1p3, hDD02_i1p4, hDD02_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD110 = SIMD_fd_function_dupD0_fdorder8(hDD11, hDD11_i0m1, hDD11_i0m2, hDD11_i0m3, hDD11_i0p1,
                                                                                          hDD11_i0p2, hDD11_i0p3, hDD11_i0p4, hDD11_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD111 = SIMD_fd_function_dupD1_fdorder8(hDD11, hDD11_i1m1, hDD11_i1m2, hDD11_i1m3, hDD11_i1p1,
                                                                                          hDD11_i1p2, hDD11_i1p3, hDD11_i1p4, hDD11_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD120 = SIMD_fd_function_dupD0_fdorder8(hDD12, hDD12_i0m1, hDD12_i0m2, hDD12_i0m3, hDD12_i0p1,
                                                                                          hDD12_i0p2, hDD12_i0p3, hDD12_i0p4, hDD12_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD121 = SIMD_fd_function_dupD1_fdorder8(hDD12, hDD12_i1m1, hDD12_i1m2, hDD12_i1m3, hDD12_i1p1,
                                                                                          hDD12_i1p2, hDD12_i1p3, hDD12_i1p4, hDD12_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD220 = SIMD_fd_function_dupD0_fdorder8(hDD22, hDD22_i0m1, hDD22_i0m2, hDD22_i0m3, hDD22_i0p1,
                                                                                          hDD22_i0p2, hDD22_i0p3, hDD22_i0p4, hDD22_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD221 = SIMD_fd_function_dupD1_fdorder8(hDD22, hDD22_i1m1, hDD22_i1m2, hDD22_i1m3, hDD22_i1p1,
                                                                                          hDD22_i1p2, hDD22_i1p3, hDD22_i1p4, hDD22_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD00 =
            SIMD_fd_function_ddnD0_fdorder8(lambdaU0, lambdaU0_i0m1, lambdaU0_i0m2, lambdaU0_i0m3, lambdaU0_i0m4, lambdaU0_i0m5, lambdaU0_i0p1,
                                            lambdaU0_i0p2, lambdaU0_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD01 =
            SIMD_fd_function_ddnD1_fdorder8(lambdaU0, lambdaU0_i1m1, lambdaU0_i1m2, lambdaU0_i1m3, lambdaU0_i1m4, lambdaU0_i1m5, lambdaU0_i1p1,
                                            lambdaU0_i1p2, lambdaU0_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD10 =
            SIMD_fd_function_ddnD0_fdorder8(lambdaU1, lambdaU1_i0m1, lambdaU1_i0m2, lambdaU1_i0m3, lambdaU1_i0m4, lambdaU1_i0m5, lambdaU1_i0p1,
                                            lambdaU1_i0p2, lambdaU1_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD11 =
            SIMD_fd_function_ddnD1_fdorder8(lambdaU1, lambdaU1_i1m1, lambdaU1_i1m2, lambdaU1_i1m3, lambdaU1_i1m4, lambdaU1_i1m5, lambdaU1_i1p1,
                                            lambdaU1_i1p2, lambdaU1_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD20 =
            SIMD_fd_function_ddnD0_fdorder8(lambdaU2, lambdaU2_i0m1, lambdaU2_i0m2, lambdaU2_i0m3, lambdaU2_i0m4, lambdaU2_i0m5, lambdaU2_i0p1,
                                            lambdaU2_i0p2, lambdaU2_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD21 =
            SIMD_fd_function_ddnD1_fdorder8(lambdaU2, lambdaU2_i1m1, lambdaU2_i1m2, lambdaU2_i1m3, lambdaU2_i1m4, lambdaU2_i1m5, lambdaU2_i1p1,
                                            lambdaU2_i1p2, lambdaU2_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD00 =
            SIMD_fd_function_dupD0_fdorder8(lambdaU0, lambdaU0_i0m1, lambdaU0_i0m2, lambdaU0_i0m3, lambdaU0_i0p1, lambdaU0_i0p2, lambdaU0_i0p3,
                                            lambdaU0_i0p4, lambdaU0_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD01 =
            SIMD_fd_function_dupD1_fdorder8(lambdaU0, lambdaU0_i1m1, lambdaU0_i1m2, lambdaU0_i1m3, lambdaU0_i1p1, lambdaU0_i1p2, lambdaU0_i1p3,
                                            lambdaU0_i1p4, lambdaU0_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD10 =
            SIMD_fd_function_dupD0_fdorder8(lambdaU1, lambdaU1_i0m1, lambdaU1_i0m2, lambdaU1_i0m3, lambdaU1_i0p1, lambdaU1_i0p2, lambdaU1_i0p3,
                                            lambdaU1_i0p4, lambdaU1_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD11 =
            SIMD_fd_function_dupD1_fdorder8(lambdaU1, lambdaU1_i1m1, lambdaU1_i1m2, lambdaU1_i1m3, lambdaU1_i1p1, lambdaU1_i1p2, lambdaU1_i1p3,
                                            lambdaU1_i1p4, lambdaU1_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD20 =
            SIMD_fd_function_dupD0_fdorder8(lambdaU2, lambdaU2_i0m1, lambdaU2_i0m2, lambdaU2_i0m3, lambdaU2_i0p1, lambdaU2_i0p2, lambdaU2_i0p3,
                                            lambdaU2_i0p4, lambdaU2_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD21 =
            SIMD_fd_function_dupD1_fdorder8(lambdaU2, lambdaU2_i1m1, lambdaU2_i1m2, lambdaU2_i1m3, lambdaU2_i1p1, lambdaU2_i1p2, lambdaU2_i1p3,
                                            lambdaU2_i1p4, lambdaU2_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputtrK_ddnD0 =
            SIMD_fd_function_ddnD0_fdorder8(trK, trK_i0m1, trK_i0m2, trK_i0m3, trK_i0m4, trK_i0m5, trK_i0p1, trK_i0p2, trK_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputtrK_ddnD1 =
            SIMD_fd_function_ddnD1_fdorder8(trK, trK_i1m1, trK_i1m2, trK_i1m3, trK_i1m4, trK_i1m5, trK_i1p1, trK_i1p2, trK_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputtrK_dupD0 =
            SIMD_fd_function_dupD0_fdorder8(trK, trK_i0m1, trK_i0m2, trK_i0m3, trK_i0p1, trK_i0p2, trK_i0p3, trK_i0p4, trK_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputtrK_dupD1 =
            SIMD_fd_function_dupD1_fdorder8(trK, trK_i1m1, trK_i1m2, trK_i1m3, trK_i1p1, trK_i1p2, trK_i1p3, trK_i1p4, trK_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD00 = SIMD_fd_function_ddnD0_fdorder8(vetU0, vetU0_i0m1, vetU0_i0m2, vetU0_i0m3, vetU0_i0m4,
                                                                                          vetU0_i0m5, vetU0_i0p1, vetU0_i0p2, vetU0_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD01 = SIMD_fd_function_ddnD1_fdorder8(vetU0, vetU0_i1m1, vetU0_i1m2, vetU0_i1m3, vetU0_i1m4,
                                                                                          vetU0_i1m5, vetU0_i1p1, vetU0_i1p2, vetU0_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD10 = SIMD_fd_function_ddnD0_fdorder8(vetU1, vetU1_i0m1, vetU1_i0m2, vetU1_i0m3, vetU1_i0m4,
                                                                                          vetU1_i0m5, vetU1_i0p1, vetU1_i0p2, vetU1_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD11 = SIMD_fd_function_ddnD1_fdorder8(vetU1, vetU1_i1m1, vetU1_i1m2, vetU1_i1m3, vetU1_i1m4,
                                                                                          vetU1_i1m5, vetU1_i1p1, vetU1_i1p2, vetU1_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD20 = SIMD_fd_function_ddnD0_fdorder8(vetU2, vetU2_i0m1, vetU2_i0m2, vetU2_i0m3, vetU2_i0m4,
                                                                                          vetU2_i0m5, vetU2_i0p1, vetU2_i0p2, vetU2_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD21 = SIMD_fd_function_ddnD1_fdorder8(vetU2, vetU2_i1m1, vetU2_i1m2, vetU2_i1m3, vetU2_i1m4,
                                                                                          vetU2_i1m5, vetU2_i1p1, vetU2_i1p2, vetU2_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD00 = SIMD_fd_function_dupD0_fdorder8(vetU0, vetU0_i0m1, vetU0_i0m2, vetU0_i0m3, vetU0_i0p1,
                                                                                          vetU0_i0p2, vetU0_i0p3, vetU0_i0p4, vetU0_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD01 = SIMD_fd_function_dupD1_fdorder8(vetU0, vetU0_i1m1, vetU0_i1m2, vetU0_i1m3, vetU0_i1p1,
                                                                                          vetU0_i1p2, vetU0_i1p3, vetU0_i1p4, vetU0_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD10 = SIMD_fd_function_dupD0_fdorder8(vetU1, vetU1_i0m1, vetU1_i0m2, vetU1_i0m3, vetU1_i0p1,
                                                                                          vetU1_i0p2, vetU1_i0p3, vetU1_i0p4, vetU1_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD11 = SIMD_fd_function_dupD1_fdorder8(vetU1, vetU1_i1m1, vetU1_i1m2, vetU1_i1m3, vetU1_i1p1,
                                                                                          vetU1_i1p2, vetU1_i1p3, vetU1_i1p4, vetU1_i1p5, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD20 = SIMD_fd_function_dupD0_fdorder8(vetU2, vetU2_i0m1, vetU2_i0m2, vetU2_i0m3, vetU2_i0p1,
                                                                                          vetU2_i0p2, vetU2_i0p3, vetU2_i0p4, vetU2_i0p5, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD21 = SIMD_fd_function_dupD1_fdorder8(vetU2, vetU2_i1m1, vetU2_i1m2, vetU2_i1m3, vetU2_i1p1,
                                                                                          vetU2_i1p2, vetU2_i1p3, vetU2_i1p4, vetU2_i1p5, invdxx1);
        const REAL_SIMD_ARRAY aDD_dKOD000 = SIMD_fd_function_dKOD0_fdorder8(aDD00, aDD00_i0m1, aDD00_i0m2, aDD00_i0m3, aDD00_i0m4, aDD00_i0m5,
                                                                            aDD00_i0p1, aDD00_i0p2, aDD00_i0p3, aDD00_i0p4, aDD00_i0p5, invdxx0);
        const REAL_SIMD_ARRAY aDD_dKOD001 = SIMD_fd_function_dKOD1_fdorder8(aDD00, aDD00_i1m1, aDD00_i1m2, aDD00_i1m3, aDD00_i1m4, aDD00_i1m5,
                                                                            aDD00_i1p1, aDD00_i1p2, aDD00_i1p3, aDD00_i1p4, aDD00_i1p5, invdxx1);
        const REAL_SIMD_ARRAY aDD_dKOD010 = SIMD_fd_function_dKOD0_fdorder8(aDD01, aDD01_i0m1, aDD01_i0m2, aDD01_i0m3, aDD01_i0m4, aDD01_i0m5,
                                                                            aDD01_i0p1, aDD01_i0p2, aDD01_i0p3, aDD01_i0p4, aDD01_i0p5, invdxx0);
        const REAL_SIMD_ARRAY aDD_dKOD011 = SIMD_fd_function_dKOD1_fdorder8(aDD01, aDD01_i1m1, aDD01_i1m2, aDD01_i1m3, aDD01_i1m4, aDD01_i1m5,
                                                                            aDD01_i1p1, aDD01_i1p2, aDD01_i1p3, aDD01_i1p4, aDD01_i1p5, invdxx1);
        const REAL_SIMD_ARRAY aDD_dKOD020 = SIMD_fd_function_dKOD0_fdorder8(aDD02, aDD02_i0m1, aDD02_i0m2, aDD02_i0m3, aDD02_i0m4, aDD02_i0m5,
                                                                            aDD02_i0p1, aDD02_i0p2, aDD02_i0p3, aDD02_i0p4, aDD02_i0p5, invdxx0);
        const REAL_SIMD_ARRAY aDD_dKOD021 = SIMD_fd_function_dKOD1_fdorder8(aDD02, aDD02_i1m1, aDD02_i1m2, aDD02_i1m3, aDD02_i1m4, aDD02_i1m5,
                                                                            aDD02_i1p1, aDD02_i1p2, aDD02_i1p3, aDD02_i1p4, aDD02_i1p5, invdxx1);
        const REAL_SIMD_ARRAY aDD_dKOD110 = SIMD_fd_function_dKOD0_fdorder8(aDD11, aDD11_i0m1, aDD11_i0m2, aDD11_i0m3, aDD11_i0m4, aDD11_i0m5,
                                                                            aDD11_i0p1, aDD11_i0p2, aDD11_i0p3, aDD11_i0p4, aDD11_i0p5, invdxx0);
        const REAL_SIMD_ARRAY aDD_dKOD111 = SIMD_fd_function_dKOD1_fdorder8(aDD11, aDD11_i1m1, aDD11_i1m2, aDD11_i1m3, aDD11_i1m4, aDD11_i1m5,
                                                                            aDD11_i1p1, aDD11_i1p2, aDD11_i1p3, aDD11_i1p4, aDD11_i1p5, invdxx1);
        const REAL_SIMD_ARRAY aDD_dKOD120 = SIMD_fd_function_dKOD0_fdorder8(aDD12, aDD12_i0m1, aDD12_i0m2, aDD12_i0m3, aDD12_i0m4, aDD12_i0m5,
                                                                            aDD12_i0p1, aDD12_i0p2, aDD12_i0p3, aDD12_i0p4, aDD12_i0p5, invdxx0);
        const REAL_SIMD_ARRAY aDD_dKOD121 = SIMD_fd_function_dKOD1_fdorder8(aDD12, aDD12_i1m1, aDD12_i1m2, aDD12_i1m3, aDD12_i1m4, aDD12_i1m5,
                                                                            aDD12_i1p1, aDD12_i1p2, aDD12_i1p3, aDD12_i1p4, aDD12_i1p5, invdxx1);
        const REAL_SIMD_ARRAY aDD_dKOD220 = SIMD_fd_function_dKOD0_fdorder8(aDD22, aDD22_i0m1, aDD22_i0m2, aDD22_i0m3, aDD22_i0m4, aDD22_i0m5,
                                                                            aDD22_i0p1, aDD22_i0p2, aDD22_i0p3, aDD22_i0p4, aDD22_i0p5, invdxx0);
        const REAL_SIMD_ARRAY aDD_dKOD221 = SIMD_fd_function_dKOD1_fdorder8(aDD22, aDD22_i1m1, aDD22_i1m2, aDD22_i1m3, aDD22_i1m4, aDD22_i1m5,
                                                                            aDD22_i1p1, aDD22_i1p2, aDD22_i1p3, aDD22_i1p4, aDD22_i1p5, invdxx1);
        const REAL_SIMD_ARRAY alpha_dD0 =
            SIMD_fd_function_dD0_fdorder8(alpha_i0m1, alpha_i0m2, alpha_i0m3, alpha_i0m4, alpha_i0p1, alpha_i0p2, alpha_i0p3, alpha_i0p4, invdxx0);
        const REAL_SIMD_ARRAY alpha_dD1 =
            SIMD_fd_function_dD1_fdorder8(alpha_i1m1, alpha_i1m2, alpha_i1m3, alpha_i1m4, alpha_i1p1, alpha_i1p2, alpha_i1p3, alpha_i1p4, invdxx1);
        const REAL_SIMD_ARRAY alpha_dDD00 = SIMD_fd_function_dDD00_fdorder8(alpha, alpha_i0m1, alpha_i0m2, alpha_i0m3, alpha_i0m4, alpha_i0p1,
                                                                            alpha_i0p2, alpha_i0p3, alpha_i0p4, invdxx0);
        const REAL_SIMD_ARRAY alpha_dDD01 = SIMD_fd_function_dDD01_fdorder8(
            alpha_i0m1_i1m1, alpha_i0m1_i1m2, alpha_i0m1_i1m3, alpha_i0m1_i1m4, alpha_i0m1_i1p1, alpha_i0m1_i1p2, alpha_i0m1_i1p3, alpha_i0m1_i1p4,
            alpha_i0m2_i1m1, alpha_i0m2_i1m2, alpha_i0m2_i1m3, alpha_i0m2_i1m4, alpha_i0m2_i1p1, alpha_i0m2_i1p2, alpha_i0m2_i1p3, alpha_i0m2_i1p4,
            alpha_i0m3_i1m1, alpha_i0m3_i1m2, alpha_i0m3_i1m3, alpha_i0m3_i1m4, alpha_i0m3_i1p1, alpha_i0m3_i1p2, alpha_i0m3_i1p3, alpha_i0m3_i1p4,
            alpha_i0m4_i1m1, alpha_i0m4_i1m2, alpha_i0m4_i1m3, alpha_i0m4_i1m4, alpha_i0m4_i1p1, alpha_i0m4_i1p2, alpha_i0m4_i1p3, alpha_i0m4_i1p4,
            alpha_i0p1_i1m1, alpha_i0p1_i1m2, alpha_i0p1_i1m3, alpha_i0p1_i1m4, alpha_i0p1_i1p1, alpha_i0p1_i1p2, alpha_i0p1_i1p3, alpha_i0p1_i1p4,
            alpha_i0p2_i1m1, alpha_i0p2_i1m2, alpha_i0p2_i1m3, alpha_i0p2_i1m4, alpha_i0p2_i1p1, alpha_i0p2_i1p2, alpha_i0p2_i1p3, alpha_i0p2_i1p4,
            alpha_i0p3_i1m1, alpha_i0p3_i1m2, alpha_i0p3_i1m3, alpha_i0p3_i1m4, alpha_i0p3_i1p1, alpha_i0p3_i1p2, alpha_i0p3_i1p3, alpha_i0p3_i1p4,
            alpha_i0p4_i1m1, alpha_i0p4_i1m2, alpha_i0p4_i1m3, alpha_i0p4_i1m4, alpha_i0p4_i1p1, alpha_i0p4_i1p2, alpha_i0p4_i1p3, alpha_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL_SIMD_ARRAY alpha_dDD11 = SIMD_fd_function_dDD11_fdorder8(alpha, alpha_i1m1, alpha_i1m2, alpha_i1m3, alpha_i1m4, alpha_i1p1,
                                                                            alpha_i1p2, alpha_i1p3, alpha_i1p4, invdxx1);
        const REAL_SIMD_ARRAY alpha_dKOD0 = SIMD_fd_function_dKOD0_fdorder8(alpha, alpha_i0m1, alpha_i0m2, alpha_i0m3, alpha_i0m4, alpha_i0m5,
                                                                            alpha_i0p1, alpha_i0p2, alpha_i0p3, alpha_i0p4, alpha_i0p5, invdxx0);
        const REAL_SIMD_ARRAY alpha_dKOD1 = SIMD_fd_function_dKOD1_fdorder8(alpha, alpha_i1m1, alpha_i1m2, alpha_i1m3, alpha_i1m4, alpha_i1m5,
                                                                            alpha_i1p1, alpha_i1p2, alpha_i1p3, alpha_i1p4, alpha_i1p5, invdxx1);
        const REAL_SIMD_ARRAY betU_dKOD00 = SIMD_fd_function_dKOD0_fdorder8(betU0, betU0_i0m1, betU0_i0m2, betU0_i0m3, betU0_i0m4, betU0_i0m5,
                                                                            betU0_i0p1, betU0_i0p2, betU0_i0p3, betU0_i0p4, betU0_i0p5, invdxx0);
        const REAL_SIMD_ARRAY betU_dKOD01 = SIMD_fd_function_dKOD1_fdorder8(betU0, betU0_i1m1, betU0_i1m2, betU0_i1m3, betU0_i1m4, betU0_i1m5,
                                                                            betU0_i1p1, betU0_i1p2, betU0_i1p3, betU0_i1p4, betU0_i1p5, invdxx1);
        const REAL_SIMD_ARRAY betU_dKOD10 = SIMD_fd_function_dKOD0_fdorder8(betU1, betU1_i0m1, betU1_i0m2, betU1_i0m3, betU1_i0m4, betU1_i0m5,
                                                                            betU1_i0p1, betU1_i0p2, betU1_i0p3, betU1_i0p4, betU1_i0p5, invdxx0);
        const REAL_SIMD_ARRAY betU_dKOD11 = SIMD_fd_function_dKOD1_fdorder8(betU1, betU1_i1m1, betU1_i1m2, betU1_i1m3, betU1_i1m4, betU1_i1m5,
                                                                            betU1_i1p1, betU1_i1p2, betU1_i1p3, betU1_i1p4, betU1_i1p5, invdxx1);
        const REAL_SIMD_ARRAY betU_dKOD20 = SIMD_fd_function_dKOD0_fdorder8(betU2, betU2_i0m1, betU2_i0m2, betU2_i0m3, betU2_i0m4, betU2_i0m5,
                                                                            betU2_i0p1, betU2_i0p2, betU2_i0p3, betU2_i0p4, betU2_i0p5, invdxx0);
        const REAL_SIMD_ARRAY betU_dKOD21 = SIMD_fd_function_dKOD1_fdorder8(betU2, betU2_i1m1, betU2_i1m2, betU2_i1m3, betU2_i1m4, betU2_i1m5,
                                                                            betU2_i1p1, betU2_i1p2, betU2_i1p3, betU2_i1p4, betU2_i1p5, invdxx1);
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
        const REAL_SIMD_ARRAY cf_dKOD0 =
            SIMD_fd_function_dKOD0_fdorder8(cf, cf_i0m1, cf_i0m2, cf_i0m3, cf_i0m4, cf_i0m5, cf_i0p1, cf_i0p2, cf_i0p3, cf_i0p4, cf_i0p5, invdxx0);
        const REAL_SIMD_ARRAY cf_dKOD1 =
            SIMD_fd_function_dKOD1_fdorder8(cf, cf_i1m1, cf_i1m2, cf_i1m3, cf_i1m4, cf_i1m5, cf_i1p1, cf_i1p2, cf_i1p3, cf_i1p4, cf_i1p5, invdxx1);
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
        const REAL_SIMD_ARRAY hDD_dKOD000 = SIMD_fd_function_dKOD0_fdorder8(hDD00, hDD00_i0m1, hDD00_i0m2, hDD00_i0m3, hDD00_i0m4, hDD00_i0m5,
                                                                            hDD00_i0p1, hDD00_i0p2, hDD00_i0p3, hDD00_i0p4, hDD00_i0p5, invdxx0);
        const REAL_SIMD_ARRAY hDD_dKOD001 = SIMD_fd_function_dKOD1_fdorder8(hDD00, hDD00_i1m1, hDD00_i1m2, hDD00_i1m3, hDD00_i1m4, hDD00_i1m5,
                                                                            hDD00_i1p1, hDD00_i1p2, hDD00_i1p3, hDD00_i1p4, hDD00_i1p5, invdxx1);
        const REAL_SIMD_ARRAY hDD_dKOD010 = SIMD_fd_function_dKOD0_fdorder8(hDD01, hDD01_i0m1, hDD01_i0m2, hDD01_i0m3, hDD01_i0m4, hDD01_i0m5,
                                                                            hDD01_i0p1, hDD01_i0p2, hDD01_i0p3, hDD01_i0p4, hDD01_i0p5, invdxx0);
        const REAL_SIMD_ARRAY hDD_dKOD011 = SIMD_fd_function_dKOD1_fdorder8(hDD01, hDD01_i1m1, hDD01_i1m2, hDD01_i1m3, hDD01_i1m4, hDD01_i1m5,
                                                                            hDD01_i1p1, hDD01_i1p2, hDD01_i1p3, hDD01_i1p4, hDD01_i1p5, invdxx1);
        const REAL_SIMD_ARRAY hDD_dKOD020 = SIMD_fd_function_dKOD0_fdorder8(hDD02, hDD02_i0m1, hDD02_i0m2, hDD02_i0m3, hDD02_i0m4, hDD02_i0m5,
                                                                            hDD02_i0p1, hDD02_i0p2, hDD02_i0p3, hDD02_i0p4, hDD02_i0p5, invdxx0);
        const REAL_SIMD_ARRAY hDD_dKOD021 = SIMD_fd_function_dKOD1_fdorder8(hDD02, hDD02_i1m1, hDD02_i1m2, hDD02_i1m3, hDD02_i1m4, hDD02_i1m5,
                                                                            hDD02_i1p1, hDD02_i1p2, hDD02_i1p3, hDD02_i1p4, hDD02_i1p5, invdxx1);
        const REAL_SIMD_ARRAY hDD_dKOD110 = SIMD_fd_function_dKOD0_fdorder8(hDD11, hDD11_i0m1, hDD11_i0m2, hDD11_i0m3, hDD11_i0m4, hDD11_i0m5,
                                                                            hDD11_i0p1, hDD11_i0p2, hDD11_i0p3, hDD11_i0p4, hDD11_i0p5, invdxx0);
        const REAL_SIMD_ARRAY hDD_dKOD111 = SIMD_fd_function_dKOD1_fdorder8(hDD11, hDD11_i1m1, hDD11_i1m2, hDD11_i1m3, hDD11_i1m4, hDD11_i1m5,
                                                                            hDD11_i1p1, hDD11_i1p2, hDD11_i1p3, hDD11_i1p4, hDD11_i1p5, invdxx1);
        const REAL_SIMD_ARRAY hDD_dKOD120 = SIMD_fd_function_dKOD0_fdorder8(hDD12, hDD12_i0m1, hDD12_i0m2, hDD12_i0m3, hDD12_i0m4, hDD12_i0m5,
                                                                            hDD12_i0p1, hDD12_i0p2, hDD12_i0p3, hDD12_i0p4, hDD12_i0p5, invdxx0);
        const REAL_SIMD_ARRAY hDD_dKOD121 = SIMD_fd_function_dKOD1_fdorder8(hDD12, hDD12_i1m1, hDD12_i1m2, hDD12_i1m3, hDD12_i1m4, hDD12_i1m5,
                                                                            hDD12_i1p1, hDD12_i1p2, hDD12_i1p3, hDD12_i1p4, hDD12_i1p5, invdxx1);
        const REAL_SIMD_ARRAY hDD_dKOD220 = SIMD_fd_function_dKOD0_fdorder8(hDD22, hDD22_i0m1, hDD22_i0m2, hDD22_i0m3, hDD22_i0m4, hDD22_i0m5,
                                                                            hDD22_i0p1, hDD22_i0p2, hDD22_i0p3, hDD22_i0p4, hDD22_i0p5, invdxx0);
        const REAL_SIMD_ARRAY hDD_dKOD221 = SIMD_fd_function_dKOD1_fdorder8(hDD22, hDD22_i1m1, hDD22_i1m2, hDD22_i1m3, hDD22_i1m4, hDD22_i1m5,
                                                                            hDD22_i1p1, hDD22_i1p2, hDD22_i1p3, hDD22_i1p4, hDD22_i1p5, invdxx1);
        const REAL_SIMD_ARRAY lambdaU_dKOD00 =
            SIMD_fd_function_dKOD0_fdorder8(lambdaU0, lambdaU0_i0m1, lambdaU0_i0m2, lambdaU0_i0m3, lambdaU0_i0m4, lambdaU0_i0m5, lambdaU0_i0p1,
                                            lambdaU0_i0p2, lambdaU0_i0p3, lambdaU0_i0p4, lambdaU0_i0p5, invdxx0);
        const REAL_SIMD_ARRAY lambdaU_dKOD01 =
            SIMD_fd_function_dKOD1_fdorder8(lambdaU0, lambdaU0_i1m1, lambdaU0_i1m2, lambdaU0_i1m3, lambdaU0_i1m4, lambdaU0_i1m5, lambdaU0_i1p1,
                                            lambdaU0_i1p2, lambdaU0_i1p3, lambdaU0_i1p4, lambdaU0_i1p5, invdxx1);
        const REAL_SIMD_ARRAY lambdaU_dKOD10 =
            SIMD_fd_function_dKOD0_fdorder8(lambdaU1, lambdaU1_i0m1, lambdaU1_i0m2, lambdaU1_i0m3, lambdaU1_i0m4, lambdaU1_i0m5, lambdaU1_i0p1,
                                            lambdaU1_i0p2, lambdaU1_i0p3, lambdaU1_i0p4, lambdaU1_i0p5, invdxx0);
        const REAL_SIMD_ARRAY lambdaU_dKOD11 =
            SIMD_fd_function_dKOD1_fdorder8(lambdaU1, lambdaU1_i1m1, lambdaU1_i1m2, lambdaU1_i1m3, lambdaU1_i1m4, lambdaU1_i1m5, lambdaU1_i1p1,
                                            lambdaU1_i1p2, lambdaU1_i1p3, lambdaU1_i1p4, lambdaU1_i1p5, invdxx1);
        const REAL_SIMD_ARRAY lambdaU_dKOD20 =
            SIMD_fd_function_dKOD0_fdorder8(lambdaU2, lambdaU2_i0m1, lambdaU2_i0m2, lambdaU2_i0m3, lambdaU2_i0m4, lambdaU2_i0m5, lambdaU2_i0p1,
                                            lambdaU2_i0p2, lambdaU2_i0p3, lambdaU2_i0p4, lambdaU2_i0p5, invdxx0);
        const REAL_SIMD_ARRAY lambdaU_dKOD21 =
            SIMD_fd_function_dKOD1_fdorder8(lambdaU2, lambdaU2_i1m1, lambdaU2_i1m2, lambdaU2_i1m3, lambdaU2_i1m4, lambdaU2_i1m5, lambdaU2_i1p1,
                                            lambdaU2_i1p2, lambdaU2_i1p3, lambdaU2_i1p4, lambdaU2_i1p5, invdxx1);
        const REAL_SIMD_ARRAY trK_dD0 =
            SIMD_fd_function_dD0_fdorder8(trK_i0m1, trK_i0m2, trK_i0m3, trK_i0m4, trK_i0p1, trK_i0p2, trK_i0p3, trK_i0p4, invdxx0);
        const REAL_SIMD_ARRAY trK_dD1 =
            SIMD_fd_function_dD1_fdorder8(trK_i1m1, trK_i1m2, trK_i1m3, trK_i1m4, trK_i1p1, trK_i1p2, trK_i1p3, trK_i1p4, invdxx1);
        const REAL_SIMD_ARRAY trK_dKOD0 = SIMD_fd_function_dKOD0_fdorder8(trK, trK_i0m1, trK_i0m2, trK_i0m3, trK_i0m4, trK_i0m5, trK_i0p1, trK_i0p2,
                                                                          trK_i0p3, trK_i0p4, trK_i0p5, invdxx0);
        const REAL_SIMD_ARRAY trK_dKOD1 = SIMD_fd_function_dKOD1_fdorder8(trK, trK_i1m1, trK_i1m2, trK_i1m3, trK_i1m4, trK_i1m5, trK_i1p1, trK_i1p2,
                                                                          trK_i1p3, trK_i1p4, trK_i1p5, invdxx1);
        const REAL_SIMD_ARRAY vetU_dD00 =
            SIMD_fd_function_dD0_fdorder8(vetU0_i0m1, vetU0_i0m2, vetU0_i0m3, vetU0_i0m4, vetU0_i0p1, vetU0_i0p2, vetU0_i0p3, vetU0_i0p4, invdxx0);
        const REAL_SIMD_ARRAY vetU_dD01 =
            SIMD_fd_function_dD1_fdorder8(vetU0_i1m1, vetU0_i1m2, vetU0_i1m3, vetU0_i1m4, vetU0_i1p1, vetU0_i1p2, vetU0_i1p3, vetU0_i1p4, invdxx1);
        const REAL_SIMD_ARRAY vetU_dD10 =
            SIMD_fd_function_dD0_fdorder8(vetU1_i0m1, vetU1_i0m2, vetU1_i0m3, vetU1_i0m4, vetU1_i0p1, vetU1_i0p2, vetU1_i0p3, vetU1_i0p4, invdxx0);
        const REAL_SIMD_ARRAY vetU_dD11 =
            SIMD_fd_function_dD1_fdorder8(vetU1_i1m1, vetU1_i1m2, vetU1_i1m3, vetU1_i1m4, vetU1_i1p1, vetU1_i1p2, vetU1_i1p3, vetU1_i1p4, invdxx1);
        const REAL_SIMD_ARRAY vetU_dD20 =
            SIMD_fd_function_dD0_fdorder8(vetU2_i0m1, vetU2_i0m2, vetU2_i0m3, vetU2_i0m4, vetU2_i0p1, vetU2_i0p2, vetU2_i0p3, vetU2_i0p4, invdxx0);
        const REAL_SIMD_ARRAY vetU_dD21 =
            SIMD_fd_function_dD1_fdorder8(vetU2_i1m1, vetU2_i1m2, vetU2_i1m3, vetU2_i1m4, vetU2_i1p1, vetU2_i1p2, vetU2_i1p3, vetU2_i1p4, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD000 = SIMD_fd_function_dDD00_fdorder8(vetU0, vetU0_i0m1, vetU0_i0m2, vetU0_i0m3, vetU0_i0m4, vetU0_i0p1,
                                                                            vetU0_i0p2, vetU0_i0p3, vetU0_i0p4, invdxx0);
        const REAL_SIMD_ARRAY vetU_dDD001 = SIMD_fd_function_dDD01_fdorder8(
            vetU0_i0m1_i1m1, vetU0_i0m1_i1m2, vetU0_i0m1_i1m3, vetU0_i0m1_i1m4, vetU0_i0m1_i1p1, vetU0_i0m1_i1p2, vetU0_i0m1_i1p3, vetU0_i0m1_i1p4,
            vetU0_i0m2_i1m1, vetU0_i0m2_i1m2, vetU0_i0m2_i1m3, vetU0_i0m2_i1m4, vetU0_i0m2_i1p1, vetU0_i0m2_i1p2, vetU0_i0m2_i1p3, vetU0_i0m2_i1p4,
            vetU0_i0m3_i1m1, vetU0_i0m3_i1m2, vetU0_i0m3_i1m3, vetU0_i0m3_i1m4, vetU0_i0m3_i1p1, vetU0_i0m3_i1p2, vetU0_i0m3_i1p3, vetU0_i0m3_i1p4,
            vetU0_i0m4_i1m1, vetU0_i0m4_i1m2, vetU0_i0m4_i1m3, vetU0_i0m4_i1m4, vetU0_i0m4_i1p1, vetU0_i0m4_i1p2, vetU0_i0m4_i1p3, vetU0_i0m4_i1p4,
            vetU0_i0p1_i1m1, vetU0_i0p1_i1m2, vetU0_i0p1_i1m3, vetU0_i0p1_i1m4, vetU0_i0p1_i1p1, vetU0_i0p1_i1p2, vetU0_i0p1_i1p3, vetU0_i0p1_i1p4,
            vetU0_i0p2_i1m1, vetU0_i0p2_i1m2, vetU0_i0p2_i1m3, vetU0_i0p2_i1m4, vetU0_i0p2_i1p1, vetU0_i0p2_i1p2, vetU0_i0p2_i1p3, vetU0_i0p2_i1p4,
            vetU0_i0p3_i1m1, vetU0_i0p3_i1m2, vetU0_i0p3_i1m3, vetU0_i0p3_i1m4, vetU0_i0p3_i1p1, vetU0_i0p3_i1p2, vetU0_i0p3_i1p3, vetU0_i0p3_i1p4,
            vetU0_i0p4_i1m1, vetU0_i0p4_i1m2, vetU0_i0p4_i1m3, vetU0_i0p4_i1m4, vetU0_i0p4_i1p1, vetU0_i0p4_i1p2, vetU0_i0p4_i1p3, vetU0_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD011 = SIMD_fd_function_dDD11_fdorder8(vetU0, vetU0_i1m1, vetU0_i1m2, vetU0_i1m3, vetU0_i1m4, vetU0_i1p1,
                                                                            vetU0_i1p2, vetU0_i1p3, vetU0_i1p4, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD100 = SIMD_fd_function_dDD00_fdorder8(vetU1, vetU1_i0m1, vetU1_i0m2, vetU1_i0m3, vetU1_i0m4, vetU1_i0p1,
                                                                            vetU1_i0p2, vetU1_i0p3, vetU1_i0p4, invdxx0);
        const REAL_SIMD_ARRAY vetU_dDD101 = SIMD_fd_function_dDD01_fdorder8(
            vetU1_i0m1_i1m1, vetU1_i0m1_i1m2, vetU1_i0m1_i1m3, vetU1_i0m1_i1m4, vetU1_i0m1_i1p1, vetU1_i0m1_i1p2, vetU1_i0m1_i1p3, vetU1_i0m1_i1p4,
            vetU1_i0m2_i1m1, vetU1_i0m2_i1m2, vetU1_i0m2_i1m3, vetU1_i0m2_i1m4, vetU1_i0m2_i1p1, vetU1_i0m2_i1p2, vetU1_i0m2_i1p3, vetU1_i0m2_i1p4,
            vetU1_i0m3_i1m1, vetU1_i0m3_i1m2, vetU1_i0m3_i1m3, vetU1_i0m3_i1m4, vetU1_i0m3_i1p1, vetU1_i0m3_i1p2, vetU1_i0m3_i1p3, vetU1_i0m3_i1p4,
            vetU1_i0m4_i1m1, vetU1_i0m4_i1m2, vetU1_i0m4_i1m3, vetU1_i0m4_i1m4, vetU1_i0m4_i1p1, vetU1_i0m4_i1p2, vetU1_i0m4_i1p3, vetU1_i0m4_i1p4,
            vetU1_i0p1_i1m1, vetU1_i0p1_i1m2, vetU1_i0p1_i1m3, vetU1_i0p1_i1m4, vetU1_i0p1_i1p1, vetU1_i0p1_i1p2, vetU1_i0p1_i1p3, vetU1_i0p1_i1p4,
            vetU1_i0p2_i1m1, vetU1_i0p2_i1m2, vetU1_i0p2_i1m3, vetU1_i0p2_i1m4, vetU1_i0p2_i1p1, vetU1_i0p2_i1p2, vetU1_i0p2_i1p3, vetU1_i0p2_i1p4,
            vetU1_i0p3_i1m1, vetU1_i0p3_i1m2, vetU1_i0p3_i1m3, vetU1_i0p3_i1m4, vetU1_i0p3_i1p1, vetU1_i0p3_i1p2, vetU1_i0p3_i1p3, vetU1_i0p3_i1p4,
            vetU1_i0p4_i1m1, vetU1_i0p4_i1m2, vetU1_i0p4_i1m3, vetU1_i0p4_i1m4, vetU1_i0p4_i1p1, vetU1_i0p4_i1p2, vetU1_i0p4_i1p3, vetU1_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD111 = SIMD_fd_function_dDD11_fdorder8(vetU1, vetU1_i1m1, vetU1_i1m2, vetU1_i1m3, vetU1_i1m4, vetU1_i1p1,
                                                                            vetU1_i1p2, vetU1_i1p3, vetU1_i1p4, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD200 = SIMD_fd_function_dDD00_fdorder8(vetU2, vetU2_i0m1, vetU2_i0m2, vetU2_i0m3, vetU2_i0m4, vetU2_i0p1,
                                                                            vetU2_i0p2, vetU2_i0p3, vetU2_i0p4, invdxx0);
        const REAL_SIMD_ARRAY vetU_dDD201 = SIMD_fd_function_dDD01_fdorder8(
            vetU2_i0m1_i1m1, vetU2_i0m1_i1m2, vetU2_i0m1_i1m3, vetU2_i0m1_i1m4, vetU2_i0m1_i1p1, vetU2_i0m1_i1p2, vetU2_i0m1_i1p3, vetU2_i0m1_i1p4,
            vetU2_i0m2_i1m1, vetU2_i0m2_i1m2, vetU2_i0m2_i1m3, vetU2_i0m2_i1m4, vetU2_i0m2_i1p1, vetU2_i0m2_i1p2, vetU2_i0m2_i1p3, vetU2_i0m2_i1p4,
            vetU2_i0m3_i1m1, vetU2_i0m3_i1m2, vetU2_i0m3_i1m3, vetU2_i0m3_i1m4, vetU2_i0m3_i1p1, vetU2_i0m3_i1p2, vetU2_i0m3_i1p3, vetU2_i0m3_i1p4,
            vetU2_i0m4_i1m1, vetU2_i0m4_i1m2, vetU2_i0m4_i1m3, vetU2_i0m4_i1m4, vetU2_i0m4_i1p1, vetU2_i0m4_i1p2, vetU2_i0m4_i1p3, vetU2_i0m4_i1p4,
            vetU2_i0p1_i1m1, vetU2_i0p1_i1m2, vetU2_i0p1_i1m3, vetU2_i0p1_i1m4, vetU2_i0p1_i1p1, vetU2_i0p1_i1p2, vetU2_i0p1_i1p3, vetU2_i0p1_i1p4,
            vetU2_i0p2_i1m1, vetU2_i0p2_i1m2, vetU2_i0p2_i1m3, vetU2_i0p2_i1m4, vetU2_i0p2_i1p1, vetU2_i0p2_i1p2, vetU2_i0p2_i1p3, vetU2_i0p2_i1p4,
            vetU2_i0p3_i1m1, vetU2_i0p3_i1m2, vetU2_i0p3_i1m3, vetU2_i0p3_i1m4, vetU2_i0p3_i1p1, vetU2_i0p3_i1p2, vetU2_i0p3_i1p3, vetU2_i0p3_i1p4,
            vetU2_i0p4_i1m1, vetU2_i0p4_i1m2, vetU2_i0p4_i1m3, vetU2_i0p4_i1m4, vetU2_i0p4_i1p1, vetU2_i0p4_i1p2, vetU2_i0p4_i1p3, vetU2_i0p4_i1p4,
            invdxx0, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD211 = SIMD_fd_function_dDD11_fdorder8(vetU2, vetU2_i1m1, vetU2_i1m2, vetU2_i1m3, vetU2_i1m4, vetU2_i1p1,
                                                                            vetU2_i1p2, vetU2_i1p3, vetU2_i1p4, invdxx1);
        const REAL_SIMD_ARRAY vetU_dKOD00 = SIMD_fd_function_dKOD0_fdorder8(vetU0, vetU0_i0m1, vetU0_i0m2, vetU0_i0m3, vetU0_i0m4, vetU0_i0m5,
                                                                            vetU0_i0p1, vetU0_i0p2, vetU0_i0p3, vetU0_i0p4, vetU0_i0p5, invdxx0);
        const REAL_SIMD_ARRAY vetU_dKOD01 = SIMD_fd_function_dKOD1_fdorder8(vetU0, vetU0_i1m1, vetU0_i1m2, vetU0_i1m3, vetU0_i1m4, vetU0_i1m5,
                                                                            vetU0_i1p1, vetU0_i1p2, vetU0_i1p3, vetU0_i1p4, vetU0_i1p5, invdxx1);
        const REAL_SIMD_ARRAY vetU_dKOD10 = SIMD_fd_function_dKOD0_fdorder8(vetU1, vetU1_i0m1, vetU1_i0m2, vetU1_i0m3, vetU1_i0m4, vetU1_i0m5,
                                                                            vetU1_i0p1, vetU1_i0p2, vetU1_i0p3, vetU1_i0p4, vetU1_i0p5, invdxx0);
        const REAL_SIMD_ARRAY vetU_dKOD11 = SIMD_fd_function_dKOD1_fdorder8(vetU1, vetU1_i1m1, vetU1_i1m2, vetU1_i1m3, vetU1_i1m4, vetU1_i1m5,
                                                                            vetU1_i1p1, vetU1_i1p2, vetU1_i1p3, vetU1_i1p4, vetU1_i1p5, invdxx1);
        const REAL_SIMD_ARRAY vetU_dKOD20 = SIMD_fd_function_dKOD0_fdorder8(vetU2, vetU2_i0m1, vetU2_i0m2, vetU2_i0m3, vetU2_i0m4, vetU2_i0m5,
                                                                            vetU2_i0p1, vetU2_i0p2, vetU2_i0p3, vetU2_i0p4, vetU2_i0p5, invdxx0);
        const REAL_SIMD_ARRAY vetU_dKOD21 = SIMD_fd_function_dKOD1_fdorder8(vetU2, vetU2_i1m1, vetU2_i1m2, vetU2_i1m3, vetU2_i1m4, vetU2_i1m5,
                                                                            vetU2_i1p1, vetU2_i1p2, vetU2_i1p3, vetU2_i1p4, vetU2_i1p5, invdxx1);
        const double dblFDPart1_Integer_1 = 1.0;
        const REAL_SIMD_ARRAY FDPart1_Integer_1 = ConstSIMD(dblFDPart1_Integer_1);

        const double dblFDPart1_NegativeOne_ = -1.0;
        const REAL_SIMD_ARRAY FDPart1_NegativeOne_ = ConstSIMD(dblFDPart1_NegativeOne_);

        const REAL_SIMD_ARRAY UpwindControlVectorU0 = DivSIMD(vetU0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY UpwindControlVectorU1 = DivSIMD(vetU1, f0_of_xx0);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 3:
         * Implement upwinding algorithm.
         */
        const double tmp_upwind_Integer_0 = 0.000000000000000000000000000000000;

        const REAL_SIMD_ARRAY upwind_Integer_0 = ConstSIMD(tmp_upwind_Integer_0);
        const double tmp_upwind_Integer_1 = 1.000000000000000000000000000000000;

        const REAL_SIMD_ARRAY upwind_Integer_1 = ConstSIMD(tmp_upwind_Integer_1);
        const REAL_SIMD_ARRAY Upwind0 = UPWIND_ALG(UpwindControlVectorU0);
        const REAL_SIMD_ARRAY Upwind1 = UPWIND_ALG(UpwindControlVectorU1);
        const double dblFDPart2_NegativeOne_ = -1.0;
        const REAL_SIMD_ARRAY FDPart2_NegativeOne_ = ConstSIMD(dblFDPart2_NegativeOne_);

        const REAL_SIMD_ARRAY aDD_dupD000 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD000, UpwindAlgInputaDD_ddnD000), UpwindAlgInputaDD_ddnD000);
        const REAL_SIMD_ARRAY aDD_dupD001 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD001, UpwindAlgInputaDD_ddnD001), UpwindAlgInputaDD_ddnD001);
        const REAL_SIMD_ARRAY aDD_dupD010 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD010, UpwindAlgInputaDD_ddnD010), UpwindAlgInputaDD_ddnD010);
        const REAL_SIMD_ARRAY aDD_dupD011 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD011, UpwindAlgInputaDD_ddnD011), UpwindAlgInputaDD_ddnD011);
        const REAL_SIMD_ARRAY aDD_dupD020 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD020, UpwindAlgInputaDD_ddnD020), UpwindAlgInputaDD_ddnD020);
        const REAL_SIMD_ARRAY aDD_dupD021 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD021, UpwindAlgInputaDD_ddnD021), UpwindAlgInputaDD_ddnD021);
        const REAL_SIMD_ARRAY aDD_dupD110 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD110, UpwindAlgInputaDD_ddnD110), UpwindAlgInputaDD_ddnD110);
        const REAL_SIMD_ARRAY aDD_dupD111 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD111, UpwindAlgInputaDD_ddnD111), UpwindAlgInputaDD_ddnD111);
        const REAL_SIMD_ARRAY aDD_dupD120 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD120, UpwindAlgInputaDD_ddnD120), UpwindAlgInputaDD_ddnD120);
        const REAL_SIMD_ARRAY aDD_dupD121 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD121, UpwindAlgInputaDD_ddnD121), UpwindAlgInputaDD_ddnD121);
        const REAL_SIMD_ARRAY aDD_dupD220 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD220, UpwindAlgInputaDD_ddnD220), UpwindAlgInputaDD_ddnD220);
        const REAL_SIMD_ARRAY aDD_dupD221 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD221, UpwindAlgInputaDD_ddnD221), UpwindAlgInputaDD_ddnD221);
        const REAL_SIMD_ARRAY alpha_dupD0 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputalpha_dupD0, UpwindAlgInputalpha_ddnD0), UpwindAlgInputalpha_ddnD0);
        const REAL_SIMD_ARRAY alpha_dupD1 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputalpha_dupD1, UpwindAlgInputalpha_ddnD1), UpwindAlgInputalpha_ddnD1);
        const REAL_SIMD_ARRAY betU_dupD00 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputbetU_dupD00, UpwindAlgInputbetU_ddnD00), UpwindAlgInputbetU_ddnD00);
        const REAL_SIMD_ARRAY betU_dupD01 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputbetU_dupD01, UpwindAlgInputbetU_ddnD01), UpwindAlgInputbetU_ddnD01);
        const REAL_SIMD_ARRAY betU_dupD10 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputbetU_dupD10, UpwindAlgInputbetU_ddnD10), UpwindAlgInputbetU_ddnD10);
        const REAL_SIMD_ARRAY betU_dupD11 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputbetU_dupD11, UpwindAlgInputbetU_ddnD11), UpwindAlgInputbetU_ddnD11);
        const REAL_SIMD_ARRAY betU_dupD20 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputbetU_dupD20, UpwindAlgInputbetU_ddnD20), UpwindAlgInputbetU_ddnD20);
        const REAL_SIMD_ARRAY betU_dupD21 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputbetU_dupD21, UpwindAlgInputbetU_ddnD21), UpwindAlgInputbetU_ddnD21);
        const REAL_SIMD_ARRAY cf_dupD0 = FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputcf_dupD0, UpwindAlgInputcf_ddnD0), UpwindAlgInputcf_ddnD0);
        const REAL_SIMD_ARRAY cf_dupD1 = FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputcf_dupD1, UpwindAlgInputcf_ddnD1), UpwindAlgInputcf_ddnD1);
        const REAL_SIMD_ARRAY hDD_dupD000 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD000, UpwindAlgInputhDD_ddnD000), UpwindAlgInputhDD_ddnD000);
        const REAL_SIMD_ARRAY hDD_dupD001 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD001, UpwindAlgInputhDD_ddnD001), UpwindAlgInputhDD_ddnD001);
        const REAL_SIMD_ARRAY hDD_dupD010 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD010, UpwindAlgInputhDD_ddnD010), UpwindAlgInputhDD_ddnD010);
        const REAL_SIMD_ARRAY hDD_dupD011 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD011, UpwindAlgInputhDD_ddnD011), UpwindAlgInputhDD_ddnD011);
        const REAL_SIMD_ARRAY hDD_dupD020 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD020, UpwindAlgInputhDD_ddnD020), UpwindAlgInputhDD_ddnD020);
        const REAL_SIMD_ARRAY hDD_dupD021 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD021, UpwindAlgInputhDD_ddnD021), UpwindAlgInputhDD_ddnD021);
        const REAL_SIMD_ARRAY hDD_dupD110 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD110, UpwindAlgInputhDD_ddnD110), UpwindAlgInputhDD_ddnD110);
        const REAL_SIMD_ARRAY hDD_dupD111 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD111, UpwindAlgInputhDD_ddnD111), UpwindAlgInputhDD_ddnD111);
        const REAL_SIMD_ARRAY hDD_dupD120 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD120, UpwindAlgInputhDD_ddnD120), UpwindAlgInputhDD_ddnD120);
        const REAL_SIMD_ARRAY hDD_dupD121 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD121, UpwindAlgInputhDD_ddnD121), UpwindAlgInputhDD_ddnD121);
        const REAL_SIMD_ARRAY hDD_dupD220 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD220, UpwindAlgInputhDD_ddnD220), UpwindAlgInputhDD_ddnD220);
        const REAL_SIMD_ARRAY hDD_dupD221 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD221, UpwindAlgInputhDD_ddnD221), UpwindAlgInputhDD_ddnD221);
        const REAL_SIMD_ARRAY lambdaU_dupD00 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputlambdaU_dupD00, UpwindAlgInputlambdaU_ddnD00), UpwindAlgInputlambdaU_ddnD00);
        const REAL_SIMD_ARRAY lambdaU_dupD01 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputlambdaU_dupD01, UpwindAlgInputlambdaU_ddnD01), UpwindAlgInputlambdaU_ddnD01);
        const REAL_SIMD_ARRAY lambdaU_dupD10 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputlambdaU_dupD10, UpwindAlgInputlambdaU_ddnD10), UpwindAlgInputlambdaU_ddnD10);
        const REAL_SIMD_ARRAY lambdaU_dupD11 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputlambdaU_dupD11, UpwindAlgInputlambdaU_ddnD11), UpwindAlgInputlambdaU_ddnD11);
        const REAL_SIMD_ARRAY lambdaU_dupD20 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputlambdaU_dupD20, UpwindAlgInputlambdaU_ddnD20), UpwindAlgInputlambdaU_ddnD20);
        const REAL_SIMD_ARRAY lambdaU_dupD21 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputlambdaU_dupD21, UpwindAlgInputlambdaU_ddnD21), UpwindAlgInputlambdaU_ddnD21);
        const REAL_SIMD_ARRAY trK_dupD0 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputtrK_dupD0, UpwindAlgInputtrK_ddnD0), UpwindAlgInputtrK_ddnD0);
        const REAL_SIMD_ARRAY trK_dupD1 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputtrK_dupD1, UpwindAlgInputtrK_ddnD1), UpwindAlgInputtrK_ddnD1);
        const REAL_SIMD_ARRAY vetU_dupD00 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputvetU_dupD00, UpwindAlgInputvetU_ddnD00), UpwindAlgInputvetU_ddnD00);
        const REAL_SIMD_ARRAY vetU_dupD01 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputvetU_dupD01, UpwindAlgInputvetU_ddnD01), UpwindAlgInputvetU_ddnD01);
        const REAL_SIMD_ARRAY vetU_dupD10 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputvetU_dupD10, UpwindAlgInputvetU_ddnD10), UpwindAlgInputvetU_ddnD10);
        const REAL_SIMD_ARRAY vetU_dupD11 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputvetU_dupD11, UpwindAlgInputvetU_ddnD11), UpwindAlgInputvetU_ddnD11);
        const REAL_SIMD_ARRAY vetU_dupD20 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputvetU_dupD20, UpwindAlgInputvetU_ddnD20), UpwindAlgInputvetU_ddnD20);
        const REAL_SIMD_ARRAY vetU_dupD21 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputvetU_dupD21, UpwindAlgInputvetU_ddnD21), UpwindAlgInputvetU_ddnD21);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 3 of 3:
         * Evaluate SymPy expressions and write to main memory.
         */
        const double dblFDPart3_Integer_1 = 1.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(dblFDPart3_Integer_1);

        const double dblFDPart3_Integer_12 = 12.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_12 = ConstSIMD(dblFDPart3_Integer_12);

        const double dblFDPart3_Integer_2 = 2.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(dblFDPart3_Integer_2);

        const double dblFDPart3_Integer_20 = 20.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_20 = ConstSIMD(dblFDPart3_Integer_20);

        const double dblFDPart3_Integer_3 = 3.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_3 = ConstSIMD(dblFDPart3_Integer_3);

        const double dblFDPart3_Integer_4 = 4.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_4 = ConstSIMD(dblFDPart3_Integer_4);

        const double dblFDPart3_Integer_8 = 8.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_8 = ConstSIMD(dblFDPart3_Integer_8);

        const double dblFDPart3_NegativeOne_ = -1.0;
        const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(dblFDPart3_NegativeOne_);

        const double dblFDPart3_Rational_1_12 = 1.0 / 12.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_12 = ConstSIMD(dblFDPart3_Rational_1_12);

        const double dblFDPart3_Rational_1_2 = 1.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_2 = ConstSIMD(dblFDPart3_Rational_1_2);

        const double dblFDPart3_Rational_1_3 = 1.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_3 = ConstSIMD(dblFDPart3_Rational_1_3);

        const double dblFDPart3_Rational_1_4 = 1.0 / 4.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_4 = ConstSIMD(dblFDPart3_Rational_1_4);

        const double dblFDPart3_Rational_1_6 = 1.0 / 6.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_6 = ConstSIMD(dblFDPart3_Rational_1_6);

        const double dblFDPart3_Rational_2_3 = 2.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_2_3 = ConstSIMD(dblFDPart3_Rational_2_3);

        const double dblFDPart3_Rational_3_2 = 3.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_3_2 = ConstSIMD(dblFDPart3_Rational_3_2);

        const double dblFDPart3_Rational_3_4 = 3.0 / 4.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_3_4 = ConstSIMD(dblFDPart3_Rational_3_4);

        const double dblFDPart3_Rational_4_3 = 4.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_4_3 = ConstSIMD(dblFDPart3_Rational_4_3);

        const REAL_SIMD_ARRAY FDPart3tmp0 = DivSIMD(FDPart3_Integer_1, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp3 = DivSIMD(FDPart3_Integer_1, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp5 = MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp8 = MulSIMD(alpha, trK);
        const REAL_SIMD_ARRAY FDPart3tmp11 = MulSIMD(FDPart3_Integer_2, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp14 = MulSIMD(f0_of_xx0__D0, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp15 = MulSIMD(f0_of_xx0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp28 = DivSIMD(FDPart3_Integer_1, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp41 = MulSIMD(f1_of_xx1, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp43 = MulSIMD(MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0), f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp48 = MulSIMD(MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp49 = MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp74 = MulSIMD(f0_of_xx0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp98 = MulSIMD(cf, cf);
        const REAL_SIMD_ARRAY FDPart3tmp115 = DivSIMD(FDPart3_Integer_1, cf);
        const REAL_SIMD_ARRAY FDPart3tmp118 = MulSIMD(f1_of_xx1, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp166 = MulSIMD(f0_of_xx0, f0_of_xx0__DD00);
        const REAL_SIMD_ARRAY FDPart3tmp182 = MulSIMD(f0_of_xx0__D0, vetU_dD11);
        const REAL_SIMD_ARRAY FDPart3tmp207 = MulSIMD(f0_of_xx0, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp248 = MulSIMD(f1_of_xx1__D1, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp249 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(f1_of_xx1, f1_of_xx1), f1_of_xx1));
        const REAL_SIMD_ARRAY FDPart3tmp285 = MulSIMD(FDPart3_Integer_3, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp308 = MulSIMD(FDPart3_Rational_3_2, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp423 = MulSIMD(FDPart3_Integer_2, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp427 = MulSIMD(FDPart3_Integer_4, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp2 = MulSIMD(FDPart3tmp0, MulSIMD(cf, KreissOliger_strength_nongauge));
        const REAL_SIMD_ARRAY FDPart3tmp4 = MulSIMD(FDPart3tmp3, MulSIMD(cf, KreissOliger_strength_nongauge));
        const REAL_SIMD_ARRAY FDPart3tmp6 = DivSIMD(FDPart3_Integer_1, FDPart3tmp5);
        const REAL_SIMD_ARRAY FDPart3tmp7 = MulSIMD(FDPart3tmp5, aDD00);
        const REAL_SIMD_ARRAY FDPart3tmp9 = MulSIMD(FDPart3tmp3, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp12 = MulSIMD(FDPart3tmp11, f0_of_xx0__DD00);
        const REAL_SIMD_ARRAY FDPart3tmp13 = MulSIMD(FDPart3tmp0, vetU0);
        const REAL_SIMD_ARRAY FDPart3tmp16 = DivSIMD(FDPart3_Integer_1, FDPart3tmp15);
        const REAL_SIMD_ARRAY FDPart3tmp20 = MulSIMD(FDPart3tmp11, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp22 = MulSIMD(FDPart3tmp0, vetU_dD00);
        const REAL_SIMD_ARRAY FDPart3tmp29 = MulSIMD(FDPart3tmp28, FDPart3tmp3);
        const REAL_SIMD_ARRAY FDPart3tmp38 = MulSIMD(FDPart3tmp3, vetU_dD11);
        const REAL_SIMD_ARRAY FDPart3tmp39 = MulSIMD(FDPart3tmp28, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp42 = DivSIMD(FDPart3_Integer_1, FDPart3tmp41);
        const REAL_SIMD_ARRAY FDPart3tmp44 = DivSIMD(FDPart3_Integer_1, FDPart3tmp43);
        const REAL_SIMD_ARRAY FDPart3tmp47 = MulSIMD(FDPart3tmp41, FDPart3tmp43);
        const REAL_SIMD_ARRAY FDPart3tmp50 = MulSIMD(FDPart3tmp41, FDPart3tmp49);
        const REAL_SIMD_ARRAY FDPart3tmp52 = DivSIMD(FDPart3_Integer_1, FDPart3tmp48);
        const REAL_SIMD_ARRAY FDPart3tmp58 = MulSIMD(FDPart3tmp5, hDD00);
        const REAL_SIMD_ARRAY FDPart3tmp62 = MulSIMD(FDPart3tmp15, hDD11);
        const REAL_SIMD_ARRAY FDPart3tmp64 = MulSIMD(FDPart3tmp15, FDPart3tmp41);
        const REAL_SIMD_ARRAY FDPart3tmp75 = MulSIMD(FDPart3tmp74, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp79 = MulSIMD(FDPart3tmp5, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp80 = MulSIMD(FDPart3tmp15, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp84 = MulSIMD(FDPart3tmp74, aDD01);
        const REAL_SIMD_ARRAY FDPart3tmp119 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp118, FDPart3tmp15));
        const REAL_SIMD_ARRAY FDPart3tmp136 = MulSIMD(FDPart3tmp115, cf_dD1);
        const REAL_SIMD_ARRAY FDPart3tmp139 = MulSIMD(FDPart3tmp115, cf_dD0);
        const REAL_SIMD_ARRAY FDPart3tmp142 = DivSIMD(alpha, FDPart3tmp98);
        const REAL_SIMD_ARRAY FDPart3tmp151 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp115);
        const REAL_SIMD_ARRAY FDPart3tmp184 = MulSIMD(FDPart3tmp15, aDD11);
        const REAL_SIMD_ARRAY FDPart3tmp185 = MulSIMD(FDPart3tmp3, vetU2);
        const REAL_SIMD_ARRAY FDPart3tmp205 = MulSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, vetU_dD01));
        const REAL_SIMD_ARRAY FDPart3tmp212 = MulSIMD(FDPart3tmp0, MulSIMD(cf, KreissOliger_strength_gauge));
        const REAL_SIMD_ARRAY FDPart3tmp213 = MulSIMD(FDPart3tmp3, MulSIMD(cf, KreissOliger_strength_gauge));
        const REAL_SIMD_ARRAY FDPart3tmp214 = MulSIMD(FDPart3tmp0, betU0);
        const REAL_SIMD_ARRAY FDPart3tmp220 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(alpha, trK_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp221 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(alpha, trK_dD0));
        const REAL_SIMD_ARRAY FDPart3tmp224 = MulSIMD(FDPart3tmp3, vetU_dD01);
        const REAL_SIMD_ARRAY FDPart3tmp226 = MulSIMD(FDPart3tmp0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp234 = MulSIMD(FDPart3tmp0, vetU_dD01);
        const REAL_SIMD_ARRAY FDPart3tmp244 = MulSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0__DD00, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp254 =
            FusedMulAddSIMD(MulSIMD(FDPart3_Integer_8, FDPart3tmp118), MulSIMD(FDPart3tmp48, FDPart3tmp49),
                            MulSIMD(FDPart3tmp43, MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3tmp118), MulSIMD(f0_of_xx0__D0, f0_of_xx0__DD00))));
        const REAL_SIMD_ARRAY FDPart3tmp274 = MulSIMD(FDPart3tmp15, hDD_dD111);
        const REAL_SIMD_ARRAY FDPart3tmp278 = MulSIMD(FDPart3tmp5, hDD_dD001);
        const REAL_SIMD_ARRAY FDPart3tmp281 = MulSIMD(FDPart3tmp3, betU1);
        const REAL_SIMD_ARRAY FDPart3tmp367 = DivSIMD(FDPart3tmp5, FDPart3tmp49);
        const REAL_SIMD_ARRAY FDPart3tmp419 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp0, lambdaU0));
        const REAL_SIMD_ARRAY FDPart3tmp30 = MulSIMD(FDPart3tmp16, FDPart3tmp28);
        const REAL_SIMD_ARRAY FDPart3tmp35 = MulSIMD(FDPart3tmp20, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp51 =
            FusedMulAddSIMD(FDPart3tmp12, FDPart3tmp47, MulSIMD(FDPart3_Integer_4, MulSIMD(FDPart3tmp48, FDPart3tmp50)));
        const REAL_SIMD_ARRAY FDPart3tmp59 = AddSIMD(FDPart3tmp5, FDPart3tmp58);
        const REAL_SIMD_ARRAY FDPart3tmp61 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp41), MulSIMD(FDPart3tmp43, MulSIMD(hDD12, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp63 = AddSIMD(FDPart3tmp15, FDPart3tmp62);
        const REAL_SIMD_ARRAY FDPart3tmp66 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp5), MulSIMD(FDPart3tmp64, MulSIMD(hDD02, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp67 = MulSIMD(FDPart3tmp64, hDD22);
        const REAL_SIMD_ARRAY FDPart3tmp69 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp15), MulSIMD(FDPart3tmp5, MulSIMD(hDD01, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp76 = MulSIMD(FDPart3tmp75, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp81 = MulSIMD(FDPart3tmp80, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp87 = MulSIMD(FDPart3tmp75, aDD02);
        const REAL_SIMD_ARRAY FDPart3tmp101 = MulSIMD(FDPart3tmp20, FDPart3tmp41);
        const REAL_SIMD_ARRAY FDPart3tmp106 = MulSIMD(FDPart3tmp74, MulSIMD(f1_of_xx1__D1, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp120 = FusedMulAddSIMD(FDPart3tmp119, hDD22, FDPart3tmp119);
        const REAL_SIMD_ARRAY FDPart3tmp138 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp136, alpha_dD0));
        const REAL_SIMD_ARRAY FDPart3tmp141 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp139, alpha_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp144 = FusedMulAddSIMD(FDPart3tmp20, hDD11, FDPart3tmp20);
        const REAL_SIMD_ARRAY FDPart3tmp155 = MulSIMD(FDPart3_Integer_2, FDPart3tmp80);
        const REAL_SIMD_ARRAY FDPart3tmp157 = MulSIMD(FDPart3tmp15, MulSIMD(f1_of_xx1__D1, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp159 =
            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, f0_of_xx0__D0),
                               FusedMulAddSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, hDD11)),
                                               FusedMulSubSIMD(FDPart3tmp20, hDD_dD011, MulSIMD(FDPart3tmp15, hDD_dD110))));
        const REAL_SIMD_ARRAY FDPart3tmp167 = FusedMulAddSIMD(FDPart3tmp166, f1_of_xx1, FDPart3tmp79);
        const REAL_SIMD_ARRAY FDPart3tmp170 = FusedMulAddSIMD(FDPart3tmp12, hDD00, FDPart3tmp12);
        const REAL_SIMD_ARRAY FDPart3tmp173 = MulSIMD(hDD01, AddSIMD(FDPart3tmp166, FDPart3tmp5));
        const REAL_SIMD_ARRAY FDPart3tmp187 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp185), MulSIMD(FDPart3tmp42, f1_of_xx1__D1));
        const REAL_SIMD_ARRAY FDPart3tmp189 = MulSIMD(FDPart3tmp80, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp197 = MulSIMD(FDPart3tmp64, aDD22);
        const REAL_SIMD_ARRAY FDPart3tmp210 = MulSIMD(FDPart3tmp16, FDPart3tmp42);
        const REAL_SIMD_ARRAY FDPart3tmp216 = MulSIMD(FDPart3tmp0, FDPart3tmp9);
        const REAL_SIMD_ARRAY FDPart3tmp235 = MulSIMD(FDPart3tmp226, FDPart3tmp41);
        const REAL_SIMD_ARRAY FDPart3tmp238 = MulSIMD(FDPart3tmp29, vetU2);
        const REAL_SIMD_ARRAY FDPart3tmp261 = MulSIMD(FDPart3tmp29, betU2);
        const REAL_SIMD_ARRAY FDPart3tmp276 = MulSIMD(FDPart3tmp16, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp280 = MulSIMD(FDPart3tmp214, FDPart3tmp9);
        const REAL_SIMD_ARRAY FDPart3tmp282 = MulSIMD(FDPart3tmp13, FDPart3tmp281);
        const REAL_SIMD_ARRAY FDPart3tmp284 = MulSIMD(FDPart3tmp6, MulSIMD(betU0, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp286 = FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp285, alpha_dD0);
        const REAL_SIMD_ARRAY FDPart3tmp299 = MulSIMD(FDPart3tmp20, aDD01);
        const REAL_SIMD_ARRAY FDPart3tmp310 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp0, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp323 = FusedMulAddSIMD(FDPart3tmp136, FDPart3tmp285, alpha_dD1);
        const REAL_SIMD_ARRAY FDPart3tmp359 = FusedMulSubSIMD(FDPart3tmp3, f0_of_xx0__DD00, MulSIMD(FDPart3tmp16, FDPart3tmp5));
        const REAL_SIMD_ARRAY FDPart3tmp383 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp3, f0_of_xx0__D0));
        const REAL_SIMD_ARRAY FDPart3tmp388 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp248, FDPart3tmp42));
        const REAL_SIMD_ARRAY FDPart3tmp389 = FusedMulAddSIMD(FDPart3tmp38, FDPart3tmp39, FDPart3tmp224);
        const REAL_SIMD_ARRAY FDPart3tmp437 = MulSIMD(FDPart3tmp16, MulSIMD(vetU1, vetU1));
        const REAL_SIMD_ARRAY FDPart3tmp439 = MulSIMD(FDPart3tmp6, MulSIMD(vetU0, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp18 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp14, FDPart3tmp16));
        const REAL_SIMD_ARRAY FDPart3tmp25 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp6), MulSIMD(f0_of_xx0__DD00, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp31 = MulSIMD(FDPart3tmp30, vetU2);
        const REAL_SIMD_ARRAY FDPart3tmp46 = MulSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp42, FDPart3tmp44));
        const REAL_SIMD_ARRAY FDPart3tmp54 = MulSIMD(FDPart3tmp51, MulSIMD(FDPart3tmp52, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp68 = AddSIMD(FDPart3tmp64, FDPart3tmp67);
        const REAL_SIMD_ARRAY FDPart3tmp77 =
            FusedMulSubSIMD(f1_of_xx1, MulSIMD(MulSIMD(FDPart3tmp49, f0_of_xx0__D0), MulSIMD(hDD01, hDD12)), MulSIMD(FDPart3tmp63, FDPart3tmp76));
        const REAL_SIMD_ARRAY FDPart3tmp82 =
            FusedMulSubSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp79), MulSIMD(hDD01, hDD02), MulSIMD(FDPart3tmp59, FDPart3tmp81));
        const REAL_SIMD_ARRAY FDPart3tmp85 = FusedMulAddSIMD(FDPart3tmp59, FDPart3tmp63, FDPart3tmp69);
        const REAL_SIMD_ARRAY FDPart3tmp102 = FusedMulAddSIMD(FDPart3tmp101, hDD22, FDPart3tmp101);
        const REAL_SIMD_ARRAY FDPart3tmp121 = FusedMulAddSIMD(FDPart3tmp64, hDD_dD221, FDPart3tmp120);
        const REAL_SIMD_ARRAY FDPart3tmp122 = FusedMulAddSIMD(FDPart3tmp75, hDD_dD021, FDPart3tmp106);
        const REAL_SIMD_ARRAY FDPart3tmp129 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp15, f1_of_xx1), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1__D1, hDD22)),
            FusedMulSubSIMD(FDPart3tmp15, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, f1_of_xx1__D1)),
                            MulSIMD(FDPart3tmp64, hDD_dD221)));
        const REAL_SIMD_ARRAY FDPart3tmp130 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp41, f0_of_xx0), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, hDD22)),
            FusedMulSubSIMD(FDPart3tmp41, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0)),
                            MulSIMD(FDPart3tmp64, hDD_dD220)));
        const REAL_SIMD_ARRAY FDPart3tmp145 = FusedMulAddSIMD(FDPart3tmp15, hDD_dD110, FDPart3tmp144);
        const REAL_SIMD_ARRAY FDPart3tmp158 = FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp157, MulSIMD(FDPart3tmp155, hDD_dD121));
        const REAL_SIMD_ARRAY FDPart3tmp171 = FusedMulAddSIMD(FDPart3tmp5, hDD_dD000, FDPart3tmp170);
        const REAL_SIMD_ARRAY FDPart3tmp174 =
            FusedMulAddSIMD(FDPart3tmp20, hDD_dD010, FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp173, MulSIMD(FDPart3tmp5, hDD_dD001)));
        const REAL_SIMD_ARRAY FDPart3tmp188 = FusedMulAddSIMD(FDPart3tmp29, vetU_dD21, FDPart3tmp187);
        const REAL_SIMD_ARRAY FDPart3tmp218 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp6, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp230 = MulSIMD(FDPart3tmp166, MulSIMD(FDPart3tmp52, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp245 = FusedMulAddSIMD(
            vetU0, FusedMulSubSIMD(FDPart3tmp244, FDPart3tmp52, MulSIMD(FDPart3tmp6, f0_of_xx0__DDD000)),
            FusedMulSubSIMD(FDPart3tmp0, vetU_dDD000, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp6), MulSIMD(f0_of_xx0__DD00, vetU_dD00))));
        const REAL_SIMD_ARRAY FDPart3tmp251 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp43, FDPart3tmp79));
        const REAL_SIMD_ARRAY FDPart3tmp260 = MulSIMD(FDPart3tmp214, FDPart3tmp238);
        const REAL_SIMD_ARRAY FDPart3tmp262 = MulSIMD(FDPart3tmp13, FDPart3tmp261);
        const REAL_SIMD_ARRAY FDPart3tmp266 = MulSIMD(FDPart3tmp30, betU2);
        const REAL_SIMD_ARRAY FDPart3tmp271 = FusedMulSubSIMD(FDPart3tmp3, vetU_dDD101, MulSIMD(FDPart3tmp16, FDPart3tmp182));
        const REAL_SIMD_ARRAY FDPart3tmp277 = MulSIMD(FDPart3tmp276, betU1);
        const REAL_SIMD_ARRAY FDPart3tmp301 = MulSIMD(FDPart3tmp35, aDD02);
        const REAL_SIMD_ARRAY FDPart3tmp304 = MulSIMD(FDPart3tmp155, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp350 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp16, f0_of_xx0__D0));
        const REAL_SIMD_ARRAY FDPart3tmp387 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp3), MulSIMD(FDPart3tmp42, f1_of_xx1__D1));
        const REAL_SIMD_ARRAY FDPart3tmp391 = FusedMulAddSIMD(FDPart3tmp28, f1_of_xx1__DD11, FDPart3tmp388);
        const REAL_SIMD_ARRAY FDPart3tmp397 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp30, f0_of_xx0__D0));
        const REAL_SIMD_ARRAY FDPart3tmp399 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp30, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp436 = MulSIMD(FDPart3tmp210, MulSIMD(vetU2, vetU2));
        const REAL_SIMD_ARRAY FDPart3tmp441 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp13, FDPart3tmp238));
        const REAL_SIMD_ARRAY FDPart3tmp442 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp13, FDPart3tmp9));
        const REAL_SIMD_ARRAY FDPart3tmp19 = FusedMulAddSIMD(FDPart3tmp3, vetU_dD10, FDPart3tmp18);
        const REAL_SIMD_ARRAY FDPart3tmp26 = AddSIMD(FDPart3tmp22, FDPart3tmp25);
        const REAL_SIMD_ARRAY FDPart3tmp33 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp31, f0_of_xx0__D0));
        const REAL_SIMD_ARRAY FDPart3tmp90 =
            FusedMulSubSIMD(MulSIMD(FDPart3tmp50, f0_of_xx0__D0), MulSIMD(hDD02, hDD12), MulSIMD(FDPart3tmp68, MulSIMD(FDPart3tmp74, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp96 = FusedMulAddSIMD(FDPart3tmp59, FDPart3tmp68, FDPart3tmp66);
        const REAL_SIMD_ARRAY FDPart3tmp103 = FusedMulAddSIMD(FDPart3tmp64, hDD_dD220, FDPart3tmp102);
        const REAL_SIMD_ARRAY FDPart3tmp109 = FusedMulAddSIMD(FDPart3tmp35, hDD12, MulSIMD(FDPart3tmp80, hDD_dD120));
        const REAL_SIMD_ARRAY FDPart3tmp123 =
            FusedMulAddSIMD(MulSIMD(f0_of_xx0, f0_of_xx0__D0), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, hDD12)),
                            NegFusedMulAddSIMD(FDPart3tmp80, hDD_dD120, FDPart3tmp122));
        const REAL_SIMD_ARRAY FDPart3tmp169 = FusedMulAddSIMD(FDPart3tmp35, hDD_dD020, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp167, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp253 = MulSIMD(FDPart3tmp46, FDPart3tmp6);
        const REAL_SIMD_ARRAY FDPart3tmp265 = MulSIMD(FDPart3tmp31, betU1);
        const REAL_SIMD_ARRAY FDPart3tmp270 = MulSIMD(FDPart3tmp210, MulSIMD(betU2, vetU2));
        const REAL_SIMD_ARRAY FDPart3tmp287 = MulSIMD(FDPart3tmp77, FDPart3tmp77);
        const REAL_SIMD_ARRAY FDPart3tmp295 = MulSIMD(FDPart3tmp82, FDPart3tmp82);
        const REAL_SIMD_ARRAY FDPart3tmp298 = MulSIMD(FDPart3tmp77, FDPart3tmp82);
        const REAL_SIMD_ARRAY FDPart3tmp34 = FusedMulAddSIMD(FDPart3tmp29, vetU_dD20, FDPart3tmp33);
        const REAL_SIMD_ARRAY FDPart3tmp55 =
            FusedMulAddSIMD(FDPart3tmp39, FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp46, FDPart3tmp54, AddSIMD(FDPart3tmp26, FDPart3tmp38)));
        const REAL_SIMD_ARRAY FDPart3tmp71 = FusedMulAddSIMD(
            FDPart3tmp68, FDPart3tmp69,
            FusedMulAddSIMD(FDPart3tmp59, MulSIMD(FDPart3tmp63, FDPart3tmp68),
                            FusedMulAddSIMD(MulSIMD(FDPart3tmp5, hDD01), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp47), MulSIMD(hDD02, hDD12)),
                                            FusedMulAddSIMD(FDPart3tmp59, FDPart3tmp61, MulSIMD(FDPart3tmp63, FDPart3tmp66)))));
        const REAL_SIMD_ARRAY FDPart3tmp93 = FusedMulAddSIMD(FDPart3tmp63, FDPart3tmp68, FDPart3tmp61);
        const REAL_SIMD_ARRAY FDPart3tmp110 = SubSIMD(SubSIMD(FDPart3tmp109, FDPart3tmp106), MulSIMD(FDPart3tmp75, hDD_dD021));
        const REAL_SIMD_ARRAY FDPart3tmp146 = AddSIMD(FDPart3tmp109, FDPart3tmp122);
        const REAL_SIMD_ARRAY FDPart3tmp255 = FusedMulAddSIMD(
            FDPart3tmp253,
            FusedMulAddSIMD(FDPart3tmp9,
                            FusedMulAddSIMD(FDPart3tmp251, f1_of_xx1__DD11,
                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp248), MulSIMD(FDPart3tmp43, FDPart3tmp5))),
                            MulSIMD(FDPart3tmp79, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp49), MulSIMD(f1_of_xx1__D1, vetU_dD11)))),
            FusedMulAddSIMD(
                FDPart3tmp3, vetU_dDD111,
                FusedMulAddSIMD(
                    FDPart3tmp218, vetU_dD01,
                    FusedMulAddSIMD(
                        FDPart3tmp253, FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp254, MulSIMD(FDPart3tmp234, FDPart3tmp51)),
                        FusedMulAddSIMD(MulSIMD(FDPart3tmp248, FDPart3tmp3),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp42, vetU1)),
                                        FusedMulSubSIMD(FDPart3tmp0, vetU_dDD001,
                                                        MulSIMD(MulSIMD(FDPart3tmp249, FDPart3tmp44), MulSIMD(FDPart3tmp54, f1_of_xx1__D1))))))));
        const REAL_SIMD_ARRAY FDPart3tmp272 = FusedMulAddSIMD(
            FDPart3tmp253, FusedMulAddSIMD(FDPart3tmp254, FDPart3tmp9, MulSIMD(FDPart3tmp19, MulSIMD(FDPart3tmp251, f1_of_xx1__D1))),
            AddSIMD(
                FusedMulAddSIMD(
                    FDPart3tmp253,
                    FusedMulAddSIMD(
                        FDPart3tmp13,
                        FusedMulAddSIMD(
                            FDPart3tmp11, MulSIMD(FDPart3tmp47, f0_of_xx0__DDD000),
                            FusedMulAddSIMD(
                                MulSIMD(FDPart3_Integer_20, FDPart3tmp5), MulSIMD(FDPart3tmp50, f0_of_xx0__DD00),
                                FusedMulAddSIMD(FDPart3tmp244, FDPart3tmp47,
                                                MulSIMD(FDPart3_Integer_12,
                                                        MulSIMD(FDPart3tmp64, MulSIMD(MulSIMD(MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0),
                                                                                      f0_of_xx0__D0)))))),
                        MulSIMD(FDPart3tmp26, FDPart3tmp51)),
                    FDPart3tmp271),
                FusedMulAddSIMD(
                    f1_of_xx1__D1,
                    MulSIMD(MulSIMD(FDPart3tmp51, FDPart3tmp6),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp249),
                                    DivSIMD(vetU1, MulSIMD(MulSIMD(MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0), f0_of_xx0), f0_of_xx0)))),
                    NegFusedMulAddSIMD(
                        DivSIMD(
                            DivSIMD(FDPart3_Integer_1,
                                    MulSIMD(MulSIMD(MulSIMD(MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0), f0_of_xx0__D0), f0_of_xx0__D0)),
                            MulSIMD(MulSIMD(MulSIMD(f1_of_xx1, f1_of_xx1), f1_of_xx1), f1_of_xx1)),
                        MulSIMD(
                            MulSIMD(FDPart3_Rational_1_2, vetU0),
                            DivSIMD(MulSIMD(FDPart3tmp51, FDPart3tmp51),
                                    MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0), f0_of_xx0), f0_of_xx0),
                                                            f0_of_xx0),
                                                    f0_of_xx0),
                                            f0_of_xx0))),
                        FDPart3tmp245))));
        const REAL_SIMD_ARRAY FDPart3tmp289 = MulSIMD(FDPart3tmp90, FDPart3tmp90);
        const REAL_SIMD_ARRAY FDPart3tmp355 = MulSIMD(FDPart3tmp26, MulSIMD(FDPart3tmp3, f0_of_xx0__D0));
        const REAL_SIMD_ARRAY FDPart3tmp364 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp118, FDPart3tmp188));
        const REAL_SIMD_ARRAY FDPart3tmp440 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp31, vetU1));
        const REAL_SIMD_ARRAY FDPart3tmp72 = DivSIMD(FDPart3_Integer_1, FDPart3tmp71);
        const REAL_SIMD_ARRAY FDPart3tmp228 = SubSIMD(FusedMulSubSIMD(FDPart3tmp0, vetU_dDD001, FDPart3tmp224), MulSIMD(FDPart3tmp19, FDPart3tmp226));
        const REAL_SIMD_ARRAY FDPart3tmp256 = MulSIMD(FDPart3_Rational_1_4, FDPart3tmp255);
        const REAL_SIMD_ARRAY FDPart3tmp273 = MulSIMD(FDPart3_Rational_1_4, FDPart3tmp272);
        const REAL_SIMD_ARRAY FDPart3tmp288 = DivSIMD(FDPart3_Integer_1, MulSIMD(FDPart3tmp71, FDPart3tmp71));
        const REAL_SIMD_ARRAY FDPart3tmp357 =
            FusedMulAddSIMD(FDPart3tmp350, vetU0, FusedMulAddSIMD(FDPart3tmp13, MulSIMD(FDPart3tmp3, f0_of_xx0__DD00), FDPart3tmp355));
        const REAL_SIMD_ARRAY FDPart3tmp360 = FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp359, FDPart3tmp355);
        const REAL_SIMD_ARRAY FDPart3tmp398 =
            FusedMulAddSIMD(FDPart3tmp387, vetU_dD20,
                            FusedMulAddSIMD(FDPart3tmp397, vetU_dD21,
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp210, f0_of_xx0__D0), MulSIMD(f1_of_xx1__D1, vetU2),
                                                            FusedMulAddSIMD(FDPart3tmp29, vetU_dDD201, MulSIMD(FDPart3tmp34, FDPart3tmp39)))));
        const REAL_SIMD_ARRAY FDPart3tmp421 = MulSIMD(FDPart3_Rational_1_3, FDPart3tmp255);
        const REAL_SIMD_ARRAY FDPart3tmp422 = MulSIMD(FDPart3_Rational_1_3, FDPart3tmp272);
        const REAL_SIMD_ARRAY FDPart3tmp78 = MulSIMD(FDPart3tmp72, FDPart3tmp77);
        const REAL_SIMD_ARRAY FDPart3tmp83 = MulSIMD(FDPart3tmp72, FDPart3tmp82);
        const REAL_SIMD_ARRAY FDPart3tmp86 = MulSIMD(FDPart3tmp72, FDPart3tmp85);
        const REAL_SIMD_ARRAY FDPart3tmp91 = MulSIMD(FDPart3tmp72, FDPart3tmp90);
        const REAL_SIMD_ARRAY FDPart3tmp94 = MulSIMD(FDPart3tmp72, FDPart3tmp93);
        const REAL_SIMD_ARRAY FDPart3tmp97 = MulSIMD(FDPart3tmp72, FDPart3tmp96);
        const REAL_SIMD_ARRAY FDPart3tmp111 =
            FusedMulSubSIMD(FDPart3tmp110, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp90)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp103), MulSIMD(FDPart3tmp72, FDPart3tmp77)));
        const REAL_SIMD_ARRAY FDPart3tmp113 =
            FusedMulSubSIMD(FDPart3tmp110, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp96)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp103), MulSIMD(FDPart3tmp72, FDPart3tmp82)));
        const REAL_SIMD_ARRAY FDPart3tmp124 =
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp90)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp72, FDPart3tmp82)));
        const REAL_SIMD_ARRAY FDPart3tmp126 =
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp93)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp72, FDPart3tmp77)));
        const REAL_SIMD_ARRAY FDPart3tmp131 =
            FusedMulSubSIMD(FDPart3tmp130, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp93)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp129), MulSIMD(FDPart3tmp72, FDPart3tmp90)));
        const REAL_SIMD_ARRAY FDPart3tmp133 =
            FusedMulSubSIMD(FDPart3tmp130, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp90)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp129), MulSIMD(FDPart3tmp72, FDPart3tmp96)));
        const REAL_SIMD_ARRAY FDPart3tmp147 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp5, FDPart3tmp72), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp93, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp146, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp77)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp145), MulSIMD(FDPart3tmp72, FDPart3tmp90))));
        const REAL_SIMD_ARRAY FDPart3tmp149 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp5, FDPart3tmp72), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp90, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp146, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp82)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp145), MulSIMD(FDPart3tmp72, FDPart3tmp96))));
        const REAL_SIMD_ARRAY FDPart3tmp160 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp15, FDPart3tmp72), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp90, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp159, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp93)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp158), MulSIMD(FDPart3tmp72, FDPart3tmp77))));
        const REAL_SIMD_ARRAY FDPart3tmp162 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp15, FDPart3tmp72), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp96, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp159, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp90)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp158), MulSIMD(FDPart3tmp72, FDPart3tmp82))));
        const REAL_SIMD_ARRAY FDPart3tmp175 = FusedMulAddSIMD(
            FDPart3tmp174, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp90)),
            FusedMulSubSIMD(FDPart3tmp171, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp93)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp169), MulSIMD(FDPart3tmp72, FDPart3tmp77))));
        const REAL_SIMD_ARRAY FDPart3tmp177 = FusedMulAddSIMD(
            FDPart3tmp174, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp96)),
            FusedMulSubSIMD(FDPart3tmp171, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp72, FDPart3tmp90)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp169), MulSIMD(FDPart3tmp72, FDPart3tmp82))));
        const REAL_SIMD_ARRAY FDPart3tmp240 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp226), MulSIMD(FDPart3tmp34, FDPart3tmp41));
        const REAL_SIMD_ARRAY FDPart3tmp294 = MulSIMD(FDPart3tmp288, FDPart3tmp7);
        const REAL_SIMD_ARRAY FDPart3tmp296 = MulSIMD(FDPart3tmp184, FDPart3tmp288);
        const REAL_SIMD_ARRAY FDPart3tmp297 = MulSIMD(FDPart3tmp197, FDPart3tmp288);
        const REAL_SIMD_ARRAY FDPart3tmp300 = MulSIMD(FDPart3tmp288, FDPart3tmp299);
        const REAL_SIMD_ARRAY FDPart3tmp324 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp96, aDD01),
            MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp93),
                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
            FusedMulAddSIMD(
                f0_of_xx0,
                MulSIMD(MulSIMD(FDPart3tmp90, aDD02),
                        MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp77),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f1_of_xx1)))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp77, FDPart3tmp96),
                    MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp82, FDPart3tmp90),
                        MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
                        FusedMulAddSIMD(
                            FDPart3tmp90,
                            MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp5),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp93, aDD00))),
                            FusedMulAddSIMD(
                                MulSIMD(FDPart3tmp41, FDPart3tmp77),
                                MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp82, aDD22))),
                                FusedMulAddSIMD(
                                    f0_of_xx0,
                                    MulSIMD(MulSIMD(FDPart3tmp93, aDD02),
                                            MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp82),
                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f1_of_xx1)))),
                                    FusedMulSubSIMD(
                                        aDD01,
                                        MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp289),
                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                        MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp90),
                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp15), MulSIMD(FDPart3tmp96, aDD11)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp326 = MulSIMD(FDPart3tmp288, FDPart3tmp87);
        const REAL_SIMD_ARRAY FDPart3tmp327 = MulSIMD(FDPart3tmp189, FDPart3tmp288);
        const REAL_SIMD_ARRAY FDPart3tmp329 = MulSIMD(FDPart3tmp288, FDPart3tmp84);
        const REAL_SIMD_ARRAY FDPart3tmp190 =
            FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp86, FusedMulAddSIMD(FDPart3tmp78, FDPart3tmp84, MulSIMD(FDPart3tmp184, FDPart3tmp83)));
        const REAL_SIMD_ARRAY FDPart3tmp191 =
            FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp78, FusedMulAddSIMD(FDPart3tmp84, FDPart3tmp94, MulSIMD(FDPart3tmp184, FDPart3tmp91)));
        const REAL_SIMD_ARRAY FDPart3tmp200 =
            FusedMulAddSIMD(FDPart3tmp197, FDPart3tmp78, FusedMulAddSIMD(FDPart3tmp87, FDPart3tmp94, MulSIMD(FDPart3tmp189, FDPart3tmp91)));
        const REAL_SIMD_ARRAY FDPart3tmp201 =
            FusedMulAddSIMD(FDPart3tmp197, FDPart3tmp83, FusedMulAddSIMD(FDPart3tmp87, FDPart3tmp91, MulSIMD(FDPart3tmp189, FDPart3tmp97)));
        const REAL_SIMD_ARRAY FDPart3tmp257 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp78);
        const REAL_SIMD_ARRAY FDPart3tmp258 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp91);
        const REAL_SIMD_ARRAY FDPart3tmp263 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp94);
        const REAL_SIMD_ARRAY FDPart3tmp303 = MulSIMD(FDPart3tmp288, MulSIMD(FDPart3tmp301, FDPart3tmp77));
        const REAL_SIMD_ARRAY FDPart3tmp306 = MulSIMD(FDPart3tmp288, MulSIMD(FDPart3tmp304, FDPart3tmp82));
        const REAL_SIMD_ARRAY FDPart3tmp343 = MulSIMD(FDPart3_Rational_4_3, FDPart3tmp78);
        const REAL_SIMD_ARRAY FDPart3tmp344 = MulSIMD(FDPart3_Rational_4_3, FDPart3tmp83);
        const REAL_SIMD_ARRAY FDPart3tmp345 = MulSIMD(FDPart3_Rational_2_3, FDPart3tmp86);
        const REAL_SIMD_ARRAY FDPart3tmp346 = MulSIMD(FDPart3_Rational_4_3, FDPart3tmp91);
        const REAL_SIMD_ARRAY FDPart3tmp347 = MulSIMD(FDPart3_Rational_2_3, FDPart3tmp97);
        const REAL_SIMD_ARRAY FDPart3tmp348 = MulSIMD(FDPart3_Rational_2_3, FDPart3tmp94);
        const REAL_SIMD_ARRAY FDPart3tmp351 =
            MulSIMD(FDPart3tmp97, FusedMulAddSIMD(FDPart3tmp19, FDPart3tmp226,
                                                  FusedMulAddSIMD(FDPart3tmp3, vetU_dDD111, MulSIMD(FDPart3_Integer_2, FDPart3tmp224))));
        const REAL_SIMD_ARRAY FDPart3tmp358 = MulSIMD(FDPart3tmp91, AddSIMD(FDPart3tmp271, FDPart3tmp357));
        const REAL_SIMD_ARRAY FDPart3tmp361 = MulSIMD(FDPart3tmp91, AddSIMD(FDPart3tmp271, FDPart3tmp360));
        const REAL_SIMD_ARRAY FDPart3tmp371 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp83);
        const REAL_SIMD_ARRAY FDPart3tmp374 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp97);
        const REAL_SIMD_ARRAY FDPart3tmp392 = MulSIMD(FDPart3tmp83, FusedMulAddSIMD(FDPart3tmp391, FDPart3tmp9, FDPart3tmp389));
        const REAL_SIMD_ARRAY FDPart3tmp394 = MulSIMD(FDPart3tmp78, FusedMulAddSIMD(FDPart3tmp19, FDPart3tmp39, FDPart3tmp357));
        const REAL_SIMD_ARRAY FDPart3tmp395 = MulSIMD(FDPart3tmp78, FusedMulAddSIMD(FDPart3tmp19, FDPart3tmp39, FDPart3tmp360));
        const REAL_SIMD_ARRAY FDPart3tmp396 =
            MulSIMD(FDPart3tmp86, FusedMulAddSIMD(FDPart3tmp118, FDPart3tmp188, MulSIMD(FDPart3tmp235, FDPart3tmp34)));
        const REAL_SIMD_ARRAY FDPart3tmp402 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp86);
        const REAL_SIMD_ARRAY FDPart3tmp433 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp91, FDPart3tmp98));
        const REAL_SIMD_ARRAY FDPart3tmp116 =
            SubSIMD(NegFusedMulAddSIMD(
                        FDPart3_Integer_2,
                        MulSIMD(alpha, FusedMulSubSIMD(FDPart3tmp113,
                                                       MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp115, cf_dD1)),
                                                       MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp111), MulSIMD(FDPart3tmp115, cf_dD0)))),
                        FusedMulSubSIMD(RbarDD02, alpha, MulSIMD(FDPart3tmp111, alpha_dD0))),
                    MulSIMD(FDPart3tmp113, alpha_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp128 =
            SubSIMD(NegFusedMulAddSIMD(
                        FDPart3_Integer_2,
                        MulSIMD(alpha, FusedMulSubSIMD(FDPart3tmp115,
                                                       MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp126, cf_dD0)),
                                                       MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp115), MulSIMD(FDPart3tmp124, cf_dD1)))),
                        FusedMulSubSIMD(RbarDD12, alpha, MulSIMD(FDPart3tmp124, alpha_dD1))),
                    MulSIMD(FDPart3tmp126, alpha_dD0));
        const REAL_SIMD_ARRAY FDPart3tmp135 =
            SubSIMD(NegFusedMulAddSIMD(
                        FDPart3_Integer_2,
                        MulSIMD(alpha, FusedMulSubSIMD(FDPart3tmp115,
                                                       MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp133, cf_dD1)),
                                                       MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp115), MulSIMD(FDPart3tmp131, cf_dD0)))),
                        FusedMulSubSIMD(RbarDD22, alpha, MulSIMD(FDPart3tmp131, alpha_dD0))),
                    MulSIMD(FDPart3tmp133, alpha_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp152 = FusedMulAddSIMD(
            FDPart3tmp142, MulSIMD(cf_dD0, cf_dD1),
            SubSIMD(
                SubSIMD(AddSIMD(FusedMulAddSIMD(RbarDD01, alpha, FDPart3tmp141),
                                NegFusedMulAddSIMD(
                                    FDPart3_Integer_2,
                                    MulSIMD(alpha, FusedMulAddSIMD(
                                                       FDPart3tmp115,
                                                       MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp149, cf_dD1)),
                                                       FusedMulSubSIMD(
                                                           FDPart3tmp151, FusedMulSubSIMD(FDPart3tmp136, cf_dD0, cf_dDD01),
                                                           MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp115), MulSIMD(FDPart3tmp147, cf_dD0))))),
                                    SubSIMD(FDPart3tmp138, alpha_dDD01))),
                        MulSIMD(FDPart3tmp149, alpha_dD1)),
                MulSIMD(FDPart3tmp147, alpha_dD0)));
        const REAL_SIMD_ARRAY FDPart3tmp164 = SubSIMD(
            NegFusedMulAddSIMD(
                FDPart3_Integer_2,
                MulSIMD(alpha,
                        FusedMulAddSIMD(FDPart3tmp115, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp162, cf_dD1)),
                                        FusedMulSubSIMD(FDPart3tmp151, FusedMulSubSIMD(FDPart3tmp115, MulSIMD(cf_dD1, cf_dD1), cf_dDD11),
                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp115), MulSIMD(FDPart3tmp160, cf_dD0))))),
                FusedMulAddSIMD(
                    RbarDD11, alpha,
                    SubSIMD(FusedMulAddSIMD(FDPart3tmp115, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha_dD1, cf_dD1)),
                                            FusedMulSubSIMD(FDPart3tmp142, MulSIMD(cf_dD1, cf_dD1), alpha_dDD11)),
                            MulSIMD(FDPart3tmp160, alpha_dD0)))),
            MulSIMD(FDPart3tmp162, alpha_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp179 = SubSIMD(
            NegFusedMulAddSIMD(
                FDPart3_Integer_2,
                MulSIMD(alpha,
                        FusedMulAddSIMD(FDPart3tmp115, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp177, cf_dD1)),
                                        FusedMulSubSIMD(FDPart3tmp151, FusedMulSubSIMD(FDPart3tmp115, MulSIMD(cf_dD0, cf_dD0), cf_dDD00),
                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp115), MulSIMD(FDPart3tmp175, cf_dD0))))),
                FusedMulAddSIMD(
                    RbarDD00, alpha,
                    SubSIMD(FusedMulAddSIMD(FDPart3tmp115, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha_dD0, cf_dD0)),
                                            FusedMulSubSIMD(FDPart3tmp142, MulSIMD(cf_dD0, cf_dD0), alpha_dDD00)),
                            MulSIMD(FDPart3tmp175, alpha_dD0)))),
            MulSIMD(FDPart3tmp177, alpha_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp194 =
            FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp83, FusedMulAddSIMD(FDPart3tmp84, FDPart3tmp91, MulSIMD(FDPart3tmp184, FDPart3tmp97)));
        const REAL_SIMD_ARRAY FDPart3tmp199 =
            FusedMulAddSIMD(FDPart3tmp197, FDPart3tmp86, FusedMulAddSIMD(FDPart3tmp78, FDPart3tmp87, MulSIMD(FDPart3tmp189, FDPart3tmp83)));
        const REAL_SIMD_ARRAY FDPart3tmp259 = FusedMulAddSIMD(FDPart3tmp103, FDPart3tmp257, MulSIMD(FDPart3tmp110, FDPart3tmp258));
        const REAL_SIMD_ARRAY FDPart3tmp264 = FusedMulAddSIMD(FDPart3tmp121, FDPart3tmp257, MulSIMD(FDPart3tmp123, FDPart3tmp263));
        const REAL_SIMD_ARRAY FDPart3tmp268 = FusedMulAddSIMD(FDPart3tmp129, FDPart3tmp258, MulSIMD(FDPart3tmp130, FDPart3tmp263));
        const REAL_SIMD_ARRAY FDPart3tmp275 =
            FusedMulAddSIMD(FDPart3tmp159, FDPart3tmp263, FusedMulAddSIMD(FDPart3tmp258, FDPart3tmp274, MulSIMD(FDPart3tmp158, FDPart3tmp257)));
        const REAL_SIMD_ARRAY FDPart3tmp279 =
            FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp257, FusedMulAddSIMD(FDPart3tmp263, FDPart3tmp278, MulSIMD(FDPart3tmp145, FDPart3tmp258)));
        const REAL_SIMD_ARRAY FDPart3tmp283 =
            FusedMulAddSIMD(FDPart3tmp171, FDPart3tmp263, FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp258, MulSIMD(FDPart3tmp169, FDPart3tmp257)));
        const REAL_SIMD_ARRAY FDPart3tmp307 = FusedMulAddSIMD(
            FDPart3tmp297, MulSIMD(FDPart3tmp85, FDPart3tmp85),
            FusedMulAddSIMD(FDPart3tmp298, FDPart3tmp300,
                            FusedMulAddSIMD(FDPart3tmp303, FDPart3tmp85,
                                            FusedMulAddSIMD(FDPart3tmp306, FDPart3tmp85,
                                                            FusedMulAddSIMD(FDPart3tmp287, FDPart3tmp294, MulSIMD(FDPart3tmp295, FDPart3tmp296))))));
        const REAL_SIMD_ARRAY FDPart3tmp314 = FusedMulAddSIMD(
            FDPart3tmp294, MulSIMD(FDPart3tmp93, FDPart3tmp93),
            FusedMulAddSIMD(FDPart3tmp303, FDPart3tmp93,
                            FusedMulAddSIMD(FDPart3tmp300, MulSIMD(FDPart3tmp90, FDPart3tmp93),
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp304), MulSIMD(FDPart3tmp77, FDPart3tmp90),
                                                            FusedMulAddSIMD(FDPart3tmp287, FDPart3tmp297, MulSIMD(FDPart3tmp289, FDPart3tmp296))))));
        const REAL_SIMD_ARRAY FDPart3tmp321 = FusedMulAddSIMD(
            FDPart3tmp306, FDPart3tmp96,
            FusedMulAddSIMD(FDPart3tmp197, MulSIMD(FDPart3tmp288, FDPart3tmp295),
                            FusedMulAddSIMD(FDPart3tmp300, MulSIMD(FDPart3tmp90, FDPart3tmp96),
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp301), MulSIMD(FDPart3tmp82, FDPart3tmp90),
                                                            FusedMulAddSIMD(FDPart3tmp289, FDPart3tmp294,
                                                                            MulSIMD(FDPart3tmp296, MulSIMD(FDPart3tmp96, FDPart3tmp96)))))));
        const REAL_SIMD_ARRAY FDPart3tmp337 = FusedMulAddSIMD(
            FDPart3tmp327, MulSIMD(FDPart3tmp85, FDPart3tmp90),
            FusedMulAddSIMD(
                FDPart3tmp326, MulSIMD(FDPart3tmp85, FDPart3tmp93),
                FusedMulAddSIMD(
                    FDPart3tmp327, MulSIMD(FDPart3tmp77, FDPart3tmp82),
                    FusedMulAddSIMD(
                        FDPart3tmp296, MulSIMD(FDPart3tmp82, FDPart3tmp90),
                        FusedMulAddSIMD(
                            FDPart3tmp297, MulSIMD(FDPart3tmp77, FDPart3tmp85),
                            FusedMulAddSIMD(FDPart3tmp329, MulSIMD(FDPart3tmp77, FDPart3tmp90),
                                            FusedMulAddSIMD(FDPart3tmp329, MulSIMD(FDPart3tmp82, FDPart3tmp93),
                                                            FusedMulAddSIMD(FDPart3tmp287, FDPart3tmp326,
                                                                            MulSIMD(FDPart3tmp294, MulSIMD(FDPart3tmp77, FDPart3tmp93))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp339 = FusedMulAddSIMD(
            FDPart3tmp327, MulSIMD(FDPart3tmp85, FDPart3tmp96),
            FusedMulAddSIMD(
                FDPart3tmp297, MulSIMD(FDPart3tmp82, FDPart3tmp85),
                FusedMulAddSIMD(
                    FDPart3tmp326, MulSIMD(FDPart3tmp85, FDPart3tmp90),
                    FusedMulAddSIMD(
                        FDPart3tmp294, MulSIMD(FDPart3tmp77, FDPart3tmp90),
                        FusedMulAddSIMD(
                            FDPart3tmp296, MulSIMD(FDPart3tmp82, FDPart3tmp96),
                            FusedMulAddSIMD(FDPart3tmp329, MulSIMD(FDPart3tmp77, FDPart3tmp96),
                                            FusedMulAddSIMD(FDPart3tmp329, MulSIMD(FDPart3tmp82, FDPart3tmp90),
                                                            FusedMulAddSIMD(FDPart3tmp298, FDPart3tmp326,
                                                                            MulSIMD(FDPart3tmp189, MulSIMD(FDPart3tmp288, FDPart3tmp295))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp341 = FusedMulAddSIMD(
            FDPart3tmp327, MulSIMD(FDPart3tmp77, FDPart3tmp96),
            FusedMulAddSIMD(
                FDPart3tmp326, MulSIMD(FDPart3tmp77, FDPart3tmp90),
                FusedMulAddSIMD(
                    FDPart3tmp326, MulSIMD(FDPart3tmp82, FDPart3tmp93),
                    FusedMulAddSIMD(FDPart3tmp294, MulSIMD(FDPart3tmp90, FDPart3tmp93),
                                    FusedMulAddSIMD(FDPart3tmp296, MulSIMD(FDPart3tmp90, FDPart3tmp96),
                                                    FusedMulAddSIMD(FDPart3tmp327, MulSIMD(FDPart3tmp82, FDPart3tmp90),
                                                                    FusedMulAddSIMD(FDPart3tmp329, MulSIMD(FDPart3tmp93, FDPart3tmp96),
                                                                                    FusedMulAddSIMD(FDPart3tmp289, FDPart3tmp329,
                                                                                                    MulSIMD(FDPart3tmp297, FDPart3tmp298)))))))));
        const REAL_SIMD_ARRAY FDPart3tmp372 = FusedMulAddSIMD(FDPart3tmp121, FDPart3tmp371, MulSIMD(FDPart3tmp123, FDPart3tmp258));
        const REAL_SIMD_ARRAY FDPart3tmp375 = FusedMulAddSIMD(FDPart3tmp103, FDPart3tmp371, MulSIMD(FDPart3tmp110, FDPart3tmp374));
        const REAL_SIMD_ARRAY FDPart3tmp376 = FusedMulAddSIMD(FDPart3tmp129, FDPart3tmp374, MulSIMD(FDPart3tmp130, FDPart3tmp258));
        const REAL_SIMD_ARRAY FDPart3tmp377 =
            FusedMulAddSIMD(FDPart3tmp159, FDPart3tmp258, FusedMulAddSIMD(FDPart3tmp274, FDPart3tmp374, MulSIMD(FDPart3tmp158, FDPart3tmp371)));
        const REAL_SIMD_ARRAY FDPart3tmp378 =
            FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp371, FusedMulAddSIMD(FDPart3tmp258, FDPart3tmp278, MulSIMD(FDPart3tmp145, FDPart3tmp374)));
        const REAL_SIMD_ARRAY FDPart3tmp379 =
            FusedMulAddSIMD(FDPart3tmp171, FDPart3tmp258, FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp374, MulSIMD(FDPart3tmp169, FDPart3tmp371)));
        const REAL_SIMD_ARRAY FDPart3tmp403 = FusedMulAddSIMD(FDPart3tmp121, FDPart3tmp402, MulSIMD(FDPart3tmp123, FDPart3tmp257));
        const REAL_SIMD_ARRAY FDPart3tmp404 = FusedMulAddSIMD(FDPart3tmp129, FDPart3tmp371, MulSIMD(FDPart3tmp130, FDPart3tmp257));
        const REAL_SIMD_ARRAY FDPart3tmp405 = FusedMulAddSIMD(FDPart3tmp103, FDPart3tmp402, MulSIMD(FDPart3tmp110, FDPart3tmp371));
        const REAL_SIMD_ARRAY FDPart3tmp406 =
            FusedMulAddSIMD(FDPart3tmp159, FDPart3tmp257, FusedMulAddSIMD(FDPart3tmp274, FDPart3tmp371, MulSIMD(FDPart3tmp158, FDPart3tmp402)));
        const REAL_SIMD_ARRAY FDPart3tmp407 =
            FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp402, FusedMulAddSIMD(FDPart3tmp257, FDPart3tmp278, MulSIMD(FDPart3tmp145, FDPart3tmp371)));
        const REAL_SIMD_ARRAY FDPart3tmp408 =
            FusedMulAddSIMD(FDPart3tmp171, FDPart3tmp257, FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp371, MulSIMD(FDPart3tmp169, FDPart3tmp402)));
        const REAL_SIMD_ARRAY FDPart3tmp416 = FusedMulAddSIMD(
            FDPart3tmp6, MulSIMD(f0_of_xx0__DD00, vetU0),
            FusedMulAddSIMD(
                alpha,
                FusedMulAddSIMD(FDPart3tmp299, FDPart3tmp91,
                                FusedMulAddSIMD(FDPart3tmp301, FDPart3tmp78,
                                                FusedMulAddSIMD(FDPart3tmp304, FDPart3tmp83,
                                                                FusedMulAddSIMD(FDPart3tmp7, FDPart3tmp94,
                                                                                FusedMulAddSIMD(FDPart3tmp184, FDPart3tmp97,
                                                                                                MulSIMD(FDPart3tmp197, FDPart3tmp86)))))),
                SubSIMD(FusedMulAddSIMD(FDPart3tmp51,
                                        MulSIMD(MulSIMD(FDPart3tmp42, FDPart3tmp44),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, vetU0))),
                                        MulSIMD(FDPart3_NegativeOne_, AddSIMD(FDPart3tmp38, FDPart3tmp22))),
                        MulSIMD(FDPart3tmp39, FDPart3tmp9))));
        const REAL_SIMD_ARRAY FDPart3tmp432 = FusedMulAddSIMD(FDPart3tmp147, alpha_dD0, FusedMulAddSIMD(FDPart3tmp149, alpha_dD1, alpha_dDD01));
        const REAL_SIMD_ARRAY FDPart3tmp180 = FusedMulAddSIMD(
            FDPart3tmp179, FDPart3tmp94,
            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp116, FDPart3tmp78),
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp128, FDPart3tmp83),
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp152, FDPart3tmp91),
                                                            FusedMulAddSIMD(FDPart3tmp135, FDPart3tmp86, MulSIMD(FDPart3tmp164, FDPart3tmp97))))));
        const REAL_SIMD_ARRAY FDPart3tmp293 = AddSIMD(FDPart3tmp235, FDPart3tmp268);
        const REAL_SIMD_ARRAY FDPart3tmp309 = MulSIMD(FDPart3tmp307, FDPart3tmp308);
        const REAL_SIMD_ARRAY FDPart3tmp311 = AddSIMD(FDPart3tmp283, FDPart3tmp310);
        const REAL_SIMD_ARRAY FDPart3tmp315 = MulSIMD(FDPart3tmp308, FDPart3tmp314);
        const REAL_SIMD_ARRAY FDPart3tmp316 = AddSIMD(FDPart3tmp226, FDPart3tmp275);
        const REAL_SIMD_ARRAY FDPart3tmp322 = MulSIMD(FDPart3tmp308, FDPart3tmp321);
        const REAL_SIMD_ARRAY FDPart3tmp338 = MulSIMD(FDPart3tmp285, FDPart3tmp337);
        const REAL_SIMD_ARRAY FDPart3tmp340 = MulSIMD(FDPart3tmp285, FDPart3tmp339);
        const REAL_SIMD_ARRAY FDPart3tmp342 = MulSIMD(FDPart3tmp285, FDPart3tmp341);
        const REAL_SIMD_ARRAY FDPart3tmp381 = AddSIMD(FDPart3tmp118, FDPart3tmp376);
        const REAL_SIMD_ARRAY FDPart3tmp384 = AddSIMD(FDPart3tmp378, FDPart3tmp383);
        const REAL_SIMD_ARRAY FDPart3tmp411 = AddSIMD(FDPart3tmp383, FDPart3tmp405);
        const REAL_SIMD_ARRAY FDPart3tmp413 = SubSIMD(FDPart3tmp403, FDPart3tmp39);
        const REAL_SIMD_ARRAY FDPart3tmp424 = MulSIMD(FDPart3tmp307, FDPart3tmp423);
        const REAL_SIMD_ARRAY FDPart3tmp425 = MulSIMD(FDPart3tmp314, FDPart3tmp423);
        const REAL_SIMD_ARRAY FDPart3tmp426 = MulSIMD(FDPart3tmp321, FDPart3tmp423);
        const REAL_SIMD_ARRAY FDPart3tmp428 = MulSIMD(FDPart3tmp337, FDPart3tmp427);
        const REAL_SIMD_ARRAY FDPart3tmp429 = MulSIMD(FDPart3tmp339, FDPart3tmp427);
        const REAL_SIMD_ARRAY FDPart3tmp430 = MulSIMD(FDPart3tmp341, FDPart3tmp427);
        const REAL_SIMD_ARRAY FDPart3tmp349 = MulSIMD(
            FDPart3tmp55, FusedMulAddSIMD(FDPart3tmp279, FDPart3tmp346,
                                          FusedMulAddSIMD(FDPart3tmp293, FDPart3tmp345,
                                                          FusedMulAddSIMD(FDPart3tmp311, FDPart3tmp348,
                                                                          FusedMulAddSIMD(FDPart3tmp316, FDPart3tmp347,
                                                                                          FusedMulAddSIMD(FDPart3tmp259, FDPart3tmp343,
                                                                                                          MulSIMD(FDPart3tmp264, FDPart3tmp344)))))));
        const REAL_SIMD_ARRAY FDPart3tmp385 = MulSIMD(
            FDPart3tmp55, FusedMulAddSIMD(FDPart3tmp345, FDPart3tmp381,
                                          FusedMulAddSIMD(FDPart3tmp346, FDPart3tmp384,
                                                          FusedMulAddSIMD(FDPart3tmp347, FDPart3tmp377,
                                                                          FusedMulAddSIMD(FDPart3tmp348, FDPart3tmp379,
                                                                                          FusedMulAddSIMD(FDPart3tmp343, FDPart3tmp375,
                                                                                                          MulSIMD(FDPart3tmp344, FDPart3tmp372)))))));
        const REAL_SIMD_ARRAY FDPart3tmp414 = MulSIMD(
            FDPart3tmp55, FusedMulAddSIMD(FDPart3tmp345, FDPart3tmp404,
                                          FusedMulAddSIMD(FDPart3tmp346, FDPart3tmp407,
                                                          FusedMulAddSIMD(FDPart3tmp347, FDPart3tmp406,
                                                                          FusedMulAddSIMD(FDPart3tmp348, FDPart3tmp408,
                                                                                          FusedMulAddSIMD(FDPart3tmp343, FDPart3tmp411,
                                                                                                          MulSIMD(FDPart3tmp344, FDPart3tmp413)))))));
        const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(
            FDPart3tmp4, aDD_dKOD001,
            FusedMulAddSIMD(
                FDPart3tmp6,
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp5, aDD00),
                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                            MulSIMD(alpha, FusedMulAddSIMD(FDPart3tmp78, FDPart3tmp87,
                                                           FusedMulAddSIMD(FDPart3tmp84, FDPart3tmp91, MulSIMD(FDPart3tmp7, FDPart3tmp94))))),
                    FusedMulAddSIMD(
                        f0_of_xx0,
                        MulSIMD(MulSIMD(aDD01, alpha),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                        MulSIMD(f0_of_xx0__D0,
                                                FusedMulAddSIMD(FDPart3tmp83, FDPart3tmp87,
                                                                FusedMulAddSIMD(FDPart3tmp84, FDPart3tmp97, MulSIMD(FDPart3tmp7, FDPart3tmp91)))))),
                        FusedMulAddSIMD(
                            FDPart3tmp34, MulSIMD(FDPart3tmp35, aDD02),
                            FusedMulAddSIMD(
                                FDPart3tmp5, MulSIMD(FDPart3tmp9, aDD_dupD001),
                                FusedMulAddSIMD(
                                    FDPart3_Integer_2, MulSIMD(FDPart3tmp26, FDPart3tmp7),
                                    FusedMulAddSIMD(
                                        FDPart3tmp19, MulSIMD(FDPart3tmp20, aDD01),
                                        FusedMulAddSIMD(
                                            FDPart3tmp7, FDPart3tmp8,
                                            FusedMulAddSIMD(
                                                FDPart3tmp98,
                                                NegFusedMulAddSIMD(
                                                    FDPart3tmp180,
                                                    FusedMulAddSIMD(FDPart3_Rational_1_3, FDPart3tmp5, MulSIMD(FDPart3_Rational_1_3, FDPart3tmp58)),
                                                    FDPart3tmp179),
                                                FusedMulAddSIMD(
                                                    MulSIMD(f0_of_xx0, f0_of_xx0__D0),
                                                    MulSIMD(MulSIMD(aDD02, alpha),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                    MulSIMD(f1_of_xx1,
                                                                            FusedMulAddSIMD(FDPart3tmp83, FDPart3tmp84,
                                                                                            FusedMulAddSIMD(FDPart3tmp86, FDPart3tmp87,
                                                                                                            MulSIMD(FDPart3tmp7, FDPart3tmp78)))))),
                                                    FusedMulSubSIMD(
                                                        FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp12, aDD00, MulSIMD(FDPart3tmp5, aDD_dupD000)),
                                                        MulSIMD(MulSIMD(FDPart3_Rational_2_3, FDPart3tmp5), MulSIMD(FDPart3tmp55, aDD00)))))))))))),
                MulSIMD(FDPart3tmp2, aDD_dKOD000)));
        const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(
            FDPart3tmp4, aDD_dKOD011,
            FusedMulAddSIMD(
                FDPart3tmp0,
                MulSIMD(FDPart3tmp3,
                        FusedMulAddSIMD(
                            alpha,
                            MulSIMD(MulSIMD(FDPart3tmp194, aDD01),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                            FusedMulAddSIMD(
                                aDD00, MulSIMD(f0_of_xx0__D0, vetU_dD01),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3tmp55, aDD01),
                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0, f0_of_xx0__D0)),
                                    FusedMulAddSIMD(
                                        FDPart3tmp98,
                                        NegFusedMulAddSIMD(f0_of_xx0,
                                                           MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp180), MulSIMD(f0_of_xx0__D0, hDD01)),
                                                           FDPart3tmp152),
                                        FusedMulAddSIMD(
                                            FDPart3tmp26, MulSIMD(FDPart3tmp74, aDD01),
                                            FusedMulAddSIMD(
                                                FDPart3tmp189, FDPart3tmp34,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp8, FDPart3tmp84,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp184, FDPart3tmp19,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp188, FDPart3tmp87,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp14, aDD_dupD011,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp182, aDD01,
                                                                    FusedMulAddSIMD(
                                                                        MulSIMD(alpha, f0_of_xx0),
                                                                        MulSIMD(MulSIMD(FDPart3tmp190, aDD02),
                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                        MulSIMD(f0_of_xx0__D0, f1_of_xx1))),
                                                                        FusedMulSubSIMD(
                                                                            FDPart3tmp13,
                                                                            FusedMulAddSIMD(FDPart3tmp74, aDD_dupD010,
                                                                                            MulSIMD(aDD01, AddSIMD(FDPart3tmp166, FDPart3tmp5))),
                                                                            MulSIMD(FDPart3tmp5, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp191),
                                                                                                         MulSIMD(aDD00, alpha))))))))))))))))),
                MulSIMD(FDPart3tmp2, aDD_dKOD010)));
        const REAL_SIMD_ARRAY __RHS_exp_2 = FusedMulAddSIMD(
            FDPart3tmp4, aDD_dKOD021,
            FusedMulAddSIMD(
                FDPart3tmp0,
                MulSIMD(FDPart3tmp29,
                        FusedMulAddSIMD(
                            alpha,
                            MulSIMD(MulSIMD(FDPart3tmp201, aDD01),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                            FusedMulAddSIMD(
                                f0_of_xx0,
                                MulSIMD(MulSIMD(FDPart3tmp55, aDD02),
                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0__D0, f1_of_xx1))),
                                FusedMulAddSIMD(
                                    FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp75, aDD_dupD021, MulSIMD(FDPart3tmp74, MulSIMD(aDD02, f1_of_xx1__D1))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp98,
                                        NegFusedMulAddSIMD(MulSIMD(f0_of_xx0, f0_of_xx0__D0),
                                                           MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp180), MulSIMD(f1_of_xx1, hDD02)),
                                                           FDPart3tmp116),
                                        FusedMulAddSIMD(
                                            FDPart3tmp26, FDPart3tmp87,
                                            FusedMulAddSIMD(
                                                FDPart3tmp8, FDPart3tmp87,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp189, FDPart3tmp19,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp197, FDPart3tmp34,
                                                        FusedMulAddSIMD(
                                                            MulSIMD(alpha, f0_of_xx0),
                                                            MulSIMD(MulSIMD(FDPart3tmp199, aDD02),
                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                            MulSIMD(f0_of_xx0__D0, f1_of_xx1))),
                                                            FusedMulSubSIMD(FDPart3tmp13,
                                                                            FusedMulAddSIMD(FDPart3tmp167, aDD02, MulSIMD(FDPart3tmp75, aDD_dupD020)),
                                                                            MulSIMD(FDPart3tmp5, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp200),
                                                                                                         MulSIMD(aDD00, alpha)))))))))))))),
                MulSIMD(FDPart3tmp2, aDD_dKOD020)));
        const REAL_SIMD_ARRAY __RHS_exp_3 = FusedMulAddSIMD(
            FDPart3tmp2, aDD_dKOD110,
            FusedMulAddSIMD(
                FDPart3tmp4, aDD_dKOD111,
                MulSIMD(
                    FDPart3tmp16,
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp15, FDPart3tmp194), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD11, alpha)),
                        FusedMulAddSIMD(
                            aDD12,
                            MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp190),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha, f1_of_xx1))),
                            FusedMulAddSIMD(
                                aDD_dupD111, MulSIMD(f0_of_xx0, vetU1),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3_Integer_2, aDD11), MulSIMD(f0_of_xx0, vetU_dD11),
                                    FusedMulAddSIMD(
                                        FDPart3tmp98,
                                        NegFusedMulAddSIMD(
                                            FDPart3tmp180,
                                            FusedMulAddSIMD(FDPart3_Rational_1_3, FDPart3tmp15, MulSIMD(FDPart3_Rational_1_3, FDPart3tmp62)),
                                            FDPart3tmp164),
                                        FusedMulAddSIMD(
                                            FDPart3tmp155, MulSIMD(FDPart3tmp188, aDD12),
                                            FusedMulAddSIMD(
                                                FDPart3tmp184, FDPart3tmp8,
                                                FusedMulAddSIMD(FDPart3tmp205, aDD01,
                                                                FusedMulAddSIMD(alpha,
                                                                                MulSIMD(MulSIMD(FDPart3tmp191, aDD01),
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                                FusedMulSubSIMD(FDPart3tmp13,
                                                                                                FusedMulAddSIMD(FDPart3tmp15, aDD_dupD110,
                                                                                                                MulSIMD(FDPart3tmp20, aDD11)),
                                                                                                MulSIMD(MulSIMD(FDPart3_Rational_2_3, FDPart3tmp15),
                                                                                                        MulSIMD(FDPart3tmp55, aDD11)))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_4 = FusedMulAddSIMD(
            FDPart3tmp30,
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp15, FDPart3tmp55), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(aDD12, f1_of_xx1)),
                FusedMulAddSIMD(
                    aDD12,
                    MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp199),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha, f1_of_xx1))),
                    FusedMulAddSIMD(
                        FDPart3tmp207, MulSIMD(aDD02, vetU_dD01),
                        FusedMulAddSIMD(
                            FDPart3tmp207, MulSIMD(aDD12, vetU_dD11),
                            FusedMulAddSIMD(
                                FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp80, aDD_dupD121, MulSIMD(FDPart3tmp15, MulSIMD(aDD12, f1_of_xx1__D1))),
                                FusedMulAddSIMD(
                                    FDPart3tmp98,
                                    NegFusedMulAddSIMD(FDPart3tmp180, MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp15), MulSIMD(f1_of_xx1, hDD12)),
                                                       FDPart3tmp128),
                                    FusedMulAddSIMD(
                                        FDPart3tmp188, FDPart3tmp197,
                                        FusedMulAddSIMD(
                                            FDPart3tmp189, FDPart3tmp8,
                                            FusedMulAddSIMD(
                                                alpha,
                                                MulSIMD(MulSIMD(FDPart3tmp200, aDD01),
                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                FusedMulSubSIMD(FDPart3tmp13,
                                                                FusedMulAddSIMD(FDPart3tmp35, aDD12, MulSIMD(FDPart3tmp80, aDD_dupD120)),
                                                                MulSIMD(FDPart3tmp201, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp15),
                                                                                               MulSIMD(aDD11, alpha))))))))))))),
            FusedMulAddSIMD(FDPart3tmp4, aDD_dKOD121, MulSIMD(FDPart3tmp2, aDD_dKOD120)));
        const REAL_SIMD_ARRAY __RHS_exp_5 = FusedMulAddSIMD(
            FDPart3tmp210,
            FusedMulAddSIMD(
                aDD12,
                MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp201), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha, f1_of_xx1))),
                FusedMulAddSIMD(
                    FDPart3tmp98,
                    NegFusedMulAddSIMD(FDPart3tmp180,
                                       FusedMulAddSIMD(FDPart3_Rational_1_3, FDPart3tmp64, MulSIMD(FDPart3_Rational_1_3, FDPart3tmp67)),
                                       FDPart3tmp135),
                    FusedMulAddSIMD(
                        FDPart3tmp41,
                        MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp199),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD22, alpha))),
                        FusedMulAddSIMD(
                            FDPart3tmp197, FDPart3tmp8,
                            FusedMulAddSIMD(
                                FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp119, aDD22, MulSIMD(FDPart3tmp64, aDD_dupD221)),
                                FusedMulAddSIMD(MulSIMD(alpha, f0_of_xx0),
                                                MulSIMD(MulSIMD(FDPart3tmp200, aDD02),
                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f1_of_xx1))),
                                                FusedMulSubSIMD(FDPart3tmp13,
                                                                FusedMulAddSIMD(FDPart3tmp101, aDD22, MulSIMD(FDPart3tmp64, aDD_dupD220)),
                                                                MulSIMD(FDPart3tmp41, MulSIMD(MulSIMD(FDPart3_Rational_2_3, FDPart3tmp15),
                                                                                              MulSIMD(FDPart3tmp55, aDD22)))))))))),
            FusedMulAddSIMD(FDPart3tmp4, aDD_dKOD221, MulSIMD(FDPart3tmp2, aDD_dKOD220)));
        const REAL_SIMD_ARRAY __RHS_exp_6 = FusedMulAddSIMD(
            FDPart3tmp212, alpha_dKOD0,
            FusedMulAddSIMD(FDPart3tmp213, alpha_dKOD1,
                            FusedMulAddSIMD(FDPart3tmp9, alpha_dupD1,
                                            FusedMulSubSIMD(FDPart3tmp13, alpha_dupD0, MulSIMD(FDPart3_Integer_2, MulSIMD(alpha, trK))))));
        const REAL_SIMD_ARRAY __RHS_exp_7 = FusedMulAddSIMD(
            FDPart3tmp213, betU_dKOD01,
            FusedMulAddSIMD(
                f0_of_xx0__D0,
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp28, FDPart3tmp3),
                    MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp259),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU0))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp41, FDPart3tmp72),
                        MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp188),
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_2), MulSIMD(FDPart3tmp82, f0_of_xx0))),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3tmp28, FDPart3tmp3),
                            MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp259),
                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU2))),
                            FusedMulAddSIMD(
                                FDPart3tmp28,
                                MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp264),
                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU1))),
                                FusedMulAddSIMD(
                                    FDPart3tmp42,
                                    MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp268),
                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU2))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp3,
                                        MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp279),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU0))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp28,
                                            MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp264),
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU2))),
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3tmp283, FDPart3tmp6),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU0)),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp3,
                                                    MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp279),
                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU1))),
                                                    FusedMulAddSIMD(
                                                        MulSIMD(FDPart3tmp0, FDPart3tmp3),
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU_dD01)),
                                                        FusedMulAddSIMD(
                                                            MulSIMD(FDPart3tmp16, FDPart3tmp275),
                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU1)),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp264, MulSIMD(FDPart3tmp266, vetU1),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp0,
                                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4),
                                                                            MulSIMD(FDPart3tmp26, lambdaU0)),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3_Rational_3_4,
                                                                        MulSIMD(
                                                                            FDPart3tmp94,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp13,
                                                                                FusedMulSubSIMD(
                                                                                    FDPart3tmp0, f0_of_xx0__DDD000,
                                                                                    MulSIMD(FDPart3tmp6, MulSIMD(f0_of_xx0__DD00, f0_of_xx0__DD00))),
                                                                                FusedMulAddSIMD(FDPart3tmp0, MulSIMD(FDPart3tmp26, f0_of_xx0__DD00),
                                                                                                FDPart3tmp245))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3_Rational_3_4,
                                                                            MulSIMD(FDPart3tmp97,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp0, vetU_dDD011,
                                                                                        FusedMulAddSIMD(FDPart3tmp226, FDPart3tmp26,
                                                                                                        NegFusedMulAddSIMD(
                                                                                                            FDPart3_Integer_2,
                                                                                                            MulSIMD(FDPart3tmp0, vetU_dD11),
                                                                                                            SubSIMD(FDPart3tmp230, FDPart3tmp13))))),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3_Rational_3_4,
                                                                                MulSIMD(FDPart3tmp86,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp230, FDPart3tmp41,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp235, FDPart3tmp26,
                                                                                                FusedMulAddSIMD(
                                                                                                    MulSIMD(FDPart3tmp0, f1_of_xx1),
                                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                            MulSIMD(f1_of_xx1__D1, vetU1)),
                                                                                                    FusedMulSubSIMD(
                                                                                                        FDPart3tmp118, FDPart3tmp234,
                                                                                                        MulSIMD(FDPart3tmp13, FDPart3tmp41)))))),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3_Rational_3_4,
                                                                                    MulSIMD(FDPart3tmp91,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp218, vetU1,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp9,
                                                                                                    FusedMulAddSIMD(FDPart3tmp166, FDPart3tmp6,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                    AddSIMD(FDPart3tmp228, FDPart3tmp9)))),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3_Rational_3_4, MulSIMD(FDPart3tmp323, FDPart3tmp324),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3_Rational_3_4,
                                                                                            MulSIMD(FDPart3tmp78,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp238,
                                                                                                        FusedMulSubSIMD(
                                                                                                            FDPart3tmp166,
                                                                                                            MulSIMD(FDPart3tmp41, FDPart3tmp6),
                                                                                                            FDPart3tmp41),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp218, MulSIMD(f1_of_xx1, vetU2),
                                                                                                            FusedMulAddSIMD(FDPart3tmp185, f1_of_xx1,
                                                                                                                            FDPart3tmp240)))),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3_Rational_3_4,
                                                                                                MulSIMD(FDPart3tmp228, FDPart3tmp91),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3_Rational_3_4,
                                                                                                    MulSIMD(
                                                                                                        FDPart3tmp286,
                                                                                                        FusedMulAddSIMD(
                                                                                                            MulSIMD(FDPart3tmp93, aDD01),
                                                                                                            MulSIMD(
                                                                                                                MulSIMD(FDPart3tmp288, FDPart3tmp90),
                                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_4,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(f0_of_xx0,
                                                                                                                                f0_of_xx0__D0))),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp288,
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3tmp15,
                                                                                                                            FDPart3tmp287),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(FDPart3tmp41,
                                                                                                                                aDD22))),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    MulSIMD(FDPart3tmp77,
                                                                                                                            FDPart3tmp90),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3tmp15,
                                                                                                                                FDPart3tmp288),
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3_Integer_4,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                            MulSIMD(aDD12,
                                                                                                                                    f1_of_xx1))),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        f0_of_xx0,
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(FDPart3tmp93,
                                                                                                                                    aDD02),
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(FDPart3tmp288,
                                                                                                                                        FDPart3tmp77),
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Integer_4,
                                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                                    MulSIMD(
                                                                                                                                        f0_of_xx0__D0,
                                                                                                                                        f1_of_xx1)))),
                                                                                                                        FusedMulSubSIMD(
                                                                                                                            MulSIMD(FDPart3tmp288,
                                                                                                                                    FDPart3tmp5),
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                    FDPart3_NegativeOne_),
                                                                                                                                MulSIMD(
                                                                                                                                    aDD00,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp93,
                                                                                                                                        FDPart3tmp93))),
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3tmp288,
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                        FDPart3tmp15),
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp289,
                                                                                                                                        aDD11))))))))),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp311, FDPart3tmp315,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp316, FDPart3tmp322,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp283, FDPart3tmp284,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp293, FDPart3tmp309,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp279, FDPart3tmp282,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp279,
                                                                                                                            FDPart3tmp342,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp275,
                                                                                                                                FDPart3tmp277,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp279,
                                                                                                                                    FDPart3tmp280,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp268,
                                                                                                                                        FDPart3tmp270,
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3tmp273,
                                                                                                                                            FDPart3tmp94,
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3tmp264,
                                                                                                                                                FDPart3tmp265,
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3tmp264,
                                                                                                                                                    FDPart3tmp340,
                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                        FDPart3tmp259,
                                                                                                                                                        FDPart3tmp262, FusedMulAddSIMD(FDPart3tmp259, FDPart3tmp338, FusedMulAddSIMD(FDPart3tmp256, FDPart3tmp91, FusedMulAddSIMD(FDPart3tmp259, FDPart3tmp260, FusedMulAddSIMD(FDPart3tmp220, FDPart3tmp91, FusedMulAddSIMD(FDPart3tmp221, FDPart3tmp94, FusedMulAddSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp0, betU_dupD00, MulSIMD(FDPart3tmp218, betU0)), FusedMulAddSIMD(FDPart3tmp216, betU_dupD01, FusedMulAddSIMD(MulSIMD(FDPart3tmp41, FDPart3tmp72), MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp34), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(FDPart3tmp77, f0_of_xx0))), FusedMulSubSIMD(FDPart3_Rational_3_4, FDPart3tmp349, MulSIMD(FDPart3tmp214, eta)))))))))))))))))))))))))))))))))))))))))))),
                MulSIMD(FDPart3tmp212, betU_dKOD00)));
        const REAL_SIMD_ARRAY __RHS_exp_8 = FusedMulAddSIMD(
            FDPart3tmp213, betU_dKOD11,
            FusedMulAddSIMD(
                f0_of_xx0,
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp3, FDPart3tmp375),
                    MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp28),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU2))),
                    FusedMulAddSIMD(
                        FDPart3tmp42,
                        MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp376),
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU2))),
                        FusedMulAddSIMD(
                            FDPart3tmp82,
                            MulSIMD(MulSIMD(FDPart3tmp188, FDPart3tmp72),
                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(f1_of_xx1, f1_of_xx1__D1))),
                            FusedMulAddSIMD(
                                FDPart3tmp372,
                                MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp28),
                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU2))),
                                FusedMulAddSIMD(
                                    FDPart3tmp372,
                                    MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp28),
                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU1))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp378,
                                        MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp3),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU1))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp378,
                                            MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp3),
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU0))),
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3tmp379, FDPart3tmp6),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU0)),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp77,
                                                    MulSIMD(MulSIMD(FDPart3tmp34, FDPart3tmp72),
                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_2),
                                                                    MulSIMD(f1_of_xx1, f1_of_xx1__D1))),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp16,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU_dD11)),
                                                        FusedMulAddSIMD(
                                                            MulSIMD(FDPart3tmp16, FDPart3tmp377),
                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU1)),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp266, MulSIMD(FDPart3tmp372, vetU1),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp0,
                                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4),
                                                                            MulSIMD(FDPart3tmp19, lambdaU0)),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3_Rational_3_4,
                                                                        MulSIMD(FDPart3tmp86,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp19, FDPart3tmp235,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp41, FDPart3tmp9,
                                                                                        FusedMulSubSIMD(FDPart3tmp118, FDPart3tmp38,
                                                                                                        MulSIMD(FDPart3tmp248, FDPart3tmp9))))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3_Rational_3_4,
                                                                            MulSIMD(
                                                                                FDPart3tmp94,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp11, MulSIMD(FDPart3tmp19, FDPart3tmp3),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp367, vetU1,
                                                                                        FusedMulAddSIMD(
                                                                                            vetU1,
                                                                                            FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp367,
                                                                                                            MulSIMD(FDPart3tmp16, f0_of_xx0__DD00)),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp3, vetU_dDD100,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp359, FDPart3tmp9,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp16,
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(f0_of_xx0__D0, vetU_dD10)),
                                                                                                        FusedMulSubSIMD(
                                                                                                            FDPart3tmp19, FDPart3tmp310,
                                                                                                            MulSIMD(FDPart3tmp276,
                                                                                                                    f0_of_xx0__DD00))))))))),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3_Rational_3_4,
                                                                                MulSIMD(
                                                                                    FDPart3tmp323,
                                                                                    FusedMulAddSIMD(
                                                                                        MulSIMD(FDPart3tmp96, aDD01),
                                                                                        MulSIMD(
                                                                                            MulSIMD(FDPart3tmp288, FDPart3tmp90),
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_),
                                                                                                    MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp295,
                                                                                            MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                            MulSIMD(FDPart3tmp41, aDD22))),
                                                                                            FusedMulAddSIMD(
                                                                                                MulSIMD(FDPart3tmp82, FDPart3tmp96),
                                                                                                MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_4,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(aDD12, f1_of_xx1))),
                                                                                                FusedMulAddSIMD(
                                                                                                    f0_of_xx0,
                                                                                                    MulSIMD(
                                                                                                        MulSIMD(FDPart3tmp90, aDD02),
                                                                                                        MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp82),
                                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_4,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(f0_of_xx0__D0,
                                                                                                                                f1_of_xx1)))),
                                                                                                    FusedMulSubSIMD(
                                                                                                        MulSIMD(FDPart3tmp288, FDPart3tmp289),
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(FDPart3tmp5, aDD00)),
                                                                                                        MulSIMD(
                                                                                                            FDPart3tmp288,
                                                                                                            MulSIMD(
                                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3tmp15),
                                                                                                                MulSIMD(
                                                                                                                    aDD11,
                                                                                                                    MulSIMD(FDPart3tmp96,
                                                                                                                            FDPart3tmp96)))))))))),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3_Rational_3_4,
                                                                                    MulSIMD(FDPart3tmp83,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp238, FDPart3tmp248,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp238,
                                                                                                    NegFusedMulSubSIMD(f1_of_xx1, f1_of_xx1__DD11,
                                                                                                                       FDPart3tmp248),
                                                                                                    NegFusedMulAddSIMD(FDPart3tmp185, f1_of_xx1,
                                                                                                                       FDPart3tmp364)))),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp342, FDPart3tmp384,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3_Rational_3_4,
                                                                                            MulSIMD(FDPart3tmp286, FDPart3tmp324),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp338, FDPart3tmp375,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp340, FDPart3tmp372,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp315, FDPart3tmp379,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp322, FDPart3tmp377,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp284, FDPart3tmp379,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp309, FDPart3tmp381,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp280, FDPart3tmp378,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp282,
                                                                                                                            FDPart3tmp378,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp276,
                                                                                                                                betU_dupD11,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp277,
                                                                                                                                    FDPart3tmp377,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp270,
                                                                                                                                        FDPart3tmp376,
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3tmp273,
                                                                                                                                            FDPart3tmp91,
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3tmp262,
                                                                                                                                                FDPart3tmp375,
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3tmp265,
                                                                                                                                                    FDPart3tmp372,
                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                        FDPart3tmp256,
                                                                                                                                                        FDPart3tmp97,
                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                            FDPart3tmp260,
                                                                                                                                                            FDPart3tmp375,
                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                FDPart3tmp220,
                                                                                                                                                                FDPart3tmp97,
                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                    FDPart3tmp221,
                                                                                                                                                                    FDPart3tmp91,
                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                        FDPart3_Rational_3_4,
                                                                                                                                                                        FDPart3tmp385,
                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                            FDPart3tmp13,
                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                FDPart3tmp3,
                                                                                                                                                                                betU_dupD10,
                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                    FDPart3tmp350,
                                                                                                                                                                                    betU1)),
                                                                                                                                                                            FusedMulAddSIMD(FDPart3_Rational_3_4, FDPart3tmp358, FusedMulAddSIMD(FDPart3_Rational_3_4, FDPart3tmp361, FusedMulAddSIMD(MulSIMD(FDPart3tmp3, FDPart3tmp375), MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp28), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU0))), FusedMulSubSIMD(FDPart3_Rational_3_4, FDPart3tmp351, MulSIMD(FDPart3tmp281, eta)))))))))))))))))))))))))))))))))))))))))))),
                MulSIMD(FDPart3tmp212, betU_dKOD10)));
        const REAL_SIMD_ARRAY __RHS_exp_9 = FusedMulAddSIMD(
            FDPart3tmp212, betU_dKOD20,
            FusedMulAddSIMD(
                FDPart3tmp213, betU_dKOD21,
                MulSIMD(
                    FDPart3tmp207,
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp3, FDPart3tmp405),
                        MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp28),
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU2))),
                        FusedMulAddSIMD(
                            FDPart3tmp403,
                            MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp28),
                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU1))),
                            FusedMulAddSIMD(
                                FDPart3tmp42,
                                MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp404),
                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU2))),
                                FusedMulAddSIMD(
                                    FDPart3tmp407,
                                    MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp3),
                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU0))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp403,
                                        MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp28),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU2))),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3tmp408, FDPart3tmp6),
                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU0)),
                                            FusedMulAddSIMD(
                                                FDPart3tmp407,
                                                MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp3),
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU1))),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp188,
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(FDPart3tmp3, lambdaU1)),
                                                    FusedMulAddSIMD(
                                                        MulSIMD(FDPart3tmp16, FDPart3tmp406),
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU1)),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp266, MulSIMD(FDPart3tmp403, vetU1),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp0,
                                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4),
                                                                        MulSIMD(FDPart3tmp34, lambdaU0)),
                                                                FusedMulAddSIMD(
                                                                    FDPart3_Rational_3_4,
                                                                    MulSIMD(FDPart3tmp94,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp11, MulSIMD(FDPart3tmp3, FDPart3tmp34),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp399, vetU2,
                                                                                    FusedMulAddSIMD(
                                                                                        vetU2,
                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp28, FDPart3tmp367),
                                                                                                        FDPart3tmp399),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp29, vetU_dDD200,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp310, FDPart3tmp34,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp28, MulSIMD(FDPart3tmp367, vetU2),
                                                                                                    FusedMulSubSIMD(
                                                                                                        FDPart3tmp238, FDPart3tmp359,
                                                                                                        MulSIMD(FDPart3tmp28,
                                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3tmp16),
                                                                                                                        MulSIMD(f0_of_xx0__D0,
                                                                                                                                vetU_dD20))))))))))),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3_Rational_3_4,
                                                                        MulSIMD(
                                                                            FDPart3tmp97,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3_Integer_2, MulSIMD(FDPart3tmp188, FDPart3tmp39),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp29, vetU_dDD211,
                                                                                    FusedMulAddSIMD(
                                                                                        vetU2,
                                                                                        FusedMulSubSIMD(
                                                                                            MulSIMD(FDPart3_Integer_2, FDPart3tmp248),
                                                                                            MulSIMD(FDPart3tmp249, FDPart3tmp3),
                                                                                            MulSIMD(FDPart3tmp3,
                                                                                                    MulSIMD(FDPart3tmp42, f1_of_xx1__DD11))),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp226, FDPart3tmp34,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp238, FDPart3tmp391,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp185,
                                                                                                    MulSIMD(FDPart3tmp248, FDPart3tmp249),
                                                                                                    NegFusedMulAddSIMD(
                                                                                                        FDPart3tmp42,
                                                                                                        MulSIMD(
                                                                                                            MulSIMD(FDPart3_Integer_2, FDPart3tmp3),
                                                                                                            MulSIMD(f1_of_xx1__D1, vetU_dD21)),
                                                                                                        FDPart3tmp238)))))))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3_Rational_3_4,
                                                                            MulSIMD(
                                                                                FDPart3tmp323,
                                                                                FusedMulAddSIMD(
                                                                                    MulSIMD(FDPart3tmp90, aDD01),
                                                                                    MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp82),
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                    MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                                    FusedMulAddSIMD(
                                                                                        f0_of_xx0,
                                                                                        MulSIMD(MulSIMD(FDPart3tmp82, aDD02),
                                                                                                MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp77),
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(f0_of_xx0__D0, f1_of_xx1)))),
                                                                                        FusedMulAddSIMD(
                                                                                            MulSIMD(FDPart3tmp85, FDPart3tmp96),
                                                                                            MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                            MulSIMD(aDD12, f1_of_xx1))),
                                                                                            FusedMulAddSIMD(
                                                                                                MulSIMD(FDPart3tmp96, aDD01),
                                                                                                MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp77),
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp77,
                                                                                                    MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp5),
                                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                            FDPart3_NegativeOne_),
                                                                                                                    MulSIMD(FDPart3tmp90, aDD00))),
                                                                                                    FusedMulAddSIMD(
                                                                                                        MulSIMD(FDPart3tmp41, FDPart3tmp82),
                                                                                                        MulSIMD(
                                                                                                            MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                            FDPart3_NegativeOne_),
                                                                                                                    MulSIMD(FDPart3tmp85, aDD22))),
                                                                                                        FusedMulAddSIMD(
                                                                                                            f0_of_xx0,
                                                                                                            MulSIMD(
                                                                                                                MulSIMD(FDPart3tmp90, aDD02),
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3tmp288,
                                                                                                                            FDPart3tmp85),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(f0_of_xx0__D0,
                                                                                                                                f1_of_xx1)))),
                                                                                                            FusedMulSubSIMD(
                                                                                                                FDPart3tmp82,
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3tmp15,
                                                                                                                            FDPart3tmp288),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(FDPart3tmp96,
                                                                                                                                aDD11))),
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3tmp288,
                                                                                                                            FDPart3tmp295),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3tmp15),
                                                                                                                        MulSIMD(aDD12,
                                                                                                                                f1_of_xx1)))))))))))),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3_Rational_3_4,
                                                                                MulSIMD(FDPart3tmp83,
                                                                                        FusedMulAddSIMD(FDPart3tmp388, FDPart3tmp9,
                                                                                                        SubSIMD(FDPart3tmp389, FDPart3tmp9))),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3_Rational_3_2, MulSIMD(FDPart3tmp398, FDPart3tmp91),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3_Rational_3_4,
                                                                                        MulSIMD(
                                                                                            FDPart3tmp286,
                                                                                            FusedMulAddSIMD(
                                                                                                MulSIMD(FDPart3tmp90, aDD01),
                                                                                                MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp77),
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                                                FusedMulAddSIMD(
                                                                                                    MulSIMD(FDPart3tmp93, aDD01),
                                                                                                    MulSIMD(
                                                                                                        MulSIMD(FDPart3tmp288, FDPart3tmp82),
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                                                    FusedMulAddSIMD(
                                                                                                        MulSIMD(FDPart3tmp85, FDPart3tmp90),
                                                                                                        MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(aDD12, f1_of_xx1))),
                                                                                                        FusedMulAddSIMD(
                                                                                                            MulSIMD(aDD02, f0_of_xx0),
                                                                                                            MulSIMD(
                                                                                                                MulSIMD(FDPart3tmp287, FDPart3tmp288),
                                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(f0_of_xx0__D0,
                                                                                                                                f1_of_xx1))),
                                                                                                            FusedMulAddSIMD(
                                                                                                                MulSIMD(FDPart3tmp41, FDPart3tmp77),
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3tmp15,
                                                                                                                            FDPart3tmp288),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(FDPart3tmp85,
                                                                                                                                aDD22))),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    MulSIMD(FDPart3tmp77,
                                                                                                                            FDPart3tmp82),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3tmp15,
                                                                                                                                FDPart3tmp288),
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                            MulSIMD(aDD12,
                                                                                                                                    f1_of_xx1))),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        f0_of_xx0,
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(FDPart3tmp93,
                                                                                                                                    aDD02),
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(FDPart3tmp288,
                                                                                                                                        FDPart3tmp85),
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                                    MulSIMD(
                                                                                                                                        f0_of_xx0__D0,
                                                                                                                                        f1_of_xx1)))),
                                                                                                                        FusedMulSubSIMD(
                                                                                                                            FDPart3tmp77,
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(FDPart3tmp288,
                                                                                                                                        FDPart3tmp5),
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp93,
                                                                                                                                        aDD00))),
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(FDPart3tmp288,
                                                                                                                                        FDPart3tmp82),
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                        FDPart3tmp15),
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp90,
                                                                                                                                        aDD11)))))))))))),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp342, FDPart3tmp407,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp9,
                                                                                                FusedMulAddSIMD(FDPart3tmp29, betU_dupD21,
                                                                                                                MulSIMD(FDPart3tmp387, betU2)),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp338, FDPart3tmp411,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp340, FDPart3tmp413,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp315, FDPart3tmp408,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp322, FDPart3tmp406,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp284, FDPart3tmp408,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp309, FDPart3tmp404,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp280,
                                                                                                                            FDPart3tmp407,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp282,
                                                                                                                                FDPart3tmp407,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp273,
                                                                                                                                    FDPart3tmp78,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp277,
                                                                                                                                        FDPart3tmp406,
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3tmp265,
                                                                                                                                            FDPart3tmp403,
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3tmp270,
                                                                                                                                                FDPart3tmp404,
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3tmp260,
                                                                                                                                                    FDPart3tmp405,
                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                        FDPart3tmp262,
                                                                                                                                                        FDPart3tmp405,
                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                            FDPart3tmp221,
                                                                                                                                                            FDPart3tmp78,
                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                FDPart3tmp256,
                                                                                                                                                                FDPart3tmp83,
                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                    FDPart3tmp13,
                                                                                                                                                                    FusedMulSubSIMD(
                                                                                                                                                                        FDPart3tmp29,
                                                                                                                                                                        betU_dupD20,
                                                                                                                                                                        MulSIMD(
                                                                                                                                                                            FDPart3tmp266,
                                                                                                                                                                            f0_of_xx0__D0)),
                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                        FDPart3tmp220,
                                                                                                                                                                        FDPart3tmp83,
                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                            FDPart3_Rational_3_4,
                                                                                                                                                                            FDPart3tmp396,
                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                FDPart3_Rational_3_4,
                                                                                                                                                                                FDPart3tmp414,
                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                    FDPart3_Rational_3_4,
                                                                                                                                                                                    FDPart3tmp394,
                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                        FDPart3_Rational_3_4, FDPart3tmp395, FusedMulAddSIMD(MulSIMD(FDPart3tmp3, FDPart3tmp405), MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp28), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU0))), FusedMulSubSIMD(FDPart3_Rational_3_4, FDPart3tmp392, MulSIMD(FDPart3tmp261, eta)))))))))))))))))))))))))))))))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_10 = FusedMulAddSIMD(
            FDPart3tmp4, cf_dKOD1,
            FusedMulSubSIMD(
                FDPart3tmp2, cf_dKOD0,
                MulSIMD(
                    FDPart3_Integer_2,
                    MulSIMD(cf,
                            FusedMulAddSIMD(
                                MulSIMD(FDPart3tmp0, FDPart3tmp115),
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(cf_dupD0, vetU0)),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3_Rational_1_12, FDPart3tmp42), MulSIMD(FDPart3tmp44, FDPart3tmp54),
                                    FusedMulAddSIMD(
                                        FDPart3tmp6, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_6), MulSIMD(f0_of_xx0__DD00, vetU0)),
                                        FusedMulAddSIMD(FDPart3_Rational_1_6, FDPart3tmp38,
                                                        FusedMulAddSIMD(FDPart3_Rational_1_6, MulSIMD(FDPart3tmp39, FDPart3tmp9),
                                                                        FusedMulAddSIMD(MulSIMD(FDPart3tmp115, FDPart3tmp3),
                                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                                MulSIMD(cf_dupD1, vetU1)),
                                                                                        FusedMulSubSIMD(FDPart3_Rational_1_6, FDPart3tmp22,
                                                                                                        MulSIMD(FDPart3_Rational_1_6,
                                                                                                                MulSIMD(alpha, trK)))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_11 = FusedMulAddSIMD(
            FDPart3tmp4, hDD_dKOD001,
            FusedMulAddSIMD(
                FDPart3tmp6,
                FusedMulAddSIMD(
                    FDPart3tmp19, MulSIMD(FDPart3tmp20, hDD01),
                    FusedMulAddSIMD(
                        FDPart3tmp34, MulSIMD(FDPart3tmp35, hDD02),
                        FusedMulAddSIMD(
                            FDPart3tmp416, FusedMulAddSIMD(FDPart3_Rational_2_3, FDPart3tmp5, MulSIMD(FDPart3_Rational_2_3, FDPart3tmp58)),
                            FusedMulAddSIMD(
                                FDPart3_Integer_2, MulSIMD(FDPart3tmp26, FDPart3tmp59),
                                FusedMulAddSIMD(FDPart3tmp5, MulSIMD(FDPart3tmp9, hDD_dupD001),
                                                FusedMulSubSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp5, hDD_dupD000, FDPart3tmp170),
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp5), MulSIMD(aDD00, alpha)))))))),
                MulSIMD(FDPart3tmp2, hDD_dKOD000)));
        const REAL_SIMD_ARRAY __RHS_exp_12 = FusedMulAddSIMD(
            FDPart3tmp4, hDD_dKOD011,
            FusedMulAddSIMD(
                FDPart3tmp0,
                MulSIMD(FDPart3tmp3,
                        FusedMulAddSIMD(
                            FDPart3tmp26, MulSIMD(FDPart3tmp74, hDD01),
                            FusedMulAddSIMD(
                                FDPart3tmp234, FDPart3tmp59,
                                FusedMulAddSIMD(
                                    FDPart3tmp34, FDPart3tmp81,
                                    FusedMulAddSIMD(
                                        FDPart3tmp188, FDPart3tmp76,
                                        FusedMulAddSIMD(
                                            FDPart3tmp19, FDPart3tmp63,
                                            FusedMulAddSIMD(
                                                FDPart3tmp14, hDD_dupD011,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp182, hDD01,
                                                    FusedMulAddSIMD(
                                                        MulSIMD(FDPart3_Rational_2_3, FDPart3tmp416), MulSIMD(FDPart3tmp74, hDD01),
                                                        FusedMulSubSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp74, hDD_dupD010, FDPart3tmp173),
                                                                        MulSIMD(alpha, MulSIMD(MulSIMD(FDPart3_Integer_2, aDD01),
                                                                                               MulSIMD(f0_of_xx0, f0_of_xx0__D0))))))))))))),
                MulSIMD(FDPart3tmp2, hDD_dKOD010)));
        const REAL_SIMD_ARRAY __RHS_exp_13 = FusedMulAddSIMD(
            FDPart3tmp4, hDD_dKOD021,
            FusedMulAddSIMD(
                FDPart3tmp0,
                MulSIMD(FDPart3tmp29,
                        FusedMulAddSIMD(
                            FDPart3tmp34, FDPart3tmp68,
                            FusedMulAddSIMD(
                                FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp75, hDD_dupD021, FDPart3tmp106),
                                FusedMulAddSIMD(
                                    FDPart3tmp19, FDPart3tmp81,
                                    FusedMulAddSIMD(
                                        FDPart3tmp26, FDPart3tmp76,
                                        FusedMulAddSIMD(
                                            FDPart3_Rational_2_3, MulSIMD(FDPart3tmp416, FDPart3tmp76),
                                            FusedMulSubSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp167, hDD02, MulSIMD(FDPart3tmp75, hDD_dupD020)),
                                                            MulSIMD(MulSIMD(alpha, f0_of_xx0), MulSIMD(MulSIMD(FDPart3_Integer_2, aDD02),
                                                                                                       MulSIMD(f0_of_xx0__D0, f1_of_xx1)))))))))),
                MulSIMD(FDPart3tmp2, hDD_dKOD020)));
        const REAL_SIMD_ARRAY __RHS_exp_14 = FusedMulAddSIMD(
            FDPart3tmp2, hDD_dKOD110,
            FusedMulAddSIMD(
                FDPart3tmp4, hDD_dKOD111,
                MulSIMD(
                    FDPart3tmp16,
                    FusedMulAddSIMD(
                        FDPart3_Integer_2, MulSIMD(FDPart3tmp38, FDPart3tmp63),
                        FusedMulAddSIMD(
                            FDPart3tmp155, MulSIMD(FDPart3tmp188, hDD12),
                            FusedMulAddSIMD(
                                FDPart3tmp205, hDD01,
                                FusedMulAddSIMD(
                                    FDPart3tmp416, FusedMulAddSIMD(FDPart3_Rational_2_3, FDPart3tmp15, MulSIMD(FDPart3_Rational_2_3, FDPart3tmp62)),
                                    FusedMulAddSIMD(f0_of_xx0, MulSIMD(hDD_dupD111, vetU1),
                                                    FusedMulSubSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp15, hDD_dupD110, FDPart3tmp144),
                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp15), MulSIMD(aDD11, alpha)))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_15 = FusedMulAddSIMD(
            FDPart3tmp30,
            FusedMulAddSIMD(
                FDPart3tmp207, MulSIMD(hDD02, vetU_dD01),
                FusedMulAddSIMD(
                    FDPart3tmp207, MulSIMD(hDD12, vetU_dD11),
                    FusedMulAddSIMD(
                        FDPart3tmp188, FDPart3tmp68,
                        FusedMulAddSIMD(
                            FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp80, hDD_dupD121, FDPart3tmp157),
                            FusedMulAddSIMD(
                                FDPart3tmp416, MulSIMD(MulSIMD(FDPart3_Rational_2_3, FDPart3tmp15), MulSIMD(f1_of_xx1, hDD12)),
                                FusedMulSubSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp35, hDD12, MulSIMD(FDPart3tmp80, hDD_dupD120)),
                                                MulSIMD(aDD12, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp15), MulSIMD(alpha, f1_of_xx1))))))))),
            FusedMulAddSIMD(FDPart3tmp4, hDD_dKOD121, MulSIMD(FDPart3tmp2, hDD_dKOD120)));
        const REAL_SIMD_ARRAY __RHS_exp_16 = FusedMulAddSIMD(
            FDPart3tmp210,
            FusedMulAddSIMD(
                FDPart3tmp416, FusedMulAddSIMD(FDPart3_Rational_2_3, FDPart3tmp64, MulSIMD(FDPart3_Rational_2_3, FDPart3tmp67)),
                FusedMulAddSIMD(FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp64, hDD_dupD221, FDPart3tmp120),
                                FusedMulSubSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp64, hDD_dupD220, FDPart3tmp102),
                                                MulSIMD(FDPart3tmp41, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp15), MulSIMD(aDD22, alpha)))))),
            FusedMulAddSIMD(FDPart3tmp4, hDD_dKOD221, MulSIMD(FDPart3tmp2, hDD_dKOD220)));
        const REAL_SIMD_ARRAY __RHS_exp_17 = FusedMulAddSIMD(
            FDPart3tmp4, lambdaU_dKOD01,
            FusedMulAddSIMD(
                f0_of_xx0__D0,
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp72, FDPart3tmp90), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(alpha, trK_dD1)),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp72, FDPart3tmp93), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(alpha, trK_dD0)),
                        FusedMulAddSIMD(
                            FDPart3tmp94,
                            FusedMulAddSIMD(
                                FDPart3tmp13,
                                FusedMulSubSIMD(FDPart3tmp0, f0_of_xx0__DDD000, MulSIMD(FDPart3tmp6, MulSIMD(f0_of_xx0__DD00, f0_of_xx0__DD00))),
                                FusedMulAddSIMD(FDPart3tmp0, MulSIMD(FDPart3tmp26, f0_of_xx0__DD00), FDPart3tmp245)),
                            FusedMulAddSIMD(
                                FDPart3tmp97,
                                FusedMulAddSIMD(FDPart3tmp0, vetU_dDD011,
                                                FusedMulAddSIMD(FDPart3tmp226, FDPart3tmp26,
                                                                NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp0, vetU_dD11),
                                                                                   SubSIMD(FDPart3tmp230, FDPart3tmp13)))),
                                FusedMulAddSIMD(
                                    FDPart3tmp86,
                                    FusedMulAddSIMD(
                                        FDPart3tmp230, FDPart3tmp41,
                                        FusedMulAddSIMD(
                                            FDPart3tmp235, FDPart3tmp26,
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp0, f1_of_xx1),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1__D1, vetU1)),
                                                            FusedMulSubSIMD(FDPart3tmp118, FDPart3tmp234, MulSIMD(FDPart3tmp13, FDPart3tmp41))))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp91,
                                        FusedMulAddSIMD(FDPart3tmp218, vetU1,
                                                        FusedMulAddSIMD(FDPart3tmp9,
                                                                        FusedMulAddSIMD(FDPart3tmp166, FDPart3tmp6, FDPart3_NegativeOne_),
                                                                        AddSIMD(FDPart3tmp228, FDPart3tmp9))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp422, FDPart3tmp94,
                                            FusedMulAddSIMD(
                                                FDPart3tmp78,
                                                FusedMulAddSIMD(FDPart3tmp238,
                                                                FusedMulSubSIMD(FDPart3tmp166, MulSIMD(FDPart3tmp41, FDPart3tmp6), FDPart3tmp41),
                                                                FusedMulAddSIMD(FDPart3tmp218, MulSIMD(f1_of_xx1, vetU2),
                                                                                FusedMulAddSIMD(FDPart3tmp185, f1_of_xx1, FDPart3tmp240))),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp323, FDPart3tmp324,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp421, FDPart3tmp91,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp311, FDPart3tmp425,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp316, FDPart3tmp426,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp286,
                                                                    FusedMulAddSIMD(
                                                                        MulSIMD(FDPart3tmp93, aDD01),
                                                                        MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp90),
                                                                                MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_),
                                                                                        MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp288,
                                                                            MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp287),
                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                            MulSIMD(FDPart3tmp41, aDD22))),
                                                                            FusedMulAddSIMD(
                                                                                MulSIMD(FDPart3tmp77, FDPart3tmp90),
                                                                                MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_),
                                                                                                MulSIMD(aDD12, f1_of_xx1))),
                                                                                FusedMulAddSIMD(
                                                                                    f0_of_xx0,
                                                                                    MulSIMD(MulSIMD(FDPart3tmp93, aDD02),
                                                                                            MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp77),
                                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_4,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                            MulSIMD(f0_of_xx0__D0, f1_of_xx1)))),
                                                                                    FusedMulSubSIMD(
                                                                                        MulSIMD(FDPart3tmp288, FDPart3tmp5),
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(aDD00, MulSIMD(FDPart3tmp93, FDPart3tmp93))),
                                                                                        MulSIMD(FDPart3tmp288,
                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp15),
                                                                                                        MulSIMD(FDPart3tmp289, aDD11)))))))),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp293, FDPart3tmp424,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp264, FDPart3tmp429,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp279, FDPart3tmp430,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp259, FDPart3tmp428,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp26, FDPart3tmp419,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp228, FDPart3tmp91,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp240, FDPart3tmp78,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp13,
                                                                                                    FusedMulAddSIMD(FDPart3tmp0, lambdaU_dupD00,
                                                                                                                    MulSIMD(FDPart3tmp218, lambdaU0)),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp216, lambdaU_dupD01,
                                                                                                        FusedMulAddSIMD(
                                                                                                            MulSIMD(FDPart3tmp41, FDPart3tmp72),
                                                                                                            MulSIMD(
                                                                                                                MulSIMD(FDPart3tmp0, FDPart3tmp188),
                                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(FDPart3tmp82,
                                                                                                                                f0_of_xx0))),
                                                                                                            NegFusedMulAddSIMD(
                                                                                                                FDPart3tmp0,
                                                                                                                MulSIMD(FDPart3tmp224, lambdaU1),
                                                                                                                FDPart3tmp349)))))))))))))))))))))))),
                MulSIMD(FDPart3tmp2, lambdaU_dKOD00)));
        const REAL_SIMD_ARRAY __RHS_exp_18 = FusedMulAddSIMD(
            FDPart3tmp4, lambdaU_dKOD11,
            FusedMulAddSIMD(
                f0_of_xx0,
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp72, FDPart3tmp90), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(alpha, trK_dD0)),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp72, FDPart3tmp96), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(alpha, trK_dD1)),
                        FusedMulAddSIMD(
                            FDPart3tmp86,
                            FusedMulAddSIMD(FDPart3tmp19, FDPart3tmp235,
                                            FusedMulAddSIMD(FDPart3tmp41, FDPart3tmp9,
                                                            FusedMulSubSIMD(FDPart3tmp118, FDPart3tmp38, MulSIMD(FDPart3tmp248, FDPart3tmp9)))),
                            FusedMulAddSIMD(
                                FDPart3tmp94,
                                FusedMulAddSIMD(
                                    FDPart3tmp11, MulSIMD(FDPart3tmp19, FDPart3tmp3),
                                    FusedMulAddSIMD(
                                        FDPart3tmp367, vetU1,
                                        FusedMulAddSIMD(
                                            vetU1, FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp367, MulSIMD(FDPart3tmp16, f0_of_xx0__DD00)),
                                            FusedMulAddSIMD(
                                                FDPart3tmp3, vetU_dDD100,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp359, FDPart3tmp9,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp16,
                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, vetU_dD10)),
                                                        FusedMulSubSIMD(FDPart3tmp19, FDPart3tmp310, MulSIMD(FDPart3tmp276, f0_of_xx0__DD00)))))))),
                                FusedMulAddSIMD(
                                    FDPart3tmp422, FDPart3tmp91,
                                    FusedMulAddSIMD(
                                        FDPart3tmp83,
                                        FusedMulAddSIMD(FDPart3tmp238, FDPart3tmp248,
                                                        FusedMulAddSIMD(FDPart3tmp238, NegFusedMulSubSIMD(f1_of_xx1, f1_of_xx1__DD11, FDPart3tmp248),
                                                                        NegFusedMulAddSIMD(FDPart3tmp185, f1_of_xx1, FDPart3tmp364))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp384, FDPart3tmp430,
                                            FusedMulAddSIMD(
                                                FDPart3tmp421, FDPart3tmp97,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp379, FDPart3tmp425,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp381, FDPart3tmp424,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp375, FDPart3tmp428,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp377, FDPart3tmp426,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp364, FDPart3tmp83,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp372, FDPart3tmp429,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp286, FDPart3tmp324,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp323,
                                                                                FusedMulAddSIMD(
                                                                                    MulSIMD(FDPart3tmp96, aDD01),
                                                                                    MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp90),
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_),
                                                                                                    MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp295,
                                                                                        MulSIMD(
                                                                                            MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                    MulSIMD(FDPart3tmp41, aDD22))),
                                                                                        FusedMulAddSIMD(
                                                                                            MulSIMD(FDPart3tmp82, FDPart3tmp96),
                                                                                            MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_4,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                            MulSIMD(aDD12, f1_of_xx1))),
                                                                                            FusedMulAddSIMD(
                                                                                                f0_of_xx0,
                                                                                                MulSIMD(MulSIMD(FDPart3tmp90, aDD02),
                                                                                                        MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp82),
                                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_4,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(f0_of_xx0__D0,
                                                                                                                                f1_of_xx1)))),
                                                                                                FusedMulSubSIMD(
                                                                                                    MulSIMD(FDPart3tmp288, FDPart3tmp289),
                                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                            MulSIMD(FDPart3tmp5, aDD00)),
                                                                                                    MulSIMD(
                                                                                                        FDPart3tmp288,
                                                                                                        MulSIMD(
                                                                                                            MulSIMD(FDPart3_Integer_2, FDPart3tmp15),
                                                                                                            MulSIMD(aDD11,
                                                                                                                    MulSIMD(FDPart3tmp96,
                                                                                                                            FDPart3tmp96))))))))),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp19, FDPart3tmp419,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp276, lambdaU_dupD11,
                                                                                        AddSIMD(
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp13,
                                                                                                FusedMulAddSIMD(FDPart3tmp3, lambdaU_dupD10,
                                                                                                                MulSIMD(FDPart3tmp350, lambdaU1)),
                                                                                                FDPart3tmp385),
                                                                                            AddSIMD(
                                                                                                AddSIMD(FDPart3tmp358, FDPart3tmp361),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp77,
                                                                                                    MulSIMD(
                                                                                                        MulSIMD(FDPart3tmp34, FDPart3tmp72),
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(f1_of_xx1, f1_of_xx1__D1))),
                                                                                                    NegFusedMulAddSIMD(
                                                                                                        FDPart3tmp16, MulSIMD(lambdaU1, vetU_dD11),
                                                                                                        FDPart3tmp351)))))))))))))))))))))),
                MulSIMD(FDPart3tmp2, lambdaU_dKOD10)));
        const REAL_SIMD_ARRAY __RHS_exp_19 = FusedMulAddSIMD(
            FDPart3tmp207,
            FusedMulAddSIMD(
                FDPart3_Integer_2, MulSIMD(FDPart3tmp398, FDPart3tmp91),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp72, FDPart3tmp77), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(alpha, trK_dD0)),
                    FusedMulAddSIMD(
                        FDPart3tmp94,
                        FusedMulAddSIMD(
                            FDPart3tmp11, MulSIMD(FDPart3tmp3, FDPart3tmp34),
                            FusedMulAddSIMD(
                                FDPart3tmp399, vetU2,
                                FusedMulAddSIMD(
                                    vetU2, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp28, FDPart3tmp367), FDPart3tmp399),
                                    FusedMulAddSIMD(
                                        FDPart3tmp29, vetU_dDD200,
                                        FusedMulAddSIMD(
                                            FDPart3tmp310, FDPart3tmp34,
                                            FusedMulAddSIMD(FDPart3tmp28, MulSIMD(FDPart3tmp367, vetU2),
                                                            FusedMulSubSIMD(FDPart3tmp238, FDPart3tmp359,
                                                                            MulSIMD(FDPart3tmp28, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp16),
                                                                                                          MulSIMD(f0_of_xx0__D0, vetU_dD20)))))))))),
                        FusedMulAddSIMD(
                            FDPart3tmp97,
                            FusedMulAddSIMD(
                                FDPart3_Integer_2, MulSIMD(FDPart3tmp188, FDPart3tmp39),
                                FusedMulAddSIMD(
                                    FDPart3tmp29, vetU_dDD211,
                                    FusedMulAddSIMD(
                                        vetU2,
                                        FusedMulSubSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp248), MulSIMD(FDPart3tmp249, FDPart3tmp3),
                                                        MulSIMD(FDPart3tmp3, MulSIMD(FDPart3tmp42, f1_of_xx1__DD11))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp226, FDPart3tmp34,
                                            FusedMulAddSIMD(FDPart3tmp238, FDPart3tmp391,
                                                            FusedMulAddSIMD(FDPart3tmp185, MulSIMD(FDPart3tmp248, FDPart3tmp249),
                                                                            NegFusedMulAddSIMD(FDPart3tmp42,
                                                                                               MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp3),
                                                                                                       MulSIMD(f1_of_xx1__D1, vetU_dD21)),
                                                                                               FDPart3tmp238))))))),
                            FusedMulAddSIMD(
                                FDPart3tmp83, FusedMulAddSIMD(FDPart3tmp388, FDPart3tmp9, SubSIMD(FDPart3tmp389, FDPart3tmp9)),
                                FusedMulAddSIMD(
                                    FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp29, lambdaU_dupD21, MulSIMD(FDPart3tmp387, lambdaU2)),
                                    FusedMulAddSIMD(
                                        FDPart3tmp421, FDPart3tmp83,
                                        FusedMulAddSIMD(
                                            FDPart3tmp422, FDPart3tmp78,
                                            FusedMulAddSIMD(
                                                FDPart3tmp411, FDPart3tmp428,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp413, FDPart3tmp429,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp407, FDPart3tmp430,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp408, FDPart3tmp425,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp404, FDPart3tmp424,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp406, FDPart3tmp426,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp323,
                                                                        FusedMulAddSIMD(
                                                                            MulSIMD(FDPart3tmp90, aDD01),
                                                                            MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp82),
                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                            MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                            FusedMulAddSIMD(
                                                                                f0_of_xx0,
                                                                                MulSIMD(
                                                                                    MulSIMD(FDPart3tmp82, aDD02),
                                                                                    MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp77),
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                    MulSIMD(f0_of_xx0__D0, f1_of_xx1)))),
                                                                                FusedMulAddSIMD(
                                                                                    MulSIMD(FDPart3tmp85, FDPart3tmp96),
                                                                                    MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                    MulSIMD(aDD12, f1_of_xx1))),
                                                                                    FusedMulAddSIMD(
                                                                                        MulSIMD(FDPart3tmp96, aDD01),
                                                                                        MulSIMD(
                                                                                            MulSIMD(FDPart3tmp288, FDPart3tmp77),
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                    MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp77,
                                                                                            MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp5),
                                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                            MulSIMD(FDPart3tmp90, aDD00))),
                                                                                            FusedMulAddSIMD(
                                                                                                MulSIMD(FDPart3tmp41, FDPart3tmp82),
                                                                                                MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(FDPart3tmp85, aDD22))),
                                                                                                FusedMulAddSIMD(
                                                                                                    f0_of_xx0,
                                                                                                    MulSIMD(
                                                                                                        MulSIMD(FDPart3tmp90, aDD02),
                                                                                                        MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp85),
                                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(f0_of_xx0__D0,
                                                                                                                                f1_of_xx1)))),
                                                                                                    FusedMulSubSIMD(
                                                                                                        FDPart3tmp82,
                                                                                                        MulSIMD(
                                                                                                            MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                            FDPart3_NegativeOne_),
                                                                                                                    MulSIMD(FDPart3tmp96, aDD11))),
                                                                                                        MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp295),
                                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3tmp15),
                                                                                                                        MulSIMD(aDD12,
                                                                                                                                f1_of_xx1))))))))))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp34, FDPart3tmp419,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp13,
                                                                                FusedMulAddSIMD(FDPart3tmp29, lambdaU_dupD20,
                                                                                                MulSIMD(FDPart3tmp397, lambdaU2)),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp286,
                                                                                    FusedMulAddSIMD(
                                                                                        MulSIMD(FDPart3tmp90, aDD01),
                                                                                        MulSIMD(
                                                                                            MulSIMD(FDPart3tmp288, FDPart3tmp77),
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                    MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                                        FusedMulAddSIMD(
                                                                                            MulSIMD(FDPart3tmp93, aDD01),
                                                                                            MulSIMD(MulSIMD(FDPart3tmp288, FDPart3tmp82),
                                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                            MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                                                            FusedMulAddSIMD(
                                                                                                MulSIMD(FDPart3tmp85, FDPart3tmp90),
                                                                                                MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(aDD12, f1_of_xx1))),
                                                                                                FusedMulAddSIMD(
                                                                                                    MulSIMD(aDD02, f0_of_xx0),
                                                                                                    MulSIMD(
                                                                                                        MulSIMD(FDPart3tmp287, FDPart3tmp288),
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                MulSIMD(f0_of_xx0__D0, f1_of_xx1))),
                                                                                                    FusedMulAddSIMD(
                                                                                                        MulSIMD(FDPart3tmp41, FDPart3tmp77),
                                                                                                        MulSIMD(
                                                                                                            MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                            FDPart3_NegativeOne_),
                                                                                                                    MulSIMD(FDPart3tmp85, aDD22))),
                                                                                                        FusedMulAddSIMD(
                                                                                                            MulSIMD(FDPart3tmp77, FDPart3tmp82),
                                                                                                            MulSIMD(
                                                                                                                MulSIMD(FDPart3tmp15, FDPart3tmp288),
                                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                        MulSIMD(aDD12, f1_of_xx1))),
                                                                                                            FusedMulAddSIMD(
                                                                                                                f0_of_xx0,
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3tmp93, aDD02),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3tmp288,
                                                                                                                                FDPart3tmp85),
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                            MulSIMD(f0_of_xx0__D0,
                                                                                                                                    f1_of_xx1)))),
                                                                                                                FusedMulSubSIMD(
                                                                                                                    FDPart3tmp77,
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3tmp288,
                                                                                                                                FDPart3tmp5),
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3_Integer_2,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                            MulSIMD(FDPart3tmp93,
                                                                                                                                    aDD00))),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3tmp288,
                                                                                                                                FDPart3tmp82),
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(FDPart3_Integer_2,
                                                                                                                                    FDPart3tmp15),
                                                                                                                            MulSIMD(FDPart3tmp90,
                                                                                                                                    aDD11))))))))))),
                                                                                    AddSIMD(AddSIMD(FDPart3tmp396, FDPart3tmp414),
                                                                                            AddSIMD(AddSIMD(FDPart3tmp394, FDPart3tmp395),
                                                                                                    FusedMulAddSIMD(
                                                                                                        MulSIMD(FDPart3tmp72, FDPart3tmp82),
                                                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_,
                                                                                                                        FDPart3_Rational_4_3),
                                                                                                                MulSIMD(alpha, trK_dD1)),
                                                                                                        NegFusedMulAddSIMD(
                                                                                                            FDPart3tmp188,
                                                                                                            MulSIMD(FDPart3tmp3, lambdaU1),
                                                                                                            FDPart3tmp392)))))))))))))))))))))),
            FusedMulAddSIMD(FDPart3tmp4, lambdaU_dKOD21, MulSIMD(FDPart3tmp2, lambdaU_dKOD20)));
        const REAL_SIMD_ARRAY __RHS_exp_20 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp72, FDPart3tmp77),
            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                    MulSIMD(FDPart3tmp98, FusedMulAddSIMD(FDPart3tmp111, alpha_dD0, MulSIMD(FDPart3tmp113, alpha_dD1)))),
            NegFusedMulAddSIMD(
                FDPart3tmp94,
                MulSIMD(FDPart3tmp98,
                        FusedMulAddSIMD(FDPart3tmp175, alpha_dD0,
                                        FusedMulAddSIMD(FDPart3tmp177, alpha_dD1, NegFusedMulAddSIMD(FDPart3tmp139, alpha_dD0, alpha_dDD00)))),
                NegFusedMulAddSIMD(
                    FDPart3tmp97,
                    MulSIMD(FDPart3tmp98,
                            FusedMulAddSIMD(FDPart3tmp160, alpha_dD0,
                                            FusedMulAddSIMD(FDPart3tmp162, alpha_dD1, NegFusedMulAddSIMD(FDPart3tmp136, alpha_dD1, alpha_dDD11)))),
                    FusedMulAddSIMD(
                        FDPart3tmp304, MulSIMD(FDPart3tmp339, alpha),
                        FusedMulAddSIMD(
                            FDPart3tmp314, MulSIMD(FDPart3tmp7, alpha),
                            FusedMulAddSIMD(
                                FDPart3tmp299, MulSIMD(FDPart3tmp341, alpha),
                                FusedMulAddSIMD(
                                    FDPart3tmp301, MulSIMD(FDPart3tmp337, alpha),
                                    FusedMulAddSIMD(
                                        FDPart3tmp184, MulSIMD(FDPart3tmp321, alpha),
                                        FusedMulAddSIMD(
                                            FDPart3tmp197, MulSIMD(FDPart3tmp307, alpha),
                                            FusedMulAddSIMD(
                                                FDPart3tmp9, trK_dupD1,
                                                FusedMulAddSIMD(
                                                    FDPart3_Rational_1_3, MulSIMD(alpha, MulSIMD(trK, trK)),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp433, AddSIMD(FDPart3tmp138, FDPart3tmp432),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp433, AddSIMD(FDPart3tmp141, FDPart3tmp432),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp2, trK_dKOD0,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp4, trK_dKOD1,
                                                                    FusedMulAddSIMD(
                                                                        MulSIMD(FDPart3tmp72, FDPart3tmp82),
                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                MulSIMD(FDPart3tmp98,
                                                                                        FusedMulAddSIMD(FDPart3tmp124, alpha_dD1,
                                                                                                        MulSIMD(FDPart3tmp126, alpha_dD0)))),
                                                                        FusedMulSubSIMD(
                                                                            FDPart3tmp13, trK_dupD0,
                                                                            MulSIMD(FDPart3tmp86,
                                                                                    MulSIMD(FDPart3tmp98,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp131, alpha_dD0,
                                                                                                MulSIMD(FDPart3tmp133, alpha_dD1)))))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_21 = FusedMulAddSIMD(
            FDPart3tmp213, vetU_dKOD01,
            FusedMulAddSIMD(
                f0_of_xx0__D0,
                FusedMulAddSIMD(
                    MulSIMD(FDPart3_Integer_2, FDPart3tmp13), MulSIMD(FDPart3tmp238, FDPart3tmp259),
                    FusedMulAddSIMD(
                        FDPart3tmp275, FDPart3tmp437,
                        FusedMulAddSIMD(
                            FDPart3tmp283, FDPart3tmp439,
                            FusedMulAddSIMD(
                                FDPart3tmp216, vetU_dupD01,
                                FusedMulAddSIMD(
                                    FDPart3tmp268, FDPart3tmp436,
                                    FusedMulAddSIMD(
                                        MulSIMD(FDPart3_Integer_2, FDPart3tmp13), MulSIMD(FDPart3tmp279, FDPart3tmp9),
                                        FusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp264), MulSIMD(FDPart3tmp31, vetU1),
                                                        FusedMulAddSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp0, vetU_dupD00, FDPart3tmp25),
                                                                        FDPart3tmp214)))))))),
                MulSIMD(FDPart3tmp212, vetU_dKOD00)));
        const REAL_SIMD_ARRAY __RHS_exp_22 = FusedMulAddSIMD(
            FDPart3tmp213, vetU_dKOD11,
            FusedMulAddSIMD(
                f0_of_xx0,
                FusedMulAddSIMD(
                    FDPart3tmp377, FDPart3tmp437,
                    FusedMulAddSIMD(
                        FDPart3tmp375, FDPart3tmp441,
                        FusedMulAddSIMD(
                            FDPart3tmp376, FDPart3tmp436,
                            FusedMulAddSIMD(FDPart3tmp276, vetU_dupD11,
                                            FusedMulAddSIMD(FDPart3tmp372, FDPart3tmp440,
                                                            FusedMulAddSIMD(FDPart3tmp378, FDPart3tmp442,
                                                                            FusedMulAddSIMD(FDPart3tmp379, FDPart3tmp439,
                                                                                            FusedMulAddSIMD(FDPart3tmp13,
                                                                                                            FusedMulAddSIMD(FDPart3tmp3, vetU_dupD10,
                                                                                                                            FDPart3tmp18),
                                                                                                            FDPart3tmp281)))))))),
                MulSIMD(FDPart3tmp212, vetU_dKOD10)));
        const REAL_SIMD_ARRAY __RHS_exp_23 = FusedMulAddSIMD(
            FDPart3tmp212, vetU_dKOD20,
            FusedMulAddSIMD(
                FDPart3tmp213, vetU_dKOD21,
                MulSIMD(FDPart3tmp207,
                        FusedMulAddSIMD(
                            FDPart3tmp407, FDPart3tmp442,
                            FusedMulAddSIMD(
                                FDPart3tmp405, FDPart3tmp441,
                                FusedMulAddSIMD(
                                    FDPart3tmp406, FDPart3tmp437,
                                    FusedMulAddSIMD(
                                        FDPart3tmp403, FDPart3tmp440,
                                        FusedMulAddSIMD(
                                            FDPart3tmp404, FDPart3tmp436,
                                            FusedMulAddSIMD(FDPart3tmp408, FDPart3tmp439,
                                                            FusedMulAddSIMD(FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp29, vetU_dupD21, FDPart3tmp187),
                                                                            FusedMulAddSIMD(FDPart3tmp13,
                                                                                            FusedMulAddSIMD(FDPart3tmp29, vetU_dupD20, FDPart3tmp33),
                                                                                            FDPart3tmp261)))))))))));

        WriteSIMD(&rhs_gfs[IDX4(ADD00GF, i0, i1, i2)], __RHS_exp_0);
        WriteSIMD(&rhs_gfs[IDX4(ADD01GF, i0, i1, i2)], __RHS_exp_1);
        WriteSIMD(&rhs_gfs[IDX4(ADD02GF, i0, i1, i2)], __RHS_exp_2);
        WriteSIMD(&rhs_gfs[IDX4(ADD11GF, i0, i1, i2)], __RHS_exp_3);
        WriteSIMD(&rhs_gfs[IDX4(ADD12GF, i0, i1, i2)], __RHS_exp_4);
        WriteSIMD(&rhs_gfs[IDX4(ADD22GF, i0, i1, i2)], __RHS_exp_5);
        WriteSIMD(&rhs_gfs[IDX4(ALPHAGF, i0, i1, i2)], __RHS_exp_6);
        WriteSIMD(&rhs_gfs[IDX4(BETU0GF, i0, i1, i2)], __RHS_exp_7);
        WriteSIMD(&rhs_gfs[IDX4(BETU1GF, i0, i1, i2)], __RHS_exp_8);
        WriteSIMD(&rhs_gfs[IDX4(BETU2GF, i0, i1, i2)], __RHS_exp_9);
        WriteSIMD(&rhs_gfs[IDX4(CFGF, i0, i1, i2)], __RHS_exp_10);
        WriteSIMD(&rhs_gfs[IDX4(HDD00GF, i0, i1, i2)], __RHS_exp_11);
        WriteSIMD(&rhs_gfs[IDX4(HDD01GF, i0, i1, i2)], __RHS_exp_12);
        WriteSIMD(&rhs_gfs[IDX4(HDD02GF, i0, i1, i2)], __RHS_exp_13);
        WriteSIMD(&rhs_gfs[IDX4(HDD11GF, i0, i1, i2)], __RHS_exp_14);
        WriteSIMD(&rhs_gfs[IDX4(HDD12GF, i0, i1, i2)], __RHS_exp_15);
        WriteSIMD(&rhs_gfs[IDX4(HDD22GF, i0, i1, i2)], __RHS_exp_16);
        WriteSIMD(&rhs_gfs[IDX4(LAMBDAU0GF, i0, i1, i2)], __RHS_exp_17);
        WriteSIMD(&rhs_gfs[IDX4(LAMBDAU1GF, i0, i1, i2)], __RHS_exp_18);
        WriteSIMD(&rhs_gfs[IDX4(LAMBDAU2GF, i0, i1, i2)], __RHS_exp_19);
        WriteSIMD(&rhs_gfs[IDX4(TRKGF, i0, i1, i2)], __RHS_exp_20);
        WriteSIMD(&rhs_gfs[IDX4(VETU0GF, i0, i1, i2)], __RHS_exp_21);
        WriteSIMD(&rhs_gfs[IDX4(VETU1GF, i0, i1, i2)], __RHS_exp_22);
        WriteSIMD(&rhs_gfs[IDX4(VETU2GF, i0, i1, i2)], __RHS_exp_23);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += simd_width)
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
