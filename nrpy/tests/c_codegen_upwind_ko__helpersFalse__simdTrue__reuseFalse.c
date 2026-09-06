/*
 * NRPy-Generated GF Access/FD Code, Step 1 of 3:
 * Read gridfunction(s) from main memory and compute FD stencils as needed.
 */
const REAL_SIMD_ARRAY u_i2m3 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1, i2 - 3)]);
const REAL_SIMD_ARRAY u_i2m2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1, i2 - 2)]);
const REAL_SIMD_ARRAY u_i2m1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1, i2 - 1)]);
const REAL_SIMD_ARRAY u_i1m3 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1 - 3, i2)]);
const REAL_SIMD_ARRAY u_i0m2_i1m2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 2, i1 - 2, i2)]);
const REAL_SIMD_ARRAY u_i0m1_i1m2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 1, i1 - 2, i2)]);
const REAL_SIMD_ARRAY u_i1m2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1 - 2, i2)]);
const REAL_SIMD_ARRAY u_i0p1_i1m2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 1, i1 - 2, i2)]);
const REAL_SIMD_ARRAY u_i0p2_i1m2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 2, i1 - 2, i2)]);
const REAL_SIMD_ARRAY u_i0m2_i1m1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 2, i1 - 1, i2)]);
const REAL_SIMD_ARRAY u_i0m1_i1m1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 1, i1 - 1, i2)]);
const REAL_SIMD_ARRAY u_i1m1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1 - 1, i2)]);
const REAL_SIMD_ARRAY u_i0p1_i1m1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 1, i1 - 1, i2)]);
const REAL_SIMD_ARRAY u_i0p2_i1m1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 2, i1 - 1, i2)]);
const REAL_SIMD_ARRAY u_i0m3 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 3, i1, i2)]);
const REAL_SIMD_ARRAY u_i0m2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 2, i1, i2)]);
const REAL_SIMD_ARRAY u_i0m1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 1, i1, i2)]);
const REAL_SIMD_ARRAY u = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1, i2)]);
const REAL_SIMD_ARRAY u_i0p1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 1, i1, i2)]);
const REAL_SIMD_ARRAY u_i0p2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 2, i1, i2)]);
const REAL_SIMD_ARRAY u_i0p3 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 3, i1, i2)]);
const REAL_SIMD_ARRAY u_i0m2_i1p1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 2, i1 + 1, i2)]);
const REAL_SIMD_ARRAY u_i0m1_i1p1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 1, i1 + 1, i2)]);
const REAL_SIMD_ARRAY u_i1p1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1 + 1, i2)]);
const REAL_SIMD_ARRAY u_i0p1_i1p1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 1, i1 + 1, i2)]);
const REAL_SIMD_ARRAY u_i0p2_i1p1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 2, i1 + 1, i2)]);
const REAL_SIMD_ARRAY u_i0m2_i1p2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 2, i1 + 2, i2)]);
const REAL_SIMD_ARRAY u_i0m1_i1p2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 1, i1 + 2, i2)]);
const REAL_SIMD_ARRAY u_i1p2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1 + 2, i2)]);
const REAL_SIMD_ARRAY u_i0p1_i1p2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 1, i1 + 2, i2)]);
const REAL_SIMD_ARRAY u_i0p2_i1p2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 2, i1 + 2, i2)]);
const REAL_SIMD_ARRAY u_i1p3 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1 + 3, i2)]);
const REAL_SIMD_ARRAY u_i2p1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1, i2 + 1)]);
const REAL_SIMD_ARRAY u_i2p2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1, i2 + 2)]);
const REAL_SIMD_ARRAY u_i2p3 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1, i2 + 3)]);
const REAL_SIMD_ARRAY v_i0m3 = ReadSIMD(&auxevol_gfs[IDX4(VGF, i0 - 3, i1, i2)]);
const REAL_SIMD_ARRAY v_i0m2 = ReadSIMD(&auxevol_gfs[IDX4(VGF, i0 - 2, i1, i2)]);
const REAL_SIMD_ARRAY v_i0m1 = ReadSIMD(&auxevol_gfs[IDX4(VGF, i0 - 1, i1, i2)]);
const REAL_SIMD_ARRAY v = ReadSIMD(&auxevol_gfs[IDX4(VGF, i0, i1, i2)]);
const REAL_SIMD_ARRAY v_i0p1 = ReadSIMD(&auxevol_gfs[IDX4(VGF, i0 + 1, i1, i2)]);
const REAL_SIMD_ARRAY v_i0p2 = ReadSIMD(&auxevol_gfs[IDX4(VGF, i0 + 2, i1, i2)]);
const REAL_SIMD_ARRAY v_i0p3 = ReadSIMD(&auxevol_gfs[IDX4(VGF, i0 + 3, i1, i2)]);
const REAL_SIMD_ARRAY w_i2m3 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1, i2 - 3)]);
const REAL_SIMD_ARRAY w_i2m2 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1, i2 - 2)]);
const REAL_SIMD_ARRAY w_i2m1 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1, i2 - 1)]);
const REAL_SIMD_ARRAY w_i1m3 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 - 3, i2)]);
const REAL_SIMD_ARRAY w_i1m2 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 - 2, i2)]);
const REAL_SIMD_ARRAY w_i1m1 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 - 1, i2)]);
const REAL_SIMD_ARRAY w = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1, i2)]);
const REAL_SIMD_ARRAY w_i1p1 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 + 1, i2)]);
const REAL_SIMD_ARRAY w_i1p2 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 + 2, i2)]);
const REAL_SIMD_ARRAY w_i1p3 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 + 3, i2)]);
const REAL_SIMD_ARRAY w_i2p1 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1, i2 + 1)]);
const REAL_SIMD_ARRAY w_i2p2 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1, i2 + 2)]);
const REAL_SIMD_ARRAY w_i2p3 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1, i2 + 3)]);
static const double dblFDPart1_NegativeOne_ = -1.0;
MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart1_NegativeOne_ = ConstSIMD(dblFDPart1_NegativeOne_);

static const double dblFDPart1_Rational_15_64 = 15.0 / 64.0;
const REAL_SIMD_ARRAY FDPart1_Rational_15_64 = ConstSIMD(dblFDPart1_Rational_15_64);

static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

static const double dblFDPart1_Rational_1_144 = 1.0 / 144.0;
const REAL_SIMD_ARRAY FDPart1_Rational_1_144 = ConstSIMD(dblFDPart1_Rational_1_144);

static const double dblFDPart1_Rational_1_18 = 1.0 / 18.0;
const REAL_SIMD_ARRAY FDPart1_Rational_1_18 = ConstSIMD(dblFDPart1_Rational_1_18);

static const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

static const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

static const double dblFDPart1_Rational_1_64 = 1.0 / 64.0;
const REAL_SIMD_ARRAY FDPart1_Rational_1_64 = ConstSIMD(dblFDPart1_Rational_1_64);

static const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

static const double dblFDPart1_Rational_3_32 = 3.0 / 32.0;
const REAL_SIMD_ARRAY FDPart1_Rational_3_32 = ConstSIMD(dblFDPart1_Rational_3_32);

static const double dblFDPart1_Rational_4_9 = 4.0 / 9.0;
const REAL_SIMD_ARRAY FDPart1_Rational_4_9 = ConstSIMD(dblFDPart1_Rational_4_9);

static const double dblFDPart1_Rational_5_16 = 5.0 / 16.0;
const REAL_SIMD_ARRAY FDPart1_Rational_5_16 = ConstSIMD(dblFDPart1_Rational_5_16);

static const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
const REAL_SIMD_ARRAY FDPart1_Rational_5_6 = ConstSIMD(dblFDPart1_Rational_5_6);

const REAL_SIMD_ARRAY FDPart1tmp0 = MulSIMD(FDPart1_Rational_5_6, u);
const REAL_SIMD_ARRAY FDPart1tmp3 = MulSIMD(FDPart1_Rational_5_16, u);
const REAL_SIMD_ARRAY UpwindAlgInputu_ddnD0 = MulSIMD(
    invdxx0,
    FusedMulAddSIMD(FDPart1_Rational_1_2, u_i0m2,
                    FusedMulAddSIMD(FDPart1_Rational_1_4, u_i0p1,
                                    SubSIMD(FDPart1tmp0, FusedMulAddSIMD(FDPart1_Rational_1_12, u_i0m3, MulSIMD(FDPart1_Rational_3_2, u_i0m1))))));
const REAL_SIMD_ARRAY UpwindAlgInputu_ddnD1 = MulSIMD(
    invdxx1,
    FusedMulAddSIMD(FDPart1_Rational_1_2, u_i1m2,
                    FusedMulAddSIMD(FDPart1_Rational_1_4, u_i1p1,
                                    SubSIMD(FDPart1tmp0, FusedMulAddSIMD(FDPart1_Rational_1_12, u_i1m3, MulSIMD(FDPart1_Rational_3_2, u_i1m1))))));
const REAL_SIMD_ARRAY UpwindAlgInputu_ddnD2 = MulSIMD(
    invdxx2,
    FusedMulAddSIMD(FDPart1_Rational_1_2, u_i2m2,
                    FusedMulAddSIMD(FDPart1_Rational_1_4, u_i2p1,
                                    SubSIMD(FDPart1tmp0, FusedMulAddSIMD(FDPart1_Rational_1_12, u_i2m3, MulSIMD(FDPart1_Rational_3_2, u_i2m1))))));
const REAL_SIMD_ARRAY UpwindAlgInputu_dupD0 = MulSIMD(
    invdxx0,
    FusedMulAddSIMD(FDPart1_Rational_3_2, u_i0p1,
                    FusedMulSubSIMD(FDPart1_Rational_1_12, u_i0p3,
                                    FusedMulAddSIMD(FDPart1_Rational_1_2, u_i0p2, FusedMulAddSIMD(FDPart1_Rational_1_4, u_i0m1, FDPart1tmp0)))));
const REAL_SIMD_ARRAY UpwindAlgInputu_dupD1 = MulSIMD(
    invdxx1,
    FusedMulAddSIMD(FDPart1_Rational_3_2, u_i1p1,
                    FusedMulSubSIMD(FDPart1_Rational_1_12, u_i1p3,
                                    FusedMulAddSIMD(FDPart1_Rational_1_2, u_i1p2, FusedMulAddSIMD(FDPart1_Rational_1_4, u_i1m1, FDPart1tmp0)))));
const REAL_SIMD_ARRAY UpwindAlgInputu_dupD2 = MulSIMD(
    invdxx2,
    FusedMulAddSIMD(FDPart1_Rational_3_2, u_i2p1,
                    FusedMulSubSIMD(FDPart1_Rational_1_12, u_i2p3,
                                    FusedMulAddSIMD(FDPart1_Rational_1_2, u_i2p2, FusedMulAddSIMD(FDPart1_Rational_1_4, u_i2m1, FDPart1tmp0)))));
const REAL_SIMD_ARRAY UpwindAlgInputv_ddnD0 =
    MulSIMD(invdxx0,
            FusedMulAddSIMD(FDPart1_Rational_1_4, v_i0p1,
                            FusedMulAddSIMD(FDPart1_Rational_5_6, v,
                                            FusedMulSubSIMD(FDPart1_Rational_1_2, v_i0m2,
                                                            FusedMulAddSIMD(FDPart1_Rational_1_12, v_i0m3, MulSIMD(FDPart1_Rational_3_2, v_i0m1))))));
const REAL_SIMD_ARRAY UpwindAlgInputv_dupD0 = MulSIMD(
    invdxx0, FusedMulAddSIMD(FDPart1_Rational_3_2, v_i0p1,
                             FusedMulSubSIMD(FDPart1_Rational_1_12, v_i0p3,
                                             FusedMulAddSIMD(FDPart1_Rational_1_4, v_i0m1,
                                                             FusedMulAddSIMD(FDPart1_Rational_5_6, v, MulSIMD(FDPart1_Rational_1_2, v_i0p2))))));
const REAL_SIMD_ARRAY UpwindAlgInputw_ddnD1 =
    MulSIMD(invdxx1,
            FusedMulAddSIMD(FDPart1_Rational_1_4, w_i1p1,
                            FusedMulAddSIMD(FDPart1_Rational_5_6, w,
                                            FusedMulSubSIMD(FDPart1_Rational_1_2, w_i1m2,
                                                            FusedMulAddSIMD(FDPart1_Rational_1_12, w_i1m3, MulSIMD(FDPart1_Rational_3_2, w_i1m1))))));
const REAL_SIMD_ARRAY UpwindAlgInputw_dupD1 = MulSIMD(
    invdxx1, FusedMulAddSIMD(FDPart1_Rational_3_2, w_i1p1,
                             FusedMulSubSIMD(FDPart1_Rational_1_12, w_i1p3,
                                             FusedMulAddSIMD(FDPart1_Rational_1_4, w_i1m1,
                                                             FusedMulAddSIMD(FDPart1_Rational_5_6, w, MulSIMD(FDPart1_Rational_1_2, w_i1p2))))));
const REAL_SIMD_ARRAY UpwindControlVectorU0 = u;
const REAL_SIMD_ARRAY UpwindControlVectorU1 = v;
const REAL_SIMD_ARRAY UpwindControlVectorU2 = w;
const REAL_SIMD_ARRAY u_dDD01 = MulSIMD(
    invdxx0,
    MulSIMD(invdxx1,
            FusedMulAddSIMD(
                FDPart1_Rational_1_18,
                AddSIMD(AddSIMD(u_i0m2_i1p1, u_i0p1_i1m2),
                        AddSIMD(u_i0p2_i1m1, SubSIMD(u_i0m1_i1p2, AddSIMD(AddSIMD(u_i0m1_i1m2, u_i0m2_i1m1), AddSIMD(u_i0p1_i1p2, u_i0p2_i1p1))))),
                FusedMulAddSIMD(FDPart1_Rational_4_9, AddSIMD(u_i0p1_i1p1, SubSIMD(u_i0m1_i1m1, AddSIMD(u_i0m1_i1p1, u_i0p1_i1m1))),
                                MulSIMD(FDPart1_Rational_1_144, AddSIMD(u_i0p2_i1p2, SubSIMD(u_i0m2_i1m2, AddSIMD(u_i0m2_i1p2, u_i0p2_i1m2))))))));
const REAL_SIMD_ARRAY u_dKOD0 =
    MulSIMD(invdxx0, FusedMulAddSIMD(FDPart1_Rational_1_64, AddSIMD(u_i0m3, u_i0p3),
                                     FusedMulSubSIMD(FDPart1_Rational_15_64, AddSIMD(u_i0m1, u_i0p1),
                                                     FusedMulAddSIMD(FDPart1_Rational_3_32, AddSIMD(u_i0m2, u_i0p2), FDPart1tmp3))));
const REAL_SIMD_ARRAY u_dKOD1 =
    MulSIMD(invdxx1, FusedMulAddSIMD(FDPart1_Rational_1_64, AddSIMD(u_i1m3, u_i1p3),
                                     FusedMulSubSIMD(FDPart1_Rational_15_64, AddSIMD(u_i1m1, u_i1p1),
                                                     FusedMulAddSIMD(FDPart1_Rational_3_32, AddSIMD(u_i1m2, u_i1p2), FDPart1tmp3))));
const REAL_SIMD_ARRAY u_dKOD2 =
    MulSIMD(invdxx2, FusedMulAddSIMD(FDPart1_Rational_1_64, AddSIMD(u_i2m3, u_i2p3),
                                     FusedMulSubSIMD(FDPart1_Rational_15_64, AddSIMD(u_i2m1, u_i2p1),
                                                     FusedMulAddSIMD(FDPart1_Rational_3_32, AddSIMD(u_i2m2, u_i2p2), FDPart1tmp3))));
const REAL_SIMD_ARRAY w_dKOD2 = MulSIMD(
    invdxx2, FusedMulAddSIMD(FDPart1_Rational_1_64, AddSIMD(w_i2m3, w_i2p3),
                             FusedMulSubSIMD(FDPart1_Rational_15_64, AddSIMD(w_i2m1, w_i2p1),
                                             FusedMulAddSIMD(FDPart1_Rational_3_32, AddSIMD(w_i2m2, w_i2p2), MulSIMD(FDPart1_Rational_5_16, w)))));

/*
 * NRPy-Generated GF Access/FD Code, Step 2 of 3:
 * Implement upwinding algorithm.
 */
MAYBE_UNUSED const double tmp_upwind_Integer_0 = 0.000000000000000000000000000000000;

MAYBE_UNUSED const REAL_SIMD_ARRAY upwind_Integer_0 = ConstSIMD(tmp_upwind_Integer_0);
MAYBE_UNUSED const double tmp_upwind_Integer_1 = 1.000000000000000000000000000000000;

MAYBE_UNUSED const REAL_SIMD_ARRAY upwind_Integer_1 = ConstSIMD(tmp_upwind_Integer_1);
const REAL_SIMD_ARRAY Upwind0 = UPWIND_ALG(UpwindControlVectorU0);
const REAL_SIMD_ARRAY Upwind1 = UPWIND_ALG(UpwindControlVectorU1);
const REAL_SIMD_ARRAY Upwind2 = UPWIND_ALG(UpwindControlVectorU2);
static const double dblFDPart2_NegativeOne_ = -1.0;
MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart2_NegativeOne_ = ConstSIMD(dblFDPart2_NegativeOne_);

const REAL_SIMD_ARRAY u_dupD0 = FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputu_dupD0, UpwindAlgInputu_ddnD0), UpwindAlgInputu_ddnD0);
const REAL_SIMD_ARRAY u_dupD1 = FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputu_dupD1, UpwindAlgInputu_ddnD1), UpwindAlgInputu_ddnD1);
const REAL_SIMD_ARRAY u_dupD2 = FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputu_dupD2, UpwindAlgInputu_ddnD2), UpwindAlgInputu_ddnD2);
const REAL_SIMD_ARRAY v_dupD0 = FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputv_dupD0, UpwindAlgInputv_ddnD0), UpwindAlgInputv_ddnD0);
const REAL_SIMD_ARRAY w_dupD1 = FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputw_dupD1, UpwindAlgInputw_ddnD1), UpwindAlgInputw_ddnD1);

/*
 * NRPy-Generated GF Access/FD Code, Step 3 of 3:
 * Evaluate SymPy expressions and write to main memory.
 */
const REAL_SIMD_ARRAY __RHS_exp_0 = AddSIMD(u_dKOD0, u_dupD0);
const REAL_SIMD_ARRAY __RHS_exp_1 = AddSIMD(u_dKOD1, u_dupD1);
const REAL_SIMD_ARRAY __RHS_exp_2 = AddSIMD(u_dKOD2, u_dupD2);
const REAL_SIMD_ARRAY __RHS_exp_3 = v_dupD0;
const REAL_SIMD_ARRAY __RHS_exp_4 = AddSIMD(w_dKOD2, w_dupD1);
const REAL_SIMD_ARRAY __RHS_exp_5 = u_dDD01;

WriteSIMD(&out[0], __RHS_exp_0);
WriteSIMD(&out[1], __RHS_exp_1);
WriteSIMD(&out[2], __RHS_exp_2);
WriteSIMD(&out[3], __RHS_exp_3);
WriteSIMD(&out[4], __RHS_exp_4);
WriteSIMD(&out[5], __RHS_exp_5);
