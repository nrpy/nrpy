/*
 * NRPy-Generated GF Access/FD Code, Step 1 of 3:
 * Read gridfunction(s) from main memory and compute FD stencils as needed.
 */
const REAL_SIMD_ARRAY u_i2m1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1, i2 - 1)]);
const REAL_SIMD_ARRAY u_i1m1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1 - 1, i2)]);
const REAL_SIMD_ARRAY u_i0m1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 - 1, i1, i2)]);
const REAL_SIMD_ARRAY u = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1, i2)]);
const REAL_SIMD_ARRAY u_i0p1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 1, i1, i2)]);
const REAL_SIMD_ARRAY u_i0p2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 2, i1, i2)]);
const REAL_SIMD_ARRAY u_i0p3 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0 + 3, i1, i2)]);
const REAL_SIMD_ARRAY u_i1p1 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1 + 1, i2)]);
const REAL_SIMD_ARRAY u_i1p2 = ReadSIMD(&auxevol_gfs[IDX4(UGF, i0, i1 + 2, i2)]);
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
const REAL_SIMD_ARRAY w_i1m3 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 - 3, i2)]);
const REAL_SIMD_ARRAY w_i1m2 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 - 2, i2)]);
const REAL_SIMD_ARRAY w_i1m1 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 - 1, i2)]);
const REAL_SIMD_ARRAY w = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1, i2)]);
const REAL_SIMD_ARRAY w_i1p1 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 + 1, i2)]);
const REAL_SIMD_ARRAY w_i1p2 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 + 2, i2)]);
const REAL_SIMD_ARRAY w_i1p3 = ReadSIMD(&auxevol_gfs[IDX4(WGF, i0, i1 + 3, i2)]);
const REAL_SIMD_ARRAY UpwindAlgInputu_dupD0 = SIMD_fd_function_dupD0_fdorder4(u, u_i0m1, u_i0p1, u_i0p2, u_i0p3, invdxx0);
const REAL_SIMD_ARRAY UpwindAlgInputu_dupD1 = SIMD_fd_function_dupD1_fdorder4(u, u_i1m1, u_i1p1, u_i1p2, u_i1p3, invdxx1);
const REAL_SIMD_ARRAY UpwindAlgInputu_dupD2 = SIMD_fd_function_dupD2_fdorder4(u, u_i2m1, u_i2p1, u_i2p2, u_i2p3, invdxx2);
const REAL_SIMD_ARRAY UpwindAlgInputv_ddnD0 = SIMD_fd_function_ddnD0_fdorder4(v, v_i0m1, v_i0m2, v_i0m3, v_i0p1, invdxx0);
const REAL_SIMD_ARRAY UpwindAlgInputv_dupD0 = SIMD_fd_function_dupD0_fdorder4(v, v_i0m1, v_i0p1, v_i0p2, v_i0p3, invdxx0);
const REAL_SIMD_ARRAY UpwindAlgInputw_ddnD1 = SIMD_fd_function_ddnD1_fdorder4(w, w_i1m1, w_i1m2, w_i1m3, w_i1p1, invdxx1);
const REAL_SIMD_ARRAY UpwindAlgInputw_dupD1 = SIMD_fd_function_dupD1_fdorder4(w, w_i1m1, w_i1p1, w_i1p2, w_i1p3, invdxx1);
const REAL_SIMD_ARRAY u_dDD01 = SIMD_fd_function_dDD01_fdorder4(&auxevol_gfs[IDX4(UGF, i0, i1, i2)], Nxx_plus_2NGHOSTS0, invdxx0, invdxx1);
const REAL_SIMD_ARRAY u_dKOD0 = SIMD_fd_function_dKOD0_fdorder4(&auxevol_gfs[IDX4(UGF, i0, i1, i2)], invdxx0);
const REAL_SIMD_ARRAY u_dKOD1 = SIMD_fd_function_dKOD1_fdorder4(&auxevol_gfs[IDX4(UGF, i0, i1, i2)], Nxx_plus_2NGHOSTS0, invdxx1);
const REAL_SIMD_ARRAY u_dKOD2 = SIMD_fd_function_dKOD2_fdorder4(&auxevol_gfs[IDX4(UGF, i0, i1, i2)], Nxx_plus_2NGHOSTS0 *Nxx_plus_2NGHOSTS1, invdxx2);
const REAL_SIMD_ARRAY w_dKOD2 = SIMD_fd_function_dKOD2_fdorder4(&auxevol_gfs[IDX4(WGF, i0, i1, i2)], Nxx_plus_2NGHOSTS0 *Nxx_plus_2NGHOSTS1, invdxx2);
const REAL_SIMD_ARRAY UpwindControlVectorU0 = u;
const REAL_SIMD_ARRAY UpwindControlVectorU1 = v;
const REAL_SIMD_ARRAY UpwindControlVectorU2 = w;

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

static const double dblFDPart2_Rational_16_3 = 16.0 / 3.0;
const REAL_SIMD_ARRAY FDPart2_Rational_16_3 = ConstSIMD(dblFDPart2_Rational_16_3);

const REAL_SIMD_ARRAY u_dupD0 =
    FusedMulAddSIMD(u_dKOD0, FusedMulSubSIMD(FDPart2_Rational_16_3, Upwind0, FDPart2_Rational_16_3), UpwindAlgInputu_dupD0);
const REAL_SIMD_ARRAY u_dupD1 =
    FusedMulAddSIMD(u_dKOD1, FusedMulSubSIMD(FDPart2_Rational_16_3, Upwind1, FDPart2_Rational_16_3), UpwindAlgInputu_dupD1);
const REAL_SIMD_ARRAY u_dupD2 =
    FusedMulAddSIMD(u_dKOD2, FusedMulSubSIMD(FDPart2_Rational_16_3, Upwind2, FDPart2_Rational_16_3), UpwindAlgInputu_dupD2);
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
