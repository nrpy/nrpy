/*
 * NRPy-Generated GF Access/FD Code, Step 1 of 3:
 * Read gridfunction(s) from main memory and compute FD stencils as needed.
 */
const REAL u_i2m3 = auxevol_gfs[IDX4(UGF, i0, i1, i2 - 3)];
const REAL u_i2m2 = auxevol_gfs[IDX4(UGF, i0, i1, i2 - 2)];
const REAL u_i2m1 = auxevol_gfs[IDX4(UGF, i0, i1, i2 - 1)];
const REAL u_i1m3 = auxevol_gfs[IDX4(UGF, i0, i1 - 3, i2)];
const REAL u_i1m2 = auxevol_gfs[IDX4(UGF, i0, i1 - 2, i2)];
const REAL u_i1m1 = auxevol_gfs[IDX4(UGF, i0, i1 - 1, i2)];
const REAL u_i0m3 = auxevol_gfs[IDX4(UGF, i0 - 3, i1, i2)];
const REAL u_i0m2 = auxevol_gfs[IDX4(UGF, i0 - 2, i1, i2)];
const REAL u_i0m1 = auxevol_gfs[IDX4(UGF, i0 - 1, i1, i2)];
const REAL u = auxevol_gfs[IDX4(UGF, i0, i1, i2)];
const REAL u_i0p1 = auxevol_gfs[IDX4(UGF, i0 + 1, i1, i2)];
const REAL u_i0p2 = auxevol_gfs[IDX4(UGF, i0 + 2, i1, i2)];
const REAL u_i0p3 = auxevol_gfs[IDX4(UGF, i0 + 3, i1, i2)];
const REAL u_i1p1 = auxevol_gfs[IDX4(UGF, i0, i1 + 1, i2)];
const REAL u_i1p2 = auxevol_gfs[IDX4(UGF, i0, i1 + 2, i2)];
const REAL u_i1p3 = auxevol_gfs[IDX4(UGF, i0, i1 + 3, i2)];
const REAL u_i2p1 = auxevol_gfs[IDX4(UGF, i0, i1, i2 + 1)];
const REAL u_i2p2 = auxevol_gfs[IDX4(UGF, i0, i1, i2 + 2)];
const REAL u_i2p3 = auxevol_gfs[IDX4(UGF, i0, i1, i2 + 3)];
const REAL v_i0m3 = auxevol_gfs[IDX4(VGF, i0 - 3, i1, i2)];
const REAL v_i0m2 = auxevol_gfs[IDX4(VGF, i0 - 2, i1, i2)];
const REAL v_i0m1 = auxevol_gfs[IDX4(VGF, i0 - 1, i1, i2)];
const REAL v = auxevol_gfs[IDX4(VGF, i0, i1, i2)];
const REAL v_i0p1 = auxevol_gfs[IDX4(VGF, i0 + 1, i1, i2)];
const REAL v_i0p2 = auxevol_gfs[IDX4(VGF, i0 + 2, i1, i2)];
const REAL v_i0p3 = auxevol_gfs[IDX4(VGF, i0 + 3, i1, i2)];
const REAL w_i2m3 = auxevol_gfs[IDX4(WGF, i0, i1, i2 - 3)];
const REAL w_i2m2 = auxevol_gfs[IDX4(WGF, i0, i1, i2 - 2)];
const REAL w_i2m1 = auxevol_gfs[IDX4(WGF, i0, i1, i2 - 1)];
const REAL w_i1m3 = auxevol_gfs[IDX4(WGF, i0, i1 - 3, i2)];
const REAL w_i1m2 = auxevol_gfs[IDX4(WGF, i0, i1 - 2, i2)];
const REAL w_i1m1 = auxevol_gfs[IDX4(WGF, i0, i1 - 1, i2)];
const REAL w = auxevol_gfs[IDX4(WGF, i0, i1, i2)];
const REAL w_i1p1 = auxevol_gfs[IDX4(WGF, i0, i1 + 1, i2)];
const REAL w_i1p2 = auxevol_gfs[IDX4(WGF, i0, i1 + 2, i2)];
const REAL w_i1p3 = auxevol_gfs[IDX4(WGF, i0, i1 + 3, i2)];
const REAL w_i2p1 = auxevol_gfs[IDX4(WGF, i0, i1, i2 + 1)];
const REAL w_i2p2 = auxevol_gfs[IDX4(WGF, i0, i1, i2 + 2)];
const REAL w_i2p3 = auxevol_gfs[IDX4(WGF, i0, i1, i2 + 3)];
const REAL UpwindAlgInputu_dupD0 = fd_function_dupD0_fdorder4(u, u_i0m1, u_i0p1, u_i0p2, u_i0p3, invdxx0);
const REAL UpwindAlgInputu_dupD1 = fd_function_dupD1_fdorder4(u, u_i1m1, u_i1p1, u_i1p2, u_i1p3, invdxx1);
const REAL UpwindAlgInputu_dupD2 = fd_function_dupD2_fdorder4(u, u_i2m1, u_i2p1, u_i2p2, u_i2p3, invdxx2);
const REAL UpwindAlgInputv_ddnD0 = fd_function_ddnD0_fdorder4(v, v_i0m1, v_i0m2, v_i0m3, v_i0p1, invdxx0);
const REAL UpwindAlgInputv_dupD0 = fd_function_dupD0_fdorder4(v, v_i0m1, v_i0p1, v_i0p2, v_i0p3, invdxx0);
const REAL UpwindAlgInputw_ddnD1 = fd_function_ddnD1_fdorder4(w, w_i1m1, w_i1m2, w_i1m3, w_i1p1, invdxx1);
const REAL UpwindAlgInputw_dupD1 = fd_function_dupD1_fdorder4(w, w_i1m1, w_i1p1, w_i1p2, w_i1p3, invdxx1);
const REAL u_dDD01 = fd_function_dDD01_fdorder4(&auxevol_gfs[IDX4(UGF, i0, i1, i2)], Nxx_plus_2NGHOSTS0, invdxx0, invdxx1);
const REAL u_dKOD0 = fd_function_dKOD0_fdorder4(u, u_i0m1, u_i0m2, u_i0m3, u_i0p1, u_i0p2, u_i0p3, invdxx0);
const REAL u_dKOD1 = fd_function_dKOD1_fdorder4(u, u_i1m1, u_i1m2, u_i1m3, u_i1p1, u_i1p2, u_i1p3, invdxx1);
const REAL u_dKOD2 = fd_function_dKOD2_fdorder4(u, u_i2m1, u_i2m2, u_i2m3, u_i2p1, u_i2p2, u_i2p3, invdxx2);
const REAL w_dKOD2 = fd_function_dKOD2_fdorder4(w, w_i2m1, w_i2m2, w_i2m3, w_i2p1, w_i2p2, w_i2p3, invdxx2);
const REAL UpwindControlVectorU0 = u;
const REAL UpwindControlVectorU1 = v;
const REAL UpwindControlVectorU2 = w;

/*
 * NRPy-Generated GF Access/FD Code, Step 2 of 3:
 * Implement upwinding algorithm.
 */
const REAL Upwind0 = UPWIND_ALG(UpwindControlVectorU0);
const REAL Upwind1 = UPWIND_ALG(UpwindControlVectorU1);
const REAL Upwind2 = UPWIND_ALG(UpwindControlVectorU2);
const REAL u_dupD0 = UpwindAlgInputu_dupD0 + u_dKOD0 * ((16.0 / 3.0) * Upwind0 - 16.0 / 3.0);
const REAL u_dupD1 = UpwindAlgInputu_dupD1 + u_dKOD1 * ((16.0 / 3.0) * Upwind1 - 16.0 / 3.0);
const REAL u_dupD2 = UpwindAlgInputu_dupD2 + u_dKOD2 * ((16.0 / 3.0) * Upwind2 - 16.0 / 3.0);
const REAL v_dupD0 = Upwind0 * (-UpwindAlgInputv_ddnD0 + UpwindAlgInputv_dupD0) + UpwindAlgInputv_ddnD0;
const REAL w_dupD1 = Upwind1 * (-UpwindAlgInputw_ddnD1 + UpwindAlgInputw_dupD1) + UpwindAlgInputw_ddnD1;

/*
 * NRPy-Generated GF Access/FD Code, Step 3 of 3:
 * Evaluate SymPy expressions and write to main memory.
 */
out[0] = u_dKOD0 + u_dupD0;
out[1] = u_dKOD1 + u_dupD1;
out[2] = u_dKOD2 + u_dupD2;
out[3] = v_dupD0;
out[4] = w_dKOD2 + w_dupD1;
out[5] = u_dDD01;
