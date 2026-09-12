static const REAL FDPart1_Rational_4_9 = 4.0 / 9.0;
static const REAL FDPart1_Rational_1_18 = 1.0 / 18.0;
static const REAL FDPart1_Rational_1_144 = 1.0 / 144.0;
static const REAL FDPart1_Rational_5_16 = 5.0 / 16.0;
static const REAL FDPart1_Rational_3_32 = 3.0 / 32.0;
static const REAL FDPart1_Rational_1_64 = 1.0 / 64.0;
static const REAL FDPart1_Rational_15_64 = 15.0 / 64.0;
static const REAL FDPart1_Rational_1_2 = 1.0 / 2.0;
static const REAL FDPart1_Rational_3_2 = 3.0 / 2.0;
static const REAL FDPart1_Rational_1_12 = 1.0 / 12.0;
static const REAL FDPart1_Rational_1_4 = 1.0 / 4.0;
static const REAL FDPart1_Rational_5_6 = 5.0 / 6.0;
const REAL u_i2m3 = auxevol_gfs[IDX4(UGF, i0, i1, i2 - 3)];
const REAL u_i2m2 = auxevol_gfs[IDX4(UGF, i0, i1, i2 - 2)];
const REAL u_i2m1 = auxevol_gfs[IDX4(UGF, i0, i1, i2 - 1)];
const REAL u_i1m3 = auxevol_gfs[IDX4(UGF, i0, i1 - 3, i2)];
const REAL u_i0m2_i1m2 = auxevol_gfs[IDX4(UGF, i0 - 2, i1 - 2, i2)];
const REAL u_i0m1_i1m2 = auxevol_gfs[IDX4(UGF, i0 - 1, i1 - 2, i2)];
const REAL u_i1m2 = auxevol_gfs[IDX4(UGF, i0, i1 - 2, i2)];
const REAL u_i0p1_i1m2 = auxevol_gfs[IDX4(UGF, i0 + 1, i1 - 2, i2)];
const REAL u_i0p2_i1m2 = auxevol_gfs[IDX4(UGF, i0 + 2, i1 - 2, i2)];
const REAL u_i0m2_i1m1 = auxevol_gfs[IDX4(UGF, i0 - 2, i1 - 1, i2)];
const REAL u_i0m1_i1m1 = auxevol_gfs[IDX4(UGF, i0 - 1, i1 - 1, i2)];
const REAL u_i1m1 = auxevol_gfs[IDX4(UGF, i0, i1 - 1, i2)];
const REAL u_i0p1_i1m1 = auxevol_gfs[IDX4(UGF, i0 + 1, i1 - 1, i2)];
const REAL u_i0p2_i1m1 = auxevol_gfs[IDX4(UGF, i0 + 2, i1 - 1, i2)];
const REAL u_i0m3 = auxevol_gfs[IDX4(UGF, i0 - 3, i1, i2)];
const REAL u_i0m2 = auxevol_gfs[IDX4(UGF, i0 - 2, i1, i2)];
const REAL u_i0m1 = auxevol_gfs[IDX4(UGF, i0 - 1, i1, i2)];
const REAL u = auxevol_gfs[IDX4(UGF, i0, i1, i2)];
const REAL u_i0p1 = auxevol_gfs[IDX4(UGF, i0 + 1, i1, i2)];
const REAL u_i0p2 = auxevol_gfs[IDX4(UGF, i0 + 2, i1, i2)];
const REAL u_i0p3 = auxevol_gfs[IDX4(UGF, i0 + 3, i1, i2)];
const REAL u_i0m2_i1p1 = auxevol_gfs[IDX4(UGF, i0 - 2, i1 + 1, i2)];
const REAL u_i0m1_i1p1 = auxevol_gfs[IDX4(UGF, i0 - 1, i1 + 1, i2)];
const REAL u_i1p1 = auxevol_gfs[IDX4(UGF, i0, i1 + 1, i2)];
const REAL u_i0p1_i1p1 = auxevol_gfs[IDX4(UGF, i0 + 1, i1 + 1, i2)];
const REAL u_i0p2_i1p1 = auxevol_gfs[IDX4(UGF, i0 + 2, i1 + 1, i2)];
const REAL u_i0m2_i1p2 = auxevol_gfs[IDX4(UGF, i0 - 2, i1 + 2, i2)];
const REAL u_i0m1_i1p2 = auxevol_gfs[IDX4(UGF, i0 - 1, i1 + 2, i2)];
const REAL u_i1p2 = auxevol_gfs[IDX4(UGF, i0, i1 + 2, i2)];
const REAL u_i0p1_i1p2 = auxevol_gfs[IDX4(UGF, i0 + 1, i1 + 2, i2)];
const REAL u_i0p2_i1p2 = auxevol_gfs[IDX4(UGF, i0 + 2, i1 + 2, i2)];
const REAL u_dDD01 =
    invdxx0 * invdxx1 *
    (FDPart1_Rational_1_144 * (u_i0m2_i1m2 - u_i0m2_i1p2 - u_i0p2_i1m2 + u_i0p2_i1p2) +
     FDPart1_Rational_1_18 * (-u_i0m1_i1m2 + u_i0m1_i1p2 - u_i0m2_i1m1 + u_i0m2_i1p1 + u_i0p1_i1m2 - u_i0p1_i1p2 + u_i0p2_i1m1 - u_i0p2_i1p1) +
     FDPart1_Rational_4_9 * (u_i0m1_i1m1 - u_i0m1_i1p1 - u_i0p1_i1m1 + u_i0p1_i1p1));
const REAL u_i1p3 = auxevol_gfs[IDX4(UGF, i0, i1 + 3, i2)];
const REAL u_i2p1 = auxevol_gfs[IDX4(UGF, i0, i1, i2 + 1)];
const REAL u_i2p2 = auxevol_gfs[IDX4(UGF, i0, i1, i2 + 2)];
const REAL u_i2p3 = auxevol_gfs[IDX4(UGF, i0, i1, i2 + 3)];
const REAL v_i0m3 = auxevol_gfs[IDX4(VGF, i0 - 3, i1, i2)];
const REAL v_i0m2 = auxevol_gfs[IDX4(VGF, i0 - 2, i1, i2)];
const REAL v_i0m1 = auxevol_gfs[IDX4(VGF, i0 - 1, i1, i2)];
const REAL v = auxevol_gfs[IDX4(VGF, i0, i1, i2)];
const REAL v_i0p1 = auxevol_gfs[IDX4(VGF, i0 + 1, i1, i2)];
const REAL UpwindAlgInputv_ddnD0 = invdxx0 * (-FDPart1_Rational_1_12 * v_i0m3 + FDPart1_Rational_1_2 * v_i0m2 + FDPart1_Rational_1_4 * v_i0p1 -
                                              FDPart1_Rational_3_2 * v_i0m1 + FDPart1_Rational_5_6 * v);
const REAL v_i0p2 = auxevol_gfs[IDX4(VGF, i0 + 2, i1, i2)];
const REAL v_i0p3 = auxevol_gfs[IDX4(VGF, i0 + 3, i1, i2)];
const REAL UpwindAlgInputv_dupD0 = invdxx0 * (FDPart1_Rational_1_12 * v_i0p3 - FDPart1_Rational_1_2 * v_i0p2 - FDPart1_Rational_1_4 * v_i0m1 +
                                              FDPart1_Rational_3_2 * v_i0p1 - FDPart1_Rational_5_6 * v);
const REAL UpwindControlVectorU1 = v;
const REAL Upwind1 = UPWIND_ALG(UpwindControlVectorU1);
const REAL w_i2m3 = auxevol_gfs[IDX4(WGF, i0, i1, i2 - 3)];
const REAL w_i2m2 = auxevol_gfs[IDX4(WGF, i0, i1, i2 - 2)];
const REAL w_i2m1 = auxevol_gfs[IDX4(WGF, i0, i1, i2 - 1)];
const REAL w_i1m3 = auxevol_gfs[IDX4(WGF, i0, i1 - 3, i2)];
const REAL w_i1m2 = auxevol_gfs[IDX4(WGF, i0, i1 - 2, i2)];
const REAL w_i1m1 = auxevol_gfs[IDX4(WGF, i0, i1 - 1, i2)];
const REAL w = auxevol_gfs[IDX4(WGF, i0, i1, i2)];
const REAL w_i1p1 = auxevol_gfs[IDX4(WGF, i0, i1 + 1, i2)];
const REAL UpwindAlgInputw_ddnD1 = invdxx1 * (-FDPart1_Rational_1_12 * w_i1m3 + FDPart1_Rational_1_2 * w_i1m2 + FDPart1_Rational_1_4 * w_i1p1 -
                                              FDPart1_Rational_3_2 * w_i1m1 + FDPart1_Rational_5_6 * w);
const REAL w_i1p2 = auxevol_gfs[IDX4(WGF, i0, i1 + 2, i2)];
const REAL w_i1p3 = auxevol_gfs[IDX4(WGF, i0, i1 + 3, i2)];
const REAL UpwindAlgInputw_dupD1 = invdxx1 * (FDPart1_Rational_1_12 * w_i1p3 - FDPart1_Rational_1_2 * w_i1p2 - FDPart1_Rational_1_4 * w_i1m1 +
                                              FDPart1_Rational_3_2 * w_i1p1 - FDPart1_Rational_5_6 * w);
const REAL w_dupD1 = Upwind1 * (-UpwindAlgInputw_ddnD1 + UpwindAlgInputw_dupD1) + UpwindAlgInputw_ddnD1;
const REAL FDPart1tmp0 = FDPart1_Rational_5_6 * u;
const REAL w_i2p1 = auxevol_gfs[IDX4(WGF, i0, i1, i2 + 1)];
const REAL w_i2p2 = auxevol_gfs[IDX4(WGF, i0, i1, i2 + 2)];
const REAL w_i2p3 = auxevol_gfs[IDX4(WGF, i0, i1, i2 + 3)];
const REAL w_dKOD2 = invdxx2 * (FDPart1_Rational_15_64 * (w_i2m1 + w_i2p1) + FDPart1_Rational_1_64 * (w_i2m3 + w_i2p3) +
                                FDPart1_Rational_3_32 * (-w_i2m2 - w_i2p2) - FDPart1_Rational_5_16 * w);
const REAL FDPart1tmp3 = -FDPart1_Rational_5_16 * u;
const REAL u_dKOD0 = invdxx0 * (FDPart1_Rational_15_64 * (u_i0m1 + u_i0p1) + FDPart1_Rational_1_64 * (u_i0m3 + u_i0p3) +
                                FDPart1_Rational_3_32 * (-u_i0m2 - u_i0p2) + FDPart1tmp3);
const REAL UpwindAlgInputu_dupD0 = invdxx0 * (FDPart1_Rational_1_12 * u_i0p3 - FDPart1_Rational_1_2 * u_i0p2 - FDPart1_Rational_1_4 * u_i0m1 +
                                              FDPart1_Rational_3_2 * u_i0p1 - FDPart1tmp0);
const REAL u_dKOD1 = invdxx1 * (FDPart1_Rational_15_64 * (u_i1m1 + u_i1p1) + FDPart1_Rational_1_64 * (u_i1m3 + u_i1p3) +
                                FDPart1_Rational_3_32 * (-u_i1m2 - u_i1p2) + FDPart1tmp3);
const REAL u_dKOD2 = invdxx2 * (FDPart1_Rational_15_64 * (u_i2m1 + u_i2p1) + FDPart1_Rational_1_64 * (u_i2m3 + u_i2p3) +
                                FDPart1_Rational_3_32 * (-u_i2m2 - u_i2p2) + FDPart1tmp3);
const REAL UpwindAlgInputu_dupD1 = invdxx1 * (FDPart1_Rational_1_12 * u_i1p3 - FDPart1_Rational_1_2 * u_i1p2 - FDPart1_Rational_1_4 * u_i1m1 +
                                              FDPart1_Rational_3_2 * u_i1p1 - FDPart1tmp0);
const REAL UpwindAlgInputu_dupD2 = invdxx2 * (FDPart1_Rational_1_12 * u_i2p3 - FDPart1_Rational_1_2 * u_i2p2 - FDPart1_Rational_1_4 * u_i2m1 +
                                              FDPart1_Rational_3_2 * u_i2p1 - FDPart1tmp0);
const REAL u_dupD1 = UpwindAlgInputu_dupD1 + u_dKOD1 * ((16.0 / 3.0) * Upwind1 - 16.0 / 3.0);
const REAL UpwindControlVectorU0 = u;
const REAL UpwindControlVectorU2 = w;
const REAL Upwind0 = UPWIND_ALG(UpwindControlVectorU0);
const REAL v_dupD0 = Upwind0 * (-UpwindAlgInputv_ddnD0 + UpwindAlgInputv_dupD0) + UpwindAlgInputv_ddnD0;
const REAL u_dupD0 = UpwindAlgInputu_dupD0 + u_dKOD0 * ((16.0 / 3.0) * Upwind0 - 16.0 / 3.0);
const REAL Upwind2 = UPWIND_ALG(UpwindControlVectorU2);
const REAL u_dupD2 = UpwindAlgInputu_dupD2 + u_dKOD2 * ((16.0 / 3.0) * Upwind2 - 16.0 / 3.0);
out[0] = u_dKOD0 + u_dupD0;
out[1] = u_dKOD1 + u_dupD1;
out[2] = u_dKOD2 + u_dupD2;
out[3] = v_dupD0;
out[4] = w_dKOD2 + w_dupD1;
out[5] = u_dDD01;
