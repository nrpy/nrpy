#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Kernel: auxevol_gfs_single_point_host.
 * Kernel to compute AUXEVOL grid functions at a single point.
 */
static void auxevol_gfs_single_point_host(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xx0,
                                          const REAL xx1, const REAL xx2, REAL *restrict psi_background, REAL *restrict ADD_times_AUU) {
  // Load necessary parameters from params_struct
  const REAL AMPL = params->AMPL;
  const REAL SINHW = params->SINHW;

  // Load necessary parameters from commondata_struct
  const REAL P0_x = commondata->P0_x;
  const REAL P0_y = commondata->P0_y;
  const REAL P0_z = commondata->P0_z;
  const REAL P1_x = commondata->P1_x;
  const REAL P1_y = commondata->P1_y;
  const REAL P1_z = commondata->P1_z;
  const REAL S0_x = commondata->S0_x;
  const REAL S0_y = commondata->S0_y;
  const REAL S0_z = commondata->S0_z;
  const REAL S1_x = commondata->S1_x;
  const REAL S1_y = commondata->S1_y;
  const REAL S1_z = commondata->S1_z;
  const REAL bare_mass_0 = commondata->bare_mass_0;
  const REAL bare_mass_1 = commondata->bare_mass_1;
  const REAL zPunc = commondata->zPunc;

  static const REAL _Integer_1 = 1;
  static const REAL _Rational_1_2 = 1.0 / 2.0;
  static const REAL _Rational_7_3 = 7.0 / 3.0;
  static const REAL _Rational_1_3 = 1.0 / 3.0;
  static const REAL _Rational_4_3 = 4.0 / 3.0;
  static const REAL _Integer_3 = 3;
  static const REAL _Rational_1_4 = 1.0 / 4.0;
  static const REAL _Integer_2 = 2;
  static const REAL _Rational_2_3 = 2.0 / 3.0;
  static const REAL _Rational_1_6 = 1.0 / 6.0;
  static const REAL _Rational_7_6 = 7.0 / 6.0;
  static const REAL _Integer_4 = 4;
  static const REAL _Rational_3_2 = 3.0 / 2.0;
  static const REAL _Integer_6 = 6;
  static const REAL _Rational_7_4 = 7.0 / 4.0;
  const REAL tmp0 = (1.0 / (SINHW));
  const REAL tmp8 = cos(xx2);
  const REAL tmp9 = sin(xx1);
  const REAL tmp12 = sin(xx2);
  const REAL tmp1 = exp(tmp0) - exp(-tmp0);
  const REAL tmp2 = (1.0 / (tmp1));
  const REAL tmp4 = exp(tmp0 * xx0) - exp(-tmp0 * xx0);
  const REAL tmp5 = AMPL * tmp2 * tmp4;
  const REAL tmp10 = ((AMPL) * (AMPL)) * ((tmp4) * (tmp4)) * ((tmp9) * (tmp9)) / ((tmp1) * (tmp1));
  const REAL tmp6 = tmp5 * cos(xx1);
  const REAL tmp11 = tmp10 * ((tmp8) * (tmp8));
  const REAL tmp13 = tmp10 * ((tmp12) * (tmp12));
  const REAL tmp76 = _Integer_6 * tmp10 * tmp12 * tmp8;
  const REAL tmp7 = tmp6 + zPunc;
  const REAL tmp16 = tmp6 - zPunc;
  const REAL tmp19 = tmp5 * tmp8 * tmp9;
  const REAL tmp21 = tmp12 * tmp5 * tmp9;
  const REAL tmp15 = tmp11 + tmp13 + ((tmp7) * (tmp7));
  const REAL tmp17 = tmp11 + tmp13 + ((tmp16) * (tmp16));
  const REAL tmp30 = S1_x * tmp7 - S1_z * tmp19;
  const REAL tmp32 = S0_x * tmp16 - S0_z * tmp19;
  const REAL tmp52 = AMPL * S1_y * tmp2 * tmp4 * tmp8 * tmp9 - S1_x * tmp21;
  const REAL tmp55 = AMPL * S0_y * tmp2 * tmp4 * tmp8 * tmp9 - S0_x * tmp21;
  const REAL tmp85 = _Rational_1_4 * tmp21;
  const REAL tmp88 = _Rational_1_4 * tmp19;
  const REAL tmp20 = pow(tmp17, -5.0 / 2.0);
  const REAL tmp23 = -S0_y * tmp16 + S0_z * tmp21;
  const REAL tmp25 = pow(tmp15, -5.0 / 2.0);
  const REAL tmp27 = -S1_y * tmp7 + S1_z * tmp21;
  const REAL tmp35 = pow(tmp15, -3.0 / 2.0);
  const REAL tmp37 = pow(tmp17, -3.0 / 2.0);
  const REAL tmp43 = -P1_x * tmp19 - P1_y * tmp21 - P1_z * tmp7;
  const REAL tmp47 = -P0_x * tmp19 - P0_y * tmp21 - P0_z * tmp16;
  const REAL tmp39 = P0_z * tmp16 * tmp37 + P1_z * tmp35 * tmp7;
  const REAL tmp44 = tmp25 * tmp43;
  const REAL tmp48 = tmp20 * tmp47;
  const REAL tmp53 = -tmp25 * tmp7;
  const REAL tmp56 = -tmp16 * tmp20;
  const REAL tmp63 = -tmp16 * tmp37;
  const REAL tmp64 = -tmp35 * tmp7;
  const REAL tmp66 = P0_y * tmp21 * tmp37;
  const REAL tmp67 = P1_y * tmp21 * tmp35;
  const REAL tmp69 = P0_x * tmp19 * tmp37;
  const REAL tmp70 = P1_x * tmp19 * tmp35;
  const REAL tmp29 = tmp19 * tmp20 * tmp23 + tmp19 * tmp25 * tmp27;
  const REAL tmp34 = tmp20 * tmp21 * tmp32 + tmp21 * tmp25 * tmp30;
  const REAL tmp49 = _Integer_2 * tmp13 * tmp44 + _Integer_2 * tmp13 * tmp48;
  const REAL tmp51 = _Integer_2 * tmp11 * tmp44 + _Integer_2 * tmp11 * tmp48;
  const REAL tmp58 = tmp43 * tmp53;
  const REAL tmp61 = tmp47 * tmp56;
  const REAL tmp65 = P0_z * tmp63 + P1_z * tmp64;
  const REAL tmp68 = -tmp66 - tmp67;
  const REAL tmp71 = -tmp69 - tmp70;
  const REAL tmp74 = _Integer_3 * (-_Rational_1_4 * _Rational_2_3 * tmp16 * tmp61 - _Rational_1_4 * _Rational_2_3 * tmp58 * tmp7 -
                                   _Rational_2_3 * tmp52 * tmp53 - _Rational_2_3 * tmp55 * tmp56) +
                     _Rational_1_6 * tmp39 + _Rational_7_6 * tmp65;
  *psi_background = _Integer_1 + _Rational_1_2 * (bare_mass_0 / sqrt(tmp17) + bare_mass_1 / sqrt(tmp15));
  *ADD_times_AUU = _Integer_2 * (((_Integer_3 * (-tmp19 * tmp20 * tmp32 - tmp19 * tmp25 * tmp30 - tmp20 * tmp21 * tmp23 - tmp21 * tmp25 * tmp27) +
                                   _Rational_1_4 * (-tmp44 * tmp76 - tmp48 * tmp76) +
                                   _Rational_3_2 * (P0_x * tmp21 * tmp37 + P0_y * tmp19 * tmp37 + P1_x * tmp21 * tmp35 + P1_y * tmp19 * tmp35)) *
                                  (_Integer_3 * (-tmp19 * tmp20 * tmp32 - tmp19 * tmp25 * tmp30 - tmp20 * tmp21 * tmp23 - tmp21 * tmp25 * tmp27) +
                                   _Rational_1_4 * (-tmp44 * tmp76 - tmp48 * tmp76) +
                                   _Rational_3_2 * (P0_x * tmp21 * tmp37 + P0_y * tmp19 * tmp37 + P1_x * tmp21 * tmp35 + P1_y * tmp19 * tmp35))) +
                                 ((_Integer_3 * (-tmp16 * tmp48 * tmp85 - tmp20 * tmp21 * tmp55 - tmp21 * tmp25 * tmp52 + tmp30 * tmp53 +
                                                 tmp32 * tmp56 - tmp44 * tmp7 * tmp85 + tmp58 * tmp85 + tmp61 * tmp85) +
                                   _Rational_1_4 * (-P0_y * tmp16 * tmp37 - P1_y * tmp35 * tmp7) +
                                   _Rational_3_2 * (P0_z * tmp21 * tmp37 + P1_z * tmp21 * tmp35) + _Rational_7_4 * (-P0_y * tmp63 - P1_y * tmp64)) *
                                  (_Integer_3 * (-tmp16 * tmp48 * tmp85 - tmp20 * tmp21 * tmp55 - tmp21 * tmp25 * tmp52 + tmp30 * tmp53 +
                                                 tmp32 * tmp56 - tmp44 * tmp7 * tmp85 + tmp58 * tmp85 + tmp61 * tmp85) +
                                   _Rational_1_4 * (-P0_y * tmp16 * tmp37 - P1_y * tmp35 * tmp7) +
                                   _Rational_3_2 * (P0_z * tmp21 * tmp37 + P1_z * tmp21 * tmp35) + _Rational_7_4 * (-P0_y * tmp63 - P1_y * tmp64))) +
                                 ((_Integer_3 * (-tmp16 * tmp48 * tmp88 - tmp19 * tmp20 * tmp55 - tmp19 * tmp25 * tmp52 + tmp23 * tmp56 +
                                                 tmp27 * tmp53 - tmp44 * tmp7 * tmp88 + tmp58 * tmp88 + tmp61 * tmp88) +
                                   _Rational_1_4 * (-P0_x * tmp16 * tmp37 - P1_x * tmp35 * tmp7) +
                                   _Rational_3_2 * (P0_z * tmp19 * tmp37 + P1_z * tmp19 * tmp35) + _Rational_7_4 * (-P0_x * tmp63 - P1_x * tmp64)) *
                                  (_Integer_3 * (-tmp16 * tmp48 * tmp88 - tmp19 * tmp20 * tmp55 - tmp19 * tmp25 * tmp52 + tmp23 * tmp56 +
                                                 tmp27 * tmp53 - tmp44 * tmp7 * tmp88 + tmp58 * tmp88 + tmp61 * tmp88) +
                                   _Rational_1_4 * (-P0_x * tmp16 * tmp37 - P1_x * tmp35 * tmp7) +
                                   _Rational_3_2 * (P0_z * tmp19 * tmp37 + P1_z * tmp19 * tmp35) + _Rational_7_4 * (-P0_x * tmp63 - P1_x * tmp64)))) +
                   ((_Integer_2 * (tmp29 + tmp66 + tmp67) - _Integer_4 * tmp34 +
                     _Rational_1_4 * (-_Integer_4 * tmp13 * tmp44 - _Integer_4 * tmp13 * tmp48 + tmp51) + tmp71 + tmp74) *
                    (_Integer_2 * (tmp29 + tmp66 + tmp67) - _Integer_4 * tmp34 +
                     _Rational_1_4 * (-_Integer_4 * tmp13 * tmp44 - _Integer_4 * tmp13 * tmp48 + tmp51) + tmp71 + tmp74)) +
                   ((_Integer_2 * (tmp34 + tmp69 + tmp70) - _Integer_4 * tmp29 +
                     _Rational_1_4 * (-_Integer_4 * tmp11 * tmp44 - _Integer_4 * tmp11 * tmp48 + tmp49) + tmp68 + tmp74) *
                    (_Integer_2 * (tmp34 + tmp69 + tmp70) - _Integer_4 * tmp29 +
                     _Rational_1_4 * (-_Integer_4 * tmp11 * tmp44 - _Integer_4 * tmp11 * tmp48 + tmp49) + tmp68 + tmp74)) +
                   ((_Integer_2 * (tmp29 + tmp34) - _Rational_1_3 * tmp39 + _Rational_1_4 * (tmp49 + tmp51) +
                     _Rational_4_3 * (_Integer_3 * _Rational_1_4 * tmp16 * tmp61 + _Integer_3 * _Rational_1_4 * tmp58 * tmp7 +
                                      _Integer_3 * tmp52 * tmp53 + _Integer_3 * tmp55 * tmp56) -
                     _Rational_7_3 * tmp65 + tmp68 + tmp71) *
                    (_Integer_2 * (tmp29 + tmp34) - _Rational_1_3 * tmp39 + _Rational_1_4 * (tmp49 + tmp51) +
                     _Rational_4_3 * (_Integer_3 * _Rational_1_4 * tmp16 * tmp61 + _Integer_3 * _Rational_1_4 * tmp58 * tmp7 +
                                      _Integer_3 * tmp52 * tmp53 + _Integer_3 * tmp55 * tmp56) -
                     _Rational_7_3 * tmp65 + tmp68 + tmp71));
} // END FUNCTION: auxevol_gfs_single_point_host

/**
 * Kernel: auxevol_gfs_all_points_host.
 * Kernel to initialize auxillary grid functions at all grid points.
 */
static void auxevol_gfs_all_points_host(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL *restrict x0,
                                        const REAL *restrict x1, const REAL *restrict x2, REAL *restrict in_gfs) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

#pragma omp parallel for
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        auxevol_gfs_single_point_host(commondata, params, xx0, xx1, xx2, &in_gfs[IDX4(PSI_BACKGROUNDGF, i0, i1, i2)],
                                      &in_gfs[IDX4(ADD_TIMES_AUUGF, i0, i1, i2)]);

      } // END LOOP: for i0 over [0, Nxx_plus_2NGHOSTS0)
    } // END LOOP: for i1 over [0, Nxx_plus_2NGHOSTS1)
  } // END LOOP: for i2 over [0, Nxx_plus_2NGHOSTS2)
} // END FUNCTION: auxevol_gfs_all_points_host

/**
 * Kernel: variable_wavespeed_gfs_all_points_host.
 * Kernel to compute variable wavespeed for all grids based on local grid spacing.
 */
static void variable_wavespeed_gfs_all_points_host(const params_struct *restrict params, const REAL *restrict x0, const REAL *restrict x1,
                                                   const REAL *restrict x2, REAL *restrict in_gfs, const REAL dt,
                                                   const REAL MINIMUM_GLOBAL_WAVESPEED) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

  // Load necessary parameters from params_struct
  const REAL AMPL = params->AMPL;
  const REAL SINHW = params->SINHW;

#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        /*
         *  Original SymPy expressions:
         *  "[const REAL dsmin0 = AMPL*params->dxx0*(exp(xx0/SINHW)/SINHW + exp(-xx0/SINHW)/SINHW)/(exp(1/SINHW) - exp(-1/SINHW))]"
         *  "[const REAL dsmin1 = AMPL*params->dxx1*(exp(xx0/SINHW) - exp(-xx0/SINHW))/(exp(1/SINHW) - exp(-1/SINHW))]"
         *  "[const REAL dsmin2 = AMPL*params->dxx2*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)/(exp(1/SINHW) - exp(-1/SINHW))]"
         */
        const REAL tmp0 = (1.0 / (SINHW));
        const REAL tmp4 = AMPL / (exp(tmp0) - exp(-tmp0));
        const REAL tmp2 = exp(tmp0 * xx0);
        const REAL tmp3 = exp(-tmp0 * xx0);
        const REAL tmp5 = tmp4 * (tmp2 - tmp3);
        const REAL dsmin0 = params->dxx0 * tmp4 * (tmp0 * tmp2 + tmp0 * tmp3);
        const REAL dsmin1 = params->dxx1 * tmp5;
        const REAL dsmin2 = params->dxx2 * tmp5 * sin(xx1);

        // Set local wavespeed
        in_gfs[IDX4(VARIABLE_WAVESPEEDGF, i0, i1, i2)] = MINIMUM_GLOBAL_WAVESPEED * NRPYMIN(dsmin0, NRPYMIN(dsmin1, dsmin2)) / dt;

      } // END LOOP: for i0 over [NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS)
    } // END LOOP: for i1 over [NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS)
  } // END LOOP: for i2 over [NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
} // END FUNCTION: variable_wavespeed_gfs_all_points_host

/**
 * Call functions that set up all AUXEVOL gridfunctions.
 */
void auxevol_gfs_set_to_constant__rfm__HoleySinhSpherical(commondata_struct *restrict commondata, params_struct *restrict params,
                                                          REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs) {
#include "set_CodeParameters.h"
  REAL *auxevol_gfs = gridfuncs->auxevol_gfs;
  REAL *x0 = xx[0];
  REAL *x1 = xx[1];
  REAL *x2 = xx[2];

  // Set up variable wavespeed
  variable_wavespeed_gfs_all_points_host(params, x0, x1, x2, auxevol_gfs, dt, MINIMUM_GLOBAL_WAVESPEED);

  // Set up all other AUXEVOL gridfunctions
  auxevol_gfs_all_points_host(commondata, params, x0, x1, x2, auxevol_gfs);

} // END FUNCTION: auxevol_gfs_set_to_constant__rfm__HoleySinhSpherical
