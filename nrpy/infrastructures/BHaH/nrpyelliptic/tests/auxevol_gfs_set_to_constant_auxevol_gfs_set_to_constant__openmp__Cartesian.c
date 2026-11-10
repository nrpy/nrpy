#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Kernel: auxevol_gfs_single_point_host.
 * Kernel to compute AUXEVOL grid functions at a single point.
 */
static void auxevol_gfs_single_point_host(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xx0,
                                          const REAL xx1, const REAL xx2, REAL *restrict psi_background, REAL *restrict ADD_times_AUU) {
  // Load necessary parameters from params_struct

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
  static const REAL _Integer_2 = 2;
  static const REAL _Rational_1_4 = 1.0 / 4.0;
  static const REAL _Rational_7_3 = 7.0 / 3.0;
  static const REAL _Rational_1_3 = 1.0 / 3.0;
  static const REAL _Rational_4_3 = 4.0 / 3.0;
  static const REAL _Integer_3 = 3;
  static const REAL _Integer_4 = 4;
  static const REAL _Rational_2_3 = 2.0 / 3.0;
  static const REAL _Rational_1_6 = 1.0 / 6.0;
  static const REAL _Rational_7_6 = 7.0 / 6.0;
  static const REAL _Rational_3_2 = 3.0 / 2.0;
  static const REAL _Integer_6 = 6;
  static const REAL _Rational_7_4 = 7.0 / 4.0;
  const REAL tmp0 = xx2 + zPunc;
  const REAL tmp1 = ((xx0) * (xx0));
  const REAL tmp2 = ((xx1) * (xx1));
  const REAL tmp5 = xx2 - zPunc;
  const REAL tmp20 = S1_x * tmp0 - S1_z * xx0;
  const REAL tmp22 = S0_x * tmp5 - S0_z * xx0;
  const REAL tmp4 = ((tmp0) * (tmp0)) + tmp1 + tmp2;
  const REAL tmp6 = tmp1 + tmp2 + ((tmp5) * (tmp5));
  const REAL tmp12 = -P1_x * xx0 - P1_y * xx1 - P1_z * tmp0;
  const REAL tmp18 = -P0_x * xx0 - P0_y * xx1 - P0_z * tmp5;
  const REAL tmp28 = -S0_y * tmp5 + S0_z * xx1;
  const REAL tmp31 = -S1_y * tmp0 + S1_z * xx1;
  const REAL tmp8 = pow(tmp4, -5.0 / 2.0);
  const REAL tmp14 = pow(tmp6, -5.0 / 2.0);
  const REAL tmp35 = pow(tmp4, -3.0 / 2.0);
  const REAL tmp36 = pow(tmp6, -3.0 / 2.0);
  const REAL tmp13 = tmp12 * tmp8;
  const REAL tmp19 = tmp14 * tmp18;
  const REAL tmp37 = P0_z * tmp36 * tmp5 + P1_z * tmp0 * tmp35;
  const REAL tmp41 = tmp8 * (-S1_x * xx1 + S1_y * xx0);
  const REAL tmp43 = tmp14 * (-S0_x * xx1 + S0_y * xx0);
  const REAL tmp45 = P0_z * tmp36;
  const REAL tmp46 = P1_z * tmp35;
  const REAL tmp48 = P0_y * tmp36 * xx1;
  const REAL tmp49 = P1_y * tmp35 * xx1;
  const REAL tmp51 = P0_x * tmp36 * xx0;
  const REAL tmp52 = P1_x * tmp35 * xx0;
  const REAL tmp62 = P1_x * tmp35;
  const REAL tmp63 = P1_y * tmp35;
  const REAL tmp64 = P0_x * tmp36;
  const REAL tmp65 = P0_y * tmp36;
  const REAL tmp24 = tmp14 * tmp22 * xx1 + tmp20 * tmp8 * xx1;
  const REAL tmp33 = tmp14 * tmp28 * xx0 + tmp31 * tmp8 * xx0;
  const REAL tmp38 = -((tmp0) * (tmp0)) * tmp13;
  const REAL tmp40 = -tmp19 * ((tmp5) * (tmp5));
  const REAL tmp42 = -tmp0 * tmp41;
  const REAL tmp44 = -tmp43 * tmp5;
  const REAL tmp47 = -tmp0 * tmp46 - tmp45 * tmp5;
  const REAL tmp50 = -tmp48 - tmp49;
  const REAL tmp53 = -tmp51 - tmp52;
  const REAL tmp25 = _Rational_1_4 * tmp13 * tmp2 + _Rational_1_4 * tmp19 * tmp2 + tmp24;
  const REAL tmp34 = _Rational_1_4 * tmp1 * tmp13 + _Rational_1_4 * tmp1 * tmp19 + tmp33;
  const REAL tmp56 = _Integer_3 * _Rational_2_3 * tmp38 + _Integer_3 * _Rational_2_3 * tmp40;
  const REAL tmp57 = _Integer_3 * (-_Rational_2_3 * tmp42 - _Rational_2_3 * tmp44) + _Rational_1_6 * tmp37 + _Rational_7_6 * tmp47;
  *psi_background = _Integer_1 + _Rational_1_2 * (bare_mass_0 / sqrt(tmp6) + bare_mass_1 / sqrt(tmp4));
  *ADD_times_AUU =
      _Integer_2 * (((_Integer_3 * (-tmp14 * tmp22 * xx0 - tmp14 * tmp28 * xx1 - tmp20 * tmp8 * xx0 - tmp31 * tmp8 * xx1) +
                      _Rational_1_4 * (-_Integer_6 * tmp13 * xx0 * xx1 - _Integer_6 * tmp19 * xx0 * xx1) +
                      _Rational_3_2 * (tmp62 * xx1 + tmp63 * xx0 + tmp64 * xx1 + tmp65 * xx0)) *
                     (_Integer_3 * (-tmp14 * tmp22 * xx0 - tmp14 * tmp28 * xx1 - tmp20 * tmp8 * xx0 - tmp31 * tmp8 * xx1) +
                      _Rational_1_4 * (-_Integer_6 * tmp13 * xx0 * xx1 - _Integer_6 * tmp19 * xx0 * xx1) +
                      _Rational_3_2 * (tmp62 * xx1 + tmp63 * xx0 + tmp64 * xx1 + tmp65 * xx0))) +
                    ((_Integer_3 * (-tmp0 * tmp20 * tmp8 - tmp14 * tmp22 * tmp5 - tmp41 * xx1 - tmp43 * xx1) +
                      _Rational_1_4 * (-_Integer_3 * tmp0 * tmp12 * tmp8 * xx1 - _Integer_3 * tmp0 * tmp13 * xx1 -
                                       _Integer_3 * tmp14 * tmp18 * tmp5 * xx1 - _Integer_3 * tmp19 * tmp5 * xx1 - tmp0 * tmp63 - tmp5 * tmp65) +
                      _Rational_3_2 * (tmp45 * xx1 + tmp46 * xx1) + _Rational_7_4 * (tmp0 * tmp63 + tmp5 * tmp65)) *
                     (_Integer_3 * (-tmp0 * tmp20 * tmp8 - tmp14 * tmp22 * tmp5 - tmp41 * xx1 - tmp43 * xx1) +
                      _Rational_1_4 * (-_Integer_3 * tmp0 * tmp12 * tmp8 * xx1 - _Integer_3 * tmp0 * tmp13 * xx1 -
                                       _Integer_3 * tmp14 * tmp18 * tmp5 * xx1 - _Integer_3 * tmp19 * tmp5 * xx1 - tmp0 * tmp63 - tmp5 * tmp65) +
                      _Rational_3_2 * (tmp45 * xx1 + tmp46 * xx1) + _Rational_7_4 * (tmp0 * tmp63 + tmp5 * tmp65))) +
                    ((_Integer_3 * (-tmp0 * tmp31 * tmp8 - tmp14 * tmp28 * tmp5 - tmp41 * xx0 - tmp43 * xx0) +
                      _Rational_1_4 * (-_Integer_3 * tmp0 * tmp12 * tmp8 * xx0 - _Integer_3 * tmp0 * tmp13 * xx0 -
                                       _Integer_3 * tmp14 * tmp18 * tmp5 * xx0 - _Integer_3 * tmp19 * tmp5 * xx0 - tmp0 * tmp62 - tmp5 * tmp64) +
                      _Rational_3_2 * (tmp45 * xx0 + tmp46 * xx0) + _Rational_7_4 * (tmp0 * tmp62 + tmp5 * tmp64)) *
                     (_Integer_3 * (-tmp0 * tmp31 * tmp8 - tmp14 * tmp28 * tmp5 - tmp41 * xx0 - tmp43 * xx0) +
                      _Rational_1_4 * (-_Integer_3 * tmp0 * tmp12 * tmp8 * xx0 - _Integer_3 * tmp0 * tmp13 * xx0 -
                                       _Integer_3 * tmp14 * tmp18 * tmp5 * xx0 - _Integer_3 * tmp19 * tmp5 * xx0 - tmp0 * tmp62 - tmp5 * tmp64) +
                      _Rational_3_2 * (tmp45 * xx0 + tmp46 * xx0) + _Rational_7_4 * (tmp0 * tmp62 + tmp5 * tmp64)))) +
      ((_Integer_2 * (tmp25 + tmp51 + tmp52) - _Integer_4 * tmp33 + _Rational_1_4 * (-_Integer_4 * tmp1 * tmp13 - _Integer_4 * tmp1 * tmp19 - tmp56) +
        tmp50 + tmp57) *
       (_Integer_2 * (tmp25 + tmp51 + tmp52) - _Integer_4 * tmp33 + _Rational_1_4 * (-_Integer_4 * tmp1 * tmp13 - _Integer_4 * tmp1 * tmp19 - tmp56) +
        tmp50 + tmp57)) +
      ((_Integer_2 * (tmp34 + tmp48 + tmp49) - _Integer_4 * tmp24 + _Rational_1_4 * (-_Integer_4 * tmp13 * tmp2 - _Integer_4 * tmp19 * tmp2 - tmp56) +
        tmp53 + tmp57) *
       (_Integer_2 * (tmp34 + tmp48 + tmp49) - _Integer_4 * tmp24 + _Rational_1_4 * (-_Integer_4 * tmp13 * tmp2 - _Integer_4 * tmp19 * tmp2 - tmp56) +
        tmp53 + tmp57)) +
      ((_Integer_2 * (tmp25 + tmp34) - _Rational_1_3 * tmp37 +
        _Rational_1_4 * (_Integer_3 * _Rational_4_3 * tmp38 + _Integer_3 * _Rational_4_3 * tmp40) +
        _Rational_4_3 * (_Integer_3 * tmp42 + _Integer_3 * tmp44) - _Rational_7_3 * tmp47 + tmp50 + tmp53) *
       (_Integer_2 * (tmp25 + tmp34) - _Rational_1_3 * tmp37 +
        _Rational_1_4 * (_Integer_3 * _Rational_4_3 * tmp38 + _Integer_3 * _Rational_4_3 * tmp40) +
        _Rational_4_3 * (_Integer_3 * tmp42 + _Integer_3 * tmp44) - _Rational_7_3 * tmp47 + tmp50 + tmp53));
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

#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        /*
         *  Original SymPy expressions:
         *  "[const REAL dsmin0 = params->dxx0]"
         *  "[const REAL dsmin1 = params->dxx1]"
         *  "[const REAL dsmin2 = params->dxx2]"
         */
        const REAL dsmin0 = params->dxx0;
        const REAL dsmin1 = params->dxx1;
        const REAL dsmin2 = params->dxx2;

        // Set local wavespeed
        in_gfs[IDX4(VARIABLE_WAVESPEEDGF, i0, i1, i2)] = MINIMUM_GLOBAL_WAVESPEED * NRPYMIN(dsmin0, NRPYMIN(dsmin1, dsmin2)) / dt;

      } // END LOOP: for i0 over [NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS)
    } // END LOOP: for i1 over [NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS)
  } // END LOOP: for i2 over [NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
} // END FUNCTION: variable_wavespeed_gfs_all_points_host

/**
 * Call functions that set up all AUXEVOL gridfunctions.
 */
void auxevol_gfs_set_to_constant__rfm__Cartesian(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                                                 MoL_gridfunctions_struct *restrict gridfuncs) {
#include "set_CodeParameters.h"
  REAL *auxevol_gfs = gridfuncs->auxevol_gfs;
  REAL *x0 = xx[0];
  REAL *x1 = xx[1];
  REAL *x2 = xx[2];

  // Set up variable wavespeed
  variable_wavespeed_gfs_all_points_host(params, x0, x1, x2, auxevol_gfs, dt, MINIMUM_GLOBAL_WAVESPEED);

  // Set up all other AUXEVOL gridfunctions
  auxevol_gfs_all_points_host(commondata, params, x0, x1, x2, auxevol_gfs);

} // END FUNCTION: auxevol_gfs_set_to_constant__rfm__Cartesian
