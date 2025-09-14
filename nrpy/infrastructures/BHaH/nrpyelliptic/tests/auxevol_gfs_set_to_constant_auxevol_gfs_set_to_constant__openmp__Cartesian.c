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

  const REAL tmp0 = xx2 + zPunc;
  const REAL tmp1 = ((xx0) * (xx0));
  const REAL tmp2 = ((xx1) * (xx1));
  const REAL tmp5 = xx2 - zPunc;
  const REAL tmp30 = 6 * xx0;
  const REAL tmp40 = -S1_x * xx1 + S1_y * xx0;
  const REAL tmp46 = -S0_x * xx1 + S0_y * xx0;
  const REAL tmp52 = 3 * xx0;
  const REAL tmp53 = 3 * xx1;
  const REAL tmp9 = S1_x * tmp0 - S1_z * xx0;
  const REAL tmp12 = S0_x * tmp5 - S0_z * xx0;
  const REAL tmp43 = -3 * xx2 - 3 * zPunc;
  const REAL tmp47 = -3 * xx2 + 3 * zPunc;
  const REAL tmp4 = ((tmp0) * (tmp0)) + tmp1 + tmp2;
  const REAL tmp6 = tmp1 + tmp2 + ((tmp5) * (tmp5));
  const REAL tmp17 = -1.0 / 4.0 * P1_x * xx0 - 1.0 / 4.0 * P1_y * xx1 - 1.0 / 4.0 * P1_z * tmp0;
  const REAL tmp21 = -1.0 / 4.0 * P0_x * xx0 - 1.0 / 4.0 * P0_y * xx1 - 1.0 / 4.0 * P0_z * tmp5;
  const REAL tmp32 = -S1_y * tmp0 + S1_z * xx1;
  const REAL tmp35 = -S0_y * tmp5 + S0_z * xx1;
  const REAL tmp8 = pow(tmp4, -5.0 / 2.0);
  const REAL tmp11 = pow(tmp6, -5.0 / 2.0);
  const REAL tmp23 = pow(tmp6, -3.0 / 2.0);
  const REAL tmp25 = pow(tmp4, -3.0 / 2.0);
  const REAL tmp18 = tmp17 * tmp8;
  const REAL tmp22 = tmp11 * tmp21;
  const REAL tmp24 = P0_y * tmp23 * xx1;
  const REAL tmp26 = P1_y * tmp25 * xx1;
  const REAL tmp29 = 2 * tmp17 * tmp25 + 2 * tmp21 * tmp23;
  const REAL tmp38 = P0_x * tmp23 * xx0;
  const REAL tmp39 = P1_x * tmp25 * xx0;
  const REAL tmp44 = tmp43 * tmp8;
  const REAL tmp48 = tmp11 * tmp47;
  const REAL tmp50 = -tmp23 * tmp5;
  const REAL tmp51 = -tmp0 * tmp25;
  const REAL tmp58 = (1.0 / 4.0) * tmp23 * tmp5;
  const REAL tmp59 = (1.0 / 4.0) * tmp0 * tmp25;
  const REAL tmp55 = P0_z * tmp23 * tmp5 + P1_z * tmp0 * tmp25 + (2.0 / 3.0) * tmp0 * tmp17 * tmp44 - 2 * tmp1 * tmp11 * tmp21 -
                     2 * tmp1 * tmp17 * tmp8 - 2 * tmp11 * tmp12 * xx1 - 2 * tmp11 * tmp2 * tmp21 - 2 * tmp11 * tmp35 * xx0 -
                     2 * tmp17 * tmp2 * tmp8 + (2.0 / 3.0) * tmp21 * tmp48 * tmp5 + tmp24 + tmp26 + tmp29 - 2 * tmp32 * tmp8 * xx0 + tmp38 + tmp39 +
                     (2.0 / 3.0) * tmp40 * tmp44 + (2.0 / 3.0) * tmp46 * tmp48 - 2 * tmp8 * tmp9 * xx1;
  *psi_background = (1.0 / 2.0) * bare_mass_0 / sqrt(tmp6) + (1.0 / 2.0) * bare_mass_1 / sqrt(tmp4) + 1;
  *ADD_times_AUU =
      ((-6 * tmp1 * tmp18 - 6 * tmp1 * tmp22 - tmp11 * tmp30 * tmp35 + tmp29 - tmp30 * tmp32 * tmp8 + 3 * tmp38 + 3 * tmp39) *
       (-6 * tmp1 * tmp18 - 6 * tmp1 * tmp22 - tmp11 * tmp30 * tmp35 + tmp29 - tmp30 * tmp32 * tmp8 + 3 * tmp38 + 3 * tmp39)) +
      ((-6 * tmp11 * tmp12 * xx1 - 6 * tmp18 * tmp2 - 6 * tmp2 * tmp22 + 3 * tmp24 + 3 * tmp26 + tmp29 - 6 * tmp8 * tmp9 * xx1) *
       (-6 * tmp11 * tmp12 * xx1 - 6 * tmp18 * tmp2 - 6 * tmp2 * tmp22 + 3 * tmp24 + 3 * tmp26 + tmp29 - 6 * tmp8 * tmp9 * xx1)) +
      ((-1.0 / 2.0 * P0_z * tmp23 * tmp5 - 7.0 / 2.0 * P0_z * tmp50 - 1.0 / 2.0 * P1_z * tmp0 * tmp25 - 7.0 / 2.0 * P1_z * tmp51 +
        2 * tmp0 * tmp17 * tmp44 + 2 * tmp21 * tmp48 * tmp5 + tmp29 + 2 * tmp40 * tmp44 + 2 * tmp46 * tmp48) *
       (-1.0 / 2.0 * P0_z * tmp23 * tmp5 - 7.0 / 2.0 * P0_z * tmp50 - 1.0 / 2.0 * P1_z * tmp0 * tmp25 - 7.0 / 2.0 * P1_z * tmp51 +
        2 * tmp0 * tmp17 * tmp44 + 2 * tmp21 * tmp48 * tmp5 + tmp29 + 2 * tmp40 * tmp44 + 2 * tmp46 * tmp48)) +
      2 * (((3.0 / 2.0) * P0_x * tmp23 * xx1 + (3.0 / 2.0) * P0_y * tmp23 * xx0 + (3.0 / 2.0) * P1_x * tmp25 * xx1 +
            (3.0 / 2.0) * P1_y * tmp25 * xx0 - tmp11 * tmp12 * tmp52 - tmp11 * tmp35 * tmp53 - tmp18 * tmp30 * xx1 - tmp22 * tmp30 * xx1 -
            tmp32 * tmp53 * tmp8 - tmp52 * tmp8 * tmp9 - tmp55) *
           ((3.0 / 2.0) * P0_x * tmp23 * xx1 + (3.0 / 2.0) * P0_y * tmp23 * xx0 + (3.0 / 2.0) * P1_x * tmp25 * xx1 +
            (3.0 / 2.0) * P1_y * tmp25 * xx0 - tmp11 * tmp12 * tmp52 - tmp11 * tmp35 * tmp53 - tmp18 * tmp30 * xx1 - tmp22 * tmp30 * xx1 -
            tmp32 * tmp53 * tmp8 - tmp52 * tmp8 * tmp9 - tmp55)) +
      2 * ((-7.0 / 4.0 * P0_x * tmp50 - P0_x * tmp58 + (3.0 / 2.0) * P0_z * tmp23 * xx0 - 7.0 / 4.0 * P1_x * tmp51 - P1_x * tmp59 +
            (3.0 / 2.0) * P1_z * tmp25 * xx0 - tmp0 * tmp18 * tmp52 + tmp11 * tmp21 * tmp47 * xx0 + tmp11 * tmp35 * tmp47 - tmp11 * tmp46 * tmp52 +
            tmp17 * tmp43 * tmp8 * xx0 - tmp22 * tmp5 * tmp52 + tmp32 * tmp43 * tmp8 - tmp40 * tmp52 * tmp8 - tmp55) *
           (-7.0 / 4.0 * P0_x * tmp50 - P0_x * tmp58 + (3.0 / 2.0) * P0_z * tmp23 * xx0 - 7.0 / 4.0 * P1_x * tmp51 - P1_x * tmp59 +
            (3.0 / 2.0) * P1_z * tmp25 * xx0 - tmp0 * tmp18 * tmp52 + tmp11 * tmp21 * tmp47 * xx0 + tmp11 * tmp35 * tmp47 - tmp11 * tmp46 * tmp52 +
            tmp17 * tmp43 * tmp8 * xx0 - tmp22 * tmp5 * tmp52 + tmp32 * tmp43 * tmp8 - tmp40 * tmp52 * tmp8 - tmp55)) +
      2 * ((-7.0 / 4.0 * P0_y * tmp50 - P0_y * tmp58 + (3.0 / 2.0) * P0_z * tmp23 * xx1 - 7.0 / 4.0 * P1_y * tmp51 - P1_y * tmp59 +
            (3.0 / 2.0) * P1_z * tmp25 * xx1 - tmp0 * tmp18 * tmp53 + tmp11 * tmp12 * tmp47 + tmp11 * tmp21 * tmp47 * xx1 - tmp11 * tmp46 * tmp53 +
            tmp17 * tmp43 * tmp8 * xx1 - tmp22 * tmp5 * tmp53 - tmp40 * tmp53 * tmp8 + tmp43 * tmp8 * tmp9 - tmp55) *
           (-7.0 / 4.0 * P0_y * tmp50 - P0_y * tmp58 + (3.0 / 2.0) * P0_z * tmp23 * xx1 - 7.0 / 4.0 * P1_y * tmp51 - P1_y * tmp59 +
            (3.0 / 2.0) * P1_z * tmp25 * xx1 - tmp0 * tmp18 * tmp53 + tmp11 * tmp12 * tmp47 + tmp11 * tmp21 * tmp47 * xx1 - tmp11 * tmp46 * tmp53 +
            tmp17 * tmp43 * tmp8 * xx1 - tmp22 * tmp5 * tmp53 - tmp40 * tmp53 * tmp8 + tmp43 * tmp8 * tmp9 - tmp55));
} // END FUNCTION auxevol_gfs_single_point_host

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

      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
} // END FUNCTION auxevol_gfs_all_points_host

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
        in_gfs[IDX4(VARIABLE_WAVESPEEDGF, i0, i1, i2)] = MINIMUM_GLOBAL_WAVESPEED * MIN(dsmin0, MIN(dsmin1, dsmin2)) / dt;

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION variable_wavespeed_gfs_all_points_host

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

} // END FUNCTION auxevol_gfs_set_to_constant__rfm__Cartesian
