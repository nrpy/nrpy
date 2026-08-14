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
  const REAL tmp12 = 3 * xx0;
  const REAL tmp17 = 3 * xx1;
  const REAL tmp28 = 6 * xx0 * xx1;
  const REAL tmp40 = -S1_x * xx1 + S1_y * xx0;
  const REAL tmp42 = -S0_x * xx1 + S0_y * xx0;
  const REAL tmp10 = S1_x * tmp0 - S1_z * xx0;
  const REAL tmp21 = S0_x * tmp5 - S0_z * xx0;
  const REAL tmp46 = -3 * xx2 - 3 * zPunc;
  const REAL tmp48 = -3 * xx2 + 3 * zPunc;
  const REAL tmp4 = ((tmp0) * (tmp0)) + tmp1 + tmp2;
  const REAL tmp6 = tmp1 + tmp2 + ((tmp5) * (tmp5));
  const REAL tmp15 = -S0_y * tmp5 + S0_z * xx1;
  const REAL tmp19 = -S1_y * tmp0 + S1_z * xx1;
  const REAL tmp26 = -1.0 / 4.0 * P1_x * xx0 - 1.0 / 4.0 * P1_y * xx1 - 1.0 / 4.0 * P1_z * tmp0;
  const REAL tmp32 = -1.0 / 4.0 * P0_x * xx0 - 1.0 / 4.0 * P0_y * xx1 - 1.0 / 4.0 * P0_z * tmp5;
  const REAL tmp7 = pow(tmp4, -3.0 / 2.0);
  const REAL tmp8 = pow(tmp6, -3.0 / 2.0);
  const REAL tmp9 = pow(tmp4, -5.0 / 2.0);
  const REAL tmp13 = pow(tmp6, -5.0 / 2.0);
  const REAL tmp11 = tmp10 * tmp9;
  const REAL tmp16 = tmp13 * tmp15;
  const REAL tmp20 = tmp19 * tmp9;
  const REAL tmp22 = tmp13 * tmp21;
  const REAL tmp27 = tmp26 * tmp9;
  const REAL tmp33 = tmp13 * tmp32;
  const REAL tmp35 = -tmp5 * tmp8;
  const REAL tmp38 = -tmp0 * tmp7;
  const REAL tmp50 = P0_x * tmp8 * xx0;
  const REAL tmp51 = P1_x * tmp7 * xx0;
  const REAL tmp56 = P0_z * tmp5 * tmp8;
  const REAL tmp57 = P1_z * tmp0 * tmp7;
  const REAL tmp59 = tmp40 * tmp46 * tmp9;
  const REAL tmp61 = tmp13 * tmp42 * tmp48;
  const REAL tmp62 = tmp0 * tmp26 * tmp46 * tmp9;
  const REAL tmp63 = tmp13 * tmp32 * tmp48 * tmp5;
  const REAL tmp66 = -P0_y * tmp8 * xx1 - P1_y * tmp7 * xx1 + 2 * tmp11 * xx1 + 2 * tmp2 * tmp27 + 2 * tmp2 * tmp33 + 2 * tmp22 * xx1;
  *psi_background = (1.0 / 2.0) * bare_mass_0 / sqrt(tmp6) + (1.0 / 2.0) * bare_mass_1 / sqrt(tmp4) + 1;
  *ADD_times_AUU =
      2 * (((3.0 / 2.0) * P0_x * tmp8 * xx1 + (3.0 / 2.0) * P0_y * tmp8 * xx0 + (3.0 / 2.0) * P1_x * tmp7 * xx1 + (3.0 / 2.0) * P1_y * tmp7 * xx0 -
            tmp11 * tmp12 - tmp12 * tmp22 - tmp16 * tmp17 - tmp17 * tmp20 - tmp27 * tmp28 - tmp28 * tmp33) *
           ((3.0 / 2.0) * P0_x * tmp8 * xx1 + (3.0 / 2.0) * P0_y * tmp8 * xx0 + (3.0 / 2.0) * P1_x * tmp7 * xx1 + (3.0 / 2.0) * P1_y * tmp7 * xx0 -
            tmp11 * tmp12 - tmp12 * tmp22 - tmp16 * tmp17 - tmp17 * tmp20 - tmp27 * tmp28 - tmp28 * tmp33)) +
      2 * ((-7.0 / 4.0 * P0_x * tmp35 - 1.0 / 4.0 * P0_x * tmp5 * tmp8 + (3.0 / 2.0) * P0_z * tmp8 * xx0 - 1.0 / 4.0 * P1_x * tmp0 * tmp7 -
            7.0 / 4.0 * P1_x * tmp38 + (3.0 / 2.0) * P1_z * tmp7 * xx0 - tmp0 * tmp12 * tmp27 - tmp12 * tmp13 * tmp42 - tmp12 * tmp33 * tmp5 -
            tmp12 * tmp40 * tmp9 + tmp13 * tmp15 * tmp48 + tmp13 * tmp32 * tmp48 * xx0 + tmp19 * tmp46 * tmp9 + tmp26 * tmp46 * tmp9 * xx0) *
           (-7.0 / 4.0 * P0_x * tmp35 - 1.0 / 4.0 * P0_x * tmp5 * tmp8 + (3.0 / 2.0) * P0_z * tmp8 * xx0 - 1.0 / 4.0 * P1_x * tmp0 * tmp7 -
            7.0 / 4.0 * P1_x * tmp38 + (3.0 / 2.0) * P1_z * tmp7 * xx0 - tmp0 * tmp12 * tmp27 - tmp12 * tmp13 * tmp42 - tmp12 * tmp33 * tmp5 -
            tmp12 * tmp40 * tmp9 + tmp13 * tmp15 * tmp48 + tmp13 * tmp32 * tmp48 * xx0 + tmp19 * tmp46 * tmp9 + tmp26 * tmp46 * tmp9 * xx0)) +
      2 * ((-7.0 / 4.0 * P0_y * tmp35 - 1.0 / 4.0 * P0_y * tmp5 * tmp8 + (3.0 / 2.0) * P0_z * tmp8 * xx1 - 1.0 / 4.0 * P1_y * tmp0 * tmp7 -
            7.0 / 4.0 * P1_y * tmp38 + (3.0 / 2.0) * P1_z * tmp7 * xx1 - tmp0 * tmp17 * tmp27 + tmp10 * tmp46 * tmp9 - tmp13 * tmp17 * tmp42 +
            tmp13 * tmp21 * tmp48 + tmp13 * tmp32 * tmp48 * xx1 - tmp17 * tmp33 * tmp5 - tmp17 * tmp40 * tmp9 + tmp26 * tmp46 * tmp9 * xx1) *
           (-7.0 / 4.0 * P0_y * tmp35 - 1.0 / 4.0 * P0_y * tmp5 * tmp8 + (3.0 / 2.0) * P0_z * tmp8 * xx1 - 1.0 / 4.0 * P1_y * tmp0 * tmp7 -
            7.0 / 4.0 * P1_y * tmp38 + (3.0 / 2.0) * P1_z * tmp7 * xx1 - tmp0 * tmp17 * tmp27 + tmp10 * tmp46 * tmp9 - tmp13 * tmp17 * tmp42 +
            tmp13 * tmp21 * tmp48 + tmp13 * tmp32 * tmp48 * xx1 - tmp17 * tmp33 * tmp5 - tmp17 * tmp40 * tmp9 + tmp26 * tmp46 * tmp9 * xx1)) +
      ((-7.0 / 3.0 * P0_z * tmp35 - 7.0 / 3.0 * P1_z * tmp38 + 2 * tmp1 * tmp27 + 2 * tmp1 * tmp33 + 2 * tmp16 * xx0 + 2 * tmp20 * xx0 - tmp50 -
        tmp51 - 1.0 / 3.0 * tmp56 - 1.0 / 3.0 * tmp57 + (4.0 / 3.0) * tmp59 + (4.0 / 3.0) * tmp61 + (4.0 / 3.0) * tmp62 + (4.0 / 3.0) * tmp63 +
        tmp66) *
       (-7.0 / 3.0 * P0_z * tmp35 - 7.0 / 3.0 * P1_z * tmp38 + 2 * tmp1 * tmp27 + 2 * tmp1 * tmp33 + 2 * tmp16 * xx0 + 2 * tmp20 * xx0 - tmp50 -
        tmp51 - 1.0 / 3.0 * tmp56 - 1.0 / 3.0 * tmp57 + (4.0 / 3.0) * tmp59 + (4.0 / 3.0) * tmp61 + (4.0 / 3.0) * tmp62 + (4.0 / 3.0) * tmp63 +
        tmp66)) +
      (((7.0 / 6.0) * P0_z * tmp35 + (7.0 / 6.0) * P1_z * tmp38 - 4 * tmp1 * tmp27 - 4 * tmp1 * tmp33 - 4 * tmp16 * xx0 - 4 * tmp20 * xx0 +
        2 * tmp50 + 2 * tmp51 + (1.0 / 6.0) * tmp56 + (1.0 / 6.0) * tmp57 - 2.0 / 3.0 * tmp59 - 2.0 / 3.0 * tmp61 - 2.0 / 3.0 * tmp62 -
        2.0 / 3.0 * tmp63 + tmp66) *
       ((7.0 / 6.0) * P0_z * tmp35 + (7.0 / 6.0) * P1_z * tmp38 - 4 * tmp1 * tmp27 - 4 * tmp1 * tmp33 - 4 * tmp16 * xx0 - 4 * tmp20 * xx0 +
        2 * tmp50 + 2 * tmp51 + (1.0 / 6.0) * tmp56 + (1.0 / 6.0) * tmp57 - 2.0 / 3.0 * tmp59 - 2.0 / 3.0 * tmp61 - 2.0 / 3.0 * tmp62 -
        2.0 / 3.0 * tmp63 + tmp66)) +
      ((2 * P0_y * tmp8 * xx1 - P0_z * tmp5 * tmp8 + 2 * P1_y * tmp7 * xx1 - P1_z * tmp0 * tmp7 + 2 * tmp1 * tmp13 * tmp32 + 2 * tmp1 * tmp26 * tmp9 -
        4 * tmp11 * xx1 + 2 * tmp13 * tmp15 * xx0 + 2 * tmp19 * tmp9 * xx0 - 4 * tmp2 * tmp27 - 4 * tmp2 * tmp33 - 4 * tmp22 * xx1 - tmp50 - tmp51 -
        2.0 / 3.0 * tmp59 - 2.0 / 3.0 * tmp61 - 2.0 / 3.0 * tmp62 - 2.0 / 3.0 * tmp63) *
       (2 * P0_y * tmp8 * xx1 - P0_z * tmp5 * tmp8 + 2 * P1_y * tmp7 * xx1 - P1_z * tmp0 * tmp7 + 2 * tmp1 * tmp13 * tmp32 + 2 * tmp1 * tmp26 * tmp9 -
        4 * tmp11 * xx1 + 2 * tmp13 * tmp15 * xx0 + 2 * tmp19 * tmp9 * xx0 - 4 * tmp2 * tmp27 - 4 * tmp2 * tmp33 - 4 * tmp22 * xx1 - tmp50 - tmp51 -
        2.0 / 3.0 * tmp59 - 2.0 / 3.0 * tmp61 - 2.0 / 3.0 * tmp62 - 2.0 / 3.0 * tmp63));
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
