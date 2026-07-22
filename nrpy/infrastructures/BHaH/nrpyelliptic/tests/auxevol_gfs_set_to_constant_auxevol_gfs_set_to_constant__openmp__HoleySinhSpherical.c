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

  const REAL tmp0 = (1.0 / (SINHW));
  const REAL tmp8 = cos(xx2);
  const REAL tmp10 = ((AMPL) * (AMPL));
  const REAL tmp13 = sin(xx1);
  const REAL tmp17 = sin(xx2);
  const REAL tmp1 = exp(tmp0) - exp(-tmp0);
  const REAL tmp9 = ((tmp8) * (tmp8));
  const REAL tmp14 = ((tmp13) * (tmp13));
  const REAL tmp2 = (1.0 / (tmp1));
  const REAL tmp4 = exp(tmp0 * xx0) - exp(-tmp0 * xx0);
  const REAL tmp11 = (1.0 / ((tmp1) * (tmp1)));
  const REAL tmp5 = AMPL * tmp2 * tmp4;
  const REAL tmp12 = ((tmp4) * (tmp4));
  const REAL tmp6 = tmp5 * cos(xx1);
  const REAL tmp15 = tmp10 * tmp11 * tmp12 * tmp14;
  const REAL tmp7 = tmp6 + zPunc;
  const REAL tmp16 = tmp15 * tmp9;
  const REAL tmp18 = tmp15 * ((tmp17) * (tmp17));
  const REAL tmp21 = tmp6 - zPunc;
  const REAL tmp26 = tmp13 * tmp5 * tmp8;
  const REAL tmp31 = tmp13 * tmp17 * tmp5;
  const REAL tmp46 = 6 * tmp15 * tmp17 * tmp8;
  const REAL tmp27 = S1_x * tmp7 - S1_z * tmp26;
  const REAL tmp40 = S0_x * tmp21 - S0_z * tmp26;
  const REAL tmp57 = AMPL * S1_y * tmp13 * tmp2 * tmp4 * tmp8 - S1_x * tmp31;
  const REAL tmp58 = AMPL * S0_y * tmp13 * tmp2 * tmp4 * tmp8 - S0_x * tmp31;
  const REAL tmp61 = -3 * tmp6 - 3 * zPunc;
  const REAL tmp62 = -3 * tmp6 + 3 * zPunc;
  const REAL tmp20 = tmp16 + tmp18 + ((tmp7) * (tmp7));
  const REAL tmp22 = tmp16 + tmp18 + ((tmp21) * (tmp21));
  const REAL tmp33 = -S0_y * tmp21 + S0_z * tmp31;
  const REAL tmp38 = -S1_y * tmp7 + S1_z * tmp31;
  const REAL tmp45 = -1.0 / 4.0 * P1_x * tmp26 - 1.0 / 4.0 * P1_y * tmp31 - 1.0 / 4.0 * P1_z * tmp7;
  const REAL tmp50 = -1.0 / 4.0 * P0_x * tmp26 - 1.0 / 4.0 * P0_y * tmp31 - 1.0 / 4.0 * P0_z * tmp21;
  const REAL tmp23 = pow(tmp20, -3.0 / 2.0);
  const REAL tmp24 = pow(tmp22, -3.0 / 2.0);
  const REAL tmp28 = pow(tmp20, -5.0 / 2.0);
  const REAL tmp34 = pow(tmp22, -5.0 / 2.0);
  const REAL tmp52 = -tmp21 * tmp24;
  const REAL tmp54 = -tmp23 * tmp7;
  const REAL tmp65 = P1_x * tmp23 * tmp26;
  const REAL tmp66 = P0_x * tmp24 * tmp26;
  const REAL tmp68 = 4 * tmp28;
  const REAL tmp73 = 4 * tmp34;
  const REAL tmp83 = P0_z * tmp21 * tmp24;
  const REAL tmp84 = P1_z * tmp23 * tmp7;
  const REAL tmp89 = 2 * tmp28;
  const REAL tmp90 = 2 * tmp34;
  const REAL tmp30 = 3 * tmp26 * tmp28;
  const REAL tmp36 = 3 * tmp31 * tmp34;
  const REAL tmp39 = 3 * tmp28 * tmp31;
  const REAL tmp41 = 3 * tmp26 * tmp34;
  const REAL tmp70 = tmp28 * tmp57 * tmp61;
  const REAL tmp75 = tmp34 * tmp58 * tmp62;
  const REAL tmp79 = tmp28 * tmp45 * tmp61 * tmp7;
  const REAL tmp81 = tmp21 * tmp34 * tmp50 * tmp62;
  const REAL tmp93 =
      -P0_y * tmp24 * tmp31 - P1_y * tmp23 * tmp31 + tmp18 * tmp45 * tmp89 + tmp18 * tmp50 * tmp90 + tmp27 * tmp31 * tmp89 + tmp31 * tmp40 * tmp90;
  *psi_background = (1.0 / 2.0) * bare_mass_0 / sqrt(tmp22) + (1.0 / 2.0) * bare_mass_1 / sqrt(tmp20) + 1;
  *ADD_times_AUU =
      2 * (((3.0 / 2.0) * AMPL * P0_x * tmp13 * tmp17 * tmp2 * tmp24 * tmp4 + (3.0 / 2.0) * AMPL * P0_y * tmp13 * tmp2 * tmp24 * tmp4 * tmp8 +
            (3.0 / 2.0) * AMPL * P1_x * tmp13 * tmp17 * tmp2 * tmp23 * tmp4 + (3.0 / 2.0) * AMPL * P1_y * tmp13 * tmp2 * tmp23 * tmp4 * tmp8 -
            tmp27 * tmp30 - tmp28 * tmp45 * tmp46 - tmp33 * tmp36 - tmp34 * tmp46 * tmp50 - tmp38 * tmp39 - tmp40 * tmp41) *
           ((3.0 / 2.0) * AMPL * P0_x * tmp13 * tmp17 * tmp2 * tmp24 * tmp4 + (3.0 / 2.0) * AMPL * P0_y * tmp13 * tmp2 * tmp24 * tmp4 * tmp8 +
            (3.0 / 2.0) * AMPL * P1_x * tmp13 * tmp17 * tmp2 * tmp23 * tmp4 + (3.0 / 2.0) * AMPL * P1_y * tmp13 * tmp2 * tmp23 * tmp4 * tmp8 -
            tmp27 * tmp30 - tmp28 * tmp45 * tmp46 - tmp33 * tmp36 - tmp34 * tmp46 * tmp50 - tmp38 * tmp39 - tmp40 * tmp41)) +
      2 * (((3.0 / 2.0) * AMPL * P0_z * tmp13 * tmp17 * tmp2 * tmp24 * tmp4 + (3.0 / 2.0) * AMPL * P1_z * tmp13 * tmp17 * tmp2 * tmp23 * tmp4 +
            AMPL * tmp13 * tmp17 * tmp2 * tmp28 * tmp4 * tmp45 * tmp61 + AMPL * tmp13 * tmp17 * tmp2 * tmp34 * tmp4 * tmp50 * tmp62 -
            1.0 / 4.0 * P0_y * tmp21 * tmp24 - 7.0 / 4.0 * P0_y * tmp52 - 1.0 / 4.0 * P1_y * tmp23 * tmp7 - 7.0 / 4.0 * P1_y * tmp54 -
            tmp21 * tmp36 * tmp50 + tmp27 * tmp28 * tmp61 + tmp34 * tmp40 * tmp62 - tmp36 * tmp58 - tmp39 * tmp45 * tmp7 - tmp39 * tmp57) *
           ((3.0 / 2.0) * AMPL * P0_z * tmp13 * tmp17 * tmp2 * tmp24 * tmp4 + (3.0 / 2.0) * AMPL * P1_z * tmp13 * tmp17 * tmp2 * tmp23 * tmp4 +
            AMPL * tmp13 * tmp17 * tmp2 * tmp28 * tmp4 * tmp45 * tmp61 + AMPL * tmp13 * tmp17 * tmp2 * tmp34 * tmp4 * tmp50 * tmp62 -
            1.0 / 4.0 * P0_y * tmp21 * tmp24 - 7.0 / 4.0 * P0_y * tmp52 - 1.0 / 4.0 * P1_y * tmp23 * tmp7 - 7.0 / 4.0 * P1_y * tmp54 -
            tmp21 * tmp36 * tmp50 + tmp27 * tmp28 * tmp61 + tmp34 * tmp40 * tmp62 - tmp36 * tmp58 - tmp39 * tmp45 * tmp7 - tmp39 * tmp57)) +
      2 * (((3.0 / 2.0) * AMPL * P0_z * tmp13 * tmp2 * tmp24 * tmp4 * tmp8 + (3.0 / 2.0) * AMPL * P1_z * tmp13 * tmp2 * tmp23 * tmp4 * tmp8 +
            AMPL * tmp13 * tmp2 * tmp28 * tmp4 * tmp45 * tmp61 * tmp8 + AMPL * tmp13 * tmp2 * tmp34 * tmp4 * tmp50 * tmp62 * tmp8 -
            1.0 / 4.0 * P0_x * tmp21 * tmp24 - 7.0 / 4.0 * P0_x * tmp52 - 1.0 / 4.0 * P1_x * tmp23 * tmp7 - 7.0 / 4.0 * P1_x * tmp54 -
            tmp21 * tmp41 * tmp50 + tmp28 * tmp38 * tmp61 - tmp30 * tmp45 * tmp7 - tmp30 * tmp57 + tmp33 * tmp34 * tmp62 - tmp41 * tmp58) *
           ((3.0 / 2.0) * AMPL * P0_z * tmp13 * tmp2 * tmp24 * tmp4 * tmp8 + (3.0 / 2.0) * AMPL * P1_z * tmp13 * tmp2 * tmp23 * tmp4 * tmp8 +
            AMPL * tmp13 * tmp2 * tmp28 * tmp4 * tmp45 * tmp61 * tmp8 + AMPL * tmp13 * tmp2 * tmp34 * tmp4 * tmp50 * tmp62 * tmp8 -
            1.0 / 4.0 * P0_x * tmp21 * tmp24 - 7.0 / 4.0 * P0_x * tmp52 - 1.0 / 4.0 * P1_x * tmp23 * tmp7 - 7.0 / 4.0 * P1_x * tmp54 -
            tmp21 * tmp41 * tmp50 + tmp28 * tmp38 * tmp61 - tmp30 * tmp45 * tmp7 - tmp30 * tmp57 + tmp33 * tmp34 * tmp62 - tmp41 * tmp58)) +
      ((-7.0 / 3.0 * P0_z * tmp52 - 7.0 / 3.0 * P1_z * tmp54 + tmp16 * tmp45 * tmp89 + tmp16 * tmp50 * tmp90 + tmp26 * tmp33 * tmp90 +
        tmp26 * tmp38 * tmp89 - tmp65 - tmp66 + (4.0 / 3.0) * tmp70 + (4.0 / 3.0) * tmp75 + (4.0 / 3.0) * tmp79 + (4.0 / 3.0) * tmp81 -
        1.0 / 3.0 * tmp83 - 1.0 / 3.0 * tmp84 + tmp93) *
       (-7.0 / 3.0 * P0_z * tmp52 - 7.0 / 3.0 * P1_z * tmp54 + tmp16 * tmp45 * tmp89 + tmp16 * tmp50 * tmp90 + tmp26 * tmp33 * tmp90 +
        tmp26 * tmp38 * tmp89 - tmp65 - tmp66 + (4.0 / 3.0) * tmp70 + (4.0 / 3.0) * tmp75 + (4.0 / 3.0) * tmp79 + (4.0 / 3.0) * tmp81 -
        1.0 / 3.0 * tmp83 - 1.0 / 3.0 * tmp84 + tmp93)) +
      (((7.0 / 6.0) * P0_z * tmp52 + (7.0 / 6.0) * P1_z * tmp54 - tmp16 * tmp45 * tmp68 - tmp16 * tmp50 * tmp73 - tmp26 * tmp33 * tmp73 -
        tmp26 * tmp38 * tmp68 + 2 * tmp65 + 2 * tmp66 - 2.0 / 3.0 * tmp70 - 2.0 / 3.0 * tmp75 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp81 +
        (1.0 / 6.0) * tmp83 + (1.0 / 6.0) * tmp84 + tmp93) *
       ((7.0 / 6.0) * P0_z * tmp52 + (7.0 / 6.0) * P1_z * tmp54 - tmp16 * tmp45 * tmp68 - tmp16 * tmp50 * tmp73 - tmp26 * tmp33 * tmp73 -
        tmp26 * tmp38 * tmp68 + 2 * tmp65 + 2 * tmp66 - 2.0 / 3.0 * tmp70 - 2.0 / 3.0 * tmp75 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp81 +
        (1.0 / 6.0) * tmp83 + (1.0 / 6.0) * tmp84 + tmp93)) +
      ((2 * AMPL * P0_y * tmp13 * tmp17 * tmp2 * tmp24 * tmp4 + 2 * AMPL * P1_y * tmp13 * tmp17 * tmp2 * tmp23 * tmp4 +
        2 * AMPL * tmp13 * tmp2 * tmp28 * tmp38 * tmp4 * tmp8 + 2 * AMPL * tmp13 * tmp2 * tmp33 * tmp34 * tmp4 * tmp8 - P0_z * tmp21 * tmp24 -
        P1_z * tmp23 * tmp7 + 2 * tmp10 * tmp11 * tmp12 * tmp14 * tmp28 * tmp45 * tmp9 + 2 * tmp10 * tmp11 * tmp12 * tmp14 * tmp34 * tmp50 * tmp9 -
        tmp18 * tmp45 * tmp68 - tmp18 * tmp50 * tmp73 - tmp27 * tmp31 * tmp68 - tmp31 * tmp40 * tmp73 - tmp65 - tmp66 - 2.0 / 3.0 * tmp70 -
        2.0 / 3.0 * tmp75 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp81) *
       (2 * AMPL * P0_y * tmp13 * tmp17 * tmp2 * tmp24 * tmp4 + 2 * AMPL * P1_y * tmp13 * tmp17 * tmp2 * tmp23 * tmp4 +
        2 * AMPL * tmp13 * tmp2 * tmp28 * tmp38 * tmp4 * tmp8 + 2 * AMPL * tmp13 * tmp2 * tmp33 * tmp34 * tmp4 * tmp8 - P0_z * tmp21 * tmp24 -
        P1_z * tmp23 * tmp7 + 2 * tmp10 * tmp11 * tmp12 * tmp14 * tmp28 * tmp45 * tmp9 + 2 * tmp10 * tmp11 * tmp12 * tmp14 * tmp34 * tmp50 * tmp9 -
        tmp18 * tmp45 * tmp68 - tmp18 * tmp50 * tmp73 - tmp27 * tmp31 * tmp68 - tmp31 * tmp40 * tmp73 - tmp65 - tmp66 - 2.0 / 3.0 * tmp70 -
        2.0 / 3.0 * tmp75 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp81));
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
