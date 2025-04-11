#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"

/**
 * Kernel: auxevol_gfs_single_point_host.
 * Kernel to compute AUXEVOL grid functions at a single point.
 */
static void auxevol_gfs_single_point_host(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xx0,
                                          const REAL xx1, const REAL xx2, REAL *restrict psi_background, REAL *restrict ADD_times_AUU) {
  // Load necessary parameters from params_struct
  const REAL AMAX = params->AMAX;
  const REAL SINHWAA = params->SINHWAA;
  const REAL bScale = params->bScale;

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

  const REAL tmp0 = (1.0 / (SINHWAA));
  const REAL tmp7 = cos(xx2);
  const REAL tmp8 = sin(xx1);
  const REAL tmp11 = sin(xx2);
  const REAL tmp1 = exp(tmp0) - exp(-tmp0);
  const REAL tmp3 = exp(tmp0 * xx0) - exp(-tmp0 * xx0);
  const REAL tmp18 = (1.0 / (tmp1));
  const REAL tmp4 = ((AMAX) * (AMAX)) * ((tmp3) * (tmp3)) / ((tmp1) * (tmp1));
  const REAL tmp21 = AMAX * tmp18 * tmp3 * tmp8;
  const REAL tmp5 = sqrt(((bScale) * (bScale)) + tmp4) * cos(xx1);
  const REAL tmp9 = tmp4 * ((tmp8) * (tmp8));
  const REAL tmp22 = tmp11 * tmp21;
  const REAL tmp25 = tmp21 * tmp7;
  const REAL tmp6 = tmp5 + zPunc;
  const REAL tmp10 = ((tmp7) * (tmp7)) * tmp9;
  const REAL tmp12 = ((tmp11) * (tmp11)) * tmp9;
  const REAL tmp15 = -tmp5 + zPunc;
  const REAL tmp60 = AMAX * S1_y * tmp18 * tmp3 * tmp7 * tmp8 - S1_x * tmp22;
  const REAL tmp65 = AMAX * S0_y * tmp18 * tmp3 * tmp7 * tmp8 - S0_x * tmp22;
  const REAL tmp70 = 3 * tmp25;
  const REAL tmp71 = 3 * tmp22;
  const REAL tmp73 = (3.0 / 2.0) * tmp22;
  const REAL tmp75 = (3.0 / 2.0) * tmp25;
  const REAL tmp78 = 6 * tmp11 * tmp7 * tmp9;
  const REAL tmp26 = S1_x * tmp6 - S1_z * tmp25;
  const REAL tmp51 = S0_y * tmp15 + S0_z * tmp22;
  const REAL tmp14 = tmp10 + tmp12 + ((tmp6) * (tmp6));
  const REAL tmp17 = tmp10 + tmp12 + ((tmp15) * (tmp15));
  const REAL tmp29 = -S0_x * tmp15 - S0_z * tmp25;
  const REAL tmp35 = -1.0 / 4.0 * P1_x * tmp25 - 1.0 / 4.0 * P1_y * tmp22 - 1.0 / 4.0 * P1_z * tmp6;
  const REAL tmp39 = -P0_z * tmp15;
  const REAL tmp49 = -S1_y * tmp6 + S1_z * tmp22;
  const REAL tmp19 = pow(tmp17, -3.0 / 2.0);
  const REAL tmp20 = pow(tmp14, -3.0 / 2.0);
  const REAL tmp24 = pow(tmp14, -5.0 / 2.0);
  const REAL tmp28 = pow(tmp17, -5.0 / 2.0);
  const REAL tmp40 = -1.0 / 4.0 * P0_x * tmp25 - 1.0 / 4.0 * P0_y * tmp22 - 1.0 / 4.0 * tmp39;
  const REAL tmp42 = 2 * tmp35;
  const REAL tmp27 = tmp24 * tmp26;
  const REAL tmp30 = tmp28 * tmp29;
  const REAL tmp36 = tmp24 * tmp35;
  const REAL tmp41 = tmp28 * tmp40;
  const REAL tmp44 = 2 * tmp40;
  const REAL tmp50 = tmp24 * tmp49;
  const REAL tmp52 = tmp28 * tmp51;
  const REAL tmp54 = P1_z * tmp20 * tmp6;
  const REAL tmp57 = -tmp20 * tmp6;
  const REAL tmp63 = tmp24 * (-3 * tmp5 - 3 * zPunc);
  const REAL tmp66 = tmp28 * (-3 * tmp5 + 3 * zPunc);
  const REAL tmp86 = -tmp15 * tmp19;
  const REAL tmp56 = P0_z * tmp15 * tmp19;
  const REAL tmp69 = -tmp15 * tmp66;
  const REAL tmp84 = (7.0 / 4.0) * tmp15 * tmp19;
  const REAL tmp95 = -tmp15 * tmp41;
  const REAL tmp46 = -tmp19 * tmp44 - tmp20 * tmp42;
  const REAL tmp83 = -P0_x * tmp19 * tmp25 - P0_y * tmp19 * tmp22 - P1_x * tmp20 * tmp25 - P1_y * tmp20 * tmp22 + (7.0 / 6.0) * P1_z * tmp57 +
                     tmp10 * tmp24 * tmp42 + tmp10 * tmp28 * tmp44 + tmp12 * tmp24 * tmp42 + tmp12 * tmp28 * tmp44 + (1.0 / 6.0) * tmp19 * tmp39 +
                     2 * tmp22 * tmp27 + 2 * tmp22 * tmp30 + 2 * tmp25 * tmp50 + 2 * tmp25 * tmp52 - 2.0 / 3.0 * tmp35 * tmp6 * tmp63 -
                     2.0 / 3.0 * tmp40 * tmp69 + tmp46 + (1.0 / 6.0) * tmp54 + (7.0 / 6.0) * tmp56 - 2.0 / 3.0 * tmp60 * tmp63 -
                     2.0 / 3.0 * tmp65 * tmp66;
  *psi_background = (1.0 / 2.0) * bare_mass_0 / sqrt(tmp17) + (1.0 / 2.0) * bare_mass_1 / sqrt(tmp14) + 1;
  *ADD_times_AUU = ((3 * AMAX * P0_x * tmp18 * tmp19 * tmp3 * tmp7 * tmp8 + 3 * AMAX * P1_x * tmp18 * tmp20 * tmp3 * tmp7 * tmp8 - 6 * tmp10 * tmp36 -
                     6 * tmp10 * tmp41 - 6 * tmp25 * tmp50 - 6 * tmp25 * tmp52 - tmp46) *
                    (3 * AMAX * P0_x * tmp18 * tmp19 * tmp3 * tmp7 * tmp8 + 3 * AMAX * P1_x * tmp18 * tmp20 * tmp3 * tmp7 * tmp8 - 6 * tmp10 * tmp36 -
                     6 * tmp10 * tmp41 - 6 * tmp25 * tmp50 - 6 * tmp25 * tmp52 - tmp46)) +
                   ((3 * AMAX * P0_y * tmp11 * tmp18 * tmp19 * tmp3 * tmp8 + 3 * AMAX * P1_y * tmp11 * tmp18 * tmp20 * tmp3 * tmp8 -
                     6 * tmp12 * tmp36 - 6 * tmp12 * tmp41 - 6 * tmp22 * tmp27 - 6 * tmp22 * tmp30 - tmp46) *
                    (3 * AMAX * P0_y * tmp11 * tmp18 * tmp19 * tmp3 * tmp8 + 3 * AMAX * P1_y * tmp11 * tmp18 * tmp20 * tmp3 * tmp8 -
                     6 * tmp12 * tmp36 - 6 * tmp12 * tmp41 - 6 * tmp22 * tmp27 - 6 * tmp22 * tmp30 - tmp46)) +
                   ((-7.0 / 2.0 * P1_z * tmp57 - 1.0 / 2.0 * tmp19 * tmp39 + tmp19 * tmp44 + tmp20 * tmp42 + tmp42 * tmp6 * tmp63 + tmp44 * tmp69 -
                     1.0 / 2.0 * tmp54 - 7.0 / 2.0 * tmp56 + 2 * tmp60 * tmp63 + 2 * tmp65 * tmp66) *
                    (-7.0 / 2.0 * P1_z * tmp57 - 1.0 / 2.0 * tmp19 * tmp39 + tmp19 * tmp44 + tmp20 * tmp42 + tmp42 * tmp6 * tmp63 + tmp44 * tmp69 -
                     1.0 / 2.0 * tmp54 - 7.0 / 2.0 * tmp56 + 2 * tmp60 * tmp63 + 2 * tmp65 * tmp66)) +
                   2 * ((P0_x * tmp19 * tmp73 + P0_y * tmp19 * tmp75 + P1_x * tmp20 * tmp73 + P1_y * tmp20 * tmp75 - tmp27 * tmp70 - tmp30 * tmp70 -
                         tmp36 * tmp78 - tmp41 * tmp78 - tmp50 * tmp71 - tmp52 * tmp71 + tmp83) *
                        (P0_x * tmp19 * tmp73 + P0_y * tmp19 * tmp75 + P1_x * tmp20 * tmp73 + P1_y * tmp20 * tmp75 - tmp27 * tmp70 - tmp30 * tmp70 -
                         tmp36 * tmp78 - tmp41 * tmp78 - tmp50 * tmp71 - tmp52 * tmp71 + tmp83)) +
                   2 * ((-P0_x * tmp84 - 1.0 / 4.0 * P0_x * tmp86 + P0_z * tmp19 * tmp75 - 1.0 / 4.0 * P1_x * tmp20 * tmp6 -
                         7.0 / 4.0 * P1_x * tmp57 + P1_z * tmp20 * tmp75 - tmp24 * tmp60 * tmp70 + tmp25 * tmp35 * tmp63 + tmp25 * tmp40 * tmp66 -
                         tmp28 * tmp65 * tmp70 - tmp36 * tmp6 * tmp70 + tmp49 * tmp63 + tmp51 * tmp66 - tmp70 * tmp95 + tmp83) *
                        (-P0_x * tmp84 - 1.0 / 4.0 * P0_x * tmp86 + P0_z * tmp19 * tmp75 - 1.0 / 4.0 * P1_x * tmp20 * tmp6 -
                         7.0 / 4.0 * P1_x * tmp57 + P1_z * tmp20 * tmp75 - tmp24 * tmp60 * tmp70 + tmp25 * tmp35 * tmp63 + tmp25 * tmp40 * tmp66 -
                         tmp28 * tmp65 * tmp70 - tmp36 * tmp6 * tmp70 + tmp49 * tmp63 + tmp51 * tmp66 - tmp70 * tmp95 + tmp83)) +
                   2 * ((-P0_y * tmp84 - 1.0 / 4.0 * P0_y * tmp86 + P0_z * tmp19 * tmp73 - 1.0 / 4.0 * P1_y * tmp20 * tmp6 -
                         7.0 / 4.0 * P1_y * tmp57 + P1_z * tmp20 * tmp73 + tmp22 * tmp35 * tmp63 + tmp22 * tmp40 * tmp66 - tmp24 * tmp60 * tmp71 +
                         tmp26 * tmp63 - tmp28 * tmp65 * tmp71 + tmp29 * tmp66 - tmp36 * tmp6 * tmp71 - tmp71 * tmp95 + tmp83) *
                        (-P0_y * tmp84 - 1.0 / 4.0 * P0_y * tmp86 + P0_z * tmp19 * tmp73 - 1.0 / 4.0 * P1_y * tmp20 * tmp6 -
                         7.0 / 4.0 * P1_y * tmp57 + P1_z * tmp20 * tmp73 + tmp22 * tmp35 * tmp63 + tmp22 * tmp40 * tmp66 - tmp24 * tmp60 * tmp71 +
                         tmp26 * tmp63 - tmp28 * tmp65 * tmp71 + tmp29 * tmp66 - tmp36 * tmp6 * tmp71 - tmp71 * tmp95 + tmp83));
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
 * Kernel to initialize auxillary grid functions at all grid points.
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
  const REAL AMAX = params->AMAX;
  const REAL SINHWAA = params->SINHWAA;
  const REAL bScale = params->bScale;

#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        /*
         *  Original SymPy expressions:
         *  "[const REAL dsmin0 = AMAX*dxx0*(exp(xx0/SINHWAA)/SINHWAA + exp(-xx0/SINHWAA)/SINHWAA)*sqrt(AMAX**2*(exp(xx0/SINHWAA) -
         * exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) -
         * exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA)))]"
         *  "[const REAL dsmin1 = dxx1*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 +
         * bScale**2*sin(xx1)**2)]"
         *  "[const REAL dsmin2 = AMAX*dxx2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)/(exp(1/SINHWAA) - exp(-1/SINHWAA))]"
         */
        const REAL tmp1 = sin(xx1);
        const REAL tmp2 = (1.0 / (SINHWAA));
        const REAL tmp3 = exp(tmp2) - exp(-tmp2);
        const REAL tmp5 = exp(tmp2 * xx0);
        const REAL tmp6 = exp(-tmp2 * xx0);
        const REAL tmp10 = AMAX / tmp3;
        const REAL tmp7 = tmp5 - tmp6;
        const REAL tmp8 = ((AMAX) * (AMAX)) * ((tmp7) * (tmp7)) / ((tmp3) * (tmp3));
        const REAL tmp9 = sqrt(((bScale) * (bScale)) * ((tmp1) * (tmp1)) + tmp8);
        const REAL dsmin0 = dxx0 * tmp10 * tmp9 * (tmp2 * tmp5 + tmp2 * tmp6) / sqrt(((bScale) * (bScale)) + tmp8);
        const REAL dsmin1 = dxx1 * tmp9;
        const REAL dsmin2 = dxx2 * tmp1 * tmp10 * tmp7;

        // Set local wavespeed
        in_gfs[IDX4(VARIABLE_WAVESPEEDGF, i0, i1, i2)] = MINIMUM_GLOBAL_WAVESPEED * MIN(dsmin0, MIN(dsmin1, dsmin2)) / dt;

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION variable_wavespeed_gfs_all_points_host

/**
 * Call functions that set up all AUXEVOL gridfunctions.
 */
void initialize_constant_auxevol__rfm__SinhSymTP(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                                                 MoL_gridfunctions_struct *restrict gridfuncs) {
#include "../set_CodeParameters.h"

  REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
  REAL *restrict x0 = xx[0];
  REAL *restrict x1 = xx[1];
  REAL *restrict x2 = xx[2];

  // Set up variable wavespeed
  variable_wavespeed_gfs_all_points_host(params, x0, x1, x2, auxevol_gfs, dt, MINIMUM_GLOBAL_WAVESPEED);

  // Set up all other AUXEVOL gridfunctions
  auxevol_gfs_all_points_host(commondata, params, x0, x1, x2, auxevol_gfs);

} // END FUNCTION initialize_constant_auxevol__rfm__SinhSymTP
