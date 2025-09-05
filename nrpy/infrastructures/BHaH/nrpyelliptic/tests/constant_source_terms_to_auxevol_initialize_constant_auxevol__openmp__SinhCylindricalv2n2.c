#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Kernel: auxevol_gfs_single_point_host.
 * Kernel to compute AUXEVOL grid functions at a single point.
 */
static void auxevol_gfs_single_point_host(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xx0,
                                          const REAL xx1, const REAL xx2, REAL *restrict psi_background, REAL *restrict ADD_times_AUU) {
  // Load necessary parameters from params_struct
  const REAL AMPLRHO = params->AMPLRHO;
  const REAL AMPLZ = params->AMPLZ;
  const REAL SINHWRHO = params->SINHWRHO;
  const REAL SINHWZ = params->SINHWZ;
  const REAL rho_slope = params->rho_slope;
  const REAL z_slope = params->z_slope;

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

  const REAL tmp1 = (1.0 / (SINHWZ));
  const REAL tmp6 = cos(xx1);
  const REAL tmp8 = (1.0 / (SINHWRHO));
  const REAL tmp13 = sin(xx1);
  const REAL tmp7 = ((tmp6) * (tmp6));
  const REAL tmp14 = ((tmp13) * (tmp13));
  const REAL tmp3 = ((xx2) * (xx2)) * (AMPLZ - z_slope) * (exp(tmp1 * xx2) - exp(-tmp1 * xx2)) / (exp(tmp1) - exp(-tmp1));
  const REAL tmp10 = rho_slope * xx0 + ((xx0) * (xx0)) * (AMPLRHO - rho_slope) * (exp(tmp8 * xx0) - exp(-tmp8 * xx0)) / (exp(tmp8) - exp(-tmp8));
  const REAL tmp4 = tmp3 + xx2 * z_slope;
  const REAL tmp11 = ((tmp10) * (tmp10));
  const REAL tmp21 = tmp10 * tmp6;
  const REAL tmp23 = tmp10 * tmp13;
  const REAL tmp63 = 3 * tmp3 + 3 * xx2 * z_slope;
  const REAL tmp5 = tmp4 + zPunc;
  const REAL tmp12 = tmp11 * tmp7;
  const REAL tmp15 = tmp11 * tmp14;
  const REAL tmp18 = tmp4 - zPunc;
  const REAL tmp61 = -S1_x * tmp23 + S1_y * tmp10 * tmp6;
  const REAL tmp64 = -tmp63 - 3 * zPunc;
  const REAL tmp67 = -S0_x * tmp23 + S0_y * tmp10 * tmp6;
  const REAL tmp68 = -tmp63 + 3 * zPunc;
  const REAL tmp82 = 3 * tmp21;
  const REAL tmp83 = 3 * tmp23;
  const REAL tmp88 = 6 * tmp11 * tmp13 * tmp6;
  const REAL tmp34 = S1_x * tmp5 - S1_z * tmp21;
  const REAL tmp36 = S0_x * tmp18 - S0_z * tmp21;
  const REAL tmp17 = tmp12 + tmp15 + ((tmp5) * (tmp5));
  const REAL tmp19 = tmp12 + tmp15 + ((tmp18) * (tmp18));
  const REAL tmp25 = -1.0 / 4.0 * P1_x * tmp21 - 1.0 / 4.0 * P1_y * tmp23 - 1.0 / 4.0 * P1_z * tmp5;
  const REAL tmp30 = -1.0 / 4.0 * P0_x * tmp21 - 1.0 / 4.0 * P0_y * tmp23 - 1.0 / 4.0 * P0_z * tmp18;
  const REAL tmp52 = -S1_y * tmp5 + S1_z * tmp23;
  const REAL tmp55 = -S0_y * tmp18 + S0_z * tmp23;
  const REAL tmp26 = pow(tmp17, -5.0 / 2.0);
  const REAL tmp31 = pow(tmp19, -5.0 / 2.0);
  const REAL tmp38 = pow(tmp19, -3.0 / 2.0);
  const REAL tmp41 = pow(tmp17, -3.0 / 2.0);
  const REAL tmp44 = 2 * tmp25;
  const REAL tmp46 = 2 * tmp30;
  const REAL tmp27 = tmp25 * tmp26;
  const REAL tmp32 = tmp30 * tmp31;
  const REAL tmp35 = tmp26 * tmp34;
  const REAL tmp37 = tmp31 * tmp36;
  const REAL tmp39 = P0_y * tmp38;
  const REAL tmp42 = P1_y * tmp41;
  const REAL tmp53 = tmp26 * tmp52;
  const REAL tmp56 = tmp31 * tmp55;
  const REAL tmp57 = P0_x * tmp38;
  const REAL tmp59 = P1_x * tmp41;
  const REAL tmp71 = -tmp18 * tmp38;
  const REAL tmp73 = -tmp41 * tmp5;
  const REAL tmp75 = P0_z * tmp18 * tmp38;
  const REAL tmp76 = P1_z * tmp41 * tmp5;
  const REAL tmp40 = tmp23 * tmp39;
  const REAL tmp43 = tmp23 * tmp42;
  const REAL tmp48 = tmp38 * tmp46 + tmp41 * tmp44;
  const REAL tmp58 = tmp21 * tmp57;
  const REAL tmp60 = tmp21 * tmp59;
  const REAL tmp66 = tmp26 * tmp61 * tmp64;
  const REAL tmp70 = tmp31 * tmp67 * tmp68;
  const REAL tmp77 = tmp26 * tmp5 * tmp64;
  const REAL tmp78 = tmp18 * tmp31 * tmp68;
  const REAL tmp90 = (2.0 / 3.0) * tmp25 * tmp77;
  const REAL tmp91 = (2.0 / 3.0) * tmp30 * tmp78;
  const REAL tmp100 = P0_z * tmp18 * tmp38 + P1_z * tmp41 * tmp5 - 2 * tmp10 * tmp13 * tmp26 * tmp34 - 2 * tmp10 * tmp13 * tmp31 * tmp36 -
                      2 * tmp10 * tmp26 * tmp52 * tmp6 - 2 * tmp10 * tmp31 * tmp55 * tmp6 - 2 * tmp11 * tmp14 * tmp25 * tmp26 -
                      2 * tmp11 * tmp14 * tmp30 * tmp31 - 2 * tmp11 * tmp25 * tmp26 * tmp7 - 2 * tmp11 * tmp30 * tmp31 * tmp7 + tmp40 + tmp43 +
                      tmp48 + tmp58 + tmp60 + (2.0 / 3.0) * tmp66 + (2.0 / 3.0) * tmp70 + tmp90 + tmp91;
  *psi_background = (1.0 / 2.0) * bare_mass_0 / sqrt(tmp19) + (1.0 / 2.0) * bare_mass_1 / sqrt(tmp17) + 1;
  *ADD_times_AUU =
      ((-6 * tmp12 * tmp27 - 6 * tmp12 * tmp32 - 6 * tmp21 * tmp53 - 6 * tmp21 * tmp56 + tmp48 + 3 * tmp58 + 3 * tmp60) *
       (-6 * tmp12 * tmp27 - 6 * tmp12 * tmp32 - 6 * tmp21 * tmp53 - 6 * tmp21 * tmp56 + tmp48 + 3 * tmp58 + 3 * tmp60)) +
      ((-6 * tmp15 * tmp27 - 6 * tmp15 * tmp32 - 6 * tmp23 * tmp35 - 6 * tmp23 * tmp37 + 3 * tmp40 + 3 * tmp43 + tmp48) *
       (-6 * tmp15 * tmp27 - 6 * tmp15 * tmp32 - 6 * tmp23 * tmp35 - 6 * tmp23 * tmp37 + 3 * tmp40 + 3 * tmp43 + tmp48)) +
      ((-7.0 / 2.0 * P0_z * tmp71 - 7.0 / 2.0 * P1_z * tmp73 + tmp44 * tmp77 + tmp46 * tmp78 + tmp48 + 2 * tmp66 + 2 * tmp70 - 1.0 / 2.0 * tmp75 -
        1.0 / 2.0 * tmp76) *
       (-7.0 / 2.0 * P0_z * tmp71 - 7.0 / 2.0 * P1_z * tmp73 + tmp44 * tmp77 + tmp46 * tmp78 + tmp48 + 2 * tmp66 + 2 * tmp70 - 1.0 / 2.0 * tmp75 -
        1.0 / 2.0 * tmp76)) +
      2 * ((-7.0 / 4.0 * P0_x * tmp71 + (3.0 / 2.0) * P0_z * tmp10 * tmp38 * tmp6 - 7.0 / 4.0 * P1_x * tmp73 +
            (3.0 / 2.0) * P1_z * tmp10 * tmp41 * tmp6 + tmp10 * tmp25 * tmp26 * tmp6 * tmp64 + tmp10 * tmp30 * tmp31 * tmp6 * tmp68 - tmp100 -
            tmp18 * tmp32 * tmp82 - 1.0 / 4.0 * tmp18 * tmp57 + tmp26 * tmp52 * tmp64 - tmp26 * tmp61 * tmp82 - tmp27 * tmp5 * tmp82 +
            tmp31 * tmp55 * tmp68 - tmp31 * tmp67 * tmp82 - 1.0 / 4.0 * tmp5 * tmp59) *
           (-7.0 / 4.0 * P0_x * tmp71 + (3.0 / 2.0) * P0_z * tmp10 * tmp38 * tmp6 - 7.0 / 4.0 * P1_x * tmp73 +
            (3.0 / 2.0) * P1_z * tmp10 * tmp41 * tmp6 + tmp10 * tmp25 * tmp26 * tmp6 * tmp64 + tmp10 * tmp30 * tmp31 * tmp6 * tmp68 - tmp100 -
            tmp18 * tmp32 * tmp82 - 1.0 / 4.0 * tmp18 * tmp57 + tmp26 * tmp52 * tmp64 - tmp26 * tmp61 * tmp82 - tmp27 * tmp5 * tmp82 +
            tmp31 * tmp55 * tmp68 - tmp31 * tmp67 * tmp82 - 1.0 / 4.0 * tmp5 * tmp59)) +
      2 * ((-7.0 / 4.0 * P0_y * tmp71 + (3.0 / 2.0) * P0_z * tmp10 * tmp13 * tmp38 - 7.0 / 4.0 * P1_y * tmp73 +
            (3.0 / 2.0) * P1_z * tmp10 * tmp13 * tmp41 + tmp10 * tmp13 * tmp25 * tmp26 * tmp64 + tmp10 * tmp13 * tmp30 * tmp31 * tmp68 - tmp100 -
            tmp18 * tmp32 * tmp83 - 1.0 / 4.0 * tmp18 * tmp39 + tmp26 * tmp34 * tmp64 - tmp26 * tmp61 * tmp83 - tmp27 * tmp5 * tmp83 +
            tmp31 * tmp36 * tmp68 - tmp31 * tmp67 * tmp83 - 1.0 / 4.0 * tmp42 * tmp5) *
           (-7.0 / 4.0 * P0_y * tmp71 + (3.0 / 2.0) * P0_z * tmp10 * tmp13 * tmp38 - 7.0 / 4.0 * P1_y * tmp73 +
            (3.0 / 2.0) * P1_z * tmp10 * tmp13 * tmp41 + tmp10 * tmp13 * tmp25 * tmp26 * tmp64 + tmp10 * tmp13 * tmp30 * tmp31 * tmp68 - tmp100 -
            tmp18 * tmp32 * tmp83 - 1.0 / 4.0 * tmp18 * tmp39 + tmp26 * tmp34 * tmp64 - tmp26 * tmp61 * tmp83 - tmp27 * tmp5 * tmp83 +
            tmp31 * tmp36 * tmp68 - tmp31 * tmp67 * tmp83 - 1.0 / 4.0 * tmp42 * tmp5)) +
      2 * (((7.0 / 6.0) * P0_z * tmp71 + (7.0 / 6.0) * P1_z * tmp73 + tmp12 * tmp26 * tmp44 + tmp12 * tmp31 * tmp46 + tmp15 * tmp26 * tmp44 +
            tmp15 * tmp31 * tmp46 + (3.0 / 2.0) * tmp21 * tmp39 + (3.0 / 2.0) * tmp21 * tmp42 + 2 * tmp21 * tmp53 + 2 * tmp21 * tmp56 +
            2 * tmp23 * tmp35 + 2 * tmp23 * tmp37 + (3.0 / 2.0) * tmp23 * tmp57 + (3.0 / 2.0) * tmp23 * tmp59 - tmp27 * tmp88 - tmp32 * tmp88 -
            tmp35 * tmp82 - tmp37 * tmp82 - tmp38 * tmp46 - tmp40 - tmp41 * tmp44 - tmp43 - tmp53 * tmp83 - tmp56 * tmp83 - tmp58 - tmp60 -
            2.0 / 3.0 * tmp66 - 2.0 / 3.0 * tmp70 + (1.0 / 6.0) * tmp75 + (1.0 / 6.0) * tmp76 - tmp90 - tmp91) *
           ((7.0 / 6.0) * P0_z * tmp71 + (7.0 / 6.0) * P1_z * tmp73 + tmp12 * tmp26 * tmp44 + tmp12 * tmp31 * tmp46 + tmp15 * tmp26 * tmp44 +
            tmp15 * tmp31 * tmp46 + (3.0 / 2.0) * tmp21 * tmp39 + (3.0 / 2.0) * tmp21 * tmp42 + 2 * tmp21 * tmp53 + 2 * tmp21 * tmp56 +
            2 * tmp23 * tmp35 + 2 * tmp23 * tmp37 + (3.0 / 2.0) * tmp23 * tmp57 + (3.0 / 2.0) * tmp23 * tmp59 - tmp27 * tmp88 - tmp32 * tmp88 -
            tmp35 * tmp82 - tmp37 * tmp82 - tmp38 * tmp46 - tmp40 - tmp41 * tmp44 - tmp43 - tmp53 * tmp83 - tmp56 * tmp83 - tmp58 - tmp60 -
            2.0 / 3.0 * tmp66 - 2.0 / 3.0 * tmp70 + (1.0 / 6.0) * tmp75 + (1.0 / 6.0) * tmp76 - tmp90 - tmp91));
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
  const REAL AMPLRHO = params->AMPLRHO;
  const REAL AMPLZ = params->AMPLZ;
  const REAL SINHWRHO = params->SINHWRHO;
  const REAL SINHWZ = params->SINHWZ;
  const REAL rho_slope = params->rho_slope;
  const REAL z_slope = params->z_slope;

#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        /*
         *  Original SymPy expressions:
         *  "[const REAL dsmin0 = dxx0*(rho_slope + xx0**2*(AMPLRHO - rho_slope)*(exp(xx0/SINHWRHO)/SINHWRHO +
         * exp(-xx0/SINHWRHO)/SINHWRHO)/(exp(1/SINHWRHO) - exp(-1/SINHWRHO)) + 2*xx0*(AMPLRHO - rho_slope)*(exp(xx0/SINHWRHO) -
         * exp(-xx0/SINHWRHO))/(exp(1/SINHWRHO) - exp(-1/SINHWRHO)))]"
         *  "[const REAL dsmin1 = dxx1*(rho_slope*xx0 + xx0**2*(AMPLRHO - rho_slope)*(exp(xx0/SINHWRHO) - exp(-xx0/SINHWRHO))/(exp(1/SINHWRHO) -
         * exp(-1/SINHWRHO)))]"
         *  "[const REAL dsmin2 = dxx2*(xx2**2*(AMPLZ - z_slope)*(exp(xx2/SINHWZ)/SINHWZ + exp(-xx2/SINHWZ)/SINHWZ)/(exp(1/SINHWZ) - exp(-1/SINHWZ)) +
         * 2*xx2*(AMPLZ - z_slope)*(exp(xx2/SINHWZ) - exp(-xx2/SINHWZ))/(exp(1/SINHWZ) - exp(-1/SINHWZ)) + z_slope)]"
         */
        const REAL tmp0 = (1.0 / (SINHWRHO));
        const REAL tmp7 = (1.0 / (SINHWZ));
        const REAL tmp5 = (AMPLRHO - rho_slope) / (exp(tmp0) - exp(-tmp0));
        const REAL tmp11 = (AMPLZ - z_slope) / (exp(tmp7) - exp(-tmp7));
        const REAL tmp2 = exp(tmp0 * xx0);
        const REAL tmp3 = exp(-tmp0 * xx0);
        const REAL tmp6 = tmp5 * ((xx0) * (xx0));
        const REAL tmp9 = exp(tmp7 * xx2);
        const REAL tmp10 = exp(-tmp7 * xx2);
        const REAL tmp4 = tmp2 - tmp3;
        const REAL dsmin0 = dxx0 * (rho_slope + 2 * tmp4 * tmp5 * xx0 + tmp6 * (tmp0 * tmp2 + tmp0 * tmp3));
        const REAL dsmin1 = dxx1 * (rho_slope * xx0 + tmp4 * tmp6);
        const REAL dsmin2 = dxx2 * (tmp11 * ((xx2) * (xx2)) * (tmp10 * tmp7 + tmp7 * tmp9) + 2 * tmp11 * xx2 * (-tmp10 + tmp9) + z_slope);

        // Set local wavespeed
        in_gfs[IDX4(VARIABLE_WAVESPEEDGF, i0, i1, i2)] = MINIMUM_GLOBAL_WAVESPEED * MIN(dsmin0, MIN(dsmin1, dsmin2)) / dt;

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION variable_wavespeed_gfs_all_points_host

/**
 * Call functions that set up all AUXEVOL gridfunctions.
 */
void initialize_constant_auxevol__rfm__SinhCylindricalv2n2(commondata_struct *restrict commondata, params_struct *restrict params,
                                                           REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs) {
#include "set_CodeParameters.h"

  REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
  REAL *restrict x0 = xx[0];
  REAL *restrict x1 = xx[1];
  REAL *restrict x2 = xx[2];

  // Set up variable wavespeed
  variable_wavespeed_gfs_all_points_host(params, x0, x1, x2, auxevol_gfs, dt, MINIMUM_GLOBAL_WAVESPEED);

  // Set up all other AUXEVOL gridfunctions
  auxevol_gfs_all_points_host(commondata, params, x0, x1, x2, auxevol_gfs);

} // END FUNCTION initialize_constant_auxevol__rfm__SinhCylindricalv2n2
