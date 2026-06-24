#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Kernel: auxevol_gfs_single_point_gpu.
 * Kernel to compute AUXEVOL grid functions at a single point.
 */
__device__ static void auxevol_gfs_single_point_gpu(const size_t streamid, const REAL xx0, const REAL xx1, const REAL xx2,
                                                    REAL *restrict psi_background, REAL *restrict ADD_times_AUU) {
  // Load necessary parameters from params_struct
  const REAL AMAX = d_params[streamid].AMAX;
  const REAL SINHWAA = d_params[streamid].SINHWAA;
  const REAL bScale = d_params[streamid].bScale;

  // Load necessary parameters from commondata_struct
  const REAL P0_x = d_commondata.P0_x;
  const REAL P0_y = d_commondata.P0_y;
  const REAL P0_z = d_commondata.P0_z;
  const REAL P1_x = d_commondata.P1_x;
  const REAL P1_y = d_commondata.P1_y;
  const REAL P1_z = d_commondata.P1_z;
  const REAL S0_x = d_commondata.S0_x;
  const REAL S0_y = d_commondata.S0_y;
  const REAL S0_z = d_commondata.S0_z;
  const REAL S1_x = d_commondata.S1_x;
  const REAL S1_y = d_commondata.S1_y;
  const REAL S1_z = d_commondata.S1_z;
  const REAL bare_mass_0 = d_commondata.bare_mass_0;
  const REAL bare_mass_1 = d_commondata.bare_mass_1;
  const REAL zPunc = d_commondata.zPunc;

  const REAL tmp0 = ((AMAX) * (AMAX));
  const REAL tmp1 = (1.0 / (SINHWAA));
  const REAL tmp10 = cos(xx2);
  const REAL tmp12 = sin(xx1);
  const REAL tmp16 = sin(xx2);
  const REAL tmp2 = exp(tmp1) - exp(-tmp1);
  const REAL tmp11 = ((tmp10) * (tmp10));
  const REAL tmp13 = ((tmp12) * (tmp12));
  const REAL tmp3 = (1.0 / ((tmp2) * (tmp2)));
  const REAL tmp5 = exp(tmp1 * xx0) - exp(-tmp1 * xx0);
  const REAL tmp23 = (1.0 / (tmp2));
  const REAL tmp6 = ((tmp5) * (tmp5));
  const REAL tmp26 = AMAX * tmp12 * tmp23 * tmp5;
  const REAL tmp7 = tmp0 * tmp3 * tmp6;
  const REAL tmp27 = tmp10 * tmp26;
  const REAL tmp32 = tmp16 * tmp26;
  const REAL tmp8 = sqrt(((bScale) * (bScale)) + tmp7) * cos(xx1);
  const REAL tmp14 = tmp13 * tmp7;
  const REAL tmp57 = AMAX * S1_y * tmp10 * tmp12 * tmp23 * tmp5 - S1_x * tmp32;
  const REAL tmp58 = AMAX * S0_y * tmp10 * tmp12 * tmp23 * tmp5 - S0_x * tmp32;
  const REAL tmp9 = tmp8 + zPunc;
  const REAL tmp15 = tmp11 * tmp14;
  const REAL tmp17 = tmp14 * ((tmp16) * (tmp16));
  const REAL tmp20 = -tmp8 + zPunc;
  const REAL tmp46 = 6 * tmp10 * tmp14 * tmp16;
  const REAL tmp28 = S1_x * tmp9 - S1_z * tmp27;
  const REAL tmp33 = S0_y * tmp20 + S0_z * tmp32;
  const REAL tmp61 = -3 * tmp8 - 3 * zPunc;
  const REAL tmp62 = -3 * tmp8 + 3 * zPunc;
  const REAL tmp19 = tmp15 + tmp17 + ((tmp9) * (tmp9));
  const REAL tmp22 = tmp15 + tmp17 + ((tmp20) * (tmp20));
  const REAL tmp38 = -S1_y * tmp9 + S1_z * tmp32;
  const REAL tmp40 = -S0_x * tmp20 - S0_z * tmp27;
  const REAL tmp45 = -1.0 / 4.0 * P1_x * tmp27 - 1.0 / 4.0 * P1_y * tmp32 - 1.0 / 4.0 * P1_z * tmp9;
  const REAL tmp49 = -P0_z * tmp20;
  const REAL tmp24 = pow(tmp19, -3.0 / 2.0);
  const REAL tmp25 = pow(tmp22, -3.0 / 2.0);
  const REAL tmp29 = pow(tmp19, -5.0 / 2.0);
  const REAL tmp34 = pow(tmp22, -5.0 / 2.0);
  const REAL tmp50 = -1.0 / 4.0 * P0_x * tmp27 - 1.0 / 4.0 * P0_y * tmp32 - 1.0 / 4.0 * tmp49;
  const REAL tmp54 = -tmp24 * tmp9;
  const REAL tmp56 = -tmp20 * tmp25;
  const REAL tmp64 = -tmp20 * tmp50;
  const REAL tmp65 = P1_x * tmp24 * tmp27;
  const REAL tmp66 = P0_x * tmp25 * tmp27;
  const REAL tmp68 = 4 * tmp29;
  const REAL tmp73 = 4 * tmp34;
  const REAL tmp84 = P1_z * tmp24 * tmp9;
  const REAL tmp89 = 2 * tmp29;
  const REAL tmp90 = 2 * tmp34;
  const REAL tmp31 = 3 * tmp27 * tmp29;
  const REAL tmp36 = 3 * tmp32 * tmp34;
  const REAL tmp39 = 3 * tmp29 * tmp32;
  const REAL tmp41 = 3 * tmp27 * tmp34;
  const REAL tmp53 = (7.0 / 4.0) * tmp20 * tmp25;
  const REAL tmp70 = tmp29 * tmp57 * tmp61;
  const REAL tmp75 = tmp34 * tmp58 * tmp62;
  const REAL tmp79 = tmp29 * tmp45 * tmp61 * tmp9;
  const REAL tmp81 = -tmp20 * tmp34 * tmp50 * tmp62;
  const REAL tmp85 = P0_z * tmp20 * tmp25;
  const REAL tmp93 =
      -P0_y * tmp25 * tmp32 - P1_y * tmp24 * tmp32 + tmp17 * tmp45 * tmp89 + tmp17 * tmp50 * tmp90 + tmp28 * tmp32 * tmp89 + tmp32 * tmp40 * tmp90;
  *psi_background = (1.0 / 2.0) * bare_mass_0 / sqrt(tmp22) + (1.0 / 2.0) * bare_mass_1 / sqrt(tmp19) + 1;
  *ADD_times_AUU =
      2 * (((3.0 / 2.0) * AMAX * P0_x * tmp12 * tmp16 * tmp23 * tmp25 * tmp5 + (3.0 / 2.0) * AMAX * P0_y * tmp10 * tmp12 * tmp23 * tmp25 * tmp5 +
            (3.0 / 2.0) * AMAX * P1_x * tmp12 * tmp16 * tmp23 * tmp24 * tmp5 + (3.0 / 2.0) * AMAX * P1_y * tmp10 * tmp12 * tmp23 * tmp24 * tmp5 -
            tmp28 * tmp31 - tmp29 * tmp45 * tmp46 - tmp33 * tmp36 - tmp34 * tmp46 * tmp50 - tmp38 * tmp39 - tmp40 * tmp41) *
           ((3.0 / 2.0) * AMAX * P0_x * tmp12 * tmp16 * tmp23 * tmp25 * tmp5 + (3.0 / 2.0) * AMAX * P0_y * tmp10 * tmp12 * tmp23 * tmp25 * tmp5 +
            (3.0 / 2.0) * AMAX * P1_x * tmp12 * tmp16 * tmp23 * tmp24 * tmp5 + (3.0 / 2.0) * AMAX * P1_y * tmp10 * tmp12 * tmp23 * tmp24 * tmp5 -
            tmp28 * tmp31 - tmp29 * tmp45 * tmp46 - tmp33 * tmp36 - tmp34 * tmp46 * tmp50 - tmp38 * tmp39 - tmp40 * tmp41)) +
      2 * (((3.0 / 2.0) * AMAX * P0_z * tmp10 * tmp12 * tmp23 * tmp25 * tmp5 + (3.0 / 2.0) * AMAX * P1_z * tmp10 * tmp12 * tmp23 * tmp24 * tmp5 +
            AMAX * tmp10 * tmp12 * tmp23 * tmp29 * tmp45 * tmp5 * tmp61 + AMAX * tmp10 * tmp12 * tmp23 * tmp34 * tmp5 * tmp50 * tmp62 - P0_x * tmp53 -
            1.0 / 4.0 * P0_x * tmp56 - 1.0 / 4.0 * P1_x * tmp24 * tmp9 - 7.0 / 4.0 * P1_x * tmp54 + tmp29 * tmp38 * tmp61 - tmp31 * tmp45 * tmp9 -
            tmp31 * tmp57 + tmp33 * tmp34 * tmp62 - tmp41 * tmp58 - tmp41 * tmp64) *
           ((3.0 / 2.0) * AMAX * P0_z * tmp10 * tmp12 * tmp23 * tmp25 * tmp5 + (3.0 / 2.0) * AMAX * P1_z * tmp10 * tmp12 * tmp23 * tmp24 * tmp5 +
            AMAX * tmp10 * tmp12 * tmp23 * tmp29 * tmp45 * tmp5 * tmp61 + AMAX * tmp10 * tmp12 * tmp23 * tmp34 * tmp5 * tmp50 * tmp62 - P0_x * tmp53 -
            1.0 / 4.0 * P0_x * tmp56 - 1.0 / 4.0 * P1_x * tmp24 * tmp9 - 7.0 / 4.0 * P1_x * tmp54 + tmp29 * tmp38 * tmp61 - tmp31 * tmp45 * tmp9 -
            tmp31 * tmp57 + tmp33 * tmp34 * tmp62 - tmp41 * tmp58 - tmp41 * tmp64)) +
      2 * (((3.0 / 2.0) * AMAX * P0_z * tmp12 * tmp16 * tmp23 * tmp25 * tmp5 + (3.0 / 2.0) * AMAX * P1_z * tmp12 * tmp16 * tmp23 * tmp24 * tmp5 +
            AMAX * tmp12 * tmp16 * tmp23 * tmp29 * tmp45 * tmp5 * tmp61 + AMAX * tmp12 * tmp16 * tmp23 * tmp34 * tmp5 * tmp50 * tmp62 - P0_y * tmp53 -
            1.0 / 4.0 * P0_y * tmp56 - 1.0 / 4.0 * P1_y * tmp24 * tmp9 - 7.0 / 4.0 * P1_y * tmp54 + tmp28 * tmp29 * tmp61 + tmp34 * tmp40 * tmp62 -
            tmp36 * tmp58 - tmp36 * tmp64 - tmp39 * tmp45 * tmp9 - tmp39 * tmp57) *
           ((3.0 / 2.0) * AMAX * P0_z * tmp12 * tmp16 * tmp23 * tmp25 * tmp5 + (3.0 / 2.0) * AMAX * P1_z * tmp12 * tmp16 * tmp23 * tmp24 * tmp5 +
            AMAX * tmp12 * tmp16 * tmp23 * tmp29 * tmp45 * tmp5 * tmp61 + AMAX * tmp12 * tmp16 * tmp23 * tmp34 * tmp5 * tmp50 * tmp62 - P0_y * tmp53 -
            1.0 / 4.0 * P0_y * tmp56 - 1.0 / 4.0 * P1_y * tmp24 * tmp9 - 7.0 / 4.0 * P1_y * tmp54 + tmp28 * tmp29 * tmp61 + tmp34 * tmp40 * tmp62 -
            tmp36 * tmp58 - tmp36 * tmp64 - tmp39 * tmp45 * tmp9 - tmp39 * tmp57)) +
      ((-7.0 / 3.0 * P1_z * tmp54 + tmp15 * tmp45 * tmp89 + tmp15 * tmp50 * tmp90 - 1.0 / 3.0 * tmp25 * tmp49 + tmp27 * tmp33 * tmp90 +
        tmp27 * tmp38 * tmp89 - tmp65 - tmp66 + (4.0 / 3.0) * tmp70 + (4.0 / 3.0) * tmp75 + (4.0 / 3.0) * tmp79 + (4.0 / 3.0) * tmp81 -
        1.0 / 3.0 * tmp84 - 7.0 / 3.0 * tmp85 + tmp93) *
       (-7.0 / 3.0 * P1_z * tmp54 + tmp15 * tmp45 * tmp89 + tmp15 * tmp50 * tmp90 - 1.0 / 3.0 * tmp25 * tmp49 + tmp27 * tmp33 * tmp90 +
        tmp27 * tmp38 * tmp89 - tmp65 - tmp66 + (4.0 / 3.0) * tmp70 + (4.0 / 3.0) * tmp75 + (4.0 / 3.0) * tmp79 + (4.0 / 3.0) * tmp81 -
        1.0 / 3.0 * tmp84 - 7.0 / 3.0 * tmp85 + tmp93)) +
      (((7.0 / 6.0) * P1_z * tmp54 - tmp15 * tmp45 * tmp68 - tmp15 * tmp50 * tmp73 + (1.0 / 6.0) * tmp25 * tmp49 - tmp27 * tmp33 * tmp73 -
        tmp27 * tmp38 * tmp68 + 2 * tmp65 + 2 * tmp66 - 2.0 / 3.0 * tmp70 - 2.0 / 3.0 * tmp75 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp81 +
        (1.0 / 6.0) * tmp84 + (7.0 / 6.0) * tmp85 + tmp93) *
       ((7.0 / 6.0) * P1_z * tmp54 - tmp15 * tmp45 * tmp68 - tmp15 * tmp50 * tmp73 + (1.0 / 6.0) * tmp25 * tmp49 - tmp27 * tmp33 * tmp73 -
        tmp27 * tmp38 * tmp68 + 2 * tmp65 + 2 * tmp66 - 2.0 / 3.0 * tmp70 - 2.0 / 3.0 * tmp75 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp81 +
        (1.0 / 6.0) * tmp84 + (7.0 / 6.0) * tmp85 + tmp93)) +
      ((2 * AMAX * P0_y * tmp12 * tmp16 * tmp23 * tmp25 * tmp5 + 2 * AMAX * P1_y * tmp12 * tmp16 * tmp23 * tmp24 * tmp5 +
        2 * AMAX * tmp10 * tmp12 * tmp23 * tmp29 * tmp38 * tmp5 + 2 * AMAX * tmp10 * tmp12 * tmp23 * tmp33 * tmp34 * tmp5 + P0_z * tmp20 * tmp25 -
        P1_z * tmp24 * tmp9 + 2 * tmp0 * tmp11 * tmp13 * tmp29 * tmp3 * tmp45 * tmp6 + 2 * tmp0 * tmp11 * tmp13 * tmp3 * tmp34 * tmp50 * tmp6 -
        tmp17 * tmp45 * tmp68 - tmp17 * tmp50 * tmp73 - tmp28 * tmp32 * tmp68 - tmp32 * tmp40 * tmp73 - tmp65 - tmp66 - 2.0 / 3.0 * tmp70 -
        2.0 / 3.0 * tmp75 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp81) *
       (2 * AMAX * P0_y * tmp12 * tmp16 * tmp23 * tmp25 * tmp5 + 2 * AMAX * P1_y * tmp12 * tmp16 * tmp23 * tmp24 * tmp5 +
        2 * AMAX * tmp10 * tmp12 * tmp23 * tmp29 * tmp38 * tmp5 + 2 * AMAX * tmp10 * tmp12 * tmp23 * tmp33 * tmp34 * tmp5 + P0_z * tmp20 * tmp25 -
        P1_z * tmp24 * tmp9 + 2 * tmp0 * tmp11 * tmp13 * tmp29 * tmp3 * tmp45 * tmp6 + 2 * tmp0 * tmp11 * tmp13 * tmp3 * tmp34 * tmp50 * tmp6 -
        tmp17 * tmp45 * tmp68 - tmp17 * tmp50 * tmp73 - tmp28 * tmp32 * tmp68 - tmp32 * tmp40 * tmp73 - tmp65 - tmp66 - 2.0 / 3.0 * tmp70 -
        2.0 / 3.0 * tmp75 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp81));
} // END FUNCTION: auxevol_gfs_single_point_gpu

/**
 * Kernel: auxevol_gfs_all_points_gpu.
 * Kernel to initialize auxillary grid functions at all grid points.
 */
__global__ static void auxevol_gfs_all_points_gpu(const size_t streamid, const REAL *restrict x0, const REAL *restrict x1, const REAL *restrict x2,
                                                  REAL *restrict in_gfs) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = d_params[streamid].Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = d_params[streamid].invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = d_params[streamid].invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = d_params[streamid].invdxx2;

  MAYBE_UNUSED const int tid0 = blockIdx.x * blockDim.x + threadIdx.x;
  MAYBE_UNUSED const int tid1 = blockIdx.y * blockDim.y + threadIdx.y;
  MAYBE_UNUSED const int tid2 = blockIdx.z * blockDim.z + threadIdx.z;

  MAYBE_UNUSED const int stride0 = blockDim.x * gridDim.x;
  MAYBE_UNUSED const int stride1 = blockDim.y * gridDim.y;
  MAYBE_UNUSED const int stride2 = blockDim.z * gridDim.z;

  for (int i2 = tid2; i2 < Nxx_plus_2NGHOSTS2; i2 += stride2) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = tid1; i1 < Nxx_plus_2NGHOSTS1; i1 += stride1) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = tid0; i0 < Nxx_plus_2NGHOSTS0; i0 += stride0) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        auxevol_gfs_single_point_gpu(streamid, xx0, xx1, xx2, &in_gfs[IDX4(PSI_BACKGROUNDGF, i0, i1, i2)],
                                     &in_gfs[IDX4(ADD_TIMES_AUUGF, i0, i1, i2)]);

      } // END LOOP: for i0 over [tid0, Nxx_plus_2NGHOSTS0)
    } // END LOOP: for i1 over [tid1, Nxx_plus_2NGHOSTS1)
  } // END LOOP: for i2 over [tid2, Nxx_plus_2NGHOSTS2)
} // END FUNCTION: auxevol_gfs_all_points_gpu

/**
 * Kernel: variable_wavespeed_gfs_all_points_gpu.
 * Kernel to compute variable wavespeed for all grids based on local grid spacing.
 */
__global__ static void variable_wavespeed_gfs_all_points_gpu(const size_t streamid, const REAL *restrict x0, const REAL *restrict x1,
                                                             const REAL *restrict x2, REAL *restrict in_gfs, const REAL dt,
                                                             const REAL MINIMUM_GLOBAL_WAVESPEED) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = d_params[streamid].Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = d_params[streamid].invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = d_params[streamid].invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = d_params[streamid].invdxx2;

  MAYBE_UNUSED const int tid0 = blockIdx.x * blockDim.x + threadIdx.x;
  MAYBE_UNUSED const int tid1 = blockIdx.y * blockDim.y + threadIdx.y;
  MAYBE_UNUSED const int tid2 = blockIdx.z * blockDim.z + threadIdx.z;

  MAYBE_UNUSED const int stride0 = blockDim.x * gridDim.x;
  MAYBE_UNUSED const int stride1 = blockDim.y * gridDim.y;
  MAYBE_UNUSED const int stride2 = blockDim.z * gridDim.z;

  // Load necessary parameters from params_struct
  const REAL AMAX = d_params[streamid].AMAX;
  const REAL SINHWAA = d_params[streamid].SINHWAA;
  const REAL bScale = d_params[streamid].bScale;

  for (int i2 = tid2 + NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2 += stride2) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = tid1 + NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1 += stride1) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = tid0 + NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += stride0) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        /*
         *  Original SymPy expressions:
         *  "[const REAL dsmin0 = AMAX*d_params[streamid].dxx0*(exp(xx0/SINHWAA)/SINHWAA + exp(-xx0/SINHWAA)/SINHWAA)*sqrt(AMAX**2*(exp(xx0/SINHWAA) -
         * exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)/(sqrt(AMAX**2*(exp(xx0/SINHWAA) -
         * exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA)))]"
         *  "[const REAL dsmin1 = d_params[streamid].dxx1*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2
         * + bScale**2*sin(xx1)**2)]"
         *  "[const REAL dsmin2 = AMAX*d_params[streamid].dxx2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)/(exp(1/SINHWAA) - exp(-1/SINHWAA))]"
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
        const REAL dsmin0 = d_params[streamid].dxx0 * tmp10 * tmp9 * (tmp2 * tmp5 + tmp2 * tmp6) / sqrt(((bScale) * (bScale)) + tmp8);
        const REAL dsmin1 = d_params[streamid].dxx1 * tmp9;
        const REAL dsmin2 = d_params[streamid].dxx2 * tmp1 * tmp10 * tmp7;

        // Set local wavespeed
        in_gfs[IDX4(VARIABLE_WAVESPEEDGF, i0, i1, i2)] = MINIMUM_GLOBAL_WAVESPEED * NRPYMIN(dsmin0, NRPYMIN(dsmin1, dsmin2)) / dt;

      } // END LOOP: for i0 over [tid0+NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS)
    } // END LOOP: for i1 over [tid1+NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS)
  } // END LOOP: for i2 over [tid2+NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
} // END FUNCTION: variable_wavespeed_gfs_all_points_gpu

/**
 * Call functions that set up all AUXEVOL gridfunctions.
 */
void auxevol_gfs_set_to_constant__rfm__SinhSymTP(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                                                 MoL_gridfunctions_struct *restrict gridfuncs) {
#include "set_CodeParameters.h"
  cpyHosttoDevice_commondata__constant(commondata);

  REAL *auxevol_gfs = gridfuncs->auxevol_gfs;
  REAL *x0 = xx[0];
  REAL *x1 = xx[1];
  REAL *x2 = xx[2];

  // Set up variable wavespeed
  {
    const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_NELL_WAVESPEED;
    const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_NELL_WAVESPEED;
    const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_NELL_WAVESPEED;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((params->Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                         (params->Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                         (params->Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    variable_wavespeed_gfs_all_points_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, x0, x1, x2, auxevol_gfs, dt,
                                                                                                         MINIMUM_GLOBAL_WAVESPEED);
    cudaCheckErrors(cudaKernel, "variable_wavespeed_gfs_all_points_gpu failure");
  }

  // Set up all other AUXEVOL gridfunctions
  {
    const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_NELL_AUX;
    const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_NELL_AUX;
    const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_NELL_AUX;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((params->Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                         (params->Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                         (params->Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    auxevol_gfs_all_points_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, x0, x1, x2, auxevol_gfs);
    cudaCheckErrors(cudaKernel, "auxevol_gfs_all_points_gpu failure");
  }

} // END FUNCTION: auxevol_gfs_set_to_constant__rfm__SinhSymTP
