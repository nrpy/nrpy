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
  const REAL tmp0 = (1.0 / (SINHWAA));
  const REAL tmp7 = cos(xx2);
  const REAL tmp8 = sin(xx1);
  const REAL tmp11 = sin(xx2);
  const REAL tmp1 = exp(tmp0) - exp(-tmp0);
  const REAL tmp3 = exp(tmp0 * xx0) - exp(-tmp0 * xx0);
  const REAL tmp18 = (1.0 / (tmp1));
  const REAL tmp4 = ((AMAX) * (AMAX)) * ((tmp3) * (tmp3)) / ((tmp1) * (tmp1));
  const REAL tmp19 = AMAX * tmp18 * tmp3 * tmp8;
  const REAL tmp5 = sqrt(((bScale) * (bScale)) + tmp4) * cos(xx1);
  const REAL tmp9 = tmp4 * ((tmp8) * (tmp8));
  const REAL tmp20 = tmp19 * tmp7;
  const REAL tmp22 = tmp11 * tmp19;
  const REAL tmp6 = tmp5 + zPunc;
  const REAL tmp10 = ((tmp7) * (tmp7)) * tmp9;
  const REAL tmp12 = ((tmp11) * (tmp11)) * tmp9;
  const REAL tmp15 = -tmp5 + zPunc;
  const REAL tmp52 = AMAX * S1_y * tmp18 * tmp3 * tmp7 * tmp8 - S1_x * tmp22;
  const REAL tmp55 = AMAX * S0_y * tmp18 * tmp3 * tmp7 * tmp8 - S0_x * tmp22;
  const REAL tmp76 = _Integer_6 * tmp11 * tmp7 * tmp9;
  const REAL tmp85 = _Rational_1_4 * tmp22;
  const REAL tmp88 = _Rational_1_4 * tmp20;
  const REAL tmp23 = S0_y * tmp15 + S0_z * tmp22;
  const REAL tmp30 = S1_x * tmp6 - S1_z * tmp20;
  const REAL tmp14 = tmp10 + tmp12 + ((tmp6) * (tmp6));
  const REAL tmp17 = tmp10 + tmp12 + ((tmp15) * (tmp15));
  const REAL tmp27 = -S1_y * tmp6 + S1_z * tmp22;
  const REAL tmp32 = -S0_x * tmp15 - S0_z * tmp20;
  const REAL tmp38 = -P0_z * tmp15;
  const REAL tmp43 = -P1_x * tmp20 - P1_y * tmp22 - P1_z * tmp6;
  const REAL tmp21 = pow(tmp17, -5.0 / 2.0);
  const REAL tmp25 = pow(tmp14, -5.0 / 2.0);
  const REAL tmp35 = pow(tmp14, -3.0 / 2.0);
  const REAL tmp37 = pow(tmp17, -3.0 / 2.0);
  const REAL tmp47 = -P0_x * tmp20 - P0_y * tmp22 - tmp38;
  const REAL tmp39 = P1_z * tmp35 * tmp6 + tmp37 * tmp38;
  const REAL tmp44 = tmp25 * tmp43;
  const REAL tmp48 = tmp21 * tmp47;
  const REAL tmp53 = -tmp25 * tmp6;
  const REAL tmp56 = tmp15 * tmp21;
  const REAL tmp63 = tmp15 * tmp37;
  const REAL tmp64 = -tmp35 * tmp6;
  const REAL tmp66 = P0_y * tmp22 * tmp37;
  const REAL tmp67 = P1_y * tmp22 * tmp35;
  const REAL tmp69 = P0_x * tmp20 * tmp37;
  const REAL tmp70 = P1_x * tmp20 * tmp35;
  const REAL tmp82 = -tmp15 * tmp37;
  const REAL tmp29 = tmp20 * tmp21 * tmp23 + tmp20 * tmp25 * tmp27;
  const REAL tmp34 = tmp21 * tmp22 * tmp32 + tmp22 * tmp25 * tmp30;
  const REAL tmp49 = _Integer_2 * tmp12 * tmp44 + _Integer_2 * tmp12 * tmp48;
  const REAL tmp51 = _Integer_2 * tmp10 * tmp44 + _Integer_2 * tmp10 * tmp48;
  const REAL tmp58 = tmp43 * tmp53;
  const REAL tmp61 = tmp47 * tmp56;
  const REAL tmp65 = P0_z * tmp63 + P1_z * tmp64;
  const REAL tmp68 = -tmp66 - tmp67;
  const REAL tmp71 = -tmp69 - tmp70;
  const REAL tmp87 = -tmp15 * tmp48;
  const REAL tmp62 = -tmp15 * tmp61;
  const REAL tmp74 = _Integer_3 * (-_Rational_1_4 * _Rational_2_3 * tmp58 * tmp6 - _Rational_1_4 * _Rational_2_3 * tmp62 -
                                   _Rational_2_3 * tmp52 * tmp53 - _Rational_2_3 * tmp55 * tmp56) +
                     _Rational_1_6 * tmp39 + _Rational_7_6 * tmp65;
  *psi_background = _Integer_1 + _Rational_1_2 * (bare_mass_0 / sqrt(tmp17) + bare_mass_1 / sqrt(tmp14));
  *ADD_times_AUU = _Integer_2 * (((_Integer_3 * (-tmp20 * tmp21 * tmp32 - tmp20 * tmp25 * tmp30 - tmp21 * tmp22 * tmp23 - tmp22 * tmp25 * tmp27) +
                                   _Rational_1_4 * (-tmp44 * tmp76 - tmp48 * tmp76) +
                                   _Rational_3_2 * (P0_x * tmp22 * tmp37 + P0_y * tmp20 * tmp37 + P1_x * tmp22 * tmp35 + P1_y * tmp20 * tmp35)) *
                                  (_Integer_3 * (-tmp20 * tmp21 * tmp32 - tmp20 * tmp25 * tmp30 - tmp21 * tmp22 * tmp23 - tmp22 * tmp25 * tmp27) +
                                   _Rational_1_4 * (-tmp44 * tmp76 - tmp48 * tmp76) +
                                   _Rational_3_2 * (P0_x * tmp22 * tmp37 + P0_y * tmp20 * tmp37 + P1_x * tmp22 * tmp35 + P1_y * tmp20 * tmp35))) +
                                 ((_Integer_3 * (-tmp20 * tmp21 * tmp55 - tmp20 * tmp25 * tmp52 + tmp23 * tmp56 + tmp27 * tmp53 -
                                                 tmp44 * tmp6 * tmp88 + tmp58 * tmp88 + tmp61 * tmp88 - tmp87 * tmp88) +
                                   _Rational_1_4 * (-P0_x * tmp82 - P1_x * tmp35 * tmp6) +
                                   _Rational_3_2 * (P0_z * tmp20 * tmp37 + P1_z * tmp20 * tmp35) + _Rational_7_4 * (-P0_x * tmp63 - P1_x * tmp64)) *
                                  (_Integer_3 * (-tmp20 * tmp21 * tmp55 - tmp20 * tmp25 * tmp52 + tmp23 * tmp56 + tmp27 * tmp53 -
                                                 tmp44 * tmp6 * tmp88 + tmp58 * tmp88 + tmp61 * tmp88 - tmp87 * tmp88) +
                                   _Rational_1_4 * (-P0_x * tmp82 - P1_x * tmp35 * tmp6) +
                                   _Rational_3_2 * (P0_z * tmp20 * tmp37 + P1_z * tmp20 * tmp35) + _Rational_7_4 * (-P0_x * tmp63 - P1_x * tmp64))) +
                                 ((_Integer_3 * (-tmp21 * tmp22 * tmp55 - tmp22 * tmp25 * tmp52 + tmp30 * tmp53 + tmp32 * tmp56 -
                                                 tmp44 * tmp6 * tmp85 + tmp58 * tmp85 + tmp61 * tmp85 - tmp85 * tmp87) +
                                   _Rational_1_4 * (-P0_y * tmp82 - P1_y * tmp35 * tmp6) +
                                   _Rational_3_2 * (P0_z * tmp22 * tmp37 + P1_z * tmp22 * tmp35) + _Rational_7_4 * (-P0_y * tmp63 - P1_y * tmp64)) *
                                  (_Integer_3 * (-tmp21 * tmp22 * tmp55 - tmp22 * tmp25 * tmp52 + tmp30 * tmp53 + tmp32 * tmp56 -
                                                 tmp44 * tmp6 * tmp85 + tmp58 * tmp85 + tmp61 * tmp85 - tmp85 * tmp87) +
                                   _Rational_1_4 * (-P0_y * tmp82 - P1_y * tmp35 * tmp6) +
                                   _Rational_3_2 * (P0_z * tmp22 * tmp37 + P1_z * tmp22 * tmp35) + _Rational_7_4 * (-P0_y * tmp63 - P1_y * tmp64)))) +
                   ((_Integer_2 * (tmp29 + tmp66 + tmp67) - _Integer_4 * tmp34 +
                     _Rational_1_4 * (-_Integer_4 * tmp12 * tmp44 - _Integer_4 * tmp12 * tmp48 + tmp51) + tmp71 + tmp74) *
                    (_Integer_2 * (tmp29 + tmp66 + tmp67) - _Integer_4 * tmp34 +
                     _Rational_1_4 * (-_Integer_4 * tmp12 * tmp44 - _Integer_4 * tmp12 * tmp48 + tmp51) + tmp71 + tmp74)) +
                   ((_Integer_2 * (tmp34 + tmp69 + tmp70) - _Integer_4 * tmp29 +
                     _Rational_1_4 * (-_Integer_4 * tmp10 * tmp44 - _Integer_4 * tmp10 * tmp48 + tmp49) + tmp68 + tmp74) *
                    (_Integer_2 * (tmp34 + tmp69 + tmp70) - _Integer_4 * tmp29 +
                     _Rational_1_4 * (-_Integer_4 * tmp10 * tmp44 - _Integer_4 * tmp10 * tmp48 + tmp49) + tmp68 + tmp74)) +
                   ((_Integer_2 * (tmp29 + tmp34) - _Rational_1_3 * tmp39 + _Rational_1_4 * (tmp49 + tmp51) +
                     _Rational_4_3 * (_Integer_3 * _Rational_1_4 * tmp58 * tmp6 + _Integer_3 * _Rational_1_4 * tmp62 + _Integer_3 * tmp52 * tmp53 +
                                      _Integer_3 * tmp55 * tmp56) -
                     _Rational_7_3 * tmp65 + tmp68 + tmp71) *
                    (_Integer_2 * (tmp29 + tmp34) - _Rational_1_3 * tmp39 + _Rational_1_4 * (tmp49 + tmp51) +
                     _Rational_4_3 * (_Integer_3 * _Rational_1_4 * tmp58 * tmp6 + _Integer_3 * _Rational_1_4 * tmp62 + _Integer_3 * tmp52 * tmp53 +
                                      _Integer_3 * tmp55 * tmp56) -
                     _Rational_7_3 * tmp65 + tmp68 + tmp71));
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
