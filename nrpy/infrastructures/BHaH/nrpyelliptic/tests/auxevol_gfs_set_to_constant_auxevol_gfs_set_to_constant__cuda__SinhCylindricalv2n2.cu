#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Kernel: auxevol_gfs_single_point_gpu.
 * Kernel to compute AUXEVOL grid functions at a single point.
 */
__device__ static void auxevol_gfs_single_point_gpu(const size_t streamid, const REAL xx0, const REAL xx1, const REAL xx2,
                                                    REAL *restrict psi_background, REAL *restrict ADD_times_AUU) {
  // Load necessary parameters from params_struct
  const REAL AMPLRHO = d_params[streamid].AMPLRHO;
  const REAL AMPLZ = d_params[streamid].AMPLZ;
  const REAL SINHWRHO = d_params[streamid].SINHWRHO;
  const REAL SINHWZ = d_params[streamid].SINHWZ;
  const REAL rho_slope = d_params[streamid].rho_slope;
  const REAL z_slope = d_params[streamid].z_slope;

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
  static const REAL _Integer_2 = 2;
  static const REAL _Rational_1_4 = 1.0 / 4.0;
  static const REAL _Rational_2_3 = 2.0 / 3.0;
  static const REAL _Rational_1_6 = 1.0 / 6.0;
  static const REAL _Rational_7_6 = 7.0 / 6.0;
  static const REAL _Integer_4 = 4;
  static const REAL _Rational_3_2 = 3.0 / 2.0;
  static const REAL _Integer_6 = 6;
  static const REAL _Rational_7_4 = 7.0 / 4.0;
  const REAL tmp0 = (1.0 / (SINHWZ));
  const REAL tmp4 = cos(xx1);
  const REAL tmp5 = (1.0 / (SINHWRHO));
  const REAL tmp10 = sin(xx1);
  const REAL tmp2 = ((xx2) * (xx2)) * (AMPLZ - z_slope) * (exp(tmp0 * xx2) - exp(-tmp0 * xx2)) / (exp(tmp0) - exp(-tmp0)) + xx2 * z_slope;
  const REAL tmp7 = rho_slope * xx0 + ((xx0) * (xx0)) * (AMPLRHO - rho_slope) * (exp(tmp5 * xx0) - exp(-tmp5 * xx0)) / (exp(tmp5) - exp(-tmp5));
  const REAL tmp3 = tmp2 + zPunc;
  const REAL tmp8 = ((tmp7) * (tmp7));
  const REAL tmp14 = tmp2 - zPunc;
  const REAL tmp19 = tmp4 * tmp7;
  const REAL tmp21 = tmp10 * tmp7;
  const REAL tmp9 = ((tmp4) * (tmp4)) * tmp8;
  const REAL tmp11 = ((tmp10) * (tmp10)) * tmp8;
  const REAL tmp31 = S1_x * tmp3 - S1_z * tmp19;
  const REAL tmp33 = S0_x * tmp14 - S0_z * tmp19;
  const REAL tmp49 = -S0_x * tmp21 + S0_y * tmp4 * tmp7;
  const REAL tmp52 = -S1_x * tmp21 + S1_y * tmp4 * tmp7;
  const REAL tmp73 = _Integer_6 * tmp10 * tmp4 * tmp8;
  const REAL tmp82 = _Rational_1_4 * tmp21;
  const REAL tmp85 = _Rational_1_4 * tmp19;
  const REAL tmp23 = -P1_x * tmp19 - P1_y * tmp21 - P1_z * tmp3;
  const REAL tmp29 = -P0_x * tmp19 - P0_y * tmp21 - P0_z * tmp14;
  const REAL tmp39 = -S0_y * tmp14 + S0_z * tmp21;
  const REAL tmp42 = -S1_y * tmp3 + S1_z * tmp21;
  const REAL tmp13 = tmp11 + ((tmp3) * (tmp3)) + tmp9;
  const REAL tmp15 = tmp11 + ((tmp14) * (tmp14)) + tmp9;
  const REAL tmp17 = pow(tmp13, -5.0 / 2.0);
  const REAL tmp25 = pow(tmp15, -5.0 / 2.0);
  const REAL tmp46 = pow(tmp13, -3.0 / 2.0);
  const REAL tmp47 = pow(tmp15, -3.0 / 2.0);
  const REAL tmp24 = tmp17 * tmp23;
  const REAL tmp30 = tmp25 * tmp29;
  const REAL tmp48 = P0_z * tmp14 * tmp47 + P1_z * tmp3 * tmp46;
  const REAL tmp50 = -tmp14 * tmp25;
  const REAL tmp53 = -tmp17 * tmp3;
  const REAL tmp60 = -tmp14 * tmp47;
  const REAL tmp61 = -tmp3 * tmp46;
  const REAL tmp63 = P0_x * tmp19 * tmp47;
  const REAL tmp64 = P1_x * tmp19 * tmp46;
  const REAL tmp66 = P0_y * tmp21 * tmp47;
  const REAL tmp67 = P1_y * tmp21 * tmp46;
  const REAL tmp35 = tmp17 * tmp21 * tmp31 + tmp21 * tmp25 * tmp33;
  const REAL tmp44 = tmp17 * tmp19 * tmp42 + tmp19 * tmp25 * tmp39;
  const REAL tmp55 = tmp23 * tmp53;
  const REAL tmp58 = tmp29 * tmp50;
  const REAL tmp62 = P0_z * tmp60 + P1_z * tmp61;
  const REAL tmp36 = _Rational_1_4 * tmp11 * tmp24 + _Rational_1_4 * tmp11 * tmp30 + tmp35;
  const REAL tmp45 = _Rational_1_4 * tmp24 * tmp9 + _Rational_1_4 * tmp30 * tmp9 + tmp44;
  const REAL tmp71 = _Integer_3 * (-_Rational_1_4 * _Rational_2_3 * tmp14 * tmp58 - _Rational_1_4 * _Rational_2_3 * tmp3 * tmp55 -
                                   _Rational_2_3 * tmp49 * tmp50 - _Rational_2_3 * tmp52 * tmp53) +
                     _Rational_1_6 * tmp48 + _Rational_7_6 * tmp62;
  *psi_background = _Integer_1 + _Rational_1_2 * (bare_mass_0 / sqrt(tmp15) + bare_mass_1 / sqrt(tmp13));
  *ADD_times_AUU = _Integer_2 * (((_Integer_3 * (-tmp17 * tmp19 * tmp31 - tmp17 * tmp21 * tmp42 - tmp19 * tmp25 * tmp33 - tmp21 * tmp25 * tmp39) +
                                   _Rational_1_4 * (-tmp24 * tmp73 - tmp30 * tmp73) +
                                   _Rational_3_2 * (P0_x * tmp21 * tmp47 + P0_y * tmp19 * tmp47 + P1_x * tmp21 * tmp46 + P1_y * tmp19 * tmp46)) *
                                  (_Integer_3 * (-tmp17 * tmp19 * tmp31 - tmp17 * tmp21 * tmp42 - tmp19 * tmp25 * tmp33 - tmp21 * tmp25 * tmp39) +
                                   _Rational_1_4 * (-tmp24 * tmp73 - tmp30 * tmp73) +
                                   _Rational_3_2 * (P0_x * tmp21 * tmp47 + P0_y * tmp19 * tmp47 + P1_x * tmp21 * tmp46 + P1_y * tmp19 * tmp46))) +
                                 ((_Integer_3 * (-tmp14 * tmp30 * tmp82 - tmp17 * tmp21 * tmp52 - tmp21 * tmp25 * tmp49 - tmp24 * tmp3 * tmp82 +
                                                 tmp31 * tmp53 + tmp33 * tmp50 + tmp55 * tmp82 + tmp58 * tmp82) +
                                   _Rational_1_4 * (-P0_y * tmp14 * tmp47 - P1_y * tmp3 * tmp46) +
                                   _Rational_3_2 * (P0_z * tmp21 * tmp47 + P1_z * tmp21 * tmp46) + _Rational_7_4 * (-P0_y * tmp60 - P1_y * tmp61)) *
                                  (_Integer_3 * (-tmp14 * tmp30 * tmp82 - tmp17 * tmp21 * tmp52 - tmp21 * tmp25 * tmp49 - tmp24 * tmp3 * tmp82 +
                                                 tmp31 * tmp53 + tmp33 * tmp50 + tmp55 * tmp82 + tmp58 * tmp82) +
                                   _Rational_1_4 * (-P0_y * tmp14 * tmp47 - P1_y * tmp3 * tmp46) +
                                   _Rational_3_2 * (P0_z * tmp21 * tmp47 + P1_z * tmp21 * tmp46) + _Rational_7_4 * (-P0_y * tmp60 - P1_y * tmp61))) +
                                 ((_Integer_3 * (-tmp14 * tmp30 * tmp85 - tmp17 * tmp19 * tmp52 - tmp19 * tmp25 * tmp49 - tmp24 * tmp3 * tmp85 +
                                                 tmp39 * tmp50 + tmp42 * tmp53 + tmp55 * tmp85 + tmp58 * tmp85) +
                                   _Rational_1_4 * (-P0_x * tmp14 * tmp47 - P1_x * tmp3 * tmp46) +
                                   _Rational_3_2 * (P0_z * tmp19 * tmp47 + P1_z * tmp19 * tmp46) + _Rational_7_4 * (-P0_x * tmp60 - P1_x * tmp61)) *
                                  (_Integer_3 * (-tmp14 * tmp30 * tmp85 - tmp17 * tmp19 * tmp52 - tmp19 * tmp25 * tmp49 - tmp24 * tmp3 * tmp85 +
                                                 tmp39 * tmp50 + tmp42 * tmp53 + tmp55 * tmp85 + tmp58 * tmp85) +
                                   _Rational_1_4 * (-P0_x * tmp14 * tmp47 - P1_x * tmp3 * tmp46) +
                                   _Rational_3_2 * (P0_z * tmp19 * tmp47 + P1_z * tmp19 * tmp46) + _Rational_7_4 * (-P0_x * tmp60 - P1_x * tmp61)))) +
                   ((_Integer_2 * (tmp36 + tmp63 + tmp64) - _Integer_4 * tmp44 +
                     _Rational_1_4 * (-_Integer_4 * tmp24 * tmp9 - _Integer_4 * tmp30 * tmp9) - tmp66 - tmp67 + tmp71) *
                    (_Integer_2 * (tmp36 + tmp63 + tmp64) - _Integer_4 * tmp44 +
                     _Rational_1_4 * (-_Integer_4 * tmp24 * tmp9 - _Integer_4 * tmp30 * tmp9) - tmp66 - tmp67 + tmp71)) +
                   ((_Integer_2 * (tmp45 + tmp66 + tmp67) - _Integer_4 * tmp35 +
                     _Rational_1_4 * (-_Integer_4 * tmp11 * tmp24 - _Integer_4 * tmp11 * tmp30) - tmp63 - tmp64 + tmp71) *
                    (_Integer_2 * (tmp45 + tmp66 + tmp67) - _Integer_4 * tmp35 +
                     _Rational_1_4 * (-_Integer_4 * tmp11 * tmp24 - _Integer_4 * tmp11 * tmp30) - tmp63 - tmp64 + tmp71)) +
                   ((_Integer_2 * (tmp36 + tmp45) - _Rational_1_3 * tmp48 +
                     _Rational_4_3 * (_Integer_3 * _Rational_1_4 * tmp14 * tmp58 + _Integer_3 * _Rational_1_4 * tmp3 * tmp55 +
                                      _Integer_3 * tmp49 * tmp50 + _Integer_3 * tmp52 * tmp53) -
                     _Rational_7_3 * tmp62 - tmp63 - tmp64 - tmp66 - tmp67) *
                    (_Integer_2 * (tmp36 + tmp45) - _Rational_1_3 * tmp48 +
                     _Rational_4_3 * (_Integer_3 * _Rational_1_4 * tmp14 * tmp58 + _Integer_3 * _Rational_1_4 * tmp3 * tmp55 +
                                      _Integer_3 * tmp49 * tmp50 + _Integer_3 * tmp52 * tmp53) -
                     _Rational_7_3 * tmp62 - tmp63 - tmp64 - tmp66 - tmp67));
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
  const REAL AMPLRHO = d_params[streamid].AMPLRHO;
  const REAL AMPLZ = d_params[streamid].AMPLZ;
  const REAL SINHWRHO = d_params[streamid].SINHWRHO;
  const REAL SINHWZ = d_params[streamid].SINHWZ;
  const REAL rho_slope = d_params[streamid].rho_slope;
  const REAL z_slope = d_params[streamid].z_slope;

  for (int i2 = tid2 + NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2 += stride2) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = tid1 + NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1 += stride1) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = tid0 + NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += stride0) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        /*
         *  Original SymPy expressions:
         *  "[const REAL dsmin0 = d_params[streamid].dxx0*(rho_slope + xx0**2*(AMPLRHO - rho_slope)*(exp(xx0/SINHWRHO)/SINHWRHO +
         * exp(-xx0/SINHWRHO)/SINHWRHO)/(exp(1/SINHWRHO) - exp(-1/SINHWRHO)) + 2*xx0*(AMPLRHO - rho_slope)*(exp(xx0/SINHWRHO) -
         * exp(-xx0/SINHWRHO))/(exp(1/SINHWRHO) - exp(-1/SINHWRHO)))]"
         *  "[const REAL dsmin1 = d_params[streamid].dxx1*(rho_slope*xx0 + xx0**2*(AMPLRHO - rho_slope)*(exp(xx0/SINHWRHO) -
         * exp(-xx0/SINHWRHO))/(exp(1/SINHWRHO) - exp(-1/SINHWRHO)))]"
         *  "[const REAL dsmin2 = d_params[streamid].dxx2*(xx2**2*(AMPLZ - z_slope)*(exp(xx2/SINHWZ)/SINHWZ + exp(-xx2/SINHWZ)/SINHWZ)/(exp(1/SINHWZ)
         * - exp(-1/SINHWZ)) + 2*xx2*(AMPLZ - z_slope)*(exp(xx2/SINHWZ) - exp(-xx2/SINHWZ))/(exp(1/SINHWZ) - exp(-1/SINHWZ)) + z_slope)]"
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
        const REAL dsmin0 = d_params[streamid].dxx0 * (rho_slope + 2 * tmp4 * tmp5 * xx0 + tmp6 * (tmp0 * tmp2 + tmp0 * tmp3));
        const REAL dsmin1 = d_params[streamid].dxx1 * (rho_slope * xx0 + tmp4 * tmp6);
        const REAL dsmin2 =
            d_params[streamid].dxx2 * (tmp11 * ((xx2) * (xx2)) * (tmp10 * tmp7 + tmp7 * tmp9) + 2 * tmp11 * xx2 * (-tmp10 + tmp9) + z_slope);

        // Set local wavespeed
        in_gfs[IDX4(VARIABLE_WAVESPEEDGF, i0, i1, i2)] = MINIMUM_GLOBAL_WAVESPEED * NRPYMIN(dsmin0, NRPYMIN(dsmin1, dsmin2)) / dt;

      } // END LOOP: for i0 over [tid0+NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS)
    } // END LOOP: for i1 over [tid1+NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS)
  } // END LOOP: for i2 over [tid2+NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
} // END FUNCTION: variable_wavespeed_gfs_all_points_gpu

/**
 * Call functions that set up all AUXEVOL gridfunctions.
 */
void auxevol_gfs_set_to_constant__rfm__SinhCylindricalv2n2(commondata_struct *restrict commondata, params_struct *restrict params,
                                                           REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs) {
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

} // END FUNCTION: auxevol_gfs_set_to_constant__rfm__SinhCylindricalv2n2
