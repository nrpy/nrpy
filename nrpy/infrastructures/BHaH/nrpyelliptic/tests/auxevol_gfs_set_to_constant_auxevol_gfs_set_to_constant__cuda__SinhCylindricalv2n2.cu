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

  const REAL tmp1 = (1.0 / (SINHWZ));
  const REAL tmp6 = cos(xx1);
  const REAL tmp8 = (1.0 / (SINHWRHO));
  const REAL tmp13 = sin(xx1);
  const REAL tmp7 = ((tmp6) * (tmp6));
  const REAL tmp3 = ((xx2) * (xx2)) * (AMPLZ - z_slope) * (exp(tmp1 * xx2) - exp(-tmp1 * xx2)) / (exp(tmp1) - exp(-tmp1));
  const REAL tmp10 = rho_slope * xx0 + ((xx0) * (xx0)) * (AMPLRHO - rho_slope) * (exp(tmp8 * xx0) - exp(-tmp8 * xx0)) / (exp(tmp8) - exp(-tmp8));
  const REAL tmp4 = tmp3 + xx2 * z_slope;
  const REAL tmp11 = ((tmp10) * (tmp10));
  const REAL tmp21 = tmp10 * tmp6;
  const REAL tmp26 = tmp10 * tmp13;
  const REAL tmp56 = 3 * tmp3 + 3 * xx2 * z_slope;
  const REAL tmp5 = tmp4 + zPunc;
  const REAL tmp12 = tmp11 * tmp7;
  const REAL tmp14 = tmp11 * ((tmp13) * (tmp13));
  const REAL tmp17 = tmp4 - zPunc;
  const REAL tmp42 = 6 * tmp11 * tmp13 * tmp6;
  const REAL tmp54 = -S1_x * tmp26 + S1_y * tmp10 * tmp6;
  const REAL tmp57 = -tmp56 - 3 * zPunc;
  const REAL tmp58 = -S0_x * tmp26 + S0_y * tmp10 * tmp6;
  const REAL tmp59 = -tmp56 + 3 * zPunc;
  const REAL tmp22 = S1_x * tmp5 - S1_z * tmp21;
  const REAL tmp35 = S0_x * tmp17 - S0_z * tmp21;
  const REAL tmp16 = tmp12 + tmp14 + ((tmp5) * (tmp5));
  const REAL tmp18 = tmp12 + tmp14 + ((tmp17) * (tmp17));
  const REAL tmp28 = -S0_y * tmp17 + S0_z * tmp26;
  const REAL tmp33 = -S1_y * tmp5 + S1_z * tmp26;
  const REAL tmp40 = -1.0 / 4.0 * P1_x * tmp21 - 1.0 / 4.0 * P1_y * tmp26 - 1.0 / 4.0 * P1_z * tmp5;
  const REAL tmp46 = -1.0 / 4.0 * P0_x * tmp21 - 1.0 / 4.0 * P0_y * tmp26 - 1.0 / 4.0 * P0_z * tmp17;
  const REAL tmp19 = pow(tmp16, -3.0 / 2.0);
  const REAL tmp20 = pow(tmp18, -3.0 / 2.0);
  const REAL tmp23 = pow(tmp16, -5.0 / 2.0);
  const REAL tmp29 = pow(tmp18, -5.0 / 2.0);
  const REAL tmp41 = tmp23 * tmp40;
  const REAL tmp47 = tmp29 * tmp46;
  const REAL tmp49 = -tmp17 * tmp20;
  const REAL tmp52 = -tmp19 * tmp5;
  const REAL tmp66 = P0_z * tmp17 * tmp20;
  const REAL tmp67 = P1_z * tmp19 * tmp5;
  const REAL tmp72 = P0_x * tmp20 * tmp21;
  const REAL tmp73 = P1_x * tmp19 * tmp21;
  const REAL tmp25 = 3 * tmp21 * tmp23;
  const REAL tmp31 = 3 * tmp26 * tmp29;
  const REAL tmp34 = 3 * tmp23 * tmp26;
  const REAL tmp36 = 3 * tmp21 * tmp29;
  const REAL tmp69 = tmp23 * tmp54 * tmp57;
  const REAL tmp71 = tmp29 * tmp58 * tmp59;
  const REAL tmp79 = tmp23 * tmp40 * tmp5 * tmp57;
  const REAL tmp80 = tmp17 * tmp29 * tmp46 * tmp59;
  const REAL tmp84 =
      -P0_y * tmp20 * tmp26 - P1_y * tmp19 * tmp26 + 2 * tmp14 * tmp41 + 2 * tmp14 * tmp47 + 2 * tmp22 * tmp23 * tmp26 + 2 * tmp26 * tmp29 * tmp35;
  *psi_background = (1.0 / 2.0) * bare_mass_0 / sqrt(tmp18) + (1.0 / 2.0) * bare_mass_1 / sqrt(tmp16) + 1;
  *ADD_times_AUU = 2 * (((3.0 / 2.0) * P0_x * tmp10 * tmp13 * tmp20 + (3.0 / 2.0) * P0_y * tmp10 * tmp20 * tmp6 +
                         (3.0 / 2.0) * P1_x * tmp10 * tmp13 * tmp19 + (3.0 / 2.0) * P1_y * tmp10 * tmp19 * tmp6 - tmp22 * tmp25 - tmp28 * tmp31 -
                         tmp33 * tmp34 - tmp35 * tmp36 - tmp41 * tmp42 - tmp42 * tmp47) *
                        ((3.0 / 2.0) * P0_x * tmp10 * tmp13 * tmp20 + (3.0 / 2.0) * P0_y * tmp10 * tmp20 * tmp6 +
                         (3.0 / 2.0) * P1_x * tmp10 * tmp13 * tmp19 + (3.0 / 2.0) * P1_y * tmp10 * tmp19 * tmp6 - tmp22 * tmp25 - tmp28 * tmp31 -
                         tmp33 * tmp34 - tmp35 * tmp36 - tmp41 * tmp42 - tmp42 * tmp47)) +
                   2 * ((-1.0 / 4.0 * P0_x * tmp17 * tmp20 - 7.0 / 4.0 * P0_x * tmp49 + (3.0 / 2.0) * P0_z * tmp10 * tmp20 * tmp6 -
                         1.0 / 4.0 * P1_x * tmp19 * tmp5 - 7.0 / 4.0 * P1_x * tmp52 + (3.0 / 2.0) * P1_z * tmp10 * tmp19 * tmp6 +
                         tmp10 * tmp23 * tmp40 * tmp57 * tmp6 + tmp10 * tmp29 * tmp46 * tmp59 * tmp6 - 3 * tmp17 * tmp21 * tmp47 -
                         3 * tmp21 * tmp41 * tmp5 + tmp23 * tmp33 * tmp57 - tmp25 * tmp54 + tmp28 * tmp29 * tmp59 - tmp36 * tmp58) *
                        (-1.0 / 4.0 * P0_x * tmp17 * tmp20 - 7.0 / 4.0 * P0_x * tmp49 + (3.0 / 2.0) * P0_z * tmp10 * tmp20 * tmp6 -
                         1.0 / 4.0 * P1_x * tmp19 * tmp5 - 7.0 / 4.0 * P1_x * tmp52 + (3.0 / 2.0) * P1_z * tmp10 * tmp19 * tmp6 +
                         tmp10 * tmp23 * tmp40 * tmp57 * tmp6 + tmp10 * tmp29 * tmp46 * tmp59 * tmp6 - 3 * tmp17 * tmp21 * tmp47 -
                         3 * tmp21 * tmp41 * tmp5 + tmp23 * tmp33 * tmp57 - tmp25 * tmp54 + tmp28 * tmp29 * tmp59 - tmp36 * tmp58)) +
                   2 * ((-1.0 / 4.0 * P0_y * tmp17 * tmp20 - 7.0 / 4.0 * P0_y * tmp49 + (3.0 / 2.0) * P0_z * tmp10 * tmp13 * tmp20 -
                         1.0 / 4.0 * P1_y * tmp19 * tmp5 - 7.0 / 4.0 * P1_y * tmp52 + (3.0 / 2.0) * P1_z * tmp10 * tmp13 * tmp19 +
                         tmp10 * tmp13 * tmp23 * tmp40 * tmp57 + tmp10 * tmp13 * tmp29 * tmp46 * tmp59 - 3 * tmp17 * tmp26 * tmp47 +
                         tmp22 * tmp23 * tmp57 - 3 * tmp26 * tmp41 * tmp5 + tmp29 * tmp35 * tmp59 - tmp31 * tmp58 - tmp34 * tmp54) *
                        (-1.0 / 4.0 * P0_y * tmp17 * tmp20 - 7.0 / 4.0 * P0_y * tmp49 + (3.0 / 2.0) * P0_z * tmp10 * tmp13 * tmp20 -
                         1.0 / 4.0 * P1_y * tmp19 * tmp5 - 7.0 / 4.0 * P1_y * tmp52 + (3.0 / 2.0) * P1_z * tmp10 * tmp13 * tmp19 +
                         tmp10 * tmp13 * tmp23 * tmp40 * tmp57 + tmp10 * tmp13 * tmp29 * tmp46 * tmp59 - 3 * tmp17 * tmp26 * tmp47 +
                         tmp22 * tmp23 * tmp57 - 3 * tmp26 * tmp41 * tmp5 + tmp29 * tmp35 * tmp59 - tmp31 * tmp58 - tmp34 * tmp54)) +
                   ((-7.0 / 3.0 * P0_z * tmp49 - 7.0 / 3.0 * P1_z * tmp52 + 2 * tmp12 * tmp41 + 2 * tmp12 * tmp47 + 2 * tmp21 * tmp23 * tmp33 +
                     2 * tmp21 * tmp28 * tmp29 - 1.0 / 3.0 * tmp66 - 1.0 / 3.0 * tmp67 + (4.0 / 3.0) * tmp69 + (4.0 / 3.0) * tmp71 - tmp72 - tmp73 +
                     (4.0 / 3.0) * tmp79 + (4.0 / 3.0) * tmp80 + tmp84) *
                    (-7.0 / 3.0 * P0_z * tmp49 - 7.0 / 3.0 * P1_z * tmp52 + 2 * tmp12 * tmp41 + 2 * tmp12 * tmp47 + 2 * tmp21 * tmp23 * tmp33 +
                     2 * tmp21 * tmp28 * tmp29 - 1.0 / 3.0 * tmp66 - 1.0 / 3.0 * tmp67 + (4.0 / 3.0) * tmp69 + (4.0 / 3.0) * tmp71 - tmp72 - tmp73 +
                     (4.0 / 3.0) * tmp79 + (4.0 / 3.0) * tmp80 + tmp84)) +
                   (((7.0 / 6.0) * P0_z * tmp49 + (7.0 / 6.0) * P1_z * tmp52 - 4 * tmp12 * tmp41 - 4 * tmp12 * tmp47 - 4 * tmp21 * tmp23 * tmp33 -
                     4 * tmp21 * tmp28 * tmp29 + (1.0 / 6.0) * tmp66 + (1.0 / 6.0) * tmp67 - 2.0 / 3.0 * tmp69 - 2.0 / 3.0 * tmp71 + 2 * tmp72 +
                     2 * tmp73 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp80 + tmp84) *
                    ((7.0 / 6.0) * P0_z * tmp49 + (7.0 / 6.0) * P1_z * tmp52 - 4 * tmp12 * tmp41 - 4 * tmp12 * tmp47 - 4 * tmp21 * tmp23 * tmp33 -
                     4 * tmp21 * tmp28 * tmp29 + (1.0 / 6.0) * tmp66 + (1.0 / 6.0) * tmp67 - 2.0 / 3.0 * tmp69 - 2.0 / 3.0 * tmp71 + 2 * tmp72 +
                     2 * tmp73 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp80 + tmp84)) +
                   ((2 * P0_y * tmp10 * tmp13 * tmp20 - P0_z * tmp17 * tmp20 + 2 * P1_y * tmp10 * tmp13 * tmp19 - P1_z * tmp19 * tmp5 +
                     2 * tmp10 * tmp23 * tmp33 * tmp6 + 2 * tmp10 * tmp28 * tmp29 * tmp6 + 2 * tmp11 * tmp23 * tmp40 * tmp7 +
                     2 * tmp11 * tmp29 * tmp46 * tmp7 - 4 * tmp14 * tmp41 - 4 * tmp14 * tmp47 - 4 * tmp22 * tmp23 * tmp26 -
                     4 * tmp26 * tmp29 * tmp35 - 2.0 / 3.0 * tmp69 - 2.0 / 3.0 * tmp71 - tmp72 - tmp73 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp80) *
                    (2 * P0_y * tmp10 * tmp13 * tmp20 - P0_z * tmp17 * tmp20 + 2 * P1_y * tmp10 * tmp13 * tmp19 - P1_z * tmp19 * tmp5 +
                     2 * tmp10 * tmp23 * tmp33 * tmp6 + 2 * tmp10 * tmp28 * tmp29 * tmp6 + 2 * tmp11 * tmp23 * tmp40 * tmp7 +
                     2 * tmp11 * tmp29 * tmp46 * tmp7 - 4 * tmp14 * tmp41 - 4 * tmp14 * tmp47 - 4 * tmp22 * tmp23 * tmp26 -
                     4 * tmp26 * tmp29 * tmp35 - 2.0 / 3.0 * tmp69 - 2.0 / 3.0 * tmp71 - tmp72 - tmp73 - 2.0 / 3.0 * tmp79 - 2.0 / 3.0 * tmp80));
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
