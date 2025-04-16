#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"

/**
 * Kernel: auxevol_gfs_single_point_gpu.
 * Kernel to compute AUXEVOL grid functions at a single point.
 */
__device__ static void auxevol_gfs_single_point_gpu(const size_t streamid, const REAL xx0, const REAL xx1, const REAL xx2,
                                                    REAL *restrict psi_background, REAL *restrict ADD_times_AUU) {
  // Load necessary parameters from params_struct

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
} // END FUNCTION auxevol_gfs_single_point_gpu

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
        cudaCheckErrors(cudaKernel, "auxevol_gfs_single_point_gpu failure");

      } // END LOOP: for (int i0 = tid0; i0 < Nxx_plus_2NGHOSTS0; i0 += stride0)
    } // END LOOP: for (int i1 = tid1; i1 < Nxx_plus_2NGHOSTS1; i1 += stride1)
  } // END LOOP: for (int i2 = tid2; i2 < Nxx_plus_2NGHOSTS2; i2 += stride2)
} // END FUNCTION auxevol_gfs_all_points_gpu

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

  for (int i2 = tid2 + NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2 += stride2) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = tid1 + NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1 += stride1) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = tid0 + NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += stride0) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        /*
         *  Original SymPy expressions:
         *  "[const REAL dsmin0 = dxx0]"
         *  "[const REAL dsmin1 = dxx1]"
         *  "[const REAL dsmin2 = dxx2]"
         */
        const REAL dsmin0 = dxx0;
        const REAL dsmin1 = dxx1;
        const REAL dsmin2 = dxx2;

        // Set local wavespeed
        in_gfs[IDX4(VARIABLE_WAVESPEEDGF, i0, i1, i2)] = MINIMUM_GLOBAL_WAVESPEED * MIN(dsmin0, MIN(dsmin1, dsmin2)) / dt;

      } // END LOOP: for (int i0 = tid0+NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += stride0)
    } // END LOOP: for (int i1 = tid1+NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1 += stride1)
  } // END LOOP: for (int i2 = tid2+NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2 += stride2)
} // END FUNCTION variable_wavespeed_gfs_all_points_gpu

/**
 * Call functions that set up all AUXEVOL gridfunctions.
 */
void initialize_constant_auxevol__rfm__Cartesian(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                                                 MoL_gridfunctions_struct *restrict gridfuncs) {
#include "../set_CodeParameters.h"
  cpyHosttoDevice_commondata__constant(commondata);

  REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
  REAL *restrict x0 = xx[0];
  REAL *restrict x1 = xx[1];
  REAL *restrict x2 = xx[2];

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
} // END FUNCTION initialize_constant_auxevol__rfm__Cartesian
