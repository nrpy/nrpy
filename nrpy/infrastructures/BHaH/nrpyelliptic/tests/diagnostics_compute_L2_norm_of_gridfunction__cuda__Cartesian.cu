#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 * Kernel: compute_L2_norm_of_gridfunction_gpu.
 * Kernel to compute L2 quantities pointwise (not summed).
 */
__global__ static void compute_L2_norm_of_gridfunction_gpu(const size_t streamid, const REAL *restrict x0, const REAL *restrict x1,
                                                           const REAL *restrict x2, const REAL *restrict in_gfs, REAL *restrict aux_gfs,
                                                           const REAL integration_radius, const int gf_index) {
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
  const REAL dxx0 = d_params[streamid].dxx0;
  const REAL dxx1 = d_params[streamid].dxx1;
  const REAL dxx2 = d_params[streamid].dxx2;

  for (int i2 = tid2 + NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2 += stride2) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = tid1 + NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1 += stride1) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = tid0 + NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += stride0) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        /*
         *  Original SymPy expressions:
         *  "[const DOUBLE r = sqrt(xx0**2 + xx1**2 + xx2**2)]"
         *  "[const DOUBLE sqrtdetgamma = 1]"
         */
        const DOUBLE r = sqrt(((xx0) * (xx0)) + ((xx1) * (xx1)) + ((xx2) * (xx2)));
        const DOUBLE sqrtdetgamma = 1;

        if (r < integration_radius) {
          const DOUBLE gf_of_x = in_gfs[IDX4(gf_index, i0, i1, i2)];
          const DOUBLE dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;

          aux_gfs[IDX4(L2_SQUARED_DVGF, i0, i1, i2)] = gf_of_x * gf_of_x * dV;
          aux_gfs[IDX4(L2_DVGF, i0, i1, i2)] = dV;
        } // END if(r < integration_radius)

      } // END LOOP: for (int i0 = tid0+NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += stride0)
    } // END LOOP: for (int i1 = tid1+NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1 += stride1)
  } // END LOOP: for (int i2 = tid2+NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2 += stride2)
} // END FUNCTION compute_L2_norm_of_gridfunction_gpu

/**
 * Compute l2-norm of a gridfunction assuming a single grid.
 */
void compute_L2_norm_of_gridfunction__rfm__Cartesian(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                                                     const REAL integration_radius, const int gf_index, REAL *l2norm, const REAL *restrict in_gfs,
                                                     REAL *restrict aux_gfs) {
#include "../set_CodeParameters.h"

  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
  REAL *restrict x0 = xx[0];
  REAL *restrict x1 = xx[1];
  REAL *restrict x2 = xx[2];

  // Since we're performing sums, make sure arrays are zero'd
  cudaMemset(aux_gfs, 0, sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);

  const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_NELL_GRIDL2;
  const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_NELL_GRIDL2;
  const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_NELL_GRIDL2;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((params->Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                       (params->Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                       (params->Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  compute_L2_norm_of_gridfunction_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, x0, x1, x2, in_gfs, aux_gfs,
                                                                                                     integration_radius, gf_index);
  cudaCheckErrors(cudaKernel, "compute_L2_norm_of_gridfunction_gpu failure");

  // Set summation variables to compute l2-norm
  REAL squared_sum = find_global__sum(&aux_gfs[IDX4(L2_SQUARED_DVGF, 0, 0, 0)], Nxx_plus_2NGHOSTS_tot);
  REAL volume_sum = find_global__sum(&aux_gfs[IDX4(L2_DVGF, 0, 0, 0)], Nxx_plus_2NGHOSTS_tot);
  // Compute and output the log of the l2-norm.
  REAL local_norm = log10(1e-16 + sqrt(squared_sum / volume_sum)); // 1e-16 + ... avoids log10(0)

  // Compute and output the log of the l2-norm.
  *l2norm = log10(1e-16 + sqrt(squared_sum / volume_sum)); // 1e-16 + ... avoids log10(0)
} // END FUNCTION compute_L2_norm_of_gridfunction__rfm__Cartesian
