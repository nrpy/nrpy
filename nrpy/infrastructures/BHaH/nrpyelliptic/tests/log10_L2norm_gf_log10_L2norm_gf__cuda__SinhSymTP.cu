#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Kernel: log10_L2norm_gf_gpu.
 * Kernel to compute L2 quantities pointwise (not summed).
 */
__global__ static void log10_L2norm_gf_gpu(const size_t streamid, const REAL *restrict x0, const REAL *restrict x1, const REAL *restrict x2,
                                           const REAL *restrict in_gfs, REAL *restrict aux_gfs, const REAL integration_radius, const int gf_index) {
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
         *  "[const DOUBLE r = sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2*sin(xx1)**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 +
         * (AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2)*cos(xx1)**2)]"
         *  "[const DOUBLE sqrtdetgamma = AMAX**4*(exp(xx0/SINHWAA)/SINHWAA + exp(-xx0/SINHWAA)/SINHWAA)**2*(AMAX**2*(exp(xx0/SINHWAA) -
         * exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)**2*(exp(xx0/SINHWAA) -
         * exp(-xx0/SINHWAA))**2*sin(xx1)**2/((AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 +
         * bScale**2)*(exp(1/SINHWAA) - exp(-1/SINHWAA))**4)]"
         */
        const REAL tmp0 = ((sin(xx1)) * (sin(xx1)));
        const REAL tmp1 = (1.0 / (SINHWAA));
        const REAL tmp2 = exp(tmp1) - exp(-tmp1);
        const REAL tmp4 = exp(tmp1 * xx0);
        const REAL tmp5 = exp(-tmp1 * xx0);
        const REAL tmp6 = ((tmp4 - tmp5) * (tmp4 - tmp5));
        const REAL tmp7 = ((AMAX) * (AMAX)) * tmp6 / ((tmp2) * (tmp2));
        const REAL tmp9 = ((bScale) * (bScale)) + tmp7;
        const DOUBLE r = sqrt(tmp0 * tmp7 + tmp9 * ((cos(xx1)) * (cos(xx1))));
        const DOUBLE sqrtdetgamma = ((AMAX) * (AMAX) * (AMAX) * (AMAX)) * tmp0 * tmp6 *
                                    ((((bScale) * (bScale)) * tmp0 + tmp7) * (((bScale) * (bScale)) * tmp0 + tmp7)) *
                                    ((tmp1 * tmp4 + tmp1 * tmp5) * (tmp1 * tmp4 + tmp1 * tmp5)) / (((tmp2) * (tmp2) * (tmp2) * (tmp2)) * tmp9);

        if (r < integration_radius) {
          const DOUBLE gf_of_x = in_gfs[IDX4(gf_index, i0, i1, i2)];
          const DOUBLE dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;

          aux_gfs[IDX4(L2_SQUARED_DVGF, i0, i1, i2)] = gf_of_x * gf_of_x * dV;
          aux_gfs[IDX4(L2_DVGF, i0, i1, i2)] = dV;
        } // END if(r < integration_radius)

      } // END LOOP: for (int i0 = tid0+NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += stride0)
    } // END LOOP: for (int i1 = tid1+NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1 += stride1)
  } // END LOOP: for (int i2 = tid2+NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2 += stride2)
} // END FUNCTION log10_L2norm_gf_gpu

/**
 * Compute l2-norm of a gridfunction assuming a single grid.
 */
void log10_L2norm_gf__rfm__SinhSymTP(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                                     const REAL integration_radius, const int gf_index, REAL *l2norm, const REAL *restrict in_gfs,
                                     REAL *restrict aux_gfs) {
#include "set_CodeParameters.h"

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
  log10_L2norm_gf_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, x0, x1, x2, in_gfs, aux_gfs, integration_radius,
                                                                                     gf_index);
  cudaCheckErrors(cudaKernel, "log10_L2norm_gf_gpu failure");

  // Set summation variables to compute l2-norm
  REAL squared_sum = find_global__sum(&aux_gfs[IDX4(L2_SQUARED_DVGF, 0, 0, 0)], Nxx_plus_2NGHOSTS_tot);
  REAL volume_sum = find_global__sum(&aux_gfs[IDX4(L2_DVGF, 0, 0, 0)], Nxx_plus_2NGHOSTS_tot);
  // Compute and output the log of the l2-norm.
  REAL local_norm = log10(1e-16 + sqrt(squared_sum / volume_sum)); // 1e-16 + ... avoids log10(0)

  // Compute and output the log of the l2-norm.
  *l2norm = log10(1e-16 + sqrt(squared_sum / volume_sum)); // 1e-16 + ... avoids log10(0)
} // END FUNCTION log10_L2norm_gf__rfm__SinhSymTP
