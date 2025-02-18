#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0_gpu.
 * Kernel to precompute metric quantity f0_of_xx0.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLRHO = d_params[streamid].AMPLRHO;
  const REAL SINHWRHO = d_params[streamid].SINHWRHO;
  const REAL rho_slope = d_params[streamid].rho_slope;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0[i0] = rho_slope * xx0 + powf(xx0, 2) * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) - expf(-xx0 / SINHWRHO)) /
                                                     (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0_gpu
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__D0_gpu.
 * Kernel to precompute metric quantity f0_of_xx0__D0.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0__D0_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLRHO = d_params[streamid].AMPLRHO;
  const REAL SINHWRHO = d_params[streamid].SINHWRHO;
  const REAL rho_slope = d_params[streamid].rho_slope;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__D0[i0] =
        rho_slope +
        powf(xx0, 2) * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / SINHWRHO + expf(-xx0 / SINHWRHO) / SINHWRHO) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO)) +
        2 * xx0 * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) - expf(-xx0 / SINHWRHO)) / (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__D0_gpu
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DD00_gpu.
 * Kernel to precompute metric quantity f0_of_xx0__DD00.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0__DD00_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLRHO = d_params[streamid].AMPLRHO;
  const REAL SINHWRHO = d_params[streamid].SINHWRHO;
  const REAL rho_slope = d_params[streamid].rho_slope;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DD00[i0] =
        powf(xx0, 2) * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / powf(SINHWRHO, 2) - expf(-xx0 / SINHWRHO) / powf(SINHWRHO, 2)) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO)) +
        4 * xx0 * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / SINHWRHO + expf(-xx0 / SINHWRHO) / SINHWRHO) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO)) +
        2 * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) - expf(-xx0 / SINHWRHO)) / (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DD00_gpu
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DDD000_gpu.
 * Kernel to precompute metric quantity f0_of_xx0__DDD000.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0__DDD000_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPLRHO = d_params[streamid].AMPLRHO;
  const REAL SINHWRHO = d_params[streamid].SINHWRHO;
  const REAL rho_slope = d_params[streamid].rho_slope;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DDD000[i0] =
        powf(xx0, 2) * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / powf(SINHWRHO, 3) + expf(-xx0 / SINHWRHO) / powf(SINHWRHO, 3)) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO)) +
        6 * xx0 * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / powf(SINHWRHO, 2) - expf(-xx0 / SINHWRHO) / powf(SINHWRHO, 2)) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO)) +
        6 * (AMPLRHO - rho_slope) * (expf(xx0 / SINHWRHO) / SINHWRHO + expf(-xx0 / SINHWRHO) / SINHWRHO) /
            (expf(1.0F / SINHWRHO) - expf(-1 / SINHWRHO));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DDD000_gpu
/**
 * Kernel: rfm_precompute_defines__f3_of_xx2_gpu.
 * Kernel to precompute metric quantity f3_of_xx2.
 */
__global__ static void rfm_precompute_defines__f3_of_xx2_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x2) {
  // Temporary parameters
  const REAL AMPLZ = d_params[streamid].AMPLZ;
  const REAL SINHWZ = d_params[streamid].SINHWZ;
  const REAL z_slope = d_params[streamid].z_slope;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i2 = tid0; i2 < d_params[streamid].Nxx_plus_2NGHOSTS2; i2 += stride0) {
    const REAL xx2 = x2[i2];
    rfmstruct->f3_of_xx2[i2] =
        powf(xx2, 2) * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / SINHWZ + expf(-xx2 / SINHWZ) / SINHWZ) / (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) +
        2 * xx2 * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) - expf(-xx2 / SINHWZ)) / (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) + z_slope;
  }
} // END FUNCTION rfm_precompute_defines__f3_of_xx2_gpu
/**
 * Kernel: rfm_precompute_defines__f3_of_xx2__D2_gpu.
 * Kernel to precompute metric quantity f3_of_xx2__D2.
 */
__global__ static void rfm_precompute_defines__f3_of_xx2__D2_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x2) {
  // Temporary parameters
  const REAL AMPLZ = d_params[streamid].AMPLZ;
  const REAL SINHWZ = d_params[streamid].SINHWZ;
  const REAL z_slope = d_params[streamid].z_slope;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i2 = tid0; i2 < d_params[streamid].Nxx_plus_2NGHOSTS2; i2 += stride0) {
    const REAL xx2 = x2[i2];
    rfmstruct->f3_of_xx2__D2[i2] =
        powf(xx2, 2) * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / powf(SINHWZ, 2) - expf(-xx2 / SINHWZ) / powf(SINHWZ, 2)) /
            (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) +
        4 * xx2 * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / SINHWZ + expf(-xx2 / SINHWZ) / SINHWZ) / (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) +
        2 * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) - expf(-xx2 / SINHWZ)) / (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ));
  }
} // END FUNCTION rfm_precompute_defines__f3_of_xx2__D2_gpu
/**
 * Kernel: rfm_precompute_defines__f3_of_xx2__DD22_gpu.
 * Kernel to precompute metric quantity f3_of_xx2__DD22.
 */
__global__ static void rfm_precompute_defines__f3_of_xx2__DD22_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x2) {
  // Temporary parameters
  const REAL AMPLZ = d_params[streamid].AMPLZ;
  const REAL SINHWZ = d_params[streamid].SINHWZ;
  const REAL z_slope = d_params[streamid].z_slope;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i2 = tid0; i2 < d_params[streamid].Nxx_plus_2NGHOSTS2; i2 += stride0) {
    const REAL xx2 = x2[i2];
    rfmstruct->f3_of_xx2__DD22[i2] =
        powf(xx2, 2) * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / powf(SINHWZ, 3) + expf(-xx2 / SINHWZ) / powf(SINHWZ, 3)) /
            (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) +
        6 * xx2 * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / powf(SINHWZ, 2) - expf(-xx2 / SINHWZ) / powf(SINHWZ, 2)) /
            (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ)) +
        6 * (AMPLZ - z_slope) * (expf(xx2 / SINHWZ) / SINHWZ + expf(-xx2 / SINHWZ) / SINHWZ) / (expf(1.0F / SINHWZ) - expf(-1 / SINHWZ));
  }
} // END FUNCTION rfm_precompute_defines__f3_of_xx2__DD22_gpu

/**
 * rfm_precompute_defines: reference metric precomputed lookup arrays: defines
 */
void rfm_precompute_defines__rfm__SinhCylindricalv2n2(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                      rfm_struct *restrict rfmstruct, REAL *restrict xx[3]) {
  MAYBE_UNUSED const REAL *restrict x0 = xx[0];
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const REAL *restrict x1 = xx[1];
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const REAL *restrict x2 = xx[2];
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                         (Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                         (Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__f0_of_xx0_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x0);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f0_of_xx0_gpu failure");
  }
  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                         (Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                         (Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__f0_of_xx0__D0_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x0);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f0_of_xx0__D0_gpu failure");
  }
  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                         (Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                         (Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__f0_of_xx0__DD00_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x0);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f0_of_xx0__DD00_gpu failure");
  }
  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                         (Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                         (Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__f0_of_xx0__DDD000_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x0);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f0_of_xx0__DDD000_gpu failure");
  }
  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                         (Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                         (Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__f3_of_xx2_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x2);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f3_of_xx2_gpu failure");
  }
  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                         (Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                         (Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__f3_of_xx2__D2_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x2);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f3_of_xx2__D2_gpu failure");
  }
  {

    const size_t threads_in_x_dir = 32;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                         (Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                         (Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_defines__f3_of_xx2__DD22_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x2);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f3_of_xx2__DD22_gpu failure");
  }
} // END FUNCTION rfm_precompute_defines__rfm__SinhCylindricalv2n2
