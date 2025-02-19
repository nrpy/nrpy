#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0_gpu.
 * Kernel to precompute metric quantity f0_of_xx0.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPL = d_params[streamid].AMPL;
  const REAL SINHW = d_params[streamid].SINHW;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0[i0] = AMPL * (expf(xx0 / SINHW) - expf(-xx0 / SINHW)) / (expf(1.0F / SINHW) - expf(-1 / SINHW));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0_gpu
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__D0_gpu.
 * Kernel to precompute metric quantity f0_of_xx0__D0.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0__D0_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPL = d_params[streamid].AMPL;
  const REAL SINHW = d_params[streamid].SINHW;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__D0[i0] = AMPL * (expf(xx0 / SINHW) / SINHW + expf(-xx0 / SINHW) / SINHW) / (expf(1.0F / SINHW) - expf(-1 / SINHW));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__D0_gpu
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DD00_gpu.
 * Kernel to precompute metric quantity f0_of_xx0__DD00.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0__DD00_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPL = d_params[streamid].AMPL;
  const REAL SINHW = d_params[streamid].SINHW;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DD00[i0] =
        AMPL * (expf(xx0 / SINHW) / powf(SINHW, 2) - expf(-xx0 / SINHW) / powf(SINHW, 2)) / (expf(1.0F / SINHW) - expf(-1 / SINHW));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DD00_gpu
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DDD000_gpu.
 * Kernel to precompute metric quantity f0_of_xx0__DDD000.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0__DDD000_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMPL = d_params[streamid].AMPL;
  const REAL SINHW = d_params[streamid].SINHW;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DDD000[i0] =
        AMPL * (expf(xx0 / SINHW) / powf(SINHW, 3) + expf(-xx0 / SINHW) / powf(SINHW, 3)) / (expf(1.0F / SINHW) - expf(-1 / SINHW));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DDD000_gpu
/**
 * Kernel: rfm_precompute_defines__f1_of_xx1_gpu.
 * Kernel to precompute metric quantity f1_of_xx1.
 */
__global__ static void rfm_precompute_defines__f1_of_xx1_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x1) {
  // Temporary parameters
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i1 = tid0; i1 < d_params[streamid].Nxx_plus_2NGHOSTS1; i1 += stride0) {
    const REAL xx1 = x1[i1];
    rfmstruct->f1_of_xx1[i1] = sinf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f1_of_xx1_gpu
/**
 * Kernel: rfm_precompute_defines__f1_of_xx1__D1_gpu.
 * Kernel to precompute metric quantity f1_of_xx1__D1.
 */
__global__ static void rfm_precompute_defines__f1_of_xx1__D1_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x1) {
  // Temporary parameters
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i1 = tid0; i1 < d_params[streamid].Nxx_plus_2NGHOSTS1; i1 += stride0) {
    const REAL xx1 = x1[i1];
    rfmstruct->f1_of_xx1__D1[i1] = cosf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f1_of_xx1__D1_gpu
/**
 * Kernel: rfm_precompute_defines__f1_of_xx1__DD11_gpu.
 * Kernel to precompute metric quantity f1_of_xx1__DD11.
 */
__global__ static void rfm_precompute_defines__f1_of_xx1__DD11_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x1) {
  // Temporary parameters
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i1 = tid0; i1 < d_params[streamid].Nxx_plus_2NGHOSTS1; i1 += stride0) {
    const REAL xx1 = x1[i1];
    rfmstruct->f1_of_xx1__DD11[i1] = -sinf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f1_of_xx1__DD11_gpu

/**
 * rfm_precompute_defines: reference metric precomputed lookup arrays: defines
 */
void rfm_precompute_defines__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
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
    rfm_precompute_defines__f1_of_xx1_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x1);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f1_of_xx1_gpu failure");
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
    rfm_precompute_defines__f1_of_xx1__D1_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x1);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f1_of_xx1__D1_gpu failure");
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
    rfm_precompute_defines__f1_of_xx1__DD11_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x1);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f1_of_xx1__DD11_gpu failure");
  }
} // END FUNCTION rfm_precompute_defines__rfm__SinhSpherical
