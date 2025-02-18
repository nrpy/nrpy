#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0_gpu.
 * Kernel to precompute metric quantity f0_of_xx0.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = d_params[streamid].AMAX;
  const REAL SINHWAA = d_params[streamid].SINHWAA;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0[i0] = AMAX * (expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA)) / (expf(1.0F / SINHWAA) - expf(-1 / SINHWAA));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0_gpu
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__D0_gpu.
 * Kernel to precompute metric quantity f0_of_xx0__D0.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0__D0_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = d_params[streamid].AMAX;
  const REAL SINHWAA = d_params[streamid].SINHWAA;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__D0[i0] =
        AMAX * (expf(xx0 / SINHWAA) / SINHWAA + expf(-xx0 / SINHWAA) / SINHWAA) / (expf(1.0F / SINHWAA) - expf(-1 / SINHWAA));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__D0_gpu
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DD00_gpu.
 * Kernel to precompute metric quantity f0_of_xx0__DD00.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0__DD00_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = d_params[streamid].AMAX;
  const REAL SINHWAA = d_params[streamid].SINHWAA;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DD00[i0] =
        AMAX * (expf(xx0 / SINHWAA) / powf(SINHWAA, 2) - expf(-xx0 / SINHWAA) / powf(SINHWAA, 2)) / (expf(1.0F / SINHWAA) - expf(-1 / SINHWAA));
  }
} // END FUNCTION rfm_precompute_defines__f0_of_xx0__DD00_gpu
/**
 * Kernel: rfm_precompute_defines__f0_of_xx0__DDD000_gpu.
 * Kernel to precompute metric quantity f0_of_xx0__DDD000.
 */
__global__ static void rfm_precompute_defines__f0_of_xx0__DDD000_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = d_params[streamid].AMAX;
  const REAL SINHWAA = d_params[streamid].SINHWAA;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f0_of_xx0__DDD000[i0] =
        AMAX * (expf(xx0 / SINHWAA) / powf(SINHWAA, 3) + expf(-xx0 / SINHWAA) / powf(SINHWAA, 3)) / (expf(1.0F / SINHWAA) - expf(-1 / SINHWAA));
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
 * Kernel: rfm_precompute_defines__f2_of_xx0_gpu.
 * Kernel to precompute metric quantity f2_of_xx0.
 */
__global__ static void rfm_precompute_defines__f2_of_xx0_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = d_params[streamid].AMAX;
  const REAL SINHWAA = d_params[streamid].SINHWAA;
  const REAL bScale = d_params[streamid].bScale;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f2_of_xx0[i0] = sqrtf(
        powf(AMAX, 2) * powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) / powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2) + powf(bScale, 2));
  }
} // END FUNCTION rfm_precompute_defines__f2_of_xx0_gpu
/**
 * Kernel: rfm_precompute_defines__f2_of_xx0__D0_gpu.
 * Kernel to precompute metric quantity f2_of_xx0__D0.
 */
__global__ static void rfm_precompute_defines__f2_of_xx0__D0_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = d_params[streamid].AMAX;
  const REAL SINHWAA = d_params[streamid].SINHWAA;
  const REAL bScale = d_params[streamid].bScale;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f2_of_xx0__D0[i0] =
        (1.0F / 2.0F) * powf(AMAX, 2) * (2 * expf(xx0 / SINHWAA) / SINHWAA + 2 * expf(-xx0 / SINHWAA) / SINHWAA) *
        (expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA)) /
        (sqrtf(powf(AMAX, 2) * powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) / powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2) +
               powf(bScale, 2)) *
         powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2));
  }
} // END FUNCTION rfm_precompute_defines__f2_of_xx0__D0_gpu
/**
 * Kernel: rfm_precompute_defines__f2_of_xx0__DD00_gpu.
 * Kernel to precompute metric quantity f2_of_xx0__DD00.
 */
__global__ static void rfm_precompute_defines__f2_of_xx0__DD00_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x0) {
  // Temporary parameters
  const REAL AMAX = d_params[streamid].AMAX;
  const REAL SINHWAA = d_params[streamid].SINHWAA;
  const REAL bScale = d_params[streamid].bScale;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i0 = tid0; i0 < d_params[streamid].Nxx_plus_2NGHOSTS0; i0 += stride0) {
    const REAL xx0 = x0[i0];
    rfmstruct->f2_of_xx0__DD00[i0] =
        -1.0F / 4.0F * powf(AMAX, 4) * powf(2 * expf(xx0 / SINHWAA) / SINHWAA + 2 * expf(-xx0 / SINHWAA) / SINHWAA, 2) *
            powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) /
            (powf(powf(AMAX, 2) * powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) / powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2) +
                      powf(bScale, 2),
                  3.0F / 2.0F) *
             powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 4)) +
        (1.0F / 2.0F) * powf(AMAX, 2) * (2 * expf(xx0 / SINHWAA) / powf(SINHWAA, 2) - 2 * expf(-xx0 / SINHWAA) / powf(SINHWAA, 2)) *
            (expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA)) /
            (sqrtf(powf(AMAX, 2) * powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) / powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2) +
                   powf(bScale, 2)) *
             powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2)) +
        (1.0F / 2.0F) * powf(AMAX, 2) * (expf(xx0 / SINHWAA) / SINHWAA + expf(-xx0 / SINHWAA) / SINHWAA) *
            (2 * expf(xx0 / SINHWAA) / SINHWAA + 2 * expf(-xx0 / SINHWAA) / SINHWAA) /
            (sqrtf(powf(AMAX, 2) * powf(expf(xx0 / SINHWAA) - expf(-xx0 / SINHWAA), 2) / powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2) +
                   powf(bScale, 2)) *
             powf(expf(1.0F / SINHWAA) - expf(-1 / SINHWAA), 2));
  }
} // END FUNCTION rfm_precompute_defines__f2_of_xx0__DD00_gpu
/**
 * Kernel: rfm_precompute_defines__f4_of_xx1_gpu.
 * Kernel to precompute metric quantity f4_of_xx1.
 */
__global__ static void rfm_precompute_defines__f4_of_xx1_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x1) {
  // Temporary parameters
  const REAL bScale = d_params[streamid].bScale;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i1 = tid0; i1 < d_params[streamid].Nxx_plus_2NGHOSTS1; i1 += stride0) {
    const REAL xx1 = x1[i1];
    rfmstruct->f4_of_xx1[i1] = bScale * sinf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f4_of_xx1_gpu
/**
 * Kernel: rfm_precompute_defines__f4_of_xx1__D1_gpu.
 * Kernel to precompute metric quantity f4_of_xx1__D1.
 */
__global__ static void rfm_precompute_defines__f4_of_xx1__D1_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x1) {
  // Temporary parameters
  const REAL bScale = d_params[streamid].bScale;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i1 = tid0; i1 < d_params[streamid].Nxx_plus_2NGHOSTS1; i1 += stride0) {
    const REAL xx1 = x1[i1];
    rfmstruct->f4_of_xx1__D1[i1] = bScale * cosf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f4_of_xx1__D1_gpu
/**
 * Kernel: rfm_precompute_defines__f4_of_xx1__DD11_gpu.
 * Kernel to precompute metric quantity f4_of_xx1__DD11.
 */
__global__ static void rfm_precompute_defines__f4_of_xx1__DD11_gpu(const size_t streamid, rfm_struct *restrict rfmstruct, const REAL *restrict x1) {
  // Temporary parameters
  const REAL bScale = d_params[streamid].bScale;
  // Kernel thread/stride setup
  const int tid0 = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride0 = blockDim.x * gridDim.x;

  for (int i1 = tid0; i1 < d_params[streamid].Nxx_plus_2NGHOSTS1; i1 += stride0) {
    const REAL xx1 = x1[i1];
    rfmstruct->f4_of_xx1__DD11[i1] = -bScale * sinf(xx1);
  }
} // END FUNCTION rfm_precompute_defines__f4_of_xx1__DD11_gpu

/**
 * rfm_precompute_defines: reference metric precomputed lookup arrays: defines
 */
void rfm_precompute_defines__rfm__SinhSymTP(const commondata_struct *restrict commondata, const params_struct *restrict params,
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
    rfm_precompute_defines__f2_of_xx0_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x0);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f2_of_xx0_gpu failure");
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
    rfm_precompute_defines__f2_of_xx0__D0_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x0);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f2_of_xx0__D0_gpu failure");
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
    rfm_precompute_defines__f2_of_xx0__DD00_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x0);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f2_of_xx0__DD00_gpu failure");
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
    rfm_precompute_defines__f4_of_xx1_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x1);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f4_of_xx1_gpu failure");
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
    rfm_precompute_defines__f4_of_xx1__D1_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x1);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f4_of_xx1__D1_gpu failure");
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
    rfm_precompute_defines__f4_of_xx1__DD11_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct, x1);
    cudaCheckErrors(cudaKernel, "rfm_precompute_defines__f4_of_xx1__DD11_gpu failure");
  }
} // END FUNCTION rfm_precompute_defines__rfm__SinhSymTP
