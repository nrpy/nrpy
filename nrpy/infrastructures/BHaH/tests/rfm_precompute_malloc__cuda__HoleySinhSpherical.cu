#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_malloc__allocate_gpu.
 * Kernel to allocate rfmstruct arrays.
 */
__global__ static void rfm_precompute_malloc__allocate_gpu(const size_t streamid, rfm_struct *restrict rfmstruct) {
  // Temporary parameters
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = d_params[streamid].Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = d_params[streamid].Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = d_params[streamid].Nxx_plus_2NGHOSTS2;
  rfmstruct->f0_of_xx0 = (REAL *)malloc(sizeof(REAL) * d_params[streamid].Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__D0 = (REAL *)malloc(sizeof(REAL) * d_params[streamid].Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__DD00 = (REAL *)malloc(sizeof(REAL) * d_params[streamid].Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__DDD000 = (REAL *)malloc(sizeof(REAL) * d_params[streamid].Nxx_plus_2NGHOSTS0);
  rfmstruct->f1_of_xx1 = (REAL *)malloc(sizeof(REAL) * d_params[streamid].Nxx_plus_2NGHOSTS1);
  rfmstruct->f1_of_xx1__D1 = (REAL *)malloc(sizeof(REAL) * d_params[streamid].Nxx_plus_2NGHOSTS1);
  rfmstruct->f1_of_xx1__DD11 = (REAL *)malloc(sizeof(REAL) * d_params[streamid].Nxx_plus_2NGHOSTS1);
} // END FUNCTION rfm_precompute_malloc__allocate_gpu

/**
 * rfm_precompute_malloc: reference metric precomputed lookup arrays: malloc
 */
void rfm_precompute_malloc__rfm__HoleySinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                    rfm_struct *restrict rfmstruct) {
  {

    const size_t threads_in_x_dir = 1;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid(1, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    rfm_precompute_malloc__allocate_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, rfmstruct);
    cudaCheckErrors(cudaKernel, "rfm_precompute_malloc__allocate_gpu failure");
  }
} // END FUNCTION rfm_precompute_malloc__rfm__HoleySinhSpherical
