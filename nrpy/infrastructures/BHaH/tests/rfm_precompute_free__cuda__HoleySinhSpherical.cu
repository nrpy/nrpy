#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_free__deallocate_gpu.
 * Kernel to deallocate rfmstruct arrays.
 */
__global__ static void rfm_precompute_free__deallocate_gpu(rfm_struct *restrict rfmstruct) {
  // Temporary parameters
  free(rfmstruct->f0_of_xx0);
  free(rfmstruct->f0_of_xx0__D0);
  free(rfmstruct->f0_of_xx0__DD00);
  free(rfmstruct->f0_of_xx0__DDD000);
  free(rfmstruct->f1_of_xx1);
  free(rfmstruct->f1_of_xx1__D1);
  free(rfmstruct->f1_of_xx1__DD11);
} // END FUNCTION rfm_precompute_free__deallocate_gpu

/**
 * rfm_precompute_free: reference metric precomputed lookup arrays: free
 */
void rfm_precompute_free__rfm__HoleySinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                  rfm_struct *restrict rfmstruct) {
  {

    const size_t threads_in_x_dir = 1;
    const size_t threads_in_y_dir = 1;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid(1, 1, 1);
    rfm_precompute_free__deallocate_gpu<<<blocks_per_grid, threads_per_block>>>(rfmstruct);
    cudaCheckErrors(cudaKernel, "rfm_precompute_free__deallocate_gpu failure");
  }
} // END FUNCTION rfm_precompute_free__rfm__HoleySinhSpherical
