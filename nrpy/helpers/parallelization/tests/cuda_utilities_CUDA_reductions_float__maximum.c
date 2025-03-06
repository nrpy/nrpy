#include "../BHaH_defines.h"

__global__ static void find_global__maximum__cuda(REAL *data, REAL *min, uint const data_length) {
  // shared data between all warps
  // Assumes one block = 32 warps = 32 * 32 threads
  // As of today, the standard maximum threads per
  // block is 1024 = 32 * 32
  __shared__ REAL shared_data[32];

  // Some initial value
  REAL REDUCTION_LIMIT = (REAL)0U;

  // Global data index - expecting a 1D dataset
  uint idx = threadIdx.x + blockDim.x * blockIdx.x;

  // thread index
  uint tid = threadIdx.x;

  // local thread reduction
  REAL local_reduced = REDUCTION_LIMIT;

  // warp mask - says all threads are involved in shuffle
  // 0xFFFFFFFFU in binary is 32 1's.
  unsigned mask = 0xFFFFFFFFU;

  // lane = which thread am I in the warp
  uint lane = threadIdx.x % warpSize;
  // warpID = which warp am I in the block
  uint warpID = threadIdx.x / warpSize;

  // Stride through data for each thread
  while (idx < data_length) {

    if (local_reduced < data[idx]) {
      local_reduced = data[idx];
    }

    // idx stride
    idx += gridDim.x * blockDim.x;
  }
  // Shuffle down kernel
  for (int offset = warpSize / 2; offset > 0; offset >>= 1) {
    REAL shfl = __shfl_down_sync(mask, local_reduced, offset);

    if (local_reduced < shfl) {
      local_reduced = shfl;
    }
  }
  // Shuffle results in lane 0 have the shuffle result
  if (lane == 0) {
    shared_data[warpID] = local_reduced;
  }

  // Make sure all warps in the block are synchronized
  __syncthreads();

  // Since there is only 32 partial reductions, we only
  // have one warp worth of work
  if (warpID == 0) {
    // Check to make sure we had 32 blocks of data
    if (tid < blockDim.x / warpSize) {
      local_reduced = shared_data[lane];
    } else {
      local_reduced = REDUCTION_LIMIT;
    }

    // Shuffle down kernel
    for (int offset = warpSize / 2; offset > 0; offset >>= 1) {
      REAL shfl = __shfl_down_sync(mask, local_reduced, offset);

      if (local_reduced < shfl) {
        local_reduced = shfl;
      }
    }
    if (tid == 0) {
      atomicMax((unsigned int *)min, *((unsigned int *)&local_reduced));
    }
  }
}

/**
 * Find array global maximum.
 */
REAL find_global__maximum(REAL *data, uint const data_length) {

  // This can be tested up to 1024
  uint threadCount = 32;

  // Number of blocks
  uint blockCount = (data_length + threadCount - 1) / threadCount;

  // CUDA atomics other than cas are only
  // compatible with (u)int.  To be generic
  // we use unsigned to be able to handle
  // float precision variables
  using ull = REAL;
  ull *h_reduced = (ull *)malloc(sizeof(ull));
  ull *d_reduced;
  *h_reduced = (ull)0U;

  cudaMalloc(&d_reduced, sizeof(ull));
  cudaCheckErrors(cudaMalloc, "cudaMalloc failure"); // error checking

  cudaMemcpy(d_reduced, h_reduced, sizeof(ull), cudaMemcpyHostToDevice);
  cudaCheckErrors(cudaMemcpy, "cudaCopyTo failure"); // error checking

  find_global__maximum__cuda<<<blockCount, threadCount>>>(data, d_reduced, data_length);
  cudaCheckErrors(find_min_cu, "cudaKernel - find_min_cu failed"); // error checking

  cudaMemcpy(h_reduced, d_reduced, sizeof(REAL), cudaMemcpyDeviceToHost);
  cudaCheckErrors(cudaMemcpy, "cudaCopyFrom failure"); // error checking

  cudaFree(d_reduced);
  cudaCheckErrors(cudaFree, "cudaFree failure"); // error checking

  // Recast back to result pointer type
  REAL *res = (REAL *)h_reduced;
  return *res;
} // END FUNCTION find_global__maximum
