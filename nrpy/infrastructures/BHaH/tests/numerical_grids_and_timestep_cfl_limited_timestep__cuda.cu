#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Kernel: compute_ds_min_gpu.
 * Kernel to find minimum grid spacing.
 */
__global__ static void compute_ds_min_gpu(const size_t streamid, const params_struct *restrict params, REAL *restrict xx0, REAL *restrict xx1,
                                          REAL *restrict xx2, REAL *restrict ds_min) {
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
    for (int i1 = tid1; i1 < Nxx_plus_2NGHOSTS1; i1 += stride1) {
      for (int i0 = tid0; i0 < Nxx_plus_2NGHOSTS0; i0 += stride0) {

        ds_min_single_pt(params, xx0[i0], xx1[i1], xx2[i2], &ds_min[IDX3(i0, i1, i2)]);

      } // END LOOP: for (int i0 = tid0; i0 < Nxx_plus_2NGHOSTS0; i0 += stride0)
    } // END LOOP: for (int i1 = tid1; i1 < Nxx_plus_2NGHOSTS1; i1 += stride1)
  } // END LOOP: for (int i2 = tid2; i2 < Nxx_plus_2NGHOSTS2; i2 += stride2)

} // END FUNCTION compute_ds_min_gpu

/**
 * Compute minimum timestep dt = CFL_FACTOR * ds_min.
 */
void cfl_limited_timestep(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]) {

  // Allocate memory for ds_min on the device
  const int Nxx_tot = params->Nxx_plus_2NGHOSTS0 * params->Nxx_plus_2NGHOSTS1 * params->Nxx_plus_2NGHOSTS2;
  REAL *ds_min_device;
  BHAH_MALLOC_DEVICE(ds_min_device, sizeof(REAL) * Nxx_tot);
  params_struct *device_params;
  BHAH_MALLOC_DEVICE(device_params, sizeof(params_struct));
  BHAH_MEMCPY_HOST_TO_DEVICE(device_params, params, sizeof(params_struct));

  const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_DS_MIN;
  const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_DS_MIN;
  const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_DS_MIN;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid((params->Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                       (params->Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                       (params->Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  compute_ds_min_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, device_params, xx[0], xx[1], xx[2], ds_min_device);
  cudaCheckErrors(cudaKernel, "compute_ds_min_gpu failure");

  REAL ds_min = find_global__minimum(ds_min_device, Nxx_tot);
  commondata->dt = MIN(commondata->dt, ds_min * commondata->CFL_FACTOR);
  BHAH_FREE_DEVICE(ds_min_device);
} // END FUNCTION cfl_limited_timestep
