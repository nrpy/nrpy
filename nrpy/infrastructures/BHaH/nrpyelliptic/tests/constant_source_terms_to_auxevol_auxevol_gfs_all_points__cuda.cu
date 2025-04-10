#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
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

        auxevol_gfs_single_point(streamid, xx0, xx1, xx2, &in_gfs[IDX4(PSI_BACKGROUNDGF, i0, i1, i2)], &in_gfs[IDX4(ADD_TIMES_AUUGF, i0, i1, i2)]);
      } // END LOOP: for (int i0 = tid0; i0 < Nxx_plus_2NGHOSTS0; i0 += stride0)
    } // END LOOP: for (int i1 = tid1; i1 < Nxx_plus_2NGHOSTS1; i1 += stride1)
  } // END LOOP: for (int i2 = tid2; i2 < Nxx_plus_2NGHOSTS2; i2 += stride2)
} // END FUNCTION auxevol_gfs_all_points_gpu

/**
 * Set AUXEVOL gridfunctions at all points.
 */
void auxevol_gfs_all_points(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  cpyHosttoDevice_commondata__constant(commondata);
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    // Unpack griddata struct:
    params_struct *restrict params = &griddata[grid].params;
    REAL *restrict x0 = griddata[grid].xx[0];
    REAL *restrict x1 = griddata[grid].xx[1];
    REAL *restrict x2 = griddata[grid].xx[2];
    REAL *restrict in_gfs = griddata[grid].gridfuncs.auxevol_gfs;

    const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_NELL_AUX;
    const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_NELL_AUX;
    const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_NELL_AUX;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((params->Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
                         (params->Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
                         (params->Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    auxevol_gfs_all_points_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, x0, x1, x2, in_gfs);
    cudaCheckErrors(cudaKernel, "auxevol_gfs_all_points_gpu failure");
  }
} // END FUNCTION auxevol_gfs_all_points
