#include "BHaH_defines.h"
/**
 * Kernel: apply_bcs_inner_only_specific_auxgfs_gpu.
 * Apply BCs to inner boundary points only for specified GFs.
 */
__global__ static void apply_bcs_inner_only_specific_auxgfs_gpu(const size_t streamid, const int num_inner_boundary_points,
                                                                const innerpt_bc_struct *restrict inner_bc_array, REAL *restrict gfs,
                                                                const int num_gfs, const int *restrict gfs_to_sync) {
  // Needed for IDX macros
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

  for (int which_gf = 0; which_gf < num_gfs; which_gf++) {
    for (int pt = tid0; pt < num_inner_boundary_points; pt += stride0) {
      const int dstpt = inner_bc_array[pt].dstpt;
      const int srcpt = inner_bc_array[pt].srcpt;
      gfs[IDX4pt(gfs_to_sync[which_gf], dstpt)] =
          inner_bc_array[pt].parity[d_aux_gf_parity[gfs_to_sync[which_gf]]] * gfs[IDX4pt(gfs_to_sync[which_gf], srcpt)];
    } // END for(int pt=0;pt<num_inner_pts;pt++)
  } // END for(int which_gf=0;which_gf<num_gfs;which_gf++)
} // END FUNCTION apply_bcs_inner_only_specific_auxgfs_gpu

/**
 *
 * Apply BCs to specific grid functions at inner boundary points only,
 * using data stored in bcstruct->inner_bc_array.
 * These structs are set in bcstruct_set_up().
 * Inner boundary points map to either the grid
 * interior ("pure inner") or to pure outer
 * boundary points ("inner maps to outer").
 *
 */
void apply_bcs_inner_only_specific_auxgfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                          const bc_struct *restrict bcstruct, REAL *restrict gfs, const int num_gfs, const int *gfs_to_sync) {

  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;
  const innerpt_bc_struct *restrict inner_bc_array = bcstruct->inner_bc_array;
  const int num_inner_boundary_points = bc_info->num_inner_boundary_points;

  // Allocate device memory for gfs_to_sync
  int *gfs_to_sync_device;
  BHAH_MALLOC_DEVICE(gfs_to_sync_device, num_gfs * sizeof(int));
  cudaMemcpy(gfs_to_sync_device, gfs_to_sync, num_gfs * sizeof(int), cudaMemcpyHostToDevice);

  const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_CURVIBC_INNER;
  const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_CURVIBC_INNER;
  const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_CURVIBC_INNER;
  dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
  dim3 blocks_per_grid(MAX(1U, (num_inner_boundary_points + threads_in_x_dir - 1) / threads_in_x_dir), 1, 1);
  size_t sm = 0;
  size_t streamid = params->grid_idx % NUM_STREAMS;
  apply_bcs_inner_only_specific_auxgfs_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, num_inner_boundary_points,
                                                                                                          inner_bc_array, gfs, num_gfs, gfs_to_sync);
  cudaCheckErrors(cudaKernel, "apply_bcs_inner_only_specific_auxgfs_gpu failure");

  BHAH_DEVICE_SYNC();
  // Free device memory for gfs_to_sync
  BHAH_FREE_DEVICE(gfs_to_sync_device);
} // END FUNCTION apply_bcs_inner_only_specific_auxgfs
