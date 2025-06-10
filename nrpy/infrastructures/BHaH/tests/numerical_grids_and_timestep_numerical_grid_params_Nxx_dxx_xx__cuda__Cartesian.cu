#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"

#define SET_XX_CELL_CENTERED_COORDS(COORD_DIR)                                                                                                       \
  const int index = blockIdx.x * blockDim.x + threadIdx.x;                                                                                           \
  const int stride = blockDim.x * gridDim.x;                                                                                                         \
  static constexpr REAL onehalf = 1.0 / 2.0;                                                                                                         \
  for (int j = index; j < d_params[streamid].Nxx_plus_2NGHOSTS##COORD_DIR; j += stride) {                                                            \
    xx##COORD_DIR[j] = d_params[streamid].xxmin##COORD_DIR + ((REAL)(j - NGHOSTS) + onehalf) * d_params[streamid].dxx##COORD_DIR;                    \
  }
/**
 * Kernel: initialize_grid_xx0_gpu.
 * Kernel to compute xx0 coordinates.
 */
__global__ static void initialize_grid_xx0_gpu(const size_t streamid, REAL *restrict xx0) {
  SET_XX_CELL_CENTERED_COORDS(0);
} // END FUNCTION initialize_grid_xx0_gpu
/**
 * Kernel: initialize_grid_xx1_gpu.
 * Kernel to compute xx1 coordinates.
 */
__global__ static void initialize_grid_xx1_gpu(const size_t streamid, REAL *restrict xx1) {
  SET_XX_CELL_CENTERED_COORDS(1);
} // END FUNCTION initialize_grid_xx1_gpu
/**
 * Kernel: initialize_grid_xx2_gpu.
 * Kernel to compute xx2 coordinates.
 */
__global__ static void initialize_grid_xx2_gpu(const size_t streamid, REAL *restrict xx2) {
  SET_XX_CELL_CENTERED_COORDS(2);
} // END FUNCTION initialize_grid_xx2_gpu

/**
 * Initializes a cell-centered grid in Cartesian coordinates based on physical dimensions (grid_physical_size).
 *
 * Inputs:
 * - Nx[] inputs: Specifies new grid dimensions, if needed.
 * - params.convergence_factor (set to 1.0 by default): Factor by which grid resolution is increased; set to 1.0 by default.
 * - apply_convergence_factor_and_set_xxminmax_defaults: Whether to set xxmin[3], xxmax[3] to default values set in reference_metric.py.
 *
 * Parameter outputs:
 * - Nxx: Number of grid points in each direction.
 * - Nxx_plus_2NGHOSTS: Total grid points including ghost zones.
 * - dxx: Grid spacing.
 * - invdxx: Inverse of grid spacing.
 *
 * Grid setup output:
 * - xx: Coordinate values for each (cell-centered) grid point.
 *
 */
void numerical_grid_params_Nxx_dxx_xx__rfm__Cartesian(const commondata_struct *restrict commondata, params_struct *restrict params, REAL *xx[3],
                                                      const int Nx[3], const bool apply_convergence_factor_and_set_xxminmax_defaults) {
  // Set default values for the grid resolution in each dimension.
  params->Nxx0 = 72;
  params->Nxx1 = 12;
  params->Nxx2 = 8;

  // If all components of Nx[] are set to valid values (i.e., not -1), override the default values with Nx[].
  if (Nx[0] != -1 && Nx[1] != -1 && Nx[2] != -1) {
    params->Nxx0 = Nx[0];
    params->Nxx1 = Nx[1];
    params->Nxx2 = Nx[2];
  }
  snprintf(params->CoordSystemName, 100, "Cartesian");

  // Resize grid by convergence_factor; used for convergence testing.
  {
    // convergence_factor does not increase resolution across an axis of symmetry (Nxx == 2):
    if (params->Nxx0 != 2)
      params->Nxx0 *= commondata->convergence_factor;
    if (params->Nxx1 != 2)
      params->Nxx1 *= commondata->convergence_factor;
    if (params->Nxx2 != 2)
      params->Nxx2 *= commondata->convergence_factor;
  }

  // Set the full grid size; including the ghostzones (of width NGHOSTS) on the boundaries.
  params->Nxx_plus_2NGHOSTS0 = params->Nxx0 + 2 * NGHOSTS;
  params->Nxx_plus_2NGHOSTS1 = params->Nxx1 + 2 * NGHOSTS;
  params->Nxx_plus_2NGHOSTS2 = params->Nxx2 + 2 * NGHOSTS;

  {
#include "../set_CodeParameters.h"
    // Set grid size to a function of grid_physical_size (grid_physical_size set in set_CodeParameters.h above):
    params->xmin = -grid_physical_size;
    params->ymin = -grid_physical_size;
    params->zmin = -grid_physical_size;
    params->xmax = grid_physical_size;
    params->ymax = grid_physical_size;
    params->zmax = grid_physical_size;
  }
  if (apply_convergence_factor_and_set_xxminmax_defaults) {
#include "../set_CodeParameters.h"
    // Set {xxmin[], xxmax[]} to default values, which could be functions of other rfm params (set in set_CodeParameters.h above):
    params->xxmin0 = xmin;
    params->xxmin1 = ymin;
    params->xxmin2 = zmin;
    params->xxmax0 = xmax;
    params->xxmax1 = ymax;
    params->xxmax2 = zmax;
  }

  // Set quantities that depend on Nxx and {xxmin, xxmax}: dxx, invdxx.
  params->dxx0 = (params->xxmax0 - params->xxmin0) / ((REAL)params->Nxx0);
  params->dxx1 = (params->xxmax1 - params->xxmin1) / ((REAL)params->Nxx1);
  params->dxx2 = (params->xxmax2 - params->xxmin2) / ((REAL)params->Nxx2);

  params->invdxx0 = ((REAL)params->Nxx0) / (params->xxmax0 - params->xxmin0);
  params->invdxx1 = ((REAL)params->Nxx1) / (params->xxmax1 - params->xxmin1);
  params->invdxx2 = ((REAL)params->Nxx2) / (params->xxmax2 - params->xxmin2);

  // Set up uniform, cell-centered, topologically Cartesian numerical grid,
  //   centered at (xxmin[i] + xxmax[i])/2 in direction i, and store
  //   {xx[0], xx[1], xx[2]} arrays.
  BHAH_MALLOC_DEVICE(xx[0], sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC_DEVICE(xx[1], sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
  BHAH_MALLOC_DEVICE(xx[2], sizeof(REAL) * params->Nxx_plus_2NGHOSTS2);
  cpyHosttoDevice_params__constant(params, params->grid_idx % NUM_STREAMS);
  {
    const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_DEFAULT;
    const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_DEFAULT;
    const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_DEFAULT;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((params->Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    initialize_grid_xx0_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, xx[0]);
    cudaCheckErrors(cudaKernel, "initialize_grid_xx0_gpu failure");
  }
  {
    const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_DEFAULT;
    const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_DEFAULT;
    const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_DEFAULT;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((params->Nxx_plus_2NGHOSTS1 + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    initialize_grid_xx1_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, xx[1]);
    cudaCheckErrors(cudaKernel, "initialize_grid_xx1_gpu failure");
  }
  {
    const size_t threads_in_x_dir = BHAH_THREADS_IN_X_DIR_DEFAULT;
    const size_t threads_in_y_dir = BHAH_THREADS_IN_Y_DIR_DEFAULT;
    const size_t threads_in_z_dir = BHAH_THREADS_IN_Z_DIR_DEFAULT;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid((params->Nxx_plus_2NGHOSTS2 + threads_in_x_dir - 1) / threads_in_x_dir, 1, 1);
    size_t sm = 0;
    size_t streamid = params->grid_idx % NUM_STREAMS;
    initialize_grid_xx2_gpu<<<blocks_per_grid, threads_per_block, sm, streams[streamid]>>>(streamid, xx[2]);
    cudaCheckErrors(cudaKernel, "initialize_grid_xx2_gpu failure");
  }
} // END FUNCTION numerical_grid_params_Nxx_dxx_xx__rfm__Cartesian
