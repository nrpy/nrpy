#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"

#define LOOP_OVER_XX(COORD_DIR)                                                                                                                      \
  static const REAL onehalf = 1.0 / 2.0;                                                                                                             \
  for (int j = 0; j < params->Nxx_plus_2NGHOSTS##COORD_DIR; j += 1) {                                                                                \
    xx##COORD_DIR[j] = params->xxmin##COORD_DIR + ((REAL)(j - NGHOSTS) + onehalf) * params->dxx##COORD_DIR;                                          \
  }
/**
 * Kernel: initialize_grid_xx0_host.
 * Kernel to compute xx0 coordinates.
 */
static void initialize_grid_xx0_host(const params_struct *restrict params, REAL *restrict xx0) {
  LOOP_OVER_XX(0);
} // END FUNCTION initialize_grid_xx0_host
/**
 * Kernel: initialize_grid_xx1_host.
 * Kernel to compute xx1 coordinates.
 */
static void initialize_grid_xx1_host(const params_struct *restrict params, REAL *restrict xx1) {
  LOOP_OVER_XX(1);
} // END FUNCTION initialize_grid_xx1_host
/**
 * Kernel: initialize_grid_xx2_host.
 * Kernel to compute xx2 coordinates.
 */
static void initialize_grid_xx2_host(const params_struct *restrict params, REAL *restrict xx2) {
  LOOP_OVER_XX(2);
} // END FUNCTION initialize_grid_xx2_host

/**
 * Initializes a cell-centered grid in SinhSymTP coordinates based on physical dimensions (grid_physical_size).
 *
 * Inputs:
 * - Nx[] inputs: Specifies new grid dimensions, if needed.
 * - params.convergence_factor (set to 1.0 by default): Factor by which grid resolution is increased; set to 1.0 by default.
 * - set_xxmin_xxmax_to_defaults: Whether to set xxmin[3], xxmax[3] to default values set in reference_metric.py.
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
void numerical_grid_params_Nxx_dxx_xx__rfm__SinhSymTP(const commondata_struct *restrict commondata, params_struct *restrict params,
                                                      REAL *restrict xx[3], const int Nx[3], const bool set_xxmin_xxmax_to_defaults) {
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
  snprintf(params->CoordSystemName, 100, "SinhSymTP");

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
    params->AMAX = grid_physical_size;
  }
  if (set_xxmin_xxmax_to_defaults) {
#include "../set_CodeParameters.h"
    // Set {xxmin[], xxmax[]} to default values, which could be functions of other rfm params (set in set_CodeParameters.h above):
    params->xxmin0 = 0;
    params->xxmin1 = 0;
    params->xxmin2 = -PI;
    params->xxmax0 = 1;
    params->xxmax1 = PI;
    params->xxmax2 = PI;
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
  BHAH_MALLOC(xx[0], sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC(xx[1], sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
  BHAH_MALLOC(xx[2], sizeof(REAL) * params->Nxx_plus_2NGHOSTS2);
  initialize_grid_xx0_host(params, xx[0]);
  initialize_grid_xx1_host(params, xx[1]);
  initialize_grid_xx2_host(params, xx[2]);
} // END FUNCTION numerical_grid_params_Nxx_dxx_xx__rfm__SinhSymTP
