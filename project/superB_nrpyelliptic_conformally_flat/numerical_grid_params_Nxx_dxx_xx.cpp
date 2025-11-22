#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Initializes a cell-centered grid in the local coordinate system based on physical dimensions (grid_physical_size).
 *
 * Inputs:
 * - Nx[] inputs: Specifies new grid dimensions, if needed.
 * - params.convergence_factor (set to 1.0 by default): Factor by which grid resolution is increased; set to 1.0 by default.
 * - apply_convergence_factor_and_set_xxminmax_defaults: Whether to apply convergence factor, and set xxmin[3], xxmax[3] to default values set in
 * reference_metric.py.
 *
 * Parameter outputs:
 * - Nxx: Number of grid points in each direction.
 * - Nxx_plus_2NGHOSTS: Total grid points including ghost zones.
 * - dxx: Grid spacing.
 * - invdxx: Inverse of grid spacing.
 *
 * Grid setup output:
 * - xx: Coordinate values for each (cell-centered) grid point.
 */
void numerical_grid_params_Nxx_dxx_xx(const commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                                      const int Nx[3], const bool apply_convergence_factor_and_set_xxminmax_defaults) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    numerical_grid_params_Nxx_dxx_xx__rfm__SinhSpherical(commondata, params, xx, Nx, apply_convergence_factor_and_set_xxminmax_defaults);
    break;
  default:
    fprintf(stderr, "ERROR in numerical_grid_params_Nxx_dxx_xx(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION numerical_grid_params_Nxx_dxx_xx
