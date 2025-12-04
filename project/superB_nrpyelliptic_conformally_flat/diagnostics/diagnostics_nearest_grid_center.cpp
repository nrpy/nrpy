#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @file diagnostics_nearest_grid_center.c
 * @brief Write 0D diagnostics for a single grid by sampling, without interpolation, the grid
 *        point nearest the physical grid center, and append one row per timestep to a persistent file.
 *
 * The function "diagnostics_nearest_grid_center" appends diagnostics for a single grid to a per-grid
 * text file whose name encodes the runtime coordinate system, grid number, and convergence factor:
 *
 *   out0d-<CoordSystemName>-grid<XX>-conv_factor-<CF>.txt
 *
 * On the first timestep, a header is written. On all timesteps, one row is appended containing the
 * simulation time followed by sampled gridfunction values. The nearest-to-center indices are chosen
 * based on the coordinate system: for Cartesian-like systems, the midpoint index is chosen in each
 * dimension; for Spherical, Cylindrical, and SymTP, the radial index is set to NGHOSTS and midpoint
 * indices are chosen in the angular directions. The "xx" coordinate arrays are accepted for API
 * compatibility but are not used by this routine.
 *
 * If a user-editable block is provided in the implementation, users may add custom logic such as
 * additional columns or filtering before rows are written.
 *
 * @param[in]  commondata           Pointer to global simulation metadata (e.g., time and step counters).
 * @param[in]  grid                 Grid index used for selecting data and for file naming.
 * @param[in]  params               Pointer to per-grid parameters (sizes, ghost zones, strides, names).
 * @param[in]  xx                   Per-grid coordinate arrays; accepted for interface compatibility, unused here.
 * @param[in]  NUM_GFS_NEAREST      Number of gridfunctions to sample at the selected index.
 * @param[in]  which_gfs            Array of indices identifying the source gridfunctions to sample.
 * @param[in]  diagnostic_gf_names  Array of human-readable names corresponding to the sampled gridfunctions.
 * @param[in]  gridfuncs_diags      Array of pointers to per-grid source gridfunction data; this routine reads
 *                                  from gridfuncs_diags[grid].
 *
 * @return     void.
 *
 */
void diagnostics_nearest_grid_center(commondata_struct *restrict commondata, const int grid, const params_struct *restrict params, const params_struct *restrict params_chare,
                                     const REAL *restrict xx[3], const int NUM_GFS_NEAREST, const int which_gfs[], const char **diagnostic_gf_names,
                                     const REAL *restrict gridfuncs_diags[], const int chare_index[3]) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    diagnostics_nearest_grid_center__rfm__SinhSpherical(commondata, grid, params, params_chare, xx, NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names,
                                                        gridfuncs_diags, chare_index);
    break;
  default:
    fprintf(stderr, "ERROR in diagnostics_nearest_grid_center(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION diagnostics_nearest_grid_center
