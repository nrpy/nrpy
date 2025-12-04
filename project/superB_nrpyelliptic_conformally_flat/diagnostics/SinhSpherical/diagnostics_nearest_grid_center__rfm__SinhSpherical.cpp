#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "diagnostics/diagnostics_nearest_common.h"

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
void diagnostics_nearest_grid_center__rfm__SinhSpherical(commondata_struct *restrict commondata, const int grid, const params_struct *restrict params,
                                                         const params_struct *restrict params_chare,
                                                         const REAL *restrict xx[3], const int NUM_GFS_NEAREST, const int which_gfs[],
                                                         const char **diagnostic_gf_names, const REAL *restrict gridfuncs_diags[],
                                                         const int chare_index[3]) {
  // Suppress unused parameter warning for xx, required for API compatibility.
  (void)xx;

  // Build coordsys string with runtime coordinate system name and grid number.
  char coordsys_with_grid[128];
  snprintf(coordsys_with_grid, sizeof(coordsys_with_grid), "grid%02d-%s", grid, params->CoordSystemName);

  // Persistent per-grid file (append across timesteps).
  FILE *out = open_outfile("out0d", coordsys_with_grid, commondata, /*include_time=*/0);

  if (commondata->nn == 0)
    diag_write_header(out, "time", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);

  // Nearest-to-center indices specialized for SinhSpherical.
  const int i0_center = NGHOSTS;
  const int i1_center = params->Nxx_plus_2NGHOSTS1 / 2;
  const int i2_center = params->Nxx_plus_2NGHOSTS2 / 2;
  
  const int Nxx0chare = params_chare->Nxx0;
  const int Nxx1chare = params_chare->Nxx1;
  const int Nxx2chare = params_chare->Nxx2;
      
  const int i0_center_local = MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0_center, Nxx0chare);
  const int i1_center_local = MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1_center, Nxx1chare);
  const int i2_center_local = MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2_center, Nxx2chare);

  const int idx3 = IDX3P(params_chare, i0_center_local, i1_center_local, i2_center_local);

  // Active grid data pointer.
  const REAL *restrict src = gridfuncs_diags[grid];

  // Emit a single output row: [time, sampled values...]
  const int NUM_COLS = 1 + NUM_GFS_NEAREST;
  REAL row[NUM_COLS];
  row[0] = commondata->time;
  for (int ii = 0; ii < NUM_GFS_NEAREST; ii++) {
    const int gf = which_gfs[ii];
    row[1 + ii] = src[IDX4Ppt(params_chare, gf, idx3)];
  } // END LOOP over gridfunctions

  diag_write_row(out, NUM_COLS, row);
  fclose(out);
} // END FUNCTION diagnostics_nearest_grid_center__rfm__SinhSpherical
