#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @file diagnostics_nearest_1d_y_and_z_axes.c
 * @brief Write 1D diagnostics for a single grid by sampling, without interpolation, axis-aligned
 *        lines nearest to the physical y and z axes, and append rows per timestep to persistent files.
 *
 * The function "diagnostics_nearest_1d_y_and_z_axes" appends diagnostics for a single grid to two
 * per-grid text files whose names encode the runtime coordinate system, grid number, and convergence
 * factor:
 *
 *   out1d-y-<CoordSystemName>-grid<XX>-conv_factor-<CF>.txt
 *   out1d-z-<CoordSystemName>-grid<XX>-conv_factor-<CF>.txt
 *
 * For each axis, the routine:
 *   1) Selects one or more index-line samples based on the coordinate family (e.g., Cartesian,
 *      Spherical, Cylindrical, SymTP, Wedge, Spherical_Ring).
 *   2) Converts logical coordinates xx to Cartesian via xx_to_Cart and extracts y (xCart[1]) or
 *      z (xCart[2]) as the axis coordinate.
 *   3) Buffers (axis_coord, idx3) pairs, then sorts them in ascending order using qsort.
 *   4) Writes a single-line time comment, an axis-specific header ("y" or "z"), and streams rows:
 *      [axis_coord, values of selected diagnostic gridfunctions].
 *
 * Gridfunction values are loaded from the flattened diagnostic array using IDX4Ppt with a 3D index
 * constructed by IDX3P. Sampling occurs at grid points only; no interpolation is performed.
 * On allocation or file-open failure, the routine prints an error message to stderr and terminates.
 * If a user-editable block is provided in the implementation, users may add custom logic such as
 * additional columns or filtering before rows are written.
 *
 * @param[in]  commondata           Pointer to global simulation metadata (e.g., time and step counters).
 * @param[in]  grid                 Zero-based grid index used for selecting data and for file naming.
 * @param[in]  params               Pointer to per-grid parameters (sizes, ghost zones, strides, names).
 * @param[in]  xx                   Array of 3 pointers to logical coordinates per dimension; used for
 *                                  coordinate conversion to Cartesian space.
 * @param[in]  NUM_GFS_NEAREST      Number of gridfunctions to include per output row.
 * @param[in]  which_gfs            Array of length NUM_GFS_NEAREST identifying the gridfunction indices.
 * @param[in]  diagnostic_gf_names  Array of NUM_GFS_NEAREST C strings used to construct column headers.
 * @param[in]  gridfuncs_diags      Array of per-grid pointers to flattened diagnostic data; this routine
 *                                  reads from gridfuncs_diags[grid].
 *
 * @return     void
 */
void diagnostics_nearest_1d_y_and_z_axes(commondata_struct *restrict commondata, const int grid,
                                         const params_struct *restrict params, const params_struct *restrict params_chare,const REAL *restrict xx[3], const REAL *restrict xx_chare[3],
                                         const int NUM_GFS_NEAREST, const int which_gfs[],const char **diagnostic_gf_names,
                                         const REAL *restrict gridfuncs_diags[],
                                         const charecomm_struct *restrict charecommstruct, diagnostic_struct *restrict diagnosticstruct,
                                         const int chare_index[3], Ck::IO::Session token, const int which_diagnostics_part) {

  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    diagnostics_nearest_1d_y_and_z_axes__rfm__SinhSpherical(commondata, grid, params, params_chare, xx, xx_chare, NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names,
                                                            gridfuncs_diags ,charecommstruct, diagnosticstruct, chare_index, token, which_diagnostics_part);
    break;
  default:
    fprintf(stderr, "ERROR in diagnostics_nearest_1d_y_and_z_axes(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION diagnostics_nearest_1d_y_and_z_axes
