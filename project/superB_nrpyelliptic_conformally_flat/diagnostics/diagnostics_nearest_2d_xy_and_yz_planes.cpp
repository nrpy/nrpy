#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @file diagnostics_nearest_2d_xy_and_yz_planes.c
 * @brief Sample and write 2D diagnostics on the nearest xy and yz planes for a single SinhSpherical grid.
 *
 * Overview:
 * For the specified grid at the current time, this routine:
 *  - Locates interior slices that are nearest to the xy and yz planes, with selection rules
 *    specialized by the runtime coordinate system.
 *  - Converts native coordinates xx to Cartesian coordinates (x, y, z) using xx_to_Cart so that
 *    the first columns of each row contain mapped coordinates.
 *  - Writes two files per call, one for the xy plane and one for the yz plane. Each file begins
 *    with a time comment and a header, followed by one row per interior point that contains the
 *    mapped coordinates and sampled diagnostic values.
 *  - Performs sampling without interpolation; values are read directly from gridfuncs_diags[grid].
 *
 * Plane selection notes (examples, not exhaustive):
 *  - Cartesian: xy fixes i2 at mid; yz fixes i0 at mid.
 *  - Cylindrical: xy fixes z at mid; yz emits two phi slices near +/- pi/2 to realize x=0.
 *  - Spherical and SymTP: xy fixes the polar-like angle at mid; yz emits two phi-like slices near +/- pi/2.
 *  - Wedge: xy may emit two z-like slices at quarter indices; yz fixes the across-wedge index at mid.
 *
 * If a user-editable block is provided in the implementation, users may insert custom logic such as
 * adding extra columns or filtering before rows are written.
 *
 * @param[in,out] commondata            Pointer to common runtime data used for time and I/O.
 * @param[in]     grid                  Grid index to process.
 * @param[in]     params                Pointer to simulation and grid parameters (sizes, names, strides).
 * @param[in]     xx                    Native grid coordinates; xx[d][i_d] gives the coordinate along dimension d.
 * @param[in]     NUM_GFS_NEAREST       Number of diagnostic gridfunctions to sample at each interior point.
 * @param[in]     which_gfs             Array of length NUM_GFS_NEAREST specifying which gridfunctions to sample.
 * @param[in]     diagnostic_gf_names   Array of length NUM_GFS_NEAREST with human-readable names for headers.
 * @param[in]     gridfuncs_diags       Array of pointers; gridfuncs_diags[grid] points to this grid's diagnostic data.
 *
 * @return        void                  No return value. On success two text files are written and closed. Fatal I/O
 *                                      or allocation failures result in program termination.
 */
void diagnostics_nearest_2d_xy_and_yz_planes(commondata_struct *restrict commondata, const int grid,
                                         const params_struct *restrict params, const params_struct *restrict params_chare,const REAL *restrict xx[3], const REAL *restrict xx_chare[3],
                                         const int NUM_GFS_NEAREST, const int which_gfs[],const char **diagnostic_gf_names,
                                         const REAL *restrict gridfuncs_diags[],
                                         const charecomm_struct *restrict charecommstruct, diagnostic_struct *restrict diagnosticstruct,
                                         const int chare_index[3], Ck::IO::Session token, const int which_diagnostics_part) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    diagnostics_nearest_2d_xy_and_yz_planes__rfm__SinhSpherical(commondata, grid, params, params_chare, xx, xx_chare, NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names,
                                                            gridfuncs_diags ,charecommstruct, diagnosticstruct, chare_index, token, which_diagnostics_part);
    break;
  default:
    fprintf(stderr, "ERROR in diagnostics_nearest_2d_xy_and_yz_planes(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION diagnostics_nearest_2d_xy_and_yz_planes
