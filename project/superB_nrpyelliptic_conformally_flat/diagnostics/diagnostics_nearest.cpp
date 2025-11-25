#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "diagnostic_gfs.h"

/**
 * @brief Dispatch nearest-sampled diagnostics by invoking specialized helper routines for 0D, 1D, and 2D outputs.
 *
 * The diagnostics_nearest() dispatcher coordinates sampling of diagnostic gridfunction data from caller-provided
 * buffers and delegates output to three helper functions:
 *   - diagnostics_nearest_grid_center(): emits 0D diagnostics at the index triplet nearest the physical center.
 *   - diagnostics_nearest_1d_y_and_z_axes(): emits 1D diagnostics along lines nearest the y and z axes.
 *   - diagnostics_nearest_2d_xy_and_yz_planes(): emits 2D diagnostics across planes nearest the xy and yz planes.
 *
 * A single USER-EDIT block appears before the per-grid loop. In that block, users select which_gfs arrays (by enum)
 * for each dimensionality. These selections are applied uniformly to all grids. The dispatcher itself performs no
 * memory allocation or deallocation; all buffers are owned by the caller. Helper routines may perform file I/O.
 *
 * @param[in] commondata
 *   Pointer to shared simulation metadata and runtime context, including NUMGRIDS and iteration/time information.
 *
 * @param[in] griddata
 *   Pointer to an array of per-grid data structures. For grid index "grid", griddata[grid] provides parameters,
 *   coordinates, and strides required by the diagnostics helper routines.
 *
 * @param[in] gridfuncs_diags
 *   Array of length MAXNUMGRIDS. For each grid index "grid", gridfuncs_diags[grid] must point to caller-owned
 *   REAL diagnostic gridfunction data that serve as the sampling source.
 *
 * @pre
 *   - For each active grid, gridfuncs_diags[grid] is non-null and points to valid diagnostic data.
 *   - which_gfs indices selected in the USER-EDIT block map to valid diagnostic gridfunctions.
 *   - Helper symbols diagnostics_nearest_grid_center(), diagnostics_nearest_1d_y_and_z_axes(), and
 *     diagnostics_nearest_2d_xy_and_yz_planes() are available at link time.
 *
 * @post
 *   - For each grid, helper routines may emit 0D, 1D (y and z), and 2D (xy and yz) diagnostic outputs.
 *   - No memory is allocated or freed by this dispatcher.
 *
 * @return void
 *
 * @note The USER-EDIT block is for selecting which diagnostic gridfunctions to sample. Keep it concise and avoid
 *       per-grid logic there, as the dispatcher handles iteration over grids.
 *
 */
void diagnostics_nearest(commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_chare,
                         const REAL *restrict gridfuncs_diags[MAXNUMGRIDS], const int chare_index[3], Ck::IO::Session token, const int which_diagnostics_part) {
  // --- USER-EDIT: Select diagnostic gridfunctions to sample (applies to all grids) ---

  // 0D diagnostics: nearest point to the grid center.
  const int which_gfs_0d[] = {DIAG_RESIDUAL, DIAG_UUGF};

  // 1D diagnostics: nearest lines to the y and z axes.
  const int which_gfs_1d[] = {DIAG_RESIDUAL, DIAG_UUGF};

  // 2D diagnostics: nearest planes to the xy and yz coordinate planes.
  const int which_gfs_2d[] = {DIAG_RESIDUAL, DIAG_UUGF};

  // --- END USER-EDIT ---

  // Loop once over all grids and call the helpers using the selections above.
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {

    const params_struct *restrict params = &griddata[grid].params;
    const params_struct *restrict params_chare = &griddata_chare[grid].params;
    diagnostic_struct *restrict diagnosticstruct = &griddata_chare[grid].diagnosticstruct;
    charecomm_struct *restrict charecommstruct = &griddata_chare[grid].charecommstruct;

#include "set_CodeParameters.h"

    const REAL *restrict xx[3] = {griddata[grid].xx[0], griddata[grid].xx[1], griddata[grid].xx[2]};
    const REAL *restrict xx_chare[3] = {griddata_chare[grid].xx[0], griddata_chare[grid].xx[1], griddata_chare[grid].xx[2]};


    switch (which_diagnostics_part) {

      case DIAGNOSTICS_SETUP_1D: {

        // 1D
        const int NUM_nearest_GFS_1d = (int)(sizeof which_gfs_1d / sizeof which_gfs_1d[0]);
        diagnostics_nearest_1d_y_and_z_axes(commondata, grid, params, params_chare, xx, xx_chare, NUM_nearest_GFS_1d, 0, NULL, NULL, charecommstruct, diagnosticstruct, chare_index, token, which_diagnostics_part);

        break;
      }

      case DIAGNOSTICS_SETUP_2D: {

        // 2D
        const int NUM_nearest_GFS_2d = (int)(sizeof which_gfs_2d / sizeof which_gfs_2d[0]);
        diagnostics_nearest_2d_xy_and_yz_planes(commondata, grid, params, params_chare, xx, xx_chare, NUM_nearest_GFS_2d, 0, NULL, NULL, charecommstruct, diagnosticstruct, chare_index, token, which_diagnostics_part);

        break;
      }

      case DIAGNOSTICS_WRITE_CENTER: {

        // 0D
        // Nearest-to-center indices specialized for SinhSpherical.
        const int i0_center = NGHOSTS;
        const int i1_center = Nxx_plus_2NGHOSTS1 / 2;
        const int i2_center = Nxx_plus_2NGHOSTS2 / 2;
        const int idx3 = IDX3(i0_center, i1_center, i2_center);

        if (charecommstruct->globalidx3pt_to_chareidx3[idx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {

          const int NUM_nearest_GFS_0d = (int)(sizeof which_gfs_0d / sizeof which_gfs_0d[0]);

          diagnostics_nearest_grid_center(commondata, grid, params, xx, NUM_nearest_GFS_0d, which_gfs_0d, diagnostic_gf_names, gridfuncs_diags);
        }

        break;
      }

      case DIAGNOSTICS_WRITE_Y:
      case DIAGNOSTICS_WRITE_Z: {


        // 1D
        const int NUM_nearest_GFS_1d = (int)(sizeof which_gfs_1d / sizeof which_gfs_1d[0]);
        diagnostics_nearest_1d_y_and_z_axes(commondata, grid, params, params_chare, xx, xx_chare, NUM_nearest_GFS_1d, which_gfs_1d, diagnostic_gf_names, gridfuncs_diags, charecommstruct, diagnosticstruct, chare_index, token, which_diagnostics_part);

        break;
      }

      case DIAGNOSTICS_WRITE_XY:
      case DIAGNOSTICS_WRITE_YZ: {


        // 2D
        const int NUM_nearest_GFS_2d = (int)(sizeof which_gfs_2d / sizeof which_gfs_2d[0]);
        diagnostics_nearest_2d_xy_and_yz_planes(commondata, grid, params, params_chare, xx, xx_chare, NUM_nearest_GFS_2d, which_gfs_2d, diagnostic_gf_names, gridfuncs_diags, charecommstruct, diagnosticstruct, chare_index, token, which_diagnostics_part);

        break;
      }
    }
  } // END loop over grids
} // END FUNCTION diagnostics_nearest
