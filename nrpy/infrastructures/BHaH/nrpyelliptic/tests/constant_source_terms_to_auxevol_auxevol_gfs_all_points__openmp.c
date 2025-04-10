#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Kernel: auxevol_gfs_all_points_host.
 * Kernel to initialize auxillary grid functions at all grid points.
 */
static void auxevol_gfs_all_points_host(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL *restrict x0,
                                        const REAL *restrict x1, const REAL *restrict x2, REAL *restrict in_gfs) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

#pragma omp parallel for
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        auxevol_gfs_single_point(commondata, params, xx0, xx1, xx2, &in_gfs[IDX4(PSI_BACKGROUNDGF, i0, i1, i2)],
                                 &in_gfs[IDX4(ADD_TIMES_AUUGF, i0, i1, i2)]);
      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
} // END FUNCTION auxevol_gfs_all_points_host

/**
 * Set AUXEVOL gridfunctions at all points.
 */
void auxevol_gfs_all_points(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    // Unpack griddata struct:
    params_struct *restrict params = &griddata[grid].params;
    REAL *restrict x0 = griddata[grid].xx[0];
    REAL *restrict x1 = griddata[grid].xx[1];
    REAL *restrict x2 = griddata[grid].xx[2];
    REAL *restrict in_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    auxevol_gfs_all_points_host(commondata, params, x0, x1, x2, in_gfs);
  }
} // END FUNCTION auxevol_gfs_all_points
