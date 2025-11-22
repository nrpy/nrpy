#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Free all memory within the griddata struct,
 * except perhaps non_y_n_gfs (e.g., after a regrid, in which non_y_n_gfs are freed first).
 */
void griddata_free(commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                   const bool free_non_y_n_gfs_and_core_griddata_pointers) {
  // Free memory allocated inside griddata[].
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    rfm_precompute_free(commondata, &griddata[grid].params, griddata[grid].rfmstruct);
    BHAH_FREE(griddata[grid].rfmstruct);

    BHAH_FREE(griddata[grid].bcstruct.inner_bc_array);
    for (int ng = 0; ng < NGHOSTS * 3; ng++) {
      BHAH_FREE(griddata[grid].bcstruct.pure_outer_bc_array[ng]);
    } // END LOOP freeing outer bc ghostzone arrays in all three directions

    BHAH_FREE(griddata[grid].gridfuncs.y_n_gfs);
    if (free_non_y_n_gfs_and_core_griddata_pointers) {
      MoL_free_intermediate_stage_gfs(&griddata[grid].gridfuncs);
      if (NUM_AUXEVOL_GFS > 0)
        BHAH_FREE(griddata[grid].gridfuncs.auxevol_gfs);
    } // END IF free_non_y_n_gfs_and_core_griddata_pointers

    for (int dirn = 0; dirn < 3; dirn++)
      BHAH_FREE(griddata[grid].xx[dirn]);
  } // END LOOP over grids

  if (free_non_y_n_gfs_and_core_griddata_pointers)
    BHAH_FREE(griddata);
} // END FUNCTION griddata_free
