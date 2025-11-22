#include "BHaH_defines.h"

/**
 * Free intermediate-level (k_i) storage for the "RK4" Method of Lines (MoL) scheme.
 *
 * This routine is registered as "MoL_free_intermediate_stage_gfs".
 * It frees intermediate-level gridfunction storage.
 *
 * @param commondata Pointer to common, read-only runtime data. Included for a uniform
 *                   function signature across BHaH routines; it is not modified here.
 * @param params     Pointer to grid parameter struct providing Nxx_plus_2NGHOSTS0/1/2 and
 *                   related metadata needed to compute allocation sizes.
 * @param gridfuncs  Pointer to the MoL gridfunctions struct whose intermediate-level
 *                   arrays will be allocated by this routine.
 *
 * @return void
 */
void MoL_free_intermediate_stage_gfs(MoL_gridfunctions_struct *restrict gridfuncs) {
  BHAH_FREE(gridfuncs->y_nplus1_running_total_gfs);
  BHAH_FREE(gridfuncs->k_odd_gfs);
  BHAH_FREE(gridfuncs->k_even_gfs);
} // END FUNCTION MoL_free_intermediate_stage_gfs
