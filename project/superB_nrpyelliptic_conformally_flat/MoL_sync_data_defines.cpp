#include "BHaH_defines.h"

/**
 * Define data needed for syncing data across chares
 */
void MoL_sync_data_defines(MoL_gridfunctions_struct *restrict gridfuncs) {
  gridfuncs->num_evol_gfs_to_sync = 2;
  gridfuncs->num_auxevol_gfs_to_sync = 0;
  gridfuncs->num_aux_gfs_to_sync = 0;
  gridfuncs->max_sync_gfs = MAX3(gridfuncs->num_evol_gfs_to_sync, gridfuncs->num_auxevol_gfs_to_sync, gridfuncs->num_aux_gfs_to_sync);
  gridfuncs->evol_gfs_to_sync[0] = UUGF;
  gridfuncs->evol_gfs_to_sync[1] = VVGF;
} // END FUNCTION MoL_sync_data_defines
