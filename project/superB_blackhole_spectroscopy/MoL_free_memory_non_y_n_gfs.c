#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Method of Lines (MoL) for "RK4" method: Free memory for "non_y_n_gfs" gridfunctions
 * - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
 * - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method
 *
 */
void MoL_free_memory_non_y_n_gfs(MoL_gridfunctions_struct *restrict gridfuncs) {
  free(gridfuncs->y_nplus1_running_total_gfs);
  free(gridfuncs->k_odd_gfs);
  free(gridfuncs->k_even_gfs);
  if (NUM_AUXEVOL_GFS > 0)
    free(gridfuncs->auxevol_gfs);
}
