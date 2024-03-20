#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Method of Lines (MoL) for "RK4" method: Allocate memory for "y_n_gfs" gridfunctions
 * - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep
 * - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method
 */
void MoL_malloc_y_n_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                        MoL_gridfunctions_struct *restrict gridfuncs) {
#include "set_CodeParameters.h"
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
  gridfuncs->y_n_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);

  gridfuncs->diagnostic_output_gfs = gridfuncs->y_nplus1_running_total_gfs;
  gridfuncs->diagnostic_output_gfs2 = gridfuncs->k_odd_gfs;
}
