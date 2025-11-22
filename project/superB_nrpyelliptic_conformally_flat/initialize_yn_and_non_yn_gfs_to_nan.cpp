#include "BHaH_defines.h"

/**
 * Initialize yn and non-yn gfs to nan
 */
void initialize_yn_and_non_yn_gfs_to_nan(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                         MoL_gridfunctions_struct *restrict gridfuncs) {
#include "set_CodeParameters.h"
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
  for (int i = 0; i < NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot; i++) {
    gridfuncs->y_nplus1_running_total_gfs[i] = NAN;
    gridfuncs->k_odd_gfs[i] = NAN;
    gridfuncs->k_even_gfs[i] = NAN;
    gridfuncs->y_n_gfs[i] = NAN;
  } // END LOOP over NUM_EVOL_GFS
} // END FUNCTION initialize_yn_and_non_yn_gfs_to_nan
