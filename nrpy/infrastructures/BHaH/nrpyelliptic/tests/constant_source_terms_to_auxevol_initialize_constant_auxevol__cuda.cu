#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Call functions that set up all AUXEVOL gridfunctions.
 */
void initialize_constant_auxevol(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {

  // Set up variable wavespeed
  variable_wavespeed_gfs_all_points(commondata, griddata);

  // Set up all other AUXEVOL gridfunctions
  auxevol_gfs_all_points(commondata, griddata);
} // END FUNCTION initialize_constant_auxevol
