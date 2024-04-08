#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Set up CoordSystem_hash for all grids
 */
void CoordSystem_hash_setup(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  // Set CoordSystem_hash, used for multi-coordinate-system evolutions. Hash is #define'd in BHaH_defines.h
  griddata[0].params.CoordSystem_hash = SINHSPHERICAL;
}
