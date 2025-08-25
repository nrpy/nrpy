#include "BHaH_defines.h"
/**
 * Set params_struct to default values specified within NRPy+.
 */
void params_struct_set_to_default(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  // Loop over params structs:
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    params_struct *restrict params = &griddata[grid].params;
    // Set params_struct variables to default
    params->Cart_originx = 0.0;   // nrpy.grid::Cart_originx
    params->Cart_originy = 0.0;   // nrpy.grid::Cart_originy
    params->Cart_originz = 0.0;   // nrpy.grid::Cart_originz
    params->grid_rotates = false; // nrpy.grid::grid_rotates
  } // END LOOP over grids
} // END FUNCTION params_struct_set_to_default
