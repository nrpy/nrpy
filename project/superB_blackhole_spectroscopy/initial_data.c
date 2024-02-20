#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Set initial data.
 */
void initial_data(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  // Attempt to read checkpoint file. If it doesn't exist, then continue. Otherwise return.
  if (read_checkpoint(commondata, griddata))
    return;
  ID_persist_struct ID_persist;

  initialize_ID_persist_struct(commondata, &ID_persist);
  TP_solve(&ID_persist);

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    // Unpack griddata struct:
    params_struct *restrict params = &griddata[grid].params;
    initial_data_reader__convert_ADM_Cartesian_to_BSSN(commondata, params, griddata[grid].xx, &griddata[grid].bcstruct, &griddata[grid].gridfuncs,
                                                       &ID_persist, TP_Interp);
    apply_bcs_outerextrap_and_inner(commondata, params, &griddata[grid].bcstruct, griddata[grid].gridfuncs.y_n_gfs);
  }

  {
    extern void free_derivs(derivs * v, int n); // <- Needed to free memory allocated by TwoPunctures.
    // <- Free memory allocated within ID_persist.
    // Now that we're finished with par.v and par.cf_v (needed in setting up ID, we can free up memory for TwoPunctures' grids...
    free_derivs(&ID_persist.v, ID_persist.npoints_A * ID_persist.npoints_B * ID_persist.npoints_phi);
    free_derivs(&ID_persist.cf_v, ID_persist.npoints_A * ID_persist.npoints_B * ID_persist.npoints_phi);
  }
}
