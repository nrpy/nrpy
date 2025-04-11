#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Diagnostics.
 */
void diagnostics(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  // Output progress to stderr
  progress_indicator(commondata, griddata);

  // Grid data output
  const int n_step = commondata->nn, outevery = commondata->diagnostics_output_every;
  const REAL time = commondata->time;

  REAL global_norm = -1e9;
  for (int grid = 0; grid < commondata->NUMGRIDS; ++grid) {

    // Set gridfunctions aliases
    REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;
    MAYBE_UNUSED REAL *restrict diagnostic_output_gfs2 = griddata[grid].gridfuncs.diagnostic_output_gfs2;
    // Set params
    params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
    // Set rfm_struct
    const rfm_struct *restrict rfmstruct = griddata[grid].rfmstruct;

    // Compute Hamiltonian constraint violation and store it at diagnostic_output_gfs
    compute_residual_all_points(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, diagnostic_output_gfs);

    // Set integration radius for l2-norm computation
    const REAL integration_radius = 1000;

    // Compute l2-norm of Hamiltonian constraint violation
    REAL residual_H;
    compute_L2_norm_of_gridfunction(commondata, &griddata[grid].params, griddata[grid].xx, integration_radius, RESIDUAL_HGF, &residual_H,
                                    diagnostic_output_gfs);
    global_norm = MAX(global_norm, residual_H);
  } // END for(grid=0; grid<commondata->NUMGRIDS; ++grid)

  // Update residual to be used in stop condition
  commondata->log10_current_residual = global_norm;

  // Output l2-norm of Hamiltonian constraint violation to file
  {
    char filename[256];
    sprintf(filename, "residual_l2_norm.txt");
    FILE *outfile = (n_step == 0) ? fopen(filename, "w") : fopen(filename, "a");
    if (!outfile) {
      fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
      exit(1);
    }
    fprintf(outfile, "%6d %10.4e %.17e\n", n_step, time, global_norm);
    fclose(outfile);
  }

  if (n_step % outevery == 0) {
    // Only consider a single grid for now.
    const int grid = 0;
    params_struct *restrict params = &griddata[grid].params;
    // Set reference metric grid xx
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];

    // Sync with device when relevant
    BHAH_DEVICE_SYNC();

    // 1D output
    diagnostics_nearest_1d_y_axis(commondata, params, xx, &griddata[grid].gridfuncs);
    diagnostics_nearest_1d_z_axis(commondata, params, xx, &griddata[grid].gridfuncs);

    // 2D output
    diagnostics_nearest_2d_xy_plane(commondata, params, xx, &griddata[grid].gridfuncs);
    diagnostics_nearest_2d_yz_plane(commondata, params, xx, &griddata[grid].gridfuncs);
  }

  if (commondata->time + commondata->dt > commondata->t_final)
    printf("\n");
} // END FUNCTION diagnostics
