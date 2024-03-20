#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Diagnostics.
 */
void diagnostics(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {

  const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->diagnostics_output_every;
  // Explanation of the if() below:
  // Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime.
  // Step 2: Multiplying by outevery yields the exact time we should output again, t_out.
  // Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible!
  if (fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt) {
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
      // Unpack griddata struct:
      const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
      REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
      REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;
      REAL *restrict xx[3];
      {
        for (int ww = 0; ww < 3; ww++)
          xx[ww] = griddata[grid].xx[ww];
      }
      const params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

      // Constraint output
      {
        Ricci_eval(commondata, params, &griddata[grid].rfmstruct, y_n_gfs, auxevol_gfs);
        constraints_eval(commondata, params, &griddata[grid].rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_output_gfs);
      }

      // 0D output
      diagnostics_nearest_grid_center(commondata, params, &griddata[grid].gridfuncs);

      // 1D output
      diagnostics_nearest_1d_y_axis(commondata, params, xx, &griddata[grid].gridfuncs);
      diagnostics_nearest_1d_z_axis(commondata, params, xx, &griddata[grid].gridfuncs);

      // 2D output
      diagnostics_nearest_2d_xy_plane(commondata, params, xx, &griddata[grid].gridfuncs);
      diagnostics_nearest_2d_yz_plane(commondata, params, xx, &griddata[grid].gridfuncs);
      // Do psi4 output, but only if the grid is spherical-like.
      if (strstr(CoordSystemName, "Spherical") != NULL) {

        // Adjusted to match Tutorial-Start_to_Finish-BSSNCurvilinear-Two_BHs_Collide-Psi4.ipynb
        const int psi4_spinweightm2_sph_harmonics_max_l = 2;
#define num_of_R_exts 24
        const REAL list_of_R_exts[num_of_R_exts] = {10.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,  30.0,
                                                    31.0, 32.0, 33.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0};

        // Set psi4.
        psi4_part0(commondata, params, xx, y_n_gfs, diagnostic_output_gfs);
        psi4_part1(commondata, params, xx, y_n_gfs, diagnostic_output_gfs);
        psi4_part2(commondata, params, xx, y_n_gfs, diagnostic_output_gfs);
        // Decompose psi4 into spin-weight -2  spherical harmonics & output to files.
        psi4_spinweightm2_decomposition_on_sphlike_grids(commondata, params, diagnostic_output_gfs, list_of_R_exts, num_of_R_exts,
                                                         psi4_spinweightm2_sph_harmonics_max_l, xx);
      }
    }
  }
  progress_indicator(commondata, griddata);
  if (commondata->time + commondata->dt > commondata->t_final)
    printf("\n");
}
