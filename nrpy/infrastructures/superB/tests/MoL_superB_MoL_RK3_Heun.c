#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Method of Lines (MoL) for "RK3 Heun" method: Step forward one full timestep.
 */
void MoL_step_forward_in_time(commondata_struct *restrict commondata, griddata_struct *restrict griddata, const REAL time_start,
                              const int which_RK_substep, const int which_MOL_part) {
  // C code implementation of -={ RK3 Heun }=- Method of Lines timestepping.

  switch (which_RK_substep) {

    // Diagonal RK3 uses three evolved-field arrays and three RHS evaluations.
    // Save the k1 contribution before projecting the first provisional stage.
  case RK_SUBSTEP_K1: {
    // -={ START k1 substep }=-
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
      commondata->time = time_start + 0.00000000000000000e+00 * commondata->dt;
      // Set gridfunction aliases from gridfuncs struct
      // y_n gridfunctions
      REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
      // Temporary timelevel & AUXEVOL gridfunctions:
      REAL *restrict k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs = griddata[grid].gridfuncs.k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs;
      REAL *restrict k2_or_y_nplus_a32_k2_gfs = griddata[grid].gridfuncs.k2_or_y_nplus_a32_k2_gfs;
      REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
      params_struct *restrict params = &griddata[grid].params;
      REAL *restrict xx[3];
      for (int ww = 0; ww < 3; ww++)
        xx[ww] = griddata[grid].xx[ww];
      const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
      const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
      const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;

      switch (which_MOL_part) {
      case MOL_PRE_RK_UPDATE: {
        rhs_eval(commondata, params, y_n_gfs, k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs);

        break;
      }
      case MOL_RK_UPDATE: {
        LOOP_ALL_GFS_GPS(i) {
          const REAL k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfsL = k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i];
          const REAL y_n_gfsL = y_n_gfs[i];
          k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i] =
              (1.0 / 3.0) * commondata->dt * k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfsL + y_n_gfsL;
          k2_or_y_nplus_a32_k2_gfs[i] = (1.0 / 4.0) * commondata->dt * k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfsL + y_n_gfsL;
        }

        break;
      }

      case MOL_POST_RK_UPDATE_APPLY_BCS: {
        apply_bcs(params, k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs);

        break;
      }

      case MOL_POST_RK_UPDATE: {
        project_constraints(params, k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs);
        break;
      }
      }
    }

    break;
  }
    // Evaluate k2 from the projected stage without changing that input.
    // The fused RHS advances y_n to the next stage and adds k2 to the accumulator.
  case RK_SUBSTEP_K2: {
    // -={ START k2 substep }=-
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
      commondata->time = time_start + 3.33333333333333315e-01 * commondata->dt;
      // Set gridfunction aliases from gridfuncs struct
      // y_n gridfunctions
      REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
      // Temporary timelevel & AUXEVOL gridfunctions:
      REAL *restrict k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs = griddata[grid].gridfuncs.k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs;
      REAL *restrict k2_or_y_nplus_a32_k2_gfs = griddata[grid].gridfuncs.k2_or_y_nplus_a32_k2_gfs;
      REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
      params_struct *restrict params = &griddata[grid].params;
      REAL *restrict xx[3];
      for (int ww = 0; ww < 3; ww++)
        xx[ww] = griddata[grid].xx[ww];
      const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
      const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
      const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;

      switch (which_MOL_part) {
      case MOL_PRE_RK_UPDATE: {
        rhs_eval_accumulate(commondata, params, k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs, y_n_gfs, k2_or_y_nplus_a32_k2_gfs,
                            ((2.0 / 3.0) * commondata->dt), (0));

        break;
      }
      case MOL_RK_UPDATE: {
        break;
      }

      case MOL_POST_RK_UPDATE_APPLY_BCS: {
        apply_bcs(params, y_n_gfs);

        break;
      }

      case MOL_POST_RK_UPDATE: {
        project_constraints(params, y_n_gfs);
        break;
      }
      }
    }

    break;
  }
    // Only the stage state is projected; the accumulator stays unprojected.
  case RK_SUBSTEP_K3: {
    // -={ START k3 substep }=-
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
      commondata->time = time_start + 6.66666666666666630e-01 * commondata->dt;
      // Set gridfunction aliases from gridfuncs struct
      // y_n gridfunctions
      REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
      // Temporary timelevel & AUXEVOL gridfunctions:
      REAL *restrict k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs = griddata[grid].gridfuncs.k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs;
      REAL *restrict k2_or_y_nplus_a32_k2_gfs = griddata[grid].gridfuncs.k2_or_y_nplus_a32_k2_gfs;
      REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
      params_struct *restrict params = &griddata[grid].params;
      REAL *restrict xx[3];
      for (int ww = 0; ww < 3; ww++)
        xx[ww] = griddata[grid].xx[ww];
      const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
      const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
      const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;

      switch (which_MOL_part) {
      case MOL_PRE_RK_UPDATE: {
        rhs_eval(commondata, params, y_n_gfs, k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs);

        break;
      }
      case MOL_RK_UPDATE: {
        LOOP_ALL_GFS_GPS(i) {
          const REAL k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfsL = k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i];
          const REAL k2_or_y_nplus_a32_k2_gfsL = k2_or_y_nplus_a32_k2_gfs[i];
          y_n_gfs[i] = (3.0 / 4.0) * commondata->dt * k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfsL + k2_or_y_nplus_a32_k2_gfsL;
        }

        break;
      }

      case MOL_POST_RK_UPDATE_APPLY_BCS: {
        apply_bcs(params, y_n_gfs);

        break;
      }

      case MOL_POST_RK_UPDATE: {
        project_constraints(params, y_n_gfs);
        break;
      }
      }
    }

    break;
  }
  }
} // END FUNCTION: MoL_step_forward_in_time
