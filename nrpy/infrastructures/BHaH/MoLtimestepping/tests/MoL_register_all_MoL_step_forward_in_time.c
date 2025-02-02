#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

#define LOOP_ALL_GFS_GPS(ii)                                                                                                                         \
  _Pragma("omp parallel for") for (int(ii) = 0;                                                                                                      \
                                   (ii) < params->Nxx_plus_2NGHOSTS0 * params->Nxx_plus_2NGHOSTS1 * params->Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS;       \
                                   (ii)++)
/**
 * Runge-Kutta function for substep 1.
 */
static void rk_substep_1(params_struct *restrict params, REAL *restrict k_odd_gfs, REAL *restrict y_n_gfs, REAL *restrict y_nplus1_running_total_gfs,
                         const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k_odd_gfsL = k_odd_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    const REAL RK_Rational_1_6 = 1.0 / 6.0;
    const REAL RK_Rational_1_2 = 1.0 / 2.0;
    y_nplus1_running_total_gfs[i] = RK_Rational_1_6 * dt * k_odd_gfsL;
    k_odd_gfs[i] = RK_Rational_1_2 * dt * k_odd_gfsL + y_n_gfsL;
  }
} // END FUNCTION rk_substep_1

/**
 * Runge-Kutta function for substep 1.
 */
static void rk_substep_1__launcher(params_struct *restrict params, REAL *restrict k_odd_gfs, REAL *restrict y_n_gfs,
                                   REAL *restrict y_nplus1_running_total_gfs, const REAL dt) {
  rk_substep_1(params, k_odd_gfs, y_n_gfs, y_nplus1_running_total_gfs, dt);
} // END FUNCTION rk_substep_1__launcher

/**
 * Runge-Kutta function for substep 2.
 */
static void rk_substep_2(params_struct *restrict params, REAL *restrict k_even_gfs, REAL *restrict y_nplus1_running_total_gfs, REAL *restrict y_n_gfs,
                         const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k_even_gfsL = k_even_gfs[i];
    const REAL y_nplus1_running_total_gfsL = y_nplus1_running_total_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    const REAL RK_Rational_1_3 = 1.0 / 3.0;
    const REAL RK_Rational_1_2 = 1.0 / 2.0;
    y_nplus1_running_total_gfs[i] = RK_Rational_1_3 * dt * k_even_gfsL + y_nplus1_running_total_gfsL;
    k_even_gfs[i] = RK_Rational_1_2 * dt * k_even_gfsL + y_n_gfsL;
  }
} // END FUNCTION rk_substep_2

/**
 * Runge-Kutta function for substep 2.
 */
static void rk_substep_2__launcher(params_struct *restrict params, REAL *restrict k_even_gfs, REAL *restrict y_nplus1_running_total_gfs,
                                   REAL *restrict y_n_gfs, const REAL dt) {
  rk_substep_2(params, k_even_gfs, y_nplus1_running_total_gfs, y_n_gfs, dt);
} // END FUNCTION rk_substep_2__launcher

/**
 * Runge-Kutta function for substep 3.
 */
static void rk_substep_3(params_struct *restrict params, REAL *restrict k_odd_gfs, REAL *restrict y_nplus1_running_total_gfs, REAL *restrict y_n_gfs,
                         const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k_odd_gfsL = k_odd_gfs[i];
    const REAL y_nplus1_running_total_gfsL = y_nplus1_running_total_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    const REAL RK_Rational_1_3 = 1.0 / 3.0;
    y_nplus1_running_total_gfs[i] = RK_Rational_1_3 * dt * k_odd_gfsL + y_nplus1_running_total_gfsL;
    k_odd_gfs[i] = dt * k_odd_gfsL + y_n_gfsL;
  }
} // END FUNCTION rk_substep_3

/**
 * Runge-Kutta function for substep 3.
 */
static void rk_substep_3__launcher(params_struct *restrict params, REAL *restrict k_odd_gfs, REAL *restrict y_nplus1_running_total_gfs,
                                   REAL *restrict y_n_gfs, const REAL dt) {
  rk_substep_3(params, k_odd_gfs, y_nplus1_running_total_gfs, y_n_gfs, dt);
} // END FUNCTION rk_substep_3__launcher

/**
 * Runge-Kutta function for substep 4.
 */
static void rk_substep_4(params_struct *restrict params, REAL *restrict k_even_gfs, REAL *restrict y_n_gfs, REAL *restrict y_nplus1_running_total_gfs,
                         const REAL dt) {
  LOOP_ALL_GFS_GPS(i) {
    const REAL k_even_gfsL = k_even_gfs[i];
    const REAL y_n_gfsL = y_n_gfs[i];
    const REAL y_nplus1_running_total_gfsL = y_nplus1_running_total_gfs[i];
    const REAL RK_Rational_1_6 = 1.0 / 6.0;
    y_n_gfs[i] = RK_Rational_1_6 * dt * k_even_gfsL + y_n_gfsL + y_nplus1_running_total_gfsL;
  }
} // END FUNCTION rk_substep_4

/**
 * Runge-Kutta function for substep 4.
 */
static void rk_substep_4__launcher(params_struct *restrict params, REAL *restrict k_even_gfs, REAL *restrict y_n_gfs,
                                   REAL *restrict y_nplus1_running_total_gfs, const REAL dt) {
  rk_substep_4(params, k_even_gfs, y_n_gfs, y_nplus1_running_total_gfs, dt);
} // END FUNCTION rk_substep_4__launcher

/**
 * Method of Lines (MoL) for "RK4" method: Step forward one full timestep.
 *
 */
void MoL_step_forward_in_time(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {

  // C code implementation of -={ RK4 }=- Method of Lines timestepping.

  // First set the initial time:
  const REAL time_start = commondata->time;
  // -={ START k1 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 0.00000000000000000e+00 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict y_nplus1_running_total_gfs = griddata[grid].gridfuncs.y_nplus1_running_total_gfs;
    MAYBE_UNUSED REAL *restrict k_odd_gfs = griddata[grid].gridfuncs.k_odd_gfs;
    MAYBE_UNUSED REAL *restrict k_even_gfs = griddata[grid].gridfuncs.k_even_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, y_n_gfs, k_odd_gfs);
    rk_substep_1__launcher(params, k_odd_gfs, y_n_gfs, y_nplus1_running_total_gfs, commondata->dt);
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, k_odd_gfs);
  }
  // -={ END k1 substep }=-

  // -={ START k2 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 5.00000000000000000e-01 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict y_nplus1_running_total_gfs = griddata[grid].gridfuncs.y_nplus1_running_total_gfs;
    MAYBE_UNUSED REAL *restrict k_odd_gfs = griddata[grid].gridfuncs.k_odd_gfs;
    MAYBE_UNUSED REAL *restrict k_even_gfs = griddata[grid].gridfuncs.k_even_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, k_odd_gfs, k_even_gfs);
    rk_substep_2__launcher(params, k_even_gfs, y_nplus1_running_total_gfs, y_n_gfs, commondata->dt);
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, k_even_gfs);
  }
  // -={ END k2 substep }=-

  // -={ START k3 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 5.00000000000000000e-01 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict y_nplus1_running_total_gfs = griddata[grid].gridfuncs.y_nplus1_running_total_gfs;
    MAYBE_UNUSED REAL *restrict k_odd_gfs = griddata[grid].gridfuncs.k_odd_gfs;
    MAYBE_UNUSED REAL *restrict k_even_gfs = griddata[grid].gridfuncs.k_even_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, k_even_gfs, k_odd_gfs);
    rk_substep_3__launcher(params, k_odd_gfs, y_nplus1_running_total_gfs, y_n_gfs, commondata->dt);
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, k_odd_gfs);
  }
  // -={ END k3 substep }=-

  // -={ START k4 substep }=-
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    commondata->time = time_start + 1.00000000000000000e+00 * commondata->dt;
    // Set gridfunction aliases, from griddata[].gridfuncs.
    MAYBE_UNUSED REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict y_nplus1_running_total_gfs = griddata[grid].gridfuncs.y_nplus1_running_total_gfs;
    MAYBE_UNUSED REAL *restrict k_odd_gfs = griddata[grid].gridfuncs.k_odd_gfs;
    MAYBE_UNUSED REAL *restrict k_even_gfs = griddata[grid].gridfuncs.k_even_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    // Set pointers to this grid's params, rfm_struct/xx, bc_struct, etc.
    MAYBE_UNUSED params_struct *restrict params = &griddata[grid].params;
    MAYBE_UNUSED REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, k_odd_gfs, k_even_gfs);
    rk_substep_4__launcher(params, k_even_gfs, y_n_gfs, y_nplus1_running_total_gfs, commondata->dt);
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, y_n_gfs);
  }
  // -={ END k4 substep }=-

  // Adding dt to commondata->time many times will induce roundoff error,
  // so here we set time based on the iteration number:
  commondata->time = (REAL)(commondata->nn + 1) * commondata->dt;

  // Increment the timestep n:
  commondata->nn++;
} // END FUNCTION MoL_step_forward_in_time
