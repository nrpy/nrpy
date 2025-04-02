#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Kernel: compute_ds_min_host.
 * Kernel to find minimum grid spacing.
 */
static void compute_ds_min_host(const params_struct *restrict params, REAL *restrict xx0, REAL *restrict xx1, REAL *restrict xx2,
                                REAL *restrict ds_min_result) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

  REAL ds_min = 1e38;
#pragma omp parallel for reduction(min : ds_min)
  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {

        REAL local_ds_min;
        ds_min_single_pt(params, xx0[i0], xx1[i1], xx2[i2], &local_ds_min);
        ds_min = MIN(ds_min, local_ds_min);

      } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
    } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
  } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)

  *ds_min_result = ds_min;
} // END FUNCTION compute_ds_min_host

/**
 * Compute minimum timestep dt = CFL_FACTOR * ds_min.
 */
void cfl_limited_timestep(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]) {
  REAL ds_min = 1e38;
  compute_ds_min_host(params, xx[0], xx[1], xx[2], &ds_min);

  commondata->dt = MIN(commondata->dt, ds_min * commondata->CFL_FACTOR);
} // END FUNCTION cfl_limited_timestep
