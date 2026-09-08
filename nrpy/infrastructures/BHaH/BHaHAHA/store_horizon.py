"""

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""


import nrpy.c_function as cfc
from nrpy.infrastructures import BHaH


def register_CFunction_store_horizon() -> None:
    """
    Register the C function for diagnostics in the simulation.

    :return: An NRPyEnv_type object if registration is successful, otherwise None.
    """

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = r"""
"""
    desc = """Stores a coarse horizon from y_n_gfs in commondata's coarse_horizon
for use as the next resolutions initial guess.
@param commondata - Pointer to the common data structure containing simulation parameters and data.   This function conducts various diagnostic checks and computations related to
@param griddata - Pointer to the grid data structure containing grid-related parameters and functions.the apparent horizon in the BHaHAHA simulation. Diagnostics are performed
at specified iteration intervals and may include interpolation, centroid
calculations, norm evaluations."""
    cfunc_type = "void"
    name = "store_horizon"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    cfunc_decorators = r"""
#ifdef __CUDACC__
__global__
#endif
"""
    body = r"""
  int grid = 0;
  const params_struct *restrict params = &griddata[grid].params;
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  // Store horizon data including ghost zones for interpolation in the next resolution.
  const int NUM_THETA = Nxx_plus_2NGHOSTS1; // NUM_THETA needed for IDX2() macro.
  PARALLEL_2D_LOOP(i1,0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      commondata->coarse_horizon[IDX2(i1, i2)] = griddata[grid].gridfuncs.y_n_gfs[IDX4(HHGF, NGHOSTS, i1, i2)];
  } END_PARALLEL_2D_LOOP // END LOOP: over i1 (theta) and i2 (phi)

  //for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
  PARALLEL_1D_LOOP(i0, 0, Nxx_plus_2NGHOSTS0) {
    commondata->coarse_horizon_r_theta_phi[0][i0] = griddata[grid].xx[0][i0];
  } END_PARALLEL_LOOP // END LOOP: radial coordinates


  //for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
  PARALLEL_1D_LOOP(i1, 0, Nxx_plus_2NGHOSTS1) {
    commondata->coarse_horizon_r_theta_phi[1][i1] = griddata[grid].xx[1][i1];
  } END_PARALLEL_LOOP // END LOOP: theta coordinates

  //for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
  PARALLEL_1D_LOOP(i2, 0, Nxx_plus_2NGHOSTS2) {
    commondata->coarse_horizon_r_theta_phi[2][i2] = griddata[grid].xx[2][i2];
  } END_PARALLEL_LOOP // END LOOP: phi coordinates
"""
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        cfunc_decorators=cfunc_decorators,
        body=body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
