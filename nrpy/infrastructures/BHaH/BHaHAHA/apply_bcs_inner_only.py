"""
Register C function for applying inner boundary conditions.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH.CurviBoundaryConditions.apply_bcs_inner_only import (
    APPLY_PARITY_BRANCHLESS_PREFUNC,
)


# apply_bcs_inner_only(): Apply inner boundary conditions.
# This function is documented below in 'desc' and 'body' variables.
def register_CFunction_apply_bcs_inner_only() -> None:
    """Register C function for inner boundary conditions."""
    includes = ["BHaH_defines.h"]
    prefunc = r"""
#ifdef __CUDACC__
__device__
#endif
""" 
    prefunc += APPLY_PARITY_BRANCHLESS_PREFUNC
    desc = r"""
Apply BCs to inner boundary points only,
using data stored in bcstruct->inner_bc_array.
These structs are set in bcstruct_set_up().
Inner boundary points map to either the grid
interior ("pure inner") or to pure outer
boundary points ("inner maps to outer").
"""
    cfunc_type = "void"
    name = "apply_bcs_inner_only"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs"
    cfunc_decorators = r"""
#ifdef __CUDACC__
__device__
#endif
"""
    body = r"""
#ifdef __CUDACC__
  static const int8_t evol_gf_parity[2] = {0,0};
#endif
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

#ifdef __CUDACC__
  PARALLEL_2D_LOOP(pt, 0, bc_info->num_inner_boundary_points, which_gf, 0, NUM_EVOL_GFS) {
    const int parity_idx = evol_gf_parity[which_gf];
    REAL *restrict gf = &gfs[IDX4pt(which_gf, 0)];
#else
#pragma omp parallel for schedule(static)
  for (int which_gf = 0; which_gf < NUM_EVOL_GFS; ++which_gf) {
    const int parity_idx = evol_gf_parity[which_gf];
    REAL *restrict gf = &gfs[IDX4pt(which_gf, 0)];

    for (int pt = 0; pt < bc_info->num_inner_boundary_points; ++pt) {
#endif
      const innerpt_bc_struct *restrict bc = &bcstruct->inner_bc_array[pt];
      const int dstpt = bc->dstpt;
      //  -> idx3 = i0 + Nx0*(i1 + Nx1*i2)
      //  -> i0 = mod(idx3, Nx0)
      // Only apply boundary condition if at the radial interior point (i0 == NGHOSTS).
      if (dstpt % Nxx_plus_2NGHOSTS0 != NGHOSTS)
        continue;

      const REAL v = gf[bc->srcpt];
      const int8_t p = bc->parity[parity_idx];
      gf[dstpt] = apply_parity_branchless(v, p);
#ifndef __CUDACC__
    } // END LOOP: for pt over inner boundary points
  } // END LOOP: for which_gf over evolution gridfunctions
#else
  } END_PARALLEL_2D_LOOP;
#endif
"""
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
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
