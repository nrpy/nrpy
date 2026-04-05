"""
Register C function for applying inner boundary conditions.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys

import nrpy.c_function as cfc


# apply_bcs_inner_only(): Apply inner boundary conditions.
# This function is documented below in 'desc' and 'body' variables.
def register_CFunction_apply_bcs_inner_only() -> None:
    """
    Register C function for filling inner boundary points on the computational grid, as prescribed by bcstruct.

    DocTests:
    >>> register_CFunction_apply_bcs_inner_only()
    """
    includes = ["BHaH_defines.h"]
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
    body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  // collapse(2) results in a nice speedup here, esp in 2D. Two_BHs_collide goes from
    //    5550 M/hr to 7264 M/hr on a Ryzen 9 5950X running on all 16 cores with core affinity.
#pragma omp parallel for collapse(2)  // spawn threads and distribute across them
  for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
    for(int pt=0;pt<bc_info->num_inner_boundary_points;pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      //  -> idx3 = i0 + Nx0*(i1 + Nx1*i2)
      //  -> i0 = mod(idx3, Nx0)
      // Only apply boundary condition if at the radial interior point (i0 == NGHOSTS).
      if (dstpt % Nxx_plus_2NGHOSTS0 != NGHOSTS)
        continue;

      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;

      gfs[IDX4pt(which_gf, dstpt)] = bcstruct->inner_bc_array[pt].parity[evol_gf_parity[which_gf]] * gfs[IDX4pt(which_gf, srcpt)];
    } // END for(int pt=0;pt<num_inner_pts;pt++)
  } // END for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++)
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
