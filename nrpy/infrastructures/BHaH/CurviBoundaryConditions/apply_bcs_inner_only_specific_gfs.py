# nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only.py
"""
Generates C functions for applying inner boundary conditions on a curvilinear grid.

The generated functions apply parity conditions to grid functions at inner boundary
points (e.g., for black hole excision). The boundary mappings are pre-computed
and stored in a `bc_struct`. This module provides routines to apply these
conditions to either all evolved grid functions or a specific list of auxiliary
grid functions.

This process is documented in the NRPy tutorial:
Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot* com
"""

import nrpy.c_function as cfc


def register_CFunction_apply_bcs_inner_only_specific_gfs() -> None:
    """
    Register C function for filling inner boundary points for specific gfs on the computational grid.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp"]
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    cfc.CFunction_dict.clear()
    ...    register_CFunction_apply_bcs_inner_only_specific_gfs()
    ...    generated_str = cfc.CFunction_dict[f'apply_bcs_inner_only_specific_gfs'].full_function
    ...    validation_desc = f"apply_bcs_inner_only_specific_gfs__{parallelization}"
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    # Presence of simd_intrinsics.h forces this code to be host (CPU)-only:
    includes = ["BHaH_defines.h", "intrinsics/simd_intrinsics.h"]
    desc = r"""
WARNING: CPU ONLY.

Apply BCs to specific grid functions at inner boundary points only,
using data stored in bcstruct->inner_bc_array.
These structs are set in bcstruct_set_up().
Inner boundary points map to either the grid
interior ("pure inner") or to pure outer
boundary points ("inner maps to outer").
"""
    cfunc_type = "void"
    name = "apply_bcs_inner_only_specific_gfs"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs, const int num_gfs, const int8_t *gf_parities, const int *gfs_to_sync"
    body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;
  const innerpt_bc_struct *restrict inner_bc_array = bcstruct->inner_bc_array;
  const int num_inner_boundary_points = bc_info->num_inner_boundary_points;
  // Needed for IDX macros
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  // collapse(2) results in a nice speedup here, esp in 2D. Two_BHs_collide goes from
  //    5550 M/hr to 7264 M/hr on a Ryzen 9 5950X running on all 16 cores with core affinity.
#pragma omp parallel for collapse(2) // spawn threads and distribute across them
  for (int which_gf = 0; which_gf < num_gfs; which_gf++) {
    for (int pt = 0; pt < num_inner_boundary_points; pt++) {
      const int dstpt = inner_bc_array[pt].dstpt;
      const int srcpt = inner_bc_array[pt].srcpt;
      gfs[IDX4pt(gfs_to_sync[which_gf], dstpt)] =
          inner_bc_array[pt].parity[gf_parities[gfs_to_sync[which_gf]]] * gfs[IDX4pt(gfs_to_sync[which_gf], srcpt)];
    } // END LOOP over inner boundary points
  } // END LOOP over specific gridfunctions
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
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
