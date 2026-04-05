"""
C function registration for a nearest-sampled diagnostics dispatcher.

This module provides a single function, register_CFunction_diagnostics_nearest(), which
constructs and registers the C function diagnostics_nearest(). The generated C routine
invokes specialized helpers to emit nearest-sampled diagnostics:
  - 0D: diagnostics_nearest_grid_center() samples the point nearest the grid center.
  - 1D: diagnostics_nearest_1d_y_and_z_axes() samples along the lines nearest the y and z axes.
  - 2D: diagnostics_nearest_2d_xy_and_yz_planes() samples across the planes nearest the xy and yz planes.

A single USER-EDIT block appears before the per-grid loop. In that block, users select
which_gfs arrays (by enum) for each dimensionality; these selections apply to all grids.
The third C parameter, gridfuncs_diags[grid], must point to caller-owned REAL diagnostic
gridfunctions that serve as the sampling source.

During code generation:
  - The routine is placed under the diagnostics/ subdirectory.
  - It includes the headers BHaH_defines.h, BHaH_function_prototypes.h, and diagnostic_gfs.h.

At runtime, diagnostics_nearest() loops over all grids and calls the helper routines that handle
coordinate selection, sampling, and file I/O for each dimensionality. This dispatcher allocates
no memory and frees none; ownership remains with the caller.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_diagnostics_nearest() -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register a C function that dispatches nearest-sampled diagnostics to helper routines.

    This function generates and registers the C routine diagnostics_nearest(), which contains:
      - One USER-EDIT block where users specify per-dimensional which_gfs arrays. These selections apply to all grids.
      - One loop over grids. For each grid, the routine calls the 0D, 1D, and 2D helper functions that handle
        coordinate selection, sampling from caller-provided diagnostic buffers, and file output.

    This function takes no Python parameters.

    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    TBD
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # --- C Function Registration ---
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "diagnostic_gfs.h"]
    desc = """
 * @brief Dispatch nearest-sampled diagnostics by invoking specialized helper routines for 0D, 1D, and 2D outputs.
 *
 * The diagnostics_nearest() dispatcher coordinates sampling of diagnostic gridfunction data from caller-provided
 * buffers and delegates output to three helper functions:
 *   - diagnostics_nearest_grid_center(): emits 0D diagnostics at the index triplet nearest the physical center.
 *   - diagnostics_nearest_1d_y_and_z_axes(): emits 1D diagnostics along lines nearest the y and z axes.
 *   - diagnostics_nearest_2d_xy_and_yz_planes(): emits 2D diagnostics across planes nearest the xy and yz planes.
 *
 * A single USER-EDIT block appears before the per-grid loop. In that block, users select which_gfs arrays (by enum)
 * for each dimensionality. These selections are applied uniformly to all grids. The dispatcher itself performs no
 * memory allocation or deallocation; all buffers are owned by the caller. Helper routines may perform file I/O.
 *
 * @param[in] commondata
 *   Pointer to shared simulation metadata and runtime context, including NUMGRIDS and iteration/time information.
 *
 * @param[in] griddata
 *   Pointer to an array of per-grid data structures. For grid index "grid", griddata[grid] provides parameters,
 *   coordinates, and strides required by the diagnostics helper routines.
 *
 * @param[in] gridfuncs_diags
 *   Array of length MAXNUMGRIDS. For each grid index "grid", gridfuncs_diags[grid] must point to caller-owned
 *   REAL diagnostic gridfunction data that serve as the sampling source.
 *
 * @pre
 *   - For each active grid, gridfuncs_diags[grid] is non-null and points to valid diagnostic data.
 *   - which_gfs indices selected in the USER-EDIT block map to valid diagnostic gridfunctions.
 *   - Helper symbols diagnostics_nearest_grid_center(), diagnostics_nearest_1d_y_and_z_axes(), and
 *     diagnostics_nearest_2d_xy_and_yz_planes() are available at link time.
 *
 * @post
 *   - For each grid, helper routines may emit 0D, 1D (y and z), and 2D (xy and yz) diagnostic outputs.
 *   - No memory is allocated or freed by this dispatcher.
 *
 * @return void
 *
 * @note The USER-EDIT block is for selecting which diagnostic gridfunctions to sample. Keep it concise and avoid
 *       per-grid logic there, as the dispatcher handles iteration over grids.
 """
    cfunc_type = "void"
    name = "diagnostics_nearest"
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                const REAL *restrict gridfuncs_diags[MAXNUMGRIDS]"""

    body = r"""
  // --- USER-EDIT: Select diagnostic gridfunctions to sample (applies to all grids) ---

  // 0D diagnostics: nearest point to the grid center.
  const int which_gfs_0d[] = { DIAG_RELERROR_UUGF, DIAG_RELERROR_VVGF, DIAG_UNUMGF, DIAG_UEXACTGF };

  // 1D diagnostics: nearest lines to the y and z axes.
  const int which_gfs_1d[] = { DIAG_RELERROR_UUGF, DIAG_RELERROR_VVGF, DIAG_UNUMGF, DIAG_UEXACTGF };

  // 2D diagnostics: nearest planes to the xy and yz coordinate planes.
  const int which_gfs_2d[] = { DIAG_RELERROR_UUGF, DIAG_RELERROR_VVGF, DIAG_UNUMGF, DIAG_UEXACTGF, DIAG_GRIDINDEXGF };

  // --- END USER-EDIT ---

  // Loop once over all grids and call the helpers using the selections above.
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    const params_struct *restrict params = &griddata[grid].params;
    const REAL *restrict xx[3] = { griddata[grid].xx[0], griddata[grid].xx[1], griddata[grid].xx[2] };

    // 0D
    const int NUM_nearest_GFS_0d = (int)(sizeof which_gfs_0d / sizeof which_gfs_0d[0]);
    diagnostics_nearest_grid_center(commondata, grid, params, xx, NUM_nearest_GFS_0d, which_gfs_0d,
                                    diagnostic_gf_names, gridfuncs_diags);

    // 1D
    const int NUM_nearest_GFS_1d = (int)(sizeof which_gfs_1d / sizeof which_gfs_1d[0]);
    diagnostics_nearest_1d_y_and_z_axes(commondata, grid, params, xx, NUM_nearest_GFS_1d, which_gfs_1d,
                                        diagnostic_gf_names, gridfuncs_diags);

    // 2D
    const int NUM_nearest_GFS_2d = (int)(sizeof which_gfs_2d / sizeof which_gfs_2d[0]);
    diagnostics_nearest_2d_xy_and_yz_planes(commondata, grid, params, xx, NUM_nearest_GFS_2d, which_gfs_2d,
                                            diagnostic_gf_names, gridfuncs_diags);
  } // END loop over grids
"""
    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
