"""
C function registration for a nearest-sampled diagnostics dispatcher.

This module provides a single function, register_CFunction_diagnostics_nearest(), which
constructs and registers the C function diagnostics_nearest(). The generated C routine
invokes specialized helpers to emit nearest-sampled diagnostics:
  - 0D: diagnostics_nearest_grid_center() samples the point nearest the grid center.
  - 1D: diagnostics_nearest_1d_y_and_z_axes() samples along the lines nearest the y and z axes.
  - 2D: diagnostics_nearest_2d_xy_and_yz_planes() samples across the planes nearest the xy and yz planes.

The C routine organizes a single USER-EDIT section that appears once, before the per-grid loop.
In that section, users define which_gfs arrays (by enum) for each dimensionality. These selections
apply to all grids. The third function parameter, gridfuncs_diags[grid], must point to caller-owned
REAL diagnostic gridfunctions that serve as the sampling source.

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
      - One USER-EDIT section: users specify per-dimensional which_gfs arrays (applies to all grids).
      - One loop over grids: for each grid, the helpers are invoked for 0D, 1D, and 2D outputs.

    Returns:
        None during registration, else the updated NRPy environment.

    Doctests:
        TBD
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # --- C Function Registration ---
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "diagnostic_gfs.h"]
    desc = """
 * @brief Dispatch nearest-sampled 0D, 1D, and 2D diagnostics to specialized helper routines.
 *
 * This routine coordinates diagnostic output by sampling nearest-point diagnostic
 * gridfunction data at:
 *   - 0D: the grid index triplet nearest to the physical center.
 *   - 1D: the lines nearest to the y and z axes.
 *   - 2D: the planes nearest to the xy and yz coordinate planes.
 *
 * A single USER-EDIT section defines which_gfs arrays for each dimensionality. These arrays list
 * diagnostic gridfunctions (by enum index) to sample and emit. The per-grid loop below then calls
 * the 0D, 1D, and 2D helpers for each grid using those same selections.
 *
 * This dispatcher does not allocate or free memory; the caller provides valid source buffers.
 *
 * @param[in] commondata
 *   Shared simulation metadata and runtime context (e.g., time, iteration counters, number of grids).
 *
 * @param[in] griddata
 *   Per-grid data including parameters, coordinates, and gridfunction strides required by helpers.
 *
 * @param[in] gridfuncs_diags
 *   Array of length MAXNUMGRIDS; for each grid index "grid", gridfuncs_diags[grid] must point to
 *   the REAL diagnostic gridfunction data used as the sampling source.
 *
 * @pre
 *   - For each active grid, gridfuncs_diags[grid] is non-null and points to valid diagnostic data.
 *   - Enum values referenced in which_gfs[] correspond to valid diagnostic gridfunctions.
 *   - diagnostics_nearest_grid_center(), diagnostics_nearest_1d_y_and_z_axes(), and
 *     diagnostics_nearest_2d_xy_and_yz_planes() are linked and available.
 *
 * @post
 *   - For each grid and timestep, helper routines may write out0d, out1d-y, out1d-z, out2d-xy, and out2d-yz files.
 *   - No memory is allocated or freed by this dispatcher.
 *
 * @return void
 """
    cfunc_type = "void"
    name = "diagnostics_nearest"
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                const REAL *restrict gridfuncs_diags[MAXNUMGRIDS]"""

    body = r"""
  // --- USER-EDIT: Select diagnostic gridfunctions to sample (applies to all grids) ---

  // 0D diagnostics: nearest point to the grid center.
  const int which_gfs_0d[] = { DIAG_RELERROR_UUGF, DIAG_RELERROR_VVGF, DIAG_UNUM, DIAG_UEXACT };

  // 1D diagnostics: nearest lines to the y and z axes.
  const int which_gfs_1d[] = { DIAG_RELERROR_UUGF, DIAG_RELERROR_VVGF, DIAG_UNUM, DIAG_UEXACT };

  // 2D diagnostics: nearest planes to the xy and yz coordinate planes.
  const int which_gfs_2d[] = { DIAG_RELERROR_UUGF, DIAG_RELERROR_VVGF, DIAG_UNUM, DIAG_UEXACT, DIAG_GRIDINDEX };

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
