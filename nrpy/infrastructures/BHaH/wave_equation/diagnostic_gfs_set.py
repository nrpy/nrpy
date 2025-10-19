"""
Register and emit the C `diagnostic_gfs_set()` routine for wave-equation diagnostics.

This module exposes a single entry point,
`register_CFunction_diagnostic_gfs_set()`, which generates and registers the
C function `diagnostic_gfs_set()`. The generated routine evaluates the
Cartesian exact solution for the selected wave type and fills per-grid,
gf-major diagnostic arrays with numerical fields, exact fields, and relative
errors.

Diagnostics produced by the generated routine (enum tokens must exist in
`diagnostic_gfs.h`).

During code generation:
  - The routine pulls the exact solution from
    `WaveEquation_solution_Cartesian(WaveType=..., default_sigma=...)`.
  - It records two dictionaries into `par.glb_extras_dict`:
      * `diagnostic_gfs_names_dict`: enum-token to short-name map
      * `diagnostic_gfs_0d1d2d_interp_dict`: selections for 0D/1D/2D interpolation
  - It includes the headers `BHaH_defines.h`, `BHaH_function_prototypes.h`,
    and `diagnostics/diagnostic_gfs.h`.

At runtime, `diagnostic_gfs_set()` loops over grids and grid points, converts
local coordinates to Cartesian, evaluates the exact fields `uexact` and
`vexact`, loads numerical fields from `y_n_gfs`, and writes diagnostics into
`diagnostic_gfs[grid][...]`. Relative errors are computed as
`(numerical - exact) / (exact + 1e-16)` to avoid divide-by-zero.

Ownership and scope:
  - The caller allocates `diagnostic_gfs[grid]` for each grid; this routine
    only writes into these arrays and does not free them.
  - Interpolation and file I/O are performed by separate drivers such as
    `diagnostics_interp()`.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.equations.wave_equation.WaveEquation_Solutions_InitialData import (
    WaveEquation_solution_Cartesian,
)


def register_CFunction_diagnostic_gfs_set(
    WaveType: str = "PlaneWave",
    default_k0: float = 1.0,
    default_k1: float = 1.0,
    default_k2: float = 1.0,
    default_sigma: float = 3.0,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Generate and register the C `diagnostic_gfs_set()` routine that fills per-grid diagnostic arrays.

    This routine:
      1) Records the call during the parallel codegen registration phase, returning None.
      2) Otherwise, constructs the C function metadata and body, including headers:
           - `BHaH_defines.h`
           - `BHaH_function_prototypes.h`
           - `diagnostics/diagnostic_gfs.h`
      3) Instantiates the exact wave solution via
         `WaveEquation_solution_Cartesian(WaveType=..., default_sigma=...)` and embeds
         its C expressions to compute `uexact` and `vexact`.
      4) Emits a grid-parallel loop that:
           - Converts local coordinates to Cartesian with `xx_to_Cart`.
           - Loads numerical fields `u` and `v` from `y_n_gfs`.
           - Stores numerical and exact fields, and relative errors:
               * DIAG_UNUM, DIAG_UEXACT, DIAG_VNUM, DIAG_VEXACT
               * DIAG_RELERROR_UUGF = (u_num - u_exact) / (u_exact + 1e-16)
               * DIAG_RELERROR_VVGF = (v_num - v_exact) / (v_exact + 1e-16)
             Also stores DIAG_GRIDINDEX = (REAL)grid.
      5) Publishes:
           - `par.glb_extras_dict["diagnostic_gfs_names_dict"]` mapping enum tokens to short names
           - `par.glb_extras_dict["diagnostic_gfs_0d1d2d_interp_dict"]` selecting 0D/1D/2D outputs
      6) Registers the function under `diagnostics/` and returns the updated NRPy environment.

    Memory and ownership:
      - The caller allocates `diagnostic_gfs[grid]` with room for all diagnostics and
        grid points. This routine does not allocate or free those buffers.

    Parallel and backend notes:
        device/host coordination; otherwise the parameter is omitted.
      - OpenMP-style loops are emitted via the `LOOP_OMP` macro when applicable.

    :param WaveType: Wave solution family to use for exact fields (for example, "PlaneWave").
    :param default_k0: Default wave vector component k0 used by the exact solution generator.
    :param default_k1: Default wave vector component k1 used by the exact solution generator.
    :param default_k2: Default wave vector component k2 used by the exact solution generator.
    :param default_sigma: Default width parameter used by the exact solution generator.

    :returns: None during the registration phase; otherwise the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "diagnostics/diagnostic_gfs.h",
    ]
    desc = """
 * @brief Populate per-grid diagnostic arrays for wave-equation fields and errors.
 *
 * @details
 * This routine fills caller-allocated diagnostic arrays with numerical fields,
 * exact fields, and relative errors for the wave equation on all grids.
 * For each grid and grid point:
 *
 *  - Convert local coordinates to Cartesian via xx_to_Cart().
 *  - Evaluate the exact solution (uexact, vexact) for the configured wave type.
 *  - Load numerical fields from y_n_gfs (UUGF and VVGF).
 *  - Store diagnostics.
 * The set of diagnostics is defined by the enum in diagnostics/diagnostic_gfs.h,
 * with short names provided by diagnostic_gf_names[]. Interpolation and output
 * are performed elsewhere by drivers such as diagnostics_interp().
 *
 * @param[in]  commondata      Global simulation constants (for example, sigma, wavespeed, time).
 * @param[in]  griddata        Per-grid data, including coordinates and y_n_gfs.
 * @param[out] diagnostic_gfs  Caller-allocated per-grid, gf-major diagnostic arrays to be filled.
 * @param[in]  griddata_host   (CUDA builds only) Host-side mirror of griddata for device access.
 *
 * @pre
 *  - diagnostic_gfs[grid] is non-null and large enough for all diagnostics and grid points.
 *  - griddata contains valid coordinates, parameters, and y_n_gfs arrays.
 *  - Enum tokens DIAG_RELERROR_UUGF, DIAG_RELERROR_VVGF, DIAG_UNUM, DIAG_UEXACT,
 *    DIAG_VNUM, DIAG_VEXACT, and DIAG_GRIDINDEX exist and are in range.
 *
 * @post
 *  - All listed diagnostics are populated for each grid point on each grid.
 *  - No memory is allocated or freed by this routine.
 *
 * @warning
 *  - Relative errors use a small denominator guard 1e-16 to avoid division by zero.
 *  - This routine only writes into diagnostic_gfs[grid]; allocation and deallocation
 *    are the caller's responsibility.
 *
 * @return void
 *
 * @see diagnostics/diagnostic_gfs.h
 * @see diagnostics_interp()
"""
    cfunc_type = "void"
    name = "diagnostic_gfs_set"
    params = "const commondata_struct *restrict commondata, const griddata_struct *restrict griddata, REAL *restrict diagnostic_gfs[MAXNUMGRIDS]"

    waveeq = WaveEquation_solution_Cartesian(
        WaveType=WaveType,
        default_k0=default_k0,
        default_k1=default_k1,
        default_k2=default_k2,
        default_sigma=default_sigma,
    )
    diagnostic_gfs_names_dict = {
        "DIAG_RELERROR_UUGF": "RelError_u",
        "DIAG_RELERROR_VVGF": "RelError_v",
        "DIAG_UNUM": "u_numerical",
        "DIAG_UEXACT": "u_exact",
        "DIAG_VNUM": "v_numerical",
        "DIAG_VEXACT": "v_exact",
        "DIAG_GRIDINDEX": "GridIndex",
    }

    body = rf"""
for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
    const params_struct *restrict params = &griddata[grid].params;
    SET_NXX_PLUS_2NGHOSTS_VARS(grid);
    const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;

    LOOP_OMP("omp parallel for", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {{
      REAL xCart[3];
      const REAL xx012[3] = {{griddata[grid].xx[0][i0], griddata[grid].xx[1][i1], griddata[grid].xx[2][i2]}};
      xx_to_Cart(params, xx012, xCart);
      const int idx3 = IDX3(i0, i1, i2);
      const REAL xCart0 = xCart[0], xCart1 = xCart[1], xCart2 = xCart[2];
      const REAL sigma = commondata->sigma;
      const REAL wavespeed = commondata->wavespeed;
      const REAL time = commondata->time;
      REAL uexact, vexact;
    { ccg.c_codegen([waveeq.uu_exactsoln, waveeq.vv_exactsoln], ["uexact", "vexact"]) }
    const REAL unum = y_n_gfs[IDX4pt(UUGF, idx3)];
    const REAL vnum = y_n_gfs[IDX4pt(VVGF, idx3)];
    diagnostic_gfs[grid][IDX4pt(DIAG_RELERROR_UUGF, idx3)] = (unum - uexact) / (uexact + 1e-16);
    diagnostic_gfs[grid][IDX4pt(DIAG_RELERROR_VVGF, idx3)] = (vnum - vexact) / (vexact + 1e-16);
    diagnostic_gfs[grid][IDX4pt(DIAG_UNUM,          idx3)] = unum;
    diagnostic_gfs[grid][IDX4pt(DIAG_UEXACT,        idx3)] = uexact;
    diagnostic_gfs[grid][IDX4pt(DIAG_VNUM,          idx3)] = vnum;
    diagnostic_gfs[grid][IDX4pt(DIAG_VEXACT,        idx3)] = vexact;
    diagnostic_gfs[grid][IDX4pt(DIAG_GRIDINDEX,     idx3)] = (REAL)grid;
  }} // END LOOP over all gridpoints to set diagnostic_gfs
}} // END LOOP over grids
"""
    par.glb_extras_dict["diagnostic_gfs_names_dict"] = diagnostic_gfs_names_dict
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
