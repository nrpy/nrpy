"""
C function registration for populating per-grid diagnostic arrays used by interpolation and integration routines.

This module constructs and registers the C routine "diagnostic_gfs_set".
The generated C function iterates over all grids to fill per-grid diagnostic arrays for
downstream interpolation and integration routines: it computes a residual-type diagnostic;
optionally applies inner boundary conditions using parity-consistent signs; copies selected
evolved gridfunctions from y_n_gfs into designated diagnostic channels; and records a
per-point grid identifier.

Function
--------
register_CFunction_diagnostic_gfs_set
    Construct and register the "diagnostic_gfs_set" C function.

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
    Construct and register a C function that populates per-grid wave-equation diagnostics, including exact fields and relative errors.

    This function generates and registers the C helper "diagnostic_gfs_set", which loops over all
    grids and, at each grid point, evaluates an analytic wave solution to obtain exact field values,
    reads the corresponding numerical fields from the current time level, computes relative-error
    diagnostics, copies both numerical and exact fields into designated diagnostic channels, and sets
    a per-point grid identifier. Coordinates are converted from the runtime coordinate system to
    Cartesian before evaluating the analytic solution. Runtime parameters such as time, wavespeed,
    and width are read from commondata; the provided defaults seed the generated code for cases where
    runtime values are not overridden. Each per-grid output buffer is assumed to hold
    TOTAL_NUM_DIAG_GFS times the number of points in that grid.

    :param WaveType: Name of the analytic wave solution used to generate exact fields (for example, "PlaneWave").
    :param default_k0: Default x-direction wavenumber used by the analytic solution.
    :param default_k1: Default y-direction wavenumber used by the analytic solution.
    :param default_k2: Default z-direction wavenumber used by the analytic solution.
    :param default_sigma: Default width parameter for localized analytic solutions.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    TBD
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
 * @file diagnostic_gfs_set.c
 * @brief Populate per-grid diagnostic arrays used by interpolation and integration routines.
 *
 * The function "diagnostic_gfs_set" loops over all grids and fills per-grid diagnostic arrays:
 *   1) Compute a residual-type diagnostic at all points using a helper that evaluates the
 *      finite-difference residual.
 *   2) If enabled at code-generation time, apply inner boundary conditions to that residual by
 *      copying from a source point to a destination point with a sign determined by the relevant
 *      parity, ensuring parity-consistent values near symmetry or excision boundaries.
 *   3) Copy selected evolved gridfunctions from the current time level (y_n_gfs) into designated
 *      diagnostic channels for downstream consumers.
 *   4) Set a per-point grid identifier channel to the grid index (converted to REAL).
 *
 * The routine assumes each per-grid output buffer is contiguous and large enough to store all
 * diagnostic channels:
 *     TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2
 * Loops over grid points may be parallelized with OpenMP if available.
 *
 * If a user-editable block is present in the implementation, users may add custom logic such as
 * extra diagnostics or filtering before finalizing values.
 *
 * @param[in]  commondata
 *     Pointer to global simulation metadata (e.g., counters and configuration) accessed by helpers
 *     and used to determine the number of grids to process.
 * @param[in]  griddata
 *     Pointer to an array of per-grid data. For each grid, this provides parameters, coordinates,
 *     boundary condition metadata, and gridfunctions (including y_n_gfs and any auxiliary data)
 *     referenced by this routine and its helpers.
 * @param[out] diagnostic_gfs
 *     Array of per-grid output buffers. For each grid, diagnostic_gfs[grid] must point to a buffer
 *     of size TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2.
 *
 * @return void.
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
    {ccg.c_codegen([waveeq.uu_exactsoln, waveeq.vv_exactsoln], ["uexact", "vexact"])}
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
