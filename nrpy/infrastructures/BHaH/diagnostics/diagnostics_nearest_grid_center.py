"""
Registers the C function for 0D diagnostics at the nearest grid-center point for a single grid.

This module provides the Python registration function for the C helper routine
`diagnostics_nearest_grid_center`. This C function samples gridfunction data
(no interpolation) at the single grid point nearest to the physical grid center
for the specified grid. The caller is responsible for looping over grids. It
appends one line per timestep to a persistent (per-grid) output file.

Functions
---------
register_CFunction_diagnostics_nearest_grid_center
    Registers the `diagnostics_nearest_grid_center` C function.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_diagnostics_nearest_grid_center(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for "0-dimensional" simulation diagnostics for a single grid: sample (no interpolation)
    at the grid's physical center and append to a per-grid file (one row per timestep).

    The generated C function takes the same inputs as the interpolation-based API used by
    `diagnostics_interp_*` helpers, enabling drop-in replacement in calling code that
    loops over grids.

    Behavior:
    - For the given grid, opens a persistent output file (same filename across timesteps) that encodes the
      runtime coordinate system name, grid number, and convergence factor.
    - At every timestep, appends a single row containing the simulation time followed by the sampled
      gridfunction values at the gridpoint nearest to the physical center (no interpolation).

    Filename format (persistent across timesteps):
      out0d-<CoordSystemName>-grid<XX>-conv_factor-<CF>.txt

    Index selection (nearest-to-center sampling):
    - Cartesian-like systems: i0=i0_mid, i1=i1_mid, i2=i2_mid (each mid ~= Nxx_plus_2NGHOSTSi/2).
    - Spherical/Cylindrical/SymTP: i0 = NGHOSTS (radial minimum), i1=i1_mid, i2=i2_mid.

    Columns:
    - "time" followed by the sampled gridfunction values (in the order given by `which_gfs`).

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :raises ValueError: If an unsupported coordinate system is specified (only systems with a defined center).
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests: TBD
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Choose expressions for the nearest-to-center grid indices based on CoordSystem.
    # Cartesian: (i0,i1,i2) = (mid,mid,mid)
    # Spherical/Cylindrical/SymTP: (i0,i1,i2) = (radial_min, mid, mid)
    if (
        ("Spherical" in CoordSystem)
        or ("Cylindrical" in CoordSystem)
        or ("SymTP" in CoordSystem)
    ):
        i0_center_expr = "NGHOSTS"
        i1_center_expr = "params->Nxx_plus_2NGHOSTS1/2"
        i2_center_expr = "params->Nxx_plus_2NGHOSTS2/2"
    elif "Cartesian" in CoordSystem:
        i0_center_expr = "params->Nxx_plus_2NGHOSTS0/2"
        i1_center_expr = "params->Nxx_plus_2NGHOSTS1/2"
        i2_center_expr = "params->Nxx_plus_2NGHOSTS2/2"
    else:
        raise ValueError(f"Unsupported CoordSystem: {CoordSystem}")

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "diagnostics/diagnostics_nearest_common.h",
    ]

    desc = f"""
 * @brief Samples and writes 0D diagnostics to a per-grid output file (one row per timestep) for a single grid.
 *
 * @details For the specified grid, this function samples (no interpolation) a set of specified gridfunctions
 *          at the single grid point nearest to the physical grid center, specialized for the
 *          coordinate system "{CoordSystem}". The results are appended to a persistent per-grid file
 *          whose filename encodes the runtime coordinate system name from params->CoordSystemName,
 *          grid number, and convergence factor.
 *
 *          Filename format (persistent across timesteps):
 *            out0d-<CoordSystemName>-grid<XX>-conv_factor-<CF>.txt
 *
 *          Index selection (nearest-to-center sampling):
 *            - Cartesian-like systems: i0=i0_mid, i1=i1_mid, i2=i2_mid (each mid ~= Nxx_plus_2NGHOSTSi/2).
 *            - Spherical/Cylindrical/SymTP: i0 = NGHOSTS (radial minimum), i1=i1_mid, i2=i2_mid.
 *
 *          Each appended row consists of: time followed by the sampled gridfunction values.
 *
 * @param[in]  commondata           Global simulation metadata (time, iteration, etc.).
 * @param[in]  grid                 Grid index used for file naming and selecting gridfuncs_diags[grid].
 * @param[in]  params               Per-grid parameters (sizes, ghost zones, strides, coordinates metadata).
 * @param[in]  xx                   Per-grid uniform coordinate arrays; xx[0], xx[1], xx[2]. Unused by this function.
 * @param[in]  NUM_GFS_NEAREST       Number of gridfunctions to sample.
 * @param[in]  which_gfs            Array of indices specifying which source gridfunctions to sample.
 * @param[in]  diagnostic_gf_names  Human-readable names for the diagnostics (used in headers).
 * @param[in]  gridfuncs_diags      Array of pointers to the source gridfunction data for each grid; this
 *                                  function samples from gridfuncs_diags[grid].
 *
 * @post Creates (on first timestep) or appends (subsequent timesteps) one line to this grid's output file.
 *
 * @return Void.
 """
    cfunc_type = "void"
    name = "diagnostics_nearest_grid_center"
    params = """commondata_struct *restrict commondata, const int grid, const params_struct *restrict params,
                const REAL *restrict xx[3], const int NUM_GFS_NEAREST, const int which_gfs[],
                const char **diagnostic_gf_names, const REAL *restrict gridfuncs_diags[]"""

    body = f"""
// Suppress unused parameter warning for xx, required for API compatibility.
(void)xx;

// Build coordsys string with runtime coordinate system name and grid number.
char coordsys_with_grid[128];
snprintf(coordsys_with_grid, sizeof(coordsys_with_grid), "grid%02d-%s", grid, params->CoordSystemName);

// Persistent per-grid file (append across timesteps).
FILE *out = open_outfile("out0d", coordsys_with_grid, commondata, /*include_time=*/0);

if (commondata->nn == 0)
  diag_write_header(out, "time", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);

// Nearest-to-center indices specialized for {CoordSystem}.
const int i0_center = {i0_center_expr};
const int i1_center = {i1_center_expr};
const int i2_center = {i2_center_expr};

const int idx3 = IDX3P(params, i0_center, i1_center, i2_center);

// Active grid data pointer.
const REAL *restrict src = gridfuncs_diags[grid];

// Emit a single output row: [time, sampled values...]
const int NUM_COLS = 1 + NUM_GFS_NEAREST;
REAL row[NUM_COLS];
row[0] = commondata->time;
for (int ii = 0; ii < NUM_GFS_NEAREST; ii++) {{
  const int gf = which_gfs[ii];
  row[1 + ii] = src[IDX4Ppt(params, gf, idx3)];
}} // END LOOP over gridfunctions

diag_write_row(out, NUM_COLS, row);
fclose(out);
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
        CoordSystem_for_wrapper_func=CoordSystem,
    )
    return pcg.NRPyEnv()
