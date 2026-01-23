"""
C function registration for 0D diagnostics at the nearest grid-center point for a single grid.

This module constructs and registers the C helper routine "diagnostics_nearest_grid_center".
The generated C function samples gridfunction data without interpolation at the single grid
point nearest the physical grid center for the specified grid. The caller is responsible for
looping over grids. One line is appended per timestep to a persistent per-grid output file
whose name encodes the coordinate system, grid number, and convergence factor.

Function
--------
register_CFunction_diagnostics_nearest_grid_center
    Construct and register the "diagnostics_nearest_grid_center" C function.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast, Tuple

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg



def get_center_index_exprs_for_coordsystem(CoordSystem: str) -> Tuple[str, str, str]:
    """
    Return C-code string expressions for the grid indices nearest the domain center.

    The returned expressions are intended to be embedded directly in generated C/C++ code.
    Convention:
      * Cartesian: (i0,i1,i2) = (mid, mid, mid)
      * Spherical/Cylindrical/SymTP: (i0,i1,i2) = (radial_min, mid, mid), with radial_min = NGHOSTS.

    :param CoordSystem: Coordinate system name (e.g., "Cartesian", "SinhSpherical", "Cylindrical", "SymTP").
    :return: Tuple (i0_center_expr, i1_center_expr, i2_center_expr) as strings.
    :raises ValueError: If CoordSystem is unsupported.
    """
      
    is_spherical_family = ("Spherical" in CoordSystem) or ("Cylindrical" in CoordSystem) or ("SymTP" in CoordSystem)
    is_cartesian = "Cartesian" in CoordSystem
    
    # Choose expressions for the nearest-to-center grid indices based on CoordSystem.
    # Cartesian: (i0,i1,i2) = (mid,mid,mid)
    # Spherical/Cylindrical/SymTP: (i0,i1,i2) = (radial_min, mid, mid)    
    if is_spherical_family:
        i0_center_expr = "NGHOSTS"
    elif is_cartesian:
        i0_center_expr = "params->Nxx_plus_2NGHOSTS0/2"
    else:
        raise ValueError(f"Unsupported CoordSystem: {CoordSystem}")

    i1_center_expr = "params->Nxx_plus_2NGHOSTS1/2"
    i2_center_expr = "params->Nxx_plus_2NGHOSTS2/2"
    return i0_center_expr, i1_center_expr, i2_center_expr



def register_CFunction_diagnostics_nearest_grid_center(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register a C function that writes 0D diagnostics by sampling (no interpolation).

    This function generates and registers the C helper "diagnostics_nearest_grid_center", which
    mirrors the interpolation-based diagnostics API so it can serve as a drop-in replacement in
    calling code that loops over grids. For the selected grid, the generated C code opens a
    persistent per-grid output file, writes a header on the first timestep, and on every timestep
    appends a single row containing the simulation time followed by sampled gridfunction values.
    The nearest-to-center index selection is specialized based on the coordinate system: for
    Cartesian-like systems the mid index is used in all directions; for Spherical, Cylindrical,
    and SymTP the radial index is set to NGHOSTS and mid indices are used in the angular directions.

    :param CoordSystem: Name of the coordinate system used to specialize index selection in the
                        generated C function and its wrapper.
    :raises ValueError: If the coordinate system is not supported by this helper.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    TBD
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
            
    i0_center_expr, i1_center_expr, i2_center_expr = get_center_index_exprs_for_coordsystem(CoordSystem)    

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "diagnostics/diagnostics_nearest_common.h",
    ]

    desc = """
 * @file diagnostics_nearest_grid_center.c
 * @brief Write 0D diagnostics for a single grid by sampling, without interpolation, the grid
 *        point nearest the physical grid center, and append one row per timestep to a persistent file.
 *
 * The function "diagnostics_nearest_grid_center" appends diagnostics for a single grid to a per-grid
 * text file whose name encodes the runtime coordinate system, grid number, and convergence factor:
 *
 *   out0d-<CoordSystemName>-grid<XX>-conv_factor-<CF>.txt
 *
 * On the first timestep, a header is written. On all timesteps, one row is appended containing the
 * simulation time followed by sampled gridfunction values. The nearest-to-center indices are chosen
 * based on the coordinate system: for Cartesian-like systems, the midpoint index is chosen in each
 * dimension; for Spherical, Cylindrical, and SymTP, the radial index is set to NGHOSTS and midpoint
 * indices are chosen in the angular directions. The "xx" coordinate arrays are accepted for API
 * compatibility but are not used by this routine.
 *
 * If a user-editable block is provided in the implementation, users may add custom logic such as
 * additional columns or filtering before rows are written.* 
 *
 * @return     void.
 """
    cfunc_type = "void"
    name = "diagnostics_nearest_grid_center"
    params = """commondata_struct *restrict commondata, const int grid, const params_struct *restrict params,
                                                         const params_struct *restrict params_chare,
                                                         const REAL *restrict xx[3], const int NUM_GFS_NEAREST, const int which_gfs[],
                                                         const char **diagnostic_gf_names, const REAL *restrict gridfuncs_diags[],
                                                         const int chare_index[3]"""

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


const int Nxx0chare = params_chare->Nxx0;
const int Nxx1chare = params_chare->Nxx1;
const int Nxx2chare = params_chare->Nxx2;
    
const int i0_center_local = MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0_center, Nxx0chare);
const int i1_center_local = MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1_center, Nxx1chare);
const int i2_center_local = MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2_center, Nxx2chare);

const int idx3 = IDX3P(params_chare, i0_center_local, i1_center_local, i2_center_local);

// Active grid data pointer.
const REAL *restrict src = gridfuncs_diags[grid];

// Emit a single output row: [time, sampled values...]
const int NUM_COLS = 1 + NUM_GFS_NEAREST;
REAL row[NUM_COLS];
row[0] = commondata->time;
for (int ii = 0; ii < NUM_GFS_NEAREST; ii++) {{
  const int gf = which_gfs[ii];
  row[1 + ii] = src[IDX4Ppt(params_chare, gf, idx3)];
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
