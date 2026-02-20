"""
C function registration for 2D diagnostics on the nearest xy and yz planes for a single grid.

This module constructs and registers the C helper routine "diagnostics_nearest_2d_xy_and_yz_planes".
The generated C function samples gridfunction data without interpolation at interior points on the
planes nearest to the xy and yz coordinate planes for the selected grid. For each call and grid,
two per-timestep files are written: one for the xy plane and one for the yz plane. The caller is
responsible for looping over grids.

Functions
---------
register_CFunction_diagnostics_nearest_2d_xy_and_yz_planes
    Construct and register the "diagnostics_nearest_2d_xy_and_yz_planes" C function.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Any, Dict, Sequence, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def bhah_plane_configs() -> Dict[str, Dict[str, Dict[str, Any]]]:
    """
    Return the plane_configs dict used to construct the 2D xy/yz plane sampling.

    Schema:
      plane_configs[PLANE][FAMILY] -> CONFIG, where:
        PLANE  : "xy" or "yz"
        FAMILY : one of {"Cartesian","Cylindrical","Spherical","SymTP","Wedge"}
        CONFIG : dict with keys:
            "fixed_dim" : str
            "fixed_val" : str or list[str]
            "loop_dims" : list[str]

    :return: Dictionary mapping PLANE ("xy" or "yz") and coordinate FAMILY to a CONFIG dict.
    """
    # -----------------------------------------------------------------------------
    # plane_configs SCHEMA (dict of dicts):
    #
    # plane_configs[PLANE][FAMILY] -> CONFIG, where:
    #   PLANE  : "xy" or "yz"
    #   FAMILY : one of {"Cartesian","Cylindrical","Spherical","SymTP","Wedge"}
    #   CONFIG : dict with keys:
    #       "fixed_dim" : str
    #           Which grid index is held fixed to realize the requested plane.
    #           One of {"i0","i1","i2"}.
    #       "fixed_val" : str or list[str]
    #           Value(s) assigned to fixed_dim. Each value is the name of a C local
    #           (defined later in the generated C body) that evaluates to an integer index,
    #           e.g. "i2_mid", "i2_q1", "i2_q3".
    #           If this is a list, multiple slices are emitted for that plane
    #           (for example two phi slices at approximately +/- pi/2 to capture x=0).
    #       "loop_dims" : list[str] of length 2
    #           The outer and inner loop indices that sweep the plane interior.
    #           Example: ["i1","i0"] means i1 is the outer loop and i0 the inner loop.
    #
    # Geometry notes:
    # - xy plane  -> typically fix z (Cartesian, Cylindrical) or fix theta to mid (Spherical, SymTP).
    # - yz plane  -> typically fix x (Cartesian) or choose phi near +/- pi/2 (Cylindrical, Spherical, SymTP).
    # - Wedge     -> limited angular extent; we choose mid or quarter indices as reasonable slices.
    # -----------------------------------------------------------------------------
    plane_configs = {
        "xy": {  # 2D plane with output columns "x y" (z or equivalent is fixed)
            "Cartesian": {
                # Fix z index to domain midpoint -> z approx 0 slice
                "fixed_dim": "i2",
                "fixed_val": "i2_mid",
                # Sweep over y (outer) and x (inner)
                "loop_dims": ["i1", "i0"],
            },
            "Cylindrical": {
                # Fix z at mid-plane; in (rho, phi, z) this makes a horizontal slice
                "fixed_dim": "i2",
                "fixed_val": "i2_mid",
                # Sweep phi (outer) and rho (inner) at fixed z
                "loop_dims": ["i1", "i0"],
            },
            "Spherical": {
                # Fix polar angle theta to equator (theta approx pi/2) -> nearest xy-plane
                "fixed_dim": "i1",
                "fixed_val": "i1_mid",
                # Sweep azimuth phi (outer) and radius r (inner)
                "loop_dims": ["i2", "i0"],
            },
            "SymTP": {
                # Symmetric toroidal/poloidal system: fix the polar-like angle at mid
                # to mimic theta = pi/2 (equatorial) -> nearest xy-plane
                "fixed_dim": "i1",
                "fixed_val": "i1_mid",
                # Sweep second angle (outer) and radial-like coordinate (inner)
                "loop_dims": ["i2", "i0"],
            },
            "Wedge": {
                # In wedge geometries, use two quarter-plane z-like slices to represent the xy-like plane
                "fixed_dim": "i2",
                "fixed_val": ["i2_q1", "i2_q3"],  # multi-slice at quarter indices
                "loop_dims": ["i1", "i0"],
            },
        },
        "yz": {  # 2D plane with output columns "y z" (x or equivalent is fixed)
            "Cartesian": {
                # Fix x index to domain midpoint -> x approx 0 slice (yz-plane)
                "fixed_dim": "i0",
                "fixed_val": "i0_mid",
                # Sweep z (outer) and y (inner)
                "loop_dims": ["i2", "i1"],
            },
            "Cylindrical": {
                # yz-plane corresponds to x=0 -> choose phi near +/- pi/2
                "fixed_dim": "i1",
                "fixed_val": ["i1_q1", "i1_q3"],  # phi approx +pi/2 and -pi/2
                # Sweep z (outer) and rho (inner)
                "loop_dims": ["i2", "i0"],
            },
            "Spherical": {
                # yz-plane (x=0) -> choose phi near +/- pi/2 at all r,theta
                "fixed_dim": "i2",
                "fixed_val": ["i2_q1", "i2_q3"],  # phi approx +pi/2 and -pi/2
                # Sweep theta (outer) and radius r (inner)
                "loop_dims": ["i1", "i0"],
            },
            "SymTP": {
                # Same reasoning as Spherical: choose azimuth-like angle near +/- pi/2
                "fixed_dim": "i2",
                "fixed_val": ["i2_q1", "i2_q3"],
                "loop_dims": ["i1", "i0"],
            },
            "Wedge": {
                # In a wedge, hold the across-wedge index at mid to approximate x=0
                "fixed_dim": "i1",
                "fixed_val": "i1_mid",
                # Sweep z-like (outer) and radial-like (inner)
                "loop_dims": ["i2", "i0"],
            },
        },
    }
    return plane_configs


def get_coord_family(cs: str) -> str:
    """
    Return the coordinate system family name from the full CoordSystem string.

    :param cs: The full CoordSystem string.
    :raises ValueError: If the coordinate system is not supported.
    :return: The coordinate family name.
    """
    plane_configs = bhah_plane_configs()
    for family in plane_configs["xy"]:
        if family in cs:
            return family
    raise ValueError(f"Unsupported CoordSystem: {cs}")


def generate_plane_loop_code(config: Dict[str, Sequence[str]], plane: str) -> str:
    """
    Generate a C code block for sampling a plane based on its config.

    :param config: A CONFIG as documented above (fixed_dim, fixed_val, loop_dims).
    :param plane: Either "xy" or "yz"; determines output mapping and file handle.
    :return: The generated C code string for the plane loop.
    """
    loop_dims, fixed_dim, fixed_val = (
        config["loop_dims"],
        config["fixed_dim"],
        config["fixed_val"],
    )

    # Map written coordinates and output file by plane:
    # xy -> write x,y into row[0],row[1] and use out_xy
    # yz -> write y, z into row[0],row[1] and use out_yz
    row_map, num_coords, out_file = (
        ("row[0]=xCart[0]; row[1]=xCart[1];", 2, "out_xy")
        if plane == "xy"
        else ("row[0]=xCart[1]; row[1]=xCart[2];", 2, "out_yz")
    )

    # Innermost loop body: transform coords, sample gfs, write row
    inner_body = f"""
const int idx3 = IDX3P(params, i0, i1, i2);
REAL xCart[3], xOrig[3] = {{xx[0][i0], xx[1][i1], xx[2][i2]}};
xx_to_Cart(params, xOrig, xCart);
{row_map}
for (int gf_idx = 0; gf_idx < NUM_GFS_NEAREST; gf_idx++) {{
  const int gf = which_gfs[gf_idx];
  row[{num_coords} + gf_idx] = src[IDX4Ppt(params, gf, idx3)];
}} // END LOOP over gridfunctions
diag_write_row({out_file}, {num_coords} + NUM_GFS_NEAREST, row);
"""

    # Build nested loops. The first loop_dims element is intended to be outermost.
    # We reverse when wrapping so the first becomes the outer loop.
    loop_code = inner_body
    for dim in reversed(loop_dims):
        loop_code = f"for (int {dim}=NGHOSTS; {dim}<{dim}_end; {dim}++) {{{loop_code}}} // END LOOP over {dim}\n"

    # Wrap with fixed_dim logic. If fixed_val is a list, emit multiple slices.
    if isinstance(fixed_val, list):
        slices = f"const int {fixed_dim}_slices[{len(fixed_val)}] = {{{', '.join(fixed_val)}}};"
        return f"""
{{
  {slices}
  for (int slice = 0; slice < {len(fixed_val)}; slice++) {{
    const int {fixed_dim} = {fixed_dim}_slices[slice];
    {loop_code}
  }} // END LOOP over slices
}} // END BLOCK {plane}-plane output"""
    return f"{{ const int {fixed_dim} = {fixed_val}; {loop_code} }} // END BLOCK {plane}-plane output"


def register_CFunction_diagnostics_nearest_2d_xy_and_yz_planes(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register a C function that writes 2D diagnostics on the nearest xy and yz planes for a single grid.

    This function generates and registers the C helper "diagnostics_nearest_2d_xy_and_yz_planes".
    The generated C code identifies interior slices nearest to the xy and yz planes based on the
    provided coordinate system, converts native grid coordinates to Cartesian (x, y, z) for output,
    and for each call writes two per-timestep files: one for the xy plane and one for the yz plane.
    Sampling is performed without interpolation; values are read directly from the selected grid's
    diagnostic gridfunctions. The plane selection and loop ordering are configured via a small
    data-driven table for Cartesian, Cylindrical, Spherical, SymTP, and Wedge families; some
    families emit multiple slices for a plane (for example, two phi slices near +/- pi/2 to realize x=0).

    :param CoordSystem: Name of the coordinate system family that specializes plane selection and the wrapper.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    TBD
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    family = get_coord_family(CoordSystem)
    plane_configs = bhah_plane_configs()
    xy_plane_code = generate_plane_loop_code(plane_configs["xy"][family], "xy")
    yz_plane_code = generate_plane_loop_code(plane_configs["yz"][family], "yz")

    includes = [
        "stdlib.h",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "diagnostics/diagnostics_nearest_common.h",
    ]
    desc = f"""
 * @file diagnostics_nearest_2d_xy_and_yz_planes.c
 * @brief Sample and write 2D diagnostics on the nearest xy and yz planes for a single {CoordSystem} grid.
 *
 * Overview:
 * For the specified grid at the current time, this routine:
 *  - Locates interior slices that are nearest to the xy and yz planes, with selection rules
 *    specialized by the runtime coordinate system.
 *  - Converts native coordinates xx to Cartesian coordinates (x, y, z) using xx_to_Cart so that
 *    the first columns of each row contain mapped coordinates.
 *  - Writes two files per call, one for the xy plane and one for the yz plane. Each file begins
 *    with a time comment and a header, followed by one row per interior point that contains the
 *    mapped coordinates and sampled diagnostic values.
 *  - Performs sampling without interpolation; values are read directly from gridfuncs_diags[grid].
 *
 * Plane selection notes (examples, not exhaustive):
 *  - Cartesian: xy fixes i2 at mid; yz fixes i0 at mid.
 *  - Cylindrical: xy fixes z at mid; yz emits two phi slices near +/- pi/2 to realize x=0.
 *  - Spherical and SymTP: xy fixes the polar-like angle at mid; yz emits two phi-like slices near +/- pi/2.
 *  - Wedge: xy may emit two z-like slices at quarter indices; yz fixes the across-wedge index at mid.
 *
 * If a user-editable block is provided in the implementation, users may insert custom logic such as
 * adding extra columns or filtering before rows are written.
 *
 * @param[in,out] commondata            Pointer to common runtime data used for time and I/O.
 * @param[in]     grid                  Grid index to process.
 * @param[in]     params                Pointer to simulation and grid parameters (sizes, names, strides).
 * @param[in]     xx                    Native grid coordinates; xx[d][i_d] gives the coordinate along dimension d.
 * @param[in]     NUM_GFS_NEAREST       Number of diagnostic gridfunctions to sample at each interior point.
 * @param[in]     which_gfs             Array of length NUM_GFS_NEAREST specifying which gridfunctions to sample.
 * @param[in]     diagnostic_gf_names   Array of length NUM_GFS_NEAREST with human-readable names for headers.
 * @param[in]     gridfuncs_diags       Array of pointers; gridfuncs_diags[grid] points to this grid's diagnostic data.
 *
 * @return        void                  No return value. On success two text files are written and closed. Fatal I/O
 *                                      or allocation failures result in program termination.
"""
    cfunc_type = "void"
    name = "diagnostics_nearest_2d_xy_and_yz_planes"
    params = """commondata_struct *restrict commondata, const int grid, const params_struct *restrict params,
                const REAL *restrict xx[3], const int NUM_GFS_NEAREST, const int which_gfs[],
                const char **diagnostic_gf_names, const REAL *restrict gridfuncs_diags[]"""

    body = rf"""
  // Build filename component with runtime coordinate system name and grid number
  char coordsys_with_grid[128];
  snprintf(coordsys_with_grid, sizeof(coordsys_with_grid), "%s-grid%02d", params->CoordSystemName, grid);

  // Open output files (one file per timestep, per plane, for this grid)
  FILE *out_xy = open_outfile("out2d-xy", coordsys_with_grid, commondata, /*include_time=*/1);
  FILE *out_yz = open_outfile("out2d-yz", coordsys_with_grid, commondata, /*include_time=*/1);

  if (!out_xy || !out_yz) {{
    if (out_xy)
      fclose(out_xy);
    if (out_yz)
      fclose(out_yz);
    fprintf(stderr, "Error: Cannot open output files for grid %d.\n", grid);
    exit(1);
  }} // END IF cannot open output files

  // Write time comment and headers
  diag_write_time_comment(out_xy, commondata->time);
  diag_write_header(out_xy, "x y", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
  diag_write_time_comment(out_yz, commondata->time);
  diag_write_header(out_yz, "y z", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);

  // Active grid data pointer and reusable row buffer
  const REAL *restrict src = gridfuncs_diags[grid];
  const int NUM_COLS = 2 + NUM_GFS_NEAREST;
  REAL *row = (REAL *)malloc(sizeof(REAL) * (size_t)NUM_COLS);
  if (!row) {{
    fprintf(stderr, "Error: Failed to allocate memory for row buffer.\n");
    exit(1);
  }} // END IF row allocation failure

  // Interior grid counts and loop bounds
  MAYBE_UNUSED const int N0int = params->Nxx_plus_2NGHOSTS0 - 2 * NGHOSTS;
  MAYBE_UNUSED const int N1int = params->Nxx_plus_2NGHOSTS1 - 2 * NGHOSTS;
  MAYBE_UNUSED const int N2int = params->Nxx_plus_2NGHOSTS2 - 2 * NGHOSTS;
  const int i0_end = params->Nxx_plus_2NGHOSTS0 - NGHOSTS;
  const int i1_end = params->Nxx_plus_2NGHOSTS1 - NGHOSTS;
  const int i2_end = params->Nxx_plus_2NGHOSTS2 - NGHOSTS;

  // Fixed-point index helpers
  MAYBE_UNUSED const int i0_mid = params->Nxx_plus_2NGHOSTS0 / 2;
  MAYBE_UNUSED const int i1_mid = params->Nxx_plus_2NGHOSTS1 / 2;
  MAYBE_UNUSED const int i2_mid = params->Nxx_plus_2NGHOSTS2 / 2;
  MAYBE_UNUSED const int i1_q1 = (int)(NGHOSTS + 0.25 * (REAL)N1int - 0.5);
  MAYBE_UNUSED const int i1_q3 = (int)(NGHOSTS + 0.75 * (REAL)N1int - 0.5);
  MAYBE_UNUSED const int i2_q1 = (int)(NGHOSTS + 0.25 * (REAL)N2int - 0.5);
  MAYBE_UNUSED const int i2_q3 = (int)(NGHOSTS + 0.75 * (REAL)N2int - 0.5);

  // --- Sample and write data for the xy-plane for {CoordSystem} ---
  {xy_plane_code}

  // --- Sample and write data for the yz-plane for {CoordSystem} ---
  {yz_plane_code}

  // Finalize
  free(row);
  fclose(out_xy);
  fclose(out_yz);
"""

    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
        CoordSystem_for_wrapper_func=CoordSystem,
    )
    return pcg.NRPyEnv()
