"""
C function registration for 1D diagnostics along the nearest y and z axes for a single grid.

This module constructs and registers the C helper routine "diagnostics_nearest_1d_y_and_z_axes".
The generated C function samples gridfunction data without interpolation along axis-aligned
lines that are nearest to the physical y- and z-axes for the specified grid. The caller is
responsible for looping over grids. For each timestep, rows are appended to two persistent
per-grid output files (one for y and one for z) whose names encode the coordinate system,
grid number, and convergence factor. Each file begins with a one-line time comment and an
axis-specific header, followed by rows whose first column is the axis coordinate.

Function
--------
register_CFunction_diagnostics_nearest_1d_y_and_z_axes
    Construct and register the "diagnostics_nearest_1d_y_and_z_axes" C function.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Any, Dict, List, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def bhah_family_from_coord(CoordSystem: str) -> str:
    """
    Return the coordinate *family* string given a CoordSystem name.
    Shared between BHaH and superB diagnostics.

    :param CoordSystem: Name of the coordinate system (e.g., Cartesian, Spherical, Cylindrical, SymTP,
                        Wedge, or Spherical_Ring variants).
    :return: Coordinate family string (e.g., "Cartesian", "Spherical", "Cylindrical", "SymTP", "Wedge",
             or "Spherical_Ring").
    :raises ValueError: If CoordSystem does not match a supported family.
    """
    if ("Spherical" in CoordSystem) and ("Ring" in CoordSystem):
        return "Spherical_Ring"
    for fam in ("Cartesian", "Spherical", "Cylindrical", "SymTP", "Wedge"):
        if fam in CoordSystem:
            return fam
    raise ValueError(f"Unsupported CoordSystem: {CoordSystem}")


def bhah_axis_configs() -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
    """
    Return axis_configs dict used to construct the 1D y/z lines.

    Structure:
      axis_configs[AXIS][FAMILY] -> list[LINE_SPEC]

    :return: Dictionary mapping AXIS ("y" or "z") and FAMILY to a list of LINE_SPEC dictionaries.
    """
    # -------------------------------------------------------------------------
    # axis_configs schema (dict of dicts):
    # axis_configs[AXIS][FAMILY] -> list[LINE_SPEC], where:
    #   AXIS   : "y" or "z" (the physical axis we will output as the first column)
    #   FAMILY : one of {"Cartesian","Spherical","Cylindrical","SymTP","Wedge","Spherical_Ring"}
    #   LINE_SPEC is a dict with:
    #       "vary" : str      -> which grid index to iterate over ("i0", "i1", or "i2")
    #       "fixed": dict     -> all *other* indices held fixed at named helper indices
    #            keys   are "i0","i1","i2"
    #            values are names of C locals defined in `body` (e.g. "i0_mid","i2_q1",...)
    #
    # Notes:
    # - We always convert (xx -> Cartesian) and then take xCart[1] for y, xCart[2] for z.
    # - The specific fixed-point choices differ by coordinate family to pick the line
    #   that is geometrically "nearest" to the physical axis without interpolation.
    # - An empty list means: no meaningful axis line exists for that (AXIS, FAMILY).
    # -------------------------------------------------------------------------
    axis_configs = {
        "y": {  # Lines used to sample along the physical y-axis (output column 0 = y)
            "Cartesian": [
                {
                    "vary": "i1",  # march along y-index to trace a y-line
                    "fixed": {
                        "i0": "i0_mid",  # hold x at domain midpoint
                        "i2": "i2_mid",  # hold z at domain midpoint
                    },
                },
            ],
            "Spherical": [
                {
                    "vary": "i0",  # march in radius to sample a radial cut at a fixed (theta,phi)
                    "fixed": {
                        "i1": "i1_mid",  # theta ~ pi/2 (equator) for nearest y-line (Q1 longitude)
                        "i2": "i2_q1",  # phi ~ +pi/2 (quadrant 1)
                    },
                },
                {
                    "vary": "i0",
                    "fixed": {
                        "i1": "i1_mid",  # theta ~ pi/2 (equator)
                        "i2": "i2_q3",  # phi ~ -pi/2 (quadrant 3) — the opposite y-direction
                    },
                },
            ],
            "Cylindrical": [
                {
                    "vary": "i0",  # march in rho at fixed (phi,z) to follow a line aligned with y
                    "fixed": {
                        "i1": "i1_q1",  # phi ~ +pi/2 (quadrant 1) -> points toward +y
                        "i2": "i2_mid",  # z mid-plane
                    },
                },
                {
                    "vary": "i0",
                    "fixed": {
                        "i1": "i1_q3",  # phi ~ -pi/2 (quadrant 3) -> points toward -y
                        "i2": "i2_mid",  # z mid-plane
                    },
                },
            ],
            "SymTP": [
                {
                    "vary": "i0",  # march in radial-like coordinate at fixed angles
                    "fixed": {
                        "i1": "i1_mid",  # mid in the first angular coordinate
                        "i2": "i2_q1",  # second angle set to quadrant 1 (y>0)
                    },
                },
                {
                    "vary": "i0",
                    "fixed": {
                        "i1": "i1_mid",
                        "i2": "i2_q3",  # second angle set to quadrant 3 (y<0)
                    },
                },
            ],
            "Wedge": [
                # No well-defined "nearest y-axis" line in wedge geometry -> no samples
            ],
            "Spherical_Ring": [
                {
                    "vary": "i0",  # march in radius within the ring
                    "fixed": {
                        "i1": "i1_mid",  # polar angle at mid (equatorial-like)
                        "i2": "i2_q1",  # azimuth at quadrant 1 (y>0)
                    },
                },
                {
                    "vary": "i0",
                    "fixed": {
                        "i1": "i1_mid",
                        "i2": "i2_q3",  # azimuth at quadrant 3 (y<0)
                    },
                },
            ],
        },
        "z": {  # Lines used to sample along the physical z-axis (output column 0 = z)
            "Cartesian": [
                {
                    "vary": "i2",  # march along z-index to trace a z-line
                    "fixed": {
                        "i0": "i0_mid",  # hold x at domain midpoint
                        "i1": "i1_mid",  # hold y at domain midpoint
                    },
                },
            ],
            "Spherical": [
                {
                    "vary": "i0",  # march in radius at north pole direction
                    "fixed": {
                        "i1": "i1_min",  # theta ~ 0   (north pole)
                        "i2": "i2_min",  # arbitrary reference phi for a single pole line
                    },
                },
                {
                    "vary": "i0",  # march in radius at south pole direction
                    "fixed": {
                        "i1": "i1_max",  # theta ~ pi  (south pole)
                        "i2": "i2_min",  # same phi reference
                    },
                },
            ],
            "Cylindrical": [
                {
                    "vary": "i2",  # march along z at rho=0, phi=-pi (on the axis)
                    "fixed": {
                        "i0": "i0_rmin",  # rho = 0 -> cylindrical axis (the z-axis)
                        "i1": "i1_pmin",  # phi = -pi (any phi is equivalent at rho=0)
                    },
                },
            ],
            "SymTP": [
                {
                    "vary": "i0",  # march in radial-like coordinate at one pole
                    "fixed": {
                        "i1": "i1_min",  # first angle min -> one pole (z>0)
                        "i2": "i2_min",  # second angle ref
                    },
                },
                {
                    "vary": "i0",  # march in radial-like coordinate at the opposite pole
                    "fixed": {
                        "i1": "i1_max",  # first angle max -> opposite pole (z<0)
                        "i2": "i2_min",  # second angle ref
                    },
                },
            ],
            "Wedge": [
                {
                    "vary": "i0",  # march in the radial-like direction that aligns with wedge's z
                    "fixed": {
                        "i1": "i1_mid",  # hold across-wedge index mid
                        "i2": "i2_mid",  # hold the other index mid
                    },
                },
            ],
            "Spherical_Ring": [
                # No unique z-axis line within a pure spherical ring configuration -> no samples
            ],
        },
    }
    return axis_configs


def count_expr(lines: List[Dict[str, Any]]) -> str:
    """
    Construct a C expression for the total number of interior sample points implied by `lines`.

    This helper maps each LINE_SPEC's varying index ("i0", "i1", "i2") to the corresponding interior
    count symbol ("N0int", "N1int", "N1int") and returns a C expression that sums these counts.

    :param lines: List of LINE_SPEC dictionaries, each containing a "vary" key.
    :return: C expression string for the total point count (e.g., "N0int + N2int" or "0").
    """
    # Map each index name to the interior count symbol used in C.
    # keys:   "i0","i1","i2" (grid indices)
    # values: "N0int","N1int","N2int" (C locals representing interior sizes)
    counts = {"i0": "N0int", "i1": "N1int", "i2": "N2int"}
    # For each LINE_SPEC, the number of points equals the interior length
    # of the varying dimension; sum across all lines to get the max buffer size.
    terms = [counts[line["vary"]] for line in lines]
    return " + ".join(terms) if terms else "0"


def gen_fill(axis_char: str, lines: List[Dict[str, Any]]) -> str:
    """
    Construct C code that fills (coord, idx3) samples for the requested physical axis from `lines`.

    The generated code loops over each LINE_SPEC, binds fixed indices, converts xx -> Cartesian, and
    stores the selected axis coordinate (y or z) plus idx3 into the appropriate buffer and counter.
    If `lines` is empty, a no-op stub is returned.

    :param axis_char: Physical axis selector ("y" or "z").
    :param lines: List of LINE_SPEC dictionaries describing the axis-line samples.
    :return: C code string that fills the per-axis sample buffer.
    """
    if not lines:
        return "(void)xx;  // no points for this axis/family\n"

    # coord_idx selects the physical axis to store in data_point_1d_struct.coord:
    #   "y" -> xCart[1], "z" -> xCart[2]
    coord_idx = "1" if axis_char == "y" else "2"

    # Names of the C arrays/counters we fill for each axis
    arr = "data_points_y" if axis_char == "y" else "data_points_z"
    cnt = "count_y" if axis_char == "y" else "count_z"

    # Map Python-side index labels to dimension numbers used in params->Nxx_plus_2NGHOSTS{d}
    # i0 -> dim 0, i1 -> dim 1, i2 -> dim 2
    dim_index = {"i0": "0", "i1": "1", "i2": "2"}

    blocks = []
    for line in lines:
        vary = line["vary"]  # e.g., "i0"
        vary_dim = dim_index[vary]  # e.g., "0"
        fixed = line["fixed"]  # dict like {"i1":"i1_mid","i2":"i2_q1"}
        blocks.append(
            "\n".join(
                [
                    # loop bounds use the correct dimension number from vary_dim
                    f"for (int {vary} = NGHOSTS; {vary} < params->Nxx_plus_2NGHOSTS{vary_dim} - NGHOSTS; {vary}++) {{",
                    *(
                        # For each non-varying index, bind that index to the named constant
                        f"  const int {i} = {fixed[i]};"
                        for i in ("i0", "i1", "i2")
                        if i != vary
                    ),
                    "  const int idx3 = IDX3P(params, i0, i1, i2);",
                    "  REAL xCart[3], xOrig[3] = { xx[0][i0], xx[1][i1], xx[2][i2] };",
                    "  xx_to_Cart(params, xOrig, xCart);",
                    # Save the physical coordinate (y or z) and the flat 3D index
                    f"  {arr}[{cnt}].coord = xCart[{coord_idx}];",
                    f"  {arr}[{cnt}].idx3  = idx3;",
                    f"  {cnt}++;",
                    f"}} // END LOOP over {vary}",
                ]
            )
        )
    return "\n".join(blocks) + "\n"


def register_CFunction_diagnostics_nearest_1d_y_and_z_axes(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct & register a C helper for 1D diagnostics along the nearest points to the y- and z-axes.

    This function generates and registers the C helper "diagnostics_nearest_1d_y_and_z_axes", which
    mirrors the interpolation-based diagnostics API so it can serve as a drop-in replacement in
    calling code that loops over grids. The generated C code (a) selects axis-line samples according
    to the coordinate family implied by `CoordSystem`, (b) converts logical coordinates to Cartesian,
    (c) buffers (axis_coord, idx3) pairs, (d) sorts by the physical axis coordinate, and (e) streams
    rows whose first column is the axis coordinate followed by selected diagnostic gridfunction values.
    Two persistent per-grid output files are produced per call, one for the y-axis and one for the z-axis.
    Each file begins with a one-line time comment and an axis-specific header. Sampling is performed at
    grid points only; no interpolation is used.

    :param CoordSystem: Name of the coordinate system used to specialize axis-line selection and
                        wrapper generation (e.g., Cartesian, Spherical, Cylindrical, SymTP, Wedge,
                        Spherical_Ring).
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    TBD
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    fam = bhah_family_from_coord(CoordSystem)
    axis_configs = bhah_axis_configs()

    y_lines = axis_configs["y"][fam]
    z_lines = axis_configs["z"][fam]

    y_count_expr = count_expr(y_lines)
    z_count_expr = count_expr(z_lines)

    fill_y = gen_fill("y", y_lines)
    fill_z = gen_fill("z", z_lines)

    prefunc = r"""
// Data point for sorting by physical axis coordinate
typedef struct {
  REAL coord; // physical y or z
  int idx3;   // 3D index
} data_point_1d_struct;

// qsort comparator (file scope to satisfy -std=c11)
/**
 * @brief Compare two data_point_1d_struct items by their coord value.
 * @param a Pointer to the left-hand data_point_1d_struct.
 * @param b Pointer to the right-hand data_point_1d_struct.
 * @return Negative if a < b, zero if equal, positive if a > b.
 */
static int compare_by_coord(const void *a, const void *b) {
  const REAL lv = ((const data_point_1d_struct *)a)->coord;
  const REAL rv = ((const data_point_1d_struct *)b)->coord;
  return (lv > rv) - (lv < rv);
} // END FUNCTION compare_by_coord
"""

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "diagnostics/diagnostics_nearest_common.h",
    ]

    desc = """
 * @file diagnostics_nearest_1d_y_and_z_axes.c
 * @brief Write 1D diagnostics for a single grid by sampling, without interpolation, axis-aligned
 *        lines nearest to the physical y and z axes, and append rows per timestep to persistent files.
 *
 * The function "diagnostics_nearest_1d_y_and_z_axes" appends diagnostics for a single grid to two
 * per-grid text files whose names encode the runtime coordinate system, grid number, and convergence
 * factor:
 *
 *   out1d-y-<CoordSystemName>-grid<XX>-conv_factor-<CF>.txt
 *   out1d-z-<CoordSystemName>-grid<XX>-conv_factor-<CF>.txt
 *
 * For each axis, the routine:
 *   1) Selects one or more index-line samples based on the coordinate family (e.g., Cartesian,
 *      Spherical, Cylindrical, SymTP, Wedge, Spherical_Ring).
 *   2) Converts logical coordinates xx to Cartesian via xx_to_Cart and extracts y (xCart[1]) or
 *      z (xCart[2]) as the axis coordinate.
 *   3) Buffers (axis_coord, idx3) pairs, then sorts them in ascending order using qsort.
 *   4) Writes a single-line time comment, an axis-specific header ("y" or "z"), and streams rows:
 *      [axis_coord, values of selected diagnostic gridfunctions].
 *
 * Gridfunction values are loaded from the flattened diagnostic array using IDX4Ppt with a 3D index
 * constructed by IDX3P. Sampling occurs at grid points only; no interpolation is performed.
 * On allocation or file-open failure, the routine prints an error message to stderr and terminates.
 * If a user-editable block is provided in the implementation, users may add custom logic such as
 * additional columns or filtering before rows are written.
 *
 * @param[in]  commondata           Pointer to global simulation metadata (e.g., time and step counters).
 * @param[in]  grid                 Zero-based grid index used for selecting data and for file naming.
 * @param[in]  params               Pointer to per-grid parameters (sizes, ghost zones, strides, names).
 * @param[in]  xx                   Array of 3 pointers to logical coordinates per dimension; used for
 *                                  coordinate conversion to Cartesian space.
 * @param[in]  NUM_GFS_NEAREST      Number of gridfunctions to include per output row.
 * @param[in]  which_gfs            Array of length NUM_GFS_NEAREST identifying the gridfunction indices.
 * @param[in]  diagnostic_gf_names  Array of NUM_GFS_NEAREST C strings used to construct column headers.
 * @param[in]  gridfuncs_diags      Array of per-grid pointers to flattened diagnostic data; this routine
 *                                  reads from gridfuncs_diags[grid].
 *
 * @return     void
"""

    cfunc_type = "void"
    name = "diagnostics_nearest_1d_y_and_z_axes"
    params = """commondata_struct *restrict commondata, const int grid, const params_struct *restrict params,
                const REAL *restrict xx[3], const int NUM_GFS_NEAREST, const int which_gfs[],
                const char **diagnostic_gf_names, const REAL *restrict gridfuncs_diags[]"""

    body = rf"""
  // Interior counts
  MAYBE_UNUSED const int N0int = params->Nxx_plus_2NGHOSTS0 - 2 * NGHOSTS;
  MAYBE_UNUSED const int N1int = params->Nxx_plus_2NGHOSTS1 - 2 * NGHOSTS;
  MAYBE_UNUSED const int N2int = params->Nxx_plus_2NGHOSTS2 - 2 * NGHOSTS;

  // Common fixed-point helpers
  MAYBE_UNUSED const int i0_mid = params->Nxx_plus_2NGHOSTS0 / 2;
  MAYBE_UNUSED const int i1_mid = params->Nxx_plus_2NGHOSTS1 / 2;
  MAYBE_UNUSED const int i2_mid = params->Nxx_plus_2NGHOSTS2 / 2;

  MAYBE_UNUSED const int i1_min = NGHOSTS;
  MAYBE_UNUSED const int i1_max = params->Nxx_plus_2NGHOSTS1 - NGHOSTS - 1;
  MAYBE_UNUSED const int i2_min = NGHOSTS;

  MAYBE_UNUSED const int i0_rmin = NGHOSTS; // rho = 0 (Cylindrical)
  MAYBE_UNUSED const int i1_pmin = NGHOSTS; // phi = -pi

  // Quarter-plane indices for cell-centered grids
  MAYBE_UNUSED const int i2_q1 = (int)(NGHOSTS + 0.25 * (REAL)N2int - 0.5);
  MAYBE_UNUSED const int i2_q3 = (int)(NGHOSTS + 0.75 * (REAL)N2int - 0.5);
  MAYBE_UNUSED const int i1_q1 = (int)(NGHOSTS + 0.25 * (REAL)N1int - 0.5);
  MAYBE_UNUSED const int i1_q3 = (int)(NGHOSTS + 0.75 * (REAL)N1int - 0.5);

  // File naming: out1d-AXIS-<CoordSystemName>-gridXX-...
  char coordsys_with_grid[128];
  snprintf(coordsys_with_grid, sizeof(coordsys_with_grid), "%s-grid%02d", params->CoordSystemName, grid);

  FILE *out_y = open_outfile("out1d-y", coordsys_with_grid, commondata, /*include_time=*/1);
  FILE *out_z = open_outfile("out1d-z", coordsys_with_grid, commondata, /*include_time=*/1);
  if (!out_y || !out_z) {{
    if (out_y)
      fclose(out_y);
    if (out_z)
      fclose(out_z);
    fprintf(stderr, "Error: Cannot open output files for grid %d.\n", grid);
    exit(1);
  }} // END IF file open failure

  // Emit time comment then axis-specific headers (no 'time' column)
  diag_write_time_comment(out_y, commondata->time);
  diag_write_time_comment(out_z, commondata->time);
  diag_write_header(out_y, "y", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
  diag_write_header(out_z, "z", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);

  // Source pointer for this grid
  const REAL *restrict src = gridfuncs_diags[grid];

  // Allocate buffers (tight upper bounds)
  const int max_y = {y_count_expr};
  const int max_z = {z_count_expr};
  data_point_1d_struct *data_points_y = max_y > 0 ? (data_point_1d_struct*)malloc(sizeof(data_point_1d_struct)*(size_t)max_y) : NULL;
  data_point_1d_struct *data_points_z = max_z > 0 ? (data_point_1d_struct*)malloc(sizeof(data_point_1d_struct)*(size_t)max_z) : NULL;

  // Row buffer: [axis_coord, gfs...]
  REAL *row = (REAL*)malloc(sizeof(REAL) * (size_t)(1 + NUM_GFS_NEAREST));
  if ((max_y > 0 && !data_points_y) || (max_z > 0 && !data_points_z) || !row) {{
    fprintf(stderr, "Error: Allocation failure in diagnostics_nearest_1d_y_and_z_axes.\\n");
    free(data_points_y); free(data_points_z); free(row);
    exit(1);
  }} // END IF allocation failure

  // ----------------------
  // Build y-axis samples
  // ----------------------
  int count_y = 0;
  {fill_y}
  if (count_y > 1) qsort(data_points_y, (size_t)count_y, sizeof(data_point_1d_struct), compare_by_coord);
  for (int p = 0; p < count_y; ++p) {{
    row[0] = data_points_y[p].coord;
    const int idx3 = data_points_y[p].idx3;
    for (int gf_i = 0; gf_i < NUM_GFS_NEAREST; ++gf_i) {{
      const int gf = which_gfs[gf_i];
      row[1 + gf_i] = src[IDX4Ppt(params, gf, idx3)];
    }} // END LOOP over gridfunctions
    diag_write_row(out_y, 1 + NUM_GFS_NEAREST, row);
  }} // END LOOP over *sorted* points closest to y-axis.

  // ----------------------
  // Build z-axis samples
  // ----------------------
  int count_z = 0;
  {fill_z}
  if (count_z > 1) qsort(data_points_z, (size_t)count_z, sizeof(data_point_1d_struct), compare_by_coord);
  for (int p = 0; p < count_z; ++p) {{
    row[0] = data_points_z[p].coord;
    const int idx3 = data_points_z[p].idx3;
    for (int gf_i = 0; gf_i < NUM_GFS_NEAREST; ++gf_i) {{
      const int gf = which_gfs[gf_i];
      row[1 + gf_i] = src[IDX4Ppt(params, gf, idx3)];
    }} // END LOOP over gridfunctions
    diag_write_row(out_z, 1 + NUM_GFS_NEAREST, row);
  }} // END LOOP over *sorted* points closest to z-axis.

  // Cleanup
  free(data_points_y);
  free(data_points_z);
  free(row);
  fclose(out_y);
  fclose(out_z);
"""

    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
        CoordSystem_for_wrapper_func=CoordSystem,
    )
    return pcg.NRPyEnv()
