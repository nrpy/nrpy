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
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Any, Dict, List, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.infrastructures.BHaH.diagnostics.diagnostics_nearest_1d_y_and_z_axes import (
    bhah_axis_configs,
    bhah_family_from_coord,
    count_expr,
)


def gen_fill(axis_char: str, lines: List[Dict[str, Any]]) -> str:
    """
    Return a C code snippet that fills and counts axis-line sample points.

    Given a list of axis-line specifications, emit C loops that populate the
    appropriate data_point_1d_struct array (y or z) with (axis coordinate, idx3, i0/i1/i2)
    entries and increment the corresponding counter.

    :param axis_char: Axis selector, "y" or "z".
    :param lines: Axis-line specifications (varying index and fixed indices).
    :return: C source code as a string.
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
    dim_index = {"i0": "0", "i1": "1", "i2": "2"}

    blocks: List[str] = []
    for line in lines:
        vary = line["vary"]  # e.g., "i0"
        vary_dim = dim_index[vary]  # e.g., "0"
        fixed = line["fixed"]  # e.g., {"i1":"i1_mid","i2":"i2_q1"}

        # Build the loop body
        body_lines: List[str] = [
            f"for (int {vary} = NGHOSTS; {vary} < params->Nxx_plus_2NGHOSTS{vary_dim} - NGHOSTS; {vary}++) {{",
        ]

        # Bind non-varying indices to named constants
        for i in ("i0", "i1", "i2"):
            if i != vary:
                body_lines.append(f"  const int {i} = {fixed[i]};")

        body_lines += [
            "  const int idx3 = IDX3P(params, i0, i1, i2);",
            "  REAL xCart[3], xOrig[3] = { xx[0][i0], xx[1][i1], xx[2][i2] };",
            "  xx_to_Cart(params, xOrig, xCart);",
            f"  {arr}[{cnt}].coord = xCart[{coord_idx}];",
            f"  {arr}[{cnt}].idx3  = idx3;",
            # store the indices too (your downstream code expects these!)
            f"  {arr}[{cnt}].i0    = i0;",
            f"  {arr}[{cnt}].i1    = i1;",
            f"  {arr}[{cnt}].i2    = i2;",
            f"  {cnt}++;",
            f"}} // END LOOP over {vary}",
        ]

        blocks.append("\n".join(body_lines))

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
  int i0;
  int i1;
  int i2;
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
 *   4) Chare (0,0,0) writes the single-line time comment and axis-specific header ("y" or "z");
 *    all chares stream their owned rows at precomputed file offsets.
 *
 * Gridfunction values are loaded from the flattened diagnostic array using IDX4Ppt with a 3D index
 * constructed by IDX3P. Sampling occurs at grid points only; no interpolation is performed.
 * On allocation or file-open failure, the routine prints an error message to stderr and terminates.
 * If a user-editable block is provided in the implementation, users may add custom logic such as
 * additional columns or filtering before rows are written.
 *
 * @param[in,out] commondata       Pointer to global simulation metadata.
 * @param[in]     grid             Grid index.
 * @param[in]     params           Pointer to per-grid parameters (global grid).
 * @param[in]     params_chare     Pointer to per-chare parameters (local grid).
 * @param[in]     xx               Global-grid logical coordinates.
 * @param[in]     xx_chare         Chare-local logical coordinates.
 * @param[in]     NUM_GFS_NEAREST  Number of diagnostic gridfunctions to output.
 * @param[in]     which_gfs        Array of length NUM_GFS_NEAREST giving gridfunction indices.
 * @param[in]     diagnostic_gf_names  Array of length NUM_GFS_NEAREST giving column names.
 * @param[in]     gridfuncs_diags  Per-grid pointers to diagnostic data arrays.
 * @param[in]     charecommstruct  Chare communication metadata/mappings.
 * @param[in,out] diagnosticstruct Diagnostic bookkeeping for point indices and file offsets.
 * @param[in]     chare_index      3D chare index.
 * @param[in]     token            Ck::IO session token.
 * @param[in]     which_diagnostics_part  Enum selecting diagnostics stage/action.
 *
 * @return     void
"""

    cfunc_type = "void"
    name = "diagnostics_nearest_1d_y_and_z_axes"
    params = """commondata_struct *restrict commondata, const int grid,
                const params_struct *restrict params, const params_struct *restrict params_chare,
                const REAL *restrict xx[3], const REAL *restrict xx_chare[3],
                const int NUM_GFS_NEAREST, const int which_gfs[],
                const char **diagnostic_gf_names, const REAL *restrict gridfuncs_diags[],
                const charecomm_struct *restrict charecommstruct,
                diagnostic_struct *restrict diagnosticstruct,
                const int chare_index[3], Ck::IO::Session token,
                const int which_diagnostics_part"""

    body = rf"""
#include "set_CodeParameters.h"
  switch (which_diagnostics_part) {{

    case DIAGNOSTICS_SETUP_1D: {{

      const int Nxx0chare = params_chare->Nxx0;
      const int Nxx1chare = params_chare->Nxx1;
      const int Nxx2chare = params_chare->Nxx2;

      // Build filename component with runtime coordinate system name and grid number
      char coordsys_with_grid[128];
      snprintf(coordsys_with_grid, sizeof(coordsys_with_grid), "%s-grid%02d", params->CoordSystemName, grid);
      strcpy(diagnosticstruct->filename_1d_y, coordsys_with_grid);
      strcpy(diagnosticstruct->filename_1d_z, coordsys_with_grid);

      diagnosticstruct->num_output_quantities = NUM_GFS_NEAREST;

      // Compute bytes common to both y and z outputs
      diagnosticstruct->sizeinbytes_per_pt_1d = 23 * (diagnosticstruct->num_output_quantities + 1);
      int time_bytes = diag_time_comment_size_bytes(commondata->time);

      // Interior counts
      const int N0int = params->Nxx_plus_2NGHOSTS0 - 2 * NGHOSTS;
      const int N1int = params->Nxx_plus_2NGHOSTS1 - 2 * NGHOSTS;
      const int N2int = params->Nxx_plus_2NGHOSTS2 - 2 * NGHOSTS;

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

      // Allocate buffers (tight upper bounds)
      const int max_y = {y_count_expr};
      const int max_z = {z_count_expr};
      data_point_1d_struct *data_points_y =
        max_y > 0 ? (data_point_1d_struct *)malloc(sizeof(data_point_1d_struct) * (size_t)max_y) : NULL;
      data_point_1d_struct *data_points_z =
        max_z > 0 ? (data_point_1d_struct *)malloc(sizeof(data_point_1d_struct) * (size_t)max_z) : NULL;

      // ----------------------
      // Build y-axis samples
      // ----------------------
      int count_y = 0;
      {fill_y}
      if (count_y > 1)
        qsort(data_points_y, (size_t)count_y, sizeof(data_point_1d_struct), compare_by_coord);

      // Set values for y
      diagnosticstruct->tot_num_diagnostic_1d_y_pts = count_y;
      int header_size_bytes_y =
        time_bytes + diag_header_size_bytes("y", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
      diagnosticstruct->totsizeinbytes_1d_y =
        header_size_bytes_y + (diagnosticstruct->sizeinbytes_per_pt_1d *
                               diagnosticstruct->tot_num_diagnostic_1d_y_pts);

      int num_diagnostics_chare = 0;
      for (int i = 0; i < count_y; i++) {{
        const int i0 = data_points_y[i].i0;
        const int i1 = data_points_y[i].i1;
        const int i2 = data_points_y[i].i2;
        const int idx3 = IDX3(i0, i1, i2);
        if (charecommstruct->globalidx3pt_to_chareidx3[idx3] ==
            IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {{
          num_diagnostics_chare++;
        }}
      }}

      diagnosticstruct->num_diagnostic_1d_y_pts = num_diagnostics_chare;

      diagnosticstruct->localidx3_diagnostic_1d_y_pt =
        (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali0_diagnostic_1d_y_pt =
        (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali1_diagnostic_1d_y_pt =
        (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali2_diagnostic_1d_y_pt =
        (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->offset_diagnostic_1d_y_pt =
        (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);

      int which_diagnostics_chare = 0;
      int which_diagnostic_global = 0;
      for (int i = 0; i < count_y; i++) {{
        const int i0 = data_points_y[i].i0;
        const int i1 = data_points_y[i].i1;
        const int i2 = data_points_y[i].i2;
        const int idx3 = IDX3(i0, i1, i2);
        if (charecommstruct->globalidx3pt_to_chareidx3[idx3] ==
            IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {{
          int localidx3 = charecommstruct->globalidx3pt_to_localidx3pt[idx3];
          diagnosticstruct->localidx3_diagnostic_1d_y_pt[which_diagnostics_chare] = localidx3;
          diagnosticstruct->locali0_diagnostic_1d_y_pt[which_diagnostics_chare] =
            MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0, Nxx0chare);
          diagnosticstruct->locali1_diagnostic_1d_y_pt[which_diagnostics_chare] =
            MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1, Nxx1chare);
          diagnosticstruct->locali2_diagnostic_1d_y_pt[which_diagnostics_chare] =
            MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2, Nxx2chare);
          diagnosticstruct->offset_diagnostic_1d_y_pt[which_diagnostics_chare] =
            header_size_bytes_y + (which_diagnostic_global *
                                   diagnosticstruct->sizeinbytes_per_pt_1d);
          which_diagnostics_chare++;
        }}
        which_diagnostic_global++;
      }}

      // ----------------------
      // Build z-axis samples
      // ----------------------
      int count_z = 0;
      {fill_z}
      if (count_z > 1)
        qsort(data_points_z, (size_t)count_z, sizeof(data_point_1d_struct), compare_by_coord);

      diagnosticstruct->tot_num_diagnostic_1d_z_pts = count_z;
      int header_size_bytes_z =
        time_bytes + diag_header_size_bytes("z", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
      diagnosticstruct->totsizeinbytes_1d_z =
        header_size_bytes_z + (diagnosticstruct->sizeinbytes_per_pt_1d *
                               diagnosticstruct->tot_num_diagnostic_1d_z_pts);

      num_diagnostics_chare = 0;
      for (int i = 0; i < count_z; i++) {{
        const int i0 = data_points_z[i].i0;
        const int i1 = data_points_z[i].i1;
        const int i2 = data_points_z[i].i2;
        const int idx3 = IDX3(i0, i1, i2);
        if (charecommstruct->globalidx3pt_to_chareidx3[idx3] ==
            IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {{
          num_diagnostics_chare++;
        }}
      }}

      diagnosticstruct->num_diagnostic_1d_z_pts = num_diagnostics_chare;

      diagnosticstruct->localidx3_diagnostic_1d_z_pt =
        (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali0_diagnostic_1d_z_pt =
        (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali1_diagnostic_1d_z_pt =
        (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->locali2_diagnostic_1d_z_pt =
        (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
      diagnosticstruct->offset_diagnostic_1d_z_pt =
        (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);

      which_diagnostics_chare = 0;
      which_diagnostic_global = 0;
      for (int i = 0; i < count_z; i++) {{
        const int i0 = data_points_z[i].i0;
        const int i1 = data_points_z[i].i1;
        const int i2 = data_points_z[i].i2;
        const int idx3 = IDX3(i0, i1, i2);
        if (charecommstruct->globalidx3pt_to_chareidx3[idx3] ==
            IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {{
          int localidx3 = charecommstruct->globalidx3pt_to_localidx3pt[idx3];
          diagnosticstruct->localidx3_diagnostic_1d_z_pt[which_diagnostics_chare] = localidx3;
          diagnosticstruct->locali0_diagnostic_1d_z_pt[which_diagnostics_chare] =
            MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0, Nxx0chare);
          diagnosticstruct->locali1_diagnostic_1d_z_pt[which_diagnostics_chare] =
            MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1, Nxx1chare);
          diagnosticstruct->locali2_diagnostic_1d_z_pt[which_diagnostics_chare] =
            MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2, Nxx2chare);
          diagnosticstruct->offset_diagnostic_1d_z_pt[which_diagnostics_chare] =
            header_size_bytes_z + (which_diagnostic_global *
                                   diagnosticstruct->sizeinbytes_per_pt_1d);
          which_diagnostics_chare++;
        }}
        which_diagnostic_global++;
      }}

      // Cleanup temporary point buffers
      free(data_points_y);
      free(data_points_z);

      break;
    }}

    case DIAGNOSTICS_WRITE_Y: {{

      // only chare (0,0,0) writes header
      if (chare_index[0] == 0 && chare_index[1] == 0 && chare_index[2] == 0) {{
        int header_bytes = diag_time_comment_size_bytes(commondata->time)
                         + diag_header_size_bytes("y", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);

        char *hdr = (char *)malloc((size_t)header_bytes + 1);
        int written = diag_ckio_build_time_comment_and_header(
            hdr, (size_t)header_bytes + 1,
            commondata->time, "y",
            NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
        if (written <= 0 || written > header_bytes) {{
           fprintf(stderr, "Error: failed to build diagnostics header.\n");
           free(hdr);
           break;
         }}
        Ck::IO::write(token, hdr, header_bytes, 0);
        free(hdr);
      }}

      // Source pointer for this grid
      const REAL *restrict src = gridfuncs_diags[grid];

      // Row buffer: [axis_coord, gfs...]
      const int NUM_COLS = 1 + NUM_GFS_NEAREST;
      REAL *row = (REAL *)malloc(sizeof(REAL) * (size_t)NUM_COLS);
      if (!row) {{
        fprintf(stderr, "Error: Failed to allocate memory for row buffer.\\n");
        exit(1);
      }}

      // Unpack diagnosticstruct
      const int num_diagnostic_pts = diagnosticstruct->num_diagnostic_1d_y_pts;
      const int *restrict idx3_diagnostic_pt = diagnosticstruct->localidx3_diagnostic_1d_y_pt;
      const int *restrict i0_diagnostic_pt   = diagnosticstruct->locali0_diagnostic_1d_y_pt;
      const int *restrict i1_diagnostic_pt   = diagnosticstruct->locali1_diagnostic_1d_y_pt;
      const int *restrict i2_diagnostic_pt   = diagnosticstruct->locali2_diagnostic_1d_y_pt;
      const int *restrict offsetpt_firstfield = diagnosticstruct->offset_diagnostic_1d_y_pt;

      for (int which_pt = 0; which_pt < num_diagnostic_pts; which_pt++) {{
        const int idx3 = idx3_diagnostic_pt[which_pt];
        const int i0   = i0_diagnostic_pt[which_pt];
        const int i1   = i1_diagnostic_pt[which_pt];
        const int i2   = i2_diagnostic_pt[which_pt];

        REAL xCart[3];
        REAL xOrig[3] = {{ xx_chare[0][i0], xx_chare[1][i1], xx_chare[2][i2] }};
        xx_to_Cart(params_chare, xOrig, xCart);

        int sizeinbytes = diagnosticstruct->sizeinbytes_per_pt_1d;
        char out[sizeinbytes + 1];

        row[0] = xCart[1];
        for (int gf_i = 0; gf_i < NUM_GFS_NEAREST; ++gf_i) {{
          const int gf = which_gfs[gf_i];
          row[1 + gf_i] = src[IDX4Ppt(params_chare, gf, idx3)];
        }}

        int n = 0;
        n += snprintf(out + n, (size_t)(sizeinbytes + 1 - n), "% .15e", row[0]);
        for (int col = 1; col < NUM_COLS; col++) {{
          n += snprintf(out + n, (size_t)(sizeinbytes + 1 - n), " % .15e", row[col]);
        }}
        out[sizeinbytes - 1] = '\n';
        Ck::IO::write(token, out, sizeinbytes, offsetpt_firstfield[which_pt]);
      }}

      free(row);
      break;
    }}

    case DIAGNOSTICS_WRITE_Z: {{

      // only chare (0,0,0) writes header
      if (chare_index[0] == 0 && chare_index[1] == 0 && chare_index[2] == 0) {{
        int header_bytes = diag_time_comment_size_bytes(commondata->time)
                         + diag_header_size_bytes("z", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);

        char *hdr = (char *)malloc((size_t)header_bytes + 1);
        int written = diag_ckio_build_time_comment_and_header(
            hdr, (size_t)header_bytes + 1,
            commondata->time, "z",
            NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
        if (written <= 0 || written > header_bytes) {{
           fprintf(stderr, "Error: failed to build diagnostics header.\n");
           free(hdr);
           break;
         }}
        Ck::IO::write(token, hdr, header_bytes, 0);
        free(hdr);
      }}

      // Source pointer for this grid
      const REAL *restrict src = gridfuncs_diags[grid];

      const int NUM_COLS = 1 + NUM_GFS_NEAREST;
      REAL *row = (REAL *)malloc(sizeof(REAL) * (size_t)NUM_COLS);
      if (!row) {{
        fprintf(stderr, "Error: Failed to allocate memory for row buffer.\\n");
        exit(1);
      }}

      // Unpack diagnosticstruct
      const int num_diagnostic_pts = diagnosticstruct->num_diagnostic_1d_z_pts;
      const int *restrict idx3_diagnostic_pt = diagnosticstruct->localidx3_diagnostic_1d_z_pt;
      const int *restrict i0_diagnostic_pt   = diagnosticstruct->locali0_diagnostic_1d_z_pt;
      const int *restrict i1_diagnostic_pt   = diagnosticstruct->locali1_diagnostic_1d_z_pt;
      const int *restrict i2_diagnostic_pt   = diagnosticstruct->locali2_diagnostic_1d_z_pt;
      const int *restrict offsetpt_firstfield = diagnosticstruct->offset_diagnostic_1d_z_pt;

      for (int which_pt = 0; which_pt < num_diagnostic_pts; which_pt++) {{
        const int idx3 = idx3_diagnostic_pt[which_pt];
        const int i0   = i0_diagnostic_pt[which_pt];
        const int i1   = i1_diagnostic_pt[which_pt];
        const int i2   = i2_diagnostic_pt[which_pt];

        REAL xCart[3];
        REAL xOrig[3] = {{ xx_chare[0][i0], xx_chare[1][i1], xx_chare[2][i2] }};
        xx_to_Cart(params_chare, xOrig, xCart);

        int sizeinbytes = diagnosticstruct->sizeinbytes_per_pt_1d;
        char out[sizeinbytes + 1];

        row[0] = xCart[2];
        for (int gf_i = 0; gf_i < NUM_GFS_NEAREST; ++gf_i) {{
          const int gf = which_gfs[gf_i];
          row[1 + gf_i] = src[IDX4Ppt(params_chare, gf, idx3)];
        }}

        int n = 0;
        n += snprintf(out + n, (size_t)(sizeinbytes + 1 - n), "% .15e", row[0]);
        for (int col = 1; col < NUM_COLS; col++) {{
          n += snprintf(out + n, (size_t)(sizeinbytes + 1 - n), " % .15e", row[col]);
        }}
        out[sizeinbytes - 1] = '\n';
        Ck::IO::write(token, out, sizeinbytes, offsetpt_firstfield[which_pt]);
      }}

      free(row);
      break;
    }}
  }}
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
