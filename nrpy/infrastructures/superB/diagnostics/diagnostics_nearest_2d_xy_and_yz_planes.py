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
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Any, Dict, Sequence, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.infrastructures.BHaH.diagnostics.diagnostics_nearest_2d_xy_and_yz_planes import (
    bhah_plane_configs,    
    get_coord_family,
)
def _indent(code: str, n: int = 8) -> str:
    pad = " " * n
    return "\n".join(pad + line if line.strip() else line for line in code.splitlines())


def _emit_setup_plane_code(cfg: Dict[str, Any], plane: str) -> tuple[str, str, str]:
    """
    Return (tot_pts_expr, count_loops_code, fill_loops_code) for DIAGNOSTICS_SETUP_2D,
    using i0/i1/i2 variable names so IDX3P(...) works unchanged.
    """
    fixed_dim: str = cfg["fixed_dim"]
    fixed_val = cfg["fixed_val"]  # str or list[str]
    loop_dims: Sequence[str] = cfg["loop_dims"]  # [outer, inner]

    Nint = {"i0": "N0int", "i1": "N1int", "i2": "N2int"}
    i_end = {"i0": "i0_end", "i1": "i1_end", "i2": "i2_end"}

    # total points expression
    nslices = len(fixed_val) if isinstance(fixed_val, list) else 1
    tot_expr = f"{nslices} * {Nint[loop_dims[0]]} * {Nint[loop_dims[1]]}"

    # loop nest builder (outer = loop_dims[0], inner = loop_dims[1])
    def nest(inner: str) -> str:
        code = inner
        # wrap inner first, then outer
        for dim in reversed(loop_dims):
            code = (
                f"for (int {dim} = NGHOSTS; {dim} < {i_end[dim]}; {dim}++) {{\n"
                f"{_indent(code, 2)}\n"
                f"}} // END LOOP over {dim}\n"
            )
        return code

    owned_check = (
        "if (charecommstruct->globalidx3pt_to_chareidx3[idx3] ==\n"
        "    IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2]))"
    )

    # counting inner body
    count_inner = (
        "const int idx3 = IDX3P(params, i0, i1, i2);\n"
        f"{owned_check} {{\n"
        "  num_diagnostics_chare++;\n"
        "}\n"
    )

    # fill inner body (plane-specific field names)
    suf = "xy" if plane == "xy" else "yz"
    fill_inner = (
        "const int idx3 = IDX3P(params, i0, i1, i2);\n"
        f"{owned_check} {{\n"
        "  const int localidx3 = charecommstruct->globalidx3pt_to_localidx3pt[idx3];\n"
        f"  diagnosticstruct->localidx3_diagnostic_2d_{suf}_pt[which_diagnostics_chare] = localidx3;\n"
        f"  diagnosticstruct->locali0_diagnostic_2d_{suf}_pt[which_diagnostics_chare] =\n"
        "      MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0, Nxx0chare);\n"
        f"  diagnosticstruct->locali1_diagnostic_2d_{suf}_pt[which_diagnostics_chare] =\n"
        "      MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1, Nxx1chare);\n"
        f"  diagnosticstruct->locali2_diagnostic_2d_{suf}_pt[which_diagnostics_chare] =\n"
        "      MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2, Nxx2chare);\n"
        f"  diagnosticstruct->offset_diagnostic_2d_{suf}_pt[which_diagnostics_chare] =\n"
        f"      header_size_bytes_{suf} + (which_diagnostic_global * diagnosticstruct->sizeinbytes_per_pt_2d);\n"
        "  which_diagnostics_chare++;\n"
        "}\n"
        "which_diagnostic_global++;\n"
    )

    count_loops = nest(count_inner)
    fill_loops = nest(fill_inner)

    # wrap fixed-dim selection (single slice vs multiple slices)
    if isinstance(fixed_val, list):
        slices_decl = f"const int {fixed_dim}_slices[{nslices}] = {{{', '.join(fixed_val)}}};"
        count_loops = (
            f"{slices_decl}\n"
            f"for (int slice = 0; slice < {nslices}; slice++) {{\n"
            f"  const int {fixed_dim} = {fixed_dim}_slices[slice];\n"
            f"{_indent(count_loops, 2)}\n"
            f"}} // END LOOP over slices\n"
        )
        fill_loops = (
            f"{slices_decl}\n"
            f"for (int slice = 0; slice < {nslices}; slice++) {{\n"
            f"  const int {fixed_dim} = {fixed_dim}_slices[slice];\n"
            f"{_indent(fill_loops, 2)}\n"
            f"}} // END LOOP over slices\n"
        )
    else:
        count_loops = f"const int {fixed_dim} = {fixed_val};\n{count_loops}"
        fill_loops = f"const int {fixed_dim} = {fixed_val};\n{fill_loops}"

    return tot_expr, count_loops, fill_loops


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
    # ~ plane_configs = bhah_plane_configs()
    # ~ xy_plane_code = generate_plane_loop_code(plane_configs["xy"][family], "xy")
    # ~ yz_plane_code = generate_plane_loop_code(plane_configs["yz"][family], "yz")
    
    # ~ print(xy_plane_code)
    
    plane_configs = bhah_plane_configs()
    xy_cfg = plane_configs["xy"][family]
    yz_cfg = plane_configs["yz"][family]

    xy_tot_expr, xy_count_loops, xy_fill_loops = _emit_setup_plane_code(xy_cfg, "xy")
    yz_tot_expr, yz_count_loops, yz_fill_loops = _emit_setup_plane_code(yz_cfg, "yz")



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
 *  - Writes to two persistent per-grid files, one for the xy plane and one for the yz plane.
 *  - Chare (0,0,0) writes the time comment and header; all chares write their owned rows at precomputed offsets.
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
 * @param[in,out] commondata            Pointer to global simulation metadata.
 * @param[in]     grid                  Grid index.
 * @param[in]     params                Pointer to per-grid parameters (global grid).
 * @param[in]     params_chare          Pointer to per-chare parameters (local grid).
 * @param[in]     xx                    Global-grid logical coordinates.
 * @param[in]     xx_chare              Chare-local logical coordinates.
 * @param[in]     NUM_GFS_NEAREST       Number of diagnostic gridfunctions to output per point.
 * @param[in]     which_gfs             Array of length NUM_GFS_NEAREST giving gridfunction indices.
 * @param[in]     diagnostic_gf_names   Array of length NUM_GFS_NEAREST giving column names.
 * @param[in]     gridfuncs_diags       Per-grid pointers to diagnostic data arrays.
 * @param[in]     charecommstruct       Chare communication metadata/mappings.
 * @param[in,out] diagnosticstruct      Diagnostic bookkeeping for point indices and file offsets.
 * @param[in]     chare_index           3D chare index.
 * @param[in]     token                Ck::IO session token.
 * @param[in]     which_diagnostics_part Enum selecting diagnostics stage/action.
 *
 * @return        void                  No return value. On success, output is written via the provided Ck::IO session token.
"""
    cfunc_type = "void"
    name = "diagnostics_nearest_2d_xy_and_yz_planes"
    params = """commondata_struct *restrict commondata, const int grid,
            const params_struct *restrict params, const params_struct *restrict params_chare,
            const REAL *restrict xx[3], const REAL *restrict xx_chare[3],
            const int NUM_GFS_NEAREST, const int which_gfs[], const char **diagnostic_gf_names,
            const REAL *restrict gridfuncs_diags[],
            const charecomm_struct *restrict charecommstruct, diagnostic_struct *restrict diagnosticstruct,
            const int chare_index[3], Ck::IO::Session token, const int which_diagnostics_part"""

    body = rf"""
#include "set_CodeParameters.h"


  switch (which_diagnostics_part) {{

    case DIAGNOSTICS_SETUP_2D: {{

      const int Nxx0chare = params_chare->Nxx0;
      const int Nxx1chare = params_chare->Nxx1;
      const int Nxx2chare = params_chare->Nxx2;

      // Build filename component with runtime coordinate system name and grid number
      char coordsys_with_grid[128];
      snprintf(coordsys_with_grid, sizeof(coordsys_with_grid), "%s-grid%02d", params->CoordSystemName, grid);
      strcpy(diagnosticstruct->filename_2d_xy, coordsys_with_grid);
      strcpy(diagnosticstruct->filename_2d_yz, coordsys_with_grid);

      diagnosticstruct->num_output_quantities = NUM_GFS_NEAREST;

      // Compute bytes common to both xy and yz outputs
      diagnosticstruct->sizeinbytes_per_pt_2d = 23 * (diagnosticstruct->num_output_quantities + 2);
      int time_bytes = diag_time_comment_size_bytes(commondata->time);

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

      // --- Sample and setup data for the xy-plane ---
      {{
        // Set values
        diagnosticstruct->tot_num_diagnostic_2d_xy_pts = {xy_tot_expr};
        int header_size_bytes_xy =
            time_bytes + diag_header_size_bytes("x y", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
        diagnosticstruct->totsizeinbytes_2d_xy =
            header_size_bytes_xy +
            (diagnosticstruct->sizeinbytes_per_pt_2d * diagnosticstruct->tot_num_diagnostic_2d_xy_pts);

        int num_diagnostics_chare = 0;
        {{
{xy_count_loops}
        }}

        diagnosticstruct->num_diagnostic_2d_xy_pts = num_diagnostics_chare;
        diagnosticstruct->locali0_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->locali1_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->locali2_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->localidx3_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->offset_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);

        int which_diagnostics_chare = 0;
        int which_diagnostic_global = 0;
        {{
{xy_fill_loops}
        }}
      }}

      // --- Sample and setup data for the yz-plane ---
      {{
        // Set values
        diagnosticstruct->tot_num_diagnostic_2d_yz_pts = {yz_tot_expr};
        int header_size_bytes_yz =
            time_bytes + diag_header_size_bytes("y z", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
        diagnosticstruct->totsizeinbytes_2d_yz =
            header_size_bytes_yz +
            (diagnosticstruct->sizeinbytes_per_pt_2d * diagnosticstruct->tot_num_diagnostic_2d_yz_pts);

        int num_diagnostics_chare = 0;
        {{
{yz_count_loops}
        }}

        diagnosticstruct->num_diagnostic_2d_yz_pts = num_diagnostics_chare;
        diagnosticstruct->locali0_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->locali1_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->locali2_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->localidx3_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
        diagnosticstruct->offset_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);

        int which_diagnostics_chare = 0;
        int which_diagnostic_global = 0;
        {{
{yz_fill_loops}
        }}
      }}

      break;
    }}

    case DIAGNOSTICS_WRITE_XY: {{

      // Only chare (0,0,0) writes header
      if (chare_index[0] == 0 && chare_index[1] == 0 && chare_index[2] == 0) {{
        int header_bytes = diag_time_comment_size_bytes(commondata->time) +
                           diag_header_size_bytes("x y", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);

        char *hdr = (char *)malloc((size_t)header_bytes + 1);
        int written = diag_ckio_build_time_comment_and_header(
            hdr, (size_t)header_bytes + 1,
            commondata->time, "x y",
            NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
        if (written <= 0 || written > header_bytes) {{
           fprintf(stderr, "Error: failed to build diagnostics header.\n");
           free(hdr);
           break;
        }}
        Ck::IO::write(token, hdr, header_bytes, 0);
        free(hdr);
      }}

      // Active grid data pointer and reusable row buffer
      const REAL *restrict src = gridfuncs_diags[grid];
      const int NUM_COLS = 2 + NUM_GFS_NEAREST;
      REAL *row = (REAL *)malloc(sizeof(REAL) * (size_t)NUM_COLS);
      if (!row) {{
        fprintf(stderr, "Error: Failed to allocate memory for row buffer.\n");
        exit(1);
      }} // END IF row allocation failure

      // Unpack diagnosticstruct
      const int num_diagnostic_pts = diagnosticstruct->num_diagnostic_2d_xy_pts;
      const int *restrict idx3_diagnostic_pt = diagnosticstruct->localidx3_diagnostic_2d_xy_pt;
      const int *restrict i0_diagnostic_pt = diagnosticstruct->locali0_diagnostic_2d_xy_pt;
      const int *restrict i1_diagnostic_pt = diagnosticstruct->locali1_diagnostic_2d_xy_pt;
      const int *restrict i2_diagnostic_pt = diagnosticstruct->locali2_diagnostic_2d_xy_pt;
      const int *restrict offsetpt_firstfield = diagnosticstruct->offset_diagnostic_2d_xy_pt;

      for (int which_pt = 0; which_pt < num_diagnostic_pts; which_pt++) {{
        const int idx3 = idx3_diagnostic_pt[which_pt];
        const int i0 = i0_diagnostic_pt[which_pt];
        const int i1 = i1_diagnostic_pt[which_pt];
        const int i2 = i2_diagnostic_pt[which_pt];
        REAL xCart[3];

        REAL xOrig[3] = {{xx_chare[0][i0], xx_chare[1][i1], xx_chare[2][i2]}};
        xx_to_Cart(params_chare, xOrig, xCart);

        int sizeinbytes = diagnosticstruct->sizeinbytes_per_pt_2d;
        char out[sizeinbytes + 1];
        row[0] = xCart[0];
        row[1] = xCart[1];
        for (int gf_idx = 0; gf_idx < NUM_GFS_NEAREST; gf_idx++) {{
          const int gf = which_gfs[gf_idx];
          row[2 + gf_idx] = src[IDX4Ppt(params_chare, gf, idx3)];
        }} // END LOOP over gridfunctions

        int n = 0;
        n += sprintf(out + n, "% .15e", row[0]);
        for (int col = 1; col < NUM_COLS; col++) {{
          n += sprintf(out + n, " % .15e", row[col]);
        }}
        out[sizeinbytes - 1] = '\n';
        Ck::IO::write(token, out, sizeinbytes, offsetpt_firstfield[which_pt]);
      }}
      // Finalize
      free(row);
      break;
    }}

    case DIAGNOSTICS_WRITE_YZ: {{

      // only chare (0,0,0) writes header
      if (chare_index[0] == 0 && chare_index[1] == 0 && chare_index[2] == 0) {{
        int header_bytes = diag_time_comment_size_bytes(commondata->time) +
                           diag_header_size_bytes("y z", NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);

        char *hdr = (char *)malloc((size_t)header_bytes + 1);
        int written = diag_ckio_build_time_comment_and_header(
            hdr, (size_t)header_bytes + 1,
            commondata->time, "y z",
            NUM_GFS_NEAREST, which_gfs, diagnostic_gf_names);
        if (written <= 0 || written > header_bytes) {{
           fprintf(stderr, "Error: failed to build diagnostics header.\n");
           free(hdr);
           break;
        }}
        Ck::IO::write(token, hdr, header_bytes, 0);
        free(hdr);
      }}

      // Active grid data pointer and reusable row buffer
      const REAL *restrict src = gridfuncs_diags[grid];
      const int NUM_COLS = 2 + NUM_GFS_NEAREST;
      REAL *row = (REAL *)malloc(sizeof(REAL) * (size_t)NUM_COLS);
      if (!row) {{
        fprintf(stderr, "Error: Failed to allocate memory for row buffer.\n");
        exit(1);
      }} // END IF row allocation failure

      // Unpack diagnosticstruct
      const int num_diagnostic_pts = diagnosticstruct->num_diagnostic_2d_yz_pts;
      const int *restrict idx3_diagnostic_pt = diagnosticstruct->localidx3_diagnostic_2d_yz_pt;
      const int *restrict i0_diagnostic_pt = diagnosticstruct->locali0_diagnostic_2d_yz_pt;
      const int *restrict i1_diagnostic_pt = diagnosticstruct->locali1_diagnostic_2d_yz_pt;
      const int *restrict i2_diagnostic_pt = diagnosticstruct->locali2_diagnostic_2d_yz_pt;
      const int *restrict offsetpt_firstfield = diagnosticstruct->offset_diagnostic_2d_yz_pt;

      for (int which_pt = 0; which_pt < num_diagnostic_pts; which_pt++) {{
        const int idx3 = idx3_diagnostic_pt[which_pt];
        const int i0 = i0_diagnostic_pt[which_pt];
        const int i1 = i1_diagnostic_pt[which_pt];
        const int i2 = i2_diagnostic_pt[which_pt];
        REAL xCart[3];
        REAL xOrig[3] = {{xx_chare[0][i0], xx_chare[1][i1], xx_chare[2][i2]}};
        xx_to_Cart(params_chare, xOrig, xCart);

        int sizeinbytes = diagnosticstruct->sizeinbytes_per_pt_2d;
        char out[sizeinbytes + 1];
        row[0] = xCart[1];
        row[1] = xCart[2];
        for (int gf_idx = 0; gf_idx < NUM_GFS_NEAREST; gf_idx++) {{
          const int gf = which_gfs[gf_idx];
          row[2 + gf_idx] = src[IDX4Ppt(params_chare, gf, idx3)];
        }} // END LOOP over gridfunctions

        int n = 0;
        n += sprintf(out + n, "% .15e", row[0]);
        for (int col = 1; col < NUM_COLS; col++) {{
          n += sprintf(out + n, " % .15e", row[col]);
        }}
        out[sizeinbytes - 1] = '\n';
        Ck::IO::write(token, out, sizeinbytes, offsetpt_firstfield[which_pt]);
      }}
      // Finalize
      free(row);
      break;
    }}
  }}
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
