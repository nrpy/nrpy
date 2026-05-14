"""
GRoovy-specific nearest-point diagnostics registration.

This module registers a diagnostics_nearest() dispatcher for GRHD runs. It
preserves the existing constraint-oriented nearest diagnostics and adds a
separate 0D hydro-center output that samples GRoovy evolution and auxiliary
evolution gridfunctions directly.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Set, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_diagnostics_nearest(
    set_of_CoordSystems: Set[str],
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
    evolving_neutrinos: bool = False,
    include_constraint_diagnostics: bool = True,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the GRoovy nearest-sampled diagnostics dispatcher.

    The generated C function has the same signature as the generic BHaH
    diagnostics_nearest() dispatcher, so it can be used by the top-level
    diagnostics() driver without changing that driver. In addition to the
    existing 0D/1D/2D nearest constraint outputs, it writes an "out0d-hydro"
    file for each grid by directly sampling GRHD primitive and conserved
    fields at the grid point nearest the physical center.

    :param set_of_CoordSystems: Coordinate systems that may appear in the generated run.
    :param evolving_temperature: If True, include Ye, temperature, and Ye_star columns.
    :param evolving_entropy: If True, include S and S_star columns.
    :param evolving_neutrinos: If True, include NRPyLeakage optical-depth, opacity, and
                               luminosity columns.
    :param include_constraint_diagnostics: If True, preserve the generic nearest
                                           constraint outputs.
    :return: None if in the parallel-codegen registration phase, otherwise the updated
             NRPy environment.
    :raises ValueError: If an unsupported coordinate system is requested.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    hydro_center_columns: List[Tuple[str, str]] = [
        ("rhob", "auxevol_gfs[IDX4Ppt(params, RHOBGF, idx3)]"),
        ("P", "auxevol_gfs[IDX4Ppt(params, PGF, idx3)]"),
        ("u4Ut", "auxevol_gfs[IDX4Ppt(params, U4UTGF, idx3)]"),
        ("rescaledvU0", "auxevol_gfs[IDX4Ppt(params, RESCALEDVU0GF, idx3)]"),
        ("rescaledvU1", "auxevol_gfs[IDX4Ppt(params, RESCALEDVU1GF, idx3)]"),
        ("rescaledvU2", "auxevol_gfs[IDX4Ppt(params, RESCALEDVU2GF, idx3)]"),
        ("rho_star", "y_n_gfs[IDX4Ppt(params, RHO_STARGF, idx3)]"),
        ("tau_tilde", "y_n_gfs[IDX4Ppt(params, TAU_TILDEGF, idx3)]"),
        (
            "rescaledstildeD0",
            "y_n_gfs[IDX4Ppt(params, RESCALEDSTILDED0GF, idx3)]",
        ),
        (
            "rescaledstildeD1",
            "y_n_gfs[IDX4Ppt(params, RESCALEDSTILDED1GF, idx3)]",
        ),
        (
            "rescaledstildeD2",
            "y_n_gfs[IDX4Ppt(params, RESCALEDSTILDED2GF, idx3)]",
        ),
    ]
    if evolving_temperature:
        hydro_center_columns += [
            ("Ye", "auxevol_gfs[IDX4Ppt(params, YEGF, idx3)]"),
            ("temperature", "auxevol_gfs[IDX4Ppt(params, TEMPERATUREGF, idx3)]"),
            ("Ye_star", "y_n_gfs[IDX4Ppt(params, YE_STARGF, idx3)]"),
        ]
    if evolving_entropy:
        hydro_center_columns += [
            ("S", "auxevol_gfs[IDX4Ppt(params, SGF, idx3)]"),
            ("S_star", "y_n_gfs[IDX4Ppt(params, S_STARGF, idx3)]"),
        ]
    if evolving_neutrinos:
        hydro_center_columns += [
            ("tau_0_nue", "auxevol_gfs[IDX4Ppt(params, TAU_0_NUEGF, idx3)]"),
            ("tau_1_nue", "auxevol_gfs[IDX4Ppt(params, TAU_1_NUEGF, idx3)]"),
            ("tau_0_anue", "auxevol_gfs[IDX4Ppt(params, TAU_0_ANUEGF, idx3)]"),
            ("tau_1_anue", "auxevol_gfs[IDX4Ppt(params, TAU_1_ANUEGF, idx3)]"),
            ("tau_0_nux", "auxevol_gfs[IDX4Ppt(params, TAU_0_NUXGF, idx3)]"),
            ("tau_1_nux", "auxevol_gfs[IDX4Ppt(params, TAU_1_NUXGF, idx3)]"),
            ("kappa_0_nue", "auxevol_gfs[IDX4Ppt(params, KAPPA_0_NUEGF, idx3)]"),
            ("kappa_1_nue", "auxevol_gfs[IDX4Ppt(params, KAPPA_1_NUEGF, idx3)]"),
            ("kappa_0_anue", "auxevol_gfs[IDX4Ppt(params, KAPPA_0_ANUEGF, idx3)]"),
            ("kappa_1_anue", "auxevol_gfs[IDX4Ppt(params, KAPPA_1_ANUEGF, idx3)]"),
            ("kappa_0_nux", "auxevol_gfs[IDX4Ppt(params, KAPPA_0_NUXGF, idx3)]"),
            ("kappa_1_nux", "auxevol_gfs[IDX4Ppt(params, KAPPA_1_NUXGF, idx3)]"),
            ("lum_nue", "auxevol_gfs[IDX4Ppt(params, LUM_NUEGF, idx3)]"),
            ("lum_anue", "auxevol_gfs[IDX4Ppt(params, LUM_ANUEGF, idx3)]"),
            ("lum_nux", "auxevol_gfs[IDX4Ppt(params, LUM_NUXGF, idx3)]"),
        ]

    center_index_cases = ""
    for CoordSystem in sorted(set_of_CoordSystems):
        if (
            ("Spherical" in CoordSystem)
            or ("Cylindrical" in CoordSystem)
            or ("SymTP" in CoordSystem)
        ):
            i0_center_expr = "NGHOSTS"
            i1_center_expr = "params->Nxx_plus_2NGHOSTS1 / 2"
            i2_center_expr = "params->Nxx_plus_2NGHOSTS2 / 2"
        elif "Cartesian" in CoordSystem or CoordSystem.startswith("GeneralRFM_fisheye"):
            i0_center_expr = "params->Nxx_plus_2NGHOSTS0 / 2"
            i1_center_expr = "params->Nxx_plus_2NGHOSTS1 / 2"
            i2_center_expr = "params->Nxx_plus_2NGHOSTS2 / 2"
        else:
            raise ValueError(f"Unsupported CoordSystem: {CoordSystem}")
        center_index_cases += f"""
    case {CoordSystem.upper()}:
      i0_center = {i0_center_expr};
      i1_center = {i1_center_expr};
      i2_center = {i2_center_expr};
      break;"""

    hydro_header = "# time " + " ".join(
        f"{name}({column})"
        for column, (name, _) in enumerate(hydro_center_columns, start=2)
    )
    num_hydro_center_cols = 1 + len(hydro_center_columns)
    hydro_row_assignments = "\n".join(
        f"      hydro_center_row[{column}] = {expr};"
        for column, (_, expr) in enumerate(hydro_center_columns, start=1)
    )

    constraint_selection = ""
    constraint_calls = ""
    if include_constraint_diagnostics:
        constraint_selection = r"""
  // 0D diagnostics: nearest point to the grid center.
  const int which_gfs_0d[] = {DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF};

  // 1D diagnostics: nearest lines to the y and z axes.
  const int which_gfs_1d[] = {DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF};

  // 2D diagnostics: nearest planes to the xy and yz coordinate planes.
  const int which_gfs_2d[] = {DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF};
"""
        constraint_calls = r"""
    const int NUM_nearest_GFS_0d = (int)(sizeof which_gfs_0d / sizeof which_gfs_0d[0]);
    diagnostics_nearest_grid_center(commondata, grid, params, xx, NUM_nearest_GFS_0d, which_gfs_0d,
                                    diagnostic_gf_names, gridfuncs_diags);

    const int NUM_nearest_GFS_1d = (int)(sizeof which_gfs_1d / sizeof which_gfs_1d[0]);
    diagnostics_nearest_1d_y_and_z_axes(commondata, grid, params, xx, NUM_nearest_GFS_1d, which_gfs_1d,
                                        diagnostic_gf_names, gridfuncs_diags);

    const int NUM_nearest_GFS_2d = (int)(sizeof which_gfs_2d / sizeof which_gfs_2d[0]);
    diagnostics_nearest_2d_xy_and_yz_planes(commondata, grid, params, xx, NUM_nearest_GFS_2d, which_gfs_2d,
                                            diagnostic_gf_names, gridfuncs_diags);
"""

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "diagnostic_gfs.h",
        "diagnostics_nearest_common.h",
    ]
    desc = """
Dispatch GRoovy nearest-sampled diagnostics and write hydro-center data.

The diagnostics_nearest() dispatcher preserves the standard nearest sampled
constraint diagnostics and adds a GRHD-specific 0D output file. The GRHD output
reads primitive fields from auxevol_gfs and conserved fields from y_n_gfs at the
grid point nearest the physical center.

@param[in] commondata Global simulation metadata.
@param[in] griddata Per-grid data structures.
@param[in] gridfuncs_diags Per-grid diagnostic gridfunction arrays used by the generic nearest constraint helpers.
"""
    cfunc_type = "void"
    name = "diagnostics_nearest"
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                const REAL *restrict gridfuncs_diags[MAXNUMGRIDS]"""
    body = f"""
  // --- USER-EDIT: Select constraint diagnostics to sample (applies to all grids) ---
{constraint_selection}
  // --- END USER-EDIT ---

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
    const params_struct *restrict params = &griddata[grid].params;
    const REAL *restrict xx[3] = {{griddata[grid].xx[0], griddata[grid].xx[1], griddata[grid].xx[2]}};
{constraint_calls}
    int i0_center = -1;
    int i1_center = -1;
    int i2_center = -1;
    switch (params->CoordSystem_hash) {{{center_index_cases}
    default:
      fprintf(stderr, "ERROR in diagnostics_nearest(): CoordSystem hash = %d not #define'd!\\n", params->CoordSystem_hash);
      exit(1);
    }} // END SWITCH: choose nearest center indices

    const int idx3 = IDX3P(params, i0_center, i1_center, i2_center);
    const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    const REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;

    char hydro_coordsys_with_grid[128];
    snprintf(hydro_coordsys_with_grid, sizeof(hydro_coordsys_with_grid), "grid%02d-%s", grid, params->CoordSystemName);

    FILE *out_hydro = open_outfile("out0d-hydro", hydro_coordsys_with_grid, commondata, /*include_time=*/0);
    if (out_hydro != NULL) {{
      if (commondata->nn == 0)
        fprintf(out_hydro, "{hydro_header}\\n");

      const int NUM_HYDRO_CENTER_COLS = {num_hydro_center_cols};
      REAL hydro_center_row[NUM_HYDRO_CENTER_COLS];
      hydro_center_row[0] = commondata->time;
{hydro_row_assignments}
      diag_write_row(out_hydro, NUM_HYDRO_CENTER_COLS, hydro_center_row);
      fclose(out_hydro);
    }} // END IF: hydro-center output file opened
  }} // END LOOP: for grid over active grids
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
