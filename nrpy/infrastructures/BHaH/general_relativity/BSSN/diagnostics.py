"""
Generate C functions for computing BSSN diagnostics in curvilinear coordinates.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, Set, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.diagnostics.output_0d_1d_2d_nearest_gridpoint_slices as out012d
import nrpy.params as par
from nrpy.infrastructures.BHaH import BHaH_defines_h, griddata_commondata
from nrpy.infrastructures.superB.CurviBoundaryConditions import (
    register_CFunction_apply_bcs_inner_only_specific_gfs,
    register_CFunction_apply_bcs_outerextrap_and_inner_specific_gfs,
)


def register_CFunction_psi4_diagnostics_set_up(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for setting up diagnostic struct.

    :param CoordSystem: Specifies the coordinate system for psi4 diagnostics set up.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If psi4 decomposition is not supported for the coordinate system.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = "Setup psi4_diagnostics."
    cfunc_type = "void"
    name = "psi4_diagnostics_set_up"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], diagnostic_struct *restrict diagnosticstruct"

    if "Spherical" in CoordSystem:
        choose_number_pts_in_theta_shell = "const int N_theta = params->Nxx1;"
        choose_number_pts_in_phi_shell = "const int N_phi = params->Nxx2;"
    elif "Cylindrical" in CoordSystem:
        choose_number_pts_in_theta_shell = "const int N_theta = params->Nxx2;"
        choose_number_pts_in_phi_shell = "const int N_phi = params->Nxx1;"
    elif "SymTP" in CoordSystem:
        choose_number_pts_in_theta_shell = "const int N_theta = params->Nxx1;"
        choose_number_pts_in_phi_shell = "const int N_phi = params->Nxx2;"
    else:
        raise ValueError(f"CoordSystem = {CoordSystem} not supported.")

    body = rf"""
  // Number of points on shell in theta direction
  {choose_number_pts_in_theta_shell}
  // Number of points on shell in phi direction
  {choose_number_pts_in_phi_shell}
"""
    body += r"""
const int psi4_spinweightm2_sph_harmonics_max_l = commondata->swm2sh_maximum_l_mode_to_compute;
#define NUM_OF_R_EXTS 24
  const REAL list_of_R_exts[NUM_OF_R_EXTS] = {10.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,  30.0,
                                              31.0, 32.0, 33.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0};


  // Set up uniform 2d grid in theta and phi at R_ext (2d shells at different R_ext)
  // phi values need to be exactly the same as phi values of  grid
  // # of pts in the theta direction can be chosen freely, here it set to the same as number of points in the theta-like direction of the grid
  const REAL PI = params->PI;
  const REAL theta_min = 0.0;
  REAL theta_max = PI;
  REAL phi_min = -PI;
  REAL phi_max = PI;
  const int N_tot_shell = N_theta * N_phi;
  REAL dtheta = (theta_max - theta_min) / N_theta;
  diagnosticstruct->dtheta = dtheta;
  const REAL dphi = (phi_max - phi_min) / N_phi;
  REAL ***restrict xx_shell_sph = (REAL * **restrict)malloc(NUM_OF_R_EXTS * sizeof(REAL **));
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {
    xx_shell_sph[which_R_ext] = (REAL * *restrict)malloc(2 * sizeof(REAL *));
    xx_shell_sph[which_R_ext][0] = (REAL *restrict)malloc(N_theta * sizeof(REAL));
    xx_shell_sph[which_R_ext][1] = (REAL *restrict)malloc(N_phi * sizeof(REAL));
  }
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {
    for (int j = 0; j < N_theta; j++) {
      xx_shell_sph[which_R_ext][0][j] = theta_min + ((REAL)j + 0.5) * dtheta;
    }
    for (int j = 0; j < N_phi; j++) {
      xx_shell_sph[which_R_ext][1][j] = phi_min + ((REAL)j + 0.5) * dphi;
    }
  }

  // Convert points on shells to Cartesian coordinates
  REAL ***restrict xx_shell_Cart = (REAL * **restrict)malloc(NUM_OF_R_EXTS * sizeof(REAL **));
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {
    xx_shell_Cart[which_R_ext] = (REAL * *restrict)malloc(3 * sizeof(REAL *));
    xx_shell_Cart[which_R_ext][0] = (REAL *restrict)malloc(N_tot_shell * sizeof(REAL));
    xx_shell_Cart[which_R_ext][1] = (REAL *restrict)malloc(N_tot_shell * sizeof(REAL));
    xx_shell_Cart[which_R_ext][2] = (REAL *restrict)malloc(N_tot_shell * sizeof(REAL));
  }
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {
    REAL r = list_of_R_exts[which_R_ext];
    for (int i_ph = 0; i_ph < N_phi; i_ph++) {
      REAL phi = xx_shell_sph[which_R_ext][1][i_ph];
      for (int i_th = 0; i_th < N_theta; i_th++) {
        REAL theta = xx_shell_sph[which_R_ext][0][i_th];
        REAL x = r * sin(theta) * cos(phi);
        REAL y = r * sin(theta) * sin(phi);
        REAL z = r * cos(theta);
        const int idx2 = IDX2GENERAL(i_th, i_ph, N_theta);
        xx_shell_Cart[which_R_ext][0][idx2] = x;
        xx_shell_Cart[which_R_ext][1][idx2] = y;
        xx_shell_Cart[which_R_ext][2][idx2] = z;
      }
    }
  }

  // For each pt on each spherical shell at R_ext find if pt lies within the grid
  //(it should for bhah single grids, but not for multipatch...)
  // set N_shell_pts_grid and xx_shell_grid in diagnostic struct
  diagnosticstruct->N_shell_pts_grid = (int *restrict)malloc(sizeof(int) * NUM_OF_R_EXTS);
  diagnosticstruct->N_theta_shell_grid = (int *restrict)malloc(sizeof(int) * NUM_OF_R_EXTS);
  diagnosticstruct->xx_shell_grid = (REAL * **restrict)malloc(NUM_OF_R_EXTS * sizeof(REAL **));
  diagnosticstruct->theta_shell_grid = (REAL * *restrict)malloc(NUM_OF_R_EXTS * sizeof(REAL *));

  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {
    // Count number of pts that lie within the extent of this grid
    int count_pt_on_grid = 0;
    bool b_theta_grid[N_theta];
    bool b_phi_grid[N_phi];
    for (int i = 0; i < N_theta; i++) {
        b_theta_grid[i] = false;
    }
    for (int j = 0; j < N_phi; j++) {
        b_phi_grid[j] = false;
    }
    for (int i_ph = 0; i_ph < N_phi; i_ph++) {
      for (int i_th = 0; i_th < N_theta; i_th++) {
        const int idx2 = IDX2GENERAL(i_th, i_ph, N_theta);
        const REAL xCart_pt_on_shell[3] = {xx_shell_Cart[which_R_ext][0][idx2], xx_shell_Cart[which_R_ext][1][idx2],
                                           xx_shell_Cart[which_R_ext][2][idx2]};
        int Cart_to_i0i1i2[3];
        REAL closest_xx[3];
        Cart_to_xx_and_nearest_i0i1i2(params, xCart_pt_on_shell, closest_xx, Cart_to_i0i1i2); // convert from Cart to xx
        if ((params->xxmin0 <= closest_xx[0] && closest_xx[0] <= params->xxmax0) &&
            (params->xxmin1 <= closest_xx[1] && closest_xx[1] <= params->xxmax1) &&
            (params->xxmin2 <= closest_xx[2] && closest_xx[2] <= params->xxmax2)) {
          count_pt_on_grid++;
          b_theta_grid[i_th] = true;
          b_phi_grid[i_ph] = true;
        }
      }
    }
    // Count the number of unique theta values on grid
    int count_theta_grid = 0;
    for (int i = 0; i < N_theta; i++) {
      if (b_theta_grid[i]) {
        count_theta_grid++;
      }
    }
    // Count the number of unique phi values on grid (should be the same as Nxx1)
    int count_phi_grid = 0;
    for (int i = 0; i < N_phi; i++) {
      if (b_phi_grid[i]) {
        count_phi_grid++;
      }
    }
    // Set values in diagnosticstruct
    diagnosticstruct->N_shell_pts_grid[which_R_ext] = count_pt_on_grid;
    diagnosticstruct->N_theta_shell_grid[which_R_ext] = count_theta_grid;

    // Allocate memory after counting
    diagnosticstruct->xx_shell_grid[which_R_ext] = (REAL * *restrict)malloc(count_pt_on_grid * sizeof(REAL *));
    diagnosticstruct->theta_shell_grid[which_R_ext] = (REAL *restrict)malloc(count_theta_grid * sizeof(REAL));
    for (int i = 0; i < count_pt_on_grid; i++) {
      diagnosticstruct->xx_shell_grid[which_R_ext][i] = (REAL *restrict)malloc(3 * sizeof(REAL));
    }

    // Now set them
    int which_pt_on_grid = 0;
    for (int i_ph = 0; i_ph < N_phi; i_ph++) {
      for (int i_th = 0; i_th < N_theta; i_th++) {
        const int idx2 = IDX2GENERAL(i_th, i_ph, N_theta);
        const REAL xCart_pt_on_shell[3] = {xx_shell_Cart[which_R_ext][0][idx2], xx_shell_Cart[which_R_ext][1][idx2],
                                           xx_shell_Cart[which_R_ext][2][idx2]};
        int Cart_to_i0i1i2[3];
        REAL closest_xx[3];
        Cart_to_xx_and_nearest_i0i1i2(params, xCart_pt_on_shell, closest_xx, Cart_to_i0i1i2);
        if ((params->xxmin0 <= closest_xx[0] && closest_xx[0] <= params->xxmax0) &&
            (params->xxmin1 <= closest_xx[1] && closest_xx[1] <= params->xxmax1) &&
            (params->xxmin2 <= closest_xx[2] && closest_xx[2] <= params->xxmax2)) {

          diagnosticstruct->xx_shell_grid[which_R_ext][which_pt_on_grid][0] = closest_xx[0];
          diagnosticstruct->xx_shell_grid[which_R_ext][which_pt_on_grid][1] = closest_xx[1];
          diagnosticstruct->xx_shell_grid[which_R_ext][which_pt_on_grid][2] = closest_xx[2];

          // also save theta values
          int i_th_grid, i_ph_grid;
          const int N_theta_shell_grid = diagnosticstruct->N_theta_shell_grid[which_R_ext];
          REVERSE_IDX2GENERAL(which_pt_on_grid, N_theta_shell_grid, i_th_grid, i_ph_grid);
          diagnosticstruct->theta_shell_grid[which_R_ext][i_th_grid] = xx_shell_sph[which_R_ext][0][i_th];
          which_pt_on_grid++;
        }
      }
    }
  } // end loop over all R_ext

  diagnosticstruct->num_of_R_exts_grid = NUM_OF_R_EXTS;
  diagnosticstruct->list_of_R_exts_grid = (REAL *restrict)malloc(sizeof(REAL) * NUM_OF_R_EXTS);
  for (int i = 0; i < NUM_OF_R_EXTS; i++) {
    diagnosticstruct->list_of_R_exts_grid[i] = list_of_R_exts[i];
  }
  diagnosticstruct->psi4_spinweightm2_sph_harmonics_max_l = psi4_spinweightm2_sph_harmonics_max_l;
  int sum = 0;
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {
    if (diagnosticstruct->N_shell_pts_grid[which_R_ext] > 0) {
      sum += diagnosticstruct->N_shell_pts_grid[which_R_ext];
    }
  }
  diagnosticstruct->tot_N_shell_pts_grid = sum;

  // Free memory
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {
    free(xx_shell_sph[which_R_ext][0]);
    free(xx_shell_sph[which_R_ext][1]);
    free(xx_shell_sph[which_R_ext]);
  }
  free(xx_shell_sph);
  for (int which_R_ext = 0; which_R_ext < NUM_OF_R_EXTS; which_R_ext++) {
    free(xx_shell_Cart[which_R_ext][0]);
    free(xx_shell_Cart[which_R_ext][1]);
    free(xx_shell_Cart[which_R_ext][2]);
    free(xx_shell_Cart[which_R_ext]);
  }
  free(xx_shell_Cart);
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_griddata() -> None:
    """Register the diagnostic_struct's contribution to the griddata_struct."""
    griddata_commondata.register_griddata_commondata(
        __name__,
        "diagnostic_struct diagnosticstruct",
        "data needed to do psi4 decomposition in cylindrical-like coordinates",
    )


def register_BHaH_defines_h() -> None:
    """Register the diagnostic_struct's contribution to the BHaH_defines.h file."""
    BHd_str = r"""
    // for psi4 decomposition
    typedef struct __diagnostic_struct__ {
      int num_of_R_exts_grid;
      int psi4_spinweightm2_sph_harmonics_max_l;
      REAL *restrict list_of_R_exts_grid;
      int tot_N_shell_pts_grid;
      REAL dtheta;
      int *restrict N_shell_pts_grid; // of shape int [num_of_R_exts_grid]
      int *restrict N_theta_shell_grid; // of shape int [num_of_R_exts_grid]
      REAL ***restrict xx_shell_grid; // of shape [num_of_R_exts_grid][N_shell_pts_grid][3]
      REAL **restrict theta_shell_grid; // of shape [num_of_R_exts_grid][N_theta_shell_grid]
    } diagnostic_struct;
    """

    BHaH_defines_h.register_BHaH_defines(
        __name__,
        BHd_str,
    )


def register_CFunction_diagnostics(
    set_of_CoordSystems: Set[str],
    default_diagnostics_out_every: float,
    enable_psi4_diagnostics: bool = False,
    enable_progress_indicator: bool = True,
    use_Ricci_eval_func: bool = True,
    grid_center_filename_tuple: Tuple[str, str] = (
        "out0d-conv_factor%.2f.txt",
        "convergence_factor",
    ),
    axis_filename_tuple: Tuple[str, str] = (
        "out1d-AXIS-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
    plane_filename_tuple: Tuple[str, str] = (
        "out2d-PLANE-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
    out_quantities_dict: Union[str, Dict[Tuple[str, str], str]] = "default",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for simulation diagnostics.

    :param set_of_CoordSystems: Sets of unique CoordSystems used.
    :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
    :param enable_psi4_diagnostics: Whether to enable psi4 diagnostics.
    :param enable_progress_indicator: Whether to enable the progress indicator.
    :param use_Ricci_eval_func: Whether to call Ricci_eval() before computing constraints.
    :param grid_center_filename_tuple: Tuple containing filename and variables for grid center output.
    :param axis_filename_tuple: Tuple containing filename and variables for axis output.
    :param plane_filename_tuple: Tuple containing filename and variables for plane output.
    :param out_quantities_dict: Dictionary or string specifying output quantities.

    :return: None if in registration phase, else the updated NRPy environment.
    :raises TypeError: If `out_quantities_dict` is not a dictionary and not set to "default".
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    _ = par.CodeParameter(
        "REAL",
        __name__,
        "diagnostics_output_every",
        default_diagnostics_out_every,
        commondata=True,
    )
    parallelization = par.parval_from_str("parallelization")

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    # fmt: off
    if out_quantities_dict == "default":
        out_quantities_dict = {
            ("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))",
            ("REAL", "log10sqrtM2L"): "log10(sqrt(diagnostic_output_gfs[IDX4pt(MSQUAREDGF, idx3)]) + 1e-16)",
            ("REAL", "cfL"): "y_n_gfs[IDX4pt(CFGF, idx3)]",
            ("REAL", "alphaL"): "y_n_gfs[IDX4pt(ALPHAGF, idx3)]",
            ("REAL", "trKL"): "y_n_gfs[IDX4pt(TRKGF, idx3)]",
        }
    if not isinstance(out_quantities_dict, dict):
        raise TypeError(f"out_quantities_dict was initialized to {out_quantities_dict}, which is not a dictionary!")
    # fmt: on
    out_quantities_gf_indexes_dict = {
        tmp_str.replace(" ", "")
        .split("[")[1]
        .split(",")[0]
        .replace("IDX4pt(", ""): gf_ptr
        for gf_ptr in ["y_n_gfs", "diagnostic_output_gfs"]
        for v in out_quantities_dict.values()
        for tmp_str in v.split(gf_ptr)
        if "IDX4pt" in tmp_str and tmp_str and tmp_str[0] == "["
    }

    for CoordSystem in set_of_CoordSystems:
        out012d.register_CFunction_diagnostics_nearest_grid_center(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            filename_tuple=grid_center_filename_tuple,
        )
        for axis in ["y", "z"]:
            out012d.register_CFunction_diagnostics_nearest_1d_axis(
                CoordSystem=CoordSystem,
                out_quantities_dict=out_quantities_dict,
                filename_tuple=axis_filename_tuple,
                axis=axis,
            )
        for plane in ["xy", "yz"]:
            out012d.register_CFunction_diagnostics_nearest_2d_plane(
                CoordSystem=CoordSystem,
                out_quantities_dict=out_quantities_dict,
                filename_tuple=plane_filename_tuple,
                plane=plane,
            )

    # Register diagnostic_struct's contribution to griddata_struct:
    register_griddata()

    # Register diagnostic_struct's contribution to BHaH_defines.h:
    register_BHaH_defines_h()

    # Register psi4_diagnostics_set_up and other functions needed for psi4 decomposition
    if enable_psi4_diagnostics:
        for CoordSystem in set_of_CoordSystems:
            register_CFunction_psi4_diagnostics_set_up(CoordSystem=CoordSystem)
        register_CFunction_apply_bcs_inner_only_specific_gfs()
        register_CFunction_apply_bcs_outerextrap_and_inner_specific_gfs()

    desc = r"""Diagnostics."""
    cfunc_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
        + (
            ", griddata_struct *restrict griddata_host"
            if parallelization == "cuda"
            else ""
        )
    )

    host_griddata = "griddata_host" if parallelization == "cuda" else "griddata"
    body = rf"""
const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->diagnostics_output_every;
// Explanation of the if() below:
// Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime.
// Step 2: Multiplying by outevery yields the exact time we should output again, t_out.
// Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible!
if(fabs(round(currtime / outevery) * outevery - currtime) < 0.5*currdt) {{
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
    // Unpack griddata struct:
    const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;
    REAL *restrict xx[3];
    {{
      for (int ww = 0; ww < 3; ww++)
        xx[ww] = {host_griddata}[grid].xx[ww];
    }}
    const params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
"""
    if parallelization == "cuda":
        body += r"""
    // This does not leverage async memory transfers using multiple streams at the moment
    // given the current intent is one cuda stream per grid. This could be leveraged
    // in the future by increasing NUM_STREAMS such that a diagnostic stream is included per grid
    size_t streamid = params->grid_idx % NUM_STREAMS;
    cpyHosttoDevice_params__constant(&griddata[grid].params, streamid);
    REAL *restrict host_y_n_gfs = griddata_host[grid].gridfuncs.y_n_gfs;
    REAL *restrict host_diagnostic_output_gfs = griddata_host[grid].gridfuncs.diagnostic_output_gfs;
"""
        for idx, gf in out_quantities_gf_indexes_dict.items():
            if "y_n_gfs" in gf:
                body += f"    cpyDevicetoHost__gf(commondata, params, host_{gf}, {gf}, {idx}, {idx}, streamid);\n"
    body += r"""
    // Constraint output
    {
"""
    if use_Ricci_eval_func:
        body += "Ricci_eval(params, griddata[grid].rfmstruct, y_n_gfs, auxevol_gfs);\n"
    body += r"""
      constraints_eval(params, griddata[grid].rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_output_gfs);
    }"""

    if parallelization == "cuda":
        for idx, gf in out_quantities_gf_indexes_dict.items():
            if "diagnostic_output_gfs" in gf:
                body += f"    cpyDevicetoHost__gf(commondata, params, host_{gf}, {gf}, {idx}, {idx}, streamid);\n"
            body += "cudaStreamSynchronize(streams[streamid]);"

    body += f"""
    // 0D output
    diagnostics_nearest_grid_center(commondata, params, &{host_griddata}[grid].gridfuncs);

    // 1D output
    diagnostics_nearest_1d_y_axis(commondata, params, xx, &{host_griddata}[grid].gridfuncs);
    diagnostics_nearest_1d_z_axis(commondata, params, xx, &{host_griddata}[grid].gridfuncs);

    // 2D output
    diagnostics_nearest_2d_xy_plane(commondata, params, xx, &{host_griddata}[grid].gridfuncs);
    diagnostics_nearest_2d_yz_plane(commondata, params, xx, &{host_griddata}[grid].gridfuncs);
"""
    if enable_psi4_diagnostics:
        body += r"""      // Do psi4 output
    // Set psi4.
    psi4(commondata, params, xx, y_n_gfs, diagnostic_output_gfs);

    // Apply outer and inner bcs to psi4 needed to do interpolation correctly
    int aux_gfs_to_sync[2] = {PSI4_REGF, PSI4_IMGF};
    apply_bcs_outerextrap_and_inner_specific_gfs(commondata, &griddata[grid].params, &griddata[grid].bcstruct, 2, griddata[grid].gridfuncs.diagnostic_output_gfs, aux_gfs_to_sync, aux_gf_parity);

    // Decompose psi4 into spin-weight -2  spherical harmonics & output to files.
    diagnostic_struct *restrict diagnosticstruct = &griddata[grid].diagnosticstruct;

    psi4_spinweightm2_decomposition(commondata, params, diagnostic_output_gfs, diagnosticstruct->list_of_R_exts_grid, diagnosticstruct->num_of_R_exts_grid,
    diagnosticstruct->psi4_spinweightm2_sph_harmonics_max_l, xx, diagnosticstruct->N_shell_pts_grid,
    diagnosticstruct->xx_shell_grid, diagnosticstruct->N_theta_shell_grid, diagnosticstruct->theta_shell_grid, diagnosticstruct->dtheta);
"""
    body += r"""
  }
}
"""
    if enable_progress_indicator:
        body += "progress_indicator(commondata, griddata);"
    body += r"""
if(commondata->time + commondata->dt > commondata->t_final) printf("\n");
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
