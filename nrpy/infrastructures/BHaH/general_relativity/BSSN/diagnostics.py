"""
Generate C functions for computing BSSN diagnostics in curvilinear coordinates.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, Set, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.diagnostics.output_0d_1d_2d_nearest_gridpoint_slices as out012d
import nrpy.params as par


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

    desc = r"""Diagnostics."""
    cfunc_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    body = r"""
const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->diagnostics_output_every;
// Explanation of the if() below:
// Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime.
// Step 2: Multiplying by outevery yields the exact time we should output again, t_out.
// Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible!
if(fabs(round(currtime / outevery) * outevery - currtime) < 0.5*currdt) {
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    // Unpack griddata struct:
    const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;
    REAL *restrict xx[3];
    {
      for (int ww = 0; ww < 3; ww++)
        xx[ww] = griddata[grid].xx[ww];
    }
    const params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

    // Constraint output
    {
"""
    if use_Ricci_eval_func:
        body += "Ricci_eval(commondata, params, griddata[grid].rfmstruct, y_n_gfs, auxevol_gfs);\n"
    body += r"""
      constraints_eval(commondata, params, griddata[grid].rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_output_gfs);
    }

    // 0D output
    diagnostics_nearest_grid_center(commondata, params, &griddata[grid].gridfuncs);

    // 1D output
    diagnostics_nearest_1d_y_axis(commondata, params, xx, &griddata[grid].gridfuncs);
    diagnostics_nearest_1d_z_axis(commondata, params, xx, &griddata[grid].gridfuncs);

    // 2D output
    diagnostics_nearest_2d_xy_plane(commondata, params, xx, &griddata[grid].gridfuncs);
    diagnostics_nearest_2d_yz_plane(commondata, params, xx, &griddata[grid].gridfuncs);
"""
    if enable_psi4_diagnostics:
        body += r"""      // Do psi4 output, but only if the grid is spherical-like.
      if (strstr(CoordSystemName, "Spherical") != NULL) {

        // Adjusted to match Tutorial-Start_to_Finish-BSSNCurvilinear-Two_BHs_Collide-Psi4.ipynb
        const int psi4_spinweightm2_sph_harmonics_max_l = 2;
#define num_of_R_exts 24
        const REAL list_of_R_exts[num_of_R_exts] =
        { 10.0, 20.0, 21.0, 22.0, 23.0,
          24.0, 25.0, 26.0, 27.0, 28.0,
          29.0, 30.0, 31.0, 32.0, 33.0,
          35.0, 40.0, 50.0, 60.0, 70.0,
          80.0, 90.0, 100.0, 150.0 };

        // Set psi4.
        psi4(commondata, params, xx, y_n_gfs, diagnostic_output_gfs);
        // Decompose psi4 into spin-weight -2  spherical harmonics & output to files.
        psi4_spinweightm2_decomposition_on_sphlike_grids(commondata, params, diagnostic_output_gfs, list_of_R_exts, num_of_R_exts, psi4_spinweightm2_sph_harmonics_max_l, xx);
      }
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
