"""
C function for wave equation evolution diagnostics.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, Set, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.diagnostics.output_0d_1d_2d_nearest_gridpoint_slices as nearest012d
import nrpy.params as par


def register_CFunction_diagnostics(
    set_of_CoordSystems: Set[str],
    default_diagnostics_out_every: float,
    enable_progress_indicator: bool = True,
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
    :param enable_progress_indicator: Whether to enable the progress indicator.
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
    desc = r"""Diagnostics."""
    cfunc_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    if "uuexact" not in gri.glb_gridfcs_dict:
        _ = gri.register_gridfunctions(
            "uuexact", group="AUX", gf_array_name="diagnostic_output_gfs"
        )
    if "vvexact" not in gri.glb_gridfcs_dict:
        _ = gri.register_gridfunctions(
            "vvexact", group="AUX", gf_array_name="diagnostic_output_gfs"
        )

    # fmt: off
    if out_quantities_dict == "default":
        out_quantities_dict = {
            ("REAL", "log10ErelUU"): "log10(fabs((y_n_gfs[IDX4pt(UUGF, idx3)]-diagnostic_output_gfs[IDX4pt(UUEXACTGF, idx3)])/(diagnostic_output_gfs[IDX4pt(UUEXACTGF, idx3)] + 1e-16)) + 1e-16)",
            ("REAL", "log10ErelVV"): "log10(fabs((y_n_gfs[IDX4pt(VVGF, idx3)]-diagnostic_output_gfs[IDX4pt(VVEXACTGF, idx3)])/(diagnostic_output_gfs[IDX4pt(VVEXACTGF, idx3)] + 1e-16)) + 1e-16)",
            ("REAL", "exactUU"): "diagnostic_output_gfs[IDX4pt(UUEXACTGF, idx3)]",
            ("REAL", "numUU"): "y_n_gfs[IDX4pt(UUGF, idx3)]",
            ("REAL", "exactVV"): "diagnostic_output_gfs[IDX4pt(VVEXACTGF, idx3)]",
            ("REAL", "numVV"): "y_n_gfs[IDX4pt(VVGF, idx3)]",
        }
    if not isinstance(out_quantities_dict, dict):
        raise TypeError(f"out_quantities_dict was initialized to {out_quantities_dict}, which is not a dictionary!")
    # fmt: on

    for CoordSystem in set_of_CoordSystems:
        nearest012d.register_CFunction_diagnostics_nearest_grid_center(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            filename_tuple=grid_center_filename_tuple,
        )
        for axis in ["y", "z"]:
            nearest012d.register_CFunction_diagnostics_nearest_1d_axis(
                CoordSystem=CoordSystem,
                out_quantities_dict=out_quantities_dict,
                filename_tuple=axis_filename_tuple,
                axis=axis,
            )
        for plane in ["xy", "yz"]:
            nearest012d.register_CFunction_diagnostics_nearest_2d_plane(
                CoordSystem=CoordSystem,
                out_quantities_dict=out_quantities_dict,
                filename_tuple=plane_filename_tuple,
                plane=plane,
            )
    # fmt: on

    body = r"""  const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->diagnostics_output_every;
  // Explanation of the if() below:
  // Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime.
  // Step 2: Multiplying by outevery yields the exact time we should output again, t_out.
  // Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible!
  if (fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt) {
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
      // Unpack griddata struct:
      MAYBE_UNUSED const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
      REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;
      REAL *restrict xx[3];
      for (int ww = 0; ww < 3; ww++)
        xx[ww] = griddata[grid].xx[ww];
      const params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

      LOOP_OMP("omp parallel for", i0, NGHOSTS, Nxx0 + NGHOSTS, i1, NGHOSTS, Nxx1 + NGHOSTS, i2, NGHOSTS, Nxx2 + NGHOSTS) {
        REAL xCart[3];
        REAL xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
        xx_to_Cart(params, xOrig, xCart);
        REAL *restrict uexact = &diagnostic_output_gfs[IDX4(UUEXACTGF, i0,i1,i2)];
        REAL *restrict vexact = &diagnostic_output_gfs[IDX4(VVEXACTGF, i0,i1,i2)];
        exact_solution_single_Cartesian_point(commondata, params, xCart[0], xCart[1], xCart[2], uexact, vexact);
      }

      // 0D output
      diagnostics_nearest_grid_center(commondata, params, &griddata[grid].gridfuncs);

      // 1D output
      diagnostics_nearest_1d_y_axis(commondata, params, xx, &griddata[grid].gridfuncs);
      diagnostics_nearest_1d_z_axis(commondata, params, xx, &griddata[grid].gridfuncs);

      // 2D output
      diagnostics_nearest_2d_xy_plane(commondata, params, xx, &griddata[grid].gridfuncs);
      diagnostics_nearest_2d_yz_plane(commondata, params, xx, &griddata[grid].gridfuncs);
    }
  }
"""
    if enable_progress_indicator:
        body += "progress_indicator(commondata, griddata);"
    body += r"""
  if (commondata->time + commondata->dt > commondata->t_final)
    printf("\n");
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
