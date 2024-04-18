"""
C functions for diagnostics for the superB infrastructure

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from typing import List, Union, cast, Tuple, Dict
from inspect import currentframe as cfr
from types import FrameType as FT

import nrpy.params as par
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.superB.output_0d_1d_2d_nearest_gridpoint_slices as out012d
from nrpy.infrastructures.BHaH import griddata_commondata


def register_CFunction_diagnostics(
    list_of_CoordSystems: List[str],
    default_diagnostics_out_every: float,
    enable_psi4_diagnostics: bool = False,
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

    :param list_of_CoordSystems: Lists of unique CoordSystems used.
    :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
    :param enable_psi4_diagnostics: Whether or not to enable psi4 diagnostics.
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

    for CoordSystem in list_of_CoordSystems:
        out012d.register_CFunction_diagnostics_nearest_grid_center(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            filename_tuple=grid_center_filename_tuple,
        )
        for axis in ["y", "z"]:
            out012d.register_CFunction_diagnostics_nearest_1d_axis(
                CoordSystem=CoordSystem,
                out_quantities_dict=out_quantities_dict,
                axis=axis,
            )
            out012d.register_CFunction_diagnostics_set_up_nearest_1d_axis(
                CoordSystem=CoordSystem,
                out_quantities_dict=out_quantities_dict,
                filename_tuple=axis_filename_tuple,
                axis=axis,
            )

        for plane in ["xy", "yz"]:
            out012d.register_CFunction_diagnostics_nearest_2d_plane(
                CoordSystem=CoordSystem,
                out_quantities_dict=out_quantities_dict,
                plane=plane,
            )
            out012d.register_CFunction_diagnostics_set_up_nearest_2d_plane(
                CoordSystem=CoordSystem,
                out_quantities_dict=out_quantities_dict,
                filename_tuple=plane_filename_tuple,
                plane=plane,
            )

    desc = r"""Diagnostics."""
    cfunc_type = "void"
    name = "diagnostics"
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata, Ck::IO::Session token, const int which_output, const int grid"

    body = r"""


const int num_diagnostic_1d_y_pts = griddata[grid].diagnosticstruct.num_diagnostic_1d_y_pts;
const int num_diagnostic_1d_z_pts = griddata[grid].diagnosticstruct.num_diagnostic_1d_z_pts;
const int num_diagnostic_2d_xy_pts = griddata[grid].diagnosticstruct.num_diagnostic_2d_xy_pts;
const int num_diagnostic_2d_yz_pts = griddata[grid].diagnosticstruct.num_diagnostic_2d_yz_pts;


const bool write_diagnostics = (which_output == OUTPUT_0D) ||
					(num_diagnostic_1d_y_pts > 0) ||
					(num_diagnostic_1d_z_pts > 0) ||
					(num_diagnostic_2d_xy_pts > 0) ||
					(num_diagnostic_2d_yz_pts > 0);

if (write_diagnostics) {
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
    Ricci_eval(commondata, params, &griddata[grid].rfmstruct, y_n_gfs, auxevol_gfs);
    constraints_eval(commondata, params, &griddata[grid].rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_output_gfs);
  }

  // // 0D, 1D and 2D outputs
  if (which_output == OUTPUT_0D) {
      diagnostics_nearest_grid_center(commondata, params, &griddata[grid].gridfuncs);
  } else if (which_output == OUTPUT_1D_Y) {
    if (num_diagnostic_1d_y_pts > 0) {
    diagnostics_nearest_1d_y_axis(commondata, params, xx, &griddata[grid].gridfuncs, &griddata[grid].diagnosticstruct, token);
    }
  } else if (which_output == OUTPUT_1D_Z) {
    if (num_diagnostic_1d_z_pts > 0) {
    diagnostics_nearest_1d_z_axis(commondata, params, xx, &griddata[grid].gridfuncs, &griddata[grid].diagnosticstruct, token);
    }
  } else if (which_output == OUTPUT_2D_XY) {
    if (num_diagnostic_2d_xy_pts > 0) {
    diagnostics_nearest_2d_xy_plane(commondata, params, xx, &griddata[grid].gridfuncs, &griddata[grid].diagnosticstruct, token);
    }
  } else if (which_output == OUTPUT_2D_YZ) {
    if (num_diagnostic_2d_yz_pts > 0) {
    diagnostics_nearest_2d_yz_plane(commondata, params, xx, &griddata[grid].gridfuncs, &griddata[grid].diagnosticstruct, token);
    }
  }
"""
    if enable_psi4_diagnostics:
        body += r"""      // Do psi4 output, but only if the grid is spherical-like.
    if (strstr(CoordSystemName, "Spherical") != NULL) {
      // todo
    }
"""
    body += r"""
}

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

    # Register diagnostic_struct's contribution to griddata_struct:
    griddata_commondata.register_griddata_commondata(
        __name__,
        "diagnostic_struct diagnosticstruct",
        "store indices of 1d and 2d diagnostic points, the offset in the output file, etc",
    )

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
