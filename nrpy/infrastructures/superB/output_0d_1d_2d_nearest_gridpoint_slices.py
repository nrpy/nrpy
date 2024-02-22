"""
Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from typing import Dict, Tuple, Union, cast
from inspect import currentframe as cfr
from types import FrameType as FT

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp


def register_CFunction_diagnostics_nearest_1d_axis(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    axis: str,
    filename_tuple: Tuple[str, str] = (
        "out1d-AXIS-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
) -> Union[None, pcg.NRPyEnv_type]:

    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    if axis not in ["y", "z"]:
        raise ValueError(
            f"Output along {axis} axis not supported. Please choose x or z axis."
        )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    prefunc, loop_1d = lp.simple_loop_1D(
        CoordSystem=CoordSystem,
        out_quantities_dict=out_quantities_dict,
        axis=axis,
    )

    desc = f"Output diagnostic quantities at gridpoints closest to {axis} axis."
    c_type = "void"
    name = f"diagnostics_nearest_1d_{axis}_axis"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs, const diagnostic_pt_offset_struct *restrict diagnosticptoffsetstruct, Ck::IO::Session token"

    body = rf"""
// Unpack gridfuncs struct:
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

const int num_diagnostic_1d_y_pts = diagnosticptoffsetstruct->num_diagnostic_1d_y_pts;

for (int which_pt = 0; which_pt < num_diagnostic_1d_y_pts; which_pt++) {{
  const int idx3 = diagnosticptoffsetstruct->localidx3_diagnostic_1d_y_pt[which_pt];
  const int i0 = diagnosticptoffsetstruct->locali0_diagnostic_1d_y_pt[which_pt];
  const int i1 = diagnosticptoffsetstruct->locali1_diagnostic_1d_y_pt[which_pt];
  const int i2 = diagnosticptoffsetstruct->locali2_diagnostic_1d_y_pt[which_pt];
  REAL xCart[3];
  xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
  const REAL xCart_axis = {'xCart[1];' if axis == "y" else 'xCart[2];'}
  int sizeinbytes = 23 * (len(out_quantities_dict) + 1);
  char out[sizeinbytes+1];
  snprintf(out, sizeof(out), " """
    output_to_file = "%.15e,"
    for key in out_quantities_dict.keys():
      printf_c_type = "%.15e" if key[0] != "int" else "%d"
      output_to_file += f"{printf_c_type} "

    output_to_file = (
      f'{output_to_file[:-1]}\\n", xCart_axis, '
    )
    for value in out_quantities_dict.values():
      output_to_file += f"{value}, "
    output_to_file = f"{output_to_file[:-2]});\n}}\n"

    body += output_to_file

    body += r"""
  const int offsetpt_firstfield = diagnosticptoffsetstruct->offset_diagnostic_1d_y_pt[which_pt];
  Ck::IO::write(token, output_to_file, sizeinbytes, offsetpt_firstfield);
}
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


