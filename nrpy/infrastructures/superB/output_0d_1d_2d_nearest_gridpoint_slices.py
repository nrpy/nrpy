"""
Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from typing import Dict, Tuple, Union, cast
from inspect import currentframe as cfr
from types import FrameType as FT
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


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
    desc = f"Output diagnostic quantities at gridpoints closest to {axis} axis."
    c_type = "void"
    name = f"diagnostics_nearest_1d_{axis}_axis"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs, const diagnostic_pt_offset_struct *restrict diagnosticptoffsetstruct, Ck::IO::Session token"

    body = rf"""
// Unpack gridfuncs struct:
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

const int num_diagnostic_pts = {'diagnosticptoffsetstruct->num_diagnostic_1d_y_pts;' if axis == "y" else 'diagnosticptoffsetstruct->num_diagnostic_1d_z_pts;'}
const int *restrict idx3_diagnostic_pt = {'diagnosticptoffsetstruct->localidx3_diagnostic_1d_y_pt;' if axis == "y" else 'diagnosticptoffsetstruct->localidx3_diagnostic_1d_z_pt;'}
const int *restrict i0_diagnostic_pt = {'diagnosticptoffsetstruct->locali0_diagnostic_1d_y_pt;' if axis == "y" else 'diagnosticptoffsetstruct->locali0_diagnostic_1d_z_pt;'}
const int *restrict i1_diagnostic_pt = {'diagnosticptoffsetstruct->locali1_diagnostic_1d_y_pt;' if axis == "y" else 'diagnosticptoffsetstruct->locali1_diagnostic_1d_z_pt;'}
const int *restrict i2_diagnostic_pt = {'diagnosticptoffsetstruct->locali2_diagnostic_1d_y_pt;' if axis == "y" else 'diagnosticptoffsetstruct->locali2_diagnostic_1d_z_pt;'}
const int *restrict offsetpt_firstfield = {'diagnosticptoffsetstruct->offset_diagnostic_1d_y_pt;' if axis == "y" else 'diagnosticptoffsetstruct->offset_diagnostic_1d_z_pt;'}
for (int which_pt = 0; which_pt < num_diagnostic_pts; which_pt++) {{
  const int idx3 = idx3_diagnostic_pt[which_pt];
  const int i0 = i0_diagnostic_pt[which_pt];
  const int i1 = i1_diagnostic_pt[which_pt];
  const int i2 = i2_diagnostic_ptt[which_pt];
  REAL xCart[3];
  xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
  const REAL xCart_axis = {'xCart[1];' if axis == "y" else 'xCart[2];'}
  int sizeinbytes = 23 * (len(out_quantities_dict) + 1);
  char out[sizeinbytes+1];
  snprintf(out, sizeof(out), " """
    output_to_file = "% .15e,"
    for key in out_quantities_dict.keys():
      printf_c_type = "% .15e" if key[0] != "int" else "%d"
      output_to_file += f"{printf_c_type} "

    output_to_file = (
      f'{output_to_file[:-1]}\\n", xCart_axis, '
    )
    for value in out_quantities_dict.values():
      output_to_file += f"{value}, "
    output_to_file = f"{output_to_file[:-2]});\n"

    body += output_to_file

    body += rf"""    
  Ck::IO::write(token, output_to_file, sizeinbytes, offsetpt_firstfield[which_pt]);
}}
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

def register_CFunction_diagnostics_nearest_2d_plane(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    plane: str,
    filename_tuple: Tuple[str, str] = (
        "out2d-PLANE-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for 2-dimensional simulation diagnostics at gridpoints closest to the specified plane.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary of output quantities.
    :param plane: Specifies the plane ("xy", "yz") for the diagnostics.
    :param filename_tuple: Tuple containing the format for filename and the replacement arguments.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If the specified plane is not supported.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> Coord = "SinhSymTP"
    >>> plane = "yz"
    >>> _ = register_CFunction_diagnostics_nearest_2d_plane("SinhSymTP", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane=plane)
    >>> diag2d = cfc.CFunction_dict[f"diagnostics_nearest_2d_{plane}_plane__rfm__{Coord}"].full_function
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AfSA3ldABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIwIlIQLEtn60ZvFpp58iRzPO8zL7jj2RkmJNOOz7kLOZ5pcKkv92RfBBV6XfvVRqJRIbNiOKE+HSFEwb7csuCUsf9zB0qgLo0SytFxhnW14mQi8WOzQtPby0lBoglttpysw/m7b1NSy56Jyd6v3ORQjM5XKs0j68j0X0IJzMCVR8Dl2FigSsCDy+t+cywlHWLwgbGz8o84kiEsQNi8vwXaFG0ILHJU4KOQHes11bMeMB40iuY3xBqqX6NPdu9FpyfZ66JRRjS85GSxBso4cKi5eGY674iGPLQBQJh/i3W3LgupqIV5Bx1X3o/450Xgj5NjMue4i7H/UfL1qp71Ge6ncAzZ/oKmsPb0h3jrM4DKNwTR39l94STUk3qNpztEXcij8XS4r5W26aH14PylTqhnT83m4d0lNIvXNKIg5tGVQnE6fEAdqPwak1KP1vcYvE9AH1jsWRuSglAlkgAO6zsNGNQvnWq0+BzWA3ywSu8sQ9UBWM9SM6tAc5AazjarASStyn6hOB+ezWZUqpZHcpzvBB7/IEl6/AZzMrbLtyodvAD4I79TJJj9zk7Zn+ZqfVrvfiMNupWLX9QigMAAX3uhZ8Q6+IjCMMIft5IsBCcVTLJKtPRhO/swcHW0FfsBU1Fthx5uRFAOCYFt1OM8Ng9w/f+9hPmolVQj9LnGbVtjafWZ0Ocbg3jxZAbOrRucb9ICXIOX63o05BYQZXCyvt6dqyWvgiBJlsDnAz7DjXCOdj+fWEHWu74+A5TbD+bn0ubpK0b7rwqJ+pbibLGbS6AF8KfE6cNGW108E78JdbM3D9wAv167SThtHXFcD2n3EaDFmnrhiQRpbWhJ1UbHrmOLPvDkrQR7YVP2CTEpWv9ZJ0WHhZheCvjkPAQehJyXxnuWPZ4ezEGgJuULcekUg8iv2i5K23Xm+00mf8YZqur7R+ymDFn7qI4yUTjk8Xo5UTV5BwYATg9vR/V05cXm3JyJy59eZLk3b/VGTn2DKIT90IqXwDV9+lkov+HAZt+hebaINtWfu733vfy1+wn2dj/fZmrrNhmkEIV4o8vFj/pDXb1IqsTTh5lQxFDRTw35DPalplJS79vz4B9rAmS1ghSR+iV+EAAAAAAn1cKhAngTtAAAZUH0w8AAHfeYwOxxGf7AgAAAAAEWVo=")
    >>> if diag2d != expected_string:
    ...     error_message = diff_strings(expected_string, diag2d)
    ...     raise ValueError(f"\n{error_message}\n base64-encoded output: {compress_string_to_base64(diag2d)}")
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    if plane not in ["xy", "yz"]:
        raise ValueError(
            f"Output along {plane} plane not supported. Please choose xy or yz plane."
        )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"Output diagnostic quantities at gridpoints closest to {plane} plane."
    c_type = "void"
    name = f"diagnostics_nearest_2d_{plane}_plane"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs, const diagnostic_pt_offset_struct *restrict diagnosticptoffsetstruct, Ck::IO::Session token"

    body = rf"""
// Unpack gridfuncs struct:
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;
// Unpack diagnosticptoffset struct:
const int num_diagnostic_pts = {'diagnosticptoffsetstruct->num_diagnostic_2d_xy_pts;' if plane == "xy" else 'diagnosticptoffsetstruct->num_diagnostic_2d_yz_pts;'}
const int *restrict idx3_diagnostic_pt = {'diagnosticptoffsetstruct->localidx3_diagnostic_2d_xy_pt;' if plane == "xy" else 'diagnosticptoffsetstruct->localidx3_diagnostic_2d_yz_pt;'}
const int *restrict i0_diagnostic_pt = {'diagnosticptoffsetstruct->locali0_diagnostic_2d_xy_pt;' if plane == "xy" else 'diagnosticptoffsetstruct->locali0_diagnostic_2d_yz_pt;'}
const int *restrict i1_diagnostic_pt = {'diagnosticptoffsetstruct->locali1_diagnostic_2d_xy_pt;' if plane == "xy" else 'diagnosticptoffsetstruct->locali1_diagnostic_2d_yz_pt;'}
const int *restrict i2_diagnostic_pt = {'diagnosticptoffsetstruct->locali2_diagnostic_2d_xy_pt;' if plane == "xy" else 'diagnosticptoffsetstruct->locali2_diagnostic_2d_yz_pt;'}
const int *restrict offsetpt_firstfield = {'diagnosticptoffsetstruct->offset_diagnostic_2d_xy_pt;' if plane == "xy" else 'diagnosticptoffsetstruct->offset_diagnostic_2d_yz_pt;'}
for (int which_pt = 0; which_pt < num_diagnostic_pts; which_pt++) {{
  const int idx3 = idx3_diagnostic_pt[which_pt];
  const int i0 = i0_diagnostic_pt[which_pt];
  const int i1 = i1_diagnostic_pt[which_pt];
  const int i2 = i2_diagnostic_ptt[which_pt];
  REAL xCart[3];
  xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
  int sizeinbytes = 23 * (len(out_quantities_dict) + 2);
  char out[sizeinbytes+1];
  snprintf(out, sizeof(out), " """
    output_to_file = "% .15e % .15e,"
    for key in out_quantities_dict.keys():
      printf_c_type = "% .15e" if key[0] != "int" else "%d"
      output_to_file += f"{printf_c_type} "

    if plane == "xy":
      output_to_file = (
        f'{output_to_file[:-1]}\\n", xCart[0], xCart[1], '
      )
    elif plane == "yz":
      output_to_file = (
        f'{output_to_file[:-1]}\\n", xCart[1], xCart[2], '
      )
    for value in out_quantities_dict.values():
      output_to_file += f"{value}, "
    output_to_file = f"{output_to_file[:-2]});\n"

    body += output_to_file

    body += rf"""
  const int offsetpt_firstfield = {'diagnosticptoffsetstruct->offset_diagnostic_2d_xy_pt[which_pt];' if plane == "xy" else 'diagnosticptoffsetstruct->offset_diagnostic_2d_yz_pt[which_pt];'}
  Ck::IO::write(token, output_to_file, sizeinbytes, offsetpt_firstfield[which_pt]);
}}
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

