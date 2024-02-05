"""
Output simulation data at points closest to grid physical center (0D), grid y or z axis (1D), and grid xy or yz plane (2D).

Functions:
----------
- register_CFunction_diagnostics_nearest_grid_center: Registers C functions for
    0-dimensional simulation diagnostics, focusing on gridpoint closest to the grid center.

- register_CFunction_diagnostics_nearest_1d_axis: Registers C functions for
    1-dimensional simulation diagnostics at gridpoints closest to a specified axis such as "x" or "z".

- register_CFunction_diagnostics_nearest_2d_plane: Registers C functions for
    2-dimensional simulation diagnostics at gridpoints closest to  a specified plane like "xy" or "yz".

Each function in this module prepares C functions for writing diagnostic
quantities of simulations to respective files based on the provided coordinate
system, output quantities dictionary, and other parameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, Tuple, Union, cast
from inspect import currentframe as cfr
from types import FrameType as FT

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp


def register_CFunction_diagnostics_nearest_grid_center(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    filename_tuple: Tuple[str, str] = (
        "out0d-conv_factor%.2f.txt",
        "convergence_factor",
    ),
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for "0-dimensional" simulation diagnostics -- output data at gridpoint closest to grid center.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary mapping (type, name) tuples to corresponding C definitions.
       Example: {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}
    :param filename_tuple: Tuple specifying the filename and its corresponding string format.
       Default: ("out0d-conv_factor%.2f.txt", "convergence_factor")
    :return: None if in registration phase, else the updated NRPy environment.

    The function will generate a C function that unpacks grid functions, defines a filename based on the provided format,
    and then computes and prints the specified diagnostic quantities to the file for the specific coordinate system provided.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> Coord = "SinhCylindrical"
    >>> axis = "y"
    >>> _ = register_CFunction_diagnostics_nearest_grid_center(Coord, out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"})
    >>> diag_gc = cfc.CFunction_dict[f"diagnostics_nearest_grid_center__rfm__{Coord}"].full_function
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4Aa+AuRdABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIwIlIQLEtn60ZvFpp58iRzPO8zL7jj2RkmJNOOz7kLOZ5pcKkv0uEf85zfc7uGHClgRs9RiDAbyNwyZiAZlTl/jf5X1BOhD4i/CNHqvMIASlwDFj3mV52zX1ZUDkoYG52BfOdNFZacCi5oHsDHgM7+ghHg469486fkW6QU/5PkBfqLgwrfD2kD9IC2/fjze4k4cO8FgwxirbR6/LauUHbkFJsFCBAZ83EfkAW5j/GtzV6YUHdBv/xEzPtaol2zX0IFQZZ7cwWXYpqvP1SvRMe/up5H/J9l0ysvNaZ7yv3AENv9B8I5DngpH9vEjcbljHUyG3agNS4mCDUiZ07JBgndzOv9CCMg55ByQV5EjBaIXbxsNVZeZ5BceUiRGM+NlZYB/FUEfiEUeOfpMCu6CAEr+oelS8jeJZ201U24vYd7Lw+3SHXmHlS9JJiXbWoTYNK20XfM69UZQzvbXqqTHOpGiSDDbBRs/bNc+aDPCbWNhqGxvY1FE8iTCXMFlwYbxpYTRdObzQ490084YvkKIUvwc5wZ0yyydPoLfrTdRmCptnib3BjSvZhDfq1lg3SJ8RBiuucho+lHgKIHk3BQEM4mUjmrbasqcQIvCHHv0azE3u8VA5J4StrwgrwV0MXwgLzePSgKyTDl6g36mnSTgXrkHJwkrzvhs6mtGayfw+J8WsZFqm4QlLgMh0tdyUatAFMb5MGOduOGQrL0oD8FWwnBjdb1ElrWL4XBheNoh5uvv/GYYv4GoW6rIZcyEz1WpXvhyMG9lk7uKBJpnNcoaUb9w5j6cyhLdtMERjJIDfxiemqqYn55zw0oxY69wxeTQAN3LCbDPXOx1URfkcuXxqQXNNXEOPkkDLO62vcRulZp0xWhRKYWy8mEwql4t94UXahNPSuyuyXrkne0Pef1AoAAJLfEC1BJxBOAAGABr8NAAA6+GwfscRn+wIAAAAABFla")
    >>> if diag_gc != expected_string:
    ...     raise ValueError(f"\n{diff_strings(expected_string, diag_gc)}\n base64-encoded output: {compress_string_to_base64(diag_gc)}")
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Output diagnostic quantities at grid's *physical* center.
For example:
In Cartesian this will be at i0_mid,i1_mid,i2_mid.
In Spherical, this will be at i0_min,i1_mid,i2_mid (i1 and i2 don't matter).
In Cylindrical, this will be at i0_min,i1_mid,i2_mid (i1 == phi doesn't matter).
In SinhSymTP, this will be at i0_min,i1_mid,i2_mid (i2 == phi doesn't matter)."""
    c_type = "void"
    name = "diagnostics_nearest_grid_center"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"

    NGHOSTS = sp.Symbol("NGHOSTS", real=True)
    Nxx_plus_2NGHOSTS = ixp.declarerank1("Nxx_plus_2NGHOSTS")

    # Default to i0_center, i1_center, i2_center to be i0_mid, i1_mid, i2_mid.
    i0_center = Nxx_plus_2NGHOSTS[0] / 2
    i1_center = Nxx_plus_2NGHOSTS[1] / 2
    i2_center = Nxx_plus_2NGHOSTS[2] / 2
    if (
        "Spherical" in CoordSystem
        or "Cylindrical" in CoordSystem
        or "SymTP" in CoordSystem
    ):
        i0_center = NGHOSTS
    elif "Cartesian" in CoordSystem:
        pass  # defaults are correct here.
    else:
        raise ValueError()

    declare_output_variables = ""
    for key, value in out_quantities_dict.items():
        declare_output_variables += f"const {key[0]} {key[1]} = {value};\n"

    fprintf = 'fprintf(outfile, "%e '
    for key in out_quantities_dict.keys():
        printf_c_type = "%.15e" if key[0] != "int" else "%d"
        fprintf += f"{printf_c_type} "
    fprintf = rf'{fprintf[:-1]}\n", time, '
    for key in out_quantities_dict.keys():
        fprintf += f"{key[1]}, "
    fprintf = f"{fprintf[:-2]});\n"

    body = rf"""
// Unpack gridfuncs struct:
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

// Output to file diagnostic quantities at grid's *physical* center.
char filename[256];
sprintf(filename, "{filename_tuple[0]}", {filename_tuple[1]});
FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
if (!outfile) {{
    fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
    exit(1);
}}

const int i0_center = {i0_center};
const int i1_center = {i1_center};
const int i2_center = {i2_center};
if(i0_center != -1 && i1_center != -1 && i2_center != -1) {{
  const int idx3 = IDX3(i0_center, i1_center, i2_center);
  {declare_output_variables}
  {fprintf}
}}
fclose(outfile);
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


def register_CFunction_diagnostics_nearest_1d_axis(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    axis: str,
    filename_tuple: Tuple[str, str] = (
        "out1d-AXIS-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for 1-dimensional simulation diagnostics at gridpoints closest to specified axis.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary of output quantities.
    :param axis: Specifies the axis ("x", "z") for the diagnostics.
    :param filename_tuple: Tuple containing the format for filename and the replacement arguments.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If the specified axis is not supported.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> Coord = "SinhSpherical"
    >>> axis = "y"
    >>> _ = register_CFunction_diagnostics_nearest_1d_axis(Coord, out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis=axis)
    >>> diag1d = cfc.CFunction_dict[f"diagnostics_nearest_1d_{axis}_axis__rfm__{Coord}"].full_function
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AotBFFdABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnSEgrtsCyPSYXzh6+ATGh7ABAkIBSF+o891I4eI9S4eGHzAvGSs7EJpzVyrmy021AfTIND2mYmaWYvVjEpYWAIWJHt6abEm12uy7LfoZTeEyiT+V7goyEruSX8ZbKfF7laQ8cZFJJGqW0hlaZ9xL1O2wbDhTDbIud4y41c+Frfqux0cJRkA9+ftib6biz4gC7XhKoBCKQDExonclhZNXjDOQ4WACA16U3zUNl7b/eEw9iimCglG79US4A+Evk2fi0APwIQ+oI7szdJNmthBLtTAs14+nlh4/DvWDefIt65GBBZC9Ob+ww2x9COzu7LqOqYSPXC+hhHFm10ZnTvBFLP7pSdh/9hi0QkxhP5GLLWRcoStm76DCEpolzaocROhXbaiF79ieeUGnKlzglpjUCal9XpefO4JnhNX54fJpxkUY/Yoowh2j3mUeMgVSFeefOW8+vBF9UISgWRR+92UO7Ylyz2U4mlSV2k7mNbAjy7KM9nbcm4uuZGIdkxV8oVTL0qQDxZ8UhSoZptXdVi28b6InehoNAScgBea46D0l9B6IMy214XJ6TFxrIH5aC0F7+8cVE7i0nAwwA24jrnChdHtkeiI4tJwzOEmURP1YUxlI/YTyeAmeToJXE9M171XKtOOphkvkVvpeEv8SQFNQ0jniSozgFFQeUeNf7HiBDacmYmu24nrTVwaH/T1ww+OOc+DfMPd2Hed/uG8ko2n7jKgNOikB/LLirMezqn6rosZVifm5GZrYbfxoUdewqXijEHFJnnVllqGoZhw4ZMbGNZtjbmBJjjogghtfHYOHrwpIhJOtuCdX6Ynr41Da58Yozyo1wVRawr6/5TaOVgmDeIAb40oKbnog0T7sdj/ZQYw9HJvnsLrU8/r1nxPcIEurm2eSSeLpYXlt9o96OYvSAzoQ0TR2GtM1gER1SEzwyhcmTmVx1sq8s0e5KrJOCVqXPvkuoO88MdMNdQcCS3R8EpHLgHgag3wCx9/7iPVSz9jc/EyJms++8Hf6YZedPjjhGlnfbQoacj2i3O3EwNKbbTJZW6JCM0efCz0nI2lGrVvZFi58mbzgv98CGTHl+BFG1FPMcZK6dqwoUS17Q5jxn+uE4fPOIco2LXg1IaUDDw/h8EvVjWEgcuGQDAHXkjJalD7Xfa4J9UrmG9x/iN8fESoi7ic2VKr9PR63zkrrXINGY/Lv73/ekcadl9YR+UIZ4kRGmRIhpu+OQgL0mJKgiLNS49yG8OfnHwOeWXo8QPRqNlyX7RdK/mUUpYPIdmw0v10r7NYwa9+Kziyei1dRBW8/qAQmBz2xxcL4UbUpkpiFB9HeKEN1QEc8saDLh4qyAJqnON18+wNbTJDqd9JeOp39aHsf0IyYJIuHgoL5QYkZHAgxrAJl6ogGrCp1egAAAAAmcwtyzpCcZsAAe0IrhQAANryZO+xxGf7AgAAAAAEWVo=")
    >>> if diag1d != expected_string:
    ...     raise ValueError(f"\n{diff_strings(expected_string, diag1d)}\n base64-encoded output: {compress_string_to_base64(diag1d)}")
    """
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
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs"

    body = rf"""
// Unpack gridfuncs struct:
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

// 1D output
char filename[256];
sprintf(filename, "{filename_tuple[0].replace('AXIS', axis)}", {filename_tuple[1]});
FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
if (!outfile) {{
    fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
    exit(1);
}}

{loop_1d}

fclose(outfile);
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
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

    loop_2d = lp.simple_loop_2D(
        CoordSystem=CoordSystem,
        out_quantities_dict=out_quantities_dict,
        plane=plane,
    )

    desc = f"Output diagnostic quantities at gridpoints closest to {plane} plane."
    c_type = "void"
    name = f"diagnostics_nearest_2d_{plane}_plane"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs"

    body = rf"""
// Unpack gridfuncs struct:
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

// 1D output
char filename[256];
sprintf(filename, "{filename_tuple[0].replace('PLANE', plane)}", {filename_tuple[1]});
FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
if (!outfile) {{
    fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
    exit(1);
}}

{loop_2d}

fclose(outfile);
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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
