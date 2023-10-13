"""
Output simulation data in 0D (grid physical center), 1D (axis), and 2D (plane) slices.

Functions:
----------
- register_CFunction_diagnostics_grid_center: Registers C functions for
    0-dimensional simulation diagnostics, focusing on data at the grid center.

- register_CFunction_diagnostics_1d_axis: Registers C functions for
    1-dimensional simulation diagnostics along a specified axis such as "x" or "z".

- register_CFunction_diagnostics_2d_plane: Registers C functions for
    2-dimensional simulation diagnostics on a specified plane like "xy" or "yz".

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


def register_CFunction_diagnostics_grid_center(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    filename_tuple: Tuple[str, str] = (
        "out0d-conv_factor%.2f.txt",
        "convergence_factor",
    ),
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for "0-dimensional" simulation diagnostics -- output data at the grid center.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary mapping (type, name) tuples to corresponding C definitions.
       Example: {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}
    :param filename_tuple: Tuple specifying the filename and its corresponding string format.
       Default: ("out0d-conv_factor%.2f.txt", "convergence_factor")
    :return: None if in registration phase, else the updated NRPy environment.

    The function will generate a C function that unpacks grid functions, defines a filename based on the provided format,
    and then computes and prints the specified diagnostic quantities to the file for the specific coordinate system provided.

    Doctests:
    >>> Coord = "SinhCylindrical"
    >>> _ = register_CFunction_diagnostics_grid_center(Coord, out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"})
    >>> print(cfc.CFunction_dict[f"diagnostics_grid_center__rfm__{Coord}"].full_function)
    #include "../BHaH_defines.h"
    #include "../BHaH_function_prototypes.h"
    /*
     * Output diagnostic quantities at grid's *physical* center.
     * For example:
     * In Cartesian this will be at i0_mid,i1_mid,i2_mid.
     * In Spherical, this will be at i0_min,i1_mid,i2_mid (i1 and i2 don't matter).
     * In Cylindrical, this will be at i0_min,i1_mid,i2_mid (i1 == phi doesn't matter).
     * In SinhSymTP, this will be at i0_min,i1_mid,i2_mid (i2 == phi doesn't matter).
     */
    void diagnostics_grid_center__rfm__SinhCylindrical(commondata_struct *restrict commondata, params_struct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs) {
    #include "../set_CodeParameters.h"
    <BLANKLINE>
      // Unpack gridfuncs struct:
      const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
      const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
      const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;
    <BLANKLINE>
      // Output to file diagnostic quantities at grid's *physical* center.
      char filename[256];
      sprintf(filename, "out0d-conv_factor%.2f.txt", convergence_factor);
      FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
      if (!outfile) {
        fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
        exit(1);
      }
    <BLANKLINE>
      const int i0_center = NGHOSTS;
      const int i1_center = Nxx_plus_2NGHOSTS1 / 2;
      const int i2_center = Nxx_plus_2NGHOSTS2 / 2;
      if (i0_center != -1 && i1_center != -1 && i2_center != -1) {
        const int idx3 = IDX3(i0_center, i1_center, i2_center);
        const REAL log10HL = log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16));
    <BLANKLINE>
        fprintf(outfile, "%e %.15e\n", time, log10HL);
      }
      fclose(outfile);
    }
    <BLANKLINE>
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
    name = "diagnostics_grid_center"
    params = "commondata_struct *restrict commondata, params_struct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"

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


def register_CFunction_diagnostics_1d_axis(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    axis: str,
    filename_tuple: Tuple[str, str] = (
        "out1d-AXIS-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for 1-dimensional simulation diagnostics along a specified axis.

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
    >>> _ = register_CFunction_diagnostics_1d_axis(Coord, out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis=axis)
    >>> diag1d = cfc.CFunction_dict[f"diagnostics_1d_{axis}_axis__rfm__{Coord}"].full_function
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AnoBEZdABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnSEgrtsCyPSYXzh6+ATGh7ABAkIBSF+o891I4eI9S4eGHzAvGSs7EJpzVyrmy021AfTIND2mYmaWYvVjEpYWAIWJHt6abEm12uy7LfoZTeEyiT+V7goyEruSX8ZbKfF7laQ8cZFJJGqW0hlaZ9xL1O2wbDhTDbIud4y41c+Frfqux0cJRkA9+ftib6biz4gC7XhKoBCKQDExonclhZNXjDOQ4WACA16U3zUNl7b/eEw9iimCglG79US4A+Evk2fi0APwIQ+oI7szdJNmthBLtTAs14+nlh4/DvWDefIt65GBBZC9OcGTxI9cdoCsBRm8wDgGABDWaUqjRZGNZ13gqA8OQN1e5HQdYE82Fhic08j4ZNni9746oPietgEetXIrgoIcd7LnwGr9N23gyGiXMrKBEC+0HbrPw1XmDdYD/ZGCXW/+h9L8+wf+ySU5fSQdWm5gFlbs1bZ7dCocqrYLnLxLUiZDuaMz9Q2VsitMbKGnDBjWXrt/lHLKdv15LQrNk/8d5lVz5xwDCKTAOdvo63lFUfJxeT7Ygxw0zLSEBvSAQuMz3bosH3BtbZ9wmkHNzby96pfnK8D0rPD4lg5TT+mqp1cEDr6GeUSrDtd7mTzWKKPfLW8stYxflfwh5UzX/qZ0hd/3dygbBDF4+zkNH0aVsHukqWhr9GwlwDBYxpSLdneuw9nOLRqLk90uM0PMteT7wE/2Xc8FEtu9QEY3qybO+sp3+D6Bc14nSu8n9m/8aU+/03prnGOxUqad8MM5tsdLfbgNWQbCSf6BNeh1E0/RtR9UUn9xbigmH6Hyyqiv1RCCekzgpzM6NShSHEeInaxPRFvzA39V4rb+MFw+XY8msEyytcp4fhinE5rkyfbUNzlvTv0xcq9gkRCuGvO/Bgn5WcwiLlSdG39d3Mfau1eHtbvNMt+L71JT+BdtyM/qDinKQVjGn2rAsPoRbXImJsCJ/EqQE6ZuYeCc6WnxHHyM+vNVD//z2zsGMXVNUf4HxFo3l5om+S65bp4wR4FuHPya8owdSFqnPEfBu4yshkX7+u53zlMAbo8wS7mj9yh4kgXHTkTZQlaKlPQ1zm5SdbhS+nepDD1Qe+p/FuwlvWbeQVdK/Y9iDVZO1MnROKM5/iLmWimACn5trCJkSD0tCEit4/DRLwQzaL7OoPims6iip/lgSW+ah+ja35UfnhWnO9oshwGorVsUPmpnqSklcigkPdkzHOsF7GgZHe54nEwpYx05MXwFC3mJ6LzdBDNYzGE2A8eMVZho7TXCSDpLagJTDgtHz7tnPoWEaC+2LM2faBWXqZPigsMfVAjF3JwZItzte+3V/oICMo3EbF1M8Feds44BnPOb+liy8o/ryNe9PyVPeD8rlS2KFBJG3q2lS9gAAAApkPCfw/UFk0AAeII6RMAAA42vh2xxGf7AgAAAAAEWVo=")
    >>> if diag1d != expected_string:
    ...     raise ValueError(f"{diff_strings(expected_string, diag1d)}\n base64-encoded output: {compress_string_to_base64(diag1d)}")
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
    name = f"diagnostics_1d_{axis}_axis"
    params = "commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs"

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


def register_CFunction_diagnostics_2d_plane(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    plane: str,
    filename_tuple: Tuple[str, str] = (
        "out2d-PLANE-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for 2-dimensional simulation diagnostics in a specified plane.

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
    >>> _ = register_CFunction_diagnostics_2d_plane("SinhSymTP", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane=plane)
    >>> diag2d = cfc.CFunction_dict[f"diagnostics_2d_{plane}_plane__rfm__{Coord}"].full_function
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AePA25dABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIwIlIQLEtn60ZvFpp58iRzPO8zL7jj2RkmJNOOz7kLOZ5pcKkv92RfBBV6XfvVRqJRIbNiOKE+HSFEwb7csuCUsf9zB0qgLjaNojvcXVzhXj5lrq9G7GuZvGWxSqMoQWYxmpol2GdvZjHEN5O4OgyeP4QzVSq9/Zu9h9Y0ty+TEFocghu8pf+y9YQkhJOYT5A4IjkO73sC/jYjkv0a/oWFPfH2VsSbaeM1EhZ50gJnJak+jJrUGuc96ZXOOKxZ1MoAeXOuZ4KJlXZJ4FYaB6vtsnUwqXUHh/yxI2Z1P7s5VdHCKPuNKOoZ8gOd9/BwxhwRLSR0mTe0guvTKlhQHzWoglERwDuqt+Yw0TtQ8YuGxxgcgFiSkVuZX2+xFqmY6D+LgSeTHZtfJr4DFW1xyeAw7xiafWFq8I2ROxOagZi+9VYaACh7ISl2iKOe2EaH2GOug8Ak9onC0w9mE5cky2JJoh6wUyLKAQUqD87fwYp/HHKfOv5SCh6+HF7Q8KqW21pywCsqIV8TdMTECcn81w7qmuPKn2nIhVOb3h0kALfaXD39UQ46z3o79PQuB4NLtmAAH0KvKxSC/PpMiz0FKJWsRv3/k1eczqEzFwWjdz4zGa3laQrGZytyLUprurs4vv3YTzGpWZ+XSSvFEppbK1A3A4/ljeCaKtGWn8dfjgIMS1VvgPqzkyWDt3jb1tJh4lv+WcEj8mssEEi5TXUy1jwVgP6JCgHxezu8dJ8ysJQbDH9q57Zazkzxcxf+60W03tuzZRrDlaqE40ooXuyIN+5YpGp+Xom6/0loxModsFdqYIHIgK1IjcXBpVjA3aVLFV/xEyBTKko8Ou5fSeI5YWr+5Af86znjCLUCzjEaS4wU+IVGqhoQtZ4dySCk4qVM5gL8QGESoeE6S+uyIqlvkrcOnwH7x0+2AIHugPhdUoqOLVbKsytDOc9YvCgOOaPrKik3YZewd7cdSMX5wt7DacCng7v4vhgwR3NAmwAL1Erw0ZTo0JxaZmBFBcZSTOLgBloNUGeC4AELXas80mre7CTWeWwd0l9h0hjlpRwrzUvw1L9kolTM5o2vIcxJFARVn8J0moTkqJOzgAAAAAA0PTK7WpKxcAAAYoHkA8AAOqZQnixxGf7AgAAAAAEWVo=")
    >>> if diag2d != expected_string:
    ...     error_message = diff_strings(expected_string, diag2d)
    ...     raise ValueError(f"{error_message}\n base64-encoded output: {compress_string_to_base64(diag2d)}")
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
    name = f"diagnostics_2d_{plane}_plane"
    params = "commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs"

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
