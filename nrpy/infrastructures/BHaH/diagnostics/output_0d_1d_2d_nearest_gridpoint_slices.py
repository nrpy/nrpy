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

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, Tuple, Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.simple_loop as lp


def register_CFunction_diagnostics_nearest_grid_center(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    filename_tuple: Tuple[str, str] = (
        "out0d-conv_factor%.2f.txt",
        "convergence_factor",
    ),
    pointer_decorator: str = "",
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for "0-dimensional" simulation diagnostics -- output data at gridpoint closest to grid center.

    This function generates a C function that computes and outputs specified diagnostic quantities at the grid point nearest to the physical center of the grid for a given coordinate system. The output is written to a file whose name is specified by `filename_tuple`.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary mapping (type, name) tuples to corresponding C definitions.
        Example: {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}
    :param filename_tuple: Tuple specifying the filename and its corresponding string format.
        Default: ("out0d-conv_factor%.2f.txt", "convergence_factor")
    :param pointer_decorator: Optional dectorates (e.g. [[maybe_unused]])

    :return: None if in registration phase, else the updated NRPy environment.

    :raises ValueError: If an unsupported coordinate system is specified, ensuring that diagnostics are only generated for coordinate systems with a defined grid center.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> Coord = "SinhCylindrical"
    >>> axis = "y"
    >>> _ = register_CFunction_diagnostics_nearest_grid_center(Coord, out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"})
    >>> diag_gc = cfc.CFunction_dict[f"diagnostics_nearest_grid_center__rfm__{Coord}"].full_function
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AbZAvFdABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIwIlIQLEtn60ZvFpp58iRzPO8zL7jj2RkmJNOOz7kLOZ5pcKkv0uEf85zfc7uGHClgRs9RiDAbyNwyZiAZlTl/jf5X1BOhD4i/CNHqvMIASlwDFj3mV52zX1ZUDkoYG52BfOdNFZacCi5oHsDHgM7+ghHg469486fkW6QU/5PkBfqLgwrfD2kD9IC2/fjze4k4cO8FgwxirbR6/LauUHbkFJsFCBAZ83EfkAW5j/GtzV6YUHdBv/xEzPtaol2zX0IFQZZ7cwWXYpqvP1SvRMe/up5H/J9l0ysvNaZ7yv3AENv9B8I5DngpH9vEjcbljHUyG3agNS4mCDUiZ07JBgndzOv9CCMg55ByQV5EjBaIXbxsNVZeZ5BceUiRGM+NlZYB/FUD3wzHj2FI8Tb0C2VcIooTKZ+Ql8sC/0uQVoaOXDqr//0g/oC9jlCeSVCNCr0KbdquKAuMhn5udBlHnrpPPWnsQfO+9MlpS2+wWUwqRTT1DbieX322/WBs8OcqXc3hUWj3Ge5ix8r8DYX2EirrnEQ/+/kyG9DFqPrePpQslc+mnjcRNSO/7HTOGgF+DSywHadHneQFmNLLOKf4ifRUfd+Gdjoi/OeP2NC/MkHXX9J4TF7pEHiNUgYrf3Itr8XfzcYl1W1EtnjnTVqA8svob8wAAfmLOAa7uMW6FSXSLvPdZ8YYBaCVE5Liw3lsKpPGOwKAKry7gWwPPFGDaLShoi98sqFoUK0Txxy1ZbPgNS8IXwknsjq/HPOqBJm/l3qT/gbqcolAx51ND+irqSJ9qFfIBTPMUFHdKoJfhMMHoWNMgAGBZtcIxfphteQeZ0a8YGLCWmva1Wlo0ZOjXHPhK4zm33skfijZeziJFQx9GsCXNQPvkCm/+/EEBmrmU9MV/1XfkcaDaCX+h9XLyyr5e0+dtfpZsJGHi8SgAAAAADl0tNlRAGGwQABjQbaDQAA1ewKr7HEZ/sCAAAAAARZWg==")
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
    cfunc_type = "void"
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
        printf_format = "%.15e" if key[0] != "int" else "%d"
        fprintf += f"{printf_format} "
    fprintf = rf'{fprintf[:-1]}\n", time, '
    for key in out_quantities_dict.keys():
        fprintf += f"{key[1]}, "
    fprintf = f"{fprintf[:-2]});\n"

    pointer_decorator = pointer_decorator if pointer_decorator == "" else f"{pointer_decorator} "

    body = rf"""
// Unpack grid function pointers from gridfuncs struct
{pointer_decorator}const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
{pointer_decorator}const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
{pointer_decorator}const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

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
        cfunc_type=cfunc_type,
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
    pointer_decorator: str = "",
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
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4Ar+BKhdABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIzLaV4Uh8eHMdkFke5ZEWqOrJB8gRVzpSQFQuNk6de9OlmLCKH/0/3hkEO+FpHWMRaKtUzxXmXROgszzTA582kYrIga+g7jqdNJxR8QmVOk3X41Z00Qr8cTShS0e4Ve/Rsfwz7AFvRk7u+g/DHsZZvuk1HMz/N+rflpetPQKKDsLx7wMAnJcaXnErnfc94NuhxGxcaLALyAW85ZpxXG+WHqDC6eXJKplcJ6yp1Ceh1tAXGQuwOP1JDjAbnb4He5ORH1qAjrPi5jbGJIaH2deJL06PMfVo0d8qpj7PrzdCzlupzrMfx9khZaPsE3OJxc2abXt+RClEd9+CSY5CPD7xMp4LoGM8YJxBNcjcCzB+ZJmUlbNi0tl4l3B+wj3BSaec3ts8ni77sbO9NA3OUODUPoVweEeC15HkEJgPuD1buhdoJvJrpX7fT9rRAwD2ZemKVcaf2LEtT/x/W9rvA2sN/dXnhlFxjnmkQwFHJ2QE2UKxzv46Y7bzpk63uag7KaxiKbzcBk6Ze691SMfQMafk3TNVI1zSfWnBtrcxpIQS1AVL5eN5P1BCU/+aPIPuapn41pqDfwuj6fv0/B7QhGgm7Q2bQSIiBrXoL1xCXAE4MachV1wbm5shIvuCiZV5HVAgkLTAFe3DkK2VcqzYYNS55sPUevoChWFMbr8EPBsXXiOl7cr7nCNUr8h7lwlE2fSEr/IZfOMlbmpCDOTRJqTzCcLb8wLOFe2bz6JFPPtDYY7UmNPrGaMdjTdctE5aZ24SeDOKnKucyhEFsADzBnEz5eWv+xg16jWH6t01zDhC+o1dNfiUC8rzmvSNFlLbVo0y/w7Pwgp/eDlOTpRtpHFHQihHwBCnpy3h5s3YUD+c+2o0oKRWVG5xbyZXj9hiChbAP7PV7d+FTFV8fiXT1GiiajndMMTAW/7fk2rV4GG2eWm604agNcKdsrGUtPJ+ZavYKiewIjZjGXLvOX8hcUcgISDaXhKcSEuLOcNPBI4Qsiv0PXH6Lj6WnyeYEdWajcBCczlvSzygS4ypXcTfd6IfELhwN/pl1YlH2mnNu/R9mUW/cH0j5d7VlMwHRr00ix8BgZ88K7wlGg5SqnKQ6mgFLf0YzuM1j6QoXx049omsorweG9v4i31GKvlNZ8a94tkXa+xCb7jBvIyXe3koejEdaE+HSzTNpxku24WV/eEgEPuA6jQWy0AKylTE6Vm0LDnzn5LyjnwKbUzHvT1Mjg3LTu3P/Oba1+Ahxelns3nh308Fn6S7QJki/QxA5Smp9J5ER24msk+4XRCNs4GY5XECTbf22KM26qj5sL74HlQhhrkCdEGdZQlVL5DBB7pR0st+dziPCKoGc3Lbr5crxC+vfIUaFkDC5Nx33wqXZlUJUix6qI1+Z7JXJ9LryYZrxFfup2sNgV9c/vodp3q56bORCD1Rf9mruUdUIGMW4aBTSXfSzmT9firvg5BOZsv6z7Nub0PkJ8dBAXZrENdelj5wj2oPbgiLrSr1EGA0W9qV1j0wkwAALlviG3fOResAAcQJ/xUAAGSN1ICxxGf7AgAAAAAEWVo=")
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
    cfunc_type = "void"
    name = f"diagnostics_nearest_1d_{axis}_axis"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs"
    pointer_decorator = pointer_decorator if pointer_decorator == "" else f"{pointer_decorator} "
    body = rf"""
// Unpack grid function pointers from gridfuncs struct
{pointer_decorator}const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
{pointer_decorator}const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
{pointer_decorator}const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

// Prepare output filename based on 1D axis
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
        cfunc_type=cfunc_type,
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
    pointer_decorator: str = "",
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
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AhtA8JdABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIwIlIQLEtn60ZvFpp58iRzPO8zL7jj2RkmJNOOz7kLOZ5pcKkv92RfBBV6XfvVRqJRIbNiOKE+HSFEwb7csuCUsf9zB0qgLo0SytFxhnW14mQi8WOzQtPby0lBoglttpysw/m7b1NSy56Jyd6v3ORQjM5XKs0j68j0X0IJzMCVR8Dl2FigSsCDy+t+cywlHWLwgbGz8o84kiEsQNi8vwXaFG0ILHJU4KOQHes11bMeMB40iuY3xBqqX6NPdu9FpyfZ66JRRjS85AfYWbHO9hol/CuH+D8Cp9SnpsenErAIgWBfkALisnuCPwn4W6SKOTgP9GrrVAb21KWhviWsb2Yd9el+R7wJOZTJKboiziPVI0fOtkOCmlDdGq3yKY2JcjGQMRZQDeowHkCPrImcerFP1FdFoVZ/Yo0clabNbpjyfpvhAbs4fG0yR7BDmpms5jfTDauLw5ltym/rTTrUK8VBBdawh0N4KhUt1IwsW5/aaqTjVZxCzCEQ7+FIGuJv+z/YJKs5DR3RYD0Mc1ADVFsvJIbemQ4E+mXX3o9NHmihcuOcMk2NsAkyM4W3b1q+lwr57pO4V7Jy8vw0XqFxQdR1CdOZmgrHBPm+vVObPmkQ/qFb9hxZBgF/5IrhXHrq694uEWdaqMSYF1SHEDeqjhFicfIfZ5CFbQfyhgsakceL4/ARNIZn+LbPCRQHPExJhFmqW459u4HoDpYQyGwOkNkUaByKJK6QncNLo/WsqQx/zxVHGCMN/lIUJ9Olh0yLsKxQ60bhJrZK55QJjtjkkeATjYBBY+zieZ1K/IOe7zh0jPNy505LuKmaoQu2aBwNiV9By2gDdvYMtXw/LPIiPGCM9Y4QVKgKjKHW/8sFEUAaYTXJuJtNSexr087UwiRY95KsZzF66dDfPC/wn6eSQBFMcPDLIy0aOZpLdxf5WALw8/htdEi7SUEcWrqabIl8OrQ6/rx2AwS26BykvURplN9327ShuEfNVllO3c3kbPpX02wOWSESzuiC4GZEg8WXypNapMIy22h2h3+fgiULWTjUtu56iqUGMDciv39XjqOTyFT+YPOL/6FrvsiwkAymiHXEerfUsno63HRCMG5sFuJKOWea2uPyfYn1+wZ5icQOc3rwHBNK1p8cqOA7mjynpE8GoQ3ZUiTCpxGzL74MOZBY8AA8o9dQNd9AStL4uLxJZCrVVBbAjyLHlMvh1/gAAAA8b2eOSp3QZoAAd4H7hAAAOngx3OxxGf7AgAAAAAEWVo=")
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
    cfunc_type = "void"
    name = f"diagnostics_nearest_2d_{plane}_plane"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs"
    pointer_decorator = pointer_decorator if pointer_decorator == "" else f"{pointer_decorator} "

    body = rf"""
// Unpack grid function pointers from gridfuncs struct
{pointer_decorator}const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
{pointer_decorator}const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
{pointer_decorator}const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

// Prepare output filename based on 2D plane
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
        cfunc_type=cfunc_type,
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
