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
    >>> from nrpy.helpers.generic import clang_format, compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> Coord = "SinhCylindrical"
    >>> axis = "y"
    >>> _ = register_CFunction_diagnostics_grid_center(Coord, out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"})
    >>> diag_gc = cfc.CFunction_dict[f"diagnostics_grid_center__rfm__{Coord}"].full_function
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AauAt9dABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIwIlIQLEtn60ZvFpp58iRzPO8zL7jj2RkmJNOOz7kLOZ5pcKkv0uEf85zfc7uGHClgRs9RiDAbyNwyZiAZlTl/jf5X1BOhD4i/CNHqvMIASlwDFj3mV52zX1ZUDkoYG52BfOdNFZacCi5oHsDHgM7+ghHg469486fkW6QU/5PkBfqLgwrfD2kD9IC2/fjze4k4cO8FgwxirbR6/LauUHbkFJsFCBAZ83EfkAW5j/GtzV6YUHdBv/xFH+frbxblZhQaPN6jo0+xmYnNNuMjEHQcKNBI0LpoumpFdjbOyUkc6hPyxdpuLQ8deT47gKgcrOr39mTGoBZMsN1oprbraRoHKIhRjjTK+arTvhlXp05RyD4+Lu5SVpYi3LEDKSxgCNdciAMw8SwVKFW9BZzd4YLJoAyaSEsGNVqUF1XuX2CTosnzO/mLvRu+wadiCfzvwcxhtWESATEBDZQhtu3WAn7u+CkC6YyJVNQ7HpPCtimChmMBqZDpztUVeLP248cuCAPPiXEEc3wQ/r5iW6NZbjG0Q+zAirX5mWG39zdjbqkBRh/FuMuZ+YS7sfQKZcHkfKjgcnFNboARLLXxwxdzDk8An4GNqlmRw8slXCCbd0MJE0gwK50o7TRtQ7rlEvnFz7fxJODEi8pRkEfsFruyke0PoLsexRRxHwN8WziNDagXJUrsTCluxVIHgctWsBsVAy9HmBPrr7tJSmluXnApgk8do8nbeabPj5cqm6dIiYjeNonb6KFkFa5EpFDYxbV8ppZ6/RKvy4dy4yhbDcjMKRKtKgst0m0cKeRRa5034zuHz2YvDcTf53RCaUhSBTEXAR6kFuiBC52LcXu3VpzBKm+AEdC1IwYTHuWgNF77ulbBnk4w/O/78q3G4KfBoxSgBSm/Bef60QDVFD2gAAAnWHu/ME4D8sAAfsFrw0AAHcoMGmxxGf7AgAAAAAEWVo=")
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
    name = "diagnostics_grid_center"
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
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AodBExdABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnSEgrtsCyPSYXzh6+ATGh7ABAkIBSF+o891I4eI9S4eGHzAvGSs7EJpzVyrmy021AfTIND2mYmaWYvVjEpYWAIWJHt6abEm12uy7LfoZTeEyiT+V7goyEruSX8ZbKfF7laQ8cZFJJGqW0hlaZ9xL1O2wbDhTDbIud4y41c+Frfqux0cJRkA9+ftib6biz4gC7XhKoBCKQDExonclhZNXjDOQ4WACA16U3zUNl7b/eEw9iimCglG79US4A+Evk2fi0APwIQ+oI7szdJNmthBLtTAs14+nlh4/DvWDefIt65GBBZC9OcGTxI9cdoCsBRm8wDgGABDWaUqjRZGNZ13gqA8OQN1e5HQdYE82EEsbVq/iiEZWsl1LAWAKJ5lVmhDqHFGwVVUuIG5Y2ioMtS66mCj5mFq4Rx0vgi93SERfg62JpcMFGmPHNQjBf6hHuUn3cC+1opMfRJNQHFZ8khM+QiEwiwxxbduZXOUw+rY2iClh/sx+RltG317V8j/X67sJ9jCsqDbEAEEyqf1SXW8aM3zp73HGNSJQaF9L9PpU4qprjEs5AVvaimJ6gCPt7e6cqv0zl47MEFqM8GOQ6CKmDFNWyxIU4pj1i8338eH934RhTtLSGun4Zs0EI8WV9hMYsBwVFh3AUJnz7QcLM1scWu68G9Erf8WZJPuyYlN3wlSEqwoF3mX2qJROLhzLpmsd8/YFgXjqzvFcXJubyqM1VHe/NIkZNj+pBQ0ly2B1dlm8NU6f6/rqP6KPdF39wbVADiBij0YjJxpV/AmFuZA+btm2cXkRt5CsqIzIvSocS2eLfFNN0m08TcDM3L7Os/qzfc8ym92Q1ZBjpQ4BRiV2k79/gm/hr04fFWK6XsQBojPtcJKChOTd5o9HdPRceg3x5hvfUqEYFffwcEMA9DggzITOQCiGqn9mDbo2nzayC7RuCNKnpyM030exRwDV540F+R0FAsXqrbwV+wh6LzX1LlV0155xuhHuepaEV0tLxRh06ztpWS29/yh53gpoABJUccjK9P6wnVcxtydnbkpPYlPcvG6TnZYDlPTMvpvm9zwhT4Jp+/kP8C1GgHj3He5/CmOBWARWukbXjTAzIQ206lsdoofW9Dul4BET81L61uB/zVR4RSC5ghmLot26HmQIyGQVDBI+ompVDjAE78E9GI93y8eANuBb62VZFopyZ3l+taAGD+HCMI/QJh2yghpdzcCQsF6PvRDYWmPlcjQ0RR4ZYMNpW6vqePi6VibNcb/8De0q3cXphevy0ywcUbhA3uyG+TC8hDHyJLYe7bVT85esZ1gDoMJ5C8bmUwZpC4EstKm8N1Od2ibGLWVyPZUUv7+qkZnbJJ5YW0+OF01o23WzJUpx6wl5Iq0u5imTEkpmquA3xNvi2wAACei1ry7tHpGAAHoCJ4UAADIm4JPscRn+wIAAAAABFla")
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
    name = f"diagnostics_1d_{axis}_axis"
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
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AfCA3RdABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIwIlIQLEtn60ZvFpp58iRzPO8zL7jj2RkmJNOOz7kLOZ5pcKkv92RfBBV6XfvVRqJRIbNiOKE+HSFEwb7csuCUsf9zB0qgLjaNojvcXVzhXj5lrq9G7GuZvGWxSqMoQWYxmpol2GdvZjHEN5O4OgyeP4QzVStm900cNC9r0y2my7Zytkn+Mx0Uv6hkLF9g7wKDwqH9TrxpV3/flnyYZhM6g8tfR1stiKkD0DKMOTo7G5ZuDDLTICjLIdu3yu7lwu/4AROxXXEQIAaeZvy/SyQxDcQX3caFyp38uJ4U/5Jl/iFqDXToYx7ponlNPbWD6wOY3hYrUaTd5SYq5O9EnGB/SbSHSuwniqf1W6XwfJk+7VAdPFlQMzbsyczUIqqWMHEF9e20Hr4iqGi8Q7QKk82tCHojLFxGGC6pswiqOQHDvemybNDUpMN/Im3ke4NyZCS32eX5MEPwnLG2l0edf3MSL+eq6WhLu/dv8afpeARUK5nzH52Gy5RDYJ14lPqiqCRsu+MgLGUl0dm/e9t8MRXG6Z4k1uq197gxhK2nS0u8L1Vo9P+UiIcUIn6DvSBOCq92qxysUaG4QPwpdX5gkWI40vh933DYo1DWXODZaEOka4a7zdzkHwn5HsNP3GYoNw9HPvJRy225KYYgeiIqLZiMW7QG+bz3l1ZuYWfsDY4Zceee1F4s6TCCUzIWF3TVg4GKOdGzy+vok58JY2PBo/jh2jhDb8Duljt7RHctcWCNdKAPAkyh9HaSri//2SatbEZymgWanvE9xdok/z6O48qNIEiTaLCpgopjmg+YD5hk+BtJqmABxZ/UtbMGSCim90r5MIWxkw2jGJ3z5dJaz7QvClOnn1EhoX3SF0BocSmcqubJtm0vuWByrGI1iEwyIbOclGOxOwpkbfkWFQwZLLvV3Jz7NWtfVqksoYsLSVevYdy3GsiaTYTrQc7MCVf6bA2Gp5VtAtUsM94yoYAZWMByfSQ1hL7U0tQAeD/dDNwXvdumTMlM0Vs1PmfmXml4p54XjCg7gHqTeTRbO007MyAe0pDqH1CF8/XpZDep20cIrbp7ait57gFK50+j+GFOiSVVe7ao39YvGOFvhP3+RAAANDLN86BEvPDAAGQB8MPAABbGLcDscRn+wIAAAAABFla")
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
    name = f"diagnostics_2d_{plane}_plane"
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
