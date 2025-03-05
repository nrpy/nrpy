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

from typing import Dict, Tuple

import sympy as sp

import nrpy.c_function as cfc
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.simple_loop as lp


def register_CFunction_diagnostics_nearest_grid_center(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    filename_tuple: Tuple[str, str] = (
        "out0d-conv_factor%.2f.txt",
        "convergence_factor",
    ),
) -> (
    None
):  # Do NOT enable parallel codegen here; this should be called from a function that should already have pcg enabled.
    r"""
    Register C function for "0-dimensional" simulation diagnostics -- output data at gridpoint closest to grid center.

    This function generates a C function that computes and outputs specified diagnostic quantities at the grid point nearest to the physical center of the grid for a given coordinate system. The output is written to a file whose name is specified by `filename_tuple`.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary mapping (type, name) tuples to corresponding C definitions.
        Example: {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}
    :param filename_tuple: Tuple specifying the filename and its corresponding string format.
        Default: ("out0d-conv_factor%.2f.txt", "convergence_factor")
    :raises ValueError: If an unsupported coordinate system is specified, ensuring that diagnostics are only generated for coordinate systems with a defined grid center.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> Coord = "SinhCylindrical"
    >>> _ = register_CFunction_diagnostics_nearest_grid_center(Coord, out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"})
    >>> validate_strings(cfc.CFunction_dict[f"diagnostics_nearest_grid_center__rfm__{Coord}"].full_function, "grid_center")
    """
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

    body = rf"""
// Unpack grid function pointers from gridfuncs struct
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
MAYBE_UNUSED const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
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
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


def register_CFunction_diagnostics_nearest_1d_axis(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    axis: str,
    filename_tuple: Tuple[str, str] = (
        "out1d-AXIS-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
) -> (
    None
):  # Do NOT enable parallel codegen here; this should be called from a function that should already have pcg enabled.
    r"""
    Register C function for 1-dimensional simulation diagnostics at gridpoints closest to specified axis.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary of output quantities.
    :param axis: Specifies the axis ("x", "z") for the diagnostics.
    :param filename_tuple: Tuple containing the format for filename and the replacement arguments.
    :raises ValueError: If the specified axis is not supported.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> Coord = "SinhSpherical"
    >>> axis = "y"
    >>> _ = register_CFunction_diagnostics_nearest_1d_axis(Coord, out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis=axis)
    >>> validate_strings(cfc.CFunction_dict[f"diagnostics_nearest_1d_{axis}_axis__rfm__{Coord}"].full_function, "SinhSpherical_y_axis")
    """
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
    body = rf"""
// Unpack grid function pointers from gridfuncs struct
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
MAYBE_UNUSED const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

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


def register_CFunction_diagnostics_nearest_2d_plane(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    plane: str,
    filename_tuple: Tuple[str, str] = (
        "out2d-PLANE-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
) -> (
    None
):  # Do NOT enable parallel codegen here; this should be called from a function that should already have pcg enabled.
    r"""
    Register C function for 2-dimensional simulation diagnostics at gridpoints closest to the specified plane.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary of output quantities.
    :param plane: Specifies the plane ("xy", "yz") for the diagnostics.
    :param filename_tuple: Tuple containing the format for filename and the replacement arguments.
    :raises ValueError: If the specified plane is not supported.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> Coord = "SinhSymTP"
    >>> plane = "yz"
    >>> _ = register_CFunction_diagnostics_nearest_2d_plane("SinhSymTP", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane=plane)
    >>> validate_strings(cfc.CFunction_dict[f"diagnostics_nearest_2d_{plane}_plane__rfm__{Coord}"].full_function, "SinhSymTP_yz_plane")
    """
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
    body = rf"""
// Unpack grid function pointers from gridfuncs struct
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
MAYBE_UNUSED const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
