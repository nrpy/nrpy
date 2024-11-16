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
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for "0-dimensional" simulation diagnostics -- output data at gridpoint closest to grid center.

    This function generates a C function that computes and outputs specified diagnostic quantities at the grid point nearest to the physical center of the grid for a given coordinate system. The output is written to a file whose name is specified by `filename_tuple`.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary mapping (type, name) tuples to corresponding C definitions.
        Example: {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}
    :param filename_tuple: Tuple specifying the filename and its corresponding string format.
        Default: ("out0d-conv_factor%.2f.txt", "convergence_factor")

    :return: None if in registration phase, else the updated NRPy environment.

    :raises ValueError: If an unsupported coordinate system is specified, ensuring that diagnostics are only generated for coordinate systems with a defined grid center.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> Coord = "SinhCylindrical"
    >>> axis = "y"
    >>> _ = register_CFunction_diagnostics_nearest_grid_center(Coord, out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"})
    >>> diag_gc = cfc.CFunction_dict[f"diagnostics_nearest_grid_center__rfm__{Coord}"].full_function
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AbnAv9dABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIwLVzEUvxTswCb6RR3dNTDk3j3S9iEPKV1hatbjHcmANRICzOBgbTAnR78qVNJJRYOak9sWdY+3cVOWMqbAFyHKE27fHWpYmU2iNU4iepIRlFITM+uMgvrnHBovCkYnRffA4WOgNtH3At+TSthGMDmZzmTLd4HbEj9JS3BEQX/tHO2F1aEFMZ/gs4e/3sbowpVFC17WK2uQ2JM3eYOh26DtXJrxXjBl6YCy2w3vlYZkAnHvZ0/7pTQvsEnZbsOkZLDy7QNsFZBgCGaFXAITp85iemhIUC3adaGrAErXDrIQV1n9cAj87lMdaxj8Gx+4MLgml+hHZDnjHSmvrHLJq8GeBQ0NBgJu6kc4hOL4HIKkS/TVD1ieKadUzpPNQMF8JObzURxIOdiMtvRH48Fw2ZcoLA9m4NjwBSk2unCTLm1r11OhORzVAwVeHbdUM1+ASHgTFwFenBHHiSknBYaquWu4BvSH6PbzvHKW4aQ4/1rx84k+v4tFb/lBwYO8xbWEK3+vL1UprmhRGPKyfDRlbCfETrot+/vTA7oFdRa+ko8U0OybedvmkRSdbXTt77urEmzwFXjTUdenhw7f9bZ0thpCFhhUadrpa+I6QbeJbYriSvctV/IZwh+TwrW7KmHQm0TTaSqma6BL8TSy3JV6GmHhQQwh0w4rwrMMAVtrr7eeemo2MPGoi/AHPJVFgWXg8aSkAzaosxttck2L4Gye1njtOtHeqISZnskDOzBzvsfvjv+iAaNkJXn/xgqz/6tvZBydsIZ/rWC5k586PU0u7r1RuskzrtKQcngx6gl/ZV/fZJGY0TGMI4tFp7afIdSnhAY0RSpY1WEzIcqJFwEaPNhhLc0qi4DClfi9KSle92x3EXm328Hsh3uIK6451CNyK6//phJ1oeAUeKWGM97DPGIh2w4qDf1ytYSIOmRRdgKH3q4hL1c4LtpuJgAAChKSyNZnC9qwABmwboDQAAeT2nILHEZ/sCAAAAAARZWg==")
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

    body = rf"""
// Unpack grid function pointers from gridfuncs struct
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const MAYBE_UNUSED REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
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
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AsMBLVdABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIzLaV4Uh8eHMdkFke5ZEWqOrJB8gRVzpSQFQuNk6de9OlmLCKH/0/3hkEO+FpHWMRaKtUzxXmXROgszzTA582kYrIga+g7jqdNJxR8QmVOk3X41Z00Qr8cTShS0e4Ve/Rsfwz7AFvRk7u+g/DHsZZvuk1HMz/N+rflpetPQKKDsLx7wMAnJcaXnErnfc94NuhxGxcaLALyAW85ZpxXG+WHqDC6eXJKplcJ6yp1Ceh1tAXGQuwOP1JDjAbnb4He5NpjZrujWlcCgshTP7DZeP2k2eVVqTkCUBVRTb086HVPM7jWvrVKuNll6qXgzmszOnXAxzIHw7wWr2syW+5eglF4s0B/ftqmfnf+np32rN/2EhRAVu9WyBH93UwDfAYYXa49UvHG3wWktQvLZjGH+ZK3jufFDjrIN22xVhLAA9ESsznhTNcla1Fuv6LywPG8ms4dtrefSV2p+SU1wnpEYCS369oLQ2Rh5o+2AAmiES4nDSLXZsCYQ2n+j9OQvouuvHaY2zHLS1aqSlBbIJYRYvkxCPNeGJUHTb+i2VSc3Q4jUJ0WFPCI4KvDRrenxUnUVnreepEW4D7/GzTfKRUnKwwOwbo6IlheaOwpxaj70KM5n0KRp6FzlGNlSa+W5QMOB8YtuklR1IMMM8pwfo3Cqj7bh+iJK2hu/kZBY5cqszwXBxtjNwr9ObQw2Pt6RA8cEcqzzxGRHzgOluDxmWF+jourCZDigYa5+13CNV5p0iD/+OPVX+MB7+SIna+O6nHw+P/rPddI7MuqjrhkGErL68RXkfxLrKPQ53JAoKJwpCxgOTWMbRRQhTzWk5yAxVYJCQnJHuJPmsqqNEc1utHk3xgPo3h95kGFGYo8ASgTG3B2A586BNeT5iLjsb6LF9SZJ9INcTAnnDXnrAMnMDreKh3zq0DcG3nsQPNayNeBykaR2hKr3VSjSWjpP4TBla1Z9ruB3rc6sD/suE2+wtAq9ZGMn99nrLRyHwUjSexADACB7XsmuJ63wQVy8TTvv/mS22/NF/V3zTDJTJh4ci/j8s3dia8RXDFqG3bmN07V4frgjYSg7ah7iqMBPAwN+sGGkeiFg62zBlijkdxkOzRnfBplkKjLa50h+JwhRC3SjGGc0IAt7Ia2UrlxRCsdAzaPqeUFhxtzawIyPvkBVCqijmi92VOU+ZkeK8pXFuTB1JwZlPZuHF3ASUg3DvbO/G/RfKnKa41Z5pWjnH/LV1GaiDZiu7jkYrWIM/ormvWi03c4gj+Ncwbs18g5/FSLqrpaSyvSTF/3msSnVLVe1b851D1I/dEFT7FCyYNKyV31hgSHmpKAQcMYBoy5CzMv6jJO71XyNDwTpa8S8edJFMHKvsgdcPOA6XtoRHJ38u56nMZrSVVeCHUMuBZgME6SndQV23D9BP1L27Ohv9zWBTmI1coNctQOEkHxas9NZ3QGvvIojoH1wqVF7nwI+foN4kRuU56A1j25zuGmM8CIHFhzxGCkPsznoGrW1bXp/ZPnKdHl7lCHZFWiZ6nbz4IzmaxIxEAAAAAAGWHCFuMaLeMAAHRCY0WAAACyL8QscRn+wIAAAAABFla")
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

    body = rf"""
// Unpack grid function pointers from gridfuncs struct
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const MAYBE_UNUSED REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
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
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4Ah7A9FdABGaScZHDxOiAHcc747vJCBLAVaWVRAnpWHX/kkr6QysfzLfXhoCqwdlhv59jJtBxk7nXnUIwLVzEUvxTswCb6RR3dNTDk3j3S9iEPKV1hatbjHcmANRICzOBge4rbcsq6HU2FEY6lkRkhN2vxwKCZMmULt85t8Ko3B4p2ztr6zGycakHm9rPjzEK5S5nwG+fnIQiwqx/SeaeG3zWsx4XNyAuxPfU33La+338O3FTQ3CwCWkE2ZYGfte96onLd1i0TuhNN5W70qN6VnFdXevqx5hJoXtSWqj86UzfWGYC+VR//6OsvehqwTisgLfQTyESmB6CpliP6Ny5GepAjdBb/I5sIspZiI+3v2eKedlqB6etSFPn0M2uNC97D+eCFQD4ns9W5g/qkBkAldqULkqcgLhV6mkA5CilTQDc0O4N4+hPuaTcgl5dFWVQntiAR7FIz18dqpo9kUYKVDT+xenKe4B6IIzbbbzIIKWDJmCdPbu1uKIdRHiUsNUBqfZ2wsOw962/hoNn4l7ps32KeNxeanz6m2eINavqRvIp1VbryEyvMf4ttj5aTcN6+9bo6ZZQghP+vHMjDPbQ9VkgJ28cmZG9QwP7c8Ve2SdGsCB5wLwKpdvgePfGp9shwPDtnJ08RJoLRZmEzfIY3r8DluXYHrWJ+Rq8ARuXrk5ratWanwQ0mMXD7R/tPZCci9Cp49k3qYzeEhVzUsFibimxw1DDdAimboifplyOUBvhgVK76tvxWcQAm34Xg2u0FlbUJvbpRh3b+Uq4D1cEdGQe/o3qwC0I/VJFn4rfha3Q5r7zxrk9sFCxHGFIahlGY2aZI022B300ovo51xaSWRWUyHtSidulN35hsbIQWIsPSTBWG8D2bf58VQ7wLxX6J1Lf8RWUTp4tvePaTOx44/sIknJ0AZQyB006dCtD2JwU92Oh5Oaq/7ug08yzVPUMpEWTySFhCf22OsyIoX+evxu4XD3HfkZdhRaNXfQNpjtBzgL5F0CIQWwHo6ojAngEmSycmYRbYDyvSAm8MRYmgnsUl669i85Hu+mWpgeeD+UWN0FTiFstBBdXcEYNSTS/IGCR159qPQq5pbDcU1oJy2xWFs3eThsR1WosBLfuhJ/bcZIAkx77znH9A8NQ9q8D1i9tGUqiYLdNC5h2ElqFx23NQiSJycLqfARi4KE2nfRpOqZjwUB1CD+K1OOrrg/j96TLX6MHC5xww9C90UsfP22vd3qQ2nqIjfoO6SsQUWEydGdLugYqMVLSyOCvlGMewKvt70+I1sySsn/WyxEo2Y1mBmzpgAAAAAAAH/XwTcegdJ+AAHtB/wQAAD+CjkLscRn+wIAAAAABFla")
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

    body = rf"""
// Unpack grid function pointers from gridfuncs struct
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const MAYBE_UNUSED REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
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
