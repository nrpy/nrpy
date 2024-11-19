"""
Output simulation data at points closest to grid physical center (0D), grid y or z axis (1D), and grid xy or yz plane (2D) for the superB infrastructure.

Functions:
----------
-register_CFunction_diagnostics_set_up_nearest_1d_axis
   set up diagnosticstruct containing the number of diagnostic points, their index and the offset in the file where to write 1-dimensional simulation diagnostics, etc

-register_CFunction_diagnostics_set_up_nearest_2d_plane
   set up diagnosticstruct containing the number of diagnostic points, their index and the offset in the file where to write 2-dimensional simulation diagnostics, etc

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
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, Tuple, Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.superB.simple_loop_diagnostic as lp


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
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> Coord = "SinhCylindrical"
    >>> _ = register_CFunction_diagnostics_nearest_grid_center(Coord, out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"})
    >>> validate_strings(cfc.CFunction_dict[f"diagnostics_nearest_grid_center__rfm__{Coord}"].full_function, "grid_center")
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
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for 1-dimensional simulation diagnostics at gridpoints closest to specified axis.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary of output quantities.
    :param axis: Specifies the axis ("x", "z") for the diagnostics.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If the specified axis is not supported.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    if axis not in ["y", "z"]:
        raise ValueError(
            f"Output along {axis} axis not supported. Please choose x or z axis."
        )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"Output diagnostic quantities at gridpoints closest to {axis} axis."
    cfunc_type = "void"
    name = f"diagnostics_nearest_1d_{axis}_axis"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs, const diagnostic_struct *restrict diagnosticstruct, Ck::IO::Session token"

    body = rf"""
// Unpack gridfuncs struct:
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

const int num_diagnostic_pts = {'diagnosticstruct->num_diagnostic_1d_y_pts;' if axis == "y" else 'diagnosticstruct->num_diagnostic_1d_z_pts;'}
const int *restrict idx3_diagnostic_pt = {'diagnosticstruct->localidx3_diagnostic_1d_y_pt;' if axis == "y" else 'diagnosticstruct->localidx3_diagnostic_1d_z_pt;'}
const int *restrict i0_diagnostic_pt = {'diagnosticstruct->locali0_diagnostic_1d_y_pt;' if axis == "y" else 'diagnosticstruct->locali0_diagnostic_1d_z_pt;'}
const int *restrict i1_diagnostic_pt = {'diagnosticstruct->locali1_diagnostic_1d_y_pt;' if axis == "y" else 'diagnosticstruct->locali1_diagnostic_1d_z_pt;'}
const int *restrict i2_diagnostic_pt = {'diagnosticstruct->locali2_diagnostic_1d_y_pt;' if axis == "y" else 'diagnosticstruct->locali2_diagnostic_1d_z_pt;'}
const int *restrict offsetpt_firstfield = {'diagnosticstruct->offset_diagnostic_1d_y_pt;' if axis == "y" else 'diagnosticstruct->offset_diagnostic_1d_z_pt;'}
for (int which_pt = 0; which_pt < num_diagnostic_pts; which_pt++) {{
  const int idx3 = idx3_diagnostic_pt[which_pt];
  const int i0 = i0_diagnostic_pt[which_pt];
  const int i1 = i1_diagnostic_pt[which_pt];
  const int i2 = i2_diagnostic_pt[which_pt];
  REAL xCart[3];
  xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
  const REAL xCart_axis = {'xCart[1];' if axis == "y" else 'xCart[2];'}
  int sizeinbytes = 23 * (diagnosticstruct->num_output_quantities + 1);
  char out[sizeinbytes+1];
  snprintf(out, sizeof(out),"""
    output_to_file = """ "% .15e """
    for key in out_quantities_dict.keys():
        printf_c_type = "% .15e" if key[0] != "int" else "%d"
        output_to_file += f"{printf_c_type} "

    output_to_file = f'{output_to_file[:-1]}\\n", xCart_axis, '
    for value in out_quantities_dict.values():
        output_to_file += f"{value}, "
    output_to_file = f"{output_to_file[:-2]});\n"

    body += output_to_file

    body += r"""
  Ck::IO::write(token, out, sizeinbytes, offsetpt_firstfield[which_pt]);
}
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


def register_CFunction_diagnostics_nearest_2d_plane(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    plane: str,
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for 2-dimensional simulation diagnostics at gridpoints closest to the specified plane.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary of output quantities.
    :param plane: Specifies the plane ("xy", "yz") for the diagnostics.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If the specified plane is not supported.
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
    cfunc_type = "void"
    name = f"diagnostics_nearest_2d_{plane}_plane"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs, const diagnostic_struct *restrict diagnosticstruct, Ck::IO::Session token"

    body = rf"""
// Unpack gridfuncs struct:
const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;
// Unpack diagnosticptoffset struct:
const int num_diagnostic_pts = {'diagnosticstruct->num_diagnostic_2d_xy_pts;' if plane == "xy" else 'diagnosticstruct->num_diagnostic_2d_yz_pts;'}
const int *restrict idx3_diagnostic_pt = {'diagnosticstruct->localidx3_diagnostic_2d_xy_pt;' if plane == "xy" else 'diagnosticstruct->localidx3_diagnostic_2d_yz_pt;'}
const int *restrict i0_diagnostic_pt = {'diagnosticstruct->locali0_diagnostic_2d_xy_pt;' if plane == "xy" else 'diagnosticstruct->locali0_diagnostic_2d_yz_pt;'}
const int *restrict i1_diagnostic_pt = {'diagnosticstruct->locali1_diagnostic_2d_xy_pt;' if plane == "xy" else 'diagnosticstruct->locali1_diagnostic_2d_yz_pt;'}
const int *restrict i2_diagnostic_pt = {'diagnosticstruct->locali2_diagnostic_2d_xy_pt;' if plane == "xy" else 'diagnosticstruct->locali2_diagnostic_2d_yz_pt;'}
const int *restrict offsetpt_firstfield = {'diagnosticstruct->offset_diagnostic_2d_xy_pt;' if plane == "xy" else 'diagnosticstruct->offset_diagnostic_2d_yz_pt;'}
for (int which_pt = 0; which_pt < num_diagnostic_pts; which_pt++) {{
  const int idx3 = idx3_diagnostic_pt[which_pt];
  const int i0 = i0_diagnostic_pt[which_pt];
  const int i1 = i1_diagnostic_pt[which_pt];
  const int i2 = i2_diagnostic_pt[which_pt];
  REAL xCart[3];
  xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
  int sizeinbytes = 23 * (diagnosticstruct->num_output_quantities + 2);
  char out[sizeinbytes+1];
  snprintf(out, sizeof(out),"""
    output_to_file = """ "% .15e % .15e """
    for key in out_quantities_dict.keys():
        printf_c_type = "% .15e" if key[0] != "int" else "%d"
        output_to_file += f"{printf_c_type} "

    if plane == "xy":
        output_to_file = f'{output_to_file[:-1]}\\n", xCart[0], xCart[1], '
    elif plane == "yz":
        output_to_file = f'{output_to_file[:-1]}\\n", xCart[1], xCart[2], '
    for value in out_quantities_dict.values():
        output_to_file += f"{value}, "
    output_to_file = f"{output_to_file[:-2]});\n"

    body += output_to_file

    body += r"""
  Ck::IO::write(token, out, sizeinbytes, offsetpt_firstfield[which_pt]);
}
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


def register_CFunction_diagnostics_set_up_nearest_1d_axis(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    axis: str,
    filename_tuple: Tuple[str, str] = (
        "out1d-AXIS-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for setting up diagnostic struct for 1-dimensional simulation diagnostics at gridpoints closest to specified axis.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary of output quantities.
    :param axis: Specifies the axis ("x", "z") for the diagnostics.
    :param filename_tuple: Tuple containing the format for filename and the replacement arguments.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If the specified axis is not supported.
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
        axis=axis,
    )

    desc = f"Setup diagnostic quantities at gridpoints closest to {axis} axis."
    cfunc_type = "void"
    name = f"diagnosticstruct_set_up_nearest_1d_{axis}_axis"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, const charecomm_struct *restrict charecommstruct, REAL *restrict xx[3], const int chare_index[3], diagnostic_struct *restrict diagnosticstruct"

    body = rf"""
const int Nchare0 = commondata->Nchare0;
const int Nchare1 = commondata->Nchare1;
const int Nchare2 = commondata->Nchare2;
const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
const int Nxx0 = params->Nxx0;
const int Nxx1 = params->Nxx1;
const int Nxx2 = params->Nxx2;
const int Nxx_plus_2NGHOSTS0chare = params_chare->Nxx_plus_2NGHOSTS0;
const int Nxx_plus_2NGHOSTS1chare = params_chare->Nxx_plus_2NGHOSTS1;
const int Nxx_plus_2NGHOSTS2chare = params_chare->Nxx_plus_2NGHOSTS2;
const int Nxx0chare = params_chare->Nxx0;
const int Nxx1chare = params_chare->Nxx1;
const int Nxx2chare = params_chare->Nxx2;

strcpy(diagnosticstruct->filename_1d_{axis}, "{filename_tuple[0].replace('AXIS', axis)}");

diagnosticstruct->num_output_quantities = {len(out_quantities_dict)};

{loop_1d}

diagnosticstruct->tot_num_diagnostic_1d_{axis}_pts = data_index;

int num_diagnostics_chare = 0;
for (int i = 0; i < data_index; i++) {{
  const int i0 = data_points[i].i0;
  const int i1 = data_points[i].i1;
  const int i2 = data_points[i].i2;
  const int idx3 = IDX3(i0, i1, i2);
  if (charecommstruct->globalidx3pt_to_chareidx3[idx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {{
    num_diagnostics_chare++;
  }}
}}
diagnosticstruct->num_diagnostic_1d_{axis}_pts = num_diagnostics_chare;
diagnosticstruct->localidx3_diagnostic_1d_{axis}_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
diagnosticstruct->locali0_diagnostic_1d_{axis}_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
diagnosticstruct->locali1_diagnostic_1d_{axis}_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
diagnosticstruct->locali2_diagnostic_1d_{axis}_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
diagnosticstruct->offset_diagnostic_1d_{axis}_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);

// compute offset in bytes for first field for each diagnostic pt
int sizeinbytes = 23 * (diagnosticstruct->num_output_quantities + 1);
int which_diagnostics_chare = 0;
int which_diagnostic_global = 0;
for (int i = 0; i < data_index; i++) {{
  const int i0 = data_points[i].i0;
  const int i1 = data_points[i].i1;
  const int i2 = data_points[i].i2;
  const int idx3 = IDX3(i0, i1, i2);
  if (charecommstruct->globalidx3pt_to_chareidx3[idx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])){{
    int localidx3 = charecommstruct->globalidx3pt_to_localidx3pt[idx3];
    diagnosticstruct->localidx3_diagnostic_1d_{axis}_pt[which_diagnostics_chare] = localidx3;
    diagnosticstruct->locali0_diagnostic_1d_{axis}_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0, Nxx0chare);
    diagnosticstruct->locali1_diagnostic_1d_{axis}_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1, Nxx1chare);
    diagnosticstruct->locali2_diagnostic_1d_{axis}_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2, Nxx2chare);
    diagnosticstruct->offset_diagnostic_1d_{axis}_pt[which_diagnostics_chare] = which_diagnostic_global * sizeinbytes;
    which_diagnostics_chare++;
  }}
  which_diagnostic_global++;
}}
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_diagnostics_set_up_nearest_2d_plane(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    plane: str,
    filename_tuple: Tuple[str, str] = (
        "out2d-PLANE-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
) -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for setting up diagnostic struct for 2-dimensional simulation diagnostics at gridpoints closest to specified plane.

    :param CoordSystem: Specifies the coordinate system for the diagnostics.
    :param out_quantities_dict: Dictionary of output quantities.
    :param plane: Specifies the plane ("xy", "yz") for the diagnostics.
    :param filename_tuple: Tuple containing the format for filename and the replacement arguments.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If the specified plane is not supported.
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
        plane=plane,
    )

    desc = f"Set up diagnostic quantities at gridpoints closest to {plane} plane."
    cfunc_type = "void"
    name = f"diagnosticstruct_set_up_nearest_2d_{plane}_plane"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, const charecomm_struct *restrict charecommstruct, REAL *restrict xx[3], const int chare_index[3], diagnostic_struct *restrict diagnosticstruct"

    body = rf"""
const int Nchare0 = commondata->Nchare0;
const int Nchare1 = commondata->Nchare1;
const int Nchare2 = commondata->Nchare2;
const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
const int Nxx0 = params->Nxx0;
const int Nxx1 = params->Nxx1;
const int Nxx2 = params->Nxx2;
const int Nxx_plus_2NGHOSTS0chare = params_chare->Nxx_plus_2NGHOSTS0;
const int Nxx_plus_2NGHOSTS1chare = params_chare->Nxx_plus_2NGHOSTS1;
const int Nxx_plus_2NGHOSTS2chare = params_chare->Nxx_plus_2NGHOSTS2;
const int Nxx0chare = params_chare->Nxx0;
const int Nxx1chare = params_chare->Nxx1;
const int Nxx2chare = params_chare->Nxx2;

strcpy(diagnosticstruct->filename_2d_{plane}, "{filename_tuple[0].replace('PLANE', plane)}");

diagnosticstruct->num_output_quantities = {len(out_quantities_dict)};

{loop_2d}

diagnosticstruct->tot_num_diagnostic_2d_{plane}_pts = numpts_i0 * numpts_i1 * numpts_i2;

int num_diagnostics_chare = 0;
LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {{
  const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
  const int idx3 = IDX3(i0, i1, i2);
  if (charecommstruct->globalidx3pt_to_chareidx3[idx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])) {{
    num_diagnostics_chare++;
  }}
}}
diagnosticstruct->num_diagnostic_2d_{plane}_pts = num_diagnostics_chare;
diagnosticstruct->locali0_diagnostic_2d_{plane}_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
diagnosticstruct->locali1_diagnostic_2d_{plane}_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
diagnosticstruct->locali2_diagnostic_2d_{plane}_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
diagnosticstruct->localidx3_diagnostic_2d_{plane}_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
diagnosticstruct->offset_diagnostic_2d_{plane}_pt = (int *restrict)malloc(sizeof(int) * num_diagnostics_chare);
// compute offset in bytes for first field for each diagnostic pt
int sizeinbytes = 23 * (diagnosticstruct->num_output_quantities + 2);

int which_diagnostics_chare = 0;
int which_diagnostic_global = 0;
LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {{
  const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
  const int idx3 = IDX3(i0, i1, i2);
  if (charecommstruct->globalidx3pt_to_chareidx3[idx3] == IDX3_OF_CHARE(chare_index[0], chare_index[1], chare_index[2])){{
    // store the local idx3 of diagnostic point
    int localidx3 = charecommstruct->globalidx3pt_to_localidx3pt[idx3];
    diagnosticstruct->localidx3_diagnostic_2d_{plane}_pt[which_diagnostics_chare] = localidx3;
    diagnosticstruct->locali0_diagnostic_2d_{plane}_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i0, Nxx0chare);
    diagnosticstruct->locali1_diagnostic_2d_{plane}_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX1(chare_index[1], i1, Nxx1chare);
    diagnosticstruct->locali2_diagnostic_2d_{plane}_pt[which_diagnostics_chare] = MAP_GLOBAL_TO_LOCAL_IDX2(chare_index[2], i2, Nxx2chare);
    diagnosticstruct->offset_diagnostic_2d_{plane}_pt[which_diagnostics_chare] = which_diagnostic_global * sizeinbytes;
    which_diagnostics_chare++;
  }}
  which_diagnostic_global++;
}}
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
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
