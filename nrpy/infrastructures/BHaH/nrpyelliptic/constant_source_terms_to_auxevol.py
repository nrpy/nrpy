"""
C function for setting constant auxiliary terms for hyperbolic relaxation.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
from nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms import (
    compute_psi_background_and_ADD_times_AUU,
)
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list


# Define functions to set AUXEVOL gridfunctions
def register_CFunction_auxevol_gfs_single_point(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for the AUXEVOL grid functions at a single point.

    :param CoordSystem: The coordinate system to use in setting up the AUXEVOL gridfunctions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Compute psi_background and ADD_times_AUU
    psi_background, ADD_times_AUU = compute_psi_background_and_ADD_times_AUU(
        CoordSystem
    )

    includes = ["BHaH_defines.h"]
    parallelization = par.parval_from_str("parallelization")

    desc = r"""Compute AUXEVOL grid functions at a single point."""
    cfunc_type = "void"
    name = "auxevol_gfs_single_point"
    params = r"""const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL xx0, const REAL xx1, const REAL xx2,  REAL *restrict psi_background, REAL *restrict ADD_times_AUU
""".replace(
        "const commondata_struct *restrict commondata, const params_struct *restrict params,",
        (
            "const size_t streamid,"
            if parallelization in ["cuda"]
            else "const commondata_struct *restrict commondata, const params_struct *restrict params,"
        ),
    )
    body = ""
    expr_list = [psi_background, ADD_times_AUU]
    params_symbols, commondata_symbols = get_params_commondata_symbols_from_expr_list(
        expr_list, exclude=[f"xx{i}" for i in range(3)]
    )
    body += "// Load necessary parameters from params_struct\n"
    for param in params_symbols:
        body += f"const REAL {param} = {parallel_utils.get_params_access(parallelization)}{param};\n"
    body += "\n// Load necessary parameters from commondata_struct\n"
    for param in commondata_symbols:
        body += f"const REAL {param} = {parallel_utils.get_commondata_access(parallelization)}{param};\n"
    body += "\n"
    body += ccg.c_codegen(
        expr_list,
        ["*psi_background", "*ADD_times_AUU"],
        verbose=False,
        include_braces=False,
    )
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        cfunc_decorators="__device__" if parallelization in ["cuda"] else "",
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_auxevol_gfs_all_points(
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for the AUXEVOL grid functions at all points.

    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Set AUXEVOL gridfunctions at all points."""
    cfunc_type = "void"
    name = "auxevol_gfs_all_points"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    psi_background_memaccess = gri.BHaHGridFunction.access_gf("psi_background")
    ADD_times_AUU_memaccess = gri.BHaHGridFunction.access_gf("ADD_times_AUU")

    body = (
        "cpyHosttoDevice_commondata__constant(commondata);\n"
        if parallelization in ["cuda"]
        else ""
    )
    body += r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
  REAL *restrict x0 = griddata[grid].xx[0];
  REAL *restrict x1 = griddata[grid].xx[1];
  REAL *restrict x2 = griddata[grid].xx[2];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.auxevol_gfs;
"""

    kernel_body = f"{parallel_utils.get_loop_parameters(parallelization)}\n"
    function_call = (
        "auxevol_gfs_single_point(streamid, xx0,xx1,xx2,"
        if parallelization == "cuda"
        else "auxevol_gfs_single_point(commondata, params, xx0,xx1,xx2,"
    )
    kernel_body += lp.simple_loop(
        f"{function_call}"
        f"&{psi_background_memaccess},"
        f"&{ADD_times_AUU_memaccess});",
        read_xxs=True,
        loop_region="all points",
        OMP_collapse=OMP_collapse,
    )
    for i in range(3):
        kernel_body = kernel_body.replace(f"xx[{i}]", f"x{i}")

    comments = "Kernel to initialize auxillary grid functions at all grid points."
    # Prepare the argument dicts
    arg_dict_cuda = {
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "REAL *restrict",
    }
    arg_dict_host = {
        "commondata": "const commondata_struct *restrict",
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }
    prefunc, new_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=comments,
        launch_dict={
            "blocks_per_grid": [],
            "threads_per_block": ["32", "NGHOSTS"],
            "stream": "default",
        },
        thread_tiling_macro_suffix="NELL_AUX",
    )

    body += f"{new_body}\n"
    body += "}\n"

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_initialize_constant_auxevol() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to call all functions that set up AUXEVOL gridfunctions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Call functions that set up all AUXEVOL gridfunctions."""
    cfunc_type = "void"
    name = "initialize_constant_auxevol"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    body = r"""
    // Set up variable wavespeed
    variable_wavespeed_gfs_all_points(commondata, griddata);

    // Set up all other AUXEVOL gridfunctions
    auxevol_gfs_all_points(commondata, griddata);
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
