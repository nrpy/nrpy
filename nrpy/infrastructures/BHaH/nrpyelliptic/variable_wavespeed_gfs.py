"""
C function for computing variable wavespeed for hyperbolic relaxation.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list


def register_CFunction_variable_wavespeed_gfs_all_points(
    CoordSystem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to compute variable wavespeed based on local grid spacing for a single coordinate system.

    :param CoordSystem: The coordinate system to use in the hyperbolic relaxation.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h"]
    desc = "Compute variable wavespeed for all grids based on local grid spacing."
    cfunc_type = "void"
    name = "variable_wavespeed_gfs_all_points"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    rfm = refmetric.reference_metric[CoordSystem]
    dxx0, dxx1, dxx2 = sp.symbols("dxx0 dxx1 dxx2", real=True)
    expr_list = [
        rfm.scalefactor_orthog[0] * dxx0,
        rfm.scalefactor_orthog[1] * dxx1,
        rfm.scalefactor_orthog[2] * dxx2,
    ]
    dsmin_computation_str = ccg.c_codegen(
        expr_list,
        ["const REAL dsmin0", "const REAL dsmin1", "const REAL dsmin2"],
        include_braces=False,
    )

    variable_wavespeed_memaccess = gri.BHaHGridFunction.access_gf("variable_wavespeed")

    dsmin_computation_str += f"""\n// Set local wavespeed
        {variable_wavespeed_memaccess} = MINIMUM_GLOBAL_WAVESPEED * MIN(dsmin0, MIN(dsmin1, dsmin2)) / dt;\n"""

    body = r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict x0 = griddata[grid].xx[0];
  REAL *restrict x1 = griddata[grid].xx[1];
  REAL *restrict x2 = griddata[grid].xx[2];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.auxevol_gfs;
"""

    loop_body = lp.simple_loop(
        loop_body="\n" + dsmin_computation_str,
        read_xxs=True,
        loop_region="interior",
        CoordSystem=CoordSystem,
    )
    for i in range(3):
        loop_body = loop_body.replace(f"xx[{i}]", f"x{i}")

    parallelization = par.parval_from_str("parallelization")
    loop_params = parallel_utils.get_loop_parameters(
        parallelization,
        enable_intrinsics=(parallelization in ["cuda"]),
    )

    params_symbols, _ = get_params_commondata_symbols_from_expr_list(
        expr_list, exclude=[f"xx{i}" for i in range(3)]
    )
    loop_params += "// Load necessary parameters from params_struct\n"
    for param in params_symbols:
        loop_params += f"const REAL {param} = {parallel_utils.get_params_access(parallelization)}{param};\n"
    loop_params += "\n"

    comments = "Kernel to initialize auxillary grid functions at all grid points."
    # Prepare the argument dicts
    arg_dict_cuda = {
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "REAL *restrict",
        "dt": "const REAL",
        "MINIMUM_GLOBAL_WAVESPEED": "const REAL",
    }
    arg_dict_host = {
        # "commondata": "const commondata_struct *restrict",
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }
    kernel_body = f"{loop_params}\n{loop_body}"
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
        thread_tiling_macro_suffix="NELL_WAVESPEED",
    )

    body += f"{new_body}\n"
    # for param in commondata_symbols:
    #     body += f"const REAL {param} = {parallel_utils.get_commondata_access(parallelization)}{param};\n"
    # body += "\n"

    # We must close the loop that was opened in the line 'for(int grid=0; grid<commondata->NUMGRIDS; grid++) {'
    body += r"""} // END LOOP for(int grid=0; grid<commondata->NUMGRIDS; grid++)
            """

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,  # set_CodeParameters.h is manually included after the declaration of params_struct *restrict params
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
