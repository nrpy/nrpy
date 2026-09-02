"""
Generate the local minimum physical grid-spacing auxiliary gridfunction.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Set, Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)
from nrpy.infrastructures import BHaH
from nrpy.infrastructures.BHaH.numerical_grids_and_timestep import (
    ds_min_single_pt_exprs,
)


def register_CFunction_dsmin_auxevol_gridfunction(
    set_of_CoordSystems: Set[str],
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a function that stores raw local minimum physical grid spacing.

    :param set_of_CoordSystems: Coordinate systems used.

    :return: None if in registration phase, else the updated NRPy environment.

    :raises ValueError: If a non-fisheye GeneralRFM coordinate system is requested.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # dxx{0,1,2} must be registered or they will not be declared in the kernel.
    # fmt: off
    for j in range(3):
        _ = par.CodeParameter("REAL", __name__, f"dxx{j}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    # fmt: on

    parallelization = par.parval_from_str("parallelization")
    for CoordSystem in set_of_CoordSystems:
        desc = "Initialize the local minimum physical grid-spacing gridfunction."
        name = "dsmin_auxevol_gridfunction"
        params = "const params_struct *restrict params, REAL *restrict xx[3], REAL *restrict auxevol_gfs"
        arg_dict_cuda = {
            "x0": "const REAL *restrict",
            "x1": "const REAL *restrict",
            "x2": "const REAL *restrict",
            "auxevol_gfs": "REAL *restrict",
        }
        arg_dict_host = {
            "params": "const params_struct *restrict",
            **arg_dict_cuda,
        }

        expr_list = ds_min_single_pt_exprs(CoordSystem)
        if expr_list is None:
            raise ValueError(
                "ds_min for non-fisheye GeneralRFM coordinate system "
                f"{CoordSystem} is not yet supported."
            )
        loop_body = "REAL dsmin0, dsmin1, dsmin2;\n"
        loop_body += ccg.c_codegen(
            expr_list,
            ["dsmin0", "dsmin1", "dsmin2"],
            include_braces=False,
        )
        loop_body += (
            "auxevol_gfs[IDX4(DSMINGF, i0, i1, i2)] = "
            "NRPYMIN(dsmin0, NRPYMIN(dsmin1, dsmin2));"
        )

        param_symbols, _ = get_params_commondata_symbols_from_expr_list(
            expr_list, exclude=[f"xx{j}" for j in range(3)]
        )
        loop_params = parallel_utils.get_loop_parameters(parallelization)
        params_definitions = generate_definition_header(
            param_symbols,
            var_access=parallel_utils.get_params_access(parallelization),
        )
        kernel_body = f"{loop_params}\n{params_definitions}\n"
        kernel_body += BHaH.simple_loop.simple_loop(
            loop_body=loop_body,
            loop_region="all points",
            read_xxs=True,
        )
        for i in range(3):
            kernel_body = kernel_body.replace(f"xx[{i}]", f"x{i}")

        prefunc, launch_body = parallel_utils.generate_kernel_and_launch_code(
            name,
            kernel_body.replace(
                "SIMD", "CUDA" if parallelization == "cuda" else "SIMD"
            ),
            arg_dict_cuda,
            arg_dict_host,
            parallelization=parallelization,
            comments=desc,
            cfunc_type="static void",
            launchblock_with_braces=False,
            thread_tiling_macro_suffix="DSMIN",
        )
        for i in range(3):
            launch_body = launch_body.replace(f"x{i}", f"xx[{i}]")

        cfc.register_CFunction(
            prefunc=prefunc,
            includes=["BHaH_defines.h"],
            desc=desc,
            name=name,
            params=params,
            include_CodeParameters_h=False,
            body=launch_body,
            CoordSystem_for_wrapper_func=CoordSystem,
        )
    return pcg.NRPyEnv()
