"""
Generate C function to set cahdprefactor gridfunction when solving the BSSN equations in curvilinear coordinates.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Set, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)


def register_CFunction_cahdprefactor_auxevol_gridfunction(
    set_of_CoordSystems: Set[str],
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Add function that sets cahdprefactor gridfunction = C_CAHD * CFL_FACTOR * dsmin to Cfunction dictionary.

    :param set_of_CoordSystems: Coordinate systems used.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = par.parval_from_str("parallelization")
    for CoordSystem in set_of_CoordSystems:
        desc = "cahdprefactor_auxevol_gridfunction(): Initialize CAHD prefactor (auxevol) gridfunction."
        name = "cahdprefactor_auxevol_gridfunction"
        params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], REAL *restrict auxevol_gfs"

        arg_dict_cuda = {
            "x0": "const REAL *restrict",
            "x1": "const REAL *restrict",
            "x2": "const REAL *restrict",
            "auxevol_gfs": "REAL *restrict",
        }

        arg_dict_host = {
            "commondata": "const commondata_struct *restrict",
            "params": "const params_struct *restrict",
            **arg_dict_cuda,
        }

        rfm = refmetric.reference_metric[CoordSystem]
        dxx0, dxx1, dxx2 = sp.symbols("dxx0 dxx1 dxx2", real=True)
        loop_body = r"""  // Compute cahdprefactor gridfunction = C_CAHD * CFL_FACTOR * dsmin.
REAL dsmin0, dsmin1, dsmin2;
"""
        expr_list = [
            sp.Abs(rfm.scalefactor_orthog[0] * dxx0),
            sp.Abs(rfm.scalefactor_orthog[1] * dxx1),
            sp.Abs(rfm.scalefactor_orthog[2] * dxx2),
        ]
        loop_body += ccg.c_codegen(
            expr_list,
            ["dsmin0", "dsmin1", "dsmin2"],
            include_braces=False,
        )
        loop_body += """auxevol_gfs[IDX4(CAHDPREFACTORGF, i0, i1, i2)] = C_CAHD * CFL_FACTOR * MIN(dsmin0, MIN(dsmin1, dsmin2));"""

        # Find symbols stored in params
        param_symbols, commondata_symbols = (
            get_params_commondata_symbols_from_expr_list(
                expr_list, exclude=[f"xx{j}" for j in range(3)]
            )
        )
        loop_params = parallel_utils.get_loop_parameters(parallelization)

        params_definitions = generate_definition_header(
            param_symbols,
            var_access=parallel_utils.get_params_access(parallelization),
        )

        commondata_definitions = generate_definition_header(
            commondata_symbols + ["C_CAHD", "CFL_FACTOR"],
            var_access=parallel_utils.get_commondata_access(parallelization),
        )

        kernel_body = f"{loop_params}\n{params_definitions}\n{commondata_definitions}\n"
        kernel_body += lp.simple_loop(
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
            cfunc_type="void",
            launchblock_with_braces=False,
            thread_tiling_macro_suffix="CAHD",
        )

        body = launch_body
        for i in range(3):
            body = body.replace(f"x{i}", f"xx[{i}]")

        cfc.register_CFunction(
            prefunc=prefunc,
            includes=["BHaH_defines.h"],
            desc=desc,
            name=name,
            params=params,
            include_CodeParameters_h=False,
            body=body,
            CoordSystem_for_wrapper_func=CoordSystem,
        )
    return pcg.NRPyEnv()
