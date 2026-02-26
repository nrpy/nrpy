"""
Library of C functions for solving the BSSN equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Tuple

import nrpy.c_codegen as ccg
import nrpy.finite_difference as fin
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy.equations.general_relativity import psi4
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)


def generate_CFunction_psi4_metric_deriv_quantities(
    CoordSystem: str,
    enable_fd_functions: bool,
) -> Tuple[str, str]:
    """
    Register C function for psi4 metric derivative quantity computations.

    :param CoordSystem: The coordinate system to be used.
    :param enable_fd_functions: Flag to enable or disable the finite difference functions.

    :return: Tuple containing the prefunction and prefunction launch code.
    """
    parallelization = par.parval_from_str("parallelization")
    # Initialize psi4 tetrad
    psi4_class = psi4.Psi4(
        CoordSystem,
        enable_rfm_precompute=False,
    )

    desc = "Compute metric derivative quantities gamma_{ij,kl}, Gamma^i_{jk}, and K_{ij,k} needed for psi4."
    name = "psi4_metric_deriv_quantities"
    cfunc_type = "void"
    cfunc_decorators = "__device__" if parallelization == "cuda" else ""

    arg_dict_cuda = {
        "in_gfs": "const REAL *restrict",
        "auxevol_gfs": "const REAL *restrict",
        "xx0": "const REAL",
        "xx1": "const REAL",
        "xx2": "const REAL",
        "i0": "const int",
        "i1": "const int",
        "i2": "const int",
        "arr_gammaDDdDD[81]": "REAL",
        "arr_GammaUDD[27]": "REAL",
        "arr_KDDdD[27]": "REAL",
    }

    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }

    # Find symbols stored in params
    param_symbols, _ = get_params_commondata_symbols_from_expr_list(
        psi4_class.metric_derivs_expr_list, exclude=[f"xx{j}" for j in range(3)]
    )
    loop_params = parallel_utils.get_loop_parameters(parallelization)

    params_definitions = generate_definition_header(
        param_symbols,
        var_access=parallel_utils.get_params_access(parallelization),
    )
    kernel_body = f"{loop_params}\n{params_definitions}\n"

    kernel_body += ccg.c_codegen(
        psi4_class.metric_derivs_expr_list,
        psi4_class.metric_derivs_varname_arr_list,
        verbose=False,
        enable_cse=True,
        include_braces=False,
        enable_fd_codegen=True,
        enable_fd_functions=enable_fd_functions,
    )

    prefunc = ""
    if parallelization == "cuda" and enable_fd_functions:
        prefunc = fin.construct_FD_functions_prefunc(
            cfunc_decorators="__device__ "
        ).replace("SIMD", "CUDA")
    elif enable_fd_functions:
        prefunc = fin.construct_FD_functions_prefunc()

    kernel, launch_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body.replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD"),
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=desc,
        cfunc_type=cfunc_type,
        launchblock_with_braces=False,
        cfunc_decorators=cfunc_decorators,
    )

    for arg in ["arr_gammaDDdDD[81]", "arr_GammaUDD[27]", "arr_KDDdD[27]"]:
        launch_body = launch_body.replace(arg, arg.split("[")[0])

    prefunc += kernel

    return prefunc, launch_body
