"""
Generate C functions to compute the Ricci tensor in curvilinear coordinates.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from pathlib import Path
from types import FrameType as FT
from typing import List, Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)
from nrpy.infrastructures import BHaH


def register_CFunction_Ricci_eval(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_intrinsics: bool,
    enable_fd_functions: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the Ricci evaluation function.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_intrinsics: Whether to enable SIMD instructions.
    :param enable_fd_functions: Whether to enable finite difference functions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = par.parval_from_str("parallelization")
    Bq = BSSN_quantities[
        (CoordSystem if CoordSystem != "Fisheye" else "Cartesian")
        + ("_rfm_precompute" if enable_rfm_precompute else "")
    ]

    includes = ["BHaH_defines.h"]
    if enable_intrinsics:
        includes += [
            str(
                Path("intrinsics") / "cuda_intrinsics.h"
                if parallelization == "cuda"
                else Path("intrinsics") / "simd_intrinsics.h"
            )
        ]
    desc = r"""Set Ricci tensor."""
    cfunc_type = "void"
    name = "Ricci_eval"
    arg_dict_cuda = {
        "in_gfs": "const REAL *restrict",
        "auxevol_gfs": "REAL *restrict",
    }
    if enable_rfm_precompute:
        arg_dict_cuda = {
            "rfmstruct": "const rfm_struct *restrict",
            **arg_dict_cuda,
        }
    else:
        arg_dict_cuda = {
            "x0": "const REAL *restrict",
            "x1": "const REAL *restrict",
            "x2": "const REAL *restrict",
            **arg_dict_cuda,
        }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }
    params = ",".join([f"{v} {k}" for k, v in arg_dict_host.items()])

    # Populate Ricci tensor
    Ricci_access_gfs: List[str] = []
    for var in Bq.Ricci_varnames:
        Ricci_access_gfs += [
            gri.BHaHGridFunction.access_gf(var, 0, 0, 0, gf_array_name="auxevol_gfs")
        ]
    kernel_body = BHaH.simple_loop.simple_loop(
        loop_body=ccg.c_codegen(
            Bq.Ricci_exprs,
            Ricci_access_gfs,
            enable_fd_codegen=True,
            enable_simd=enable_intrinsics,
            enable_fd_functions=enable_fd_functions,
            rational_const_alias=(
                "static constexpr" if parallelization == "cuda" else "static const"
            ),
        ).replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD"),
        loop_region="interior",
        enable_intrinsics=enable_intrinsics,
        CoordSystem=(CoordSystem if CoordSystem != "Fisheye" else "Cartesian"),
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )
    loop_params = parallel_utils.get_loop_parameters(
        parallelization, enable_intrinsics=enable_intrinsics
    )

    param_symbols, _ = get_params_commondata_symbols_from_expr_list(Bq.Ricci_exprs)
    params_definitions = params_definitions = generate_definition_header(
        param_symbols,
        enable_intrinsics=enable_intrinsics,
        var_access=parallel_utils.get_params_access(parallelization),
    )
    kernel_body = f"{loop_params}\n{params_definitions}\n{kernel_body}"

    kernel, launch_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body.replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD"),
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=desc,
        cfunc_type=cfunc_type,
        launchblock_with_braces=False,
        thread_tiling_macro_suffix="RICCI_EVAL",
    )

    prefunc = ""
    if parallelization == "cuda" and enable_fd_functions:
        prefunc = fin.construct_FD_functions_prefunc(
            cfunc_decorators="__device__ "
        ).replace("SIMD", "CUDA")
    elif enable_fd_functions:
        prefunc = fin.construct_FD_functions_prefunc()

    cfc.register_CFunction(
        include_CodeParameters_h=False,
        prefunc=prefunc + kernel,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=launch_body,
        enable_simd=enable_intrinsics,
    )
    return pcg.NRPyEnv()
