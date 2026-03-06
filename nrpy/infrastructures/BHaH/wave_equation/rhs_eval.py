"""
C function for evaluating the RHS of the wave equation in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from pathlib import Path
from types import FrameType as FT
from typing import List, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.equations.wave_equation.WaveEquationCurvilinear_RHSs import (
    WaveEquationCurvilinear_RHSs,
)
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)
from nrpy.infrastructures import BHaH


def register_CFunction_rhs_eval(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_intrinsics: bool,
    enable_KreissOliger_dissipation: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side (RHS) evaluation function for the wave equation.

    This function sets the right-hand side of the wave equation according to the
    selected coordinate system and specified parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_intrinsics: Whether to enable hardware intrinsics, e.g. SIMD.
    :param enable_KreissOliger_dissipation: Whether to enable Kreiss-Oliger dissipation, to damp high-frequency noise.
    :param OMP_collapse: Level of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h"]
    if enable_intrinsics:
        includes += [
            str(
                Path("intrinsics") / "cuda_intrinsics.h"
                if parallelization == "cuda"
                else Path("intrinsics") / "simd_intrinsics.h"
            )
        ]
    desc = r"""Set RHSs for wave equation."""
    cfunc_type = "void"
    name = "rhs_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict in_gfs, REAL *restrict rhs_gfs"
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]",
            "const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs",
        )
    # Populate uu_rhs, vv_rhs
    rhs = WaveEquationCurvilinear_RHSs(CoordSystem, enable_rfm_precompute)
    if enable_KreissOliger_dissipation:
        diss_strength = par.register_CodeParameter(
            "REAL",
            __name__,
            "KreissOliger_diss_strength",
            0.9,
            commondata=True,
        )
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        uu_dKOD = ixp.declarerank1("uu_dKOD")
        vv_dKOD = ixp.declarerank1("vv_dKOD")
        for k in range(3):
            rhs.uu_rhs += (
                diss_strength * uu_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.vv_rhs += (
                diss_strength * vv_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
    expr_list: List[sp.Expr] = [rhs.uu_rhs, rhs.vv_rhs]

    # Find symbols stored in params
    param_symbols, commondata_symbols = get_params_commondata_symbols_from_expr_list(
        expr_list, exclude=[f"xx{j}" for j in range(3)]
    )
    params_definitions = generate_definition_header(
        param_symbols,
        var_access=parallel_utils.get_params_access(parallelization),
    )

    kernel_body = parallel_utils.get_loop_parameters(
        parallelization, enable_intrinsics=enable_intrinsics
    )

    if enable_intrinsics:
        for symbol in commondata_symbols:
            kernel_body += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ConstSIMD(NOSIMD{symbol});\n"

    kernel_body += f"{params_definitions}\n"
    kernel_body += BHaH.simple_loop.simple_loop(
        loop_body=ccg.c_codegen(
            expr_list,
            [
                gri.BHaHGridFunction.access_gf("uu", gf_array_name="rhs_gfs"),
                gri.BHaHGridFunction.access_gf("vv", gf_array_name="rhs_gfs"),
            ],
            enable_fd_codegen=True,
            enable_simd=enable_intrinsics,
        ),
        loop_region="interior",
        enable_intrinsics=enable_intrinsics,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    arg_dict_cuda = {
        "in_gfs": "const REAL *restrict",
        "rhs_gfs": "REAL *restrict",
    }

    if enable_rfm_precompute:
        arg_dict_cuda = {
            "rfmstruct": "const rfm_struct *restrict",
            "auxevol_gfs": "const REAL *restrict",
            **arg_dict_cuda,
        }
    else:
        arg_dict_cuda = {
            "x0": "const REAL *restrict",
            "x1": "const REAL *restrict",
            "x2": "const REAL *restrict",
            **arg_dict_cuda,
        }

    arg_dict_cuda = {
        **arg_dict_cuda,
        **{
            f"NOCUDA{k}" if enable_intrinsics else k: "const REAL"
            for k in commondata_symbols
        },
    }

    arg_dict_host = {
        "params": "const params_struct *restrict",
        **{k.replace("CUDA", "SIMD"): v for k, v in arg_dict_cuda.items()},
    }

    prefunc, launch_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body.replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD"),
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=desc,
        cfunc_type=f"static {cfunc_type}",
        launchblock_with_braces=False,
        thread_tiling_macro_suffix="WAVE_RHS",
    )

    for symbol in commondata_symbols:
        tmp_sym = (
            f"NOCUDA{symbol}"
            if parallelization == "cuda" and enable_intrinsics
            else (
                f"NOSIMD{symbol}"
                if enable_intrinsics
                else symbol if enable_intrinsics else symbol
            )
        )
        launch_body = launch_body.replace(tmp_sym, f"commondata->{symbol}")

    cfc.register_CFunction(
        prefunc=prefunc,
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
