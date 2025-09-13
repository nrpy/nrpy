"""
Register C functions for evaluating RHSs and applying KO dissipation in BHaHAHA.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.bhahaha.ExpansionFunctionTheta import (
    ExpansionFunctionTheta,
)
from nrpy.infrastructures import BHaH


def register_CFunction_rhs_eval(
    CoordSystem: str = "Spherical",
    enable_rfm_precompute: bool = True,
    enable_simd: bool = True,
    enable_fd_functions: bool = True,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for evaluating RHSs.

    :param CoordSystem: Coordinate system to use.
    :param enable_rfm_precompute: Enable reference metric precomputation.
    :param enable_simd: Enable SIMD optimizations.
    :param enable_fd_functions: Enable finite difference functions.
    :param OMP_collapse: Number of loops to collapse in OpenMP.
    :return: The NRPyEnv_type environment object.

    >>> env = register_CFunction_rhs_eval()
    >>> env is not None
    True
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = "Evaluate RHSs"
    cfunc_type = "void"
    name = "rhs_eval"
    params = (
        "const commondata_struct *restrict commondata, "
        "const params_struct *restrict params, "
        "REAL *restrict xx[3], "
        "const REAL *restrict auxevol_gfs, "
        "const REAL *restrict in_gfs, "
        "REAL *restrict rhs_gfs"
    )
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )

    Th = ExpansionFunctionTheta[
        CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    ]
    eta_damping = par.register_CodeParameter(
        "REAL",
        __name__,
        name="eta_damping",
        defaultvalue=2.0,
        commondata=True,
    )

    # vv and hh are gridfunctions:
    hh_rhs = sp.Symbol("vv", real=True) - eta_damping * sp.Symbol("hh", real=True)
    # wavespeed = sp.Symbol("variable_wavespeed", real=True)
    wavespeed = sp.sympify(1)  # Disabling variable wavespeed, as it does not help.
    # The minus sign here is because F_{,ij} ~ -h_{,ij}
    vv_rhs = -(wavespeed**2) * Th.Theta

    body = BHaH.simple_loop.simple_loop(
        loop_body=ccg.c_codegen(
            [hh_rhs, vv_rhs],
            [
                gri.BHaHGridFunction.access_gf("hh", gf_array_name="rhs_gfs"),
                gri.BHaHGridFunction.access_gf("vv", gf_array_name="rhs_gfs"),
            ],
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            enable_fd_functions=enable_fd_functions,
        ),
        loop_region="interior",
        enable_intrinsics=enable_simd,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return pcg.NRPyEnv()


def register_CFunction_KO_apply(
    CoordSystem: str = "Spherical",
    enable_rfm_precompute: bool = True,
    enable_simd: bool = True,
    enable_fd_functions: bool = True,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for applying KO dissipation.

    :param CoordSystem: Coordinate system to use.
    :param enable_rfm_precompute: Enable reference metric precomputation.
    :param enable_simd: Enable SIMD optimizations.
    :param enable_fd_functions: Enable finite difference functions.
    :param OMP_collapse: Number of loops to collapse in OpenMP.
    :return: The NRPyEnv_type environment object.

    >>> env = register_CFunction_KO_apply()
    >>> env is not None
    True
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h"]
    desc = "Apply KO"
    cfunc_type = "void"
    name = "KO_apply"
    params = (
        "const commondata_struct *restrict commondata, "
        "const params_struct *restrict params, "
        "REAL *restrict xx[3], "
        "const REAL *restrict auxevol_gfs, "
        "const REAL *restrict in_gfs, "
        "REAL *restrict rhs_gfs"
    )
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    orig_fd_order = par.parval_from_str("fd_order")
    # Reduce the order of KO dissipation to fit within an NGHOSTS-sized stencil.
    par.set_parval_from_str("fd_order", orig_fd_order - 2)

    KO_diss_strength = par.register_CodeParameter(
        "REAL",
        __name__,
        "KO_diss_strength",
        0.5,
        commondata=True,
    )

    rfm = refmetric.reference_metric["Spherical_rfm_precompute"]
    hh_rhs = sp.Symbol(
        gri.BHaHGridFunction.access_gf("hh", gf_array_name="rhs_gfs"), real=True
    )
    vv_rhs = sp.Symbol(
        gri.BHaHGridFunction.access_gf("vv", gf_array_name="rhs_gfs"), real=True
    )
    hh_dKOD = ixp.declarerank1("hh_dKOD")
    vv_dKOD = ixp.declarerank1("vv_dKOD")
    for k in range(1, 3):
        hh_rhs += KO_diss_strength * hh_dKOD[k] * rfm.ReU[k]
        vv_rhs += KO_diss_strength * vv_dKOD[k] * rfm.ReU[k]

    body = "if(commondata->KO_diss_strength == 0.0) return;\n"
    h = sp.Symbol("hh", real=True)
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            [
                hh_rhs.subs(rfm.xx[0], h).subs(sp.sympify("f0_of_xx0"), h),
                vv_rhs.subs(rfm.xx[0], h).subs(sp.sympify("f0_of_xx0"), h),
            ],
            [
                gri.BHaHGridFunction.access_gf("hh", gf_array_name="rhs_gfs"),
                gri.BHaHGridFunction.access_gf("vv", gf_array_name="rhs_gfs"),
            ],
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            enable_fd_functions=enable_fd_functions,
        ),
        loop_region="interior",
        enable_intrinsics=enable_simd,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )
    par.set_parval_from_str("fd_order", orig_fd_order)

    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return pcg.NRPyEnv()
