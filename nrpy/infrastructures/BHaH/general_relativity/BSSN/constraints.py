"""
Generate C functions for computing the BSSN constraint equations in curvilinear coordinates, using a reference-metric formalism.

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
import nrpy.infrastructures.BHaH.simple_loop as lp
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints


def register_CFunction_constraints(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_RbarDD_gridfunctions: bool,
    enable_T4munu: bool,
    enable_simd: bool,
    enable_fd_functions: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the BSSN constraints evaluation function.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_RbarDD_gridfunctions: Whether to enable RbarDD gridfunctions.
    :param enable_T4munu: Whether to enable T4munu (stress-energy terms).
    :param enable_simd: Whether to enable SIMD instructions.
    :param enable_fd_functions: Whether to enable finite difference functions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    Bcon = BSSN_constraints[
        CoordSystem
        + ("_rfm_precompute" if enable_rfm_precompute else "")
        + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
        + ("_T4munu" if enable_T4munu else "")
    ]

    includes = ["BHaH_defines.h"]
    if enable_simd:
        includes += [str(Path("intrinsics") / "simd_intrinsics.h")]
    desc = r"""Evaluate BSSN constraints."""
    cfunc_type = "void"
    name = "constraints_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict in_gfs, const REAL *restrict auxevol_gfs, REAL *restrict diagnostic_output_gfs"
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    Constraints_access_gfs: List[str] = []
    for var in ["H", "MSQUARED"]:
        Constraints_access_gfs += [
            gri.BHaHGridFunction.access_gf(
                var, 0, 0, 0, gf_array_name="diagnostic_output_gfs"
            )
        ]
    body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            [Bcon.H, Bcon.Msquared],
            Constraints_access_gfs,
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
        prefunc=fin.construct_FD_functions_prefunc() if enable_fd_functions else "",
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
        enable_simd=enable_simd,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
