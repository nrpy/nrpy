"""
Generate C functions for computing the BSSN constraint equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.infrastructures import BHaH


def register_CFunction_constraints_eval(
    CoordSystem: str,
    enable_T4munu: bool,
    enable_fd_functions: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the BSSN constraints evaluation function.

    :param CoordSystem: The coordinate system to be used.
    :param enable_T4munu: Whether to enable T4munu (stress-energy terms).
    :param enable_fd_functions: Whether to enable finite difference functions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        # register_CFunction_Ricci_onept()
        print("registering0!", f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = [
        "BHaH_defines.h",
        "diagnostics/diagnostic_gfs.h",
        "intrinsics/simd_intrinsics.h",
    ]

    orig_parallelization = par.parval_from_str("parallelization")
    par.set_parval_from_str("parallelization", "openmp")
    desc = r"""Evaluate BSSN constraints."""
    cfunc_type = "void"
    name = "constraints_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const rfm_struct *restrict rfmstruct, const REAL *restrict in_gfs, REAL *restrict diagnostic_gfs"
    Bcon = BSSN_constraints[
        CoordSystem
        + "_rfm_precompute_RbarDD_gridfunctions"
        + ("_T4munu" if enable_T4munu else "")
    ]
    # The above dictionary instantiation just populated gri.glb_gridfcs_dict
    #  with RbarDD gridfunctions. As such, c_codegen will recognize them
    #  as needing to be read from memory in the below codegen. We want to
    #  1. make codegen of diagnostics fast (MUST enable RbarDD_gridfunctions)
    #  2. avoid allocating additional diagnostics gridfunctions
    #  So to these ends, this function will call another that sets RbarDD
    #  *pointwise* (or more accurately, SIMD-wise).
    if "RbarDD00" in gri.glb_gridfcs_dict:
        for i in range(3):
            for j in range(i, 3):
                gri.glb_gridfcs_dict.pop(f"RbarDD{i}{j}")
    expr_list = [Bcon.H, Bcon.Msquared]
    loop_body = ccg.c_codegen(
        expr_list,
        [
            "diagnostic_gfs[IDX4(DIAG_HAMILTONIAN, i0, i1, i2)]",
            "diagnostic_gfs[IDX4(DIAG_MSQUARED, i0, i1, i2)]",
        ],
        enable_fd_codegen=True,
        enable_simd=True,
        enable_fd_functions=enable_fd_functions,
        rational_const_alias="static const",
    )
    body = BHaH.simple_loop.simple_loop(
        loop_body=loop_body,
        loop_region="interior",
        enable_intrinsics=True,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=True,
        read_xxs=False,
        OMP_collapse=OMP_collapse,
    )
    _ = gri.register_gridfunctions_for_single_rank2(
        "RbarDD",
        symmetry="sym01",
        group="AUXEVOL",
        gf_array_name="auxevol_gfs",
    )
    prefunc = ""
    if enable_fd_functions:
        prefunc += fin.construct_FD_functions_prefunc()
    par.set_parval_from_str("parallelization", orig_parallelization)
    cfc.register_CFunction(
        subdirectory="diagnostics",
        enable_simd=True,
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
        CoordSystem_for_wrapper_func=CoordSystem,
    )
    return pcg.NRPyEnv()
