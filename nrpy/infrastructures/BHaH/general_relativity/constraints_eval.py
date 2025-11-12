"""
Generate C functions for computing the BSSN constraint equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import re
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
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
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const rfm_struct *restrict rfmstruct, const REAL *restrict in_gfs, const REAL *restrict auxevol_gfs, REAL *restrict diagnostic_gfs"""
    Bcon = BSSN_constraints[
        CoordSystem
        + "_rfm_precompute_RbarDD_gridfunctions"
        + ("_T4munu" if enable_T4munu else "")
    ]
    expr_list = [Bcon.H, Bcon.Msquared]
    prefunc = """
// Map gridfunction enums to the correct backing storage on CPU vs CUDA.
// CUDA: diags entirely on CPU; RbarDD output to diagnostic_gfs: 
//      diagnostic_gfs[IDX4(DIAG_RBARDD00GF, i0, i1, i2)]
// CPU: auxevol_gfs[IDX4(RBARDD00GF, i0, i1, i2)]
#ifdef __CUDACC__
  #define GF_IN(gf, i0, i1, i2) diagnostic_gfs[IDX4(DIAG_##gf, i0, i1, i2)]
#else
  #define GF_IN(gf, i0, i1, i2) auxevol_gfs[IDX4(gf, i0, i1, i2)]
#endif // __CUDACC__
// These all resolve through GF_IN so call sites stay uniform.
#define RBARDD_GF GF_IN
#define T4UU_GF   GF_IN

"""
    loop_body = ccg.c_codegen(
        expr_list,
        [
            "diagnostic_gfs[IDX4(DIAG_HAMILTONIANGF, i0, i1, i2)]",
            "diagnostic_gfs[IDX4(DIAG_MSQUAREDGF, i0, i1, i2)]",
        ],
        enable_fd_codegen=True,
        enable_simd=True,
        enable_fd_functions=enable_fd_functions,
        rational_const_alias="static const",
    )
    # RBARDD: auxevol_gfs[IDX4(RBARDDxyGF, i0, i1, i2)] -> RBARDD_GF(RBARDDxyGF, i0, i1, i2)
    loop_body = re.sub(
        r"auxevol_gfs\[IDX4\((RBARDD[0-9][0-9]GF), i0, i1, i2\)\]",
        r"RBARDD_GF(\1, i0, i1, i2)",
        loop_body,
    )
    # T4UU: auxevol_gfs[IDX4(T4UUxyGF, i0, i1, i2)] -> T4UU_GF(T4UUxyGF, i0, i1, i2)
    loop_body = re.sub(
        r"auxevol_gfs\[IDX4\((T4UU[0-9][0-9]GF), i0, i1, i2\)\]",
        r"T4UU_GF(\1, i0, i1, i2)",
        loop_body,
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
    par.set_parval_from_str("parallelization", orig_parallelization)
    cfc.register_CFunction(
        subdirectory="diagnostics",
        enable_simd=True,
        includes=includes,
        prefunc=(
            prefunc + fin.construct_FD_functions_prefunc()
            if enable_fd_functions
            else ""
        ),
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
        CoordSystem_for_wrapper_func=CoordSystem,
    )
    return pcg.NRPyEnv()
