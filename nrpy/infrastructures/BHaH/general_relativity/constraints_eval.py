"""
Generate C functions for computing the BSSN constraint equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.helpers.parallel_codegen as pcg
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)
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

    includes = ["BHaH_defines.h", "diagnostics/diagnostic_gfs.h"]
    desc = r"""Evaluate BSSN constraints."""
    cfunc_type = "void"
    name = "constraints_eval"
    params = "const params_struct *restrict params, const REAL *restrict xx[], const REAL *restrict in_gfs, REAL *restrict diagnostic_gfs"
    Bcon = BSSN_constraints[CoordSystem + ("_T4munu" if enable_T4munu else "")]
    body = """
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;
"""
    expr_list = [Bcon.H, Bcon.Msquared]
    param_symbols, _ = get_params_commondata_symbols_from_expr_list(expr_list)
    params_definitions = (
        generate_definition_header(
            param_symbols,
            enable_intrinsics=False,
            var_access="params->",
        )
        + "\n"
    )
    body += params_definitions
    body += BHaH.simple_loop.simple_loop(
        loop_body=ccg.c_codegen(
            expr_list,
            [
                "diagnostic_gfs[IDX4(DIAG_HAMILTONIAN, i0, i1, i2)]",
                "diagnostic_gfs[IDX4(DIAG_MSQUARED, i0, i1, i2)]",
            ],
            enable_fd_codegen=True,
            enable_simd=False,
            enable_fd_functions=enable_fd_functions,
            rational_const_alias="static const",
        ),
        loop_region="interior",
        enable_intrinsics=False,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=False,
        read_xxs=True,
        OMP_collapse=OMP_collapse,
    )
    prefunc = ""
    if enable_fd_functions:
        prefunc = fin.construct_FD_functions_prefunc()
    cfc.register_CFunction(
        subdirectory="diagnostics",
        enable_simd=False,
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        CoordSystem_for_wrapper_func=CoordSystem,
    )
    return pcg.NRPyEnv()
