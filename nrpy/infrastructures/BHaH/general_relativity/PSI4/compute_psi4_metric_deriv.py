"""
Library of C functions for solving the BSSN equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.helpers.parallel_codegen as pcg
from nrpy.equations.general_relativity import psi4

def register_CFunction_psi4_metric_deriv_quantities(
    CoordSystem: str,
    enable_fd_functions: bool,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for psi4 metric derivative quantity computations.

    :param CoordSystem: The coordinate system to be used.
    :param enable_fd_functions: Flag to enable or disable the finite difference functions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Initialize psi4 tetrad
    psi4_class = psi4.Psi4(
        CoordSystem,
        enable_rfm_precompute=False,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Compute metric derivative quantities gamma_{ij,kl}, Gamma^i_{jk}, and K_{ij,k} needed for psi4."
    name = "psi4_metric_deriv_quantities"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL *restrict in_gfs, const REAL xx0, const REAL xx1, const REAL xx2, const int i0, const int i1, const int i2, REAL arr_gammaDDdDD[81], REAL arr_GammaUDD[27], REAL arr_KDDdD[27]"""

    body = ccg.c_codegen(
        psi4_class.metric_derivs_expr_list,
        psi4_class.metric_derivs_varname_arr_list,
        verbose=False,
        enable_cse=True,
        include_braces=False,
        enable_fd_codegen=True,
        enable_fd_functions=enable_fd_functions,
    )

    cfc.register_CFunction(
        includes=includes,
        prefunc=fin.construct_FD_functions_prefunc() if enable_fd_functions else "",
        desc=desc,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
