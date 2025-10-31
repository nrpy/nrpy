"""
Set up C function library for SEOBNR initial conditions.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative_nodf() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5's conservative initial conditions using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
"""
    desc = """
Evaluates the conservative initial conditions for the SEOBNRv5 ODE integration.
The conservative initial conditions are given by a circular orbit (p_{r_*} = 0) such that the 
the time derivative of the tortoise momentum is zero and the time derivative of the orbital phase equals the input orbital frequency.
@params commondata - The Common data structure containing the model parameters.
@returns - GSL_SUCCESS (0) upon success.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_initial_conditions_conservative"
    params = "commondata_struct *restrict commondata"
    body = """
const size_t n = 2;
gsl_multiroot_function f = {&SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit,
                                n,
                                commondata};
REAL omega = commondata->initial_omega;
REAL pphi = pow(omega,-1./3.);
REAL r = pphi*pphi;
const REAL x_guess[2] = {r,pphi};
REAL x_result[2] = {0.,0.};
root_finding_multidimensional(2, x_guess, &f, x_result);
commondata->r = x_result[0];
commondata->pphi = x_result[1];
"""
    cfc.register_CFunction(
        subdirectory="initial_conditions",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
