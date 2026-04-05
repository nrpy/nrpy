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


def register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_RHSdRHS() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5 Hamiltonian circular orbit and derivative conditions.
    The CFunction is constructed for use with gsl/multiroots.

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
Combines the SEOBNRv5 circular orbit condition and it's Jacobian to pass to gsl/multiroots.

@params x - The gsl_vector object containing the radial separation and angular momentum.
@params params - The Common data structure containing the model parameters.
@params f - The gsl_vector object to store the respective derivatives of the Hamiltonian.
@params J - The gsl_matrix object to store the Jacobian of the circular orbit conditions.
@returns - GSL_SUCCESS (0) as required by GSL.
"""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_RHSdRHS"
    params = "const gsl_vector *restrict x, void *restrict params, gsl_vector *restrict f, gsl_matrix *restrict J"  # This convention is from the gsl multi-dimensional root finder example
    body = r"""
SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit(x , params , f) ;
SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_dRHS(x , params , J) ;
return GSL_SUCCESS;
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
