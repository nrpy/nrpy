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

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_Hamiltonian as SEOBNRv5_Ham
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5 Hamiltonian circular orbit conditions.
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
Computes the derivative of the SEOBNRv5 Hamiltonian with respect to the radial separation and angular momentum.
The above values are used to evaluate the conservative initial conditions given by a circular orbit at input orbital frequency.

@params x - The gsl_vector object containing the radial separation and angular momentum.
@params params - The Common data structure containing the model parameters.
@params f - The gsl_vector object to store the respective derivatives of the Hamiltonian.
@returns - GSL_SUCCESS (0) as required by GSL.
"""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit"
    params = "const gsl_vector *restrict x , void *restrict params , gsl_vector *restrict f"  # This convention is from the gsl multi-dimensional root finder example
    Hq = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    dHreal_dr_circ = Hq.dHreal_dr_circ
    dHreal_dpphi_circ = Hq.dHreal_dpphi_circ
    body = r"""
const REAL r = gsl_vector_get(x , 0);
const REAL pphi = gsl_vector_get(x , 1);
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL chi1 = ((commondata_struct *restrict) params)->chi1;
const REAL chi2 = ((commondata_struct *restrict) params)->chi2;
const REAL a6 = ((commondata_struct *restrict) params)->a6;
const REAL dSO = ((commondata_struct *restrict) params)->dSO;
const REAL initial_omega = ((commondata_struct *restrict) params)->initial_omega;
"""
    body += ccg.c_codegen(
        [
            dHreal_dr_circ,
            dHreal_dpphi_circ,
        ],
        [
            "const REAL dHreal_dr_circ",
            "const REAL dHreal_dpphi_circ",
        ],
        verbose=False,
        include_braces=False,
    )
    body += r"""
gsl_vector_set(f , 0 , dHreal_dr_circ);
gsl_vector_set(f , 1 , dHreal_dpphi_circ - initial_omega);
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
