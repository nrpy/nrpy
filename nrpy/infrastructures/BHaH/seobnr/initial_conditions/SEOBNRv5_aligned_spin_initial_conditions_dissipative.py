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


def register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_dissipative() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5's dissipative initial conditions using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluates the dissipative initial conditions for the SEOBNRv5 ODE integration.

@params commondata - The Common data structure containing the model parameters.
@returns - GSL_SUCCESS (0) upon success.
"""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_initial_conditions_dissipative"
    params = "commondata_struct *restrict commondata"
    body = """
int status;
int iter = 0;
const int max_iter = 100;
const gsl_root_fsolver_type *T;
gsl_root_fsolver *s;
REAL prstar;
REAL x_lo = -3e-2;
REAL x_hi = 0.0;
REAL xtol = 1e-12;
REAL rtol = 1e-10;
gsl_function F;
F.function = &SEOBNRv5_aligned_spin_radial_momentum_conditions;
F.params = commondata;

T = gsl_root_fsolver_brent;
s = gsl_root_fsolver_alloc(T);
gsl_root_fsolver_set(s, &F, x_lo, x_hi);

do {
  iter++;
  status = gsl_root_fsolver_iterate(s);
  int fsolver_status[1] = {GSL_SUCCESS};
  char fsolver_name[] = "gsl_root_fsolver_iterate";
  handle_gsl_return_status(status,fsolver_status,1,fsolver_name);
  prstar = gsl_root_fsolver_root(s);

  x_lo = gsl_root_fsolver_x_lower(s);
  x_hi = gsl_root_fsolver_x_upper(s);
  status = gsl_root_test_interval(x_lo, x_hi, xtol, rtol);
  int test_interval_status[2] = {GSL_SUCCESS,GSL_CONTINUE};
  char root_test_name[] = "gsl_root_test_interval";
  handle_gsl_return_status(status,test_interval_status,2,root_test_name);


}
while (status == GSL_CONTINUE && iter < max_iter);

gsl_root_fsolver_free (s);
commondata->prstar = prstar;
return status;
"""
    cfc.register_CFunction(
        subdirectory="initial_conditions",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
