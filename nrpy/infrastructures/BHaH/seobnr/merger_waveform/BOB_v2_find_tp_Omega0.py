"""
Set up C function library for BOBv2 peak news time and reference orbital frequency.

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
import nrpy.helpers.parallel_codegen as pcg
from nrpy.equations.seobnr.BOB_v2_waveform_quantities import BOB_v2_waveform_quantities


def register_CFunction_BOB_v2_find_tp_Omega0() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for finding the BOBv2 peak news time and reference orbital frequency using GSL.

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
Evaluates the BOBv2 peak news time and reference orbital frequency using GSL.

@params commondata - The Common data structure containing the model parameters.
@returns - GSL_SUCCESS (0) upon success.
"""
    cfunc_type = "int"
    name = "BOB_v2_find_tp_Omega0"
    params = "commondata_struct *restrict commondata"
    wf = BOB_v2_waveform_quantities()
    initial_guess_code = (
        ccg.c_codegen(
            [wf.t_p_guess, wf.Omega_0_guess],
            ["const REAL t_p_guess", "const REAL Omega_0_guess"],
            verbose=False,
            include_braces=False,
        )
        .replace("REAL", "double complex")
        .replace("exp", "cexp")
        .replace("sqrt", "csqrt")
        .replace("pow", "cpow")
        .replace("fabs", "cabs")
        .replace("tanh", "ctanh")
        .replace("sinh", "csinh")
        .replace("cosh", "ccosh")
        .replace("actanh", "catanh")
    )
    body = """
const size_t n = 2;
gsl_multiroot_function f = {&BOB_v2_peak_strain_conditions,
                                n,
                                commondata};
// Compute the initial guess for the peak news time and reference orbital frequency
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL t_0 = commondata->t_attach;
const REAL omega_qnm = commondata->omega_qnm;
const REAL tau_qnm = commondata->tau_qnm;
"""
    body += initial_guess_code
    body += """
const REAL x_guess[2] = {creal(t_p_guess),creal(Omega_0_guess)};
REAL *restrict x_result = malloc(2*sizeof(REAL));
const gsl_multiroot_fsolver_type *restrict T = gsl_multiroot_fsolver_hybrids;
gsl_multiroot_fsolver *restrict s = gsl_multiroot_fsolver_alloc(T , n);
gsl_vector *restrict x = gsl_vector_alloc(n);
size_t i , iter = 0;
int status;
const int maxiter = 100;
for (i = 0; i < n; i++){
gsl_vector_set(x , i , x_guess[i]);
}
gsl_multiroot_fsolver_set(s , &f , x);
do {
  iter++;
  status = gsl_multiroot_fsolver_iterate (s);
  int f_solver_status[1] = {GSL_SUCCESS};
  char fsolver_name[] = "gsl_multiroot_fsolver_iterate";
  handle_gsl_return_status(status,f_solver_status,1,fsolver_name);
  status = gsl_multiroot_test_residual (s->f, 6e-12);
  int test_residual_status[2] = {GSL_SUCCESS,GSL_CONTINUE};
  char residual_name[] = "gsl_multiroot_test_residual";
  handle_gsl_return_status(status,test_residual_status,2,residual_name);
}
while(status == GSL_CONTINUE && iter < maxiter);
// if the peak conditions are not met, exit out.
if (status == GSL_CONTINUE) {
  const REAL t_p_error = gsl_vector_get(s->f,0);
  const REAL Omega_0_error = gsl_vector_get(s->f,1);
  const REAL total_error = sqrt(t_p_error*t_p_error + Omega_0_error*Omega_0_error);
  printf("In function BOB_v2_find_tp_Omega0, the peak strain conditions were not solved within maximum iterations.\\n");
  printf("t_p_error = %.5e\\nOmega_0_error = %.5e\\ntotal_error = %.5e\\ntolerance = %.5e\\n",t_p_error,Omega_0_error,total_error,6.e-12);
  exit(EXIT_FAILURE);
}
for (i = 0; i < n; i++){
x_result[i] = gsl_vector_get(s->x , i);
}
gsl_multiroot_fsolver_free (s);
gsl_vector_free (x);
commondata->t_p_BOB = x_result[0];
commondata->Omega_0_BOB = x_result[1];
free(x_result);
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        subdirectory="merger_waveform",
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
