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
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_constants as SEOBNRv5_const
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_Hamiltonian as SEOBNRv5_Ham
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_coefficients() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the masses and SEOBNRv5 coefficients.
    The inputs needed to generate an EOB approximant are mass ratio, spins and
    an initial orbital frequency. Therefore, one needs to compute the individual
    masses from the mass ratio and the Hamiltonian coeffients which are a function
    of mass ratio and spins.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """Evaluate the SEOBNRv5 coefficients."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_coefficients"
    params = "commondata_struct *restrict commondata"
    v5_const = SEOBNRv5_const.SEOBNR_aligned_spin_constants()
    body = """
const REAL q = commondata->mass_ratio;
commondata->m1 = q / (1.0 + q);
commondata->m2 = 1.0 / (1.0 + q);
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
commondata->dT = commondata->dt / commondata->total_mass / 4.925490947641266978197229498498379006e-6;
"""
    body += ccg.c_codegen(
        [
            v5_const.pyseobnr_a6,
            v5_const.pyseobnr_dSO,
            v5_const.Delta_t,
            v5_const.M_f,
            v5_const.a_f,
            v5_const.rISCO,
            v5_const.rstop,
        ],
        [
            "commondata->a6",
            "commondata->dSO",
            "commondata->Delta_t",
            "commondata->M_f",
            "commondata->a_f",
            "commondata->r_ISCO",
            "commondata->r_stop",
        ],
        verbose=False,
        include_braces=False,
    )
    body += """
const REAL afinallist[107] = { -0.9996, -0.9995, -0.9994, -0.9992, -0.999, -0.9989, -0.9988,
  -0.9987, -0.9986, -0.9985, -0.998, -0.9975, -0.997, -0.996, -0.995, -0.994, -0.992, -0.99, -0.988,
  -0.986, -0.984, -0.982, -0.98, -0.975, -0.97, -0.96, -0.95, -0.94, -0.92, -0.9, -0.88, -0.86, -0.84,
  -0.82, -0.8, -0.78, -0.76, -0.74, -0.72, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3,
  -0.25, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
  0.65, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97,
  0.975, 0.98, 0.982, 0.984, 0.986, 0.988, 0.99, 0.992, 0.994, 0.995, 0.996, 0.997, 0.9975, 0.998,
  0.9985, 0.9986, 0.9987, 0.9988, 0.9989, 0.999, 0.9992, 0.9994, 0.9995, 0.9996
  };

const REAL reomegaqnm22[107] = { 0.2915755,0.291581,0.2915866,0.2915976,0.2916086,0.2916142,0.2916197,0.2916252,
  0.2916307,0.2916362,0.2916638,0.2916915,0.2917191,0.2917744,0.2918297,0.291885,0.2919958,0.2921067,0.2922178,0.2923289,
  0.2924403,0.2925517,0.2926633,0.292943,0.2932235,0.2937871,0.2943542,0.2949249,0.2960772,0.2972442,0.2984264,0.299624,
  0.3008375,0.3020672,0.3033134,0.3045767,0.3058573,0.3071558,0.3084726,0.3098081,0.3132321,0.316784,0.3204726,0.3243073,
  0.3282986,0.3324579,0.336798,0.3413329,0.3460786,0.3510526,0.3562748,0.3617677,0.3675569,0.3736717,0.3801456,0.3870175,
  0.394333,0.4021453,0.4105179,0.4195267,0.4292637,0.4398419,0.4514022,0.464123,0.4782352,0.4940448,0.5119692,0.5326002,
  0.5417937,0.5516303,0.5622007,0.5736164,0.586017,0.5995803,0.6145391,0.631206,0.6500179,0.6716143,0.6969947,0.7278753,
  0.74632,0.7676741,0.7932082,0.8082349,0.8254294,0.8331,0.8413426,0.8502722,0.8600456,0.8708927,0.8830905,0.8969183,0.9045305
  ,0.912655,0.9213264,0.9258781,0.9305797,0.9354355,0.9364255,0.937422,0.9384248,0.9394341,0.9404498,0.9425009,0.9445784,
  0.9456271,0.9466825 };

const REAL imomegaqnm22[107] = { 0.0880269,0.0880272,0.0880274,0.088028,0.0880285,0.0880288,0.088029,0.0880293,
  0.0880296,0.0880298,0.0880311,0.0880325,0.0880338,0.0880364,0.0880391,0.0880417,0.088047,0.0880523,0.0880575,0.0880628,0.088068,
  0.0880733,0.0880785,0.0880915,0.0881045,0.0881304,0.088156,0.0881813,0.0882315,0.0882807,0.0883289,0.0883763,0.0884226,0.0884679,
  0.0885122,0.0885555,0.0885976,0.0886386,0.0886785,0.0887172,0.0888085,0.0888917,0.0889663,0.0890315,0.0890868,0.0891313,0.0891643,
  0.0891846,0.0891911,0.0891825,0.0891574,0.0891138,0.0890496,0.0889623,0.0888489,0.0887057,0.0885283,0.0883112,0.0880477,0.0877293,
  0.0873453,0.086882,0.0863212,0.0856388,0.0848021,0.0837652,0.0824618,0.0807929,0.0799908,0.0790927,0.0780817,0.0769364,0.0756296,
  0.0741258,0.072378,0.0703215,0.0678642,0.0648692,0.0611186,0.0562313,0.053149,0.0494336,0.0447904,0.0419586,0.0386302,0.0371155,
  0.0354677,0.033659,0.0316517,0.0293904,0.0268082,0.0238377,0.0221857,0.0204114,0.0185063,0.0175021,0.016462,0.015385,0.0151651,
  0.0149437,0.0147207,0.0144962,0.0142701,0.0138132,0.0133501,0.0131161,0.0128806 
  };

gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, 107);
gsl_interp_accel *acc = gsl_interp_accel_alloc();

gsl_spline_init(spline, afinallist, reomegaqnm22, 107);
commondata->omega_qnm = gsl_spline_eval(spline, commondata->a_f, acc) / commondata->M_f;

gsl_spline_init(spline, afinallist, imomegaqnm22, 107);
gsl_interp_accel_reset(acc);
commondata->tau_qnm = 1./(gsl_spline_eval(spline, commondata->a_f, acc) / commondata->M_f) ;

gsl_spline_free(spline);
gsl_interp_accel_free(acc);
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative() -> (
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
    desc = """Evaluate the SEOBNRv5 conservative initial conditions."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_initial_conditions_conservative"
    params = "commondata_struct *restrict commondata"
    body = """
const size_t n = 2;
gsl_multiroot_function_fdf f = {&SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit,
                                &SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_dRHS,
                                &SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_RHSdRHS,
                                n,
                                commondata};
REAL omega = commondata->initial_omega;
REAL pphi = pow(omega,-1./3.);
REAL r = pphi*pphi;
const REAL x_guess[2] = {r,pphi};
REAL *restrict x_result = malloc(2*sizeof(REAL));
SEOBNRv5_aligned_multidimensional_root_wrapper(f,x_guess,n,x_result);
commondata->r = x_result[0];
commondata->pphi = x_result[1];
free(x_result);
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


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
    desc = """Evaluate the SEOBNRv5 conservative initial conditions."""
    cfunc_type = "int"
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
for (i = 0; i < n; i++){
x_result[i] = gsl_vector_get(s->x , i);
}
gsl_multiroot_fsolver_free (s);
gsl_vector_free (x);
commondata->r = x_result[0];
commondata->pphi = x_result[1];
free(x_result);
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


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
    desc = """Evaluate SEOBNRv5 Hamiltonian's circular orbit conditions to compute conservative initial conditions."""
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
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_dRHS() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5 Hamiltonian circular orbit conditions' derivative.
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
    desc = """Evaluate SEOBNRv5 Hamiltonian's circular orbit conditions' derivative to compute conservative initial conditions."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_dRHS"
    params = "const gsl_vector *restrict x, void *restrict params, gsl_matrix *restrict J"  # This convention is from the gsl multi-dimensional root finder example
    Hq = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    dHreal_dr_dr_circ = Hq.dHreal_dr_dr_circ
    dHreal_dr_dpphi_circ = Hq.dHreal_dr_dpphi_circ
    dHreal_dpphi_dpphi_circ = Hq.dHreal_dpphi_dpphi_circ
    body = r"""
const REAL r = gsl_vector_get(x , 0);
const REAL pphi = gsl_vector_get(x , 1);
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL chi1 = ((commondata_struct *restrict) params)->chi1;
const REAL chi2 = ((commondata_struct *restrict) params)->chi2;
const REAL a6 = ((commondata_struct *restrict) params)->a6;
const REAL dSO = ((commondata_struct *restrict) params)->dSO;
"""
    body += ccg.c_codegen(
        [dHreal_dr_dr_circ, dHreal_dr_dpphi_circ, dHreal_dpphi_dpphi_circ],
        [
            "const REAL dHreal_dr_dr_circ",
            "const REAL dHreal_dr_dpphi_circ",
            "const REAL dHreal_dpphi_dpphi_circ",
        ],
        verbose=False,
        include_braces=False,
    )
    body += r"""
gsl_matrix_set(J , 0 , 0 , dHreal_dr_dr_circ);
gsl_matrix_set(J , 0 , 1 , dHreal_dr_dpphi_circ);
gsl_matrix_set(J , 1 , 0 , dHreal_dr_dpphi_circ);
gsl_matrix_set(J , 1 , 1 , dHreal_dpphi_dpphi_circ);
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


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
    desc = """Evaluate SEOBNRv5 Hamiltonian's circular orbit and derivative conditions to compute conservative initial conditions."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_RHSdRHS"
    params = "const gsl_vector *restrict x, void *restrict params, gsl_vector *restrict f, gsl_matrix *restrict J"  # This convention is from the gsl multi-dimensional root finder example
    body = r"""
SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit(x , params , f) ;
SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_dRHS(x , params , J) ;
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


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
    desc = """Evaluate the SEOBNRv5 dissipative initial conditions."""
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
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_SEOBNRv5_aligned_spin_radial_momentum_condition() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5 initial radial tortoise momentum.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Evaluate SEOBNRv5 radial momentum condition."""
    cfunc_type = "REAL"
    name = "SEOBNRv5_aligned_spin_radial_momentum_conditions"
    params = "REAL x, void *restrict params"
    body = """
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL chi1 = ((commondata_struct *restrict) params)->chi1;
const REAL chi2 = ((commondata_struct *restrict) params)->chi2;
const REAL a6 = ((commondata_struct *restrict) params)->a6;
const REAL dSO = ((commondata_struct *restrict) params)->dSO; 
const REAL r = ((commondata_struct *restrict) params)->r; 
const REAL pphi = ((commondata_struct *restrict) params)->pphi; 
REAL prstar = x; 
"""
    H = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    body += ccg.c_codegen(
        [
            H.dHreal_dr_dr,
            H.dHreal_dr_dpphi,
            H.xi,
            H.dHreal_dprstar,
            H.Hreal,
            H.dHreal_dpphi,
            H.dHreal_dpphi_circ,
        ],
        [
            "const REAL dHreal_dr_dr",
            "const REAL dHreal_dr_dpphi",
            "const REAL xi",
            "const REAL dHreal_dprstar",
            "const REAL Hreal",
            "const REAL Omega",
            "const REAL Omega_circ",
        ],
        verbose=False,
        include_braces=False,
    )
    body += """
const REAL dLdr = - dHreal_dr_dr / dHreal_dr_dpphi;
const REAL rdot_dyn = xi * dHreal_dprstar;
REAL flux[2];
const REAL y[4] = {r , 0. , prstar , pphi};
SEOBNRv5_aligned_spin_flux(y,Hreal,Omega,Omega_circ,flux,params);
const REAL dLdt = flux[1];
const REAL rdot_rad = dLdt / dLdr;
const REAL prstar_condition = rdot_dyn - rdot_rad;    
return prstar_condition;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
