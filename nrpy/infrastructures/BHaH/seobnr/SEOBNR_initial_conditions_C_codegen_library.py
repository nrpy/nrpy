"""
Set up C function library for SEOBNR inspiral integrations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_Hamiltonian as SEOBNRv5_Ham
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_waveform_quantities as SEOBNRv5_wf
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_coefficients() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the masses and Hamiltonian coefficients.
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
    desc = """Evaluate the SEOBNRv5 Hamiltonian coefficients."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_Hamiltonian_coefficients"
    params = "commondata_struct *restrict commondata"
    H = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    body = """
const REAL q = commondata->mass_ratio;
commondata->m1 = q / (1.0 + q);
commondata->m2 = 1.0 / (1.0 + q);
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
"""
    body += ccg.c_codegen(
        [H.pyseobnr_a6, H.pyseobnr_dSO],
        ["commondata->a6", "commondata->dSO"],
        verbose=False,
        include_braces=False,
    )
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
const gsl_multiroot_fdfsolver_type *restrict T;
gsl_multiroot_fdfsolver *restrict s;
int status;
size_t i , iter = 0;
const size_t n = 2;
const int maxiter = 100;
gsl_multiroot_function_fdf f = {&SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit,
                                &SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_dRHS,
                                &SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit_RHSdRHS,
                                n,
                                commondata};
gsl_vector *restrict x = gsl_vector_alloc(n);
REAL omega = commondata->initial_omega;
REAL pphi = pow(omega,-1./3.);
REAL r = pphi*pphi;
gsl_vector_set(x , 0 , r);
gsl_vector_set(x , 1 , pphi);
T = gsl_multiroot_fdfsolver_hybridsj;
s = gsl_multiroot_fdfsolver_alloc(T , n);
gsl_multiroot_fdfsolver_set(s , &f , x);
do {
  iter++;
  status = gsl_multiroot_fdfsolver_iterate (s);
  if (status){
    break;
  }
  status = gsl_multiroot_test_residual (s->f, 6e-12);
}
while(status == GSL_CONTINUE && iter < maxiter);
commondata->r = gsl_vector_get(s->x , 0);
commondata->pphi = gsl_vector_get(s->x, 1);
gsl_multiroot_fdfsolver_free (s);
gsl_vector_free (x);
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
  prstar = gsl_root_fsolver_root(s);
  x_lo = gsl_root_fsolver_x_lower(s);
  x_hi = gsl_root_fsolver_x_upper(s);
  status = gsl_root_test_interval(x_lo, x_hi, xtol, rtol);

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

    includes = ["BHaH_defines.h"]
    desc = """Evaluate SEOBNRv5 radial momentum condition."""
    cfunc_type = "REAL"
    name = "SEOBNRv5_aligned_spin_radial_momentum_conditions"
    params = "REAL x, void *restrict params"
    body = """
REAL m1 = ((commondata_struct *restrict) params)->m1;
REAL m2 = ((commondata_struct *restrict) params)->m2;
REAL chi1 = ((commondata_struct *restrict) params)->chi1;
REAL chi2 = ((commondata_struct *restrict) params)->chi2;
REAL a6 = ((commondata_struct *restrict) params)->a6;
REAL dSO = ((commondata_struct *restrict) params)->dSO; 
REAL r = ((commondata_struct *restrict) params)->r; 
REAL pphi = ((commondata_struct *restrict) params)->pphi; 
REAL prstar = x; 
"""
    H = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    wf = SEOBNRv5_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    flux = wf.flux()
    flux = (
        flux.subs(wf.Hreal, H.Hreal)
        .subs(wf.Omega, H.dHreal_dpphi)
        .subs(wf.Omega_circ, H.dHreal_dpphi_circ)
    )
    dLdr = -H.dHreal_dr_dr / H.dHreal_dr_dpphi
    dLdt = flux / H.nu
    rdot_rad = dLdt / dLdr
    rdot_dyn = H.xi * H.dHreal_dprstar
    body += ccg.c_codegen(
        [rdot_dyn - rdot_rad],
        ["REAL prstar_condition"],
        verbose=False,
        include_braces=False,
    )
    body += r"""
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
