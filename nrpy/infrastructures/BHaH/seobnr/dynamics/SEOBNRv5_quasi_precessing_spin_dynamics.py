"""
Register CFunction for evolve the SEOBNRv5 spins and angular momenta using GSL.

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
import nrpy.params as par


def register_CFunction_SEOBNRv5_quasi_precessing_spin_dynamics() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for integrating the SEOBNRv5 spin evolution equations using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # register spline parameters needed for orbital dynamics
    par.register_CodeParameters(
        "spline_data",
        __name__,
        [
            "chi1_lnhat",
            "chi2_lnhat",
            "chi1_l",
            "chi2_l",
            "lnhat_x",
            "lnhat_y",
            "lnhat_z",
            "L_x",
            "L_y",
            "L_z",
        ],
        commondata=True,
        add_to_parfile=False,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
"""
    desc = """
Integrates the SEOBNRv5 quasi-precessing spin evolution equations of motion.
Creates splines that are accessed by the orbital dynamics.

@param commondata - Common data struct containing the model parameters.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_quasi_precessing_spin_dynamics"
    params = "commondata_struct *restrict commondata"
    body = """
const gsl_odeiv2_step_type *restrict T = gsl_odeiv2_step_rk8pd;
gsl_odeiv2_step *restrict s = gsl_odeiv2_step_alloc(T, NUMVARS_SPIN);
gsl_odeiv2_control *restrict c = gsl_odeiv2_control_standard_new(1e-12, 1e-11, 1.0, 1.0);
gsl_odeiv2_system sys = {SEOBNRv5_quasi_precessing_spin_equations, NULL, NUMVARS_SPIN, commondata};
gsl_odeiv2_evolve *restrict e = gsl_odeiv2_evolve_alloc (NUMVARS_SPIN);

REAL t = 0.0;
const REAL tmax = 2e9;
REAL z[NUMVARS_SPIN], dzdt[NUMVARS_SPIN], L[3];
int status = 0;
int stop = 0;
// Initial conditions
z[LN_X] = 0.;
z[LN_Y] = 0.;
z[LN_Z] = 1.;
z[CHI1_X] = commondata->chi1_x;
z[CHI1_Y] = commondata->chi1_y;
z[CHI1_Z] = commondata->chi1_z;
z[CHI2_X] = commondata->chi2_x;
z[CHI2_Y] = commondata->chi2_y;
z[CHI2_Z] = commondata->chi2_z;
z[OMEGA_PN] = commondata->initial_omega;
SEOBNRv5_quasi_precessing_spin_angular_momentum(z, L, commondata);
int rhs_status[1] = {GSL_SUCCESS};
char rhs_name[] = "gsl_odeiv2_evolve_apply";
// set the step size according to the orbital frequency
REAL h = 2 * M_PI / commondata->initial_omega / 5;
// declare relevant pointers
size_t bufferlength = (size_t)(tmax / h); // runs up to 0.01x maximum time (we should not ideally run that long)
REAL *chi1_lnhat = (REAL *)malloc(bufferlength * sizeof(REAL));
REAL *chi2_lnhat = (REAL *)malloc(bufferlength * sizeof(REAL));
REAL *chi1_l = (REAL *)malloc(bufferlength * sizeof(REAL));
REAL *chi2_l = (REAL *)malloc(bufferlength * sizeof(REAL));
REAL *lnhat_x = (REAL *)malloc(bufferlength * sizeof(REAL));
REAL *lnhat_y = (REAL *)malloc(bufferlength * sizeof(REAL));
REAL *lnhat_z = (REAL *)malloc(bufferlength * sizeof(REAL));
REAL *L_x = (REAL *)malloc(bufferlength * sizeof(REAL));
REAL *L_y = (REAL *)malloc(bufferlength * sizeof(REAL));
REAL *L_z = (REAL *)malloc(bufferlength * sizeof(REAL));
REAL *omega = (REAL *)malloc(bufferlength * sizeof(REAL));
if (chi1_lnhat == NULL || chi2_lnhat == NULL || chi1_l == NULL || chi2_l == NULL || lnhat_x == NULL || lnhat_y == NULL || lnhat_z == NULL || L_x == NULL || L_y == NULL || L_z == NULL || omega == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_quasi_precessing_spin_dynamics(), malloc() failed for spin_dynamics\\n");
  exit(1);
}
size_t nsteps = 0;

// store
chi1_lnhat[nsteps] = z[LN_X]*z[CHI1_X] + z[LN_Y]*z[CHI1_Y] + z[LN_Z]*z[CHI1_Z];
REAL L_mod_inv = 1./sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
status = SEOBNRv5_quasi_precessing_spin_equations(t, z, dzdt, commondata);
chi2_lnhat[nsteps] = z[LN_X]*z[CHI2_X] + z[LN_Y]*z[CHI2_Y] + z[LN_Z]*z[CHI2_Z];
chi1_l[nsteps] = (z[CHI1_X]*L[0] + z[CHI1_Y]*L[1] + z[CHI1_Z]*L[2])*L_mod_inv;
chi2_l[nsteps] = (z[CHI2_X]*L[0] + z[CHI2_Y]*L[1] + z[CHI2_Z]*L[2])*L_mod_inv;
REAL ln_mod_inv = 1./sqrt(z[LN_X]*z[LN_X] + z[LN_Y]*z[LN_Y] + z[LN_Z]*z[LN_Z]);
lnhat_x[nsteps] = z[LN_X]*ln_mod_inv;
lnhat_y[nsteps] = z[LN_Y]*ln_mod_inv;
lnhat_z[nsteps] = z[LN_Z]*ln_mod_inv;
L_x[nsteps] = L[0];
L_y[nsteps] = L[1];
L_z[nsteps] = L[2];
omega[nsteps] = z[OMEGA_PN];
nsteps++;
REAL Omega_previous = commondata->initial_omega;
REAL time_previous = t;
REAL r , v, Omega;

while (t < tmax && stop == 0) {
  // integrate
  status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, tmax, &h, z);
  handle_gsl_return_status(status,rhs_status,1,rhs_name);

  // buffercheck
if (nsteps >= bufferlength) {
    bufferlength = 2 * bufferlength;
    chi1_lnhat = (REAL *)realloc(chi1_lnhat, bufferlength * sizeof(REAL));
    chi2_lnhat = (REAL *)realloc(chi2_lnhat, bufferlength * sizeof(REAL));
    chi1_l = (REAL *)realloc(chi1_l, bufferlength * sizeof(REAL));
    chi2_l = (REAL *)realloc(chi2_l, bufferlength * sizeof(REAL));
    lnhat_x = (REAL *)realloc(lnhat_x, bufferlength * sizeof(REAL));
    lnhat_y = (REAL *)realloc(lnhat_y, bufferlength * sizeof(REAL));
    lnhat_z = (REAL *)realloc(lnhat_z, bufferlength * sizeof(REAL));
    L_x = (REAL *)realloc(L_x, bufferlength * sizeof(REAL));
    L_y = (REAL *)realloc(L_y, bufferlength * sizeof(REAL));
    L_z = (REAL *)realloc(L_z, bufferlength * sizeof(REAL));
    omega = (REAL *)realloc(omega, bufferlength * sizeof(REAL));
    if (chi1_lnhat == NULL || chi2_lnhat == NULL || chi1_l == NULL || chi2_l == NULL || lnhat_x == NULL || lnhat_y == NULL || lnhat_z == NULL || L_x == NULL || L_y == NULL || L_z == NULL || omega == NULL){
      fprintf(stderr,"Error: in SEOBNRv5_quasi_precessing_spin_dynamics(), realloc() failed for spin_dynamics\\n");
      exit(1);
    }
  }

  // update
  // store
  SEOBNRv5_quasi_precessing_spin_angular_momentum(z, L, commondata);
  L_mod_inv = 1./sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
  L_x[nsteps] = L[0];
  L_y[nsteps] = L[1];
  L_z[nsteps] = L[2];
  ln_mod_inv = 1./sqrt(z[LN_X]*z[LN_X] + z[LN_Y]*z[LN_Y] + z[LN_Z]*z[LN_Z]);
  lnhat_x[nsteps] = z[LN_X]*ln_mod_inv;
  lnhat_y[nsteps] = z[LN_Y]*ln_mod_inv;
  lnhat_z[nsteps] = z[LN_Z]*ln_mod_inv;
  chi1_lnhat[nsteps] = lnhat_x[nsteps]*z[CHI1_X] + lnhat_y[nsteps]*z[CHI1_Y] + lnhat_z[nsteps]*z[CHI1_Z];
  chi2_lnhat[nsteps] = lnhat_x[nsteps]*z[CHI2_X] + lnhat_y[nsteps]*z[CHI2_Y] + lnhat_z[nsteps]*z[CHI2_Z];
  chi1_l[nsteps] = (z[CHI1_X]*L[0] + z[CHI1_Y]*L[1] + z[CHI1_Z]*L[2])*L_mod_inv;
  chi2_l[nsteps] = (z[CHI2_X]*L[0] + z[CHI2_Y]*L[1] + z[CHI2_Z]*L[2])*L_mod_inv;
  omega[nsteps] = z[OMEGA_PN];
  nsteps++;

  // stopcheck
  Omega = omega[nsteps - 1];
  r = pow(Omega,-2./3.);
  v = 1./sqrt(r);
  // superluminal velocity
  if (v > 1.0){
    break;
  }
  // decrease in frequency: Omega peak stop = index of omega
  if (Omega < Omega_previous || (Omega - Omega_previous < 1e-9 && r < 6.)) {
    stop = OMEGA_PN;
    break;
  }
  // empirical upper bound on the frequency
  if (Omega > 0.35){
    break;
  }
  Omega_previous = Omega;
  time_previous = t;
}

// free up gsl ode solver
gsl_odeiv2_control_free(c);
gsl_odeiv2_step_free(s);
gsl_odeiv2_evolve_free(e);

// remove the last point if:
// there is a peak in omega
// the last timestep is too small
REAL dt_last = t - time_previous;
if (stop == OMEGA_PN || dt_last < 1e-12){
  nsteps--;
}

// create the splines
commondata->chi1_lnhat.spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
commondata->chi1_lnhat.acc = gsl_interp_accel_alloc();
gsl_spline_init(commondata->chi1_lnhat.spline, omega, chi1_lnhat, nsteps);

commondata->chi2_lnhat.spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
commondata->chi2_lnhat.acc = gsl_interp_accel_alloc();
gsl_spline_init(commondata->chi2_lnhat.spline, omega, chi2_lnhat, nsteps);

commondata->chi1_l.spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
commondata->chi1_l.acc = gsl_interp_accel_alloc();
gsl_spline_init(commondata->chi1_l.spline, omega, chi1_l, nsteps);

commondata->chi2_l.spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
commondata->chi2_l.acc = gsl_interp_accel_alloc();
gsl_spline_init(commondata->chi2_l.spline, omega, chi2_l, nsteps);

commondata->lnhat_x.spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
commondata->lnhat_x.acc = gsl_interp_accel_alloc();
gsl_spline_init(commondata->lnhat_x.spline, omega, lnhat_x, nsteps);

commondata->lnhat_y.spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
commondata->lnhat_y.acc = gsl_interp_accel_alloc();
gsl_spline_init(commondata->lnhat_y.spline, omega, lnhat_y, nsteps);

commondata->lnhat_z.spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
commondata->lnhat_z.acc = gsl_interp_accel_alloc();
gsl_spline_init(commondata->lnhat_z.spline, omega, lnhat_z, nsteps);

commondata->L_x.spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
commondata->L_x.acc = gsl_interp_accel_alloc();
gsl_spline_init(commondata->L_x.spline, omega, L_x, nsteps);

commondata->L_y.spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
commondata->L_y.acc = gsl_interp_accel_alloc();
gsl_spline_init(commondata->L_y.spline, omega, L_y, nsteps);

commondata->L_z.spline = gsl_spline_alloc(gsl_interp_cspline, nsteps);
commondata->L_z.acc = gsl_interp_accel_alloc();
gsl_spline_init(commondata->L_z.spline, omega, L_z, nsteps);

// free up memory
free(chi1_lnhat);
free(chi2_lnhat);
free(chi1_l);
free(chi2_l);
free(lnhat_x);
free(lnhat_y);
free(lnhat_z);
free(L_x);
free(L_y);
free(L_z);
free(omega);
"""
    cfc.register_CFunction(
        subdirectory="dynamics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        prefunc=prefunc,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
