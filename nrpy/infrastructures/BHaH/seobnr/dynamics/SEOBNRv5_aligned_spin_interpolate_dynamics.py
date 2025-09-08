"""
Register CFunction for interpolating the SEOBNRv5 dynamics.

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


def register_CFunction_SEOBNRv5_aligned_spin_interpolate_dynamics() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for interpolating the SEOBNRv5 dynamics.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Interpolates and fine samples the SEOBNRv5 dynamics close to the end point for Non Quasi-Circular (NQC) corrections.

@param commondata - Struct containing the common data for the SEOBNRv5 dynamics.
@param dynamics_fine_prelim - Array containing the portion of the dynamics that needs fine sampling.
@param nsteps_fine_prelim - Lenght of dynamics_fine_prelim.
@param t_peak - The estimated end time of the trajectory (can be the last element or the orbital frequency/tortoise momentum peak).
@param stop - Stop condition of the ODE integration; used to determine if iterative refinement was performed.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_interpolate_dynamics"
    params = "commondata_struct *restrict commondata, REAL *restrict dynamics_fine_prelim, const size_t nsteps_fine_prelim, const REAL t_peak, const int stop"
    body = """
int i;
// Intepolate the high sampled dynamics for NQCs.
REAL time_start = dynamics_fine_prelim[IDX(0,TIME)];
REAL time_end = dynamics_fine_prelim[IDX(nsteps_fine_prelim-1,TIME)];
if (stop != 0){
  time_start = MAX(t_peak - commondata->t_stepback,dynamics_fine_prelim[IDX(0,TIME)]);
  time_end = MIN(t_peak , time_end);
}

REAL *restrict ts = (REAL *)malloc(nsteps_fine_prelim*sizeof(REAL));
REAL *restrict rs = (REAL *)malloc(nsteps_fine_prelim*sizeof(REAL));
REAL *restrict phis = (REAL *)malloc(nsteps_fine_prelim*sizeof(REAL));
REAL *restrict prs = (REAL *)malloc(nsteps_fine_prelim*sizeof(REAL));
REAL *restrict pphis = (REAL *)malloc(nsteps_fine_prelim*sizeof(REAL));
for (i = 0; i < nsteps_fine_prelim; i++) {
  ts[i] = dynamics_fine_prelim[IDX(i,TIME)];
  rs[i] = dynamics_fine_prelim[IDX(i,R)];
  phis[i] = dynamics_fine_prelim[IDX(i,PHI)];
  prs[i] = dynamics_fine_prelim[IDX(i,PRSTAR)];
  pphis[i] = dynamics_fine_prelim[IDX(i,PPHI)];
}
gsl_interp_accel *restrict r_acc = gsl_interp_accel_alloc();
gsl_spline *restrict r_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
gsl_spline_init(r_spline, ts, rs, nsteps_fine_prelim);
gsl_interp_accel *restrict phi_acc = gsl_interp_accel_alloc();
gsl_spline *restrict phi_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
gsl_spline_init(phi_spline, ts, phis, nsteps_fine_prelim);
gsl_interp_accel *restrict pr_acc = gsl_interp_accel_alloc();
gsl_spline *restrict pr_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
gsl_spline_init(pr_spline, ts, prs, nsteps_fine_prelim);
gsl_interp_accel *restrict pphi_acc = gsl_interp_accel_alloc();
gsl_spline *restrict pphi_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
gsl_spline_init(pphi_spline, ts, pphis, nsteps_fine_prelim);

const REAL dt = 0.1;
commondata->nsteps_fine = (size_t)((time_end - time_start) / dt + 1);
commondata->dynamics_fine = (REAL *)malloc(NUMVARS * commondata->nsteps_fine * sizeof(REAL));
REAL t;
for (i = 0; i < commondata->nsteps_fine; i++) {
  t = time_start + i * dt;
  commondata->dynamics_fine[IDX(i , TIME)] = t;
  commondata->dynamics_fine[IDX(i , R)] = gsl_spline_eval(r_spline, t, r_acc);
  commondata->dynamics_fine[IDX(i , PHI)] = gsl_spline_eval(phi_spline, t, phi_acc);
  commondata->dynamics_fine[IDX(i , PRSTAR)] = gsl_spline_eval(pr_spline, t, pr_acc);
  commondata->dynamics_fine[IDX(i , PPHI)] = gsl_spline_eval(pphi_spline, t, pphi_acc);
  commondata->r = commondata->dynamics_fine[IDX(i , R)];
  commondata->phi = commondata->dynamics_fine[IDX(i , PHI)];
  commondata->prstar = commondata->dynamics_fine[IDX(i , PRSTAR)];
  commondata->pphi = commondata->dynamics_fine[IDX(i , PPHI)];
  SEOBNRv5_aligned_spin_augments(commondata);
  commondata->dynamics_fine[IDX(i , H)] = commondata->Hreal;
  commondata->dynamics_fine[IDX(i , OMEGA)] = commondata->dHreal_dpphi;
  commondata->dynamics_fine[IDX(i , OMEGA_CIRC)] = commondata->Omega_circ;
}

gsl_spline_free(r_spline);
gsl_interp_accel_free(r_acc);
gsl_spline_free(phi_spline);
gsl_interp_accel_free(phi_acc);
gsl_spline_free(pr_spline);
gsl_interp_accel_free(pr_acc);
gsl_spline_free(pphi_spline);
gsl_interp_accel_free(pphi_acc);

free(ts);
free(rs);
free(phis);
free(prs);
free(pphis);
"""
    cfc.register_CFunction(
        subdirectory="dynamics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
