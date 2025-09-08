"""
Set up C function library for SEOBNR and BOB attachment routines.

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


def register_CFunction_SEOBNRv5_aligned_spin_interpolate_modes() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for interpolating the (2,2) inspiral mode for the SEOBNRv5 waveform.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Interpolates the (2,2) inspiral modes and stores them in the waveform_inspiral array with constant spacing.

@param commondata - Common data structure containing the model parameters.
@param dT - Time step for interpolation.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_interpolate_modes"
    params = "commondata_struct *restrict commondata, const REAL dT"
    body = """
size_t i;
const size_t nsteps_inspiral_old = commondata->nsteps_inspiral;
const size_t nsteps_new = (size_t) ((commondata->waveform_inspiral[IDX_WF(commondata->nsteps_inspiral - 1,TIME)] - commondata->waveform_inspiral[IDX_WF(0,TIME)]) / dT) + 1; 
// build the complex inspiral modes
REAL *restrict times_old = (REAL *)malloc(nsteps_inspiral_old * sizeof(REAL));
REAL *restrict orbital_phases = (REAL *)malloc(nsteps_inspiral_old * sizeof(REAL));
REAL *restrict h22_nophase_real = (REAL *)malloc(nsteps_inspiral_old * sizeof(REAL));
REAL *restrict h22_nophase_imag = (REAL *)malloc(nsteps_inspiral_old * sizeof(REAL));
double complex h22_nophase , h22_rescaled;
for (i = 0; i < commondata->nsteps_low; i++){
  times_old[i] = commondata->waveform_inspiral[IDX_WF(i,TIME)];
  orbital_phases[i] = commondata->dynamics_low[IDX(i,PHI)];
  h22_nophase = cexp(2 * I * orbital_phases[i])*(commondata->waveform_inspiral[IDX_WF(i,STRAIN)]);
  h22_nophase_real[i] = creal(h22_nophase);
  h22_nophase_imag[i] = cimag(h22_nophase);
}
for (i = 0; i < commondata->nsteps_fine; i++){
  times_old[i + commondata->nsteps_low] = commondata->waveform_inspiral[IDX_WF(i + commondata->nsteps_low,TIME)];
  orbital_phases[i + commondata->nsteps_low] = commondata->dynamics_fine[IDX(i,PHI)];
  h22_nophase = cexp(2 * I * orbital_phases[i + commondata->nsteps_low])*(commondata->waveform_inspiral[IDX_WF(i + commondata->nsteps_low,STRAIN)]);
  h22_nophase_real[i + commondata->nsteps_low] = creal(h22_nophase);
  h22_nophase_imag[i + commondata->nsteps_low] = cimag(h22_nophase);
}

//interpolate and set the inspiral modes
REAL orbital_phase, h22_real , h22_imag, time;
const REAL tstart = commondata->waveform_inspiral[IDX_WF(0,TIME)];
gsl_interp_accel *restrict acc_real = gsl_interp_accel_alloc();
gsl_interp_accel *restrict acc_imag = gsl_interp_accel_alloc();
gsl_interp_accel *restrict acc = gsl_interp_accel_alloc();
gsl_spline *restrict spline_real = gsl_spline_alloc(gsl_interp_cspline, nsteps_inspiral_old);
gsl_spline *restrict spline_imag = gsl_spline_alloc(gsl_interp_cspline, nsteps_inspiral_old);
gsl_spline *restrict spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_inspiral_old);
gsl_spline_init(spline,times_old,orbital_phases, nsteps_inspiral_old);
gsl_spline_init(spline_real,times_old,h22_nophase_real, nsteps_inspiral_old);
gsl_spline_init(spline_imag,times_old,h22_nophase_imag, nsteps_inspiral_old);
//realloc the amount of memory needed to store the interpolated modes
commondata->nsteps_inspiral = nsteps_new;
commondata->waveform_inspiral = (double complex *)realloc(commondata->waveform_inspiral,commondata->nsteps_inspiral * NUMMODES * sizeof(double complex));
for (i = 0; i < commondata->nsteps_inspiral; i++){
  time = tstart + i * dT;
  commondata->waveform_inspiral[IDX_WF(i,TIME)] = time;
  orbital_phase = gsl_spline_eval(spline,time,acc);
  h22_real = gsl_spline_eval(spline_real,time,acc_real);
  h22_imag = gsl_spline_eval(spline_imag,time,acc_imag);
  h22_rescaled = cexp(-2. * I * orbital_phase) * (h22_real + I * h22_imag);
  commondata->waveform_inspiral[IDX_WF(i,STRAIN)] = h22_rescaled;
}
gsl_interp_accel_free(acc);
gsl_interp_accel_free(acc_real);
gsl_interp_accel_free(acc_imag);
gsl_spline_free(spline_real);
gsl_spline_free(spline_imag);
gsl_spline_free(spline);
free(times_old);
free(orbital_phases);
free(h22_nophase_real);
free(h22_nophase_imag);
"""
    cfc.register_CFunction(
        subdirectory="inspiral_waveform",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
