"""
Set up C function library for the SEOBNR aligned spin expressions.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import ast
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_waveform_quantities as SEOBNRv5_wf
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_constants as SEOBNRv5_const
import nrpy.helpers.parallel_codegen as pcg


def register_Cfunction_SEOBNRv5_aligned_spin_special_amplitude_coefficients_rholm() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register C function for computing and applying "special" amplitude coefficients needed for (2,1), (4,3), and (5,5) modes.

    :return: None if in registration phase, else the updated NRPy environment.
    """

    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = SEOBNRv5_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    const = SEOBNRv5_const.SEOBNR_aligned_spin_constants()

    rho_dict = wf.rho
    rholm = []
    rholm_labels = []

    hNR_fits = const.hNR
    hNR = []
    hNR_labels = []

    modes = [(2, 1), (4, 3), (5, 5)]

    for key in rho_dict.keys():
        mode = ast.literal_eval(key)
        l, m = mode
        for mode in modes:
            rholm.append(wf.rho[f"({l} , {m})"])
            rholm_labels.append(f"const REAL rho{l}{m}")

    for key in hNR_fits():
        mode = ast.literal_eval(key)
        l, m = mode
        for mode in modes:
            hNR.append(hNR_fits.hNR[f"({l} , {m})"])
            hNR_labels.append(f"const REAL hNR{l}{m}")

    rholm_code = ccg.c_codegen(
        rholm,
        rholm_labels,
        verbose=False,
        include_braces=False,
        cse_varprefix="rho",
    )

    hNR_code = ccg.c_codegen(
        hNR,
        hNR_labels,
        verbose=False,
        include_braces=False,
        cse_varprefix="hNR",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Computes and applies the special amplitude coefficients to inspiral waveform modes (2,1), (4,3), and (5,5).

@param commondata - Common data structure containing the model parameters.
"""
    cfunc_type = "void"
    prefunc = "#include<complex.h>"
    name = "SEOBNRv5_aligned_spin_special_coefficients_rholm_calculation"
    params = "commondata_struct *restrict commondata, REAL *rhos, REAL *restrict dynamics, REAL *hNR"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL c_21 = commondata->c_21;
const REAL c_43 = commondata->c_43;
const REAL c_55 = commondata->c_55;
const REAL Omega = dynamics[OMEGA];

"""

    body += rholm_code
    body += """
rhos[RHO21 - 1] = rho21;
rhos[RHO43 - 1] = rho43;
rhos[RHO55 - 1] = rho55;
"""

    body += hNR_code
    body += """
hNR[hNR21 - 1] = hNR21;
hNR[hNR43 - 1] = hNR43;
hNR[hNR55 - 1] = hNR55;
"""

    cfc.register_CFunction(
        subdirectory="inspiral_waveform",
        includes=includes,
        desc=desc,
        prefunc=prefunc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


def register_Cfunction_SEOBNRv5_aligned_spin_special_amplitude_coefficients() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register C function for computing and applying "special" amplitude coefficients needed for (2,1), (4,3), and (5,5) modes.

    :return: None if in registration phase, else the updated NRPy environment.
    """

    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Computes and applies the special amplitude coefficients to inspiral waveform modes (2,1), (4,3), and (5,5).

@param commondata - Common data structure containing the model parameters.
"""
    cfunc_type = "void"
    prefunc = "#include<complex.h>"
    name = "SEOBNRv5_aligned_spin_special_coefficients"
    params = "commondata_struct *restrict commondata, REAL *rhos, double complex *inspiral_modes, REAL *restrict dynamics"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL c_21 = commondata->c_21;
const REAL c_43 = commondata->c_43;
const REAL c_55 = commondata->c_55;
const REAL Omega = dynamics[OMEGA];
const REAL rho21 = rhos[RHO21];
const REAL rho43 = rhos[RHO43];
const REAL rho55 = rhos[RHO55];
const REAL hNR21 = hNR[hNR21];
const REAL hNR43 = hNR[hNR43];
const REAL hNR55 = hNR[hNR55];
double complex h21 = inspiral_modes[STRAIN21];
double complex h43 = inspiral_modes[STRAIN43];
double complex h55 = inspiral_modes[STRAIN55];

size_t i;

for (i = 0; i < commondata->nsteps_fine; i++){
  Omega[i] = commondata->dynamics_fine[IDX(i,OMEGA)];
  h21[i] = commondata->waveform_fine[IDX_WF(i,STRAIN21)];
  h43[i] = commondata->waveform_fine[IDX_WF(i,STRAIN43)];
  times[i] = commondata->dynamics_fine[IDX(i,TIME)];
  K21[i] = h21[i] / rho21;
  K43[i] = h43[i] / rho43;
  v[i] = cbrt(Omega[i])
  
}

SEOBNRv5_aligned_spin_unwrap(phase,phase_unwrapped,commondata->nsteps_fine);
// Find t_ISCO:

if (commondata->r_ISCO < r[commondata->nsteps_fine - 1]){
  commondata->t_ISCO = times[commondata->nsteps_fine - 1];
}
else{
  const REAL dt_ISCO = 0.001;
  const size_t N_zoom = (size_t) ((times[commondata->nsteps_fine - 1] - times[0]) / dt_ISCO);
  REAL *restrict t_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  REAL *restrict minus_r_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  gsl_interp_accel *restrict acc_r = gsl_interp_accel_alloc();
  gsl_spline *restrict spline_r = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  gsl_spline_init(spline_r,times,r,commondata->nsteps_fine);
  for (i = 0; i < N_zoom; i++){
    t_zoom[i] = times[0] + i * dt_ISCO;
    minus_r_zoom[i] = -1.0*gsl_spline_eval(spline_r,t_zoom[i],acc_r);
  }
  const size_t ISCO_zoom_idx = gsl_interp_bsearch(minus_r_zoom, -commondata->r_ISCO, 0 , N_zoom);
  commondata->t_ISCO = t_zoom[ISCO_zoom_idx];

  gsl_interp_accel_free(acc_r);
  gsl_spline_free(spline_r);
  free(t_zoom);
  free(minus_r_zoom);
}


for (i = 0; i < commondata->nsteps_fine; i++){
  Omega[i] = commondata->dynamics_fine[IDX(i,OMEGA)];
  h55[i] = commondata->waveform_fine[IDX_WF(i,STRAIN55)];
  times[i] = commondata->dynamics_fine[IDX(i,TIME)];
  K55[i] = h55[i] / rho55;
  v[i] = cbrt(Omega[i])
  
}

SEOBNRv5_aligned_spin_unwrap(phase,phase_unwrapped,commondata->nsteps_fine);
// Find t_ISCO:

if (commondata->r_ISCO < r[commondata->nsteps_fine - 1]){
  commondata->t_ISCO = times[commondata->nsteps_fine - 1];
}
else{
  const REAL dt_ISCO = 0.001;
  const size_t N_zoom = (size_t) ((times[commondata->nsteps_fine - 1] - times[0]) / dt_ISCO);
  REAL *restrict t_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  REAL *restrict minus_r_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  gsl_interp_accel *restrict acc_r = gsl_interp_accel_alloc();
  gsl_spline *restrict spline_r = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  gsl_spline_init(spline_r,times,r,commondata->nsteps_fine);
  for (i = 0; i < N_zoom; i++){
    t_zoom[i] = times[0] + i * dt_ISCO;
    minus_r_zoom[i] = -1.0*gsl_spline_eval(spline_r,t_zoom[i],acc_r);
  }
  const size_t ISCO_zoom_idx = gsl_interp_bsearch(minus_r_zoom, -commondata->r_ISCO, 0 , N_zoom);
  commondata->t_ISCO = t_zoom[ISCO_zoom_idx];

  gsl_interp_accel_free(acc_r);
  gsl_spline_free(spline_r);
  free(t_zoom);
  free(minus_r_zoom);
}
"""
