"""
Set up C function library for the SEOBNR aligned spin expressions.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_constants as SEOBNRv5_const
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_waveform_quantities as SEOBNRv5_wf
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

    rholm: List[sp.Expr] = []
    rholm_labels: List[str] = []

    hNR: List[sp.Expr] = []
    hNR_labels: List[str] = []

    modes = [(2, 1), (4, 3), (5, 5)]

    for l, m in modes:
        rholm.append(wf.rho[f"({l} , {m})"])
        rholm_labels.append(f"REAL rho{l}{m}")

    for i in range(len(rholm)):
        cast(sp.Expr, rholm[i])

    for l, m in modes:
        hNR.append(const.hNR[f"({l} , {m})"])
        hNR_labels.append(f"const REAL hNR{l}{m}")

    for i in range(len(hNR)):
        cast(sp.Expr, hNR[i])

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
    name = "SEOBNRv5_aligned_spin_special_coefficients_rholm"
    params = "commondata_struct *restrict commondata,  REAL *restrict dynamics, REAL *rhos, REAL *hNR"
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
rhos[RHO21] = rho21;
rhos[RHO43] = rho43;
rhos[RHO55] = rho55;
"""

    body += hNR_code
    body += """
hNR[HNR21] = hNR21;
hNR[HNR43] = hNR43;
hNR[HNR55] = hNR55;
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
    params = "commondata_struct *restrict commondata"
    body = """
    
REAL *restrict times = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (times == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed to for times\\n");
  exit(1);
}

REAL *restrict r = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (r == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed to for times\\n");
  exit(1);
}

REAL *restrict phi = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (phi == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed to for times\\n");
  exit(1);
}

REAL *restrict pphi = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (pphi == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed to for times\\n");
  exit(1);
}

REAL *restrict prstar = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (prstar == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed to for times\\n");
  exit(1);
}

REAL *restrict Omega = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Omega == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed to for times\\n");
  exit(1);
}

REAL *restrict Hreal = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Hreal == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed to for times\\n");
  exit(1);
}

REAL *restrict Omega_circ = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Omega_circ == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed to for times\\n");
  exit(1);
}


REAL c_21 = c_21;
REAL c_43 = c_43;
REAL c_55 = c_55;
REAL rhos[NUMVARS_COEFFICIENTS];
REAL hNR[NUMVARS_COEFFICIENTS];
double complex inspiral_modes[NUMMODES];

size_t i;

for (i = 0; i < commondata->nsteps_fine; i++){
  Omega[i] = commondata->dynamics_fine[IDX(i,OMEGA)];
  Omega_circ[i] = commondata->dynamics_fine[IDX(i,OMEGA_CIRC)];
  Hreal[i] = commondata->dynamics_fine[IDX(i,H)];
  r[i] = commondata->dynamics_fine[IDX(i,R)];
  phi[i] = commondata->dynamics_fine[IDX(i,PHI)];
  prstar[i] = commondata->dynamics_fine[IDX(i,PRSTAR)];
  pphi[i] = commondata->dynamics_fine[IDX(i,PPHI)];
  times[i] = commondata->dynamics_fine[IDX(i,TIME)];
}


// construct splines of dynamical variables

gsl_interp_accel *restrict acc_r = gsl_interp_accel_alloc();
  if (acc_r == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_r = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_r == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_r,times,r,commondata->nsteps_fine);
 

gsl_interp_accel *restrict acc_phi = gsl_interp_accel_alloc();
  if (acc_phi == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_phi = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_phi == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_phi,times,phi,commondata->nsteps_fine);
  

gsl_interp_accel *restrict acc_prstar = gsl_interp_accel_alloc();
  if (acc_prstar == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_prstar = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_prstar == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_prstar,times,prstar,commondata->nsteps_fine);
  


gsl_interp_accel *restrict acc_pphi = gsl_interp_accel_alloc();
  if (acc_pphi == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_pphi = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_pphi == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_pphi,times,pphi,commondata->nsteps_fine);
 

gsl_interp_accel *restrict acc_Omega = gsl_interp_accel_alloc();
  if (acc_Omega == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_Omega = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_Omega == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_Omega,times,Omega,commondata->nsteps_fine);
  

gsl_interp_accel *restrict acc_Omega_circ = gsl_interp_accel_alloc();
  if (acc_Omega_circ == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_Omega_circ = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_Omega_circ == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_Omega_circ,times,Omega_circ,commondata->nsteps_fine); 
  
gsl_interp_accel *restrict acc_Hreal = gsl_interp_accel_alloc();
  if (acc_Hreal == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_Hreal = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_Hreal == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_Hreal,times,Hreal,commondata->nsteps_fine);

// Find t_ISCO:

if (commondata->r_ISCO < r[commondata->nsteps_fine - 1]){
  commondata->t_ISCO = times[commondata->nsteps_fine - 1];
}
else{
  const REAL dt_ISCO = 0.001;
  const size_t N_zoom = (size_t) ((times[commondata->nsteps_fine - 1] - times[0]) / dt_ISCO);
  REAL *restrict t_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  if (t_zoom == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed for t_zoom\\n");
      exit(1);
    }
  REAL *restrict minus_r_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  if (minus_r_zoom == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), malloc() failed for times\\n");
      exit(1);
    }
    for (i = 0; i < N_zoom; i++){
    t_zoom[i] = times[0] + i * dt_ISCO;
    minus_r_zoom[i] = -1.0*gsl_spline_eval(spline_r,t_zoom[i],acc_r);
  }
  const size_t ISCO_zoom_idx = gsl_interp_bsearch(minus_r_zoom, -commondata->r_ISCO, 0 , N_zoom);
  commondata->t_ISCO = t_zoom[ISCO_zoom_idx];
  
  free(t_zoom);
  free(minus_r_zoom);
}

REAL t_peak_22 = commondata->t_ISCO - commondata->Delta_t;
REAL t_peak_55 = t_peak_22 + 10;

if (t_peak_22 > times[commondata->nsteps_fine - 1]){
  t_peak_22 = times[commondata->nsteps_fine - 2];
  t_peak_55 = t_peak_22;
}
if (t_peak_55 > times[commondata->nsteps_fine - 1]){
  t_peak_55 = times[commondata->nsteps_fine - 2];  
}
commondata->t_attach = t_peak_22;

REAL dynamics_22[NUMVARS];
REAL dynamics_55[NUMVARS];

dynamics_22[TIME] = t_peak_22;
dynamics_22[R] = gsl_spline_eval(spline_r,t_peak_22,acc_r);
dynamics_22[PHI] = gsl_spline_eval(spline_phi,t_peak_22,acc_phi);
dynamics_22[PPHI] = gsl_spline_eval(spline_pphi,t_peak_22,acc_pphi);
dynamics_22[PRSTAR] = gsl_spline_eval(spline_prstar,t_peak_22,acc_prstar);
dynamics_22[OMEGA] = gsl_spline_eval(spline_Omega,t_peak_22,acc_Omega);
dynamics_22[H] = gsl_spline_eval(spline_Hreal,t_peak_22,acc_Hreal);
dynamics_22[OMEGA_CIRC] = gsl_spline_eval(spline_Omega_circ,t_peak_22,acc_Omega_circ);


dynamics_55[TIME] = t_peak_55;
dynamics_55[R] = gsl_spline_eval(spline_r,t_peak_55,acc_r);
dynamics_55[PHI] = gsl_spline_eval(spline_phi,t_peak_55,acc_phi);
dynamics_55[PPHI] = gsl_spline_eval(spline_pphi,t_peak_55,acc_pphi);
dynamics_55[PRSTAR] = gsl_spline_eval(spline_prstar,t_peak_55,acc_prstar);
dynamics_55[OMEGA] = gsl_spline_eval(spline_Omega,t_peak_55,acc_Omega);
dynamics_55[H] = gsl_spline_eval(spline_Hreal,t_peak_55,acc_Hreal);
dynamics_55[OMEGA_CIRC] = gsl_spline_eval(spline_Omega_circ,t_peak_55,acc_Omega_circ);


SEOBNRv5_aligned_spin_special_coefficients_rholm(commondata, dynamics_22, rhos, hNR);
REAL rho21 = rhos[RHO21];
REAL rho43 = rhos[RHO43];
const REAL hNR21 = hNR[HNR21];
const REAL hNR43 = hNR[HNR43];

SEOBNRv5_aligned_spin_special_coefficients_rholm(commondata, dynamics_55, rhos, hNR);
REAL rho55 = rhos[RHO55];
const REAL hNR55 = hNR[HNR55];

SEOBNRv5_aligned_spin_waveform(dynamics_22, commondata, inspiral_modes);
double complex h21 = inspiral_modes[STRAIN21 - 1];
double complex h43 = inspiral_modes[STRAIN43 - 1];

SEOBNRv5_aligned_spin_waveform(dynamics_55, commondata, inspiral_modes);
double complex h55 = inspiral_modes[STRAIN55 - 1];

const REAL v22 = pow(dynamics_22[OMEGA], 1.0/3.0);
const double complex K21 = h21 / rho21;
const double complex K43 = h43 / rho43;
const REAL v55 = pow(dynamics_55[OMEGA], 1.0/3.0);
const double complex K55 = h55/rho55;

const REAL c21 = cabs((hNR21/K21)/v22);
const REAL c43 = cabs((hNR43/K43)/v22);
const REAL c55 = cabs((hNR55/K55)/v55);

//save to commondata

commondata->c_21 = c21;
commondata->c_43 = c43;
commondata->c_55 = c55;


free(times);
free(r);
free(phi);
free(pphi);
free(prstar);
free(Omega);
free(Hreal);
free(Omega_circ);


gsl_spline_free(spline_r);
gsl_interp_accel_free(acc_r);
gsl_spline_free(spline_phi);
gsl_interp_accel_free(acc_phi);
gsl_spline_free(spline_prstar);
gsl_interp_accel_free(acc_prstar);
gsl_spline_free(spline_pphi);
gsl_interp_accel_free(acc_pphi);
gsl_spline_free(spline_Omega);
gsl_interp_accel_free(acc_Omega);
gsl_spline_free(spline_Hreal);
gsl_interp_accel_free(acc_Hreal);
gsl_spline_free(spline_Omega_circ);
gsl_interp_accel_free(acc_Omega_circ);

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
