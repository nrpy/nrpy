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

    for key in hNR_fits.keys():
        mode = ast.literal_eval(key)
        l, m = mode
        for mode in modes:
            hNR.append(hNR_fits[f"({l} , {m})"])
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
const REAL Omega = dynamics[OMEGA];
REAL c_21 = c_21;
REAL c_43 = c_43;
REAL c_55 = c_55;
REAL rhos[3];
REAL hNR[3];
REAL inspiral_modes[NUMMODES];

size_t i;

for (i = 0; i < commondata->nsteps_fine; i++){
  Omega[i] = commondata->dynamics_fine[IDX(i,OMEGA)];
  Omega_circ[i] = commondata->dynamics_fine[IDX(i,OMEGA_CIRC)];
  Hreal[i] = commondata->dynamics_fine[IDX(i,HREAL)];
  R[i] = commondata->dynamics_fine[IDX(i,R)];
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
  gsl_spline *restrict spline_Omegcirc = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_Omega_circ == NULL) {
      fprintf(stderr, "Error: in SEBOBv2_NQC_corrections(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_Omega_circ,times,Omega_circ,commondata->nsteps_fine); 

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
  
  const size_t ISCO_zoom_idx = gsl_interp_bsearch(minus_r_zoom, -commondata->r_ISCO, 0 , N_zoom);
  commondata->t_ISCO = t_zoom[ISCO_zoom_idx];
}

REAL t_peak_22 = commondata->t_ISCO - commondata->Delta_t;
REAL t_peak_55 = t_peak_22 + 10;
size_t_22 peak_idx; 

if (t_peak_22 > times[commondata->nsteps_fine - 1]){
  t_peak_22 = times[commondata->nsteps_fine - 2];
  t_peak_55 = t_peak_22;
  peak_idx = commondata->nsteps_fine - 2;
}
else{
  peak_idx = gsl_interp_bsearch(times, t_peak_22, 0, commondata->nsteps_fine);
}  
commondata->t_attach = t_peak_22;

REAL dynamics_22[NUMVARS];
REAL dynamics_55[NUMVARS];

//for(i = 0; i < size_t_22; i++){
//    dynamics_22[R] = 


SEOBNRv5_aligned_spin_special_amplitude_coefficients_rholm(commondata, dynamics, rhos, hNR);
rho21 = rhos[RHO21 - 1];
rho43 = rhos[RHO43 - 1];
rho55 = rhos[RHO55 - 1];
hNR21 = hNR[hNR21 - 1];
hNR43 = hNR[hNR43 - 1];
hNR55 = hNR[hNR55 - 1];

SEOBNRv5_aligned_spin_waveform(commondata, dynamics, inspiral_modes);
h21 = inspiral_modes[STRAIN21 - 1];
h43 = inspiral_modes[STRAIN43 - 1];
h55 = inspiral_modes[STRAIN55 - 1];

v22 = cpow(dynamics_22[Omega], 1, 3);
K21 = h21 / rho21;
K43 = h43 / rho43;
v55 = cpow(dynamics_55[Omega], 1, 3);
K55 = h55/rho55;

c21 = (hNR21/K21)/v22;
c43 = (hNR43/K43)/v22;
c55 = (hNR55/K55)/v55;

//save to commondata

commondata->c_21 = c21;
commondata->c_43 = c43;
commondata->c_55 = c55;

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
