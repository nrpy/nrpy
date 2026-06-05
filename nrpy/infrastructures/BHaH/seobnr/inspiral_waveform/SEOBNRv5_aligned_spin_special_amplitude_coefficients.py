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

    wf = SEOBNRv5_wf.SEOBNRv5_aligned_spin_waveform_quantities(
        apply_special_amplitude_coefficients=True
    )

    rholm: List[sp.Expr] = []
    rholm_labels: List[str] = []

    modes = [(2, 1), (4, 3), (5, 5)]

    for l, m in modes:
        rholm.append(cast(sp.Expr, wf.pn_contribution_f[f"({l} , {m})"]))
        rholm_labels.append(f"REAL rho{l}{m}")

    rholm_code = ccg.c_codegen(
        rholm,
        rholm_labels,
        verbose=False,
        include_braces=False,
        cse_varprefix="rho",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Computes and applies the special amplitude coefficients to inspiral waveform modes (2,1), (4,3), and (5,5).

@param commondata - Common data structure containing the model parameters.
"""
    cfunc_type = "void"
    prefunc = "#include<complex.h>"
    name = "SEOBNRv5_aligned_spin_special_coefficients_rholm"
    params = (
        "commondata_struct *restrict commondata,  REAL *restrict dynamics, REAL *rhos"
    )
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

    v5_const = SEOBNRv5_const.SEOBNR_aligned_spin_constants()
    projected_risco_code = ccg.c_codegen(
        [v5_const.rISCO],
        ["const REAL projected_r_ISCO"],
        verbose=False,
        include_braces=False,
        cse_varprefix="projected_risco",
    )
    projected_delta_t_code = ccg.c_codegen(
        [v5_const.Delta_t],
        ["const REAL projected_Delta_t"],
        verbose=False,
        include_braces=False,
        cse_varprefix="projected_delta_t",
    )

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


const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL nu = m1 * m2/((m1 + m2) * (m1 + m2));
const REAL chiA = (chi1 - chi2) / 2;
const int use_projected_attachment =
    fabs(commondata->chi1 - commondata->chi1_z) <= 1e-14 &&
    fabs(commondata->chi2 - commondata->chi2_z) <= 1e-14;
REAL rhos[NUMVARS_COEFFICIENTS];
REAL hNR[NUMVARS_HNRFITS];
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

// Step 1: Build combined time-frequency samples for projected-spin attachment fits.
if (use_projected_attachment) {
  const size_t nsteps_combined = commondata->nsteps_low + commondata->nsteps_fine;
  if (commondata->nsteps_low < 2 || commondata->nsteps_fine < 2) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), insufficient dynamics samples for projected attachment inputs\\n");
    exit(1);
  }
if (commondata->chi1_lnhat.spline == NULL || commondata->chi1_lnhat.acc == NULL ||
    commondata->chi2_lnhat.spline == NULL || commondata->chi2_lnhat.acc == NULL) {
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), projected-spin splines are unavailable\\n");
  exit(1);
}

REAL *restrict times_combined = (REAL *)malloc(nsteps_combined*sizeof(REAL));
if (times_combined == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed for times_combined\\n");
  exit(1);
}
REAL *restrict Omega_combined = (REAL *)malloc(nsteps_combined*sizeof(REAL));
if (Omega_combined == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed for Omega_combined\\n");
  exit(1);
}
REAL *restrict u_rlow = (REAL *)malloc(commondata->nsteps_low*sizeof(REAL));
if (u_rlow == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed for u_rlow\\n");
  exit(1);
}
REAL *restrict t_rlow = (REAL *)malloc(commondata->nsteps_low*sizeof(REAL));
if (t_rlow == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed for t_rlow\\n");
  exit(1);
}
for (i = 0; i < commondata->nsteps_low; i++){
  const REAL r_low_i = commondata->dynamics_low[IDX(i,R)];
  if (r_low_i <= 0.0) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), nonpositive low-dynamics radius in projected attachment inputs\\n");
    exit(1);
  }
  u_rlow[i] = 1.0 / r_low_i;
  t_rlow[i] = commondata->dynamics_low[IDX(i,TIME)];
  times_combined[i] = t_rlow[i];
  Omega_combined[i] = commondata->dynamics_low[IDX(i,OMEGA)];
} // END LOOP: for i over low-dynamics projected attachment samples
for (i = 0; i < commondata->nsteps_fine; i++){
  const size_t dst_idx = commondata->nsteps_low + i;
  times_combined[dst_idx] = times[i];
  Omega_combined[dst_idx] = Omega[i];
} // END LOOP: for i over fine-dynamics projected attachment samples
REAL u_r10M = 0.1;
if (u_r10M < u_rlow[0]) {
  fprintf(stderr,"Warning: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), r=10M lies before the low-dynamics projected-spin reference domain; using first available low-dynamics point\\n");
  u_r10M = u_rlow[0];
}
if (u_r10M > u_rlow[commondata->nsteps_low - 1]) {
  fprintf(stderr,"Warning: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), r=10M lies after the low-dynamics projected-spin reference domain; using last available low-dynamics point\\n");
  u_r10M = u_rlow[commondata->nsteps_low - 1];
}

gsl_interp_accel *restrict acc_t_of_u = gsl_interp_accel_alloc();
if (acc_t_of_u == NULL) {
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_interp_accel_alloc() failed for acc_t_of_u\\n");
  exit(1);
}
gsl_spline *restrict spline_t_of_u = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_low);
if (spline_t_of_u == NULL) {
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_spline_alloc() failed for spline_t_of_u\\n");
  exit(1);
}
gsl_spline_init(spline_t_of_u,u_rlow,t_rlow,commondata->nsteps_low);

gsl_interp_accel *restrict acc_Omega_combined = gsl_interp_accel_alloc();
if (acc_Omega_combined == NULL) {
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_interp_accel_alloc() failed for acc_Omega_combined\\n");
  exit(1);
}
gsl_spline *restrict spline_Omega_combined = gsl_spline_alloc(gsl_interp_cspline, nsteps_combined);
if (spline_Omega_combined == NULL) {
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_spline_alloc() failed for spline_Omega_combined\\n");
  exit(1);
}
gsl_spline_init(spline_Omega_combined,times_combined,Omega_combined,nsteps_combined);

const REAL t_r10M = gsl_spline_eval(spline_t_of_u, u_r10M, acc_t_of_u);
REAL omega_r10M = gsl_spline_eval(spline_Omega_combined, t_r10M, acc_Omega_combined);
omega_r10M = fmin(commondata->omega_spin_max, fmax(commondata->omega_spin_min, omega_r10M));
const REAL chi1_projected_r10 = gsl_spline_eval(commondata->chi1_lnhat.spline, omega_r10M, commondata->chi1_lnhat.acc);
const REAL chi2_projected_r10 = gsl_spline_eval(commondata->chi2_lnhat.spline, omega_r10M, commondata->chi2_lnhat.acc);
{
  const REAL chi1 = chi1_projected_r10;
  const REAL chi2 = chi2_projected_r10;
"""
    body += projected_risco_code
    body += """
  commondata->r_ISCO = projected_r_ISCO;
} // END BLOCK: projected-spin r_ISCO evaluation at r=10M
gsl_spline_free(spline_t_of_u);
gsl_interp_accel_free(acc_t_of_u);
gsl_spline_free(spline_Omega_combined);
gsl_interp_accel_free(acc_Omega_combined);
free(times_combined);
free(Omega_combined);
free(u_rlow);
free(t_rlow);
} // END IF: projected-spin attachment inputs are self-consistent

// construct splines of dynamical variables

gsl_interp_accel *restrict acc_r = gsl_interp_accel_alloc();
  if (acc_r == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_r = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_r == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_r,times,r,commondata->nsteps_fine);
 

gsl_interp_accel *restrict acc_phi = gsl_interp_accel_alloc();
  if (acc_phi == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_phi = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_phi == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_phi,times,phi,commondata->nsteps_fine);
  

gsl_interp_accel *restrict acc_prstar = gsl_interp_accel_alloc();
  if (acc_prstar == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_prstar = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_prstar == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_prstar,times,prstar,commondata->nsteps_fine);
  


gsl_interp_accel *restrict acc_pphi = gsl_interp_accel_alloc();
  if (acc_pphi == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_pphi = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_pphi == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_pphi,times,pphi,commondata->nsteps_fine);
 

gsl_interp_accel *restrict acc_Omega = gsl_interp_accel_alloc();
  if (acc_Omega == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_Omega = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_Omega == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_Omega,times,Omega,commondata->nsteps_fine);
  

gsl_interp_accel *restrict acc_Omega_circ = gsl_interp_accel_alloc();
  if (acc_Omega_circ == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_Omega_circ = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_Omega_circ == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_Omega_circ,times,Omega_circ,commondata->nsteps_fine); 
  
gsl_interp_accel *restrict acc_Hreal = gsl_interp_accel_alloc();
  if (acc_Hreal == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_interp_accel_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline *restrict spline_Hreal = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_Hreal == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_spline_alloc() failed to initialize\\n");
      exit(1);
    }
  gsl_spline_init(spline_Hreal,times,Hreal,commondata->nsteps_fine);

// Step 2: Find t_ISCO with the projected-spin r_ISCO.

if (commondata->r_ISCO < r[commondata->nsteps_fine - 1]){
  commondata->t_ISCO = times[commondata->nsteps_fine - 1];
}
else{
  const REAL dt_ISCO = 0.001;
  const size_t N_zoom = (size_t) ((times[commondata->nsteps_fine - 1] - times[0]) / dt_ISCO);
  if (N_zoom == 0) {
    fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), fine dynamics time interval is too short for t_ISCO search\\n");
    exit(1);
  }
  REAL *restrict t_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  if (t_zoom == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed for t_zoom\\n");
      exit(1);
    }
  REAL *restrict minus_r_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  if (minus_r_zoom == NULL) {
      fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), malloc() failed for times\\n");
      exit(1);
    }
  for (i = 0; i < N_zoom; i++){
    t_zoom[i] = times[0] + i * dt_ISCO;
    minus_r_zoom[i] = -1.0*gsl_spline_eval(spline_r,t_zoom[i],acc_r);
  }
  size_t ISCO_zoom_idx = 0;
  if (use_projected_attachment) {
    REAL min_abs_r_gap = fabs(minus_r_zoom[0] + commondata->r_ISCO);
    for (i = 1; i < N_zoom; i++){
      const REAL abs_r_gap = fabs(minus_r_zoom[i] + commondata->r_ISCO);
      if (abs_r_gap < min_abs_r_gap) {
        min_abs_r_gap = abs_r_gap;
        ISCO_zoom_idx = (size_t)i;
      }
    } // END LOOP: for i over fine-grid projected r_ISCO search samples
  } else {
    ISCO_zoom_idx = gsl_interp_bsearch(minus_r_zoom, -commondata->r_ISCO, 0 , N_zoom - 1);
  } // END ELSE: legacy scalar aligned-spin r_ISCO search
  commondata->t_ISCO = t_zoom[ISCO_zoom_idx];
  
  free(t_zoom);
  free(minus_r_zoom);
}

// Step 3: Evaluate the attachment-time shift from projected spins at r_ISCO.
if (use_projected_attachment) {
  REAL omega_rISCO = gsl_spline_eval(spline_Omega, commondata->t_ISCO, acc_Omega);
  omega_rISCO = fmin(commondata->omega_spin_max, fmax(commondata->omega_spin_min, omega_rISCO));
  const REAL chi1_projected_rISCO = gsl_spline_eval(commondata->chi1_lnhat.spline, omega_rISCO, commondata->chi1_lnhat.acc);
  const REAL chi2_projected_rISCO = gsl_spline_eval(commondata->chi2_lnhat.spline, omega_rISCO, commondata->chi2_lnhat.acc);
  {
    const REAL chi1 = chi1_projected_rISCO;
    const REAL chi2 = chi2_projected_rISCO;
"""
    body += projected_delta_t_code
    body += """
    commondata->Delta_t = projected_Delta_t;
  } // END BLOCK: projected-spin Delta_t evaluation at r_ISCO
} // END IF: projected-spin Delta_t inputs are self-consistent

REAL t_peak_22 = commondata->t_ISCO - commondata->Delta_t;
REAL t_peak_55 = t_peak_22 - 10;
const size_t attachment_end_idx = use_projected_attachment ? commondata->nsteps_fine - 1 : commondata->nsteps_fine - 2;

if (t_peak_22 > times[commondata->nsteps_fine - 1]){
  t_peak_22 = times[attachment_end_idx];
  t_peak_55 = t_peak_22;
}
if (t_peak_55 > times[commondata->nsteps_fine - 1]){
  t_peak_55 = times[attachment_end_idx];
}
commondata->t_attach = t_peak_22;

if (use_projected_attachment) {
  REAL omega_attach = gsl_spline_eval(spline_Omega, t_peak_22, acc_Omega);
  omega_attach = fmin(commondata->omega_spin_max, fmax(commondata->omega_spin_min, omega_attach));
  commondata->chi1 = gsl_spline_eval(commondata->chi1_lnhat.spline, omega_attach, commondata->chi1_lnhat.acc);
  commondata->chi2 = gsl_spline_eval(commondata->chi2_lnhat.spline, omega_attach, commondata->chi2_lnhat.acc);
} // END IF: projected-spin waveform inputs use attachment-time spin projections

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

const REAL hNR21_threshold = 300;
const REAL hNR43_threshold = 200 * nu * (1 - 0.8 * chiA);
const REAL hNR55_threshold = 2000;

SEOBNRv5_aligned_spin_special_coefficients_rholm(commondata, dynamics_22, rhos);
REAL rho21 = rhos[RHO21];
REAL rho43 = rhos[RHO43];

SEOBNRv5_aligned_spin_hNR_fits_at_t_attach(commondata, hNR);
REAL hNR21 = hNR[HNR21] * nu;
REAL hNR43 = hNR[HNR43] * nu;
REAL hNR55 = hNR[HNR55] * nu;
REAL hNR22 = hNR[HNR22] * nu;

if (fabs(hNR21) < hNR22 / hNR21_threshold) {
    hNR21 = copysign(hNR22/hNR21_threshold, hNR21);
}

if (fabs(hNR43) < hNR22 / hNR43_threshold) {
    hNR43 = copysign(hNR22/hNR43_threshold, hNR43);
}

SEOBNRv5_aligned_spin_special_coefficients_rholm(commondata, dynamics_55, rhos);
REAL rho55 = rhos[RHO55];

if (fabs(hNR55) < hNR22 / hNR55_threshold) {
    hNR55 = copysign(hNR22/hNR55_threshold, hNR55);
}

SEOBNRv5_aligned_spin_waveform(dynamics_22, commondata, inspiral_modes);
double complex h21 = inspiral_modes[STRAIN21 - 1];
double complex h43 = inspiral_modes[STRAIN43 - 1];

SEOBNRv5_aligned_spin_waveform(dynamics_55, commondata, inspiral_modes);
double complex h55 = inspiral_modes[STRAIN55 - 1];

const REAL v22 = pow(dynamics_22[OMEGA], 1.0/3.0);
const REAL K21 = cabs(h21) / fabs(rho21);
const REAL K43 = cabs(h43) / fabs(rho43);
const REAL v55 = pow(dynamics_55[OMEGA], 1.0/3.0);
const REAL K55 = cabs(h55) / fabs(rho55);

const REAL vpow21 = pow(v22, 7.0);
const REAL vpow43 = pow(v22, 7.0);
const REAL vpow55 = pow(v55, 5.0);

const REAL c21 = (hNR21/K21 - rho21)/vpow21;
const REAL c43 = (hNR43/K43 - rho43)/vpow43;
const REAL c55 = (hNR55/K55 - rho55)/vpow55;

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
