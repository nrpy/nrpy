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
    projected_remnant_code = ccg.c_codegen(
        [v5_const.M_f, v5_const.a_f, v5_const.rISCO],
        ["commondata->M_f", "commondata->a_f", "commondata->r_ISCO"],
        verbose=False,
        include_braces=False,
        cse_varprefix="projected_remnant",
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
commondata->projected_attachment_active = false;
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
  commondata->projected_attachment_active = true;
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
    body += projected_remnant_code
    body += """
  } // END BLOCK: projected-spin remnant evaluation at r=10M
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
    0.74632,0.7676741,0.7932082,0.8082349,0.8254294,0.8331,0.8413426,0.8502722,0.8600456,0.8708927,0.8830905,0.8969183,0.9045305,
    0.912655,0.9213264,0.9258781,0.9305797,0.9354355,0.9364255,0.937422,0.9384248,0.9394341,0.9404498,0.9425009,0.9445784,
    0.9456271,0.9466825 };
  const REAL imomegaqnm22[107] = { 0.0880269,0.0880272,0.0880274,0.088028,0.0880285,0.0880288,0.088029,0.0880293,
    0.0880296,0.0880298,0.0880311,0.0880325,0.0880338,0.0880364,0.0880391,0.0880417,0.088047,0.0880523,0.0880575,0.0880628,0.088068,
    0.0880733,0.0880785,0.0880915,0.0881045,0.0881304,0.088156,0.0881813,0.0882315,0.0882807,0.0883289,0.0883763,0.0884226,0.0884679,
    0.0885122,0.0885555,0.0885976,0.0886386,0.0886785,0.0887172,0.0888085,0.0888917,0.0889663,0.0890315,0.0890868,0.0891313,0.0891643,
    0.0891846,0.0891911,0.0891825,0.0891574,0.0891138,0.0890496,0.0889623,0.0888489,0.0887057,0.0885283,0.0883112,0.0880477,0.0877293,
    0.0873453,0.086882,0.0863212,0.0856388,0.0848021,0.0837652,0.0824618,0.0807929,0.0799908,0.0790927,0.0780817,0.0769364,0.0756296,
    0.0741258,0.072378,0.0703215,0.0678642,0.0648692,0.0611186,0.0562313,0.053149,0.0494336,0.0447904,0.0419586,0.0386302,0.0371155,
    0.0354677,0.033659,0.0316517,0.0293904,0.0268082,0.0238377,0.0221857,0.0204114,0.0185063,0.0175021,0.016462,0.015385,0.0151651,
    0.0149437,0.0147207,0.0144962,0.0142701,0.0138132,0.0133501,0.0131161,0.0128806 };
  gsl_interp_accel *restrict acc_qnm = gsl_interp_accel_alloc();
  if (acc_qnm == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_interp_accel_alloc() failed for acc_qnm\\n");
    exit(1);
  }
  gsl_spline *restrict spline_qnm = gsl_spline_alloc(gsl_interp_cspline, 107);
  if (spline_qnm == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_special_amplitude_coefficients(), gsl_spline_alloc() failed for spline_qnm\\n");
    exit(1);
  }
  gsl_spline_init(spline_qnm, afinallist, reomegaqnm22, 107);
  commondata->omega_qnm = gsl_spline_eval(spline_qnm, commondata->a_f, acc_qnm) / commondata->M_f;
  gsl_spline_init(spline_qnm, afinallist, imomegaqnm22, 107);
  gsl_interp_accel_reset(acc_qnm);
  commondata->tau_qnm = 1.0 / (gsl_spline_eval(spline_qnm, commondata->a_f, acc_qnm) / commondata->M_f);
  gsl_spline_free(spline_qnm);
  gsl_interp_accel_free(acc_qnm);
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
