"""
Set up C function to force BOBv2 strain peak to align with attachment time.

Authors: Siddharth Mahesh
sm0193 at mix dot wvu dot edu
Authors: Siddharth Mahesh
sm0193 at mix dot wvu dot edu
Zachariah B. Etienne
zachetie at gmail *dot com
Anuj Kankani
aak00009 at mix dot wvu dot edu

"""
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.equations.seobnr.BOB_v2_waveform_quantities_kankani_etal import BOB_v2_waveform_quantities

def register_CFunction_BOB_v2_setup_peak_attachment() -> Union[None, pcg.NRPyEnv_type]:
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
"""
    # --------------------------------------------------------
    # 1. Helper Function: Evaluate d|h|/dt
    # --------------------------------------------------------
    # This assumes t_p (News Peak) = 0.0 to find the relative lag.
    desc_helper = "Helper: Evaluates d|h|/dt assuming t_p_news=0."
    name_helper = "BOB_v2_strain_deriv_lag_helper"
    params_helper = "double t_val, void *params"
    
    wf = BOB_v2_waveform_quantities()
    
    body_helper = """
    commondata_struct *commondata = (commondata_struct *)params;
    const REAL m1 = commondata->m1;
    const REAL m2 = commondata->m2;
    const REAL chi1 = commondata->chi1;
    const REAL chi2 = commondata->chi2;
    const REAL omega_qnm = commondata->omega_qnm;
    const REAL tau_qnm = commondata->tau_qnm;
    const REAL M_f = commondata->M_f;
    const REAL a_f = commondata->a_f;
    
    // Assume News Peak is at 0
    const REAL t_p = 0.0; 
    const REAL t = t_val;
    
    """
    body_helper += ccg.c_codegen(
        [wf.strain_amp_deriv], ["REAL deriv_val"], 
        verbose=False, include_braces=False
    )
    body_helper += "return deriv_val;\n"
    
    cfc.register_CFunction(
        subdirectory="merger_waveform", includes=includes,prefunc=prefunc, desc=desc_helper,
        cfunc_type="double", name=name_helper, params=params_helper, body=body_helper
    )


    desc = """
    1. Finds the Lag (time delay) between News Peak and Strain Peak.
    2. Calculates |h|, omega, dot_omega at that peak for NQC.
    """
    cfunc_type = "int"
    name = "BOB_v2_setup_peak_attachment"
    params = "commondata_struct *restrict commondata"
    
    body = """
    // Step 1: Find the Lag
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    gsl_function F;
    F.function = &BOB_v2_strain_deriv_lag_helper;
    F.params = commondata;
    
    // Extremely conservative search window. Generally the strain peak will only be a few M before the news peak.
    double x_lo = -20.0;
    double x_hi = 20.0;
    
    if(gsl_root_fsolver_set(s, &F, x_lo, x_hi) != GSL_SUCCESS) {
        fprintf(stderr, "Error: BOB_v2 lag solver initialization failed.\\n");
        exit(EXIT_FAILURE);
    }
    
    int status;
    int iter = 0;
    double lag = 0.0;
    
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        lag = gsl_root_fsolver_root(s);
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lo, x_hi, 0, 1e-9);
    } while (status == GSL_CONTINUE && iter < 100);
    
    gsl_root_fsolver_free(s);

    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Warning: BOB_v2 lag solver did not converge. Using lag=0.\\n");
        lag = 0.0;
    }

    // Step 2: Set the News Peak Time
    // We set the time of the news peak such that the time of the strain peak equals t_attach.
    // In the construction of BOB, t_p_news is assumed to = 0. Here based on t_attach, we perform a time shift to allow for inspiral attachment.
    // We want: t_strain_peak = t_attach
    // We know: t_strain_peak = t_p_news + lag
    // Therefore: t_p_news = t_attach - lag
    commondata->t_p_BOB = commondata->t_attach - lag;
    return GSL_SUCCESS;
    """

    cfc.register_CFunction(
        subdirectory="merger_waveform",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )
    return pcg.NRPyEnv()