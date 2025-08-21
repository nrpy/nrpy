"""
Set up C function library for BOB related routines.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_merger_quantities as SEOBNRv5_mr
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_NQC_rhs() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for calculating the NQC amplitudes and phase from BOB.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = SEOBNRv5_mr.SEOBNRv5_aligned_spin_merger_quantities()
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    SEOBNRv5_code = ccg.c_codegen(
        [
            wf.h_t_attach,
            wf.hdot_t_attach,
            wf.hddot_t_attach,
            wf.w_t_attach,
            wf.wdot_t_attach,
        ],
        [
            "const REAL h_t_attach",
            "const REAL hdot_t_attach",
            "const REAL hddot_t_attach",
            "const REAL w_t_attach",
            "const REAL wdot_t_attach",
        ],
        verbose=False,
        include_braces=False,
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Calculate the SEOBNRv5 informed NQC amplitudes and phases."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_NQC_rhs"
    params = "commondata_struct *restrict commondata , REAL *restrict amps , REAL *restrict omegas"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
//compute
"""
    body += SEOBNRv5_code
    body += """
amps[0] = h_t_attach;
amps[1] = hdot_t_attach;
amps[2] = hddot_t_attach;
omegas[0] = fabs(w_t_attach);
omegas[1] = fabs(wdot_t_attach);
commondata->nr_amp_1 = amps[0];
commondata->nr_amp_2 = amps[1];
commondata->nr_amp_3 = amps[2];
commondata->nr_omega_1 = omegas[0];
commondata->nr_omega_2 = omegas[1];
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
    return pcg.NRPyEnv()


def register_CFunction_SEOBNRv5_aligned_spin_merger_waveform() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for calculating the (2,2) mode of the SEOBNRv5 merger.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = SEOBNRv5_mr.SEOBNRv5_aligned_spin_merger_quantities()
    SEOBNRv5_code = ccg.c_codegen(
        [wf.h, wf.phi],
        ["const REAL h", "const REAL phi"],
        verbose=False,
        include_braces=False,
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Calculate the BOB 22 mode."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_merger_waveform"
    params = "const REAL t , const REAL t_0, const REAL h_0 , const REAL hdot_0 , const REAL phi_0 , const REAL phidot_0, commondata_struct *restrict commondata , REAL *restrict waveform"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL omega_qnm = commondata->omega_qnm;
const REAL tau_qnm = commondata->tau_qnm;
//compute
"""
    body += SEOBNRv5_code
    body += """
waveform[0] = h;
waveform[1] = phi;
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
    return pcg.NRPyEnv()


def register_CFunction_SEOBNRv5_aligned_spin_merger_waveform_from_times() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for computing the (2,2) mode of the SEOBNRv5 merger ringdown waveform given the evaluation times.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Calculate the BOB 22 mode."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_merger_waveform_from_times"
    params = "REAL *restrict times , REAL *restrict amps , REAL *restrict phases , const REAL t_0, const REAL h_0 , const REAL hdot_0 , const REAL phi_0 , const REAL phidot_0, const size_t nsteps_MR , commondata_struct *restrict commondata"
    body = """
size_t i;
REAL waveform[2];
for (i = 0; i < nsteps_MR; i++) {
  //compute
  SEOBNRv5_aligned_spin_merger_waveform(times[i], t_0, h_0 , hdot_0 , phi_0 , phidot_0 , commondata , waveform);
  //store
  amps[i] = waveform[0];
  phases[i] = waveform[1];
}
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
    return pcg.NRPyEnv()
