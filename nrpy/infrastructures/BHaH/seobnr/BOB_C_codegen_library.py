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
import nrpy.equations.seobnr.BOB_aligned_spin_waveform_quantities as BOB_wf
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_BOB_aligned_spin_NQC_rhs() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for calculating the NQC amplitudes and phase from BOB.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = BOB_wf.BOB_aligned_spin_waveform_quantities()
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    BOB_code = ccg.c_codegen(
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
    desc = """Calculate the BOB informed NQC amplitudes and phases."""
    cfunc_type = "void"
    name = "BOB_aligned_spin_NQC_rhs"
    params = "commondata_struct *restrict commondata , REAL *restrict amps , REAL *restrict omegas"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL omega_qnm = commondata->omega_qnm;
const REAL tau_qnm = commondata->tau_qnm;
//compute
"""
    body += BOB_code
    body += """
amps[0] = h_t_attach;
amps[1] = hdot_t_attach;
amps[2] = hddot_t_attach;
omegas[0] = w_t_attach;
omegas[1] = wdot_t_attach;
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
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_BOB_aligned_spin_waveform() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for calculating the (2,2) mode of BOB.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = BOB_wf.BOB_aligned_spin_waveform_quantities()
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    BOB_code = ccg.c_codegen(
        [wf.h, wf.phi],
        ["const REAL h", "const REAL phi"],
        verbose=False,
        include_braces=False,
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Calculate the BOB 22 mode."""
    cfunc_type = "void"
    name = "BOB_aligned_spin_waveform"
    params = "const REAL t , commondata_struct *restrict commondata , REAL *restrict waveform"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL omega_qnm = commondata->omega_qnm;
const REAL tau_qnm = commondata->tau_qnm;
const REAL t_0 = commondata->t_attach;
//compute
"""
    body += BOB_code
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
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_BOB_aligned_spin_waveform_from_times() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for computing the (2,2) mode of BOB given the evaluation times.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Calculate the BOB 22 mode."""
    cfunc_type = "void"
    name = "BOB_aligned_spin_waveform_from_times"
    params = "REAL *restrict times , REAL *restrict amps , REAL *restrict phases , const size_t nsteps_BOB , commondata_struct *restrict commondata"
    body = """
size_t i;
REAL *restrict wrapped_phases = (REAL *)malloc(nsteps_BOB*sizeof(REAL));
REAL waveform[2];
for (i = 0; i < nsteps_BOB; i++) {
  //compute
  BOB_aligned_spin_waveform(times[i], commondata , waveform);
  //store
  amps[i] = waveform[0];
  wrapped_phases[i] = waveform[1];
}
//unwrap the phase
SEOBNRv5_aligned_spin_unwrap(wrapped_phases,phases,nsteps_BOB);

//shift and take absolute value of the phase
const REAL pshift = phases[0];
for (i = 0; i < nsteps_BOB; i++){
  phases[i] = fabs(phases[i] - pshift);
}
free(wrapped_phases);
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
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
