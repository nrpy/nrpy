"""
Register C function for BOB-informed merger waveform.

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
    desc = """
Calculates the BOB (2,2) mode for a given array of times.

@param times - Array of times at which to evaluate the waveform.
@param amps - Array to store the calculated amplitudes.
@param phases - Array to store the calculated phases.
@param nsteps_BOB - length of the times array.
@param commondata - Common data structure containing the model parameters."""
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
const REAL pshift = fabs(phases[0]);
for (i = 0; i < nsteps_BOB; i++){
  phases[i] = fabs(phases[i]) - pshift;
}
free(wrapped_phases);
"""
    cfc.register_CFunction(
        subdirectory="merger_waveform",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
