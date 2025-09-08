"""
Set up C function library for native SEOBNRv5 merger-related routines.

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
    desc = """
Calculates the (2,2) mode of the native SEOBNRv5 merger-ringdown model for a given array of times.

@params times - Array of times at which to evaluate the waveform.
@params amps - Array to store the calculated amplitudes.
@params phases - Array to store the calculated phases.
@params t_0 - Attachment time.
@params h_0 - Amplitude at attachment time.
@params hdot_0 - Amplitude derivative at attachment time.
@params phi_0 - Phase at attachment time.
@params phidot_0 - Angular frequency at attachment time.
@params nsteps_MR - length of the times array.
@params commondata - Common data structure containing the model parameters.
"""
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
