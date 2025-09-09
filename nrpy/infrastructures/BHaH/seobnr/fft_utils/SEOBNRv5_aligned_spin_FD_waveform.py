"""
Register CFunction for calculating the (2,2) mode of SEOBNRv5 in frequency domain.

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


def register_CFunction_SEOBNRv5_aligned_spin_FD_waveform() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for calculating the (2,2) mode of SEOBNRv5 in frequency domain.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Calculate the SEOBNRv5 22 mode in frequency domain."""
    cfunc_type = "void "
    prefunc = "#include<fftw3.h>"
    name = "SEOBNRv5_aligned_spin_FD_waveform"
    params = "const char *wisdom_file, commondata_struct *restrict commondata"
    body = r"""
// window and zero-pad the waveform
SEOBNRv5_aligned_spin_process_waveform(commondata);

size_t i;
commondata->nsteps_IMR_FD = commondata->nsteps_IMR;
const REAL dF = 1. / (commondata->nsteps_IMR_FD * commondata->dT);
commondata->waveform_IMR_FD = (double complex *)malloc(commondata->nsteps_IMR_FD * NUMMODES * sizeof(double complex));

fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * commondata->nsteps_IMR_FD);
fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * commondata->nsteps_IMR_FD);

if (!in || !out || !commondata->waveform_IMR_FD) {
  fprintf(stderr, "Error allocating memory for FFTW arrays.\n");
  exit(EXIT_FAILURE);
}

// Copy the processed TD waveform to the FFTW input array
for (i = 0; i < commondata->nsteps_IMR_FD; i++){
  in[i] = commondata->waveform_IMR[IDX_WF(i,STRAIN)];
}
// Load FFTW wisdom if available
if (wisdom_file) {
  FILE *file = fopen(wisdom_file, "r");
  if (file) {
    fftw_import_wisdom_from_file(file);
    fclose(file);
  } else {
    fprintf(stderr, "Could not open wisdom file '%s'.\n", wisdom_file);
  }
}

// Create FFTW plan
fftw_plan plan = fftw_plan_dft_1d((int)commondata->nsteps_IMR_FD, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

if (!plan) {
  fprintf(stderr, "Error creating FFTW plan.\n");
  fftw_free(in);
  fftw_free(out);
  return;
}

// Execute the FFT
fftw_execute(plan);

// Save FFTW wisdom for future runs
if (wisdom_file) {
  FILE *file = fopen(wisdom_file, "w");
  if (file) {
    fftw_export_wisdom_to_file(file);
    fclose(file);
  } else {
    fprintf(stderr, "Could not save wisdom to file '%s'.\n", wisdom_file);
  }
}

// Store the results (real and imaginary parts of the output)
for (i = 0; i < commondata->nsteps_IMR_FD; i++) {
  if (i <= commondata->nsteps_IMR_FD/2){
    commondata->waveform_IMR_FD[IDX_WF(i,FREQ)] = i * dF;
  }
  else{
    commondata->waveform_IMR_FD[IDX_WF(i,FREQ)] = ((REAL)i - (REAL)commondata->nsteps_IMR_FD) * dF;
  }
  commondata->waveform_IMR_FD[IDX_WF(i,STRAIN)] = out[i];
}
// Clean up
fftw_destroy_plan(plan);
fftw_free(in);
fftw_free(out);
fftw_cleanup();
"""
    cfc.register_CFunction(
        subdirectory="fft_utils",
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
