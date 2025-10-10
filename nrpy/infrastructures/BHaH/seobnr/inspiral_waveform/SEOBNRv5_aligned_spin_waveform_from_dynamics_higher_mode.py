"""
Set up C function library for the SEOBNR aligned spin expressions.

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


def register_CFunction_SEOBNRv5_aligned_spin_waveform_from_dynamics() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for computing the (2,2) mode of SEOBNRv5 given the dynamics.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Calculates the (2,2) mode of the SEOBNRv5 inspiral waveform for the low- and fine-sampled ODE trajectory.

@param commondata - Common data structure containing the model parameters.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_waveform_from_dynamics"
    params = "commondata_struct *restrict commondata"
    body = """
int i;
REAL dynamics[NUMVARS];
COMPLEX *strain_modes_single_timestep = (COMPLEX *) malloc((NUMMODES- 1)*sizeof(COMPLEX));
if (strain_modes_single_timestep == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_waveform_from_dynamics(), malloc() failed to for strain_modes_single_timestep\\n");
  exit(1);
}
commondata->waveform_low = (double complex *)malloc(commondata->nsteps_low*NUMMODESSTORED*sizeof(double complex));
commondata->waveform_low = (double complex *)malloc(commondata->nsteps_low*NUMMODES*sizeof(double complex)); //t , h_+ , h_x
if (commondata->waveform_low == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_waveform_from_dynamics(), malloc() failed to for commondata->waveform_low\\n");
  exit(1);
}
commondata->waveform_fine = (double complex *)malloc(commondata->nsteps_fine*NUMMODESSTORED*sizeof(double complex));
commondata->waveform_fine = (double complex *)malloc(commondata->nsteps_fine*NUMMODES*sizeof(double complex)); //t , h_+ , h_x
if (commondata->waveform_fine == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_waveform_from_dynamics(), malloc() failed to for commondata->waveform_fine\\n");
  exit(1);
}

//low sampling
for (i = 0; i < commondata->nsteps_low; i++) {
  //assign
  dynamics[TIME] = commondata->dynamics_low[IDX(i,TIME)];
  dynamics[R] = commondata->dynamics_low[IDX(i,R)];
  dynamics[PHI] = commondata->dynamics_low[IDX(i,PHI)];
  dynamics[PRSTAR] = commondata->dynamics_low[IDX(i,PRSTAR)];
  dynamics[PPHI] = commondata->dynamics_low[IDX(i,PPHI)];
  dynamics[H] = commondata->dynamics_low[IDX(i,H)];
  dynamics[OMEGA] = commondata->dynamics_low[IDX(i,OMEGA)];
  dynamics[OMEGA_CIRC] = commondata->dynamics_low[IDX(i,OMEGA_CIRC)];

  //compute
  //store
  SEOBNRv5_aligned_spin_waveform(dynamics, commondata, strain_modes_single_timestep);
  commondata->waveform_low[IDX_WF(i,TIME)] = dynamics[TIME];
  commondata->waveform_low[IDX_WF(i,STRAIN)] = strain_modes_single_timestep[STRAIN22-1];
}
//high sampling
for (i = 0; i < commondata->nsteps_fine; i++) {
  //assign
  dynamics[TIME] = commondata->dynamics_fine[IDX(i,TIME)];
  dynamics[R] = commondata->dynamics_fine[IDX(i,R)];
  dynamics[PHI] = commondata->dynamics_fine[IDX(i,PHI)];
  dynamics[PRSTAR] = commondata->dynamics_fine[IDX(i,PRSTAR)];
  dynamics[PPHI] = commondata->dynamics_fine[IDX(i,PPHI)];
  dynamics[H] = commondata->dynamics_fine[IDX(i,H)];
  dynamics[OMEGA] = commondata->dynamics_fine[IDX(i,OMEGA)];
  dynamics[OMEGA_CIRC] = commondata->dynamics_fine[IDX(i,OMEGA_CIRC)];

  //compute
  //store
  SEOBNRv5_aligned_spin_waveform(dynamics, commondata, strain_modes_single_timestep);
  commondata->waveform_fine[IDX_WF(i,TIME)] = dynamics[TIME];
  commondata->waveform_fine[IDX_WF(i,STRAIN)] = strain_modes_single_timestep[STRAIN22-1];
  

}
free(strain_modes_single_timestep);
"""
    cfc.register_CFunction(
        subdirectory="inspiral_waveform",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
