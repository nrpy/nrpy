#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * -={ main() function }=-
 * Step 1.a: Set each commondata CodeParameter to default.
 * Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
 * Step 2: Compute Hamiltonian and derivatives.
 */
int main(int argc, const char *argv[]) {
  commondata_struct commondata; // commondata contains parameters common to all grids.
  // Step 0: Initialize a loop parameter for outputs
  //  Step 1.a: Set each commondata CodeParameter to default.
  commondata_struct_set_to_default(&commondata);
  // Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
  cmdline_input_and_parfile_parser(&commondata, argc, argv);
  // Step 1.c: Overwrite default values of m1, m2, a6, and dSO.
  SEOBNRv5_aligned_spin_coefficients(&commondata);
  // Step 2.a: Compute SEOBNRv5 conservative initial conditions.
  SEOBNRv5_aligned_spin_initial_conditions_conservative(&commondata);
  // Step 2.b: Print out the conservative initial conditions.
  // printf("r = %.15e\n",commondata.r);
  // printf("pphi = %.15e\n",commondata.pphi);
  // Step 3.a: Compute SEOBNRv5 dissipative initial conditions.
  SEOBNRv5_aligned_spin_initial_conditions_dissipative(&commondata);
  // Step 3.b: Print out the dissipative initial conditions.
  // printf("prstar = %.15e\n",commondata.prstar);
  // Step 4: Run the ODE integration.
  SEOBNRv5_aligned_spin_ode_integration(&commondata);
  // Step 5. Generate the waveform.
  SEOBNRv5_aligned_spin_waveform_from_dynamics(&commondata);
  // Step 6. Compute and apply the NQC corrections
  SEOBNRv5_aligned_spin_NQC_corrections(&commondata);
  // Step 7.a Compute the IMR waveform
  SEOBNRv5_aligned_spin_IMR_waveform(&commondata);

  // Step 6.b: Print the resulting waveform.
  for (size_t i = 0; i < commondata.nsteps_IMR; i++) {
    printf("%.15e %.15e %.15e\n", creal(commondata.waveform_IMR[IDX_WF(i, TIME)]), creal(commondata.waveform_IMR[IDX_WF(i, STRAIN)]),
           cimag(commondata.waveform_IMR[IDX_WF(i, STRAIN)]));
  }

  free(commondata.dynamics_low);
  free(commondata.dynamics_fine);
  free(commondata.waveform_low);
  free(commondata.waveform_fine);
  free(commondata.waveform_inspiral);
  free(commondata.waveform_IMR);

  return 0;
} // END FUNCTION main
