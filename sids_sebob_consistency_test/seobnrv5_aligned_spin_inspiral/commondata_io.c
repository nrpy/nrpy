#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Store the commondata struct in a binary file.
 */
int commondata_io(commondata_struct *restrict data, const char *restrict filename) {

  FILE *fp = fopen(filename, "wb");
  if (fp == NULL) {
    perror("Error opening file for writing");
    return 1;
  }

  // Write non-pointer members
  fwrite(&data->Delta_t, sizeof(REAL), 1, fp);
  fwrite(&data->Hreal, sizeof(REAL), 1, fp);
  fwrite(&data->M_f, sizeof(REAL), 1, fp);
  fwrite(&data->Omega_circ, sizeof(REAL), 1, fp);
  fwrite(&data->a6, sizeof(REAL), 1, fp);
  fwrite(&data->a_1_NQC, sizeof(REAL), 1, fp);
  fwrite(&data->a_2_NQC, sizeof(REAL), 1, fp);
  fwrite(&data->a_3_NQC, sizeof(REAL), 1, fp);
  fwrite(&data->nr_amp_1, sizeof(REAL), 1, fp);
  fwrite(&data->nr_amp_2, sizeof(REAL), 1, fp);
  fwrite(&data->nr_amp_3, sizeof(REAL), 1, fp);
  fwrite(&data->a_f, sizeof(REAL), 1, fp);
  fwrite(&data->b_1_NQC, sizeof(REAL), 1, fp);
  fwrite(&data->b_2_NQC, sizeof(REAL), 1, fp);
  fwrite(&data->nr_omega_1, sizeof(REAL), 1, fp);
  fwrite(&data->nr_omega_2, sizeof(REAL), 1, fp);
  fwrite(&data->chi1, sizeof(REAL), 1, fp);
  fwrite(&data->chi2, sizeof(REAL), 1, fp);
  fwrite(&data->dHreal_dpphi, sizeof(REAL), 1, fp);
  fwrite(&data->dHreal_dpphi_circ, sizeof(REAL), 1, fp);
  fwrite(&data->dHreal_dprstar, sizeof(REAL), 1, fp);
  fwrite(&data->dHreal_dr, sizeof(REAL), 1, fp);
  fwrite(&data->dHreal_dr_circ, sizeof(REAL), 1, fp);
  fwrite(&data->dHreal_dr_dpphi, sizeof(REAL), 1, fp);
  fwrite(&data->dHreal_dr_dr, sizeof(REAL), 1, fp);
  fwrite(&data->dSO, sizeof(REAL), 1, fp);
  fwrite(&data->dT, sizeof(REAL), 1, fp);
  fwrite(&data->dt, sizeof(REAL), 1, fp);
  fwrite(&data->flux, sizeof(REAL), 1, fp);
  fwrite(&data->initial_omega, sizeof(REAL), 1, fp);
  fwrite(&data->m1, sizeof(REAL), 1, fp);
  fwrite(&data->m2, sizeof(REAL), 1, fp);
  fwrite(&data->mass_ratio, sizeof(REAL), 1, fp);
  fwrite(&data->omega_qnm, sizeof(REAL), 1, fp);
  fwrite(&data->phi, sizeof(REAL), 1, fp);
  fwrite(&data->pphi, sizeof(REAL), 1, fp);
  fwrite(&data->prstar, sizeof(REAL), 1, fp);
  fwrite(&data->r, sizeof(REAL), 1, fp);
  fwrite(&data->r_ISCO, sizeof(REAL), 1, fp);
  fwrite(&data->r_stop, sizeof(REAL), 1, fp);
  fwrite(&data->t_ISCO, sizeof(REAL), 1, fp);
  fwrite(&data->t_attach, sizeof(REAL), 1, fp);
  fwrite(&data->t_stepback, sizeof(REAL), 1, fp);
  fwrite(&data->tau_qnm, sizeof(REAL), 1, fp);
  fwrite(&data->total_mass, sizeof(REAL), 1, fp);
  fwrite(&data->xi, sizeof(REAL), 1, fp);
  fwrite(&data->NUMGRIDS, sizeof(int), 1, fp);
  fwrite(&data->nsteps_IMR, sizeof(size_t), 1, fp);
  fwrite(&data->nsteps_fine, sizeof(size_t), 1, fp);
  fwrite(&data->nsteps_inspiral, sizeof(size_t), 1, fp);
  fwrite(&data->nsteps_low, sizeof(size_t), 1, fp);
  fwrite(&data->nsteps_raw, sizeof(size_t), 1, fp);
  int num_dynamics_vars = 8;

  // Write pointer members (arrays)
  if (data->dynamics_fine != NULL) {
    fwrite(data->dynamics_fine, sizeof(REAL), data->nsteps_fine * num_dynamics_vars, fp);
  }
  if (data->dynamics_low != NULL) {
    fwrite(data->dynamics_low, sizeof(REAL), data->nsteps_low * num_dynamics_vars, fp);
  }
  if (data->dynamics_raw != NULL) {
    fwrite(data->dynamics_raw, sizeof(REAL), data->nsteps_raw * num_dynamics_vars, fp);
  }
  if (data->waveform_IMR != NULL) {
    fwrite(data->waveform_IMR, sizeof(double complex), data->nsteps_IMR * 2, fp);
  }
  if (data->waveform_fine != NULL) {
    fwrite(data->waveform_fine, sizeof(double complex), data->nsteps_fine * 2, fp);
  }
  if (data->waveform_inspiral != NULL) {
    fwrite(data->waveform_inspiral, sizeof(double complex), data->nsteps_inspiral * 2, fp);
  }
  if (data->waveform_low != NULL) {
    fwrite(data->waveform_low, sizeof(double complex), data->nsteps_low * 2, fp);
  }

  fclose(fp);
  return GSL_SUCCESS;
} // END FUNCTION commondata_io
