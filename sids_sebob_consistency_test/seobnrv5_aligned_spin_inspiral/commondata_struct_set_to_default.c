#include "BHaH_defines.h"
/**
 * Set commondata_struct to default values specified within NRPy+.
 */
void commondata_struct_set_to_default(commondata_struct *restrict commondata) {

  // Set commondata_struct variables to default
  commondata->Delta_t = 0.0;               // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::Delta_t
  commondata->M_f = 0.0;                   // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::M_f
  commondata->NUMGRIDS = 1;                // nrpy.grid::NUMGRIDS
  commondata->a6 = 0.0;                    // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::a6
  commondata->a_1_NQC = 0.0;               // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::a_1_NQC
  commondata->a_2_NQC = 0.0;               // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::a_2_NQC
  commondata->a_3_NQC = 0.0;               // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::a_3_NQC
  commondata->a_f = 0.0;                   // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::a_f
  commondata->b_1_NQC = 0.0;               // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::b_1_NQC
  commondata->b_2_NQC = 0.0;               // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::b_2_NQC
  commondata->chi1 = 0.4;                  // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::chi1
  commondata->chi2 = -0.3;                 // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::chi2
  commondata->dSO = 0.0;                   // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::dSO
  commondata->dT = 0.0;                    // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::dT
  commondata->dt = 2.4627455127717882e-05; // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::dt
  commondata->initial_omega = 0.01118;     // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::initial_omega
  commondata->m1 = 0.5;                    // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::m1
  commondata->m2 = 0.5;                    // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::m2
  commondata->mass_ratio = 1;              // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::mass_ratio
  commondata->nr_amp_1 = 0.0;              // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::nr_amp_1
  commondata->nr_amp_2 = 0.0;              // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::nr_amp_2
  commondata->nr_amp_3 = 0.0;              // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::nr_amp_3
  commondata->nr_omega_1 = 0.0;            // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::nr_omega_1
  commondata->nr_omega_2 = 0.0;            // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::nr_omega_2
  commondata->omega_qnm = 0.0;             // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::omega_qnm
  commondata->phi = 0;                     // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::phi
  commondata->pphi = 3.3;                  // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::pphi
  commondata->prstar = 0.0;                // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::prstar
  commondata->r = 20;                      // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::r
  commondata->r_ISCO = 0.0;                // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::r_ISCO
  commondata->r_stop = 0.0;                // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::r_stop
  commondata->t_ISCO = 0.0;                // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::t_ISCO
  commondata->t_attach = 0.0;              // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::t_attach
  commondata->t_stepback = 250.0;          // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::t_stepback
  commondata->tau_qnm = 0.0;               // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::tau_qnm
  commondata->total_mass = 50;             // nrpy.infrastructures.BHaH.seobnr.SEOBNR_C_codegen_library::total_mass
} // END FUNCTION commondata_struct_set_to_default
