#include "BHaH_defines.h"
/*
 * Set commondata_struct to default values specified within NRPy+.
 */
void commondata_struct_set_to_default(commondata_struct *restrict commondata) {

  // Set commondata_struct variables to default
  commondata->CFL_FACTOR = 1.0;                   // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::CFL_FACTOR
  commondata->KreissOliger_strength_gauge = 0.99; // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::KreissOliger_strength_gauge
  commondata->KreissOliger_strength_nongauge =
      0.1;                                    // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::KreissOliger_strength_nongauge
  commondata->NUMGRIDS = 1;                   // nrpy.grid::NUMGRIDS
  commondata->TP_bare_mass_M = 0.5;           // TwoPunctures::TP_bare_mass_M
  commondata->TP_bare_mass_m = 0.5;           // TwoPunctures::TP_bare_mass_m
  commondata->TP_npoints_A = 48;              // TwoPunctures::TP_npoints_A
  commondata->TP_npoints_B = 48;              // TwoPunctures::TP_npoints_B
  commondata->TP_npoints_phi = 4;             // TwoPunctures::TP_npoints_phi
  commondata->bbhxy_BH_M_chix = 0.0;          // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_M_chix
  commondata->bbhxy_BH_M_chiy = 0;            // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_M_chiy
  commondata->bbhxy_BH_M_chiz = 0;            // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_M_chiz
  commondata->bbhxy_BH_m_chix = 0.0;          // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_m_chix
  commondata->bbhxy_BH_m_chiy = 0;            // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_m_chiy
  commondata->bbhxy_BH_m_chiz = 0;            // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_m_chiz
  commondata->checkpoint_every = 2.0;         // nrpy.infrastructures.BHaH.checkpointing::checkpoint_every
  commondata->convergence_factor = 1.0;       // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::convergence_factor
  commondata->diagnostics_output_every = 0.5; // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::diagnostics_output_every
  commondata->eta = 2.0;                      // nrpy.equations.general_relativity.BSSN_gauge_RHSs::eta
  commondata->initial_p_r = 0.0;              // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::initial_p_r
  commondata->initial_p_t = 0.0;              // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::initial_p_t
  commondata->initial_sep = 0.5;              // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::initial_sep
  commondata->mass_ratio = 1.0;               // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::mass_ratio
  commondata->swm2sh_maximum_l_mode_to_compute =
      2;                       // nrpy.infrastructures.BHaH.special_functions.spin_weight_minus2_spherical_harmonics::swm2sh_maximum_l_mode_to_compute
  commondata->t_final = 450.0; // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_final
  snprintf(commondata->outer_bc_type, 50, "radiation"); // nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions::outer_bc_type
}
