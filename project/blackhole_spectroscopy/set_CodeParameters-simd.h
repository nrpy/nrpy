const REAL NOSIMDAMPL = params->AMPL;               // nrpy.reference_metric_SinhSpherical::AMPL
const REAL_SIMD_ARRAY AMPL = ConstSIMD(NOSIMDAMPL); // nrpy.reference_metric_SinhSpherical::AMPL
const REAL NOSIMDbbhxy_BH_M_chix =
    commondata->bbhxy_BH_M_chix; // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_M_chix
const REAL_SIMD_ARRAY bbhxy_BH_M_chix =
    ConstSIMD(NOSIMDbbhxy_BH_M_chix); // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_M_chix
const REAL NOSIMDbbhxy_BH_m_chix =
    commondata->bbhxy_BH_m_chix; // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_m_chix
const REAL_SIMD_ARRAY bbhxy_BH_m_chix =
    ConstSIMD(NOSIMDbbhxy_BH_m_chix); // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_m_chix
const REAL NOSIMDbbhxy_BH_M_chiy =
    commondata->bbhxy_BH_M_chiy; // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_M_chiy
const REAL_SIMD_ARRAY bbhxy_BH_M_chiy =
    ConstSIMD(NOSIMDbbhxy_BH_M_chiy); // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_M_chiy
const REAL NOSIMDbbhxy_BH_m_chiy =
    commondata->bbhxy_BH_m_chiy; // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_m_chiy
const REAL_SIMD_ARRAY bbhxy_BH_m_chiy =
    ConstSIMD(NOSIMDbbhxy_BH_m_chiy); // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_m_chiy
const REAL NOSIMDbbhxy_BH_M_chiz =
    commondata->bbhxy_BH_M_chiz; // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_M_chiz
const REAL_SIMD_ARRAY bbhxy_BH_M_chiz =
    ConstSIMD(NOSIMDbbhxy_BH_M_chiz); // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_M_chiz
const REAL NOSIMDbbhxy_BH_m_chiz =
    commondata->bbhxy_BH_m_chiz; // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_m_chiz
const REAL_SIMD_ARRAY bbhxy_BH_m_chiz =
    ConstSIMD(NOSIMDbbhxy_BH_m_chiz);                 // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::bbhxy_BH_m_chiz
const REAL NOSIMDCart_originx = params->Cart_originx; // nrpy.grid::Cart_originx
const REAL_SIMD_ARRAY Cart_originx = ConstSIMD(NOSIMDCart_originx);         // nrpy.grid::Cart_originx
const REAL NOSIMDCart_originy = params->Cart_originy;                       // nrpy.grid::Cart_originy
const REAL_SIMD_ARRAY Cart_originy = ConstSIMD(NOSIMDCart_originy);         // nrpy.grid::Cart_originy
const REAL NOSIMDCart_originz = params->Cart_originz;                       // nrpy.grid::Cart_originz
const REAL_SIMD_ARRAY Cart_originz = ConstSIMD(NOSIMDCart_originz);         // nrpy.grid::Cart_originz
const REAL NOSIMDCFL_FACTOR = commondata->CFL_FACTOR;                       // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::CFL_FACTOR
const REAL_SIMD_ARRAY CFL_FACTOR = ConstSIMD(NOSIMDCFL_FACTOR);             // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::CFL_FACTOR
const REAL NOSIMDcheckpoint_every = commondata->checkpoint_every;           // nrpy.infrastructures.BHaH.checkpointing::checkpoint_every
const REAL_SIMD_ARRAY checkpoint_every = ConstSIMD(NOSIMDcheckpoint_every); // nrpy.infrastructures.BHaH.checkpointing::checkpoint_every
const REAL NOSIMDconvergence_factor = commondata->convergence_factor; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::convergence_factor
const REAL_SIMD_ARRAY convergence_factor =
    ConstSIMD(NOSIMDconvergence_factor);               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::convergence_factor
const int CoordSystem_hash = params->CoordSystem_hash; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::CoordSystem_hash
const REAL NOSIMDdiagnostics_output_every =
    commondata->diagnostics_output_every; // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::diagnostics_output_every
const REAL_SIMD_ARRAY diagnostics_output_every =
    ConstSIMD(NOSIMDdiagnostics_output_every);      // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::diagnostics_output_every
const REAL NOSIMDdt = commondata->dt;               // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::dt
const REAL_SIMD_ARRAY dt = ConstSIMD(NOSIMDdt);     // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::dt
const REAL NOSIMDdxx0 = params->dxx0;               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx0
const REAL_SIMD_ARRAY dxx0 = ConstSIMD(NOSIMDdxx0); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx0
const REAL NOSIMDdxx1 = params->dxx1;               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx1
const REAL_SIMD_ARRAY dxx1 = ConstSIMD(NOSIMDdxx1); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx1
const REAL NOSIMDdxx2 = params->dxx2;               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx2
const REAL_SIMD_ARRAY dxx2 = ConstSIMD(NOSIMDdxx2); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx2
const REAL NOSIMDeta = commondata->eta;             // nrpy.equations.general_relativity.BSSN_gauge_RHSs::eta
const REAL_SIMD_ARRAY eta = ConstSIMD(NOSIMDeta);   // nrpy.equations.general_relativity.BSSN_gauge_RHSs::eta
const REAL NOSIMDf0_of_xx0 = params->f0_of_xx0;     // nrpy.reference_metric_SinhSpherical::f0_of_xx0
const REAL_SIMD_ARRAY f0_of_xx0 = ConstSIMD(NOSIMDf0_of_xx0);                   // nrpy.reference_metric_SinhSpherical::f0_of_xx0
const REAL NOSIMDf1_of_xx1 = params->f1_of_xx1;                                 // nrpy.reference_metric_SinhSpherical::f1_of_xx1
const REAL_SIMD_ARRAY f1_of_xx1 = ConstSIMD(NOSIMDf1_of_xx1);                   // nrpy.reference_metric_SinhSpherical::f1_of_xx1
const REAL NOSIMDf2_of_xx0 = params->f2_of_xx0;                                 // nrpy.reference_metric_SinhSpherical::f2_of_xx0
const REAL_SIMD_ARRAY f2_of_xx0 = ConstSIMD(NOSIMDf2_of_xx0);                   // nrpy.reference_metric_SinhSpherical::f2_of_xx0
const REAL NOSIMDf3_of_xx2 = params->f3_of_xx2;                                 // nrpy.reference_metric_SinhSpherical::f3_of_xx2
const REAL_SIMD_ARRAY f3_of_xx2 = ConstSIMD(NOSIMDf3_of_xx2);                   // nrpy.reference_metric_SinhSpherical::f3_of_xx2
const REAL NOSIMDf4_of_xx1 = params->f4_of_xx1;                                 // nrpy.reference_metric_SinhSpherical::f4_of_xx1
const REAL_SIMD_ARRAY f4_of_xx1 = ConstSIMD(NOSIMDf4_of_xx1);                   // nrpy.reference_metric_SinhSpherical::f4_of_xx1
const REAL NOSIMDgrid_physical_size = params->grid_physical_size;               // nrpy.reference_metric::grid_physical_size
const REAL_SIMD_ARRAY grid_physical_size = ConstSIMD(NOSIMDgrid_physical_size); // nrpy.reference_metric::grid_physical_size
const REAL NOSIMDinitial_p_r = commondata->initial_p_r; // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::initial_p_r
const REAL_SIMD_ARRAY initial_p_r =
    ConstSIMD(NOSIMDinitial_p_r);                       // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::initial_p_r
const REAL NOSIMDinitial_p_t = commondata->initial_p_t; // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::initial_p_t
const REAL_SIMD_ARRAY initial_p_t =
    ConstSIMD(NOSIMDinitial_p_t);                       // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::initial_p_t
const REAL NOSIMDinitial_sep = commondata->initial_sep; // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::initial_sep
const REAL_SIMD_ARRAY initial_sep =
    ConstSIMD(NOSIMDinitial_sep);                         // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::initial_sep
const REAL NOSIMDinvdxx0 = params->invdxx0;               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx0
const REAL_SIMD_ARRAY invdxx0 = ConstSIMD(NOSIMDinvdxx0); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx0
const REAL NOSIMDinvdxx1 = params->invdxx1;               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx1
const REAL_SIMD_ARRAY invdxx1 = ConstSIMD(NOSIMDinvdxx1); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx1
const REAL NOSIMDinvdxx2 = params->invdxx2;               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx2
const REAL_SIMD_ARRAY invdxx2 = ConstSIMD(NOSIMDinvdxx2); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx2
const REAL NOSIMDKreissOliger_strength_gauge =
    commondata->KreissOliger_strength_gauge; // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::KreissOliger_strength_gauge
const REAL_SIMD_ARRAY KreissOliger_strength_gauge =
    ConstSIMD(NOSIMDKreissOliger_strength_gauge); // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::KreissOliger_strength_gauge
const REAL NOSIMDKreissOliger_strength_nongauge =
    commondata->KreissOliger_strength_nongauge; // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::KreissOliger_strength_nongauge
const REAL_SIMD_ARRAY KreissOliger_strength_nongauge = ConstSIMD(
    NOSIMDKreissOliger_strength_nongauge);    // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::KreissOliger_strength_nongauge
const REAL NOSIMDmass_M = commondata->mass_M; // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::mass_M
const REAL_SIMD_ARRAY mass_M = ConstSIMD(NOSIMDmass_M); // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::mass_M
const REAL NOSIMDmass_m = commondata->mass_m;           // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::mass_m
const REAL_SIMD_ARRAY mass_m = ConstSIMD(NOSIMDmass_m); // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::mass_m
const REAL NOSIMDmass_ratio = commondata->mass_ratio;   // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::mass_ratio
const REAL_SIMD_ARRAY mass_ratio =
    ConstSIMD(NOSIMDmass_ratio);                           // nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta::mass_ratio
const int nn = commondata->nn;                             // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::nn
const int nn_0 = commondata->nn_0;                         // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::nn_0
const int NUMGRIDS = commondata->NUMGRIDS;                 // nrpy.grid::NUMGRIDS
const int Nxx0 = params->Nxx0;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx0
const int Nxx1 = params->Nxx1;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx1
const int Nxx2 = params->Nxx2;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx2
const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS0
const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS1
const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS2
const REAL NOSIMDSINHW = params->SINHW;                    // nrpy.reference_metric_SinhSpherical::SINHW
const REAL_SIMD_ARRAY SINHW = ConstSIMD(NOSIMDSINHW);      // nrpy.reference_metric_SinhSpherical::SINHW
const int swm2sh_maximum_l_mode_to_compute =
    commondata
        ->swm2sh_maximum_l_mode_to_compute; // nrpy.infrastructures.BHaH.special_functions.spin_weight_minus2_spherical_harmonics::swm2sh_maximum_l_mode_to_compute
const REAL NOSIMDt_0 = commondata->t_0;                                 // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_0
const REAL_SIMD_ARRAY t_0 = ConstSIMD(NOSIMDt_0);                       // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_0
const REAL NOSIMDt_final = commondata->t_final;                         // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_final
const REAL_SIMD_ARRAY t_final = ConstSIMD(NOSIMDt_final);               // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_final
const REAL NOSIMDtime = commondata->time;                               // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::time
const REAL_SIMD_ARRAY time = ConstSIMD(NOSIMDtime);                     // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::time
const REAL NOSIMDTP_bare_mass_m = commondata->TP_bare_mass_m;           // TwoPunctures::TP_bare_mass_m
const REAL_SIMD_ARRAY TP_bare_mass_m = ConstSIMD(NOSIMDTP_bare_mass_m); // TwoPunctures::TP_bare_mass_m
const REAL NOSIMDTP_bare_mass_M = commondata->TP_bare_mass_M;           // TwoPunctures::TP_bare_mass_M
const REAL_SIMD_ARRAY TP_bare_mass_M = ConstSIMD(NOSIMDTP_bare_mass_M); // TwoPunctures::TP_bare_mass_M
const int TP_npoints_A = commondata->TP_npoints_A;                      // TwoPunctures::TP_npoints_A
const int TP_npoints_B = commondata->TP_npoints_B;                      // TwoPunctures::TP_npoints_B
const int TP_npoints_phi = commondata->TP_npoints_phi;                  // TwoPunctures::TP_npoints_phi
const REAL NOSIMDxxmax0 = params->xxmax0;                               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax0
const REAL_SIMD_ARRAY xxmax0 = ConstSIMD(NOSIMDxxmax0);                 // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax0
const REAL NOSIMDxxmax1 = params->xxmax1;                               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax1
const REAL_SIMD_ARRAY xxmax1 = ConstSIMD(NOSIMDxxmax1);                 // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax1
const REAL NOSIMDxxmax2 = params->xxmax2;                               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax2
const REAL_SIMD_ARRAY xxmax2 = ConstSIMD(NOSIMDxxmax2);                 // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax2
const REAL NOSIMDxxmin0 = params->xxmin0;                               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin0
const REAL_SIMD_ARRAY xxmin0 = ConstSIMD(NOSIMDxxmin0);                 // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin0
const REAL NOSIMDxxmin1 = params->xxmin1;                               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin1
const REAL_SIMD_ARRAY xxmin1 = ConstSIMD(NOSIMDxxmin1);                 // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin1
const REAL NOSIMDxxmin2 = params->xxmin2;                               // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin2
const REAL_SIMD_ARRAY xxmin2 = ConstSIMD(NOSIMDxxmin2);                 // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin2
