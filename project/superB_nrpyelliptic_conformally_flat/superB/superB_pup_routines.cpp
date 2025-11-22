#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

// PUP routine for struct commondata_struct
void pup_commondata_struct(PUP::er &p, commondata_struct &commondata) {
  // PUP commondata struct
  PUParray(p, commondata.outer_bc_type, 50);   // nrpy.infrastructures.BHaH.CurviBoundaryConditions.BHaH_defines::outer_bc_type
  p | commondata.CFL_FACTOR;                   // nrpy.infrastructures.superB.MoL::CFL_FACTOR
  p | commondata.MINIMUM_GLOBAL_WAVESPEED;     // nrpy.infrastructures.BHaH.nrpyelliptic.auxevol_gfs_set_to_constant::MINIMUM_GLOBAL_WAVESPEED
  p | commondata.NUMGRIDS;                     // nrpy.grid::NUMGRIDS
  p | commondata.Nchare0;                      // nrpy.infrastructures.superB.main_chare::Nchare0
  p | commondata.Nchare1;                      // nrpy.infrastructures.superB.main_chare::Nchare1
  p | commondata.Nchare2;                      // nrpy.infrastructures.superB.main_chare::Nchare2
  p | commondata.P0_x;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::P0_x
  p | commondata.P0_y;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::P0_y
  p | commondata.P0_z;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::P0_z
  p | commondata.P1_x;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::P1_x
  p | commondata.P1_y;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::P1_y
  p | commondata.P1_z;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::P1_z
  p | commondata.S0_x;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::S0_x
  p | commondata.S0_y;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::S0_y
  p | commondata.S0_z;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::S0_z
  p | commondata.S1_x;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::S1_x
  p | commondata.S1_y;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::S1_y
  p | commondata.S1_z;                         // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::S1_z
  p | commondata.bare_mass_0;                  // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::bare_mass_0
  p | commondata.bare_mass_1;                  // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::bare_mass_1
  p | commondata.convergence_factor;           // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::convergence_factor
  p | commondata.diagnostics_output_every;     // nrpy.infrastructures.BHaH.diagnostics.diagnostics::diagnostics_output_every
  p | commondata.dt;                           // nrpy.infrastructures.superB.MoL::dt
  p | commondata.eta_damping;                  // nrpy.equations.nrpyelliptic.ConformallyFlat_RHSs::eta_damping
  p | commondata.log10_current_residual;       // nrpy.infrastructures.BHaH.nrpyelliptic.stop_conditions_check::log10_current_residual
  p | commondata.log10_residual_tolerance;     // nrpy.infrastructures.BHaH.nrpyelliptic.stop_conditions_check::log10_residual_tolerance
  p | commondata.nn;                           // nrpy.infrastructures.superB.MoL::nn
  p | commondata.nn_0;                         // nrpy.infrastructures.superB.MoL::nn_0
  p | commondata.nn_max;                       // nrpy.infrastructures.BHaH.nrpyelliptic.stop_conditions_check::nn_max
  p | commondata.output_progress_every;        // nrpy.infrastructures.BHaH.diagnostics.progress_indicator::output_progress_every
  p | commondata.start_wallclock_time.tv_nsec; // nrpy.infrastructures.BHaH.diagnostics.progress_indicator::start_wallclock_time
  p | commondata.start_wallclock_time.tv_sec;  // nrpy.infrastructures.BHaH.diagnostics.progress_indicator::start_wallclock_time
  p | commondata.stop_relaxation;              // nrpy.infrastructures.BHaH.nrpyelliptic.stop_conditions_check::stop_relaxation
  p | commondata.t_0;                          // nrpy.infrastructures.superB.MoL::t_0
  p | commondata.t_final;                      // nrpy.infrastructures.superB.MoL::t_final
  p | commondata.time;                         // nrpy.infrastructures.superB.MoL::time
  p | commondata.zPunc;                        // nrpy.equations.nrpyelliptic.ConformallyFlat_SourceTerms::zPunc
}
// PUP routine for struct params_struct
void pup_params_struct(PUP::er &p, params_struct &params) {
  // PUP params struct
  PUParray(p, params.CoordSystemName, 100); // nrpy.reference_metric::CoordSystemName
  PUParray(p, params.gridname, 100);        // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::gridname
  p | params.AMPL;                          // nrpy.reference_metric::AMPL
  p | params.Cart_originx;                  // nrpy.grid::Cart_originx
  p | params.Cart_originy;                  // nrpy.grid::Cart_originy
  p | params.Cart_originz;                  // nrpy.grid::Cart_originz
  p | params.CoordSystem_hash;              // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::CoordSystem_hash
  p | params.Nxx0;                          // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx0
  p | params.Nxx1;                          // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx1
  p | params.Nxx2;                          // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx2
  p | params.Nxx_plus_2NGHOSTS0;            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS0
  p | params.Nxx_plus_2NGHOSTS1;            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS1
  p | params.Nxx_plus_2NGHOSTS2;            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS2
  p | params.PI;                            // nrpy.reference_metric::PI
  p | params.SINHW;                         // nrpy.reference_metric::SINHW
  p | params.dxx0;                          // nrpy.infrastructures.BHaH.diagnostics.sqrt_detgammahat_d3xx_volume_element::dxx0
  p | params.dxx1;                          // nrpy.infrastructures.BHaH.diagnostics.sqrt_detgammahat_d3xx_volume_element::dxx1
  p | params.dxx2;                          // nrpy.infrastructures.BHaH.diagnostics.sqrt_detgammahat_d3xx_volume_element::dxx2
  p | params.grid_hole_radius;              // nrpy.reference_metric::grid_hole_radius
  p | params.grid_idx;                      // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::grid_idx
  p | params.grid_physical_size;            // nrpy.reference_metric::grid_physical_size
  p | params.grid_rotates;                  // nrpy.grid::grid_rotates
  p | params.invdxx0;                       // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx0
  p | params.invdxx1;                       // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx1
  p | params.invdxx2;                       // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx2
  p | params.xxmax0;                        // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax0
  p | params.xxmax1;                        // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax1
  p | params.xxmax2;                        // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax2
  p | params.xxmin0;                        // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin0
  p | params.xxmin1;                        // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin1
  p | params.xxmin2;                        // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin2
}

// PUP routine for struct innerpt_bc_struct
void pup_innerpt_bc_struct(PUP::er &p, innerpt_bc_struct &ibc) {
  p | ibc.dstpt;
  p | ibc.srcpt;
  for (int i = 0; i < 10; i++) {
    p | ibc.parity[i];
  }
}

// PUP routine for struct outerpt_bc_struct
void pup_outerpt_bc_struct(PUP::er &p, outerpt_bc_struct &obc) {
  p | obc.i0;
  p | obc.i1;
  p | obc.i2;
  p | obc.FACEX0;
  p | obc.FACEX1;
  p | obc.FACEX2;
}

// PUP routine for struct bc_info_struct
void pup_bc_info_struct(PUP::er &p, bc_info_struct &bci) {
  p | bci.num_inner_boundary_points;
  p | bci.num_inner_boundary_points_nonlocal;
  for (int i = 0; i < NGHOSTS; i++) {
    for (int j = 0; j < 3; j++) {
      p | bci.num_pure_outer_boundary_points[i][j];
    }
  }
  for (int i = 0; i < NGHOSTS; i++) {
    for (int j = 0; j < 6; j++) {
      for (int k = 0; k < 6; k++) {
        p | bci.bc_loop_bounds[i][j][k];
      }
    }
  }
}

// PUP routine for struct bc_struct
void pup_bc_struct(PUP::er &p, bc_struct &bc) {

  pup_bc_info_struct(p, bc.bc_info);

  if (p.isUnpacking()) {

    bc.inner_bc_array = (innerpt_bc_struct *restrict)malloc(sizeof(innerpt_bc_struct) * bc.bc_info.num_inner_boundary_points);
    bc.inner_bc_array_nonlocal = (innerpt_bc_struct *restrict)malloc(sizeof(innerpt_bc_struct) * bc.bc_info.num_inner_boundary_points_nonlocal);

    for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      for (int dirn = 0; dirn < 3; dirn++) {
        bc.pure_outer_bc_array[dirn + (3 * which_gz)] =
            (outerpt_bc_struct *restrict)malloc(bc.bc_info.num_pure_outer_boundary_points[which_gz][dirn] * sizeof(outerpt_bc_struct));
      }
    }
  }

  for (int i = 0; i < bc.bc_info.num_inner_boundary_points; i++) {
    pup_innerpt_bc_struct(p, bc.inner_bc_array[i]);
  }

  for (int i = 0; i < bc.bc_info.num_inner_boundary_points_nonlocal; i++) {
    pup_innerpt_bc_struct(p, bc.inner_bc_array_nonlocal[i]);
  }

  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    for (int dirn = 0; dirn < 3; dirn++) {
      int num_points = bc.bc_info.num_pure_outer_boundary_points[which_gz][dirn];
      for (int pt = 0; pt < num_points; pt++) {
        pup_outerpt_bc_struct(p, bc.pure_outer_bc_array[dirn + (3 * which_gz)][pt]);
      }
    }
  }
}

// PUP routine for struct MoL_gridfunctions_struct
void pup_MoL_gridfunctions_struct(PUP::er &p, MoL_gridfunctions_struct &gridfuncs, const params_struct &params, const commondata_struct &commondata) {
  p | gridfuncs.num_evol_gfs_to_sync;
  p | gridfuncs.num_auxevol_gfs_to_sync;
  p | gridfuncs.num_aux_gfs_to_sync;
  p | gridfuncs.max_sync_gfs;
  PUParray(p, gridfuncs.evol_gfs_to_sync, gridfuncs.num_evol_gfs_to_sync);
  PUParray(p, gridfuncs.auxevol_gfs_to_sync, gridfuncs.num_auxevol_gfs_to_sync);
  PUParray(p, gridfuncs.aux_gfs_to_sync, gridfuncs.num_aux_gfs_to_sync);

  const int Nxx_plus_2NGHOSTS_tot = params.Nxx_plus_2NGHOSTS0 * params.Nxx_plus_2NGHOSTS1 * params.Nxx_plus_2NGHOSTS2;
  if (p.isUnpacking()) {
    gridfuncs.y_n_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    gridfuncs.y_nplus1_running_total_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    gridfuncs.k_odd_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    gridfuncs.k_even_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    if (NUM_AUXEVOL_GFS > 0)
      gridfuncs.auxevol_gfs = (REAL *restrict)malloc(sizeof(REAL) * NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot);

    initialize_yn_and_non_yn_gfs_to_nan(&commondata, &params, &gridfuncs);
  }
  PUParray(p, gridfuncs.y_n_gfs, NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);

  if (strstr(params.CoordSystemName, "Spherical") != NULL) {
    if (NUM_AUXEVOL_GFS > 0) {
      PUParray(p, gridfuncs.auxevol_gfs, NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    }
  }
}

// PUP routine for struct charecomm_struct
void pup_charecomm_struct(PUP::er &p, charecomm_struct &cc, const params_struct &params, const params_struct &params_chare) {
  const int ntot = params.Nxx_plus_2NGHOSTS0 * params.Nxx_plus_2NGHOSTS1 * params.Nxx_plus_2NGHOSTS2;
  const int ntotchare = params_chare.Nxx_plus_2NGHOSTS0 * params_chare.Nxx_plus_2NGHOSTS1 * params_chare.Nxx_plus_2NGHOSTS2;
  if (p.isUnpacking()) {
    cc.globalidx3pt_to_chareidx3 = (int *restrict)malloc(sizeof(int) * ntot);
    cc.globalidx3pt_to_localidx3pt = (int *restrict)malloc(sizeof(int) * ntot);
    cc.localidx3pt_to_globalidx3pt = (int *restrict)malloc(sizeof(int) * ntotchare);
  }
  PUParray(p, cc.globalidx3pt_to_chareidx3, ntot);
  PUParray(p, cc.globalidx3pt_to_localidx3pt, ntot);
  PUParray(p, cc.localidx3pt_to_globalidx3pt, ntotchare);
}
// PUP routine for struct diagnostic_struct
void pup_diagnostic_struct(PUP::er &p, diagnostic_struct &ds, const params_struct &params_chare) {
  p | ds.num_output_quantities;
  p | ds.tot_num_diagnostic_1d_y_pts;
  p | ds.tot_num_diagnostic_1d_z_pts;
  p | ds.tot_num_diagnostic_2d_xy_pts;
  p | ds.tot_num_diagnostic_2d_yz_pts;
  p | ds.num_diagnostic_1d_y_pts;
  p | ds.num_diagnostic_1d_z_pts;
  p | ds.num_diagnostic_2d_xy_pts;
  p | ds.num_diagnostic_2d_yz_pts;

  if (p.isUnpacking()) {
    ds.localidx3_diagnostic_1d_y_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_y_pts);
    ds.locali0_diagnostic_1d_y_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_y_pts);
    ds.locali1_diagnostic_1d_y_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_y_pts);
    ds.locali2_diagnostic_1d_y_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_y_pts);
    ds.offset_diagnostic_1d_y_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_y_pts);

    ds.localidx3_diagnostic_1d_z_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_z_pts);
    ds.locali0_diagnostic_1d_z_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_z_pts);
    ds.locali1_diagnostic_1d_z_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_z_pts);
    ds.locali2_diagnostic_1d_z_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_z_pts);
    ds.offset_diagnostic_1d_z_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_z_pts);

    ds.localidx3_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_xy_pts);
    ds.locali0_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_xy_pts);
    ds.locali1_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_xy_pts);
    ds.locali2_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_xy_pts);
    ds.offset_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_xy_pts);

    ds.localidx3_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_yz_pts);
    ds.locali0_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_yz_pts);
    ds.locali1_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_yz_pts);
    ds.locali2_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_yz_pts);
    ds.offset_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_yz_pts);
  }

  PUParray(p, ds.localidx3_diagnostic_1d_y_pt, ds.num_diagnostic_1d_y_pts);
  PUParray(p, ds.locali0_diagnostic_1d_y_pt, ds.num_diagnostic_1d_y_pts);
  PUParray(p, ds.locali1_diagnostic_1d_y_pt, ds.num_diagnostic_1d_y_pts);
  PUParray(p, ds.locali2_diagnostic_1d_y_pt, ds.num_diagnostic_1d_y_pts);
  PUParray(p, ds.offset_diagnostic_1d_y_pt, ds.num_diagnostic_1d_y_pts);

  PUParray(p, ds.localidx3_diagnostic_1d_z_pt, ds.num_diagnostic_1d_z_pts);
  PUParray(p, ds.locali0_diagnostic_1d_z_pt, ds.num_diagnostic_1d_z_pts);
  PUParray(p, ds.locali1_diagnostic_1d_z_pt, ds.num_diagnostic_1d_z_pts);
  PUParray(p, ds.locali2_diagnostic_1d_z_pt, ds.num_diagnostic_1d_z_pts);
  PUParray(p, ds.offset_diagnostic_1d_z_pt, ds.num_diagnostic_1d_z_pts);

  PUParray(p, ds.localidx3_diagnostic_2d_xy_pt, ds.num_diagnostic_2d_xy_pts);
  PUParray(p, ds.locali0_diagnostic_2d_xy_pt, ds.num_diagnostic_2d_xy_pts);
  PUParray(p, ds.locali1_diagnostic_2d_xy_pt, ds.num_diagnostic_2d_xy_pts);
  PUParray(p, ds.locali2_diagnostic_2d_xy_pt, ds.num_diagnostic_2d_xy_pts);
  PUParray(p, ds.offset_diagnostic_2d_xy_pt, ds.num_diagnostic_2d_xy_pts);

  PUParray(p, ds.localidx3_diagnostic_2d_yz_pt, ds.num_diagnostic_2d_yz_pts);
  PUParray(p, ds.locali0_diagnostic_2d_yz_pt, ds.num_diagnostic_2d_yz_pts);
  PUParray(p, ds.locali1_diagnostic_2d_yz_pt, ds.num_diagnostic_2d_yz_pts);
  PUParray(p, ds.locali2_diagnostic_2d_yz_pt, ds.num_diagnostic_2d_yz_pts);
  PUParray(p, ds.offset_diagnostic_2d_yz_pt, ds.num_diagnostic_2d_yz_pts);

  PUParray(p, ds.filename_1d_y, 256);
  PUParray(p, ds.filename_1d_z, 256);
  PUParray(p, ds.filename_2d_xy, 256);
  PUParray(p, ds.filename_2d_yz, 256);
}
// PUP routine for struct tmpBuffers_struct
void pup_tmpBuffers_struct(PUP::er &p, tmpBuffers_struct &tmpBuffers, const params_struct &params, const nonlocalinnerbc_struct &nonlocalinnerbc,
                           const MoL_gridfunctions_struct &gridfuncs) {
  const int Nxx_plus_2NGHOSTS_face0 = params.Nxx_plus_2NGHOSTS1 * params.Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS_face1 = params.Nxx_plus_2NGHOSTS0 * params.Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS_face2 = params.Nxx_plus_2NGHOSTS0 * params.Nxx_plus_2NGHOSTS1;
  const int max_sync_gfs = gridfuncs.max_sync_gfs;
  size_t size_EW = static_cast<size_t>(max_sync_gfs) * NGHOSTS * Nxx_plus_2NGHOSTS_face0;
  size_t size_NS = static_cast<size_t>(max_sync_gfs) * NGHOSTS * Nxx_plus_2NGHOSTS_face1;
  size_t size_TB = static_cast<size_t>(max_sync_gfs) * NGHOSTS * Nxx_plus_2NGHOSTS_face2;
  if (p.isUnpacking()) {
    tmpBuffers.tmpBuffer_EW = (REAL *restrict)malloc(sizeof(REAL) * size_EW);
    tmpBuffers.tmpBuffer_NS = (REAL *restrict)malloc(sizeof(REAL) * size_NS);
    tmpBuffers.tmpBuffer_TB = (REAL *restrict)malloc(sizeof(REAL) * size_TB);
  }
  const int tot_num_dst_chares = nonlocalinnerbc.tot_num_dst_chares;
  const int tot_num_src_chares = nonlocalinnerbc.tot_num_src_chares;
  const int *num_srcpts_tosend_each_chare = nonlocalinnerbc.num_srcpts_tosend_each_chare;
  const int *num_srcpts_each_chare = nonlocalinnerbc.num_srcpts_each_chare;
  if (p.isUnpacking()) {
    tmpBuffers.tmpBuffer_innerbc_send = (REAL **)malloc(tot_num_dst_chares * sizeof(REAL *));
    for (int which_chare = 0; which_chare < tot_num_dst_chares; which_chare++) {
      tmpBuffers.tmpBuffer_innerbc_send[which_chare] =
          (REAL *restrict)malloc(sizeof(REAL) * max_sync_gfs * num_srcpts_tosend_each_chare[which_chare]);
    }
    tmpBuffers.tmpBuffer_innerbc_receiv = (REAL **)malloc(tot_num_src_chares * sizeof(REAL *));
    for (int which_chare = 0; which_chare < tot_num_src_chares; which_chare++) {
      tmpBuffers.tmpBuffer_innerbc_receiv[which_chare] = (REAL *restrict)malloc(sizeof(REAL) * max_sync_gfs * num_srcpts_each_chare[which_chare]);
    }
  }
}
// PUP routine for struct nonlocalinnerbc_struct
void pup_nonlocalinnerbc_struct(PUP::er &p, nonlocalinnerbc_struct &nonlocal, const commondata_struct &commondata) {
  const int Nchare0 = commondata.Nchare0;
  const int Nchare1 = commondata.Nchare1;
  const int Nchare2 = commondata.Nchare2;
  const int tot_num_chares = Nchare0 * Nchare1 * Nchare2;

  p | nonlocal.tot_num_src_chares;
  p | nonlocal.tot_num_dst_chares;
  if (p.isUnpacking()) {
    nonlocal.idx3_of_src_chares = (int *restrict)malloc(sizeof(int) * nonlocal.tot_num_src_chares);
    nonlocal.idx3chare_to_src_chare_id = (int *restrict)malloc(sizeof(int) * tot_num_chares);
    nonlocal.num_srcpts_each_chare = (int *restrict)malloc(sizeof(int) * nonlocal.tot_num_src_chares);

    nonlocal.idx3_of_dst_chares = (int *restrict)malloc(sizeof(int) * nonlocal.tot_num_dst_chares);
    nonlocal.idx3chare_to_dst_chare_id = (int *restrict)malloc(sizeof(int) * nonlocal.tot_num_dst_chares);
    nonlocal.num_srcpts_tosend_each_chare = (int *restrict)malloc(sizeof(int) * nonlocal.tot_num_dst_chares);

    nonlocal.map_srcchare_and_srcpt_id_to_linear_id = (int **)malloc(sizeof(int *) * nonlocal.tot_num_src_chares);
    nonlocal.globalidx3_srcpts = (int **)malloc(sizeof(int *) * nonlocal.tot_num_src_chares);

    nonlocal.globalidx3_srcpts_tosend = (int **)malloc(sizeof(int *) * nonlocal.tot_num_dst_chares);
  }
  PUParray(p, nonlocal.idx3_of_src_chares, nonlocal.tot_num_src_chares);
  PUParray(p, nonlocal.idx3chare_to_src_chare_id, tot_num_chares);
  PUParray(p, nonlocal.num_srcpts_each_chare, nonlocal.tot_num_src_chares);

  PUParray(p, nonlocal.idx3_of_dst_chares, nonlocal.tot_num_dst_chares);
  PUParray(p, nonlocal.idx3chare_to_dst_chare_id, nonlocal.tot_num_dst_chares);
  PUParray(p, nonlocal.num_srcpts_tosend_each_chare, nonlocal.tot_num_dst_chares);

  for (int src_chare = 0; src_chare < nonlocal.tot_num_src_chares; src_chare++) {
    if (p.isUnpacking()) {
      nonlocal.map_srcchare_and_srcpt_id_to_linear_id[src_chare] = (int *restrict)malloc(sizeof(int) * nonlocal.num_srcpts_each_chare[src_chare]);
      nonlocal.globalidx3_srcpts[src_chare] = (int *restrict)malloc(sizeof(int) * nonlocal.num_srcpts_each_chare[src_chare]);
    }
    PUParray(p, nonlocal.map_srcchare_and_srcpt_id_to_linear_id[src_chare], nonlocal.num_srcpts_each_chare[src_chare]);
    PUParray(p, nonlocal.globalidx3_srcpts[src_chare], nonlocal.num_srcpts_each_chare[src_chare]);
  }
  for (int dst_chare = 0; dst_chare < nonlocal.tot_num_dst_chares; dst_chare++) {
    if (p.isUnpacking()) {
      nonlocal.globalidx3_srcpts_tosend[dst_chare] = (int *restrict)malloc(sizeof(int) * nonlocal.num_srcpts_tosend_each_chare[dst_chare]);
    }
    PUParray(p, nonlocal.globalidx3_srcpts_tosend[dst_chare], nonlocal.num_srcpts_tosend_each_chare[dst_chare]);
  }
}
// PUP routine for struct griddata
// During time evolution, need params from griddata which is used to unpack charecomm_struct in griddata_chare and xx for diagnostics.
void pup_griddata(PUP::er &p, griddata_struct &gd) {
  pup_params_struct(p, gd.params);
  if (p.isUnpacking()) {
    gd.xx[0] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS0);
    gd.xx[1] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS1);
    gd.xx[2] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS2);
  }
  PUParray(p, gd.xx[0], gd.params.Nxx_plus_2NGHOSTS0);
  PUParray(p, gd.xx[1], gd.params.Nxx_plus_2NGHOSTS1);
  PUParray(p, gd.xx[2], gd.params.Nxx_plus_2NGHOSTS2);
}
// PUP routine for struct griddata_chare
// For unpacking order is important; unpacked structs are used for unpacking the subsequent structs.
void pup_griddata_chare(PUP::er &p, griddata_struct &gd, const params_struct &params, const commondata_struct &commondata) {

  pup_params_struct(p, gd.params);

  if (p.isUnpacking()) {
    gd.xx[0] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS0);
    gd.xx[1] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS1);
    gd.xx[2] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS2);
  }
  PUParray(p, gd.xx[0], gd.params.Nxx_plus_2NGHOSTS0);
  PUParray(p, gd.xx[1], gd.params.Nxx_plus_2NGHOSTS1);
  PUParray(p, gd.xx[2], gd.params.Nxx_plus_2NGHOSTS2);

  pup_diagnostic_struct(p, gd.diagnosticstruct, gd.params);

  pup_charecomm_struct(p, gd.charecommstruct, params, gd.params);

  pup_bc_struct(p, gd.bcstruct);

  pup_nonlocalinnerbc_struct(p, gd.nonlocalinnerbcstruct, commondata);

  pup_MoL_gridfunctions_struct(p, gd.gridfuncs, gd.params, commondata);

  pup_tmpBuffers_struct(p, gd.tmpBuffers, gd.params, gd.nonlocalinnerbcstruct, gd.gridfuncs);

  if (p.isUnpacking()) {
    gd.rfmstruct = (rfm_struct *)malloc(sizeof(rfm_struct));
    rfm_precompute_malloc(&commondata, &gd.params, gd.rfmstruct);
    rfm_precompute_defines(&commondata, &gd.params, gd.rfmstruct, gd.xx);
  }
}

/**
 * This file implements a collection of Pack-Unpack (PUP) routines used in Charm++ for checkpointing,
 * and load balancing in the superB framework.
 * It includes routines for serializing and deserializing:
 *     - commondata_struct and params_struct,
 *     - rfm_struct with reference metric precomputation and memory allocation,
 *     - inner and outer boundary condition structures (innerpt_bc_struct, outerpt_bc_struct, bc_info_struct, bc_struct),
 *     - MoL grid functions (MoL_gridfunctions_struct),
 *     - chare communication structures (charecomm_struct),
 *     - diagnostic information (diagnostic_struct),
 *     - temporary buffers and nonlocal inner boundary conditions (tmpBuffers_struct, nonlocalinnerbc_struct),
 *     - and grid data structures (griddata_struct, griddata_chare).
 * This comprehensive set of routines is crucial for efficient data management and communication in high-performance, parallel simulations.
 */
void superB_pup_routines() {
  // This space intentionally left blank:
  // NRPy requires C files share the same name as the primary function within the C file.
} // END FUNCTION superB_pup_routines
