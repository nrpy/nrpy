void MoL_free_intermediate_stage_gfs(MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_malloc_intermediate_stage_gfs(const commondata_struct *restrict commondata, const params_struct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void MoL_step_forward_in_time(commondata_struct *restrict commondata, griddata_struct *restrict griddata, const REAL time_start, const int which_RK_substep, const int which_MOL_part);
void MoL_sync_data_defines(MoL_gridfunctions_struct *restrict gridfuncs);
void apply_bcs_inner_only(const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs);
void apply_bcs_inner_only_nonlocal(const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, const nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, const int NUM_GFS, REAL *restrict gfs, const int *gfs_to_sync, const int8_t* gf_parity_types, REAL **restrict tmpBuffer_innerbc_receiv);
void apply_bcs_inner_only_specific_auxgfs(const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs, const int num_gfs, const int *gfs_to_sync);
void apply_bcs_outerextrap_and_inner(const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs);
void apply_bcs_outerradiation_and_inner(const commondata_struct *restrict commondata, const params_struct *restrict params,
    const bc_struct *restrict bcstruct, REAL *restrict xx[3],
    const REAL custom_wavespeed[NUM_EVOL_GFS],
    const REAL custom_f_infinity[NUM_EVOL_GFS],
    REAL *restrict gfs, REAL *restrict rhs_gfs);
void apply_bcs_outerradiation_and_inner__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
    const bc_struct *restrict bcstruct, REAL *restrict xx[3],
    const REAL custom_wavespeed[NUM_EVOL_GFS],
    const REAL custom_f_infinity[NUM_EVOL_GFS],
    REAL *restrict gfs, REAL *restrict rhs_gfs);
void auxevol_gfs_set_to_constant(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs);
void auxevol_gfs_set_to_constant__rfm__SinhSpherical(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs);
void bcstruct_chare_set_up(const commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, const charecomm_struct *restrict charecommstruct, REAL *restrict xx[3], const bc_struct *restrict bcstruct, bc_struct *restrict bcstruct_chare, nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, const int chare_index[3]);
void bcstruct_chare_set_up__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, const charecomm_struct *restrict charecommstruct, REAL *restrict xx[3], const bc_struct *restrict bcstruct, bc_struct *restrict bcstruct_chare, nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, const int chare_index[3]);
void bcstruct_set_up(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct);
void bcstruct_set_up__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct);
void cfl_limited_timestep(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]);
void charecommstruct_set_up(const commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, charecomm_struct* charecommstruct, const int thischareindex[3]);
void charecommstruct_set_up__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, charecomm_struct* charecommstruct, const int thischareindex[3]);
void cmdline_input_and_parfile_parser(commondata_struct *restrict commondata, int argc, const char *argv[]);
void commondata_struct_set_to_default(commondata_struct *restrict commondata);
void diagnostic_gfs_set(const commondata_struct *restrict commondata, const griddata_struct *restrict griddata, REAL *restrict diagnostic_gfs[MAXNUMGRIDS]);
void diagnostics(commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_chare,
                 const REAL *restrict gridfuncs_diags[MAXNUMGRIDS],
                 const int chare_index[3], const int grid, Ck::IO::Session token, const int which_diagnostics_part);
void diagnostics_nearest(commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_chare,
                         const REAL *restrict gridfuncs_diags[MAXNUMGRIDS], const int chare_index[3], Ck::IO::Session token, const int which_diagnostics_part);
void diagnostics_nearest_1d_y_and_z_axes(commondata_struct *restrict commondata, const int grid,
                                         const params_struct *restrict params, const params_struct *restrict params_chare,const REAL *restrict xx[3], const REAL *restrict xx_chare[3],
                                         const int NUM_GFS_NEAREST, const int which_gfs[],const char **diagnostic_gf_names,
                                         const REAL *restrict gridfuncs_diags[],
                                         const charecomm_struct *restrict charecommstruct, diagnostic_struct *restrict diagnosticstruct,
                                         const int chare_index[3], Ck::IO::Session token, const int which_diagnostics_part);
void diagnostics_nearest_1d_y_and_z_axes__rfm__SinhSpherical(commondata_struct *restrict commondata, const int grid,
                                                             const params_struct *restrict params, const params_struct *restrict params_chare, const REAL *restrict xx[3], const REAL *restrict xx_chare[3],
                                                             const int NUM_GFS_NEAREST, const int which_gfs[], const char **diagnostic_gf_names,
                                                             const REAL *restrict gridfuncs_diags[],
                                                             const charecomm_struct *restrict charecommstruct, diagnostic_struct *restrict diagnosticstruct,
                                                             const int chare_index[3], Ck::IO::Session token, const int which_diagnostics_part);
void diagnostics_nearest_2d_xy_and_yz_planes(commondata_struct *restrict commondata, const int grid,
                                         const params_struct *restrict params, const params_struct *restrict params_chare,const REAL *restrict xx[3], const REAL *restrict xx_chare[3],
                                         const int NUM_GFS_NEAREST, const int which_gfs[],const char **diagnostic_gf_names,
                                         const REAL *restrict gridfuncs_diags[],
                                         const charecomm_struct *restrict charecommstruct, diagnostic_struct *restrict diagnosticstruct,
                                         const int chare_index[3], Ck::IO::Session token, const int which_diagnostics_part);
void diagnostics_nearest_2d_xy_and_yz_planes__rfm__SinhSpherical(commondata_struct *restrict commondata, const int grid,
                                                                 const params_struct *restrict params, const params_struct *restrict params_chare, const REAL *restrict xx[3], const REAL *restrict xx_chare[3],
                                                                 const int NUM_GFS_NEAREST, const int which_gfs[], const char **diagnostic_gf_names,
                                                                 const REAL *restrict gridfuncs_diags[],
                                                                 const charecomm_struct *restrict charecommstruct, diagnostic_struct *restrict diagnosticstruct,
                                                                 const int chare_index[3], Ck::IO::Session token, const int which_diagnostics_part);
void diagnostics_nearest_grid_center(commondata_struct *restrict commondata, const int grid, const params_struct *restrict params, const params_struct *restrict params_chare,
                                     const REAL *restrict xx[3], const int NUM_GFS_NEAREST, const int which_gfs[], const char **diagnostic_gf_names,
                                     const REAL *restrict gridfuncs_diags[], const int chare_index[3]);
void diagnostics_nearest_grid_center__rfm__SinhSpherical(commondata_struct *restrict commondata, const int grid, const params_struct *restrict params,
                                                        const params_struct *restrict params_chare,
                                                         const REAL *restrict xx[3], const int NUM_GFS_NEAREST, const int which_gfs[],
                                                         const char **diagnostic_gf_names, const REAL *restrict gridfuncs_diags[],
                                                         const int chare_index[3]);
void diagnostics_volume_integration(commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_chare,
                                    const REAL *restrict gridfuncs_diags[MAXNUMGRIDS],
                                    const int chare_index[3],
                                    const int which_diagnostics_part);
void ds_min_single_pt(const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict ds_min);
void ds_min_single_pt__rfm__SinhSpherical(const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict ds_min);
void griddata_free(commondata_struct *restrict commondata, griddata_struct *restrict griddata, const bool free_non_y_n_gfs_and_core_griddata_pointers);
void initial_data(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
void initial_guess_single_point(const REAL xx0, const REAL xx1, const REAL xx2,  REAL *restrict uu_ID, REAL *restrict vv_ID);
void initialize_yn_and_non_yn_gfs_to_nan(const commondata_struct *restrict commondata, const params_struct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs);
void log10_L2norm_gf(commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                const REAL integration_radius, const int gf_index, const REAL *restrict in_gf, REAL localsums_for_residualH[2]);
void numerical_grid_params_Nxx_dxx_xx(const commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], const int Nx[3], const bool apply_convergence_factor_and_set_xxminmax_defaults);
void numerical_grid_params_Nxx_dxx_xx__rfm__SinhSpherical(const commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3], const int Nx[3], const bool apply_convergence_factor_and_set_xxminmax_defaults);
void numerical_grid_params_Nxx_dxx_xx_chare(commondata_struct *restrict commondata, const params_struct *restrict params, params_struct *restrict params_chare, REAL *restrict xx[3], const int chare_index[3]);
void numerical_grid_params_Nxx_dxx_xx_chare__rfm__SinhSpherical(commondata_struct *restrict commondata, const params_struct *restrict params, params_struct *restrict params_chare, REAL *restrict xx[3], const int chare_index[3]);
void numerical_grids_and_timestep(commondata_struct *restrict commondata, griddata_struct *restrict griddata, bool calling_for_first_time);
void numerical_grids_chare(commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_chare, const int chare_index[3]);
void params_struct_set_to_default(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
void progress_indicator(commondata_struct *restrict commondata, const griddata_struct *restrict griddata);
void residual_H_compute_all_points(const commondata_struct *restrict commondata, const params_struct *restrict params,
                REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                REAL *restrict dest_gf_address);
void residual_H_compute_all_points__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                REAL *restrict dest_gf_address);
void rfm_precompute_defines(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct, REAL *restrict xx[3]);
void rfm_precompute_defines__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct, REAL *restrict xx[3]);
void rfm_precompute_free(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct);
void rfm_precompute_free__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct);
void rfm_precompute_malloc(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct);
void rfm_precompute_malloc__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct);
void rhs_eval(const commondata_struct *restrict commondata, const params_struct *restrict params, const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs);
void rhs_eval__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs);
void sqrt_detgammahat_d3xx_volume_element(const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict dV);
void sqrt_detgammahat_d3xx_volume_element__rfm__SinhSpherical(const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict dV);
void stop_conditions_check(commondata_struct *restrict commondata);
void superB_pup_routines( );
void timestepping_free_memory_tmpBuffer(const nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, tmpBuffers_struct *restrict tmpBuffers);
void timestepping_malloc_tmpBuffer(const commondata_struct *restrict commondata, const params_struct *restrict params, const MoL_gridfunctions_struct *restrict gridfuncs, const nonlocalinnerbc_struct *restrict nonlocalinnerbcstruct, tmpBuffers_struct *restrict tmpBuffers);
void xx_to_Cart(const params_struct *restrict params, const REAL xx[3], REAL xCart[3]);
void xx_to_Cart__rfm__SinhSpherical(const params_struct *restrict params, const REAL xx[3], REAL xCart[3]);
