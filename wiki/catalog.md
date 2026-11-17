# Global Catalog

> Keyword and scope index for normal leaves across `wiki/`. · Status: confirmed

Search this table with `rg -n -i -F '<term>' wiki/catalog.md` when router titles do not resolve a topic; never read the full table.

| Page | Scope | Query terms |
| --- | --- | --- |
| [Build And Run](architecture/build-and-run.md) | Install, direct-run setup, backend toolchains, and nrpyinline execution. | PYTHONPATH, requirements-dev.txt, nvcc, --cuda, Makefile |
| [C And Embedded C Style](architecture/c-and-embedded-c-style.md) | Handwritten C/CUDA/H formatting, Doxygen, and END markers. | clang-format, snake_case, UPPER_SNAKE_CASE, restrict |
| [Contribution Style And Static Analysis](architecture/contribution-style-and-static-analysis.md) | Commit attribution, coordinate admission, and guard rules. | coding_style.md, Co-authored-by:, single_file_static_analysis.sh |
| [Generated Output Boundaries](architecture/generated-output-boundaries.md) | Boundary between handwritten NRPy code and generated files. | project/, Makefile, .par, dendro_bssn.py |
| [Overview](architecture/overview.md) | Source tree layout, package metadata, and SymPy codegen stages. | c_codegen, py_codegen, CFunction_dict, PyFunction_dict |
| [Python Coding Style](architecture/python-coding-style.md) | Python formatting, naming, imports, docstrings, and typing. | Black, isort, Sphinx, __main__ |
| [C Codegen](core/c-codegen.md) | SymPy-to-C emission, CSE options, FD stencils, and SIMD C. | c_codegen, CCodeGen, GoldenKernels, enable_fd_codegen, enable_simd |
| [C Function Registry](core/c-function-registry.md) | CFunction objects, CFunction_dict, and coordinate wrappers. | register_CFunction, generate_full_function, CoordSystem_for_wrapper_func |
| [Finite Difference](core/finite-difference.md) | FD coefficients, derivative extraction, reads, and KO stencils. | fd_order, FDFunction, compute_fdcoeffs_fdstencl, dD, dDD, dKOD |
| [Gridfunctions And Parameters](core/gridfunctions-and-parameters.md) | GridFunction, NRPyParameter, CodeParameter, and runtime structs. | register_gridfunctions, glb_params_dict, griddata_struct, commondata_struct |
| [CSE And Printer Support](core/helpers/cse-and-printer-support.md) | CSE preprocessing, deterministic ordering, and C/JAX printers. | cse_preprocess, cse_postprocess, custom_functions_for_SymPy_ccode, NRPyJaxPrinter |
| [Loop Kernel And Device Helpers](core/helpers/loop-kernel-and-device-helpers.md) | Loop emission, GPUKernel wrappers, and CUDA launch helpers. | simple_loop, loop1D, tile_size, streamid, __global__ |
| [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md) | clang-format, validate_strings, .nrpycache, and file updates. | clang_format, cached_simplify, conditional_file_updater |
| [Parallel Codegen Orchestration](core/helpers/parallel-codegen-orchestration.md) | Multiprocessing code generation and global registry merges. | parallel_codegen, pcg_registration_phase, register_CFunctions |
| [SIMD And Intrinsic Support](core/helpers/simd-and-intrinsic-support.md) | Symbolic SIMD rewrites and CPU/CUDA intrinsic macro headers. | expr_convert_to_simd_intrins, simd_intrinsics.h, cuda_intrinsics.h |
| [Symbolic Expression Utilities](core/helpers/symbolic-expression-utilities.md) | Free-symbol extraction, float-to-rational f2r, and products. | expression_utils, float_to_rational, ExprTree, expr_tree.py |
| [Indexed Expressions](core/indexed-expressions.md) | Symbolic rank-1..4 tensor declarations, symmetries, and inverses. | indexedexp, declarerank1, declarerank2, symmetry_axes, inverter |
| [Python Codegen](core/python-codegen.md) | SymPy-to-JAX Python assignment emission and PyCodeGen. | py_codegen, jnp, postproc_substitution_dict |
| [Python Function Registry](core/python-function-registry.md) | PyFunction text assembly and PyFunction_dict for JAX modules. | register_PyFunction, full_function, subdirectory |
| [Reference Metrics](core/reference-metrics.md) | Coordinate reference metrics, rfm_dict, and rfm_precompute. | ReferenceMetric, CoordSystem, enable_rfm_precompute, ReU |
| [Conformally Flat Elliptic](equations/conformally-flat-elliptic.md) | NRPyElliptic conformally flat relaxation RHSs and sources. | ConformallyFlat_RHSs, ConformallyFlat_SourceTerms, psi_background |
| [Equation Setup Style](equations/equation-setup-style.md) | SymPy equation conventions, derivative names, and trusted dicts. | cached_simplify, dD, dDD, dupD, trusted_dict |
| [BSSN Family](equations/general-relativity/bssn-family.md) | Reference-metric BSSN quantities, RHSs, gauge, constraints, KO. | BSSN_quantities, BSSN_RHSs, BSSN_gauge_RHSs, BSSN_constraints, nrpylatex |
| [fCCZ4 Formulation](equations/general-relativity/fccz4-formulation.md) | Reference-metric fCCZ4 variables, Z4 damping, and gauge. | FCCZ4RHSs, FCCZ4Constraints, kappa1, kappa2, kappa3, Mewes |
| [fCCZ4 Implementation And Validation](equations/general-relativity/fccz4-implementation-and-validation.md) | fCCZ4 modules, YBS-MOM, assembly, and tests. | fCCZ4_RHSs, fCCZ4_constraints, build_fccz4_expression_bundle, H_Z4 |
| [Fishbone-Moncrief](equations/general-relativity/fishbone-moncrief.md) | Equilibrium relativistic torus initial data and B-field seed. | FishboneMoncriefID, Boyer-Lindquist, BmagU, hm1 |
| [Geodesics](equations/general-relativity/geodesics.md) | Analytic Kerr/Schwarzschild metrics, geodesic RHSs, invariants. | AnalyticSpacetimes, GeodesicEquations, geodesic_eom_rhs_photon, geodesic_eom_rhs_massive |
| [Horizon Diagnostics](equations/general-relativity/horizon-diagnostics.md) | Apparent-horizon expansion Theta, area, centroid, and spin. | ExpansionFunctionThetaClass, SpECTRESpinEstimateClass, spin_NewtonRaphson |
| [Initial Data](equations/general-relativity/initial-data.md) | Analytic Cartesian and spherical ADM initial-data builders. | InitialData_Cartesian, InitialData_Spherical, BrillLindquist, UIUCBlackHole |
| [Metric Conversions And Matter](equations/general-relativity/metric-conversions-and-matter.md) | ADM/BSSN/4-metric conversions, T4munu projections, and boosts. | ADM_to_BSSN, BSSN_to_ADM, g4DD, LorentzBoost |
| [Psi4 And Tetrads](equations/general-relativity/psi4-and-tetrads.md) | Weyl scalar Psi4 evaluation and BCL or symbolic tetrads. | Psi4Tetrads, psi4_re, psi4_im, leave_symbolic |
| [YBS Gamma-Constraint Adjustment](equations/general-relativity/ybs-gamma-constraint-adjustment.md) | Opt-in YBS_chi Gamma-constraint term in BSSN LambdaBar RHS. | enable_YBS_Gamma_constraint_adjustment, xi, BSSN_RHSs |
| [YBS-MOM Timestep-Scaled Momentum Adjustment](equations/general-relativity/ybs-momentum-damping.md) | Yo-Lin-Cao momentum-gradient adjustment with CAHD scaling. | C_YBS_mom, enable_YBS_momentum_constraint_adjustment, DSMINGF, dsmin |
| [Geometry And Special-Function Support](equations/geometry-and-special-function-support.md) | Basis transforms, fisheye, SO(3)/quaternion rotations, _s Y_{lm}. | BasisTransforms, GeneralRFMFisheye, SO3_rotations, spherical harmonics |
| [GRHD](equations/grhd.md) | Reference-metric Valencia GRHD conservatives, fluxes, and speeds. | GRHD_Equations, T4UU, u4U, HLL_fluxes |
| [GRMHD](equations/grmhd.md) | Magnetic stress-energy in GRHD, magnetosonic bounds, HLL fluxes. | GRMHD_Equations, BmagU, smallb4U, Alfven, cmin, cmax |
| [Aligned-Spin Calibration And Remnant](equations/seobnr/aligned-spin-calibration-and-remnant.md) | SEOBNRv5 calibration constants, ISCO helpers, and remnant fits. | SEOBNRv5_aligned_spin_constants, HBR2016, UIB2016 |
| [Aligned-Spin Hamiltonian](equations/seobnr/aligned-spin-hamiltonian.md) | Conservative SEOBNRv5 aligned-spin Hamiltonian and derivatives. | SEOBNRv5_aligned_spin_Hamiltonian_quantities, Hreal, dHreal_dr, dHreal_dpphi |
| [Aligned-Spin Waveforms](equations/seobnr/aligned-spin-waveforms.md) | SEOBNRv5 factorized modes, radiation flux, strain, and NQC. | SEOBNRv5_aligned_spin_waveform_quantities, hlm, Omega_circ, vphi |
| [BOB Waveforms](equations/seobnr/bob-waveforms.md) | Backwards-One-Body merger-ringdown strain, NQC, and BOBv2 modes. | BOB_aligned_spin_waveform_quantities, BOB_v2_waveform_quantities, QNM |
| [SEOBNR Precessing Dynamics](equations/seobnr/precessing-dynamics.md) | Quasi-precessing SEOBNRv5 Hamiltonian and spin/orbital RHSs. | SEOBNRv5_quasi_precessing_spin_Hamiltonian, spin_evolution, ln_dot |
| [SEOBNR Precessing Rotations And Ringdown](equations/seobnr/precessing-rotations-and-ringdown.md) | Co-precessing Euler/Wigner rotations, inertial polarizations, omega_prec. | SEOBNRv5_Coprecessing_Rotations, wigner_d_small, J_f |
| [TOV Equations](equations/tov-equations.md) | Tolman-Oppenheimer-Volkoff stellar-equilibrium ODE RHSs. | TOV_Equations, r_Schw, rho_energy, dP_dr, dM_dr |
| [Wave Equation](equations/wave-equation.md) | Cartesian/curvilinear scalar wave RHSs and analytic solutions. | WaveEquation_RHSs, WaveEquationCurvilinear_RHSs, initial data |
| [Apparent Horizon Library](examples/apparent-horizon-library.md) | Standalone libBHaHAHA.a generator and horizon-enabled examples. | bhahaha.py, BHaHAHA_header.h, find_horizons, bah_find_horizon |
| [Einstein Toolkit Thorn Generators](examples/einstein-toolkit-thorn-generators.md) | ETLegacy and CarpetX Baikal/WaveToy thorn generators. | carpet_baikal_thorns, carpetx_baikal_thorns, carpet_wavetoy_thorns |
| [Elliptic Initial Data](examples/elliptic-initial-data.md) | Standalone BHaH and superB conformally flat NRPyElliptic codes. | nrpyelliptic_conformally_flat, superB_nrpyelliptic_conformally_flat |
| [Example Generator Catalog](examples/example-generator-catalog.md) | Table of runnable nrpy/examples/*.py scripts and dependencies. | project/, main.yml, sebob_consistency_check.py |
| [Geodesic Raytracing](examples/geodesic-raytracing.md) | Standalone massive/photon geodesic integrators and raytracer. | mass_geodesic_integrator, photon_geodesic_integrator, photon_geodesic_batch_integrator |
| [Matter TOV Workflows](examples/matter-tov-workflows.md) | TOVola star, Cowling hydro, GRoovy TOV, and MANGA C library. | tovola_neutron_star, hydro_without_hydro, groovy_TOV_BSSN, manga_bhah_lib |
| [Standalone GR/BHaH](examples/standalone-gr-bhah.md) | Single-patch BHaH GR generators, initial data, diagnostics, and build modes. | blackhole_spectroscopy, two_blackholes_collide, spinning_blackhole, fisheye |
| [superB Charm++ Workflows](examples/superb-charm-workflows.md) | Charm++ superB BBH, spectroscopy, and NRPyElliptic examples. | superB_two_blackholes_collide, superB_blackhole_spectroscopy, charmrun |
| [Wave Equation Generators](examples/wave-equation-generators.md) | First-run Cartesian wave walkthrough and curvilinear wave codes. | wave_equation_cartesian, wave_equation_curvilinear, wave_equation_multicoordinates |
| [Waveform JAX PN Generators](examples/waveform-jax-pn-generators.md) | SEOBNRv5, SEBOB, sebobv1_jax, and Post-Newtonian generators. | seobnrv5_aligned_spin_inspiral, sebobv2, nrpypn_quasicircular_momenta |
| [Generated Backend Comparison](generated-backend-comparison.md) | Side-by-side comparison of BHaH, ET, superB, JAX, and Dendro. | ETLegacy, CarpetX, CFunction_dict, BHaH_defines.h |
| [BHaHAHA Grids And Interpolation](infrastructures/bhah/bhahaha/bhahaha-grid-interpolation-and-boundaries.md) | Spherical grids, Lagrange interpolation, and parity BCs. | bah_numgrid, bah_interpolation, bah_bcstruct_set_up |
| [BHaHAHA Public C Header And Inputs](infrastructures/bhah/bhahaha/bhahaha-public-api-and-input-contract.md) | BHaHAHA_header.h structs, bah_find_horizon, and Nphi parity. | bhahaha_find_horizons, bhahaha_params_and_data_struct |
| [BHaHAHA Relaxation And State](infrastructures/bhah/bhahaha/bhahaha-relaxation-diagnostics-and-state.md) | Multiresolution KO relaxation, BBH mode, and .gp files. | bah_MoL_step_forward_in_time, bah_diagnostics, BBH_mode |
| [Curvilinear Boundaries And Stencils](infrastructures/bhah/curvilinear-boundaries-and-stencils.md) | CurviBoundaryConditions parity/Sommerfeld BCs and Nphi parity. | bcstruct_set_up, apply_bcs_inner_only, SymTP, SinhSymTP |
| [Diagnostics And Checkpointing](infrastructures/bhah/diagnostics-output-and-checkpointing.md) | BHaH 0D/1D/2D/volume diagnostics, raytracing, and checkpoints. | register_all_diagnostics, read_checkpoint, write_checkpoint |
| [Geodesics And Raytracing Runtime](infrastructures/bhah/geodesics-and-raytracing-runtime.md) | BHaH RKF45 geodesic integration and azimuthal metric stencils. | output_raytracing_data, ode_gsl_wrapper_massive, Lagrange |
| [GR Evolution And Constraints](infrastructures/bhah/gr-evolution-and-constraints.md) | BHaH BSSN/fCCZ4 RHS, Ricci, YBS/CAHD options, and constraints. | rhs_eval, Ricci_eval, enforce_detgbar_equals_detghat_trAzero, constraints_eval |
| [GR Initial Data And Psi4](infrastructures/bhah/gr-initial-data-and-psi4.md) | BHaH analytic/TwoPunctures/TOVola initial data and Psi4 modes. | ID_persist_struct, register_CFunction_psi4, psi4_spinweightm2_decomposition |
| [Grids And Reference Metrics](infrastructures/bhah/grids-coordinates-and-reference-metrics.md) | BHaH numerical grids, fisheye setup, Cart_to_xx, and CFL dt. | numerical_grids_and_timestep, rfm_precompute, CoordSystem_hash, rfm_struct |
| [GRoovy GRHD Runtime](infrastructures/bhah/groovy-grhd-runtime.md) | BHaH GRoovy GRHD gridfunctions, GRHayL EOS, PPM, and HLL. | GRHD_Equations, reconstruction_loop, calculate_HLL_flux_dirn_i, grhd_rhs_eval |
| [BHaH Project Assembly](infrastructures/bhah/lifecycle-and-project-assembly.md) | BHaH main.c registration, Makefile emission, and library mode. | register_CFunction_main_c, Makefile_helpers, create_lib |
| [MoL Time Integration](infrastructures/bhah/mol-time-integration.md) | BHaH Butcher tables, RK stage buffers, and MoL stepping. | MoLtimestepping, Butcher_dict, MoL_malloc_intermediate_stage_gfs, MoL_step_forward_in_time |
| [NRPyElliptic Workflow](infrastructures/bhah/nrpyelliptic-workflow.md) | BHaH conformally flat relaxation loop, residual_H, and stops. | nrpyelliptic_conformally_flat, residual_H_compute_all_points, stop_conditions_check |
| [Runtime Data, Parameters, Headers, And CLI](infrastructures/bhah/runtime-data-parameters-headers-and-cli.md) | BHaH_defines.h structs, parfile parser, and CUDA structs. | BHaH_defines_h, cmdline_input_and_parfiles, commondata_struct |
| [SEOBNR BOB Generated Library](infrastructures/bhah/seobnr-bob-generated-library.md) | BHaH C waveform libraries, GSL ODE integration, and resampling. | SEBOBv2_IMR_waveform, gsl_odeiv2_step_rk8pd, ringdown, dT |
| [SO(3) Rotation Helpers](infrastructures/bhah/so3-rotation-helpers.md) | BHaH generated SO(3) axis-angle C helpers and string tests. | so3_axis_angle_to_R, so3_build_R_from_hats, so3_apply_R_to_vector, validate_strings |
| [CarpetX Boundaries And RHS Init](infrastructures/carpetx/boundaries-and-rhs-initialization.md) | CarpetX NewRadX radiative boundaries and zero_rhss. | BoundaryConditions, schedule_ccl |
| [CarpetX Code Params And Loops](infrastructures/carpetx/code-parameters-includes-and-loops.md) | Cactus parameter reads, includes, and grid.loop_all_device. | cctk_Arguments, loop_region, DECLARE_CCTK_PARAMETERS, SIMD |
| [CarpetX GR ADM/BSSN And Matter](infrastructures/carpetx/gr-adm-bssn-and-matter-coupling.md) | CarpetX ADMBaseX import/export and TmunuBaseX coupling. | ADM_to_BSSN, BSSN_to_ADM, T4DD_to_T4UU |
| [CarpetX GR BSSN RHS And Constraints](infrastructures/carpetx/gr-bssn-rhs-ricci-constraints-and-validation.md) | CarpetX BSSN Ricci, RHS, constraints, and trusted dicts. | ODESolvers_RHS, ODESolvers_PostStep, CCTK_DELTA_SPACE |
| [CarpetX GR Poststep Repairs](infrastructures/carpetx/gr-poststep-repairs.md) | CarpetX lapse floor and determinant/trace-free Abar repair. | enforce_detgbar_equals_detghat_trAzero, floor_the_lapse |
| [CarpetX Thorn Assembly And CCL](infrastructures/carpetx/thorn-assembly-configuration-and-ccl-files.md) | CarpetX CCL files, configuration.ccl, and make.code.defn. | USES_INCLUDEs, evol_variables, ET_thorn_name |
| [BSSN Application Wiring](infrastructures/dendro/bssn-application-wiring.md) | Dendro BSSN W/chi evolution, SSL, CAHD, Brown Lambda, and q1. | Dendro_NRPy_BSSN, C_CAHD, YBS_chi, C_YBS_mom, lambdaU |
| [Constraints And Diagnostic Norms](infrastructures/dendro/constraints-and-diagnostic-norms.md) | Dendro covariant M_i, H_Z4/M/Lambda norms, and node RMS. | M_CONSTRAINT, LAMBDA_CONSTRAINT, Ctx::diagnostic_output |
| [fCCZ4 Application Wiring](infrastructures/dendro/fccz4-application-wiring.md) | Dendro fCCZ4 W/chi generation, CAHD, lambdaU init, and DIAG. | Dendro_NRPy_fCCZ4, initial_data_lambdaU, kappa1, kappa2 |
| [FD Profiles And Conformance](infrastructures/dendro/finite-difference-profiles-and-dendro-conformance.md) | Centered advection, KO order, padding, and BSSN_ELE_ORDER. | KO_FD_ORDER, REQUIRED_PADDING, NGHOSTS |
| [Octree Grid, AMR, And Time Stepping](infrastructures/dendro/grid-amr-and-time-stepping.md) | Dendro zip/unzip, ghosts, BH_WAMR remeshing, and RK4. | ot::Block, Ctx::rhs, physical_boundary_ghosts, BH_WAMR |
| [Gridfunctions, Naming, And Loops](infrastructures/dendro/gridfunctions-naming-and-loops.md) | DendroGridFunction ordering, pointer prefixes, and block loops. | cf_W_or_chi, state_h.dendro_state_name, in_, rhs_ |
| [Initial Data And Wave Diagnostics](infrastructures/dendro/initial-data-and-wave-diagnostics.md) | Dendro TwoPunctures setup, ADM integrals, and Psi4 modes. | --tpid, _ADM.dat, Ctx::adm_output, BSSN_GW_RADAII |
| [Production Validation And Checks](infrastructures/dendro/production-validation-and-deferred-checks.md) | Dendro MPI remesh/restart checks and dendro-validation CI. | dendro_application_check, dendrolib_capability_test, dendrolib-canary |
| [Project Assembly And Emitters](infrastructures/dendro/project-assembly-and-emitters.md) | Dendro generated C++ layout, CMake sources, and TOML parser. | CMakeLists.py, FetchContent, Dendrolib |
| [Runtime Parameter Keys](infrastructures/dendro/runtime-parameters.md) | Dendro TOML parameter keys, fallbacks, q1, and guards. | BSSN_CAHD_C, KO_DISS_SIGMA, ETA_CONST, AEH_SOLVER_FREQ |
| [ETLegacy Code Params And Loops](infrastructures/etlegacy/code-parameters-includes-and-loops.md) | ETLegacy Cactus parameter reads, symbols, and CCTK loops. | DECLARE_CCTK_PARAMETERS, CCTK_ParameterGet, cctk_lsh, cctk_nghostzones, SIMD |
| [GR ADM/BSSN, Slicing, And Matter](infrastructures/etlegacy/gr-adm-bssn-slicing-and-matter-coupling.md) | ETLegacy ADMBase import/export, slicing, and TmunuBase. | Baikal, ADM_to_BSSN, CCTK_INITIAL, ADMBase_PostInitial |
| [ETLegacy GR BSSN RHS And Constraints](infrastructures/etlegacy/gr-bssn-rhs-ricci-constraints-and-validation.md) | ETLegacy BSSN Ricci, RHS, constraints, and trusted dicts. | MoL_CalcRHS, MoL_PseudoEvolution, RbarDD |
| [GR Poststep Repairs](infrastructures/etlegacy/gr-poststep-repairs.md) | ETLegacy lapse floor and determinant/trace-free Abar repair. | enforce_detgbar_equals_detghat_trAzero, floor_the_lapse |
| [MoL, Boundaries, And Symmetries](infrastructures/etlegacy/mol-boundaries-symmetry-and-rhs-initialization.md) | ETLegacy MoL registration, symmetries, NewRad, and zero_rhss. | MoLRegisterEvolvedGroup, MoL_Register, CartGrid3D |
| [Thorn Assembly And CCL Files](infrastructures/etlegacy/thorn-assembly-and-ccl-files.md) | ETLegacy interface.ccl, param.ccl, schedule.ccl, make.code.defn. | USES_INCLUDEs, is_evol_thorn, enable_NewRad |
| [Infrastructure Code Style](infrastructures/infrastructure-code-style.md) | Infrastructure naming, CodeParameter scope, and BHaH loops. | register_CFunction_, validate_strings, parallel_codegen |
| [Commondata And PyFunction Registry](infrastructures/jax/commondata-and-pyfunction-registry.md) | JAX Commondata dataclass emission and PyFunction_dict modules. | register_commondata_params, commondata_params_dict |
| [JAX Project Layout And Emission](infrastructures/jax/project-generation-lifecycle.md) | JAX package layout, pyproject.toml emission, and module writers. | output_PyFunction_files_and_construct_project, jax_project_generator, src/ |
| [SEBOBv1 JAX Workflow](infrastructures/jax/sebobv1-jax-workflow.md) | sebobv1_jax generator, SEOBNRv5 JAX kernels, and CI call. | SEOBNRv5_aligned_spin_coefficients, PyFunction, Commondata |
| [New Infrastructure Conformance](infrastructures/new-infrastructure-conformance.md) | Conformance rules, examples, and checks for new backends. | glb_code_params_dict, CFunction_dict, BHaH_defines_h.py |
| [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md) | superB Main/Timestepping chares, proxies, SDAG, PUP, and maps. | mainchare, CProxy, CkReduction, RRMap_with_offset |
| [Diagnostics And Observables](infrastructures/superb/diagnostics-and-observables.md) | superB diagnostics, CkIO 1D/2D sessions, and reductions. | CkReductionMsg, report_sums_for_volume, out0d, out1d, out2d |
| [GR, BHaHAHA, Psi4, And Interpolation](infrastructures/superb/gr-bhahaha-psi4-and-interpolation.md) | superB InterpBufMsg messages, interpolation, BHaHAHA, and Psi4. | Interpolator3d, Horizon_finder, psi4_spinweightm2 |
| [Grids, Boundaries, MoL, And Initial Data](infrastructures/superb/grids-boundaries-mol-and-initial-data.md) | superB chare subgrids, ghost/nonlocal exchange, and MoL. | Nchare0, NGHOSTS, numerical_grids_chare |
| [superB Project Assembly](infrastructures/superb/lifecycle-and-project-assembly.md) | superB Charm++ .ci/header emission, Makefile, and charmrun. | .decl.h, .def.h, charmc |
| [Code Test Policy](validation/code-test-policy.md) | Test placement, numerical assertions, and integration rules. | pytest, doctest, coding_style.md |
| [Narrow Examples Product Exception](validation/examples-product-exception.md) | Twelve rules for nrpy/examples/** standalone test scripts. | sebob_consistency_check, sebobv2_consistency_check, dendro_application_check |
| [Expression Validation Helpers](validation/expression-validation-helpers.md) | assert_equal, check_zero, and trusted_dict module validation. | process_dictionary_of_expressions, compare_or_generate_trusted_results |
| [Generated Project CI](validation/generated-project-ci.md) | GitHub Actions workflows building and running generated codes. | main.yml, dendrolib-canary.yml, sebobv2-projected-attachment.yml, full_nrpy_local_ci.sh |
| [Static Analysis](validation/static-analysis.md) | Black, isort, mypy, Pylint, pydocstyle, darglint, and doctests. | single_file_static_analysis.sh, .pylintrc, coding_style.md, trusted_dict |
| [Test Oracles And Safe Updates](validation/test-oracles-and-safe-updates.md) | Expected-output storage, two-process review, and isolation. | validate_strings, trusted_dict, compare_or_generate_trusted_results, subprocess |
