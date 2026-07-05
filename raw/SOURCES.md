# Source Manifest

> Pointer manifest for source material. Root-level documentation sources live
> under `raw/source-docs/` so `AGENTS.md` is the only root KB document. Code,
> config, fixtures, selected logs, and build inputs stay in place. Status is
> `frozen` when a source is meant not to change and `living` when drift must
> trigger re-ingest. Last audited: 07-02-2026.

## Aggregate Sources

| Source | Provenance | Status | Ingest |
| --- | --- | --- | --- |
| `core-top-level-package-modules` | Core top-level package modules from `find nrpy -maxdepth 1 -type f \( -name '*.py' -o -name '*.txt' -o -name 'py.typed' \)`, 12 files. | living | partial |
| `helpers-package-modules` | Helper package files from `find nrpy/helpers -type f \( -name '*.py' -o -name '*.h' \)`, 21 files. | living | ingested |
| `helpers-validation-and-reference-metric-tests` | Helpers, validation helpers, and reference metric tests from `nrpy/helpers`, `nrpy/validate_expressions`, and `nrpy/tests`, 38 files. | living | partial |
| `equation-modules-and-trusted-values` | Equation modules and generated trusted-value files from `nrpy/equations`, 311 files. | living | partial |
| `infrastructure-modules-and-embedded-headers` | Infrastructure modules and embedded headers from `nrpy/infrastructures`, 365 files. | living | partial |
| `carpetx-package-inventory` | CarpetX Python package inventory from `find nrpy/infrastructures/CarpetX -type f -name '*.py'`, 26 files. | living | ingested |
| `example-generators-and-companion-scripts` | Example generators and companion scripts from `nrpy/examples`, 40 files. | living | partial |
| `ci-and-local-automation` | CI and local automation files from `.github`, 3 files. | living | partial |

## Source Documents Moved Below Root

| Source | Provenance | Status | Ingest |
| --- | --- | --- | --- |
| `raw/source-docs/original-agents.md` | Previous root agent instructions, preserved byte-for-byte before replacing `AGENTS.md`. | frozen | ingested |
| `raw/source-docs/kb-instructions.md` | KB schema and governance source moved from the repository root. | frozen | ingested |

## Cited Code And Config Sources

Exact cited files are registered below; they may also be covered by an aggregate
row above. Rows here abbreviate to source and status per `wiki/SCHEMA.md`: the
repository path is the provenance, and ingest state is tracked through covering
aggregate rows and `wiki/source-map.md`.

| Source | Status |
| --- | --- |
| `README.md` | living |
| `CITATION.md` | living |
| `coding_style.md` | living |
| `setup.py` | living |
| `requirements.txt` | living |
| `requirements-dev.txt` | living |
| `pyproject.toml` | living |
| `.mypy.ini` | living |
| `.pydocstyle` | living |
| `.pylintrc` | living |
| `.pylintrc_python36` | living |
| `.darglint` | living |
| `.github/workflows/main.yml` | living |
| `.github/single_file_static_analysis.sh` | living |
| `.github/full_nrpy_local_ci.sh` | living |
| `nrpy/c_codegen.py` | living |
| `nrpy/py_codegen.py` | living |
| `nrpy/py_function.py` | living |
| `nrpy/c_function.py` | living |
| `nrpy/grid.py` | living |
| `nrpy/params.py` | living |
| `nrpy/indexedexp.py` | living |
| `nrpy/reference_metric.py` | living |
| `nrpy/finite_difference.py` | living |
| `nrpy/helpers/expr_tree.py` | living |
| `nrpy/helpers/expression_utils.py` | living |
| `nrpy/helpers/float_to_rational.py` | living |
| `nrpy/helpers/functional.py` | living |
| `nrpy/helpers/cse_preprocess_postprocess.py` | living |
| `nrpy/helpers/custom_c_codegen_functions.py` | living |
| `nrpy/helpers/jax_printer.py` | living |
| `nrpy/helpers/simd.py` | living |
| `nrpy/helpers/simd_intrinsics.h` | living |
| `nrpy/helpers/cuda_intrinsics.h` | living |
| `nrpy/helpers/loop.py` | living |
| `nrpy/helpers/parallelization/gpu_kernel.py` | living |
| `nrpy/helpers/parallelization/utilities.py` | living |
| `nrpy/helpers/parallel_codegen.py` | living |
| `nrpy/helpers/parallelization/__init__.py` | living |
| `nrpy/helpers/generic.py` | living |
| `nrpy/helpers/cached_functions.py` | living |
| `nrpy/helpers/conditional_file_updater.py` | living |
| `nrpy/helpers/colorize_text.py` | living |
| `nrpy/helpers/type_annotation_utilities.py` | living |
| `nrpy/validate_expressions/validate_expressions.py` | living |
| `nrpy/equations/general_relativity/BSSN_RHSs.py` | living |
| `nrpy/equations/general_relativity/BSSN_quantities.py` | living |
| `nrpy/equations/general_relativity/BSSN_gauge_RHSs.py` | living |
| `nrpy/equations/general_relativity/BSSN_constraints.py` | living |
| `nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py` | living |
| `nrpy/equations/general_relativity/tests/BSSN_quantities_Cartesian.py` | living |
| `nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py` | living |
| `nrpy/infrastructures/BHaH/main_c.py` | living |
| `nrpy/infrastructures/BHaH/bhah_lib.py` | living |
| `nrpy/infrastructures/BHaH/simple_loop.py` | living |
| `nrpy/infrastructures/BHaH/griddata_commondata.py` | living |
| `nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py` | living |
| `nrpy/infrastructures/BHaH/rfm_precompute.py` | living |
| `nrpy/infrastructures/BHaH/rotation/__init__.py` | living |
| `nrpy/infrastructures/BHaH/rotation/register_all.py` | living |
| `nrpy/infrastructures/BHaH/rotation/so3_apply_R_to_vector.py` | living |
| `nrpy/infrastructures/BHaH/rotation/so3_apply_RT_to_vector.py` | living |
| `nrpy/infrastructures/BHaH/rotation/so3_apply_R_to_tensorDD.py` | living |
| `nrpy/infrastructures/BHaH/rotation/so3_axis_angle_to_R.py` | living |
| `nrpy/infrastructures/BHaH/rotation/so3_build_R_from_hats.py` | living |
| `nrpy/infrastructures/BHaH/rotation/so3_find_nU_and_dphi_from_unit_vectors.py` | living |
| `nrpy/infrastructures/BHaH/rotation/so3_left_multiply_hats_with_R.py` | living |
| `nrpy/infrastructures/BHaH/rotation/so3_relative_R_dst_from_src.py` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_apply_R_to_vector_so3_apply_R_to_vector__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_apply_RT_to_vector_so3_apply_RT_to_vector__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_apply_R_to_tensorDD_so3_apply_R_to_tensorDD__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_axis_angle_to_R_so3_axis_angle_to_R__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_build_R_from_hats_so3_build_R_from_hats__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_find_nU_and_dphi_from_unit_vectors_so3_find_nU_and_dphi_from_unit_vectors__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_left_multiply_hats_with_R_so3_left_multiply_hats_with_R__openmp.c` | living |
| `nrpy/infrastructures/BHaH/rotation/tests/so3_relative_R_dst_from_src_so3_relative_R_dst_from_src__openmp.c` | living |
| `nrpy/infrastructures/ETLegacy/interface_ccl.py` | living |
| `nrpy/infrastructures/ETLegacy/param_ccl.py` | living |
| `nrpy/infrastructures/ETLegacy/schedule_ccl.py` | living |
| `nrpy/infrastructures/ETLegacy/make_code_defn.py` | living |
| `nrpy/infrastructures/ETLegacy/CodeParameters.py` | living |
| `nrpy/infrastructures/ETLegacy/ETLegacy_include_header.py` | living |
| `nrpy/infrastructures/ETLegacy/simple_loop.py` | living |
| `nrpy/infrastructures/ETLegacy/MoL_registration.py` | living |
| `nrpy/infrastructures/ETLegacy/Symmetry_registration.py` | living |
| `nrpy/infrastructures/ETLegacy/boundary_conditions.py` | living |
| `nrpy/infrastructures/ETLegacy/zero_rhss.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/rhs_eval.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/Ricci_eval.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/ADM_to_BSSN.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/BSSN_to_ADM.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/T4DD_to_T4UU.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/RegisterSlicing.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/floor_the_lapse.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgammahat_constraint.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsFalse.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsTrue.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsFalse.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsTrue.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuFalse_KOTrue_improvementsFalse.py` | living |
| `nrpy/infrastructures/ETLegacy/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuTrue_KOTrue_improvementsFalse.py` | living |
| `nrpy/infrastructures/CarpetX/interface_ccl.py` | living |
| `nrpy/infrastructures/CarpetX/param_ccl.py` | living |
| `nrpy/infrastructures/CarpetX/schedule_ccl.py` | living |
| `nrpy/infrastructures/CarpetX/configuration_ccl.py` | living |
| `nrpy/infrastructures/CarpetX/make_code_defn.py` | living |
| `nrpy/infrastructures/CarpetX/CodeParameters.py` | living |
| `nrpy/infrastructures/CarpetX/CarpetX_include_header.py` | living |
| `nrpy/infrastructures/CarpetX/simple_loop.py` | living |
| `nrpy/infrastructures/CarpetX/boundary_conditions.py` | living |
| `nrpy/infrastructures/CarpetX/zero_rhss.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/Ricci_eval.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/rhs_eval.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/BSSN_constraints.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/ADM_to_BSSN.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/BSSN_to_ADM.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/T4DD_to_T4UU.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/floor_the_lapse.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/enforce_detgammahat_constraint.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuFalse_KOTrue_improvementsFalse.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsFalse.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_NoCovariant_Cartesian_T4munuTrue_KOTrue_improvementsFalse.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuTrue_improvementsTrue.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsTrue.py` | living |
| `nrpy/infrastructures/CarpetX/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_Cartesian_T4munuFalse_improvementsFalse.py` | living |
| `nrpy/infrastructures/JAX/commondata.py` | living |
| `nrpy/infrastructures/JAX/jax_project_generator.py` | living |
| `nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py` | living |
| `nrpy/infrastructures/superB/main_chare.py` | living |
| `nrpy/infrastructures/superB/timestepping_chare.py` | living |
| `nrpy/infrastructures/superB/Makefile_helpers.py` | living |
| `nrpy/infrastructures/superB/numerical_grids.py` | living |
| `nrpy/infrastructures/superB/chare_communication_maps.py` | living |
| `nrpy/infrastructures/superB/CurviBoundaryConditions.py` | living |
| `nrpy/infrastructures/superB/MoL.py` | living |
| `nrpy/infrastructures/superB/initial_data.py` | living |
| `nrpy/infrastructures/superB/BHaH_implementation.py` | living |
| `nrpy/infrastructures/superB/horizon_finder_chare.py` | living |
| `nrpy/infrastructures/superB/interpolator3d_chare.py` | living |
| `nrpy/infrastructures/superB/diagnostics/diagnostics.py` | living |
| `nrpy/infrastructures/superB/diagnostics/diagnostics_nearest_grid_center.py` | living |
| `nrpy/infrastructures/superB/diagnostics/diagnostics_nearest_1d_y_and_z_axes.py` | living |
| `nrpy/infrastructures/superB/diagnostics/diagnostics_nearest_2d_xy_and_yz_planes.py` | living |
| `nrpy/infrastructures/superB/general_relativity/diagnostics_nearest.py` | living |
| `nrpy/infrastructures/superB/general_relativity/psi4_spinweightm2_decomposition.py` | living |
| `nrpy/infrastructures/superB/nrpyelliptic/diagnostics_nearest.py` | living |
| `nrpy/infrastructures/superB/superB/superB_pup.py` | living |
| `nrpy/infrastructures/superB/superB/superB.h` | living |
| `nrpy/infrastructures/superB/superB/superB_pup_function_prototypes.h` | living |
| `nrpy/tests/reference_metric_Cartesian.py` | living |
| `nrpy/tests/reference_metric_Spherical.py` | living |
| `nrpy/tests/reference_metric_GeneralRFM_fisheyeN2.py` | living |
| `nrpy/examples/wave_equation_cartesian.py` | living |
| `nrpy/examples/wave_equation_curvilinear.py` | living |
| `nrpy/examples/wave_equation_multicoordinates.py` | living |
| `nrpy/examples/sebobv1_jax.py` | living |
| `nrpy/examples/sebobv2.py` | living |
| `nrpy/examples/seobnrv5_aligned_spin_inspiral.py` | living |
| `nrpy/examples/nrpypn_quasicircular_momenta.py` | living |
| `nrpy/examples/two_blackholes_collide.py` | living |
| `nrpy/examples/blackhole_spectroscopy.py` | living |
| `nrpy/examples/spinning_blackhole.py` | living |
| `nrpy/examples/kasner_exact_evolution.py` | living |
| `nrpy/examples/nrpyelliptic_conformally_flat.py` | living |
| `nrpy/examples/carpet_wavetoy_thorns.py` | living |
| `nrpy/examples/carpet_baikal_thorns.py` | living |
| `nrpy/examples/carpetx_wavetoy_thorns.py` | living |
| `nrpy/examples/carpetx_baikal_thorns.py` | living |
| `nrpy/examples/superB_two_blackholes_collide.py` | living |
| `nrpy/examples/superB_blackhole_spectroscopy.py` | living |
| `nrpy/examples/superB_nrpyelliptic_conformally_flat.py` | living |
| `nrpy/examples/mass_geodesic_integrator.py` | living |
| `nrpy/examples/photon_geodesic_integrator.py` | living |
| `nrpy/examples/photon_geodesic_batch_integrator.py` | living |
| `nrpy/examples/tovola_neutron_star.py` | living |
| `nrpy/examples/hydro_without_hydro.py` | living |
| `nrpy/examples/groovy_TOV_BSSN.py` | living |
| `nrpy/examples/manga_bhah_lib.py` | living |
| `nrpy/examples/bhahaha.py` | living |
| `nrpy/examples/et_WaveToyfiles/ThornList` | living |
| `nrpy/examples/et_WaveToyfiles/WaveToyNRPy.par` | living |
| `nrpy/examples/et_WaveToyfiles/test/test.ccl` | living |
| `nrpy/examples/et_WaveToyfiles/test/WaveToyNRPy_test.par` | living |
| `nrpy/examples/et_WaveToyfiles/test/WaveToyNRPy_test/uuGF.x.asc` | living |
| `nrpy/examples/geodesic_visualizations/visualize_trajectory.py` | living |
| `nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py` | living |
| `nrpy/examples/geodesic_visualizations/render_lensed_image.py` | living |
| `nrpy/examples/geodesic_visualizations/visualize_lensed_image.py` | living |
| `nrpy/examples/geodesic_visualizations/blueprint_analysis.py` | living |
| `nrpy/examples/tests/sebob_consistency_check.py` | living |
| `nrpy/examples/tests/sebobv2_consistency_check.py` | living |
| `nrpy/equations/basis_transforms/jacobians.py` | living |
| `nrpy/equations/basis_transforms/tests/jacobians_Cartesian.py` | living |
| `nrpy/equations/basis_transforms/tests/jacobians_GeneralRFM_fisheyeN2.py` | living |
| `nrpy/equations/basis_transforms/tests/jacobians_Spherical.py` | living |
| `nrpy/equations/general_relativity/ADM_to_BSSN.py` | living |
| `nrpy/equations/general_relativity/BSSN_to_ADM.py` | living |
| `nrpy/equations/general_relativity/BSSN_to_g4Christoffel.py` | living |
| `nrpy/equations/general_relativity/InitialData_Cartesian.py` | living |
| `nrpy/equations/general_relativity/InitialData_Spherical.py` | living |
| `nrpy/equations/general_relativity/LorentzBoost.py` | living |
| `nrpy/equations/general_relativity/T4munu.py` | living |
| `nrpy/equations/general_relativity/bhahaha/ExpansionFunctionTheta.py` | living |
| `nrpy/equations/general_relativity/bhahaha/HorizonSpinVorticityDipole.py` | living |
| `nrpy/equations/general_relativity/bhahaha/SpECTRESpinEstimate.py` | living |
| `nrpy/equations/general_relativity/bhahaha/approx_killing_vector_spin.py` | living |
| `nrpy/equations/general_relativity/bhahaha/area.py` | living |
| `nrpy/equations/general_relativity/bhahaha/tests/ExpansionFunctionTheta_Spherical.py` | living |
| `nrpy/equations/general_relativity/bhahaha/tests/HorizonSpinVorticityDipole_Spherical.py` | living |
| `nrpy/equations/general_relativity/bhahaha/tests/SpECTRESpinEstimate_Spherical.py` | living |
| `nrpy/equations/general_relativity/bhahaha/tests/approx_killing_vector_spin_Spherical.py` | living |
| `nrpy/equations/general_relativity/bhahaha/tests/area_Spherical.py` | living |
| `nrpy/equations/general_relativity/fishbone_moncrief/fishbone_moncrief.py` | living |
| `nrpy/equations/general_relativity/fishbone_moncrief/tests/fishbone_moncrief.py` | living |
| `nrpy/equations/general_relativity/g4munu_conversions.py` | living |
| `nrpy/equations/general_relativity/geodesics/analytic_spacetimes.py` | living |
| `nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/conserved_quantities.py` | living |
| `nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/tests/conserved_quantities_KerrSchild_Cartesian_massive.py` | living |
| `nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/tests/conserved_quantities_KerrSchild_Cartesian_photon.py` | living |
| `nrpy/equations/general_relativity/geodesics/geodesics.py` | living |
| `nrpy/equations/general_relativity/geodesics/tests/analytic_spacetimes_KerrSchild_Cartesian.py` | living |
| `nrpy/equations/general_relativity/geodesics/tests/geodesics_KerrSchild_Cartesian_massive.py` | living |
| `nrpy/equations/general_relativity/geodesics/tests/geodesics_KerrSchild_Cartesian_photon.py` | living |
| `nrpy/equations/general_relativity/psi4.py` | living |
| `nrpy/equations/general_relativity/psi4_tetrads.py` | living |
| `nrpy/equations/general_relativity/tests/ADM_to_BSSN_StaticTrumpet.py` | living |
| `nrpy/equations/general_relativity/tests/BSSN_to_ADM_Cartesian.py` | living |
| `nrpy/equations/general_relativity/tests/BSSN_to_g4Christoffel_Cartesian.py` | living |
| `nrpy/equations/general_relativity/tests/InitialData_Cartesian_BrillLindquist.py` | living |
| `nrpy/equations/general_relativity/tests/InitialData_Cartesian_Kasner.py` | living |
| `nrpy/equations/general_relativity/tests/InitialData_Spherical_OffsetKerrSchild.py` | living |
| `nrpy/equations/general_relativity/tests/InitialData_Spherical_StaticTrumpet.py` | living |
| `nrpy/equations/general_relativity/tests/InitialData_Spherical_UIUCBlackHole.py` | living |
| `nrpy/equations/general_relativity/tests/LorentzBoost.py` | living |
| `nrpy/equations/general_relativity/tests/T4munu.py` | living |
| `nrpy/equations/general_relativity/tests/g4munu_conversions.py` | living |
| `nrpy/equations/general_relativity/tests/psi4_leave_symbolic_Spherical.py` | living |
| `nrpy/equations/general_relativity/tests/psi4_quasiKinnersley_SinhSpherical_rfm_precompute.py` | living |
| `nrpy/equations/general_relativity/tests/psi4_quasiKinnersley_Spherical.py` | living |
| `nrpy/equations/general_relativity/tests/psi4_tetrads_Spherical.py` | living |
| `nrpy/equations/generalrfm/fisheye.py` | living |
| `nrpy/equations/generalrfm/tests/fisheye_N1.py` | living |
| `nrpy/equations/generalrfm/tests/fisheye_N2.py` | living |
| `nrpy/equations/grhd/GRHD_equations.py` | living |
| `nrpy/equations/grhd/HLL_fluxes.py` | living |
| `nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py` | living |
| `nrpy/equations/grhd/characteristic_speeds.py` | living |
| `nrpy/equations/grhd/tests/GRHD_equations_Cartesian.py` | living |
| `nrpy/equations/grhd/tests/GRHD_equations_SinhSpherical_rfm_precompute.py` | living |
| `nrpy/equations/grhd/tests/GRHD_equations_Spherical.py` | living |
| `nrpy/equations/grhd/tests/HLL_fluxes.py` | living |
| `nrpy/equations/grhd/tests/Min_Max_and_Piecewise_Expressions.py` | living |
| `nrpy/equations/grhd/tests/characteristic_speeds.py` | living |
| `nrpy/equations/nrpyelliptic/ConformallyFlat_RHSs.py` | living |
| `nrpy/equations/nrpyelliptic/ConformallyFlat_SourceTerms.py` | living |
| `nrpy/equations/nrpyelliptic/tests/ConformallyFlat_RHSs_Cartesian.py` | living |
| `nrpy/equations/nrpyelliptic/tests/ConformallyFlat_RHSs_Spherical.py` | living |
| `nrpy/equations/nrpyelliptic/tests/ConformallyFlat_SourceTerms_Cartesian.py` | living |
| `nrpy/equations/nrpyelliptic/tests/ConformallyFlat_SourceTerms_Spherical.py` | living |
| `nrpy/equations/quaternion_rotations/tensor_rotation.py` | living |
| `nrpy/equations/rotation/SO3_rotations.py` | living |
| `nrpy/equations/rotation/tests/SO3_rotations.py` | living |
| `nrpy/equations/seobnr/BOB_aligned_spin_waveform_quantities.py` | living |
| `nrpy/equations/seobnr/BOB_aligned_spin_waveform_quantities_higher_modes.py` | living |
| `nrpy/equations/seobnr/BOB_v2_waveform_quantities_kankani_etal.py` | living |
| `nrpy/equations/seobnr/SEOBNRv5_aligned_spin_Hamiltonian.py` | living |
| `nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py` | living |
| `nrpy/equations/seobnr/SEOBNRv5_aligned_spin_merger_quantities.py` | living |
| `nrpy/equations/seobnr/SEOBNRv5_aligned_spin_waveform_quantities.py` | living |
| `nrpy/equations/seobnr/SEOBNRv5_coprecessing_rotations_quantities.py` | living |
| `nrpy/equations/seobnr/SEOBNRv5_merger_ringdown.py` | living |
| `nrpy/equations/seobnr/SEOBNRv5_quasi_precessing_spin_Hamiltonian.py` | living |
| `nrpy/equations/seobnr/SEOBNRv5_spin_evolution_equations.py` | living |
| `nrpy/equations/seobnr/tests/BOB_aligned_spin_waveform_quantities.py` | living |
| `nrpy/equations/seobnr/tests/BOB_aligned_spin_waveform_quantities_higher_modes.py` | living |
| `nrpy/equations/seobnr/tests/BOB_v2_waveform_quantities_kankani_etal.py` | living |
| `nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_Hamiltonian.py` | living |
| `nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_constants.py` | living |
| `nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_merger_quantities.py` | living |
| `nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_waveform_quantities.py` | living |
| `nrpy/equations/seobnr/tests/SEOBNRv5_coprecessing_rotations_quantities.py` | living |
| `nrpy/equations/seobnr/tests/SEOBNRv5_merger_ringdown.py` | living |
| `nrpy/equations/seobnr/tests/SEOBNRv5_quasi_precessing_spin_Hamiltonian.py` | living |
| `nrpy/equations/seobnr/tests/SEOBNRv5_spin_evolution_equations.py` | living |
| `nrpy/equations/special_functions/spin_weighted_spherical_harmonics.py` | living |
| `nrpy/equations/special_functions/tests/spin_weighted_spherical_harmonics.py` | living |
| `nrpy/equations/tov/TOV_equations.py` | living |
| `nrpy/equations/tov/tests/TOV_equations.py` | living |
| `nrpy/equations/wave_equation/WaveEquationCurvilinear_RHSs.py` | living |
| `nrpy/equations/wave_equation/WaveEquation_RHSs.py` | living |
| `nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py` | living |
| `nrpy/equations/wave_equation/tests/WaveEquationCurvilinear_RHSs_Spherical.py` | living |
| `nrpy/equations/wave_equation/tests/WaveEquation_RHSs_WaveEquation.py` | living |
| `nrpy/equations/wave_equation/tests/WaveEquation_Solutions_InitialData_PlaneWave.py` | living |
| `nrpy/equations/wave_equation/tests/WaveEquation_Solutions_InitialData_SphericalGaussian.py` | living |

## External Background Sources

| Source | Provenance | Status | Accessed | Ingest | Notes |
| --- | --- | --- | --- | --- | --- |
| `https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md` | Andrej Karpathy gist raw note, `LLM Wiki`, pinned to revision `ac46de1ad27f92b28ac95459c782c07f6b8c964a`. | frozen | 06-30-2026 | partial | Background approach source for persistent LLM-maintained wiki governance: raw/wiki/schema layers, index/log, and ingest/query/lint workflows. |
| `https://arxiv.org/abs/2111.02424` | arXiv abstract page for arXiv:2111.02424. | living | 06-29-2026 | partial | NRPyElliptic background for hyperbolic relaxation and conformally flat binary-puncture initial data. |
| `https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html` | Intel Intrinsics Guide landing page. | living | 06-29-2026 | partial | SSE, AVX, AVX512, and intrinsic-family terminology for SIMD helper documentation. |
| `https://docs.nvidia.com/cuda/cuda-programming-guide/index.html` | NVIDIA CUDA C++ Programming Guide. | living | 06-29-2026 | partial | CUDA C++ terminology for SIMD, CUDA header, and GPU kernel helper documentation. |
| `https://link.aps.org/doi/10.1103/PhysRev.55.374` | DOI landing page for `10.1103/PhysRev.55.374`. | living | 06-29-2026 | partial | Original Oppenheimer-Volkoff stellar-equilibrium background. |
| `https://web2.ph.utexas.edu/~gsudama/pub/1967_008.pdf` | PDF URL for Goldberg et al. spin-weighted spherical-harmonic formula reference. | living | 06-29-2026 | partial | Goldberg-formula background for spin-weighted spherical harmonics. |
| `https://pubs.aip.org/aip/jmp/article/57/9/092504/648118/How-should-spin-weighted-spherical-functions-be` | Journal of Mathematical Physics article landing page. | living | 06-29-2026 | partial | Background on spin-weighted functions and quaternion viewpoints. |
| `https://rotations.berkeley.edu/geodesics-of-the-rotation-group-so3/` | Berkeley rotations course page. | living | 06-29-2026 | partial | Background on SO(3) and quaternion rotation geometry. |
| `https://github.com/charmplusplus/charm/blob/main/doc/quickstart.rst` | Charm++ quickstart from the Charm++ repository. | living | 06-29-2026 | partial | Background for `.ci` files, generated `.decl.h` and `.def.h` files, `charmc`, and `charmrun`. |
| `https://github.com/charmplusplus/charm/blob/main/doc/charm%2B%2B/manual.rst` | Charm++ language manual from the Charm++ repository. | living | 06-29-2026 | partial | Background for chares, entry methods, proxies, SDAG, PUP, checkpoint/restart, reductions, and chare arrays. |
| `https://github.com/charmplusplus/charm/blob/main/doc/libraries/manual.rst` | Charm++ and Converse libraries manual from the Charm++ repository. | living | 06-29-2026 | partial | Background for CkIO. |
| `https://einsteintoolkit.org/usersguide/UsersGuide.html` | Cactus 4.20 Users Guide page from the Einstein Toolkit site. | living | 06-30-2026 | partial | Background for Cactus thorn-writing and Cactus terminology emitted by ETLegacy code. |
| `https://einsteintoolkit.org/referencemanual/ReferenceManual.html` | Cactus 4.20 Reference Manual page from the Einstein Toolkit site. | living | 06-30-2026 | partial | Background for `CCTK_*` thorn-writer function terminology emitted by ETLegacy code. |
| `https://www.cactuscode.org/documentation/usersguide/UsersGuidech9.html` | Cactus users guide chapter C1, `Application thorns`. | living | 06-30-2026 | partial | Background for CCL and thorn-file terminology emitted by ETLegacy code. |
| `https://www.cactuscode.org/documentation/usersguide/UsersGuidech12.html` | Cactus users guide chapter D2, CCL reference. | living | 06-30-2026 | partial | Background for CCL configuration-file terminology emitted by CarpetX thorn-assembly code. |
| `https://einsteintoolkit.org/thornguide/CarpetX/CarpetX/documentation.html` | Einstein Toolkit thorn guide page, `CarpetX`. | living | 06-30-2026 | partial | Background for CarpetX CCL and dependency terminology emitted by CarpetX thorn-assembly code. |
| `https://einsteintoolkit.org/thornguide/CactusNumerical/MoL/documentation.html` | Einstein Toolkit thorn guide page, `Method of Lines`. | living | 06-30-2026 | partial | Background for MoL terminology emitted by ETLegacy registration code. |
| `https://einsteintoolkit.org/thornguide/CactusBase/Boundary/documentation.html` | Einstein Toolkit thorn guide page, `Boundary Conditions`. | living | 06-30-2026 | partial | Background for Boundary terminology emitted by ETLegacy boundary-condition code. |
| `https://einsteintoolkit.org/thornguide/EinsteinEvolve/NewRad/documentation.html` | Einstein Toolkit thorn guide page, `NewRad`. | living | 06-30-2026 | partial | Background for NewRad terminology emitted by ETLegacy boundary-condition code. |
| `https://einsteintoolkit.org/thornguide/CactusBase/CartGrid3D/documentation.html` | Einstein Toolkit thorn guide page, `CartGrid3D`. | living | 06-30-2026 | partial | Background for CartGrid3D terminology emitted by ETLegacy symmetry code. |
| `https://einsteintoolkit.org/thornguide/EinsteinBase/ADMBase/documentation.html` | Einstein Toolkit thorn guide page, `ADMBase`. | living | 06-30-2026 | partial | Background for ADMBase terminology emitted by ETLegacy GR coupling code. |
| `https://einsteintoolkit.org/thornguide/EinsteinBase/TmunuBase/documentation.html` | Einstein Toolkit thorn guide page, `TmunuBase`. | living | 06-30-2026 | partial | Background for TmunuBase terminology emitted by ETLegacy matter-coupling code. |

## Exclusions

- `project/`, `build/`, and `dist/`.
- Caches such as `__pycache__/`, `.mypy_cache/`, and `.pytest_cache/`.
- Compiled files, object files, executables, generated PDFs, archives, images,
  rendered artifacts, and other non-text outputs.
- Generated C/CUDA projects and generated thorns unless registered as selected
  frozen evidence.
- `optimal-plan-runs/`, `plan*.md`, `synth_plan*.md`, and `rank*.md`.
- Scratch logs, prompt transcripts, token-count reports, and latest-snapshot
  reports.
