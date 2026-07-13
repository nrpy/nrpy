# Source Manifest

> Pointer manifest for source material. Root-level documentation sources live
> under `raw/source-docs/` so `AGENTS.md` is the only root KB document. Code,
> config, fixtures, selected logs, and build inputs stay in place. Status is
> `frozen` when a source is meant not to change and `living` when drift must
> trigger re-ingest. Last audited: 07-12-2026.

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
| `raw/source-docs/original-agents.md` | Previous root agent instructions, preserved byte-for-byte before replacing `AGENTS.md`. Its root style-guide reference is superseded; see the source map row for current routing. | frozen | ingested |
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
| `bin/nrpyinline.py` | living |
| `coding_style.md` | living |
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
| `nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py` | living |
| `nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py` | living |
| `nrpy/equations/general_relativity/tests/BSSN_quantities_Cartesian.py` | living |
| `nrpy/equations/general_relativity/tests/BSSN_constraints_Cartesian.py` | living |
| `nrpy/infrastructures/BHaH/main_c.py` | living |
| `nrpy/infrastructures/BHaH/bhah_lib.py` | living |
| `nrpy/infrastructures/BHaH/Makefile_helpers.py` | living |
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

## External Sources

| Source | Provenance | Status | Accessed | Ingest | Notes |
| --- | --- | --- | --- | --- | --- |
| `https://gist.githubusercontent.com/karpathy/442a6bf555914893e9891c11519de94f/raw/ac46de1ad27f92b28ac95459c782c07f6b8c964a/llm-wiki.md` | Andrej Karpathy gist raw note, `LLM Wiki`, pinned to revision `ac46de1ad27f92b28ac95459c782c07f6b8c964a`. | frozen | 06-30-2026 | partial | Background approach source for persistent LLM-maintained wiki governance: raw/wiki/schema layers and ingest/query/lint workflows. |
| `https://arxiv.org/abs/1605.01938` | arXiv abstract page for HBR2016 final-spin paper. | living | 07-13-2026 | partial | Mutable latest-revision background page; exact revision-to-fit mapping is not yet audited. |
| `https://arxiv.org/abs/1611.00332` | arXiv abstract page for UIB2016 final-state paper and ancillary implementation. | living | 07-13-2026 | partial | Mutable latest-revision background page; exact revision-to-fit/ancillary mapping is not yet audited. |
| `https://arxiv.org/abs/2111.02424` | arXiv abstract page for arXiv:2111.02424. | living | 07-13-2026 | partial | Mutable latest-revision background page; exact revision-to-claim mapping is not yet audited. |
| `https://arxiv.org/abs/2303.18039` | arXiv abstract page for SEOBNRv5HM aligned-spin waveform model. | living | 07-13-2026 | partial | Mutable latest-revision background page; exact revision-to-equation mapping is not yet audited. |
| `https://arxiv.org/abs/2303.18046` | arXiv abstract page for SEOBNRv5PHM precessing waveform model. | living | 07-13-2026 | partial | Mutable latest-revision background page; exact revision-to-equation mapping is not yet audited. |
| `https://arxiv.org/abs/2303.18143` | arXiv abstract page for SEOBNRv5 precessing dynamics groundwork. | living | 07-13-2026 | partial | Mutable latest-revision background page; exact revision-to-equation mapping is not yet audited. |
| `https://arxiv.org/abs/2508.20418` | arXiv abstract page for SEBOB aligned-spin attachment model. | living | 07-13-2026 | partial | Mutable latest-revision background page; exact revision-to-equation mapping is not yet audited. |
| `https://arxiv.org/abs/2510.25012` | arXiv abstract page for BOBv2 waveform model. | living | 07-13-2026 | partial | Mutable latest-revision background page; exact revision-to-equation/validation mapping is not yet audited. |
| `https://pypi.org/pypi/black/26.5.1/json` | PyPI release metadata for Black 26.5.1. | frozen | 07-12-2026 | partial | Official release metadata for Black's Python 3.10+ runtime requirement and installation surface. |
| `https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html` | Intel Intrinsics Guide landing page. | living | 07-12-2026 | partial | Official SSE, AVX, AVX512, and intrinsic-family terminology for SIMD helper documentation; not authority for NRPy macro behavior. |
| `https://docs.nvidia.com/cuda/cuda-programming-guide/index.html` | NVIDIA CUDA Programming Guide. | living | 07-12-2026 | partial | Official CUDA programming-model terminology for helper documentation; not authority for NRPy kernel behavior. |
| `https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/nvcc.html` | NVIDIA CUDA Programming Guide, NVCC chapter. | living | 07-12-2026 | partial | Official NVCC compiler-driver terminology used by CUDA example prerequisites. |
| `https://docs.einsteintoolkit.org/et-docs/Adding_a_test_case` | Einstein Toolkit documentation, `Adding a test case`. | living | 07-12-2026 | partial | Official test-suite terminology; local workflow configuration decides NRPy CI behavior. |
| `https://docs.jax.dev/en/latest/installation.html` | JAX installation documentation. | living | 07-12-2026 | partial | Official supported-platform and installation context; local requirements and generator code decide NRPy behavior. |
| `https://www.gnu.org/software/gsl/doc/html/usage.html` | GNU GSL 2.8 documentation, `Using the Library`. | living | 07-12-2026 | partial | Official compiling/linking context for GSL-dependent generated examples. |
| `https://link.aps.org/doi/10.1103/PhysRev.55.374` | DOI landing page for `10.1103/PhysRev.55.374`. | living | 06-29-2026 | partial | Original Oppenheimer-Volkoff stellar-equilibrium background. |
| `https://web2.ph.utexas.edu/~gsudama/pub/1967_008.pdf` | PDF URL for Goldberg et al. spin-weighted spherical-harmonic formula reference. | living | 06-29-2026 | partial | Goldberg-formula background for spin-weighted spherical harmonics. |
| `https://pubs.aip.org/aip/jmp/article/57/9/092504/648118/How-should-spin-weighted-spherical-functions-be` | Journal of Mathematical Physics article landing page. | living | 06-29-2026 | partial | Background on spin-weighted functions and quaternion viewpoints. |
| `https://rotations.berkeley.edu/geodesics-of-the-rotation-group-so3/` | Berkeley rotations course page. | living | 06-29-2026 | partial | Background on SO(3) and quaternion rotation geometry. |
| `https://charm.readthedocs.io/en/v8.0.0/quickstart.html` | Charm++ quickstart from the Charm++ repository. | frozen | 07-12-2026 | partial | External spec context for `.ci` files, generated `.decl.h` and `.def.h` files, `charmc`, and `charmrun`. |
| `https://charm.readthedocs.io/en/v8.0.0/charm%2B%2B/manual.html` | Charm++ language manual from the Charm++ repository. | frozen | 07-12-2026 | partial | External spec context for chares, entry methods, proxies, SDAG, PUP, checkpoint/restart, reductions, chare arrays, and generic message-transport caveats. |
| `https://charm.readthedocs.io/en/v8.0.0/libraries/manual.html` | Charm++ and Converse libraries manual from the Charm++ repository. | frozen | 07-12-2026 | partial | External spec context for CkIO and Charm++ library APIs. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/scripts/charmc` | Charm++ pinned Charm++ 8.0.0 compiler-driver script. | frozen | 07-12-2026 | partial | Background implementation context for `charmc`; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/scripts/Makefile` | Charm++ pinned Charm++ 8.0.0 build-script Makefile. | frozen | 07-12-2026 | partial | Background implementation context for Charm++ build scripts; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/scripts/Make.cidepends` | Charm++ pinned Charm++ 8.0.0 `.ci` dependency helper. | frozen | 07-12-2026 | partial | Background implementation context for `.ci` dependency generation; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/xi-main.C` | Charm++ pinned Charm++ 8.0.0 charmxi translator main source. | frozen | 07-12-2026 | partial | Background implementation context for interface translation; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/xi-Module.C` | Charm++ pinned Charm++ 8.0.0 charmxi module translation source. | frozen | 07-12-2026 | partial | Background implementation context for modules and generated interface files; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/libs/ck-libs/io/Makefile` | Charm++ pinned Charm++ 8.0.0 CkIO library Makefile. | frozen | 07-12-2026 | partial | Background implementation context for CkIO linking; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/tests/charm%2B%2B/io/Makefile` | Charm++ pinned Charm++ 8.0.0 CkIO test Makefile. | frozen | 07-12-2026 | partial | Background test/example context for `-module CkIO`; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/tests/charm%2B%2B/io_read/Makefile` | Charm++ pinned Charm++ 8.0.0 CkIO read-test Makefile. | frozen | 07-12-2026 | partial | Background test/example context for CkIO test builds; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/tests/charm%2B%2B/charmxi_parsing/Makefile` | Charm++ pinned Charm++ 8.0.0 charmxi parsing test Makefile. | frozen | 07-12-2026 | partial | Background test/example context for charmxi parsing builds; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/xi-grammar.y` | Charm++ pinned Charm++ 8.0.0 charmxi grammar source. | frozen | 07-12-2026 | partial | Background implementation context for `.ci` grammar terms; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/xi-Chare.C` | Charm++ pinned Charm++ 8.0.0 charmxi chare translation source. | frozen | 07-12-2026 | partial | Background implementation context for chare and chare-array generation; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/xi-Entry.C` | Charm++ pinned Charm++ 8.0.0 charmxi entry-method translation source. | frozen | 07-12-2026 | partial | Background implementation context for entry methods and generated indexes; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/xi-Parameter.C` | Charm++ pinned Charm++ 8.0.0 charmxi parameter-marshalling source. | frozen | 07-12-2026 | partial | Background implementation context for parameter-marshaled entry methods; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/xi-Member.C` | Charm++ pinned Charm++ 8.0.0 charmxi member translation source. | frozen | 07-12-2026 | partial | Background implementation context for generated class/member support; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckarray.h` | Charm++ pinned Charm++ 8.0.0 chare-array header. | frozen | 07-12-2026 | partial | Background implementation context for chare arrays, proxies, and maps; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckarray.C` | Charm++ pinned Charm++ 8.0.0 chare-array implementation. | frozen | 07-12-2026 | partial | Background implementation context for chare-array runtime behavior; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckarrayindex.h` | Charm++ pinned Charm++ 8.0.0 chare-array index header. | frozen | 07-12-2026 | partial | Background implementation context for `CkArrayIndex3D` and index objects; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckcallback.h` | Charm++ pinned Charm++ 8.0.0 callback header. | frozen | 07-12-2026 | partial | Background implementation context for `CkCallback`; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckcallback.C` | Charm++ pinned Charm++ 8.0.0 callback implementation. | frozen | 07-12-2026 | partial | Background implementation context for callback dispatch; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/readonly.h` | Charm++ pinned Charm++ 8.0.0 readonly support header. | frozen | 07-12-2026 | partial | Background implementation context for readonly globals; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/sdag/CParsedFile.C` | Charm++ pinned Charm++ 8.0.0 SDAG parsed-file source. | frozen | 07-12-2026 | partial | Background implementation context for SDAG translation; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/sdag/CEntry.C` | Charm++ pinned Charm++ 8.0.0 SDAG entry translation source. | frozen | 07-12-2026 | partial | Background implementation context for SDAG entry behavior; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/sdag/constructs/When.C` | Charm++ pinned Charm++ 8.0.0 SDAG `when` construct source. | frozen | 07-12-2026 | partial | Background implementation context for SDAG `when` matching; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/xlat-i/sdag/constructs/SdagEntry.C` | Charm++ pinned Charm++ 8.0.0 SDAG entry construct source. | frozen | 07-12-2026 | partial | Background implementation context for SDAG entry generation; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/sdag.h` | Charm++ pinned Charm++ 8.0.0 SDAG runtime header. | frozen | 07-12-2026 | partial | Background implementation context for SDAG runtime buffering; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/sdag.C` | Charm++ pinned Charm++ 8.0.0 SDAG runtime implementation. | frozen | 07-12-2026 | partial | Background implementation context for SDAG runtime state; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/charm%2B%2B.h` | Charm++ pinned Charm++ 8.0.0 core Charm++ API header. | frozen | 07-12-2026 | partial | Background implementation context for core Charm++ runtime APIs; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/util/pup.h` | Charm++ pinned Charm++ 8.0.0 PUP header. | frozen | 07-12-2026 | partial | Background implementation context for PUP serialization; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckcheckpoint.h` | Charm++ pinned Charm++ 8.0.0 checkpoint header. | frozen | 07-12-2026 | partial | Background implementation context for checkpoint/restart APIs; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckcheckpoint.C` | Charm++ pinned Charm++ 8.0.0 checkpoint implementation. | frozen | 07-12-2026 | partial | Background implementation context for checkpoint/restart behavior; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckcheckpoint.ci` | Charm++ pinned Charm++ 8.0.0 checkpoint interface file. | frozen | 07-12-2026 | partial | Background implementation context for checkpoint entry methods; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/init.C` | Charm++ pinned Charm++ 8.0.0 runtime initialization source. | frozen | 07-12-2026 | partial | Background implementation context for runtime startup and restart paths; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/cklocation.h` | Charm++ pinned Charm++ 8.0.0 location-manager header. | frozen | 07-12-2026 | partial | Background implementation context for location, migration, and placement caveats; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/cklocation.C` | Charm++ pinned Charm++ 8.0.0 location-manager implementation. | frozen | 07-12-2026 | partial | Background implementation context for location-manager PUP and migration behavior; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckmemcheckpoint.h` | Charm++ pinned Charm++ 8.0.0 memory-checkpoint header. | frozen | 07-12-2026 | partial | Background implementation context for memory checkpointing; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckmemcheckpoint.C` | Charm++ pinned Charm++ 8.0.0 memory-checkpoint implementation. | frozen | 07-12-2026 | partial | Background implementation context for memory checkpointing; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckmemcheckpoint.ci` | Charm++ pinned Charm++ 8.0.0 memory-checkpoint interface file. | frozen | 07-12-2026 | partial | Background implementation context for memory-checkpoint entry methods; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-cp/controlPoints.h` | Charm++ pinned Charm++ 8.0.0 control-points header. | frozen | 07-12-2026 | partial | Background implementation context for ck-cp caveats; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-cp/controlPoints.C` | Charm++ pinned Charm++ 8.0.0 control-points implementation. | frozen | 07-12-2026 | partial | Background implementation context for ck-cp behavior; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-cp/pathHistory.h` | Charm++ pinned Charm++ 8.0.0 control-points path-history header. | frozen | 07-12-2026 | partial | Background implementation context for ck-cp path history; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-cp/pathHistory.C` | Charm++ pinned Charm++ 8.0.0 control-points path-history implementation. | frozen | 07-12-2026 | partial | Background implementation context for ck-cp path history; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckreduction.h` | Charm++ pinned Charm++ 8.0.0 reduction header. | frozen | 07-12-2026 | partial | Background implementation context for `CkReduction` and `CkReductionMsg`; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckreduction.C` | Charm++ pinned Charm++ 8.0.0 reduction implementation. | frozen | 07-12-2026 | partial | Background implementation context for built-in reducers; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/libs/ck-libs/io/ckio.h` | Charm++ pinned Charm++ 8.0.0 CkIO header. | frozen | 07-12-2026 | partial | Background implementation context for CkIO API shape; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/libs/ck-libs/io/ckio.C` | Charm++ pinned Charm++ 8.0.0 CkIO implementation. | frozen | 07-12-2026 | partial | Background implementation context for CkIO sessions and callbacks; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/libs/ck-libs/io/ckio.ci` | Charm++ pinned Charm++ 8.0.0 CkIO interface file. | frozen | 07-12-2026 | partial | Background implementation context for CkIO entry methods; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckarrayoptions.h` | Charm++ pinned Charm++ 8.0.0 chare-array options header. | frozen | 07-12-2026 | partial | Background implementation context for array-map and placement options; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckarrayoptions.C` | Charm++ pinned Charm++ 8.0.0 chare-array options implementation. | frozen | 07-12-2026 | partial | Background implementation context for array options; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/cklocation.ci` | Charm++ pinned Charm++ 8.0.0 location-manager interface file. | frozen | 07-12-2026 | partial | Background implementation context for location-manager entry methods; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-core/ckmigratable.h` | Charm++ pinned Charm++ 8.0.0 migratable-object header. | frozen | 07-12-2026 | partial | Background implementation context for migration caveats; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-ldb/LBManager.h` | Charm++ pinned Charm++ 8.0.0 load-balancer manager header. | frozen | 07-12-2026 | partial | Background implementation context for ck-ldb and AtSync caveats; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-ldb/LBManager.C` | Charm++ pinned Charm++ 8.0.0 load-balancer manager implementation. | frozen | 07-12-2026 | partial | Background implementation context for ck-ldb and AtSync behavior; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-ldb/CommonLBs.ci` | Charm++ pinned Charm++ 8.0.0 common load-balancers interface file. | frozen | 07-12-2026 | partial | Background implementation context for load-balancer entry methods; pinned to Charm++ 8.0.0. |
| `https://github.com/charmplusplus/charm/blob/v8.0.0/src/ck-ldb/EveryLB.ci` | Charm++ pinned Charm++ 8.0.0 EveryLB interface file. | frozen | 07-12-2026 | partial | Background implementation context for load-balancer entry methods; pinned to Charm++ 8.0.0. |
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
