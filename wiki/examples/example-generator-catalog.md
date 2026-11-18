# Example Generator Catalog

> Inventory of runnable `nrpy/examples/*.py` generators, CLI flags, configured CI coverage, companion scripts, and detail leaves. · Status: confirmed
> Up: [Examples](index.md)

## Summary

Every runnable generator module under `nrpy/examples/*.py` (`__init__.py` is package import only) is listed below with its CLI flags, configured CI coverage (`manual` denotes source-supported scripts outside GitHub Actions), and defining detail leaf, followed by checked-in companion directories. Generated `project/**` trees, plots, ZIP archives, and compiled binaries remain disposable outputs ([Generated Output Boundaries](../architecture/generated-output-boundaries.md)).

## Detail

### Top-Level `nrpy/examples/*.py` Generators

Run `python -m nrpy.examples.<module> [flags]` from the repository root ([Build And Run](../architecture/build-and-run.md)); most generators recreate `project/<project_name>/` cleanly on each run.

| Generator | Key CLI flags | CI check | Detail leaf |
| --- | --- | --- | --- |
| `bhahaha.py` | `[--fdorder N] [--outrootdir DIR] [--cpp] [--no-openmp]` | Ubuntu/macOS C static-library `make` build only | [Apparent Horizon Library](apparent-horizon-library.md) |
| `blackhole_spectroscopy.py` | `[--cuda] [--fccz4] [--raytracing-time]` | Ubuntu/macOS default BSSN/OpenMP `make` build only | [Standalone GR/BHaH](standalone-gr-bhah.md) |
| `carpet_baikal_thorns.py` | none | `einsteintoolkit-validation` ET build + `Baikal`/`BaikalVacuum` testsuites | [Einstein Toolkit Thorn Generators](einstein-toolkit-thorn-generators.md) |
| `carpet_wavetoy_thorns.py` | none | `einsteintoolkit-validation` ET build + `WaveToyNRPy` testsuite | [Einstein Toolkit Thorn Generators](einstein-toolkit-thorn-generators.md) |
| `carpetx_baikal_thorns.py` | none | `full_nrpy_local_ci.sh` generation only; no CI build/run | [Einstein Toolkit Thorn Generators](einstein-toolkit-thorn-generators.md) |
| `carpetx_wavetoy_thorns.py` | none | `full_nrpy_local_ci.sh` generation only; no CI build/run | [Einstein Toolkit Thorn Generators](einstein-toolkit-thorn-generators.md) |
| `dendro_bssn.py` | `[--project-dir D] [--fd-order {4,6,8}] [--conformal-factor {W,chi}] [--ybs-{gamma,momentum}]` | `dendro-validation`, `dendrolib-canary` | [Dendro BSSN](../infrastructures/dendro/bssn-application-wiring.md) |
| `dendro_fccz4.py` | `[--project-dir D] [--fd-order {4,6,8}] [--conformal-factor {W,chi}] [--ybs-{gamma,momentum}]` | `dendro-validation`, `dendrolib-canary` | [Dendro fCCZ4](../infrastructures/dendro/fccz4-application-wiring.md) |
| `groovy_TOV_BSSN.py` | none (clones/builds GRHayL) | manual only; not invoked in CI | [Matter TOV Workflows](matter-tov-workflows.md) |
| `hydro_without_hydro.py` | `[--cuda] [--floating_point_precision TYPE]` | Ubuntu/macOS OpenMP `make` build (`--cuda` in local script only) | [Matter TOV Workflows](matter-tov-workflows.md) |
| `kasner_exact_evolution.py` | `[--floating_point_precision TYPE]` (`--cuda` rejected by `GeneralRFM`) | manual OpenMP route; not invoked in CI | [Standalone GR/BHaH](standalone-gr-bhah.md) |
| `manga_bhah_lib.py` | none (`create_lib=True`) | commented out in `main.yml` | [Matter TOV Workflows](matter-tov-workflows.md) |
| `massive_single_geodesic_integrator_analytical.py` | none (requires GSL) | manual single-ray route | [Geodesic Raytracing](geodesic-raytracing.md) |
| `nrpyelliptic_conformally_flat.py` | `[--cuda] [--floating_point_precision TYPE]` | Ubuntu/macOS OpenMP `make` build (`--cuda` in local script only) | [Elliptic Initial Data](elliptic-initial-data.md) |
| `nrpypn_quasicircular_momenta.py` | none | Ubuntu/macOS `make` build only | [Waveform JAX PN Generators](waveform-jax-pn-generators.md) |

### Photon, Waveform, And Other Evolution Generators

| Generator | Key CLI flags | CI check | Detail leaf |
| --- | --- | --- | --- |
| `photon_batch_geodesic_integrator_analytical.py` | observer geometry, `--scan-density`, optional `--cuda` | manual tiled batch route | [Geodesic Raytracing](geodesic-raytracing.md) |
| `photon_batch_geodesic_integrator_numerical.py` | spacetime `.bin`, `--scan-density`, `--eom` | manual tiled batch | [Geodesic Raytracing](geodesic-raytracing.md) |
| `photon_single_geodesic_integrator_analytical.py` | `--spacetime`, `--eom`, observer geometry | manual single-ray route | [Geodesic Raytracing](geodesic-raytracing.md) |
| `photon_single_geodesic_integrator_numerical.py` | spacetime `.bin`, observer geometry, `--eom` | manual single-ray | [Geodesic Raytracing](geodesic-raytracing.md) |
| `sebobv1_jax.py` | none | Ubuntu/macOS `pip install .` + 1 scalar call ([`sebobv1_jax` step](../validation/generated-project-ci.md#sebobv1_jax-install-and-call-step)) | [Waveform JAX PN Generators](waveform-jax-pn-generators.md) |
| `sebobv2.py` | `[--projected-attachment-diagnostics]` (GSL; `use_projected_attachment` in `.par`) | Ubuntu/macOS build + `sebobv2` & projected-attachment CI | [Waveform JAX PN Generators](waveform-jax-pn-generators.md) |
| `seobnrv5_aligned_spin_inspiral.py` | `[-seobnrv5_{bob,nrnqc_bob,nrpy}] [-{calibration_no_spin,calibration_spin,nrpy_calibrated}]` | Ubuntu/macOS; `sebob-consistency-test` | [Waveform JAX PN](waveform-jax-pn-generators.md) |
| `spinning_blackhole.py` | `[--cuda] [--floating_point_precision TYPE]` | Ubuntu/macOS OpenMP `make` build (`--cuda` in local script only) | [Standalone GR/BHaH](standalone-gr-bhah.md) |
| `superB_blackhole_spectroscopy.py` | `[--paper]` | `charmpp-validation` Charm++ build only | [superB Charm++ Workflows](superb-charm-workflows.md) |
| `superB_nrpyelliptic_conformally_flat.py` | `[--floating_point_precision TYPE]` | `charmpp-validation` Charm++ build only | [superB Charm++ Workflows](superb-charm-workflows.md) |
| `superB_two_blackholes_collide.py` | none (links BHaHAHA) | `charmpp-validation` build + `charmrun +p2` exit check | [superB Charm++ Workflows](superb-charm-workflows.md) |
| `tovola_neutron_star.py` | none (ignores `--cuda` in local script) | Ubuntu/macOS C `make` build only | [Matter TOV Workflows](matter-tov-workflows.md) |
| `two_blackholes_collide.py` | `[--cuda] [--floating_point_precision TYPE] [--raytracing-time]` | Ubuntu/macOS OpenMP `make` build (`--cuda` in local script only) | [Standalone GR/BHaH](standalone-gr-bhah.md) |
| `wave_equation_cartesian.py` | none | Ubuntu/macOS `make` build only | [Wave Equation Generators](wave-equation-generators.md) |
| `wave_equation_curvilinear.py` | `[--cuda] [--floating_point_precision P] [--disable_intrinsics] [--disable_rfm_precompute]` | Ubuntu/macOS OpenMP `make` build | [Wave Equation Generators](wave-equation-generators.md) |
| `wave_equation_multicoordinates.py` | `[--cuda] [--floating_point_precision P] [--disable_intrinsics] [--disable_rfm_precompute]` | Ubuntu/macOS OpenMP `make` build | [Wave Equation Generators](wave-equation-generators.md) |

Evidence (CI behavior): `.github/workflows/main.yml`; `.github/workflows/dendrolib-canary.yml`; `.github/workflows/sebobv2-projected-attachment.yml`; `.github/full_nrpy_local_ci.sh`; `dendro_bssn.py`; `dendro_fccz4.py`; `sebobv2.py`.

### Companion Script And Test Directories

| Companion group | Checked-in files | Generated-output boundary | Detail leaf |
| --- | --- | --- | --- |
| `nrpy/examples/et_WaveToyfiles/**` | `ThornList`, `WaveToyNRPy*.par`, `test.ccl`, `uuGF.x.asc` | Checked-in ET inputs; `project/**` disposable | [ET Thorn Generators](einstein-toolkit-thorn-generators.md) |
| `nrpy/examples/geodesic_visualizations/*.py` | Trajectory plots, lensed-image rendering, blueprint analysis | PNGs, `trajectory.txt`, and blueprint binaries are generated | [Geodesic Raytracing](geodesic-raytracing.md) |
| `nrpy/examples/tests/*.py` | `sebob*_consistency_check.py`, `sebobv2_projected_attachment_check.py`, `dendro_application_check*.py` | Temp worktrees disposable | [Generated Project CI](../validation/generated-project-ci.md) |
| `nrpy/examples/superb_paper2/*.py` | `superB_blackhole_spectroscopy_{8Mseparation,last_orbit}.py` | Generated Charm++ `project/**` disposable | [superB Charm++ Workflows](superb-charm-workflows.md) |
| `nrpy/examples/q1.par.lowres.toml` | Checked-in low-res Dendro BSSN/fCCZ4 TOML parameter file | Copied into generated `pars/` alongside `<stem>.toml` | [Dendro Parameters](../infrastructures/dendro/runtime-parameters.md) |

See [SEBOBv1 JAX Workflow](../infrastructures/jax/sebobv1-jax-workflow.md) for `SEOBNRv5_aligned_spin_coefficients` and `Commondata` fields.

## Sources

- [main.yml](../../.github/workflows/main.yml) - CI build and validation matrix
- [dendrolib-canary.yml](../../.github/workflows/dendrolib-canary.yml) - `dendro-validation-dendrolib-master`
- [sebobv2-projected-attachment.yml](../../.github/workflows/sebobv2-projected-attachment.yml) - `projected-reference-and-matching`
- [full_nrpy_local_ci.sh](../../.github/full_nrpy_local_ci.sh) - manual local build and `--cuda` configuration script
- [bhahaha.py](../../nrpy/examples/bhahaha.py) - `parser`, `project_name`
- [blackhole_spectroscopy.py](../../nrpy/examples/blackhole_spectroscopy.py) - `parser`, `project_name`
- [carpet_baikal_thorns.py](../../nrpy/examples/carpet_baikal_thorns.py) - `project_name`, `thorn_names`
- [carpet_wavetoy_thorns.py](../../nrpy/examples/carpet_wavetoy_thorns.py) - `project_name`, thorn names
- [carpetx_baikal_thorns.py](../../nrpy/examples/carpetx_baikal_thorns.py) - `project_name`, `thorn_names`
- [carpetx_wavetoy_thorns.py](../../nrpy/examples/carpetx_wavetoy_thorns.py) - `project_name`, thorn names
- [dendro_bssn.py](../../nrpy/examples/dendro_bssn.py) - `parse_args`, `SOLVER_NAME`
- [dendro_fccz4.py](../../nrpy/examples/dendro_fccz4.py) - `parse_args`, `SOLVER_NAME`
- [groovy_TOV_BSSN.py](../../nrpy/examples/groovy_TOV_BSSN.py) - `project_name`, `repo_url`
- [hydro_without_hydro.py](../../nrpy/examples/hydro_without_hydro.py) - `parser`, `project_name`
- [kasner_exact_evolution.py](../../nrpy/examples/kasner_exact_evolution.py) - `parser`, `project_name`
- [manga_bhah_lib.py](../../nrpy/examples/manga_bhah_lib.py) - `project_name`, `create_lib=True`
- [massive_single_geodesic_integrator_analytical.py](../../nrpy/examples/massive_single_geodesic_integrator_analytical.py) - `project_name`
- [nrpyelliptic_conformally_flat.py](../../nrpy/examples/nrpyelliptic_conformally_flat.py) - `parser`, `project_name`
- [nrpypn_quasicircular_momenta.py](../../nrpy/examples/nrpypn_quasicircular_momenta.py) - `project_name`
- [photon_batch_geodesic_integrator_analytical.py](../../nrpy/examples/photon_batch_geodesic_integrator_analytical.py) - `--outdir`, `--cuda`
- [photon_single_geodesic_integrator_analytical.py](../../nrpy/examples/photon_single_geodesic_integrator_analytical.py) - `project_name`
- [sebobv1_jax.py](../../nrpy/examples/sebobv1_jax.py) - `project_name`
- [sebobv2.py](../../nrpy/examples/sebobv2.py) - `--projected-attachment-diagnostics`, `project_name`
- [seobnrv5_aligned_spin_inspiral.py](../../nrpy/examples/seobnrv5_aligned_spin_inspiral.py) - `argparse`, `project_name`
- [photon_single_geodesic_integrator_numerical.py](../../nrpy/examples/photon_single_geodesic_integrator_numerical.py) - numerical photon generation options
- [photon_batch_geodesic_integrator_numerical.py](../../nrpy/examples/photon_batch_geodesic_integrator_numerical.py) - numerical photon generation options
- [spinning_blackhole.py](../../nrpy/examples/spinning_blackhole.py) - `parser`, `project_name`
- [superB_blackhole_spectroscopy.py](../../nrpy/examples/superB_blackhole_spectroscopy.py) - `parser`, `project_name`
- [superB_nrpyelliptic_conformally_flat.py](../../nrpy/examples/superB_nrpyelliptic_conformally_flat.py) - `parser`, `project_name`
- [superB_two_blackholes_collide.py](../../nrpy/examples/superB_two_blackholes_collide.py) - `project_name`
- [tovola_neutron_star.py](../../nrpy/examples/tovola_neutron_star.py) - `project_name`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `parser`, `project_name`
- [wave_equation_cartesian.py](../../nrpy/examples/wave_equation_cartesian.py) - `project_name`
- [wave_equation_curvilinear.py](../../nrpy/examples/wave_equation_curvilinear.py) - `parser`, `project_name`
- [wave_equation_multicoordinates.py](../../nrpy/examples/wave_equation_multicoordinates.py) - `parser`, `project_name`
- [q1.par.lowres.toml](../../nrpy/examples/q1.par.lowres.toml) - checked-in low-resolution Dendro BBH TOML
- [ThornList](../../nrpy/examples/et_WaveToyfiles/ThornList) - ET ThornList
- [WaveToyNRPy.par](../../nrpy/examples/et_WaveToyfiles/WaveToyNRPy.par) - WaveToy parfile
- [test.ccl](../../nrpy/examples/et_WaveToyfiles/test/test.ccl) - `TEST WaveToyNRPy_test`
- [WaveToyNRPy_test.par](../../nrpy/examples/et_WaveToyfiles/test/WaveToyNRPy_test.par) - testsuite parfile
- [uuGF.x.asc](../../nrpy/examples/et_WaveToyfiles/test/WaveToyNRPy_test/uuGF.x.asc) - CarpetIOASCII reference output
- [visualize_trajectory.py](../../nrpy/examples/geodesic_visualizations/visualize_trajectory.py) - `visualize_trajectory`
- [blueprint_config_and_schema.py](../../nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py) - `BLUEPRINT_DTYPE`
- [render_lensed_image.py](../../nrpy/examples/geodesic_visualizations/render_lensed_image.py) - `generate_static_lensed_image`
- [visualize_lensed_image.py](../../nrpy/examples/geodesic_visualizations/visualize_lensed_image.py) - `main`
- [blueprint_analysis.py](../../nrpy/examples/geodesic_visualizations/blueprint_analysis.py) - `diagnose_blueprint`
- [superB_blackhole_spectroscopy_8Mseparation.py](../../nrpy/examples/superb_paper2/superB_blackhole_spectroscopy_8Mseparation.py) - 8M separation superB spectroscopy setup
- [superB_blackhole_spectroscopy_last_orbit.py](../../nrpy/examples/superb_paper2/superB_blackhole_spectroscopy_last_orbit.py) - last-orbit superB spectroscopy setup
- [sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `__main__`
- [sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `__main__`
- [sebobv2_projected_attachment_check.py](../../nrpy/examples/tests/sebobv2_projected_attachment_check.py) - `main`
- [dendro_application_check.py](../../nrpy/examples/tests/dendro_application_check.py) - `Leg.run`
- [dendro_application_check_reference.py](../../nrpy/examples/tests/dendro_application_check_reference.py) - `trusted_dict`

## See Also

- Depends on: [Generated Output Boundaries](../architecture/generated-output-boundaries.md)
- Depends on: [Build And Run](../architecture/build-and-run.md)
- Validated by: [Generated Project CI](../validation/generated-project-ci.md)
