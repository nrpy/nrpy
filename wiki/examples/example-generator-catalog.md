# Example Generator Catalog

> Inventory the runnable example generators, companion scripts, output families, prerequisites, validation status, and owning pages. · Status: confirmed · Last reconciled: 06-30-2026
> Up: [Examples](index.md)

## Summary

This catalog is the inventory leaf for `nrpy/examples`. Each non-empty
top-level generator under `nrpy/examples/*.py` appears once below with its
command shape, output family, prerequisite class, validation route or manual
status, and owning detail page. Companion files under
`nrpy/examples/geodesic_visualizations/`, `nrpy/examples/tests/`, and
`nrpy/examples/et_WaveToyfiles/` are cataloged by group after the generator
inventory.

Generated `project/**` trees, rendered plots, images, ZIPs, binaries, and
scratch logs remain artifacts. Use the checked-in generator and companion files
as source evidence, then follow the owning page for details.

## Detail

Run module commands from the repository root after installing NRPy or setting
`PYTHONPATH` to include `.` as described in [Build And Run](../architecture/build-and-run.md).

| Generator | Command shape | Output family | Prerequisites | Validation route or status | Owning page |
| --- | --- | --- | --- | --- | --- |
| `bhahaha.py` | `python -m nrpy.examples.bhahaha [--fdorder N] [--outrootdir DIR] [--cpp] [--no-openmp]` | Static BHaHAHA apparent-horizon library under the chosen output root | Python, C compiler, `make`; OpenMP optional | Ubuntu/macOS CI generate and build default library path; horizon-enabled examples exercise linked use | [Apparent Horizon Library](apparent-horizon-library.md) |
| `blackhole_spectroscopy.py` | `python -m nrpy.examples.blackhole_spectroscopy [--cuda] [--floating_point_precision TYPE]` | Standalone BHaH binary-black-hole spectroscopy project | Python, C or CUDA toolchain, `make`, GSL | Generation/build status from standalone examples and CI context; runtime is heavier/manual | [Standalone GR/BHaH](standalone-gr-bhah.md) |
| `carpet_baikal_thorns.py` | `python -m nrpy.examples.carpet_baikal_thorns` | ETLegacy/Carpet Baikal and BaikalVacuum thorns | Python for generation; Einstein Toolkit checkout for build/test | GitHub Actions ET validation builds and runs Baikal testsuites | [Einstein Toolkit Thorn Generators](einstein-toolkit-thorn-generators.md) |
| `carpet_wavetoy_thorns.py` | `python -m nrpy.examples.carpet_wavetoy_thorns` | ETLegacy/Carpet WaveToyNRPy, IDWaveToyNRPy, and diagWaveToyNRPy thorns | Python for generation; Einstein Toolkit checkout for build/test | GitHub Actions ET validation builds and runs WaveToyNRPy testsuite | [Einstein Toolkit Thorn Generators](einstein-toolkit-thorn-generators.md) |
| `carpetx_baikal_thorns.py` | `python -m nrpy.examples.carpetx_baikal_thorns` | CarpetX BaikalX and BaikalVacuumX thorns | Python for generation; CarpetX/Einstein Toolkit environment for build/test | Local full-CI generation route; ET compile skipped in local script | [Einstein Toolkit Thorn Generators](einstein-toolkit-thorn-generators.md) |
| `carpetx_wavetoy_thorns.py` | `python -m nrpy.examples.carpetx_wavetoy_thorns` | CarpetX WaveToyNRPyX, IDWaveToyNRPyX, and diagWaveToyNRPyX thorns | Python for generation; CarpetX/Einstein Toolkit environment for build/test | Local full-CI generation route; ET compile skipped in local script | [Einstein Toolkit Thorn Generators](einstein-toolkit-thorn-generators.md) |
| `groovy_TOV_BSSN.py` | `python -m nrpy.examples.groovy_TOV_BSSN` | Standalone BHaH/GRoovy TOV GRHD evolution project | Python, C compiler, `make`, GSL, GRHayL clone/build | Manual source-supported route; inspected CI does not build it | [Matter TOV Workflows](matter-tov-workflows.md) |
| `hydro_without_hydro.py` | `python -m nrpy.examples.hydro_without_hydro [--cuda] [--floating_point_precision TYPE]` | Standalone BHaH static-fluid spacetime evolution project | Python, C or CUDA toolchain, `make`, GSL | Ubuntu/macOS CI generate and build default OpenMP path | [Matter TOV Workflows](matter-tov-workflows.md) |
| `kasner_exact_evolution.py` | `python -m nrpy.examples.kasner_exact_evolution [--cuda] [--floating_point_precision TYPE]` | Standalone BHaH Kasner benchmark project | Python, C or CUDA toolchain, `make` | Manual/source-supported benchmark route | [Standalone GR/BHaH](standalone-gr-bhah.md) |
| `manga_bhah_lib.py` | `python -m nrpy.examples.manga_bhah_lib` | MANGA-facing `bhah_lib` library project | Python, C compiler, `make`, GSL | Source-supported library route; CI commands are present but commented out | [Matter TOV Workflows](matter-tov-workflows.md) |
| `mass_geodesic_integrator.py` | `python -m nrpy.examples.mass_geodesic_integrator` | Standalone massive-particle geodesic C project plus trajectory visualization copy | Python, C compiler, `make`, GSL, NumPy/Matplotlib for visualization | Manual/source-supported single-ray route | [Geodesic Raytracing](geodesic-raytracing.md) |
| `nrpyelliptic_conformally_flat.py` | `python -m nrpy.examples.nrpyelliptic_conformally_flat [--cuda] [--floating_point_precision TYPE]` | Standalone BHaH NRPyElliptic conformally flat project | Python, C or CUDA toolchain, `make` | Ubuntu/macOS CI generate and build default path | [Elliptic Initial Data](elliptic-initial-data.md) |
| `nrpypn_quasicircular_momenta.py` | `python -m nrpy.examples.nrpypn_quasicircular_momenta` | Standalone BHaH PN momentum utility project | Python, C compiler, `make` | Ubuntu/macOS CI generate and build | [Waveform JAX PN Generators](waveform-jax-pn-generators.md) |
| `photon_geodesic_batch_integrator.py` | `python -m nrpy.examples.photon_geodesic_batch_integrator [--cuda] [--outdir DIR]` | Standalone tiled photon raytracing project emitting light-blueprint ZIP artifacts | Python, C or CUDA toolchain, `make`, NumPy/Matplotlib/Numba/Pillow for visualization | Manual/source-supported batch route | [Geodesic Raytracing](geodesic-raytracing.md) |
| `photon_geodesic_integrator.py` | `python -m nrpy.examples.photon_geodesic_integrator` | Standalone single-photon geodesic project plus trajectory visualization copy | Python, C compiler, `make`, NumPy/Matplotlib for visualization | Manual/source-supported single-ray route | [Geodesic Raytracing](geodesic-raytracing.md) |
| `sebobv1_jax.py` | `python -m nrpy.examples.sebobv1_jax` | Python/JAX package project for SEOBNRv5 coefficient initialization | Python and JAX stack | Ubuntu/macOS CI run generation; no generated-package build step in inspected workflow | [Waveform JAX PN Generators](waveform-jax-pn-generators.md) |
| `sebobv2.py` | `python -m nrpy.examples.sebobv2` | GSL-backed BHaH C waveform project | Python, C compiler, `make`, GSL | Example consistency helper covers current/trusted executable comparisons | [Waveform JAX PN Generators](waveform-jax-pn-generators.md) |
| `seobnrv5_aligned_spin_inspiral.py` | `python -m nrpy.examples.seobnrv5_aligned_spin_inspiral [-seobnrv5_bob|-seobnrv5_nrnqc_bob|-seobnrv5_nrpy] [-calibration_no_spin|-calibration_spin]` | GSL-backed BHaH C SEOBNRv5 waveform project family | Python, C compiler, `make`, GSL | Example consistency helper covers current/trusted executable comparisons | [Waveform JAX PN Generators](waveform-jax-pn-generators.md) |
| `spinning_blackhole.py` | `python -m nrpy.examples.spinning_blackhole [--cuda] [--floating_point_precision TYPE]` | Standalone BHaH spinning black-hole project | Python, C or CUDA toolchain, `make`; BHaHAHA in supported OpenMP/double path | Ubuntu/macOS CI generate and build default OpenMP path; local full-CI also exercises CUDA generation/build | [Standalone GR/BHaH](standalone-gr-bhah.md) |
| `superB_blackhole_spectroscopy.py` | `python -m nrpy.examples.superB_blackhole_spectroscopy [--paper]` | Charm++/superB black-hole spectroscopy project | Python, Charm++ toolchain, `make`, GSL | Charm++ CI generates and builds; executable not run in inspected workflow | [superB Charm++ Workflows](superb-charm-workflows.md) |
| `superB_nrpyelliptic_conformally_flat.py` | `python -m nrpy.examples.superB_nrpyelliptic_conformally_flat [--floating_point_precision TYPE]` | Charm++/superB NRPyElliptic conformally flat project | Python, Charm++ toolchain, `make` | Charm++ CI generates and builds | [superB Charm++ Workflows](superb-charm-workflows.md) |
| `superB_two_blackholes_collide.py` | `python -m nrpy.examples.superB_two_blackholes_collide` | Charm++/superB Brill-Lindquist binary black-hole project | Python, Charm++ toolchain, `make`; generated BHaHAHA library | Charm++ CI generates, builds, and runs with `+p2` | [superB Charm++ Workflows](superb-charm-workflows.md) |
| `tovola_neutron_star.py` | `python -m nrpy.examples.tovola_neutron_star` | Standalone BHaH TOVola initial-data and constraint project | Python, C compiler, `make`, GSL | Ubuntu/macOS CI generate and build | [Matter TOV Workflows](matter-tov-workflows.md) |
| `two_blackholes_collide.py` | `python -m nrpy.examples.two_blackholes_collide [--cuda] [--floating_point_precision TYPE] [--raytracing-outputs]` | Standalone BHaH Brill-Lindquist binary black-hole project | Python, C or CUDA toolchain, `make`; BHaHAHA in supported OpenMP/double path | Ubuntu/macOS CI generate and build default OpenMP path; local full-CI also exercises CUDA generation/build; raytracing output is OpenMP-only | [Standalone GR/BHaH](standalone-gr-bhah.md) |
| `wave_equation_cartesian.py` | `python -m nrpy.examples.wave_equation_cartesian` | Minimal standalone BHaH Cartesian wave project | Python, C compiler, `make` | First-run quickstart route; Ubuntu/macOS CI generate and build | [Wave Equation Generators](wave-equation-generators.md) |
| `wave_equation_curvilinear.py` | `python -m nrpy.examples.wave_equation_curvilinear [--cuda] [--floating_point_precision TYPE] [--disable_intrinsics] [--disable_rfm_precompute]` | Standalone BHaH single-coordinate curvilinear wave project | Python, C or CUDA toolchain, `make` | Ubuntu/macOS CI generate and build default OpenMP path; local full-CI also exercises CUDA generation/build | [Wave Equation Generators](wave-equation-generators.md) |
| `wave_equation_multicoordinates.py` | `python -m nrpy.examples.wave_equation_multicoordinates [--cuda] [--floating_point_precision TYPE] [--disable_intrinsics] [--disable_rfm_precompute]` | Standalone BHaH multicoordinate wave project | Python, C or CUDA toolchain, `make` | Ubuntu/macOS CI generate and build default OpenMP path; local full-CI also exercises CUDA generation/build | [Wave Equation Generators](wave-equation-generators.md) |

Companion groups:

| Companion group | Checked-in source shape | Artifact boundary | Owning page |
| --- | --- | --- | --- |
| `nrpy/examples/et_WaveToyfiles/**` | ET ThornList, Carpet WaveToy parfiles, testsuite descriptor, and checked-in ASCII fixture | Generated thorn files under `project/**` stay artifacts; checked-in fixtures are source evidence | [Einstein Toolkit Thorn Generators](einstein-toolkit-thorn-generators.md) |
| `nrpy/examples/geodesic_visualizations/*.py` | Trajectory plotter, light-blueprint schema, lensed-image renderer, image visualizer, and blueprint diagnostics | `trajectory.txt`, PNGs, downloaded textures, and `light_blueprint_*.zip` files stay artifacts | [Geodesic Raytracing](geodesic-raytracing.md) |
| `nrpy/examples/tests/sebob*_consistency_check.py` | Current-vs-trusted waveform consistency scripts | Trusted/current generated executable directories and waveform stdout are run artifacts | [Waveform JAX PN Generators](waveform-jax-pn-generators.md) |

The direct source rows in [Sources](../../raw/SOURCES.md) register the cited
example files. The aggregate row `example-generators-and-companion-scripts`
remains partial because the examples tree is broader than this catalog and
aggregate ingestion is not proven complete for all future files.

## Sources

- [wave_equation_cartesian.py](../../nrpy/examples/wave_equation_cartesian.py) - `project_name`
- [wave_equation_curvilinear.py](../../nrpy/examples/wave_equation_curvilinear.py) - `parser`, `project_name`
- [wave_equation_multicoordinates.py](../../nrpy/examples/wave_equation_multicoordinates.py) - `parser`, `project_name`
- [two_blackholes_collide.py](../../nrpy/examples/two_blackholes_collide.py) - `parser`, `project_name`, `--raytracing-outputs`
- [blackhole_spectroscopy.py](../../nrpy/examples/blackhole_spectroscopy.py) - `parser`, `project_name`
- [spinning_blackhole.py](../../nrpy/examples/spinning_blackhole.py) - `parser`, `project_name`
- [kasner_exact_evolution.py](../../nrpy/examples/kasner_exact_evolution.py) - `parser`, `project_name`
- [nrpyelliptic_conformally_flat.py](../../nrpy/examples/nrpyelliptic_conformally_flat.py) - `parser`, `project_name`
- [superB_nrpyelliptic_conformally_flat.py](../../nrpy/examples/superB_nrpyelliptic_conformally_flat.py) - `parser`, `project_name`
- [carpet_wavetoy_thorns.py](../../nrpy/examples/carpet_wavetoy_thorns.py) - `project_name`, `ID_thorn_name`, `evol_thorn_name`, `diag_thorn_name`
- [carpet_baikal_thorns.py](../../nrpy/examples/carpet_baikal_thorns.py) - `project_name`, `thorn_names`
- [carpetx_wavetoy_thorns.py](../../nrpy/examples/carpetx_wavetoy_thorns.py) - `project_name`, `ID_thorn_name`, `evol_thorn_name`, `diag_thorn_name`
- [carpetx_baikal_thorns.py](../../nrpy/examples/carpetx_baikal_thorns.py) - `project_name`, `thorn_names`
- [superB_two_blackholes_collide.py](../../nrpy/examples/superB_two_blackholes_collide.py) - `project_name`
- [superB_blackhole_spectroscopy.py](../../nrpy/examples/superB_blackhole_spectroscopy.py) - `parser`, `project_name`
- [seobnrv5_aligned_spin_inspiral.py](../../nrpy/examples/seobnrv5_aligned_spin_inspiral.py) - `argparse`, `project_name`
- [sebobv2.py](../../nrpy/examples/sebobv2.py) - `project_name`
- [sebobv1_jax.py](../../nrpy/examples/sebobv1_jax.py) - `project_name`, `Infrastructure`
- [nrpypn_quasicircular_momenta.py](../../nrpy/examples/nrpypn_quasicircular_momenta.py) - `project_name`
- [mass_geodesic_integrator.py](../../nrpy/examples/mass_geodesic_integrator.py) - `project_name`, `gsl-config`
- [photon_geodesic_integrator.py](../../nrpy/examples/photon_geodesic_integrator.py) - `project_name`
- [photon_geodesic_batch_integrator.py](../../nrpy/examples/photon_geodesic_batch_integrator.py) - `--outdir`, `--cuda`, `project_name`
- [tovola_neutron_star.py](../../nrpy/examples/tovola_neutron_star.py) - `project_name`
- [hydro_without_hydro.py](../../nrpy/examples/hydro_without_hydro.py) - `parser`, `project_name`
- [groovy_TOV_BSSN.py](../../nrpy/examples/groovy_TOV_BSSN.py) - `project_name`, `repo_url`
- [manga_bhah_lib.py](../../nrpy/examples/manga_bhah_lib.py) - `project_name`, `create_lib=True`
- [bhahaha.py](../../nrpy/examples/bhahaha.py) - `parser`, `project_name`
- [ThornList](../../nrpy/examples/et_WaveToyfiles/ThornList) - Einstein Toolkit component-list fixture
- [WaveToyNRPy.par](../../nrpy/examples/et_WaveToyfiles/WaveToyNRPy.par) - Carpet WaveToy parfile fixture
- [test.ccl](../../nrpy/examples/et_WaveToyfiles/test/test.ccl) - `TEST WaveToyNRPy_test`
- [WaveToyNRPy_test.par](../../nrpy/examples/et_WaveToyfiles/test/WaveToyNRPy_test.par) - Carpet WaveToy test parfile
- [uuGF.x.asc](../../nrpy/examples/et_WaveToyfiles/test/WaveToyNRPy_test/uuGF.x.asc) - checked-in CarpetIOASCII fixture
- [visualize_trajectory.py](../../nrpy/examples/geodesic_visualizations/visualize_trajectory.py) - `visualize_trajectory`
- [blueprint_config_and_schema.py](../../nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py) - `BLUEPRINT_DTYPE`
- [render_lensed_image.py](../../nrpy/examples/geodesic_visualizations/render_lensed_image.py) - `generate_static_lensed_image`
- [visualize_lensed_image.py](../../nrpy/examples/geodesic_visualizations/visualize_lensed_image.py) - `main`
- [blueprint_analysis.py](../../nrpy/examples/geodesic_visualizations/blueprint_analysis.py) - `diagnose_blueprint`
- [sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `__main__`
- [sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `__main__`
- [Generated Output Boundaries](../architecture/generated-output-boundaries.md) - generated-project source boundary
- [Generated Project CI](../validation/generated-project-ci.md) - compiled CI validation route summary

## See Also

- Parent: [Examples](index.md)
- Depends on: [Generated Output Boundaries](../architecture/generated-output-boundaries.md)
- Validated by: [Generated Project CI](../validation/generated-project-ci.md)
- Example: [First Wave Equation Run](first-wave-equation-run.md)
- Example: [Black Hole Evolution](black-hole-evolution.md)
