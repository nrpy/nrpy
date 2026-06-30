# Glossary

> Canonical spellings for recurring code, domain, and experiment entities. · Status: confirmed · Last reconciled: 2026-06-29

## Summary

Use these spellings when routing, citing, or adding leaves. Terms link to owning
pages where a first-pass owner exists.

## Terms

| Term | Meaning |
| --- | --- |
| NRPy | Python/SymPy-based symbolic code generation toolkit for numerical relativity and relativistic astrophysics; see [Architecture Overview](architecture/overview.md). |
| SymPy | Python symbolic algebra package used to build expressions before C generation; see [Symbolic Codegen Lifecycle](architecture/symbolic-codegen-lifecycle.md). |
| generated C | C emitted by NRPy generators and infrastructure rather than handwritten source; see [Generated Output Boundaries](architecture/generated-output-boundaries.md). |
| BHaH | NRPy's standalone generated-application infrastructure for BH@H-style workflows; see [BHaH Lifecycle](infrastructures/bhah-lifecycle.md). |
| BH@H | BlackHoles@Home, the single-patch numerical relativity workflow represented by current BHaH examples. |
| BHaHAHA | Apparent-horizon finding library integrated by BHaH black-hole workflows; see [BHaHAHA Horizon Runtime](infrastructures/bhah/bhahaha-horizon-runtime.md). |
| CFunction | Python representation of a generated C function; see [C Function Registry](core/c-function-registry.md). |
| CFunction_dict | Global registry mapping generated function names to `CFunction` objects; see [C Function Registry](core/c-function-registry.md). |
| register_CFunction | API that registers generated C functions; see [C Function Registry](core/c-function-registry.md). |
| c_codegen | Main SymPy-to-C emission function; see [C Codegen](core/c-codegen.md). |
| CCodeGen | Codegen class that backs lower-level C code emission; see [C Codegen](core/c-codegen.md). |
| CSE | Common-subexpression elimination; NRPy helper pages use this spelling for SymPy CSE preprocessing, postprocessing, and deterministic temporary ordering; see [CSE And Printer Support](core/helpers/cse-and-printer-support.md). |
| SIMD | Single instruction, multiple data; NRPy helper pages use this spelling for symbolic `*SIMD` rewrites and handwritten SIMD intrinsic macro headers; see [SIMD And Intrinsic Support](core/helpers/simd-and-intrinsic-support.md). |
| CUDA | NVIDIA GPU programming target used by NRPy helper and infrastructure code for generated kernels, launch calls, and device-oriented macro headers; see [Loop Kernel And Device Helpers](core/helpers/loop-kernel-and-device-helpers.md). |
| intrinsics | Low-level operation spellings exposed through NRPy's SIMD and CUDA helper headers, not generated project copies; see [SIMD And Intrinsic Support](core/helpers/simd-and-intrinsic-support.md). |
| clang-format | External formatter invoked by NRPy helper code when formatting generated C-like text; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| NRPyParameter | Runtime-independent NRPy parameter object; see [Gridfunctions And Parameters](core/gridfunctions-and-parameters.md). |
| CodeParameter | Code parameter emitted into generated projects; see [Gridfunctions And Parameters](core/gridfunctions-and-parameters.md). |
| GridFunction | Registered evolved, auxiliary, or diagnostic grid quantity; see [Gridfunctions And Parameters](core/gridfunctions-and-parameters.md). |
| commondata | Generated data that is common across grids, usually represented through `commondata_struct`; see [Gridfunctions And Parameters](core/gridfunctions-and-parameters.md). |
| reference metric | Coordinate-system-aware metric context used by symbolic equations and infrastructures; see [Reference Metrics](core/reference-metrics.md). |
| indexed expression | Tensor or vector structure built with `nrpy.indexedexp`; see [Indexed Expressions](core/indexed-expressions.md). |
| ADM | Arnowitt-Deser-Misner 3+1 quantities used by GR conversions and initial-data pages; see [Metric Conversions And Matter](equations/general-relativity/metric-conversions-and-matter.md). |
| BSSN | Baumgarte-Shapiro-Shibata-Nakamura formulation family implemented by NRPy equation modules; see [BSSN Family](equations/general-relativity/bssn-family.md). |
| GRHD | General relativistic hydrodynamics equation family for conserved variables, fluxes, sources, speeds, and HLL helpers; see [GRHD](equations/grhd.md). |
| TOV | Tolman-Oppenheimer-Volkoff stellar-equilibrium ODE system; see [TOV Equations](equations/tov-equations.md). |
| Psi4 | Weyl radiation scalar implementation and tetrad construction path; see [Psi4 And Tetrads](equations/general-relativity/psi4-and-tetrads.md). |
| SEOBNR | Effective-One-Body binary-black-hole dynamics and waveform family under `nrpy/equations/seobnr`; see [SEOBNR And BOB](equations/seobnr/index.md). |
| BOB | Backwards-One-Body waveform and attachment quantity family; see [BOB Waveforms](equations/seobnr/bob-waveforms.md). |
| GeneralRFM | General reference-metric support including fisheye coordinate maps; see [Geometry And Special-Function Support](equations/geometry-and-special-function-support.md). |
| Fishbone-Moncrief | Relativistic torus initial-data model implemented under GR equation modules; see [Fishbone-Moncrief](equations/general-relativity/fishbone-moncrief.md). |
| spin-weighted spherical harmonic | Special-function helper implemented with the Goldberg formula; see [Geometry And Special-Function Support](equations/geometry-and-special-function-support.md). |
| Method of Lines | Timestepping pattern used by NRPy infrastructures; ETLegacy registers evolved and RHS groups with the Einstein Toolkit MoL thorn, while BHaH examples use MoL infrastructure for standalone workflows. |
| trusted values | Generated reference dictionaries used to validate symbolic expressions; see [Trusted Expression Pipeline](equations/trusted-expression-pipeline.md). |
| trusted string files | Caller-derived text files written or compared by `validate_strings()` for generated string output; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| trusted expression dictionaries | Generated symbolic-expression dictionaries processed by `validate_expressions`, distinct from trusted string files; see [Trusted Expression Pipeline](equations/trusted-expression-pipeline.md). |
| validate_expressions | Module that processes and compares symbolic-expression validation dictionaries; see [Trusted Expression Pipeline](equations/trusted-expression-pipeline.md). |
| user cache directory | Platform-specific cache directory used by `cached_functions.py` for `.nrpycache` pickle files; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| .nrpycache | Hashed pickle cache file suffix used by `cached_functions.py`; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| PYTHONPATH | Environment variable that can include `.` for direct source-tree example runs; see [Build And Run](architecture/build-and-run.md). |
| Einstein Toolkit | External toolkit targeted by NRPy ETLegacy and CarpetX generators. |
| ETLegacy | Classic Einstein Toolkit / Carpet-style generated thorn target; see [ETLegacy](infrastructures/etlegacy/index.md). |
| CarpetX | Newer Einstein Toolkit driver-stack target generated by NRPy examples. |
| superB | Charm++-based distributed-memory infrastructure generated by NRPy examples; see [superB](infrastructures/superb/index.md). |
| Charm++ | Runtime/toolchain used by superB workflows; see [superB](infrastructures/superb/index.md). |
| chare | Charm++ object type used by superB generated runtime components; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| PUP | Charm++ pack/unpack serialization mechanism used for migration, checkpointing, and restart contracts; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| SDAG | Charm++ structured dagger notation for coordinating asynchronous entry-method control flow; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkIO | Charm++ parallel I/O library used by superB diagnostic output paths; see [Diagnostics And Observables](infrastructures/superb/diagnostics-and-observables.md). |
| JAX | Python array/accelerator ecosystem targeted by the `sebobv1_jax` generator. |

## Sources

- [README.md](../README.md) - `# NRPy`
- [README.md](../README.md) - `## Project Families and Example Generators`
- [README.md](../README.md) - `## What Gets Generated?`
- [CITATION.md](../CITATION.md) - `# How to Cite NRPy`
- [CITATION.md](../CITATION.md) - `## Quick reference`
- [c_function.py](../nrpy/c_function.py) - `CFunction`, `CFunction_dict`, `register_CFunction`.
- [c_codegen.py](../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`.
- [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py) - `cse_preprocess`, `cse_postprocess`, `sort_cse_output_deterministically`.
- [simd.py](../nrpy/helpers/simd.py) - `expr_convert_to_simd_intrins`.
- [simd_intrinsics.h](../nrpy/helpers/simd_intrinsics.h) - `REAL_SIMD_ARRAY`, `SIMD_WIDTH`.
- [cuda_intrinsics.h](../nrpy/helpers/cuda_intrinsics.h) - `REAL_CUDA_ARRAY`, `CUDA_WIDTH`.
- [generic.py](../nrpy/helpers/generic.py) - `clang_format`, `validate_strings`.
- [cached_functions.py](../nrpy/helpers/cached_functions.py) - `cache_file`, `cached_simplify`.
- [grid.py](../nrpy/grid.py) - `GridFunction`.
- [params.py](../nrpy/params.py) - `NRPyParameter`, `CodeParameter`.
- [reference_metric.py](../nrpy/reference_metric.py) - `ReferenceMetric`.
- [indexedexp.py](../nrpy/indexedexp.py) - `declarerank1`
- [validate_expressions.py](../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`
- [validate_expressions.py](../nrpy/validate_expressions/validate_expressions.py) - `compare_or_generate_trusted_results`
- [ADM_to_BSSN.py](../nrpy/equations/general_relativity/ADM_to_BSSN.py) - `ADM_to_BSSN`
- [GRHD_equations.py](../nrpy/equations/grhd/GRHD_equations.py) - `GRHD_Equations`
- [TOV_equations.py](../nrpy/equations/tov/TOV_equations.py) - `TOV_Equations`
- [psi4.py](../nrpy/equations/general_relativity/psi4.py) - `Psi4`
- [SEOBNRv5_aligned_spin_Hamiltonian.py](../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_Hamiltonian.py) - `SEOBNRv5_aligned_spin_Hamiltonian_quantities`
- [BOB_aligned_spin_waveform_quantities.py](../nrpy/equations/seobnr/BOB_aligned_spin_waveform_quantities.py) - `BOB_aligned_spin_waveform_quantities`
- [fisheye.py](../nrpy/equations/generalrfm/fisheye.py) - `GeneralRFMFisheye`
- [fishbone_moncrief.py](../nrpy/equations/general_relativity/fishbone_moncrief/fishbone_moncrief.py) - `FishboneMoncriefID`
- [spin_weighted_spherical_harmonics.py](../nrpy/equations/special_functions/spin_weighted_spherical_harmonics.py) - `Y`
- [Charm++ language manual](https://github.com/charmplusplus/charm/blob/main/doc/charm%2B%2B/manual.rst) - `Execution Model`, `Structured Control Flow: Structured Dagger`, `Serialization Using the PUP Framework`
- [Charm++ and Converse libraries manual](https://github.com/charmplusplus/charm/blob/main/doc/libraries/manual.rst) - `CkIO`, `Using CkIO`, `Parallel Output API`

## See Also

- [Schema](SCHEMA.md)
- [Workflows](workflows.md)
- [Architecture Overview](architecture/overview.md)
