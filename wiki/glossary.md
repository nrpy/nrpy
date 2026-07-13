# Glossary

> Canonical spellings for recurring code, domain, and experiment entities. · Status: confirmed · Last reconciled: 07-12-2026

## Summary

Use these spellings when routing, citing, or adding leaves. Terms link to owning
pages where a first-pass owner exists.

## Terms

| Term | Meaning |
| --- | --- |
| NRPy | Python/SymPy-based symbolic code generation toolkit for numerical relativity and relativistic astrophysics; see [Architecture Overview](architecture/overview.md). |
| SymPy | Python symbolic algebra package used to build expressions before C generation; see [Symbolic Codegen Lifecycle](architecture/symbolic-codegen-lifecycle.md). |
| generated C | C emitted by NRPy generators and infrastructure rather than handwritten source; see [Generated Output Boundaries](architecture/generated-output-boundaries.md). |
| BHaH | NRPy's standalone generated-application infrastructure for BH@H-style workflows; see [BHaH](infrastructures/bhah/index.md). |
| BH@H | BlackHoles@Home, the single-patch numerical relativity workflow represented by current BHaH examples; see [Black Hole Evolution](examples/black-hole-evolution.md). |
| BHaHAHA | Apparent-horizon finding library integrated by BHaH black-hole workflows; see [BHaHAHA Horizon Runtime](infrastructures/bhah/bhahaha-horizon-runtime.md). |
| CFunction | Python representation of a generated C function; see [C Function Registry](core/c-function-registry.md). |
| CFunction_dict | Global registry mapping generated function names to `CFunction` objects; see [C Function Registry](core/c-function-registry.md). |
| register_CFunction | API that registers generated C functions; see [C Function Registry](core/c-function-registry.md). |
| c_codegen | Main SymPy-to-C emission function; see [C Codegen](core/c-codegen.md). |
| CCodeGen | Codegen class that backs lower-level C code emission; see [C Codegen](core/c-codegen.md). |
| PyCodeGen | Option object behind JAX-compatible Python assignment emission; see [Python Codegen](core/python-codegen.md). |
| py_codegen | Main SymPy-to-JAX-compatible Python assignment emission function; see [Python Codegen](core/python-codegen.md). |
| PyFunction | Python representation of generated Python/JAX-compatible function text; see [Python Function Registry](core/python-function-registry.md). |
| PyFunction_dict | Global registry mapping generated Python function names to `PyFunction` objects; see [Python Function Registry](core/python-function-registry.md). |
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
| Method of Lines | Timestepping pattern used by NRPy infrastructures; ETLegacy registers evolved and RHS groups with the Einstein Toolkit MoL thorn, while BHaH examples use MoL infrastructure for standalone workflows; see [MoL Time Integration](infrastructures/bhah/mol-time-integration.md). |
| test oracle | Reviewed data or expected artifact used by an executable owner validation path; presence or missing-file creation is not executable coverage; see [Test Oracles And Safe Updates](validation/test-oracles-and-safe-updates.md) and [Code Test Policy](validation/code-test-policy.md). |
| executable test | Executed assertion or comparison path that checks a durable production contract, distinct from stored oracle data; see [Code Test Policy](validation/code-test-policy.md) and [Test Oracles And Safe Updates](validation/test-oracles-and-safe-updates.md). |
| trusted values | Generated reference dictionaries used to validate symbolic expressions; see [Trusted Expression Pipeline](equations/trusted-expression-pipeline.md). |
| trusted string files | Caller-derived text files written or compared by `validate_strings()` for generated string output; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| trusted expression dictionaries | Generated symbolic-expression dictionaries processed by `validate_expressions`, distinct from trusted string files; see [Trusted Expression Pipeline](equations/trusted-expression-pipeline.md). |
| Expression validation helpers | Functions that process, compare, generate, or numerically check trusted symbolic-expression dictionaries; see [Expression Validation Helpers](validation/expression-validation-helpers.md). |
| validate_expressions | Module that processes and compares symbolic-expression validation dictionaries; see [Trusted Expression Pipeline](equations/trusted-expression-pipeline.md). |
| user cache directory | Platform-specific cache directory used by `cached_functions.py` for `.nrpycache` pickle files; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| .nrpycache | Hashed pickle cache file suffix used by `cached_functions.py`; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| PYTHONPATH | Environment variable that can include `.` for direct source-tree example runs; see [Build And Run](architecture/build-and-run.md). |
| nrpyinline | Installed console script backed by `bin/nrpyinline.py`; it executes trusted embedded Python snippets from a file in one shared namespace and provides no sandbox; see [Build And Run](architecture/build-and-run.md). |
| Einstein Toolkit | External toolkit targeted by NRPy ETLegacy and CarpetX generators. |
| ETLegacy | Classic Einstein Toolkit / Carpet-style generated thorn target; see [ETLegacy](infrastructures/etlegacy/index.md). |
| CarpetX | Newer Einstein Toolkit driver-stack target generated by NRPy examples; see [CarpetX](infrastructures/carpetx/index.md). |
| superB | Charm++-based distributed-memory infrastructure generated by NRPy examples; see [superB](infrastructures/superb/index.md). |
| Charm++ | Runtime/toolchain used by superB workflows; see [superB](infrastructures/superb/index.md). |
| .ci | Charm++ interface file emitted by superB generators; see [Lifecycle And Project Assembly](infrastructures/superb/lifecycle-and-project-assembly.md). |
| mainchare | Charm++ startup chare role used by superB's generated `Main`; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| chare | Charm++ object type used by superB generated runtime components; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| readonly | Charm++ startup-published global or handle used by superB for proxy publication; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CProxy_* | Generated Charm++ proxy-handle class family used to invoke chare entry methods; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkIndex_* | Generated Charm++ entry-point index family used in callbacks and reduction targets; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkArrayIndex3D | Charm++ three-dimensional array-coordinate wrapper used to address superB `Timestepping` and `Interpolator3d` elements; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkReduction | Charm++ reduction API namespace/type family used by superB volume data and completion reductions; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkReductionMsg | Charm++ callback message used for reduction payloads and CkIO completion callbacks; see [Diagnostics And Observables](infrastructures/superb/diagnostics-and-observables.md). |
| CkCallback | Charm++ callback wrapper used by superB reductions, checkpoint coordination, and CkIO completion paths; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| InterpBufMsg | superB variable-size interpolation message carrying packed BHaHAHA or Psi4 interpolation buffers; see [GR, BHaHAHA, Psi4, And Interpolation](infrastructures/superb/gr-bhahaha-psi4-and-interpolation.md). |
| CkArrayMap | Charm++ array-placement map abstraction used by superB's optional `RRMap_with_offset`; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| PUP | Charm++ pack/unpack serialization mechanism used for migration, checkpointing, and restart contracts; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| SDAG | Charm++ structured dagger notation for coordinating asynchronous entry-method control flow; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkIO | Charm++ parallel I/O library used by superB diagnostic output paths; see [Diagnostics And Observables](infrastructures/superb/diagnostics-and-observables.md). |
| JAX | Generated Python/JAX infrastructure target used by `sebobv1_jax`; see [JAX](infrastructures/jax/index.md). |
| NRPyElliptic | Hyperbolic-relaxation elliptic initial-data solver family spanning equation modules, BHaH workflows, and superB examples; see [NRPyElliptic Workflow](infrastructures/bhah/nrpyelliptic-workflow.md). |
| SEBOB | SEOBNR-plus-BOB waveform example family (SEBOBv1, SEBOBv2) spanning generated C library and JAX routes; see [SEOBNR BOB Generated Library](infrastructures/bhah/seobnr-bob-generated-library.md). |
| Kreiss-Oliger | Numerical dissipation scheme applied through `dKOD` finite-difference operators in evolution RHS wiring; see [Finite Difference](core/finite-difference.md). |
| GRoovy | GRHD evolution code family generated through BHaH matter workflows; see [GRoovy GRHD Runtime](infrastructures/bhah/groovy-grhd-runtime.md). |
| TOVola | TOV initial-data solver used by matter workflow examples; see [Matter TOV Workflows](examples/matter-tov-workflows.md). |
| GRHayL | External GRHD library that is a hard generated-project build dependency for GRoovy examples; see [GRoovy GRHD Runtime](infrastructures/bhah/groovy-grhd-runtime.md). |
| SO(3) | Canonical spelling for the rotation-group helpers (`SO3_rotations.py` symbolic expressions and generated BHaH helper CFunctions); prefer `SO(3)` in prose, `SO3` only in code symbols; see [SO(3) Rotation Helpers](infrastructures/bhah/so3-rotation-helpers.md). |

## Sources

- [README.md](../README.md) - `# NRPy`
- [README.md](../README.md) - `## Project Families and Example Generators`
- [README.md](../README.md) - `## What Gets Generated?`
- [CITATION.md](../CITATION.md) - `# How to Cite NRPy`
- [CITATION.md](../CITATION.md) - `## Quick reference`
- [c_function.py](../nrpy/c_function.py) - `CFunction`, `CFunction_dict`, `register_CFunction`.
- [c_codegen.py](../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`.
- [py_codegen.py](../nrpy/py_codegen.py) - `PyCodeGen`, `py_codegen`.
- [py_function.py](../nrpy/py_function.py) - `PyFunction`, `PyFunction_dict`, `register_PyFunction`.
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
- [main_chare.py](../nrpy/infrastructures/superB/main_chare.py) - `mainchare Main`, readonly proxies, `RRMap_with_offset`.
- [timestepping_chare.py](../nrpy/infrastructures/superB/timestepping_chare.py) - `CkIndex_Timestepping::report_sums_for_volume`, `CkReduction::sum_double`, CkIO callbacks.
- [interpolator3d_chare.py](../nrpy/infrastructures/superB/interpolator3d_chare.py) - `InterpBufMsg`, `CkArrayIndex3D`.
- [Charm++ language manual](https://charm.readthedocs.io/en/v8.0.0/charm%2B%2B/manual.html) - `Charm++ Interface (.ci) Files`, `Execution Model`, readonly variables, reductions, array maps, structured dagger, PUP.
- [Charm++ and Converse libraries manual](https://charm.readthedocs.io/en/v8.0.0/libraries/manual.html) - `CkIO`, `Using CkIO`, `Parallel Output API`

## See Also

- [Schema](SCHEMA.md)
- [Workflows](workflows.md)
- [Architecture Overview](architecture/overview.md)
