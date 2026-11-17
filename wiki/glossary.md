# Glossary

> Canonical spellings for recurring code, domain, and experiment entities. · Status: confirmed

## Summary

Use these spellings when routing, citing, or writing leaves. Terms link to their defining pages where NRPy defines them.

## Terms

| Term | Meaning |
| --- | --- |
| NRPy | Python/SymPy symbolic code-generation package for numerical relativity and relativistic astrophysics; see [Architecture Overview](architecture/overview.md). |
| SymPy | Python symbolic algebra package used to build expressions before code emission; see [Architecture Overview](architecture/overview.md). |
| generated C | C emitted by NRPy codegen and infrastructure modules rather than handwritten source; see [Generated Output Boundaries](architecture/generated-output-boundaries.md). |
| BHaH | NRPy's standalone generated-application infrastructure for BH@H-style executables and libraries; see [BHaH](infrastructures/bhah/index.md). |
| BH@H | BlackHoles@Home single-patch numerical-relativity workflows represented by BHaH examples; see [Standalone GR/BHaH](examples/standalone-gr-bhah.md). |
| BHaHAHA | Apparent-horizon finding library integrated into BHaH, superB, and Dendro black-hole workflows; see [BHaHAHA](infrastructures/bhah/bhahaha/index.md). |
| CFunction | Python representation of a generated C function; see [C Function Registry](core/c-function-registry.md). |
| CFunction_dict | Global registry mapping generated function names to `CFunction` objects; see [C Function Registry](core/c-function-registry.md). |
| register_CFunction | Function that registers generated C functions in `CFunction_dict`; see [C Function Registry](core/c-function-registry.md). |
| c_codegen | Main SymPy-to-C emission function; see [C Codegen](core/c-codegen.md). |
| CCodeGen | Option class configuring lower-level C code emission; see [C Codegen](core/c-codegen.md). |
| PyCodeGen | Option class configuring JAX-compatible Python assignment emission; see [Python Codegen](core/python-codegen.md). |
| py_codegen | Main SymPy-to-JAX-compatible Python assignment emission function; see [Python Codegen](core/python-codegen.md). |
| PyFunction | Python representation of generated Python/JAX-compatible function text; see [Python Function Registry](core/python-function-registry.md). |
| PyFunction_dict | Global registry mapping generated Python function names to `PyFunction` objects; see [Python Function Registry](core/python-function-registry.md). |
| CSE | Common-subexpression elimination preprocessing, postprocessing, and deterministic temporary ordering; see [CSE And Printer Support](core/helpers/cse-and-printer-support.md). |
| SIMD | Single instruction, multiple data symbolic `*SIMD` rewrites and handwritten intrinsic macro headers; see [SIMD And Intrinsic Support](core/helpers/simd-and-intrinsic-support.md). |
| CUDA | NVIDIA GPU target for generated kernels, launches, and device intrinsics; see [Loop Kernel And Device Helpers](core/helpers/loop-kernel-and-device-helpers.md). |
| intrinsics | Low-level vector/device operation macros in NRPy's SIMD and CUDA helper headers; see [SIMD And Intrinsic Support](core/helpers/simd-and-intrinsic-support.md). |
| clang-format | External C/C++ formatter invoked by `nrpy.helpers.generic.clang_format`; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| NRPyParameter | Generation-time NRPy parameter object; see [Gridfunctions And Parameters](core/gridfunctions-and-parameters.md). |
| CodeParameter | Runtime parameter emitted into generated project structs and parfiles; see [Gridfunctions And Parameters](core/gridfunctions-and-parameters.md). |
| GridFunction | Registered evolved, auxiliary, or diagnostic grid quantity; see [Gridfunctions And Parameters](core/gridfunctions-and-parameters.md). |
| commondata | Grid-shared runtime data emitted in `commondata_struct` or `Commondata`; see [Gridfunctions And Parameters](core/gridfunctions-and-parameters.md). |
| reference metric | Coordinate-system metric, scale-factor, and Jacobian quantities in `rfm_dict`; see [Reference Metrics](core/reference-metrics.md). |
| inner boundary point | Ghost point mapped by `bcstruct_set_up()` to an interior source point with parity sign; see [Curvilinear Boundaries](infrastructures/bhah/curvilinear-boundaries-and-stencils.md). |
| indexed expression | Tensor or vector array built with `nrpy.indexedexp`; see [Indexed Expressions](core/indexed-expressions.md). |
| ADM | Arnowitt-Deser-Misner 3+1 spatial metric and extrinsic curvature quantities; see [Metric Conversions And Matter](equations/general-relativity/metric-conversions-and-matter.md). |
| BSSN | Baumgarte-Shapiro-Shibata-Nakamura reference-metric evolution formulation; see [BSSN Family](equations/general-relativity/bssn-family.md). |
| YBS-MOM | Opt-in Yo--Lin--Cao momentum-gradient term (`C_YBS_mom`, `W` weight, timestep-scaled); see [YBS-MOM Adjustment](equations/general-relativity/ybs-momentum-damping.md). |
| GRHD | General relativistic hydrodynamics conserved variables, fluxes, sources, speeds, and HLL helpers; see [GRHD](equations/grhd.md). |
| TOV | Tolman-Oppenheimer-Volkoff stellar-equilibrium ODE system; see [TOV Equations](equations/tov-equations.md). |
| Psi4 | Weyl gravitational-radiation scalar and tetrad construction modules; see [Psi4 And Tetrads](equations/general-relativity/psi4-and-tetrads.md). |
| SEOBNR | Effective-One-Body binary-black-hole dynamics and waveform modules under `nrpy/equations/seobnr`; see [SEOBNR And BOB](equations/seobnr/index.md). |
| BOB | Backwards-One-Body merger-ringdown waveform and attachment quantities; see [BOB Waveforms](equations/seobnr/bob-waveforms.md). |
| GeneralRFM | General reference-metric coordinate maps including fisheye stretching; see [Geometry And Special-Function Support](equations/geometry-and-special-function-support.md). |
| Fishbone-Moncrief | Relativistic equilibrium torus initial-data model; see [Fishbone-Moncrief](equations/general-relativity/fishbone-moncrief.md). |
| spin-weighted spherical harmonic | Special-function helper `_s Y_{lm}` implemented via the Goldberg formula; see [Geometry And Special-Function Support](equations/geometry-and-special-function-support.md). |
| Method of Lines | Spatial-discretization plus ODE-integrator time stepping across BHaH, ETLegacy, superB, and Dendro; see [MoL Time Integration](infrastructures/bhah/mol-time-integration.md). |
| test oracle | Reviewed reference values or expected generated output checked by an executable test; see [Test Oracles And Safe Updates](validation/test-oracles-and-safe-updates.md). |
| executable test | Executed assertion or numerical comparison verifying required module behavior; see [Code Test Policy](validation/code-test-policy.md). |
| trusted values | Reference numerical dictionaries (`trusted_dict`) validating symbolic expressions; see [Expression Validation Helpers](validation/expression-validation-helpers.md). |
| trusted string files | Caller-named reference text files compared or updated by `validate_strings()`; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| trusted expression dictionaries | Symbolic-expression reference dictionaries processed by `validate_expressions`; see [Expression Validation Helpers](validation/expression-validation-helpers.md). |
| Expression validation helpers | Numerical comparison, zero-check, and `trusted_dict` comparison utilities; see [Expression Validation Helpers](validation/expression-validation-helpers.md). |
| validate_expressions | Module that evaluates and compares symbolic-expression validation dictionaries; see [Expression Validation Helpers](validation/expression-validation-helpers.md). |
| user cache directory | Platform directory storing `.nrpycache` pickle files for `cached_functions.py`; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| .nrpycache | Hashed pickle cache file suffix written by `cached_functions.py`; see [Maintenance And Validation Helpers](core/helpers/maintenance-and-validation-helpers.md). |
| PYTHONPATH | Environment variable (`PYTHONPATH=.`) enabling direct repository module execution; see [Build And Run](architecture/build-and-run.md). |
| nrpyinline | Console script (`bin/nrpyinline.py`) executing trusted embedded Python blocks in one shared namespace; see [Build And Run](architecture/build-and-run.md). |
| Einstein Toolkit | External numerical-relativity framework targeted by NRPy ETLegacy and CarpetX thorn generators. |
| ETLegacy | Classic Einstein Toolkit / Carpet generated-thorn backend; see [ETLegacy](infrastructures/etlegacy/index.md). |
| CarpetX | AMReX-based Einstein Toolkit generated-thorn backend; see [CarpetX](infrastructures/carpetx/index.md). |
| superB | Charm++-based distributed-memory multi-patch infrastructure; see [superB](infrastructures/superb/index.md). |
| Charm++ | Parallel C++ runtime and compiler toolchain used by superB applications; see [superB](infrastructures/superb/index.md). |
| .ci | Charm++ interface definition file emitted by superB generators; see [superB Project Assembly](infrastructures/superb/lifecycle-and-project-assembly.md). |
| mainchare | Charm++ startup chare class (`Main`) in superB executables; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| chare | Charm++ message-driven object class in superB executables; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| readonly | Charm++ startup-broadcast global handle used for proxy publication in superB; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CProxy_* | Generated Charm++ proxy handle class invoking chare entry methods; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkIndex_* | Generated Charm++ entry-method index namespace used in callbacks and reductions; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkArrayIndex3D | Charm++ 3D chare-array index addressing `Timestepping` and `Interpolator3d` elements; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkReduction | Charm++ reduction namespace for volume integrals and completion barriers; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkReductionMsg | Charm++ message delivering reduced scalar arrays and CkIO callbacks; see [Diagnostics And Observables](infrastructures/superb/diagnostics-and-observables.md). |
| CkCallback | Charm++ callback descriptor for reductions, checkpoints, and CkIO sessions; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| InterpBufMsg | superB Charm++ message carrying BHaHAHA or `Psi4` interpolation buffers; see [GR, BHaHAHA, Psi4, And Interpolation](infrastructures/superb/gr-bhahaha-psi4-and-interpolation.md). |
| CkArrayMap | Charm++ chare-placement base class used by superB's `RRMap_with_offset`; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| PUP | Charm++ pack/unpack serialization interface for migration, checkpointing, and restart; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| SDAG | Charm++ structured dagger control-flow blocks inside `.ci` entry methods; see [Chare Entrypoints And Runtime](infrastructures/superb/chare-entrypoints-and-runtime.md). |
| CkIO | Charm++ parallel file-output library used by superB 1D/2D diagnostics; see [Diagnostics And Observables](infrastructures/superb/diagnostics-and-observables.md). |
| JAX | Python array/JIT library targeted by `nrpy/infrastructures/JAX`; see [JAX](infrastructures/jax/index.md). |
| Dendro | Standalone octree-AMR C++/MPI backend emitting `Dendro_NRPy_BSSN` and `Dendro_NRPy_fCCZ4` on Dendrolib; see [Dendro](infrastructures/dendro/index.md). |
| NRPyElliptic | Hyperbolic-relaxation conformally flat elliptic initial-data solver across equations, BHaH, and superB; see [NRPyElliptic Workflow](infrastructures/bhah/nrpyelliptic-workflow.md). |
| SEBOB | SEOBNR-plus-BOB binary-black-hole waveform generators across BHaH C libraries and JAX; see [SEOBNR BOB Generated Library](infrastructures/bhah/seobnr-bob-generated-library.md). |
| Kreiss-Oliger | Artificial dissipation via `dKOD` stencils in evolution RHSs; see [Finite Difference](core/finite-difference.md) and [BSSN Family](equations/general-relativity/bssn-family.md). |
| GRoovy | GRHD evolution executables generated through BHaH matter modules; see [GRoovy GRHD Runtime](infrastructures/bhah/groovy-grhd-runtime.md). |
| TOVola | Neutron-star TOV initial-data solver integrated into BHaH matter examples; see [Matter TOV Workflows](examples/matter-tov-workflows.md). |
| GRHayL | External C GRHD library required to build generated GRoovy executables; see [GRoovy GRHD Runtime](infrastructures/bhah/groovy-grhd-runtime.md). |
| SO(3) | 3D rotation helpers (`SO3_rotations.py`, `BHaH/rotation/`); `SO(3)` in prose, `SO3` in code; see [SO(3) Rotation Helpers](infrastructures/bhah/so3-rotation-helpers.md). |

## See Also

- See also: [KB Rules](SCHEMA.md)
- See also: [Workflows](workflows.md)
- See also: [Architecture Overview](architecture/overview.md)
