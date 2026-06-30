# ETLegacy Code Parameters, Includes, And Loops

> ETLegacy-local helpers for Cactus parameter reads, standard includes, and generated grid loops. · Status: confirmed · Last reconciled: 2026-06-30
> Up: [ETLegacy](index.md)

## Summary

ETLegacy generated C functions use small infrastructure-local helpers for three
recurring code fragments: reading Cactus runtime parameters, attaching standard
Cactus headers, and wrapping kernel bodies in Cactus grid loops. These helpers
stay ETLegacy-specific: generic symbolic expression emission belongs to
[C Codegen](../../core/c-codegen.md), CodeParameter registration belongs to
[Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md),
SIMD mechanics belong to
[SIMD And Intrinsic Support](../../core/helpers/simd-and-intrinsic-support.md),
and generic loop text belongs to
[Loop Kernel And Device Helpers](../../core/helpers/loop-kernel-and-device-helpers.md).

## Detail

`read_CodeParameters()` emits C text that gets Cactus runtime parameter values
with `CCTK_ParameterGet`. Callers pass `(thorn, code-parameter-name)` tuples;
the helper sorts them, looks up each parameter name in
`par.glb_code_params_dict`, and raises if the parameter was not registered. For
`CCTK_REAL` and `REAL`, generated Cactus lookups use the caller-provided tuple
thorn. For other non-char types, generated lookups use the registered
`CodeParameter.module`, while emitted comments still record the tuple thorn.

For `CCTK_REAL` and `REAL` code parameters, scalar mode emits a typed pointer
from `CCTK_ParameterGet` and then dereferences it into a `const CCTK_REAL`
local. SIMD mode emits a `NOSIMD<name>` pointer and wraps the dereferenced
value with `ConstSIMD` into a `REAL_SIMD_ARRAY`. Parameters named like
`invdxx0`, `invdxx1`, or `invdxx2` are special-cased as inverse grid-spacing
constants derived from `CCTK_DELTA_SPACE(<dir>)` instead of
`CCTK_ParameterGet`. When `declare_invdxxs` is enabled, the helper appends all
three inverse spacing declarations; in SIMD mode it emits both scalar
`NOSIMD*` values and vector constants. Code-parameter types containing `char`
are rejected before any scalar or SIMD declaration is emitted.

Other registered code-parameter types do not use the numeric SIMD wrapper path.
They are emitted as `const <type> <name> = CCTK_ParameterGet(...)` using the
registered `CodeParameter` metadata. This keeps ETLegacy code generation tied
to registered code parameters and Cactus thorn/module parameter namespaces,
while leaving the core CodeParameter object model documented in
[Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md).

`define_standard_includes()` centralizes the default include list for
NRPy-generated ETLegacy C functions. It returns `math.h`, `cctk.h`,
`cctk_Arguments.h`, and `cctk_Parameters.h`, so ETLegacy kernels can request
the same Cactus/standard C headers without repeating that list in every
registration site.

`simple_loop()` wraps a kernel body in a three-dimensional Cactus grid loop by
delegating final loop-string construction to the generic `nrpy.helpers.loop`
helper. With `loop_region="all points"`, bounds run from zero to
`cctk_lsh[2]`, `cctk_lsh[1]`, and `cctk_lsh[0]`. With
`loop_region="interior"`, bounds start at `cctk_nghostzones[dir]` and stop at
`cctk_lsh[dir]-cctk_nghostzones[dir]`, excluding ghost zones. Unsupported loop
regions raise `ValueError`.

OpenMP behavior is selected locally before the generic loop helper is called.
Default OpenMP emits `#pragma omp parallel for`; `OMP_collapse > 1` adds a
`collapse(<n>)` clause; a nonempty `OMP_custom_pragma` replaces the default
pragma; and disabling OpenMP with no custom pragma emits no pragma. SIMD mode
changes only the innermost `i0` loop stride from `1` to `SIMD_WIDTH`; the
outer `i2` and `i1` loops keep unit stride.

## Sources

- [CodeParameters.py](../../../nrpy/infrastructures/ETLegacy/CodeParameters.py) - `read_CodeParameters`
- [ETLegacy_include_header.py](../../../nrpy/infrastructures/ETLegacy/ETLegacy_include_header.py) - `define_standard_includes`
- [simple_loop.py](../../../nrpy/infrastructures/ETLegacy/simple_loop.py) - `simple_loop`
- [params.py](../../../nrpy/params.py) - `CodeParameter`, `glb_code_params_dict`
- [loop.py](../../../nrpy/helpers/loop.py) - `loop`

## See Also

- [ETLegacy](index.md)
- [Thorn Assembly And CCL Files](thorn-assembly-and-ccl-files.md)
- [MoL, Boundaries, Symmetry, And RHS Initialization](mol-boundaries-symmetry-and-rhs-initialization.md)
- [C Codegen](../../core/c-codegen.md)
- [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- [SIMD And Intrinsic Support](../../core/helpers/simd-and-intrinsic-support.md)
- [Loop Kernel And Device Helpers](../../core/helpers/loop-kernel-and-device-helpers.md)
