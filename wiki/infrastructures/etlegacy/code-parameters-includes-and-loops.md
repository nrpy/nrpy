# ETLegacy Code Parameters, Includes, And Loops

> ETLegacy-local helpers for Cactus parameter reads, standard includes, and generated grid loops. · Status: confirmed · Last reconciled: 08-18-2026
> Up: [ETLegacy](index.md)

## Summary

ETLegacy generated C functions use small infrastructure-local helpers for
reading Cactus runtime parameters, handing expression-derived parameters to
generated kernels, attaching standard Cactus headers, and wrapping kernel
bodies in Cactus grid loops. These helpers stay ETLegacy-specific: generic
symbolic expression emission belongs to
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

The BSSN RHS and constraint registrations pass their final expression lists to
`get_params_commondata_symbols_from_expr_list()` and use only its
non-commondata parameter result. In SIMD mode, that list is converted to
current-thorn `(thorn, parameter)` tuples for `read_CodeParameters()`; in
scalar mode, the generated function instead uses `DECLARE_CCTK_PARAMETERS`.
The same non-commondata list is passed into current-thorn CFunction parameter
metadata, with the RHS retaining its separate manually assembled metadata for
parameters used outside the final expression list.

No parameter name or matter flag is special-cased in this handoff. `PI` appears
only when the selected T4munu expressions contain the registered,
non-commondata `PI` CodeParameter; vacuum expressions do not contribute it.

Claim evidence:
- Claim: ETLegacy BSSN RHS and constraint registrations derive non-commondata CodeParameters from their final expressions, use that same list for current-thorn metadata and SIMD `read_CodeParameters()` declarations, use `DECLARE_CCTK_PARAMETERS` in scalar mode, and do not select `PI` by parameter name or matter flag.
- Role: descriptive behavior
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`; [BSSN_constraints.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`
- Corroboration: [expression_utils.py](../../../nrpy/helpers/expression_utils.py), `get_params_commondata_symbols_from_expr_list`; [CodeParameters.py](../../../nrpy/infrastructures/ETLegacy/CodeParameters.py), `read_CodeParameters`
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3; backend=ETLegacy; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=Cartesian, fd_order 4, RHS KO enabled, reference-metric precompute disabled, T4munu False/True, SIMD False/True; date=08-18-2026`

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
- [expression_utils.py](../../../nrpy/helpers/expression_utils.py) - `get_params_commondata_symbols_from_expr_list`
- [rhs_eval.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`
- [BSSN_constraints.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_constraints.py) - `register_CFunction_BSSN_constraints`
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
- Depends on: [Symbolic Expression Utilities](../../core/helpers/symbolic-expression-utilities.md)
- [SIMD And Intrinsic Support](../../core/helpers/simd-and-intrinsic-support.md)
- [Loop Kernel And Device Helpers](../../core/helpers/loop-kernel-and-device-helpers.md)
