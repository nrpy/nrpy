# Gridfunctions And Parameters

> Core route for symbolic gridfunctions, NRPy parameters, code parameters, and generated data structs. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Core APIs](index.md)

## Summary

Gridfunctions and parameters are the main bridge between symbolic expressions and generated infrastructure data. `nrpy.grid` registers symbolic gridfunction names and infrastructure-specific metadata. `nrpy.params` stores Python-side `NRPyParameter` values and C-side `CodeParameter` symbols. BHaH grid/common data registration records declarations that later become members of generated `griddata_struct` and `commondata_struct` support code.

## Detail

### NRPy Parameters

`NRPyParameter` stores Python-level configuration in `glb_params_dict`. Valid parameter Python types are exactly `bool`, `int`, `float`, and `str`; other types raise during object construction. `register_param()` is the public add-if-missing wrapper: construction inserts a parameter only when its bare `name` is absent from `glb_params_dict`; duplicate registration leaves the existing object in place, with only an optional verbose warning. `parval_from_str()` and `set_parval_from_str()` both accept either `name` or `module::name` and strip any module prefix by using the final token after `::`.

### Default NRPy Parameters

`params.py` registers four default Python-side NRPy parameters: `Infrastructure` defaults to `BHaH`, `fp_type` defaults to `double`, `parallelization` defaults to `openmp`, and `clang_format_options` defaults to `-style={BasedOnStyle: LLVM, ColumnLimit: 150}`. Infrastructure selection matters directly to gridfunction registration because `register_gridfunctions()` reads `Infrastructure` before selecting a gridfunction subclass.

### CodeParameter Core Rules

`CodeParameter` stores runtime C parameters in `glb_code_params_dict` and creates the SymPy symbol returned by `register_CodeParameter()` or `register_CodeParameters()`. Supported symbolic assumptions are `Real`, which creates a real symbol, and `RealPositive`, which creates a real positive symbol. Core routing flags are `commondata` for common-across-grids storage, `add_to_parfile` for parameter-file exposure, `add_to_set_CodeParameters_h` for BHaH `set_CodeParameters*.h` support, and `add_to_glb_code_params_dict` for global registry insertion. After parameter-type processing, non-array parameters with `add_to_parfile=True` and final default value `"unset"` raise `ValueError`; array defaults are first normalized to lists, so the current source does not enforce that exact `"unset"` check for array entries.

### `#define` CodeParameters

A `CodeParameter` whose `cparam_type` is `#define` is treated as a C-preprocessor constant, not a runtime-settable parameter. Construction forces both `add_to_parfile` and `add_to_set_CodeParameters_h` to `False`, regardless of caller input, before normal symbol creation and optional global registration.

### Array CodeParameters

Array `CodeParameter` types are parsed only for C-style `REAL[N]` and `int[N]` strings. `_parse_array_spec()` extracts `N`; invalid or unsupported strings do not enter the array branch. Array parameters must use `add_to_set_CodeParameters_h=False`. If the default value is a list, its length must equal `N`; if the default is scalar, construction broadcasts that scalar to a list of length `N`. `adjust_CodeParam_default("name[i]", value)` updates one indexed default in an existing list-valued array parameter, and `adjust_CodeParam_default("name", value, new_cparam_type)` updates the whole default and optionally the C type.

### Batch Versus Single Registration

`register_CodeParameter()` accepts one simple string name and returns one symbol. Direct `CodeParameter` construction and this single-name wrapper insert into `glb_code_params_dict` only if the name is absent. `register_CodeParameters()` requires a non-empty list of names, broadcasts a scalar `defaultvalues` input to all names, validates list lengths for defaults and descriptions, rejects a non-empty scalar description string, and returns one symbol per name. When `add_to_glb_code_params_dict=True`, the batch helper explicitly assigns each newly built object into `glb_code_params_dict`, so it overwrites an existing global entry.

### GridFunction Base Contract

`GridFunction` is the base object for registered gridfunctions. It stores name, group, description, rank, dimension, C type, asymptotic value `f_infinity`, characteristic `wavespeed`, and whether the name is a basename. If `is_basename=True`, the basename must be a non-empty string and must not end in a digit; this keeps scalar basenames distinct from component names such as `betU2` or `gDD12`.

### Registration Behavior

`register_gridfunctions()` normalizes one name or a list of names, reads `Infrastructure`, selects the subclass from `GF_CLASS_MAP`, stores new gridfunctions in `glb_gridfcs_dict`, and returns real SymPy symbols. Supported infrastructure keys are `BHaH`, `ETLegacy`, and `CarpetX`; an unknown key raises `ValueError`. Duplicate registration prints `Warning: Gridfunction <name> is already registered.`, leaves the existing registry object unchanged, and still returns a real SymPy symbol for the requested name. Per-gridfunction list-valued `f_infinity` and `wavespeed` inputs are indexed by position during registration; the helper does not prevalidate those list lengths.

### Parity

`get_parity_type()` returns parity code `0` for rank-0 scalars. For rank-1 components, it reads the final component digit and returns component index plus one. For rank-2 components, it reads the final two component digits and consults fixed lookup tables for dimension 3 or dimension 4. `set_parity_types()` looks up each requested name in `glb_gridfcs_dict`, computes its parity from rank and dimension, and raises if a name is missing, rank or dimension is unsupported, or the resulting parity list length does not match the input list.

### BHaH Access

`BHaHGridFunction` accepts groups `EVOL`, `AUXEVOL`, `DIAG`, and `AUX`, uses C type `REAL`, and maps default arrays by group: `EVOL` to `in_gfs`, `AUXEVOL` to `auxevol_gfs`, `DIAG` to `diagnostic_gfs`, and `AUX` to `aux_gfs`. Memory access has spelling `array[IDX4(<NAME>GF, i0+offset, i1+offset, i2+offset)]`, for example `in_gfs[IDX4(ABCGF, i0+1, i1+2, i2+3)]`; `enable_simd=True` wraps it as `ReadSIMD(&...)`. If `sync_gf_in_superB` is not supplied, it defaults to `True` for `EVOL` and `AUX`. `gridfunction_defines()` emits group `#define` blocks and, for evolved gridfunctions, `gridfunctions_f_infinity` and `gridfunctions_wavespeed` arrays.

### ETLegacy Access

`ETLegacyGridFunction` accepts groups `EVOL`, `AUX`, and `AUXEVOL`, uses C type `CCTK_REAL`, and emits memory access through `CCTK_GFINDEX3D`, for example `abcGF[CCTK_GFINDEX3D(cctkGH, i0+1, i1+2, i2+3)]`. `use_GF_suffix=False` removes the `GF` suffix from the callable name, and `enable_simd=True` wraps the access as `ReadSIMD(&...)`.

### CarpetX Access

`CarpetXGridFunction` accepts groups `EVOL`, `AUX`, `AUXEVOL`, `EXTERNAL`, `CORE`, `TILE_TMP`, and `SCALAR_TMP`, uses C type `CCTK_REAL`, and validates centering as a three-character string containing only `C` or `V`. During construction it appends `_ext` for `EXTERNAL`, `_core` for `CORE`, and `_tile_tmp` for `TILE_TMP`; other groups keep the original name. Memory access calls the gridfunction with an index expression based on `p.I + offset*p.DI[...]`, for example `aaGF(p.I + 1*p.DI[0] + 2*p.DI[1] + 3*p.DI[2])` or `defgGF(p.I - 1*p.DI[1])`. `reuse_index=True` uses the provided `index_name` instead of recomputing the `p.I` expression, `use_GF_suffix=False` removes the `GF` suffix, and `enable_simd=True` wraps the access as `ReadSIMD(&...)`.

### Rank-N Registration

Rank helpers such as `register_gridfunctions_for_single_rank1()`, `register_gridfunctions_for_single_rank2()`, and `register_gridfunctions_for_single_rankN()` use `nrpy.indexedexp` declaration functions to build nested SymPy component arrays. For rank greater than one, the optional symmetry string is passed to the indexed-expression declaration. The nested component list is flattened and deduplicated before registration, so symmetric aliases such as `gDD01` and `gDD10` register only once. Each registered component receives a component-specific description of the form `<base_desc>_<component>`, is registered with `is_basename=False`, and carries rank metadata.

For BHaH generated structs, `register_griddata_commondata()` stores `GridCommonData` declarations under `par.glb_extras_dict["griddata_struct"]` or `par.glb_extras_dict["commondata_struct"]`, depending on whether the data is grid-local or common across grids.

Generation-time parallel codegen snapshots and merges the parameter, code-parameter, C-function, Python-function, gridfunction, and extras dictionaries. Use [Parallel Codegen Orchestration](helpers/parallel-codegen-orchestration.md) for that multiprocessing merge contract; this page owns the objects and registries being merged.

Split into separate parameter and gridfunction leaves is deferred unless this leaf becomes hard to scan.

## Sources

- [nrpy/params.py](../../nrpy/params.py) - `NRPyParameter`, `CodeParameter`, `_parse_array_spec`, `register_param`, `parval_from_str`, `set_parval_from_str`, `register_CodeParameters`, `register_CodeParameter`, `adjust_CodeParam_default`, `glb_params_dict`, `glb_code_params_dict`, `glb_extras_dict`
- [nrpy/grid.py](../../nrpy/grid.py) - `GridFunction`, `BHaHGridFunction`, `ETLegacyGridFunction`, `CarpetXGridFunction`, `gridfunction_lists`, `get_parity_type`, `set_parity_types`, `define_gfs_group`, `gridfunction_defines`, `register_gridfunctions`, `register_gridfunctions_for_single_rankN`, `glb_gridfcs_dict`, `GF_CLASS_MAP`
- [nrpy/infrastructures/BHaH/griddata_commondata.py](../../nrpy/infrastructures/BHaH/griddata_commondata.py) - `GridCommonData`, `register_griddata_commondata`

## See Also

- Parent: [Core APIs](index.md)
- Depends on: [Indexed Expressions](indexed-expressions.md)
- Depends on: [Symbolic Expression Utilities](helpers/symbolic-expression-utilities.md)
- See also: [C Codegen](c-codegen.md)
- See also: [Finite Difference](finite-difference.md)
- See also: [Reference Metrics](reference-metrics.md)
- See also: [Parallel Codegen Orchestration](helpers/parallel-codegen-orchestration.md)
