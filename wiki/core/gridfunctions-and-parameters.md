# Gridfunctions And Parameters

> Core route for symbolic gridfunctions, NRPy parameters, code parameters, and generated data structs. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Core APIs](index.md)

## Summary

Gridfunctions and parameters are the main bridge between symbolic expressions and generated infrastructure data. `nrpy.grid` registers symbolic gridfunction names and infrastructure-specific metadata. `nrpy.params` stores Python-side `NRPyParameter` values and C-side `CodeParameter` symbols. BHaH grid/common data registration records declarations that later become members of generated `griddata_struct` and `commondata_struct` support code.

## Detail

`GridFunction` is the base object for registered gridfunctions. It stores the name, group, description, rank, dimension, C type, asymptotic value, wave speed, and whether the name is a basename. Base names must be non-empty strings and may not end in an integer, because rank components such as `betU2` rely on numeric suffixes being unambiguous.

Infrastructure subclasses choose storage details. `BHaHGridFunction` accepts groups `EVOL`, `AUXEVOL`, `DIAG`, and `AUX`; it maps default array names such as `in_gfs`, `auxevol_gfs`, `diagnostic_gfs`, and `aux_gfs`, and emits `IDX4(<NAME>GF, i0, i1, i2)` access strings. ETLegacy and CarpetX have their own subclasses selected by the active `Infrastructure` parameter.

`register_gridfunctions()` selects the correct infrastructure class, stores each gridfunction in `glb_gridfcs_dict`, and returns SymPy symbols. Rank helpers such as `register_gridfunctions_for_single_rank1()`, `register_gridfunctions_for_single_rank2()`, and `register_gridfunctions_for_single_rankN()` use `nrpy.indexedexp` declarations to build component symbols and register only unique components after symmetry is applied.

`NRPyParameter` stores Python-level configuration in `glb_params_dict`. `CodeParameter` stores runtime C parameters in `glb_code_params_dict`, creates a SymPy symbol, validates defaults, handles `#define` and array cases, and controls whether the parameter belongs in parfiles or generated `set_CodeParameters*.h` support. The registration helpers `register_CodeParameter()` and `register_CodeParameters()` return symbolic parameters for use in expression construction.

For BHaH generated structs, `register_griddata_commondata()` stores `GridCommonData` declarations under `par.glb_extras_dict["griddata_struct"]` or `par.glb_extras_dict["commondata_struct"]`, depending on whether the data is grid-local or common across grids.

Generation-time parallel codegen snapshots and merges the parameter, code-parameter, C-function, Python-function, gridfunction, and extras dictionaries. Use [Parallel Codegen Orchestration](helpers/parallel-codegen-orchestration.md) for that multiprocessing merge contract; this page owns the objects and registries being merged.

## Sources

- [nrpy/grid.py](../../nrpy/grid.py) - `GridFunction`, `BHaHGridFunction`, `register_gridfunctions`, `register_gridfunctions_for_single_rankN`
- [nrpy/params.py](../../nrpy/params.py) - `NRPyParameter`, `CodeParameter`, `register_CodeParameter`, `register_CodeParameters`
- [nrpy/infrastructures/BHaH/griddata_commondata.py](../../nrpy/infrastructures/BHaH/griddata_commondata.py) - `GridCommonData`, `register_griddata_commondata`
- [coding_style.md](../../coding_style.md) - gridfunction naming and code-parameter registration rules

## See Also

- [Core APIs](index.md)
- [Indexed Expressions](indexed-expressions.md)
- [Finite Difference](finite-difference.md)
- [Symbolic Expression Utilities](helpers/symbolic-expression-utilities.md)
- [Parallel Codegen Orchestration](helpers/parallel-codegen-orchestration.md)
