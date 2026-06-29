# Finite Difference

> Core route for finite-difference operators in generated C kernels. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Core APIs](index.md)

## Summary

`nrpy.finite_difference` supplies the finite-difference support used by `c_codegen()` when expressions contain derivative symbols. It computes finite-difference coefficients and stencils, extracts derivative variables from SymPy expressions, maps derivative names back to base gridfunctions and operators, emits memory reads for stencil points, and can package derivative formulas as reusable C helper functions through `FDFunction`.

## Detail

Finite-difference order is an NRPy parameter registered as `finite_difference::fd_order`, with default value `4`. `compute_fdcoeffs_fdstencl()` takes a derivative-operator string and order, chooses a stencil width and upwind/downwind shift, fetches the cached inverse finite-difference matrix, and returns coefficient and stencil lists. It supports first derivatives, unmixed and mixed second derivatives, Kreiss-Oliger derivatives, and upwind/downwind variants such as `dupD`, `ddnD`, `dfullupD`, and `dfulldnD`; third and higher derivative strings containing `DDD` are rejected.

The codegen pipeline starts from expression free symbols. `extract_list_of_deriv_var_strings_from_sympyexpr_list()` treats registered gridfunctions and registered C parameters as known quantities; derivative-like symbols with `_dD`, `_dKOD`, `_dupD`, or `_ddnD` are collected as derivative variables. When an upwind control vector is supplied, each `_dupD` also adds the corresponding `_ddnD` variable needed by the upwind algorithm.

`extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars()` parses derivative names into base gridfunction names and derivative operators, enforcing the convention that the number of final numeric indices matches the number of `U` and `D` markers. `read_gfs_from_memory()` then builds C declarations for the unique gridfunction stencil points needed by the derivatives and by any zeroth-derivative gridfunctions appearing directly in the expressions. It removes duplicate reads and sorts accesses according to the requested memory allocation style.

`FDFunction` represents a generated finite-difference helper function. It records the floating-point alias, order, operator, rational constants, symbolic formula, SIMD mode, and generated C function name. `construct_FD_functions_prefunc()` emits the registered helper functions for use in a C function `prefunc`.

Finite-difference pages own derivative naming, stencil construction, and memory-read ordering. Helper CSE and SIMD pages own shared expression preprocessing, deterministic CSE handling, and symbolic SIMD rewrite mechanics that can be used by the generated finite-difference code path.

## Sources

- [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `compute_fdcoeffs_fdstencl`, `extract_list_of_deriv_var_strings_from_sympyexpr_list`, `extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars`
- [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `read_gfs_from_memory`, `FDFunction`, `construct_FD_functions_prefunc`
- [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `gridfunction_management_and_FD_codegen`
- [coding_style.md](../../coding_style.md) - derivative naming conventions

## See Also

- [Core APIs](index.md)
- [C Codegen](c-codegen.md)
- [Gridfunctions And Parameters](gridfunctions-and-parameters.md)
- [CSE And Printer Support](helpers/cse-and-printer-support.md)
- [SIMD And Intrinsic Support](helpers/simd-and-intrinsic-support.md)
