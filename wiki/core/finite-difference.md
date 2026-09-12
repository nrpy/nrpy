# Finite Difference

> Core route for finite-difference operators in generated C kernels. · Status: confirmed
> Up: [Core APIs](index.md)

## Summary

`nrpy.finite_difference` supplies the finite-difference support used by `c_codegen()` when expressions contain derivative symbols. It computes finite-difference coefficients and stencils, extracts derivative variables from SymPy expressions, maps derivative names back to base gridfunctions and operators, emits memory reads for stencil points, and can package derivative formulas as reusable C helper functions through `FDFunction`.

## Detail

Finite-difference order is an NRPy parameter registered as `finite_difference::fd_order`, with default value `4`. `compute_fdcoeffs_fdstencl()` takes a derivative-operator string and order, chooses `stencil_width = fd_order + 1`, chooses an `UPDOWNWIND_stencil_shift`, fetches the inverse finite-difference matrix, and returns coefficient and stencil lists. Matrix inversion is cached by `FD_Matrix_dict`, a `Matrix_dict` keyed by `(stencil_width, UPDOWNWIND_stencil_shift)`. The cache computes missing entries lazily through `setup_FD_matrix__return_inverse_lowlevel()`, whose low-level inverse uses SymPy LU decomposition and multiplies `U_inv * L_inv`.

Operator support is string-driven. Centered first derivatives use the first-derivative row; unmixed second derivatives use the second-derivative row when the final two directions match; mixed second derivatives use products of first-derivative rows when the two final directions differ. `dKOD` bumps `fd_order` by two and uses the last matrix row with Kreiss-Oliger scaling. `dupD` and `ddnD` shift by one grid point, while `dfullupD` and `dfulldnD` shift by half the finite-difference order. Strings containing `DDD` are rejected because this implementation supports derivatives only through second order.

The codegen pipeline starts from expression free symbols. `symbol_is_gridfunction_Cparameter_or_other()` classifies a symbol as `gridfunction`, `Cparameter`, or `other` by checking `gri.glb_gridfcs_dict` and `par.glb_code_params_dict`; it raises if the same name is registered both ways. `extract_list_of_deriv_var_strings_from_sympyexpr_list()` treats registered gridfunctions and registered C parameters as known quantities, collects derivative-like `other` symbols containing `_dD`, `_dKOD`, `_dupD`, or `_ddnD`, removes duplicates, and returns a deterministic SymPy sort.

The helper's current `_ddnD` expansion condition is broader than the parameter name suggests: whenever `upwind_control_vec` is not a Python `str`, each `_dupD` also adds the corresponding `_ddnD`. Thus `sp.Symbol("unset")` triggers expansion, while the literal string `"unset"` does not. `c_codegen()` currently passes `sp.Symbol("unset")` to this extraction helper independently of the caller's `CCodeGen.upwind_control_vec`; only a later list-valued control vector activates final upwind selection.

`extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars()` parses derivative names into base gridfunction names and derivative operators. It enforces the naming contract that the number of final numeric suffix digits equals the total number of `U` and `D` characters in the symbol name. It then finds the last underscore, counts derivative `D` markers after it, and takes the differentiated gridfunction's rank as the remaining trailing digits. That subtraction, rather than a count of rank markers adjacent to the underscore, is what lets a gridfunction whose own name records a derivative be differentiated again, as `hDDdD_dD0011` is. It splits names such as `hDD_dDD0112` into the base gridfunction component and operator suffix expected by the finite-difference routines.

`select_stored_first_derivatives()` returns aligned source-gridfunction and operator lists without changing derivative symbols. It is called by FD codegen before prototype construction and memory-read planning, only when the consumer supplies a nonempty `stored_first_derivatives` selection. For selected, registered storage, `hDD_dD011` keeps its temporary name but reads `hDDdD011`; `hDD_dDD0112` keeps its temporary name but applies `dD2` to `hDDdD011`. An empty operator marks a pointwise read internally. Only centered first derivatives and canonical mixed pairs (`01`, `02`, `12`) select storage. Unregistered selections raise `ValueError`; registered but unselected fields have no effect. Equation symmetry zeros remain zeros, and the caller owns storage validity, array availability, and halos.

Claim evidence:
- Claim: Explicit per-kernel storage selection changes derivative evaluation sources while preserving mathematical temporary names; an empty selection retains ordinary lowering.
- Role: descriptive behavior
- Deciding authority: [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `select_stored_first_derivatives` and its selection/error doctests.
- Corroboration: [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `CCodeGen`, `gridfunction_management_and_FD_codegen`, and the polynomial evaluation doctest; [BHaH GR Application Wiring](../infrastructures/bhah/gr-application-wiring.md) records consumer integration and bounded evolution validation.

`fd_temp_variable_name()` constructs temporary names for stencil-point reads. Offsets become suffixes such as `i0m2`, `i1p4`, or `i2m1`, joined after the gridfunction basename; zero offsets contribute no suffix, so the center point keeps the bare gridfunction name. An empty basename is rejected.

`read_gfs_from_memory()` plans memory reads for the base gridfunctions needed by derivative stencils and for zeroth-derivative gridfunctions that appear directly in the expressions. It records offsets as `i0,i1,i2` triples, removes duplicate reads with `superfast_uniq()`, and sorts each gridfunction's accesses by a synthetic memory address matching `mem_alloc_style` `210` or `012`. It raises a runtime error if an expression contains derivatives but no gridfunctions are registered. In SIMD mode, supported infrastructures substitute `REAL_SIMD_ARRAY` for the gridfunction C type and emit SIMD-aware one-point reads.

Prototype operators are built through `proto_FD_operators_to_sympy_expressions()`. It creates `FDPROTO_*` derivative symbols, `FDPROTO` stencil-point symbols, and `invdxx0`, `invdxx1`, and `invdxx2` metric-spacing factors. Each prototype expression is preprocessed with `cse_preprocess(prefix="FDPart1", factor=True)` so rational constants become symbols such as `FDPart1_Rational_1_4`; SIMD mode also asks preprocessing to declare negative-one symbols. The function populates `FDFunctions_dict` by operator, storing the prototype formula and its rational-symbol map for later helper generation.

`FDFunction` represents a generated finite-difference helper function. It records the floating-point alias, finite-difference order, operator, rational constants, symbolic formula, SIMD mode, generated C function name, and static helper `CFunction`. For BHaH operators `dDD01`, `dDD02`, and `dDD12`, `FDFunction.c_function_call()` passes a pointer to the current point in the registered field's array (including AUXEVOL and custom arrays), the required non-unit strides, and both inverse spacings. Direction 0 has implicit unit stride; directions 1 and 2 use `Nxx_plus_2NGHOSTS0` and `Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1`. `FDFunction.CFunction_fd_function()` delegates these helpers to `FDFunction.pointer_stride_mixed_params_body()`, which loads its own stencil using scalar indexing or `ReadSIMD()` and evaluates a tensor product of centered first-derivative lines. For BHaH CUDA generation, `dKOD0`, `dKOD1`, and `dKOD2` also use a pointer, the direction's non-unit stride if needed, and its inverse spacing. These KO helpers load the existing stencil locally and retain the supplied expression code unchanged. Selection reads the existing `parallelization` setting; CPU KO helpers keep their value arguments. Other operators and infrastructures retain calls using the prototype's non-rational free symbols as parameters, with `FDPROTO` replaced by the concrete gridfunction name; their helpers return `FD_result` from the supplied expression code. Upwind/downwind derivative outputs retain the `UpwindAlgInput` prefix. CarpetX helpers retain host/device decorators.

Claim evidence:
- Claim: BHaH `dDD01`, `dDD02`, and `dDD12` helpers accept a field pointer, non-unit strides, and inverse spacings and own their scalar/SIMD stencil loads. BHaH CUDA KO helpers likewise own their stencil loads while preserving the existing KO expression; CPU KO and other operators/infrastructures retain the prototype-value interface.
- Role: descriptive behavior
- Deciding authority: [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `use_pointer_stride_mixed_derivative`, `use_pointer_stride_ko_derivative`, `FDFunction.c_function_call`, `FDFunction.CFunction_fd_function`, `FDFunction.pointer_stride_mixed_params_body`
- Corroboration: [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `c_codegen` scalar/SIMD mixed-helper doctest examples, including custom-array and ordinary-read cases (inspected, not executed).

`construct_FD_functions_prefunc()` walks `FDFunctions_dict`, applies any requested C-function decorators, regenerates each helper's full function text, and concatenates those functions for use in a generated C function `prefunc`.

Finite-difference pages own derivative naming, stencil construction, memory-read ordering, prototype FD operators, and FD helper-function lifecycle. Helper CSE and SIMD pages own shared expression preprocessing, deterministic CSE handling, and symbolic SIMD rewrite mechanics that can be used by the generated finite-difference code path.

## Sources

- [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `compute_fdcoeffs_fdstencl`, `extract_list_of_deriv_var_strings_from_sympyexpr_list`, `extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars`
- [nrpy/finite_difference.py](../../nrpy/finite_difference.py) - `read_gfs_from_memory`, `proto_FD_operators_to_sympy_expressions`, `FDFunction`, `construct_FD_functions_prefunc`
- [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `gridfunction_management_and_FD_codegen`

## See Also

- Parent: [Core APIs](index.md)
- See also: [C Codegen](c-codegen.md)
- Depends on: [Gridfunctions And Parameters](gridfunctions-and-parameters.md)
- Depends on: [CSE And Printer Support](helpers/cse-and-printer-support.md)
- Depends on: [SIMD And Intrinsic Support](helpers/simd-and-intrinsic-support.md)
