# MoL Time Integration

> BHaH route for Method of Lines parameters, Runge-Kutta routing, stage storage, and time-step hooks. Status: confirmed.
> Up: [BHaH](index.md)

## Summary

BHaH Method of Lines support registers timestep CodeParameters, allocates
integrator-specific stage arrays, emits a `MoL_gridfunctions_struct`, and
generates `MoL_step_forward_in_time()`. The stepping function dispatches through
explicit Runge-Kutta Butcher tables, calls user-provided RHS and post-RHS hook
strings, updates gridfunctions, advances `commondata->time`, and increments
`commondata->nn`.

## Detail

`MoLtimestepping/register_all.py` is the one-stop registration entrypoint. It
registers commondata CodeParameters `nn_0`, `nn`, `CFL_FACTOR`, `dt`, `t_0`,
`time`, and `t_final`; checks that the selected parallelization is supported;
builds the Butcher-table dictionary; registers stage malloc/free functions;
optionally registers `MoL_step_forward_in_time()`; adds
`MoL_gridfunctions_struct gridfuncs` to `griddata_struct`; and registers MoL's
`BHaH_defines.h` contribution.

`rk_butcher_table_dictionary.py` owns the available explicit integrators. Its
`generate_Butcher_tables()` dictionary includes Euler, several RK2/RK3 methods,
SSP methods, classic `RK4`, Dormand-Prince and other high-order tables, adaptive
embedded tables, and optional Adams-Bashforth generation. `validate()` checks
selected tables against simple ODEs. `is_diagonal_Butcher()` distinguishes
diagonal storage/update paths. `intermediate_stage_gf_names_list()` determines
which temporary gridfunction arrays are needed: non-diagonal methods allocate
`next_y_input_gfs` and `k*_gfs` arrays, diagonal RK3 methods use
`k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs` and
`k2_or_y_nplus_a32_k2_gfs`, and diagonal methods such as `RK4` use
`y_nplus1_running_total_gfs` plus alternating `k_odd_gfs` and `k_even_gfs`
except for Euler.

`BHaH_defines.py` emits `MoL_gridfunctions_struct`. It always contains
`REAL *y_n_gfs`, appends method-specific intermediate stage pointers, and ends
with `REAL *auxevol_gfs`. The malloc routine computes
`Nxx_plus_2NGHOSTS_tot` from `params`, then allocates each intermediate stage
array as `sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot`, using
`BHAH_MALLOC_DEVICE` for CUDA and `BHAH_MALLOC` otherwise. The free routine
inspects the same method-specific list and frees those arrays with
`BHAH_FREE_DEVICE` or `BHAH_FREE`.

`MoL_step_forward_in_time.py` generates the full step driver. It rejects adaptive
Butcher tables because adaptive RK is not implemented in BHaH MoL. For each
grid, it creates aliases for `y_n_gfs`, method-specific stage arrays,
`auxevol_gfs`, `params`, and either `rfmstruct` or `xx[3]`; if curvilinear BCs
are enabled, it also aliases `bcstruct`. It chooses algorithm-specific substep
construction for diagonal RK3, non-diagonal RK, Euler, or other diagonal RK
methods. Each substep is generated through `single_RK_substep_input_symbolic()`
with a time offset from the Butcher table.

The hook strings are textual C fragments with placeholder substitution.
`rhs_string` is evaluated with `RK_INPUT_GFS` and `RK_OUTPUT_GFS` replaced by
the substep input and output arrays. `post_rhs_string` is applied after each RK
update, usually for boundary conditions on the output array.
`post_post_rhs_string` runs after the per-grid loop with the same output
replacement, allowing global work after all grids complete a substep. The
registration flags `enable_rfm_precompute`, `enable_curviBCs`, and
`enable_intrinsics` determine whether generated aliases include `rfmstruct`,
`bcstruct`, and SIMD/CUDA intrinsic headers and update loops.

`RK3 Ralston` and `RK3 Heun` require `rhs_accumulation_string` in addition
to the ordinary RHS hook. Its placeholders `RK_INPUT_GFS`, `RK_STAGE_GFS`,
`RK_ACCUM_GFS`, `RK_STAGE_DT`, and `RK_ACCUM_DT` specify an unchanged input
stage and two distinct arrays that receive weighted derivatives. The first
substep saves the weighted first derivative before projection changes its
stage. The second RHS adds directly to `y_n_gfs` and the final running sum.
The third RHS supplies the last derivative, and the final result remains in
the original `y_n_gfs` allocation. Both methods retain their Butcher
coefficients, use three evolved-field arrays including `y_n_gfs`, and evaluate
the RHS three times. The post-RHS hook must prepare or project only the next
input stage or final state, and never the running sum.

Claim evidence:
- Claim: Diagonal RK3 requires the accumulation hook above, preserves the original Ralston and Heun coefficients, and uses three evolved arrays and three RHS evaluations when its hooks follow the stated input and projection requirements.
- Role: descriptive behavior
- Deciding authority: [MoL_step_forward_in_time.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/MoL_step_forward_in_time.py), `register_CFunction_MoL_step_forward_in_time`; [rk_butcher_table_dictionary.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/rk_butcher_table_dictionary.py), `generate_Butcher_tables` and `intermediate_stage_gf_names_list`.
- Corroboration: [runtime.c](../../../.github/tests/rk3_cpu/runtime.c), `dense_step`, the linked evaluation counters, and the allocation checks in `main`; [generate.py](../../../.github/tests/rk3_cpu/generate.py), `generate`, builds the complete CPU drivers. The configured comparison covers BSSN and does not establish fCCZ4, CUDA, distributed execution, or long-evolution convergence.

`general_relativity/rhs_eval.py::MoL_rhs_strings()` builds matching ordinary
and accumulated BSSN/fCCZ4 calls.
`register_CFunction_rhs_eval(enable_RK3_accumulation=True)` registers the
accumulated RHS, and
`CurviBoundaryConditions.register_all.register_C_functions(enable_RK3_accumulation=True)`
registers its radiation boundary function. With three distinct arrays, the
accumulated kernel reads the projected stage, computes one point's derivatives,
and immediately adds them to the two output arrays. Its partial SIMD group
stores only interior lanes, and an output with a zero coefficient remains
unchanged. Ricci and other pre-RHS calculations must leave both live output
arrays unchanged.

Claim evidence:
- Claim: The separate RHS and boundary registration flags enable the matching accumulated calls; the RHS reads a fixed stage, adds only to interior output lanes, and leaves an output with a zero multiplier unchanged.
- Role: descriptive behavior
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `MoL_rhs_strings` and `register_CFunction_rhs_eval`; [boundary register_all.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/register_all.py), `register_C_functions`.
- Corroboration: [runtime.c](../../../.github/tests/rk3_cpu/runtime.c), `check_accumulation`, compares independent initial outputs with ordinary RHS updates, checks input and zero-weight identity, and includes partial SIMD rows. [check.py](../../../.github/tests/rk3_cpu/check.py), `main`, configures scalar and AVX2 CPU execution with the generated BSSN functions; this comparison does not execute every coordinate system or formulation.

For radiation, preparation subtracts each old boundary-source contribution
in reverse dependency order; completion adds the updated source and radiation
term in forward order. With distinct arrays and a fixed input stage, this
preserves the original ghost values plus their weighted RHS increments
without another complete evolved array.

Claim evidence:
- Claim: Accumulated radiation and inner boundaries preserve arbitrary old output values plus the weighted boundary RHS by reversing source subtraction during preparation and restoring sources and radiation terms in forward order during completion.
- Role: descriptive behavior
- Deciding authority: [apply_bcs_outerradiation_and_inner.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerradiation_and_inner.py), `register_CFunction_apply_bcs_outerradiation_and_inner_accumulate`, `setup_Cfunction_apply_bcs_inner_accumulate`, and `setup_Cfunction_apply_bcs_pure_accumulate`.
- Corroboration: [runtime.c](../../../.github/tests/rk3_cpu/runtime.c), `make_boundary` and `check_accumulation`, compare the ordinary generated radiation RHS with the accumulated result across faces, edges, corners, and inner sources. The boundary owner's doctests compare complete emitted CPU and CUDA source; source comparison alone does not establish GPU execution.

`rk_substep.py` provides the low-level RK code generator. `RKFunction` turns
symbolic RK left-hand and right-hand sides into a launchable update function,
using OpenMP or CUDA loop generation and optional SIMD/CUDA intrinsic reads and
writes. `MoL_Functions_dict` stores generated substep functions so
`construct_RK_functions_prefunc()` can prepend them to the final C function.
`single_RK_substep_input_symbolic()` builds the substep order: set
`commondata->time` to `time_start + c_i * dt`, update CUDA constant params when
needed, run RHS, launch the RK update, run post-RHS hooks, close the grid loop,
then run post-post-RHS hooks. In the diagonal RK3 k2 substep, the post-RHS hook
targets the next stage state, not the running final-state accumulator; the
accumulator is only combined pointwise into the final `y_n_gfs`.

After all substeps, `MoL_step_forward_in_time()` updates
`commondata->time = commondata->t_0 + (commondata->nn - commondata->nn_0 + 1) *
commondata->dt` to reduce roundoff accumulation across many steps and across
regrids where `dt` can change. It then increments `commondata->nn`.

## Sources

- [register_all.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/register_all.py) - `register_CFunctions`
- [MoL_step_forward_in_time.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/MoL_step_forward_in_time.py) - `register_CFunction_MoL_step_forward_in_time`
- [MoL_malloc_intermediate_stage_gfs.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/MoL_malloc_intermediate_stage_gfs.py) - `register_CFunction_MoL_malloc_intermediate_stage_gfs`
- [MoL_free_intermediate_stage_gfs.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/MoL_free_intermediate_stage_gfs.py) - `register_CFunction_MoL_free_intermediate_levels`
- [apply_bcs_outerradiation_and_inner.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerradiation_and_inner.py) - radiation and inner accumulation functions
- [boundary register_all.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/register_all.py) - `register_C_functions`
- [check.py](../../../.github/tests/rk3_cpu/check.py) - scoped generated CPU comparison
- [generate.py](../../../.github/tests/rk3_cpu/generate.py) - `generate`
- [runtime.c](../../../.github/tests/rk3_cpu/runtime.c) - `dense_step`, `check_accumulation`, linked evaluation counters
- [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py) - `MoL_rhs_strings`, `register_CFunction_rhs_eval`
- [rk_substep.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/rk_substep.py) - `RKFunction`, `single_RK_substep_input_symbolic`, `construct_RK_functions_prefunc`, `MoL_Functions_dict`, `check_supported_parallelization`
- [rk_butcher_table_dictionary.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/rk_butcher_table_dictionary.py) - `generate_Butcher_tables`, `validate`, `is_diagonal_Butcher`, `intermediate_stage_gf_names_list`
- [BHaH_defines.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/BHaH_defines.py) - `register_BHaH_defines_h`

## See Also

- [BHaH](index.md)
- [Grids, Coordinates, Reference Metrics, And Boundaries](grids-coordinates-reference-metrics-and-boundaries.md)
- [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
