# MoL Time Integration

> BHaH route for Method of Lines parameters, Runge-Kutta routing, stage storage, and time-step hooks. Status: confirmed. Last reconciled: 06-29-2026
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
- [rk_substep.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/rk_substep.py) - `RKFunction`, `single_RK_substep_input_symbolic`, `construct_RK_functions_prefunc`, `MoL_Functions_dict`, `check_supported_parallelization`
- [rk_butcher_table_dictionary.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/rk_butcher_table_dictionary.py) - `generate_Butcher_tables`, `validate`, `is_diagonal_Butcher`, `intermediate_stage_gf_names_list`
- [BHaH_defines.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/BHaH_defines.py) - `register_BHaH_defines_h`

## See Also

- [BHaH](index.md)
- [Grids, Coordinates, Reference Metrics, And Boundaries](grids-coordinates-reference-metrics-and-boundaries.md)
- [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
