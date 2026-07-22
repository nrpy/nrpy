# NRPyElliptic Workflow

> BHaH workflow for conformally flat NRPyElliptic hyperbolic relaxation. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [BHaH](index.md)

## Summary

The standalone BHaH NRPyElliptic example assembles a hyperbolic-relaxation
solver around the standard grid, reference-metric, Method-of-Lines, boundary,
diagnostics, checkpointing, and generated-project build paths. Its generated
runtime initializes `uu` and `vv`, fills auxiliary-evolution background/source
gridfunctions plus a local variable wavespeed, evolves the relaxation system,
computes a Hamiltonian-residual diagnostic, and stops when either `nn_max` or
the configured log residual tolerance is reached.

## Detail

`register_CFunction_initial_guess_single_point` emits the pointwise initial
guess and currently sets both `uu_ID` and `vv_ID` to zero.
`register_CFunction_initial_guess_all_points` registers the generated
`initial_data` function: it optionally returns early from checkpoint input, then
loops over every grid and every point, reading the BHaH coordinate arrays and
writing `UUGF` and `VVGF` in `y_n_gfs`.

Auxiliary-evolution setup is separated from the initial guess.
`register_CFunction_auxevol_gfs_set_to_constant` registers
`MINIMUM_GLOBAL_WAVESPEED`, launches `variable_wavespeed_gfs_all_points`, and
then launches `auxevol_gfs_all_points`. The wavespeed kernel uses the selected
reference metric's orthogonal scale factors and local `dxx*` values to store
`VARIABLE_WAVESPEEDGF`. The auxevol kernel calls
`compute_psi_background_and_ADD_times_AUU` through
`generate_prefunc_auxevol_gfs_single_point` and stores `PSI_BACKGROUNDGF` and
`ADD_TIMES_AUUGF`. In the standalone example, this hook runs after MoL storage
allocation through `post_non_y_n_auxevol_mallocs`; CUDA builds copy the
auxiliary data back to the host after synchronizing the per-grid stream.

`register_CFunction_rhs_eval` is the evolved-system RHS bridge. It constructs
`HyperbolicRelaxationCurvilinearRHSs` with reference-metric precompute enabled,
generates C for `uu_rhs` and `vv_rhs`, and emits an interior-point BHaH loop
whose signature consumes `commondata`, `params`, `rfmstruct`, `auxevol_gfs`,
current `in_gfs`, and destination `rhs_gfs`. The example inserts that RHS into
the MoL registration and then applies radiation boundary conditions with a
custom per-grid wavespeed read from `VARIABLE_WAVESPEEDGF`; extrapolation
boundaries are applied in the post-RHS hook.

Residual evaluation uses the same symbolic owner but a separate generated
diagnostic function. `register_CFunction_residual_H_compute_all_points`
registers `residual_H_compute_all_points` under `diagnostics/`, loops over the
interior, and writes the `residual` expression into a caller-supplied
destination buffer. `register_CFunction_diagnostic_gfs_set` registers DIAG
gridfunctions for `DIAG_RESIDUAL`, `DIAG_UU`, `DIAG_VV`, and
`DIAG_GRIDINDEX`; it calls the residual helper, optionally applies inner BCs to
the residual diagnostic for interpolation users, copies `uu` and `vv` from the
current MoL state, and records the grid index at every point.

Nearest diagnostics and volume diagnostics consume those diagnostic buffers.
`register_CFunction_diagnostics_nearest` emits a dispatcher that samples
`DIAG_RESIDUALGF` and `DIAG_UUGF` at the grid center, nearest y/z axes, and
nearest xy/yz planes through the common BHaH nearest-output helpers.
`register_CFunction_diagnostics_volume_integration` builds integration recipes:
one over the whole domain and one inside `sphere_R_80`, both using the squared
`DIAG_RESIDUALGF` integrand. After executing the recipes, it extracts the RMS
for `sphere_R_80` and updates `commondata->log10_current_residual`.

`register_CFunction_stop_conditions_check` owns the relaxation exit policy. It
registers `stop_relaxation`, `nn_max`, `log10_residual_tolerance`, and
`log10_current_residual` in `commondata`, then sets `stop_relaxation` when the
iteration count reaches `nn_max` or when the current log residual drops below
the tolerance. The standalone example invokes this hook after each MoL step and
forces a checkpoint before breaking out of the main evolution loop.

The generated-project relation is ordinary BHaH assembly. The example registers
BHaH numerical grids and timestep with reference-metric precompute and CurviBCs,
generates `xx_to_Cart`, registers all diagnostics, registers rfm-precompute
wrappers, writes the `set_CodeParameters*.h` parameter-access headers, emits the
parfile/parser, writes `BHaH_defines.h`, emits `main`, registers griddata cleanup,
and constructs the Makefile. CI covers the generated elliptic project as part of
the generated project build matrix; generated `project/` outputs remain products,
not KB sources.

## Sources

- [initial_data.py](../../../nrpy/infrastructures/BHaH/nrpyelliptic/initial_data.py) - `register_CFunction_initial_guess_single_point`, `register_CFunction_initial_guess_all_points`
- [auxevol_gfs_set_to_constant.py](../../../nrpy/infrastructures/BHaH/nrpyelliptic/auxevol_gfs_set_to_constant.py) - `register_CFunction_auxevol_gfs_set_to_constant`, `generate_prefunc_variable_wavespeed_gfs_all_points`, `generate_prefunc_auxevol_gfs_single_point`
- [rhs_eval.py](../../../nrpy/infrastructures/BHaH/nrpyelliptic/rhs_eval.py) - `register_CFunction_rhs_eval`
- [residual_H_compute_all_points.py](../../../nrpy/infrastructures/BHaH/nrpyelliptic/residual_H_compute_all_points.py) - `register_CFunction_residual_H_compute_all_points`
- [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/nrpyelliptic/diagnostic_gfs_set.py) - `register_CFunction_diagnostic_gfs_set`
- [diagnostics_nearest.py](../../../nrpy/infrastructures/BHaH/nrpyelliptic/diagnostics_nearest.py) - `register_CFunction_diagnostics_nearest`
- [diagnostics_volume_integration.py](../../../nrpy/infrastructures/BHaH/nrpyelliptic/diagnostics_volume_integration.py) - `register_CFunction_diagnostics_volume_integration`
- [stop_conditions_check.py](../../../nrpy/infrastructures/BHaH/nrpyelliptic/stop_conditions_check.py) - `register_CFunction_stop_conditions_check`
- [nrpyelliptic_conformally_flat.py](../../../nrpy/examples/nrpyelliptic_conformally_flat.py) - `rhs_string`, `post_non_y_n_auxevol_mallocs`, `post_MoL_step_forward_in_time`

## See Also

- [BHaH](index.md)
- [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
- [Conformally Flat Elliptic](../../equations/conformally-flat-elliptic.md)
- [Generated Project CI](../../validation/generated-project-ci.md)
