# BHaHAHA Public API And Input Contract

> Define the BHaHAHA public C boundary, caller-owned inputs, generated-app controls, poisoning checks, and error vocabulary. · Status: confirmed · Last reconciled: 06-30-2026
> Up: [BHaH](index.md)

## Summary

`BHaHAHA_header.h` is the public C boundary for BHaHAHA. It defines the input
and diagnostics structs, the Cartesian ADM input ordering, the horizon-shape
indexing macro, the multiresolution storage limit, and the public `bah_*`
entry points for poisoning, radial-grid setup, center/radius extrapolation,
horizon finding, and diagnostic-file output.

The generated BHaH application has a separate input layer. Its
`CodeParameter` registrations, `commondata` arrays, BBH controls, and
multigrid checks configure how BHaH drives BHaHAHA, but they are not public C
API. The generated app adapts BHaH BSSN data into the public Cartesian ADM
input contract before calling `bah_find_horizon`; the runtime flow is covered
in [BHaHAHA Horizon Runtime](bhahaha-horizon-runtime.md).

## Detail

The public metric input is `input_metric_data` inside
`bhahaha_params_and_data_struct`. It is caller-owned storage for
`NUM_EXT_INPUT_CARTESIAN_GFS == 12` Cartesian ADM components on a flattened
`Nr x Ntheta x Nphi x gf` layout: six independent `gamma_{ij}` components
followed by six independent `K_{ij}` components. The public enum in
`BHaHAHA_header.h` names the order as `INTERP_GAMMADDXXGF`,
`INTERP_GAMMADDXYGF`, `INTERP_GAMMADDXZGF`, `INTERP_GAMMADDYYGF`,
`INTERP_GAMMADDYZGF`, `INTERP_GAMMADDZZGF`, then the matching `INTERP_KDD*`
slots. `IDX2(itheta, iphi)` is the public 2D horizon-shape indexing macro,
with theta as the fast physical coordinate in the argument list and phi
varying in the outer loop. `MAX_RESOLUTIONS == 16` fixes public storage for
multigrid angular-resolution arrays.

`bhahaha_params_and_data_struct` groups the caller-visible input contract into
metric/radial-grid data, external iteration/time metadata, multigrid angular
resolution, hyperbolic-relaxation controls, diagnostics settings, and previous
horizon history. Callers set the public input values before `bah_find_horizon`.
The previous-shape pointers `prev_horizon_m1`, `prev_horizon_m2`, and
`prev_horizon_m3` must point to caller-owned storage sized for the finest
`max(Ntheta) x max(Nphi)` surface. BHaHAHA owns the history values it writes
there; the public caller owns allocation lifetime, while the generated BHaH
orchestrator allocates and frees those arrays for its own commondata storage.
The struct comments label `t_m1/t_m2/t_m3`, previous min/max radii, and
previous centers as persistent quantities set by BHaHAHA.

`bhahaha_diagnostics_struct` is the public output container. It carries
convergence norms, area, centroid coordinates, centroid-relative coordinate
radii, coordinate-plane proper circumferences, spin estimates derived from
circumference ratios, general spin-axis circumference fields, and
`Theta_eval_points_counter`.

The public prototypes in `BHaHAHA_header.h` are:

- `bah_poisoning_set_inputs(bhahaha_params_and_data_struct *params)`
- `bah_poisoning_check_inputs(const bhahaha_params_and_data_struct *params)`
- `bah_radial_grid_cell_centered_set_up(...)`
- `bah_xyz_center_r_minmax(...)`
- `bah_find_horizon(bhahaha_params_and_data_struct *params, bhahaha_diagnostics_struct *diags)`
- `bah_diagnostics_file_output(...)`

The poisoning helpers are defensive caller-side checks. Public
`bah_poisoning_set_inputs` rejects a null params pointer, then poisons selected
required inputs by setting `REAL` scalars to `NaN`, pointers to `NULL`, scalar
ints and int arrays to `-100`, and the multigrid arrays to the same int
sentinel. Public `bah_poisoning_check_inputs` rejects a null params pointer,
checks the selected fields for those poisoned values, prints field-specific
messages, and exits if any selected input remains poisoned. The selected fields
include `input_metric_data`, external time/iteration, radial-grid metadata,
`num_resolutions_multigrid`, active multigrid entries, fixed-radius flag, CFL,
mass scale, eta, KO strength, max iterations, `Theta_Linf_times_M_tolerance`,
horizon index/count, verbosity, and the eta-varying flag. The poisoning helpers
do not check every public struct field; notably, BHaHAHA writes persistent
previous-horizon history values into caller-owned storage, and
`Theta_L2_times_M_tolerance` is not covered by the current poisoning set/check
functions.

Error-code vocabulary is generated from `error_code_msg_tuples_list` in
`error_message.py`. That list is the source of truth for symbolic names and
user-facing messages, including `BHAHAHA_SUCCESS`, `bah_find_horizon` timing,
iteration-limit, and small-horizon failures, boundary-structure failures,
initial-data and numerical-grid allocation failures, interpolation null-pointer
and interpolation-order failures, horizon out-of-bounds failures, and
diagnostic-circumference allocation failure. `register_CFunction_error_message`
registers the C helper that maps a `bhahaha_error_codes` value to its message
string and prints a fallback error if the code is unknown.

Generated BHaH application controls are registered in
`register_CFunction_bhahaha_find_horizons`. The generated app creates
`commondata->bhahaha_params_and_data[max_horizons]` and
`commondata->bhahaha_diagnostics[max_horizons]`, then registers BHaHAHA
`CodeParameter` inputs. Per-horizon `REAL[max_horizons]` controls are
`bah_max_search_radius`, `bah_cfl_factor`,
`bah_Theta_Linf_times_M_tolerance`, `bah_Theta_L2_times_M_tolerance`,
`bah_eta_damping_times_M`, `bah_M_scale`, `bah_KO_strength`, and initial
grid-center arrays `bah_initial_grid_x_center`,
`bah_initial_grid_y_center`, and `bah_initial_grid_z_center`.

Generated scalar int controls are `bah_Nr_interp_max`,
`bah_verbosity_level`, `bah_max_iterations`, `bah_max_num_horizons`,
`bah_num_resolutions_multigrid`, `bah_enable_BBH_mode`,
`bah_BBH_mode_common_horizon_idx`, and
`bah_enable_eta_varying_alg_for_precision_common_horizon`. Generated int-array
controls are `bah_Ntheta_array_multigrid`, `bah_Nphi_array_multigrid`, and
`bah_BBH_mode_inspiral_BH_idxs`; `bah_BBH_mode_horizon_active` is generated
commondata state and is not added to the parfile. These generated controls
feed BHaH orchestration and should not be described as public C ABI.

The generated app checks that `bah_num_resolutions_multigrid` is positive and
that active `bah_Ntheta_array_multigrid` and `bah_Nphi_array_multigrid`
entries are positive before running; the current check does not enforce an
upper bound against the generated array length or `MAX_RESOLUTIONS`.
BBH mode is also generated-app policy: when enabled, it requires exactly three
horizons, requires valid non-overlapping individual/common horizon indices,
starts the two individual horizons active, and starts the common horizon
inactive. These controls shape BHaH's multi-horizon orchestration around the
public single-horizon `bah_find_horizon` call.

## Sources

- [BHaHAHA_header.h](../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h) - `bhahaha_params_and_data_struct`, `bhahaha_diagnostics_struct`, `NUM_EXT_INPUT_CARTESIAN_GFS`, `IDX2`, `MAX_RESOLUTIONS`, public `bah_*` prototypes
- [BHaH_implementation.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaH_implementation.py) - `register_CFunction_bhahaha_find_horizons`, `check_multigrid_resolution_inputs`, `initialize_bhahaha_solver_params_and_shapes`, generated `CodeParameter` registrations and BHaH `commondata` arrays
- [poisoning_set_inputs.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/poisoning_set_inputs.py) - `register_CFunction_poisoning_set_inputs`, `poisoning_set_inputs`
- [poisoning_check_inputs.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/poisoning_check_inputs.py) - `register_CFunction_poisoning_check_inputs`, `poisoning_check_inputs`
- [error_message.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/error_message.py) - `error_code_msg_tuples_list`, `register_CFunction_error_message`

## See Also

- [BHaH](index.md)
- [BHaHAHA Horizon Runtime](bhahaha-horizon-runtime.md)
- [Runtime Data, Parameters, Headers, And CLI](runtime-data-parameters-headers-and-cli.md)
- [Diagnostics Output And Checkpointing](diagnostics-output-and-checkpointing.md)
