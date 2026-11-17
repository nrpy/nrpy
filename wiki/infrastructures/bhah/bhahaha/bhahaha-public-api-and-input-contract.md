# BHaHAHA Public C Header And Input Fields

> Define the BHaHAHA public C header, caller-allocated ADM inputs, generated BHaH controls, poisoning checks, and error codes. · Status: confirmed
> Up: [BHaHAHA](index.md)

## Summary

BHaHAHA integrates apparent-horizon finding in two layers: `bhahaha_find_horizons` orchestrates multi-horizon history, generated BHaH controls, and BSSN-to-ADM sampling, while `bah_find_horizon` (backed by `find_horizon`) solves a single horizon via multiresolution hyperbolic relaxation through the public C interface in `BHaHAHA_header.h`.

## Detail

### Public C Structs, ADM Ordering, And Entry Points

`BHaHAHA_header.h` defines the public C interface:
- `input_metric_data` inside `bhahaha_params_and_data_struct` points to caller-allocated memory for `NUM_EXT_INPUT_CARTESIAN_GFS == 12` Cartesian ADM components in `Nr x Ntheta x Nphi x gf` layout: `INTERP_GAMMADDXXGF`..`INTERP_GAMMADDZZGF` followed by `INTERP_KDDXXGF`..`INTERP_KDDZZGF`.
- `IDX2(itheta, iphi)` indexes 2D horizon shapes with `theta` fast and `phi` slow; `MAX_RESOLUTIONS == 16` bounds multigrid angular resolution arrays.
- `bhahaha_params_and_data_struct` groups radial-grid inputs, external iteration/time, multigrid angular counts, relaxation parameters, diagnostic flags, and previous-horizon arrays (`prev_horizon_m1`, `prev_horizon_m2`, `prev_horizon_m3` sized to `max(Ntheta) x max(Nphi)`). The caller manages buffer allocation (handled by `bhahaha_find_horizons` in generated BHaH apps), while BHaHAHA updates `t_m1/m2/m3`, previous radii, centers, and shape history.
- `bhahaha_diagnostics_struct` outputs `Theta` norms, horizon area, centroid, centroid-relative radii, `xy`/`xz`/`yz` and general spin-axis proper circumferences, circumference-ratio spin magnitudes, and `Theta_eval_points_counter`.
- Public prototypes: `bah_poisoning_set_inputs`, `bah_poisoning_check_inputs`, `bah_radial_grid_cell_centered_set_up`, `bah_xyz_center_r_minmax`, `bah_find_horizon`, and `bah_diagnostics_file_output`.

### Poisoning Helpers And Error Vocabulary

`bah_poisoning_set_inputs` rejects a null `params` pointer and initializes selected required fields to `NaN` (`REAL`), `NULL` (pointers), or `-100` (ints and multigrid arrays). `bah_poisoning_check_inputs` verifies those fields were populated before solving and exits with field-specific diagnostics if any sentinel remains; persistent history fields written by BHaHAHA and `Theta_L2_times_M_tolerance` are not poisoned.

`error_code_msg_tuples_list` in `error_message.py` defines `bhahaha_error_codes` (`BHAHAHA_SUCCESS`, iteration/horizon-size limits, boundary parity errors, grid/circumference allocation failures, and interpolation pointer/order/out-of-bounds codes), emitted via `register_CFunction_error_message`.

### Generated BHaH Controls And Multigrid Parity

`register_CFunction_bhahaha_find_horizons` allocates `commondata->bhahaha_params_and_data[max_horizons]` and `commondata->bhahaha_diagnostics[max_horizons]` and registers BHaH `CodeParameter` controls (distinct from the public C header):
- Per-horizon `REAL[max_horizons]`: `bah_max_search_radius`, `bah_cfl_factor`, `bah_Theta_Linf_times_M_tolerance`, `bah_Theta_L2_times_M_tolerance`, `bah_eta_damping_times_M`, `bah_M_scale`, `bah_KO_strength`, and `bah_initial_grid_{x,y,z}_center`.
- Scalar `int`: `bah_Nr_interp_max`, `bah_verbosity_level`, `bah_max_iterations`, `bah_max_num_horizons`, `bah_num_resolutions_multigrid`, `bah_enable_BBH_mode`, `bah_BBH_mode_common_horizon_idx`, and `bah_enable_eta_varying_alg_for_precision_common_horizon`.
- `int` arrays: `bah_Ntheta_array_multigrid`, `bah_Nphi_array_multigrid`, `bah_BBH_mode_inspiral_BH_idxs`, plus non-parfile `commondata` array `bah_BBH_mode_horizon_active` (BBH mode requires `3` horizons with disjoint inspiral/common indices, starting individual horizons active and the common horizon inactive).

`check_multigrid_resolution_inputs` checks positivity of `bah_num_resolutions_multigrid` and active `(Ntheta, Nphi)` entries without checking `MAX_RESOLUTIONS` or parity. Because BHaHAHA spherical grids span `theta in [0, pi]` and `phi in [-pi, pi]`, every active `bah_Nphi_array_multigrid` entry must be even per [angular inner-boundary requirements](../curvilinear-boundaries-and-stencils.md#hard-requirements-for-angular-inner-boundaries), whereas `Ntheta` may be odd or even (defaults are even).

Evidence (scientific/interface requirement): `BHaHAHA/numgrid__evol_set_up.py::register_CFunction_numgrid__evol_set_up` and `BHaHAHA/bcstruct_set_up.py`; corroboration: `BHaH_implementation.py::check_multigrid_resolution_inputs` (odd `Nphi` is unguarded at runtime).

## Sources

- [BHaHAHA_header.h](../../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h) - `bhahaha_params_and_data_struct`, `bhahaha_diagnostics_struct`, `NUM_EXT_INPUT_CARTESIAN_GFS`, `IDX2`, `MAX_RESOLUTIONS`, `bah_*` prototypes
- [BHaH_implementation.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaH_implementation.py) - `register_CFunction_bhahaha_find_horizons`, `bhahaha_find_horizons`, `check_multigrid_resolution_inputs`, `initialize_bhahaha_solver_params_and_shapes`
- [find_horizon.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/find_horizon.py) - `register_CFunction_find_horizon`, `find_horizon`
- [numgrid__evol_set_up.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/numgrid__evol_set_up.py) - `register_CFunction_numgrid__evol_set_up`
- [bcstruct_set_up.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/bcstruct_set_up.py) - `Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt`
- [poisoning_set_inputs.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/poisoning_set_inputs.py) - `register_CFunction_poisoning_set_inputs`
- [poisoning_check_inputs.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/poisoning_check_inputs.py) - `register_CFunction_poisoning_check_inputs`
- [error_message.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/error_message.py) - `error_code_msg_tuples_list`, `register_CFunction_error_message`

## See Also

- Depends on: [Curvilinear Boundaries And Stencils](../curvilinear-boundaries-and-stencils.md)
- Depends on: [Runtime Data, Parameters, Headers, And CLI](../runtime-data-parameters-headers-and-cli.md)
- Example: [Apparent Horizon Library](../../../examples/apparent-horizon-library.md)
