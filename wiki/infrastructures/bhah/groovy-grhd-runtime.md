# GRoovy GRHD Runtime

> Explain how BHaH GRoovy registers and orchestrates generated GRHD runtime kernels. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [BHaH](index.md)

## Summary

GRoovy is the BHaH runtime leaf for general relativistic hydrodynamics. It
registers the metric, primitive, conserved, reconstructed, face-centered, flux,
temperature, entropy, and optional neutrino gridfunctions needed by generated
GRHD applications, then emits C kernels that connect NRPy's symbolic GRHD
expressions to GRHayL primitive recovery, EOS calls, reconstruction, and
interface flux routines.

This page covers orchestration and generated runtime plumbing. The Valencia
equations, characteristic speeds, and HLL expressions themselves belong to the
equation-side GRHD page; GRoovy consumes them to generate BHaH kernels.
GRHayL and EOS algorithms are an external dependency boundary: GRoovy passes
`ghl_parameters`, `ghl_eos_parameters`, `ghl_primitive_quantities`,
`ghl_conservative_quantities`, and `ghl_metric_quantities` to `ghl_*` calls but
does not define GRHayL recovery or EOS internals.

## Detail

`register_all_grhd_gridfunctions` is the usual first GRoovy registration call.
When `evolving_spacetime=True`, it imports the full BSSN set through
`BSSN_quantities[CoordSystem]`; otherwise it registers fixed-metric lapse,
conformal factor, trace curvature, conformal metric, traceless curvature, and
shift in `auxevol_gfs`. It always registers face-centered metric fields
`alpha_face`, `cf_face`, `h_faceDD`, and `vet_faceU`, velocity fields
`u4Ut`, reconstructed `u4rUt/u4lUt`, `rescaledvU`, and reconstructed
`rescaledvrU/rescaledvlU`.

The evolved hydrodynamic variables are `rho_star`, `tau_tilde`, and
`rescaledstildeD`. GRoovy also registers HLL flux storage
`rho_star_HLL_flux`, `tau_tilde_HLL_flux`, and
`rescaledStilde_flux_HLLD`, primitive storage `rhob`, `P`, and their right/left
reconstructed states. Temperature mode adds `Ye`, `temperature`, their
right/left states, `Ye_star`, and `Ye_star_HLL_flux`; entropy mode adds `S`,
right/left entropy states, `S_star`, and `S_star_HLL_flux`; neutrino mode adds
the NRPyLeakage optical-depth and opacity fields.

`grhd_rhs_eval` is the top-level RHS orchestrator for
`partial_t U = -partial_i F^i + S`. It first calls
`calculate_all_source_terms` to initialize RHS gridfunctions with source and
connection terms. It then loops over `flux_dir = 0, 1, 2`, selects the matching
`calculate_HLL_fluxes_direction_i` function pointer, interpolates metric data
to cell faces, reconstructs primitive variables, optionally reconstructs
entropy, computes right/left `u^0` on faces with
`compute_up_index_velocity_time_component_pointwise`, computes HLL fluxes, and
adds flux divergences for that direction.

`calculate_all_source_terms` instantiates equation-side `GRHD_Equations`,
calls `construct_all_equations`, and emits the geometric source plus
reference-metric connection contributions into the RHS arrays. It initializes a
GRHayL primitive state from the current primitive gridfunctions and asks
`ghl_compute_h_and_cs2` for the enthalpy and sound speed used by the generated
symbolic source expressions. Optional `Ye_star` and `S_star` connection terms
are included when the corresponding variables are evolved.

`interpolate_metric_gfs_to_cell_faces` computes third-order face values from a
four-point stencil with coefficients `AM2`, `AM1`, `A0`, and `A1`. It covers
the BSSN metric variables needed at interfaces: `hDD`, `vetU`, `alpha`, and
`cf`. If spacetime is fixed, the same generated loop reads metric data from
`auxevol_gfs` instead of the evolved gridfunction array.

`reconstruction_loop` builds the primitive right/left states consumed by the
Riemann solver. `_scheme_settings` supports `wenoz` and `mc`, and
`_build_reconstruction_body` emits the stencil loop over the selected flux
direction. Temperature reconstruction reconstructs `rho`, velocity, `Ye`, and
`temperature`, then uses tabulated EOS calls to enforce bounds and compute
right/left pressure. Entropy mode can also register
`reconstruction_entropy_loop`, which reconstructs entropy and enforces
tabulated `rho/Ye/S` bounds on the interface states.

`calculate_HLL_flux_dirn_i` registers one kernel per flux direction:
`calculate_HLL_fluxes_direction_0`, `calculate_HLL_fluxes_direction_1`, and
`calculate_HLL_fluxes_direction_2`. Each kernel initializes right and left
GRHayL primitive states from reconstructed interface gridfunctions, calls
`ghl_compute_h_and_cs2` for each side, builds Cartesian symbolic face data,
uses equation-side `calculate_HLL_fluxes`, and stores HLL fluxes for density,
energy, momentum, and optional `Ye`/entropy channels.

`calculate_flux_divergences` adds the flux-difference contribution back to the
RHS. It computes the reference-metric rescaling factors and their derivatives,
uses equation-side `GRHD_Equations` to build the rescaled physical flux
expressions needed for the divergence formula, then combines those with
neighboring HLL flux gridfunctions at `index` and `indexp1`. Momentum RHSs are
rescaled back into the evolved `rescaledstildeD` basis.

The primitive/conservative conversion path is split by direction of data
movement. `primitives_to_conservatives_routine` initializes GRHayL primitive
states from `rhob`, `P`, velocity, and optional temperature/entropy fields,
updates tabulated temperature when needed, obtains `h` and `cs2`, and emits
equation-side conserved variables into `rho_star`, `tau_tilde`,
`rescaledstildeD`, and optional `Ye_star`/`S_star`. The inverse
`conservatives_to_primitives_routine` selects one of `Hybrid`,
`HybridEntropy`, `Tabulated`, or `TabulatedEntropy`, rebuilds GRHayL
conservative and metric structs via basis-transform helpers, calls the
matching GRHayL recovery routine, falls back to weighted neighbor averages when
needed, and stores recovered primitives back in the reference-metric basis. In
robust tabulated-entropy mode it tries Palenzuela and Newman energy/entropy
recoveries, compares recomputed-conservative mismatch after primitive limiting,
and falls back to atmosphere values when all methods fail.

`basis_transform_rfm_basis_to_Cartesian` and
`basis_transform_rfm_basis_to_Cartesian__read_cons_only` are the inbound GRHayL
bridge: they read NRPy gridfunctions in the current reference-metric basis,
transform velocities and conserved momentum into Cartesian components, and
initialize GRHayL metric/conservative/primitive structs. The
`basis_transform_Cartesian_to_rfm_basis` bridge writes GRHayL primitive and
conservative data back into the NRPy reference-metric basis. These helpers mark
the boundary between generated NRPy storage conventions and GRHayL's Cartesian
struct API.

`compute_up_index_velocity_time_component_pointwise` enforces the configured
maximum Lorentz factor through branchless min/max helpers, rescales velocity
when it exceeds the speed limit, and writes both corrected `rescaledvU` and
`u4Ut`. `adjust_hydrodynamic_data` re-runs GRHayL primitive limiting and stores
updated `u^0` after initial metric data are available. `compute_entropy`
recomputes entropy from tabulated temperature data or from hybrid EOS density
and pressure, depending on the selected thermodynamics mode.

`apply_copy_and_outflow_bcs` owns hydrodynamic primitive boundary conditions.
Outer ghost zones copy scalar primitive data from the nearest interior point,
zero optional leakage variables, project the velocity against the outward
radial unit vector to remove inflow, and recompute `u^0`. Inner boundary points
are filled from the `bcstruct` parity map using `auxevol_gf_parity` for
primitive and optional leakage fields.

`hybrid_EoS_TOV_initial_data` registers BHaH `initial_data` for a hybrid-EOS
TOV run. It executes caller-provided GRHayL setup code, solves the TOV model
with `TOVola_solve`, converts spherical ADM data to BSSN through
`initial_data_reader__convert_ADM_Spherical_to_BSSN`, interpolates pointwise
TOV data through `TOVola_interp`, recovers pressure and baryon density with
`ghl_hybrid_compute_rho_cold_from_P_cold`, applies primitive limiting, stores
the data back in the rfm basis, writes `u4Ut`, and frees the TOV persistent
arrays.

## Sources

- [register_all_grhd_gridfunctions.py](../../../nrpy/infrastructures/BHaH/GRoovy/register_all_grhd_gridfunctions.py) - `register_all_grhd_gridfunctions`
- [grhd_rhs_eval.py](../../../nrpy/infrastructures/BHaH/GRoovy/grhd_rhs_eval.py) - `register_CFunction_grhd_rhs_eval`, `grhd_rhs_eval`
- [calculate_all_source_terms.py](../../../nrpy/infrastructures/BHaH/GRoovy/calculate_all_source_terms.py) - `register_CFunction_calculate_all_source_terms`, `calculate_all_source_terms`
- [interpolate_metric_gfs_to_cell_faces.py](../../../nrpy/infrastructures/BHaH/GRoovy/interpolate_metric_gfs_to_cell_faces.py) - `register_CFunction_interpolate_metric_gfs_to_cell_faces`, `interpolate_metric_gfs_to_cell_faces`
- [reconstruction_loop.py](../../../nrpy/infrastructures/BHaH/GRoovy/reconstruction_loop.py) - `_scheme_settings`, `_build_reconstruction_body`, `reconstruction_loop`, `reconstruction_entropy_loop`
- [calculate_HLL_flux_dirn_i.py](../../../nrpy/infrastructures/BHaH/GRoovy/calculate_HLL_flux_dirn_i.py) - `register_CFunction_calculate_HLL_flux_dirn_i`, `calculate_HLL_fluxes`
- [calculate_flux_divergences.py](../../../nrpy/infrastructures/BHaH/GRoovy/calculate_flux_divergences.py) - `register_CFunction_calculate_flux_divergences`, `calculate_flux_divergences`
- [primitives_to_conservatives_routine.py](../../../nrpy/infrastructures/BHaH/GRoovy/primitives_to_conservatives_routine.py) - `register_CFunction_primitives_to_conservatives_routine`, `primitives_to_conservatives_routine`
- [conservatives_to_primitives_routine.py](../../../nrpy/infrastructures/BHaH/GRoovy/conservatives_to_primitives_routine.py) - `_select_recovery_mode`, `_build_tabulated_entropy_robust_block`, `register_CFunction_conservatives_to_primitives_routine`, `conservatives_to_primitives_routine`
- [basis_transform_rfm_basis_to_Cartesian.py](../../../nrpy/infrastructures/BHaH/GRoovy/basis_transform_rfm_basis_to_Cartesian.py) - `register_CFunction_basis_transform_rfm_basis_to_Cartesian`, `basis_transform_rfm_basis_to_Cartesian`
- [basis_transform_rfm_basis_to_Cartesian__read_cons_only.py](../../../nrpy/infrastructures/BHaH/GRoovy/basis_transform_rfm_basis_to_Cartesian__read_cons_only.py) - `register_CFunction_basis_transform_rfm_basis_to_Cartesian__read_cons_only`, `basis_transform_rfm_basis_to_Cartesian__read_cons_only`
- [basis_transform_Cartesian_to_rfm_basis.py](../../../nrpy/infrastructures/BHaH/GRoovy/basis_transform_Cartesian_to_rfm_basis.py) - `register_CFunction_basis_transform_Cartesian_to_rfm_basis`, `basis_transform_Cartesian_to_rfm_basis`
- [compute_up_index_velocity_time_component_pointwise.py](../../../nrpy/infrastructures/BHaH/GRoovy/compute_up_index_velocity_time_component_pointwise.py) - `register_CFunction_compute_up_index_velocity_time_component_pointwise`, `compute_up_index_velocity_time_component_pointwise`
- [adjust_hydrodynamic_data.py](../../../nrpy/infrastructures/BHaH/GRoovy/adjust_hydrodynamic_data.py) - `register_CFunction_adjust_hydrodynamic_data`, `adjust_hydrodynamic_data`
- [compute_entropy.py](../../../nrpy/infrastructures/BHaH/GRoovy/compute_entropy.py) - `register_CFunction_compute_entropy`, `compute_entropy`
- [apply_copy_and_outflow_bcs.py](../../../nrpy/infrastructures/BHaH/GRoovy/apply_copy_and_outflow_bcs.py) - `register_CFunction_apply_copy_and_outflow_bcs`, `apply_copy_and_outflow_bcs`
- [hybrid_EoS_TOV_initial_data.py](../../../nrpy/infrastructures/BHaH/GRoovy/hybrid_EoS_TOV_initial_data.py) - `register_CFunction_hybrid_EoS_TOV_initial_data`, `initial_data`
- [GRHD_equations.py](../../../nrpy/equations/grhd/GRHD_equations.py) - `GRHD_Equations`, `construct_all_equations`
- [HLL_fluxes.py](../../../nrpy/equations/grhd/HLL_fluxes.py) - `calculate_HLL_fluxes`

## See Also

- [BHaH](index.md)
- [GR Application Wiring](gr-application-wiring.md)
- [GRHD](../../equations/grhd.md)
- [Metric Conversions And Matter](../../equations/general-relativity/metric-conversions-and-matter.md)
