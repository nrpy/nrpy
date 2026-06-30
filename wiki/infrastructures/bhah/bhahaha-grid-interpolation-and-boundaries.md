# BHaHAHA Grid, Interpolation, And Boundaries

> Explain how BHaHAHA builds its spherical grids, moves metric data through interpolation layers, and applies local boundary metadata. · Status: confirmed · Last reconciled: 2026-06-30
> Up: [BHaH](index.md)

## Summary

BHaHAHA has two metric-ingest paths that meet at the same public
single-horizon solver contract. In a generated BHaH evolution,
`BHaHAHA_interpolate_metric_data_nrpy` samples BSSN gridfunctions from the
evolution grid, transforms each sampled point to Cartesian-basis BSSN data, and
converts that point to the ADM Cartesian components expected by BHaHAHA. At the
public C boundary, callers instead provide Cartesian ADM `gamma_{ij}` and
`K_{ij}` components directly in the `NUM_EXT_INPUT_CARTESIAN_GFS` order.

Inside the solver, the input data is reorganized onto three spherical grids:
the external input grid, the interpolation-source grid, and the evolved surface
grid. The external grid stores Cartesian ADM input converted to rescaled
spherical BSSN-like fields. The interpolation-source grid shares the external
radial sampling but uses the current evolved angular resolution and adds
derivatives. The evolved grid is a two-angular-dimensional spherical surface
grid with one radial index, used for `h(theta, phi)` and surface-local metric
interpolation.

## Detail

The public Cartesian metric order is defined in `BHaHAHA_header.h`: six
independent `gammaDD` components followed by six independent `KDD` components,
with names from `INTERP_GAMMADDXXGF` through `INTERP_KDDZZGF`. The generated
BHaH adapter does not use that enum while sampling BSSN data. It first uses the
private `INTERP_BSSN_GF_INDICES` order for `aDD`, `cf`, `hDD`, and `trK`, then
stores final ADM output in the private `FINAL_ADM_METRIC_INDICES` order. Both
private orders are implementation details of `BHaHAHA_interpolate_metric_data_nrpy`;
the stable public input layout remains the Cartesian ADM order in the header.

`BHaHAHA_interpolate_metric_data_nrpy` builds a spherical destination grid
around the current horizon center with cell-centered theta and phi and the
radial array already chosen for the horizon search. Each destination point is
converted to the source reference-metric coordinates with
`Cart_to_xx_and_nearest_i0i1i2`, then
`interpolation_3d_general__uniform_src_grid` samples the selected BSSN
gridfunctions from the generated BHaH evolution grid. For each sampled point,
`basis_transform_BSSN_rfm_to_Cartesian_single_point` moves native-basis BSSN
data to Cartesian basis, and `BHaHAHA_BSSN_to_ADM_Cartesian` writes ADM
`gammaDD` and `KDD` into the BHaHAHA input buffer.

`bah_radial_grid_cell_centered_set_up` chooses the radial sampling used for
the external input data. It floors `input_r_min` at `0`, caps `input_r_max` at
`max_search_radius`, enforces at least four interior radial points, adds
ghost-zone coverage, returns `Nr_external_input`, `r_min_external_input`, and
`dr_external_input`, and fills the `radii` array. If the requested inner radius
is positive but the implied ghost-zone grid would cross zero, setup falls back
to a zero-origin radial grid.

`bah_numgrid__external_input_set_up` allocates `external_input_gfs`, copies the
caller-owned no-ghost Cartesian input into ghost-zone storage, and creates
uniform cell-centered coordinate arrays in `(r, theta, phi)`. Its angular
resolution is the finest configured multigrid level; its radial interior count
depends on whether `r_min_external_input` is zero, because a zero lower radius
shifts active radial data by `NGHOSTS`. The setup then converts Cartesian ADM
input to spherical rescaled fields: `external_spherical_WW`, `external_spherical_trK`,
`external_spherical_hDD*`, and `external_spherical_aDD*`.

`bah_numgrid__interp_src_set_up` builds the 3D interpolation-source grid. It
reuses the external radial grid, takes the current evolved-grid angular counts,
allocates `interp_src_gfs`, and creates matching `(r, theta, phi)` coordinates.
`bah_interpolation_2d_external_input_to_interp_src_grid` performs angular
Lagrange interpolation from external input to this grid at every overlapping
radial point. The setup then copies the interpolated external spherical fields
into the `src_*` gridfunctions used by later radial interpolation.

Derivative setup happens on the interpolation-source grid, not on the evolved
surface. `bah_hDD_dD_and_W_dD_in_interp_src_grid_interior` computes angular
derivatives `partial_D_hDD` and `partial_D_WW` over active radial points and
radial derivatives over the radial interior. Afterward,
`bah_apply_bcs_r_maxmin_partial_r_hDD_upwinding` fills radial derivative data
at the outer radial ghost layers with finite-difference stencils selected by
boundary offset. It also fills inner radial ghost layers when
`r_min_external_input != 0`; zero-origin grids rely on the inner-boundary
mapping path instead.

`bah_numgrid__evol_set_up` creates the evolved spherical grid for horizon
relaxation. It registers `hh` and `vv` as evolved fields, registers metric and
derivative fields as auxiliary evolved fields, sets `NUMGRIDS = 1`, sets
`Nxx0` from the caller-provided evolved grid shape, and uses spherical
cell-centered coordinates with theta in `[0, pi]` and phi in `[-pi, pi]`. This
grid stores the BHaHAHA-local `bcstruct` in `griddata`, because evolved-grid
inner boundary application uses the same parity metadata as other local grids.

The 1D interpolation layer maps metric data from the interpolation-source grid
onto the current evolved surface. `bah_interpolation_1d_radial_spokes_on_3d_src_grid`
interpolates along radial spokes only, at destination radial index `NGHOSTS`,
using `hh` as the destination radius field. It returns distinct error codes for
null pointers, too-large interpolation order, horizons beyond the input grid,
and horizons below `r_min_external_input`. `bah_interpolation_2d_general__uniform_src_grid`
is the separate angular interpolation helper for arbitrary theta/phi
destinations, used by other surface-transfer paths such as coarse-to-fine
surface data.

The interpolation helpers expose their own guard surfaces. The 3D BHaH-grid
interpolator checks coordinate/output pointers, interpolation order versus each
source dimension, and out-of-bounds destination stencils. The external-to-source
2D interpolator checks required grid pointers, angular interpolation order, and
theta/phi stencil bounds. The general 2D and radial 1D helpers expose matching
null-pointer, order, and out-of-bounds checks through BHaHAHA error codes. The
grid setup functions separately return allocation errors for external and
interpolation-source gridfunction and coordinate arrays.

Boundary support is BHaHAHA-local. `bah_bcstruct_set_up` emits the BHaHAHA
`innerpt_bc_struct`, `outerpt_bc_struct`, `bc_info_struct`, and `bc_struct`
definitions, plus local origin-free coordinate helpers. Inner boundary points
map through eigencoordinates to an in-bounds point and carry both ten base
parity signs and a signed-permutation derivative Jacobian. If the parity values
or derivative Jacobian cannot be represented by those signs, setup returns a
boundary error instead of applying an invalid boundary condition.

Gridfunction parity arrays are generated from local gridfunction names and
ranks. Scalars get scalar parity, vectors get component parity, rank-2 tensor
components get tensor parity, and stored derivative gridfunctions on the
interpolation-source grid carry base-field parity plus derivative-direction
metadata. `bah_apply_bcs_inner_only` applies evolved-grid inner boundary
conditions only at the active radial index `NGHOSTS`, while external input and
interpolation-source setup apply inner boundary copies directly from their
local `bcstruct` arrays. `bcstruct_set_up.py` rejects `GeneralRFM` coordinate
systems for BHaHAHA boundary generation at code-generation time.

## Sources

- [BHaHAHA_header.h](../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h) - `NUM_EXT_INPUT_CARTESIAN_GFS`, `INTERP_GAMMADDXXGF`, `INTERP_KDDZZGF`
- [BHaH_implementation.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaH_implementation.py) - `FINAL_ADM_METRIC_INDICES`, `INTERP_BSSN_GF_INDICES`, `bhahaha_gf_interp_indices`, `BHaHAHA_BSSN_to_ADM_Cartesian`, `BHaHAHA_interpolate_metric_data_nrpy`
- [radial_grid_cell_centered_set_up.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/radial_grid_cell_centered_set_up.py) - `register_CFunction_radial_grid_cell_centered_set_up`, `radial_grid_cell_centered_set_up`
- [numgrid__external_input_set_up.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/numgrid__external_input_set_up.py) - `register_CFunction_numgrid__external_input_set_up`, `numgrid__external_input_set_up`
- [numgrid__interp_src_set_up.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/numgrid__interp_src_set_up.py) - `register_CFunction_numgrid__interp_src_set_up`, `numgrid__interp_src_set_up`
- [numgrid__evol_set_up.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/numgrid__evol_set_up.py) - `register_CFunction_numgrid__evol_set_up`, `numgrid__evol_set_up`
- [interpolation_3d_general__uniform_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_3d_general__uniform_src_grid.py) - `register_CFunction_interpolation_3d_general__uniform_src_grid`, `interpolation_3d_general__uniform_src_grid`
- [interpolation_2d_external_input_to_interp_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_2d_external_input_to_interp_src_grid.py) - `register_CFunction_interpolation_2d_external_input_to_interp_src_grid`, `interpolation_2d_external_input_to_interp_src_grid`
- [interpolation_2d_general__uniform_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_2d_general__uniform_src_grid.py) - `register_CFunction_interpolation_2d_general__uniform_src_grid`, `interpolation_2d_general__uniform_src_grid`
- [interpolation_1d_radial_spokes_on_3d_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_1d_radial_spokes_on_3d_src_grid.py) - `register_CFunction_interpolation_1d_radial_spokes_on_3d_src_grid`, `interpolation_1d_radial_spokes_on_3d_src_grid`
- [hDD_dD_and_W_dD_in_interp_src_grid_interior.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/hDD_dD_and_W_dD_in_interp_src_grid_interior.py) - `register_CFunction_hDD_dD_and_W_dD_in_interp_src_grid_interior`, `hDD_dD_and_W_dD_in_interp_src_grid_interior`
- [bcstruct_set_up.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/bcstruct_set_up.py) - `register_CFunction_bcstruct_set_up`, `BHaH_defines_set_gridfunction_defines_with_parity_types`, `parity_conditions_symbolic_dot_products`
- [apply_bcs_inner_only.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/apply_bcs_inner_only.py) - `register_CFunction_apply_bcs_inner_only`, `apply_bcs_inner_only`
- [apply_bcs_r_maxmin_partial_r_hDD_upwinding.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/apply_bcs_r_maxmin_partial_r_hDD_upwinding.py) - `register_CFunction_apply_bcs_r_maxmin_partial_r_hDD_upwinding`, `setup_Cfunction_FD1_arbitrary_upwind`

## See Also

- [BHaH](index.md)
- [BHaHAHA Horizon Runtime](bhahaha-horizon-runtime.md)
- [GR Application Wiring](gr-application-wiring.md)
- [Reference Metrics](../../core/reference-metrics.md)
