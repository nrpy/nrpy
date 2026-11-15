# Curvilinear Boundaries And Stencils

> BHaH eigencoordinate ghost-zone classification, parity-signed inner and radiation/extrapolation outer boundary kernels, angular grid-count invariants, and Lagrange stencil ghost-zone reads. · Status: confirmed
> Up: [BHaH](index.md)

## Summary

`CurviBoundaryConditions` classifies every BHaH ghost point via an eigencoordinate round-trip into `pure_outer_bc_array` or parity-signed `inner_bc_array`, then fills pure outer faces (radiation or second-order inward extrapolation) before copying inner ghost zones with 10 orthonormal-basis parity signs. Angular coordinates must span their `ReferenceMetric` defining bounds with even `Nphi` wherever `phi -> phi + pi` applies (`SymTP`/`SinhSymTP` also require even `xx1`, and `RingHoleySinhSpherical` requires `Ntheta >= 2 NGHOSTS`), enabling Lagrange interpolation stencils to cross angular seams and poles directly through filled inner ghost points.

## Detail

### Boundary Structs, Eigencoordinate Setup, And Application Kernels

`CurviBoundaryConditions/register_all.py::register_C_functions` registers `outer_bc_type` (default `"radiation"`), `bcstruct_set_up()`, `apply_bcs_outerradiation_and_inner()`, `apply_bcs_outerextrap_and_inner()` (plus specific-gridfunction variant), `apply_bcs_inner_only()`, and CPU-only `apply_bcs_inner_only_specific_gfs()`. `CurviBoundaryConditions/BHaH_defines.py` emits `innerpt_bc_struct`, `outerpt_bc_struct`, `bc_info_struct`, and `bc_struct` in `griddata_struct` (`inner_bc_array`, `pure_outer_bc_array[NGHOSTS*3]`, `bc_info`), parity tables inferred from registered gridfunction names (types `0` scalar, `1-3` vector, `4-9` symmetric rank-2 tensor; placeholder parity on high-rank `GeneralRFM` `auxevol` metric fields), and CUDA `cpyHosttoDevice_bc_struct()`.

`bcstruct_set_up()` maps each ghost point's logical `xx` through its eigencoordinate (`Spherical`, `Cylindrical`, `SymTP`, `Cartesian` including wedges, or `GeneralRFM` numerical inverse) `xx_to_Cart` -> `Cart_to_xx` -> nearest index. Points mapping to themselves enter `pure_outer_bc_array` by ghost layer and face direction; all others enter `inner_bc_array` (`dst` ghost point, `src` interior or pure-outer point, plus 10 parity signs computed as dot products of dst/src orthonormal basis vectors). `bcstruct_set_up()` aborts if dst and src Cartesian positions or unit-vector dot products (`+/-1`) exceed tolerance; optional masks skip unset points (CUDA masking is rejected).

Boundary kernels update stored gridfunction values without altering `xx` or indices. Because an inner `src` can be a pure outer point, outer fills must precede inner copies: `apply_bcs_outerradiation_and_inner()` fills outer RHS points using wavespeed and `f_infinity` before inner BCs, while `apply_bcs_outerextrap_and_inner()` uses 3-point inward second-order extrapolation along each outer face normal (requiring `>= 3` interior points along every axis with an outer face) before inner BCs.

Evidence (descriptive behavior): `CurviBoundaryConditions/bcstruct_set_up.py`; `CurviBoundaryConditions/apply_bcs_inner_only.py`; `CurviBoundaryConditions/apply_bcs_outerextrap_and_inner.py`.

### Hard Requirements For Angular Inner Boundaries

Angular coordinates must span their `ReferenceMetric` cell-face intervals (`xxmin*`, `xxmax*`, and `PI` are fixed by code generation, not parfile inputs), and the post-`convergence_factor` integer counts `Nxx*` (`convergence_factor` scales every axis with `Nxx != 2`) must satisfy:
- Even `Nphi` wherever an inner identification uses `phi -> phi + pi` (poles, origin, cylindrical axis, `SymTP` focal segment; `Ntheta = Nphi = 2` 1D spherical and `Nphi = 2` axisymmetric setups satisfy this).
- Even `xx1` count on `SymTP` and `SinhSymTP` (an odd count places cell centers at `xx1 = pi/2`, where roundoff under `sqrt` in `SymTP` `Cart_to_xx[1]` produces `NaN` indices in `bcstruct_set_up()`).
- `Ntheta >= 2 NGHOSTS` on `RingHoleySinhSpherical` so theta ghost zones stay strictly away from poles (coarser grids reach poles where single-pass `apply_bcs_inner_only()` can chain inner-to-inner points or abort). All other theta and wedge counts may be odd or even; Cartesian-like and `GeneralRFM` grids have no angular coordinates (`SymTP` lists `xx1` in `radial_like_dirns` near the focal segment while keeping `[0, pi]` bounds).

| Coordinate systems | Angular coordinates and defining bounds | Identifications used by inner boundaries |
| --- | --- | --- |
| `Spherical`, `SinhSpherical`, `SinhSphericalv2n*` | `theta = xx1` in `[0, pi]`; `phi = xx2` in `[-pi, pi]` | seam `phi -> phi +/- 2 pi`; poles `theta -> -theta` or `2 pi - theta` with `phi -> phi + pi`; origin `r -> -r`, `theta -> pi - theta`, `phi -> phi + pi` |
| `HoleySinhSpherical` | `theta = xx1` in `[0, pi]`; `phi = xx2` in `[-pi, pi]` | seam and poles; inner radial face is an outer boundary |
| `RingHoleySinhSpherical` | `theta = xx1` in `[pi/4, 3 pi/4]`; `phi = xx2` in `[-pi, pi]` | seam only; theta and inner radial faces are outer boundaries |
| `Cylindrical`, `SinhCylindrical`, `SinhCylindricalv2n*` | `phi = xx1` in `[-pi, pi]` | seam; axis `rho -> -rho` with `phi -> phi + pi` |
| `SymTP`, `SinhSymTP` | `xx1` in `[0, pi]`; `phi = xx2` in `[-pi, pi]` | seam; poles `xx1 -> -xx1` or `2 pi - xx1` with `phi -> phi + pi`; focal segment `xx0 -> -xx0` with `phi -> phi + pi` |
| `LWedgeHSinhSph`, `UWedgeHSinhSph` | `theta = xx1` in `[pi/4, 3 pi/4]`; `phi = xx2` in `[-pi/4, pi/4]` | none; eigencoordinate `Cartesian` makes all ghost points pure outer |

Cell centers sit at `xx[d][i] = xxmin_d + (i - NGHOSTS + 1/2) dxx_d` (`dxx_d = (xxmax_d - xxmin_d) / Nxx_d`): seam shifts move `k` by `+/- Nphi`, `phi -> phi + pi` moves `k` by `Nphi / 2` (integer iff `Nphi` is even), and `theta` reflections map `j` to `-1 - j`, `2 Ntheta - 1 - j`, or `Ntheta - 1 - j` for any `Ntheta`. Shifted/shortened angular intervals, odd counts where even is required, and `RingHoleySinhSpherical` with `Ntheta < 2 NGHOSTS` are excluded configurations and never justify coordinate bounds checks or stencil clamping; see [Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md).

Evidence (scientific/interface requirement): `reference_metric.py::ReferenceMetric.__init__`; `numerical_grids_and_timestep.py`; `CurviBoundaryConditions/bcstruct_set_up.py`; `reference_metric_Spherical.py`.

### Inner Ghost Points In Interpolation Stencils

On valid grids ([Hard Requirements For Angular Inner Boundaries](#hard-requirements-for-angular-inner-boundaries)), filled inner ghost points store the field at the ghost point's Cartesian location in orthonormal basis vectors at `xx[d][i]` (possibly outside `[xxmin, xxmax]`), so Lagrange stencils cross seams, poles, axes, and origins without wrapping or clamping (helpers without ghost zones, such as geodesic spatial tables, map extended indices to stored records instead).

Distinguish three quantities: the logical target `xx` from `Cart_to_xx_and_nearest_i0i1i2_assume_valid` (`acos` in `[0, pi]`, `atan2` in `(-pi, pi]`), its unclamped nearest index `(int)((xx - xxmin)/dxx + NGHOSTS)`, and the stencil center computed independently by `interpolation_3d_general__uniform_src_grid` using the same rounding and node differences `xx_target - xx[d][i]`.
Targets in `[xxmin, xxmax)` land on interior centers `[NGHOSTS, Nxx + NGHOSTS - 1]` whose half-width `NGHOSTS` stencils stay in bounds; upper-endpoint targets (`theta = pi` on `-z` or `phi = +pi` from `atan2(+0, x < 0)`) get center `Nxx + NGHOSTS` (first upper ghost point), overhang the array, and return `INTERP3D_GENERAL_HORIZON_OUT_OF_BOUNDS` without writing output. Remapping only the center near `phi = +pi` to `phi = -pi` would create `~2 pi` Lagrange coordinate differences (extrapolation), whereas shifting the center is governed by [Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md).

Verify actual BC fills per channel: evolved fields receive outer + inner fills at initial data and every MoL substep ([MoL Time Integration](mol-time-integration.md)), whereas `diagnostic_gfs_set` applies `apply_bcs_inner_only_specific_gfs()` only to `DIAG_PSI4_REGF`/`DIAG_PSI4_IMGF` (when `enable_psi4=True`) and `DIAG_HAMILTONIANGF`, `DIAG_MGF`, `DIAG_LAMBDA_CONSTRAINTGF` (only when `enable_interp_diagnostics=True`), leaving pure outer ghost layers (and inner points sourced from them) at unset `BHAH_MALLOC` values.
Finally, parity signs act only on curvilinear-basis components; Cartesian-basis tensors at two logical representatives of the same physical point require neither sign flips nor rotations (rotations apply only across symmetry-related distinct Cartesian points, as in [Geodesics And Raytracing Runtime](geodesics-and-raytracing-runtime.md)).

Evidence (descriptive behavior): `xx_tofrom_Cart.py`; `BHaHAHA/interpolation_3d_general__uniform_src_grid.py`; `BHaH/general_relativity/diagnostic_gfs_set.py`; `BHaH/diagnostics/diagnostics.py`.

## Sources

- [coding_style.md](../../../coding_style.md) - `## Coordinate Bounds-Check Prohibition`
- [register_all.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/register_all.py) - `register_C_functions`
- [bcstruct_set_up.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/bcstruct_set_up.py) - `register_CFunction_bcstruct_set_up`, `Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt`
- [BHaH_defines.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/BHaH_defines.py) - `register_BHaH_defines_h`, `BHaH_defines_set_gridfunction_defines_with_parity_types`
- [apply_bcs_outerradiation_and_inner.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerradiation_and_inner.py) - `register_CFunction_apply_bcs_outerradiation_and_inner`
- [apply_bcs_outerextrap_and_inner.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerextrap_and_inner.py) - `register_CFunction_apply_bcs_outerextrap_and_inner`
- [apply_bcs_inner_only.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only.py) - `register_CFunction_apply_bcs_inner_only`
- [apply_bcs_inner_only_specific_gfs.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only_specific_gfs.py) - `register_CFunction_apply_bcs_inner_only_specific_gfs`
- [reference_metric.py](../../../nrpy/reference_metric.py) - `ReferenceMetric.__init__`, `spherical_like`, `cylindrical_like`, `prolate_spheroidal_like`, `spherical_wedge_like`
- [reference_metric_Spherical.py](../../../nrpy/tests/reference_metric_Spherical.py) - `trusted_dict`
- [numerical_grids_and_timestep.py](../../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py) - `register_CFunction_numerical_grid_params_Nxx_dxx_xx`
- [xx_tofrom_Cart.py](../../../nrpy/infrastructures/BHaH/xx_tofrom_Cart.py) - `register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid`
- [interpolation_3d_general__uniform_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_3d_general__uniform_src_grid.py) - `register_CFunction_interpolation_3d_general__uniform_src_grid`
- [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostic_gfs_set.py) - `register_CFunction_diagnostic_gfs_set`
- [diagnostics.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics.py) - `register_all_diagnostics`

## See Also

- Depends on: [Grids, Coordinates, And Reference Metrics](grids-coordinates-and-reference-metrics.md)
- Depends on: [Reference Metrics](../../core/reference-metrics.md)
- Depends on: [Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md)
