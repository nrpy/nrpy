# Grids, Coordinates, Reference Metrics, And Boundaries

> BHaH route for numerical grid setup, coordinate wrappers, reference-metric precompute, GeneralRFM, fisheye hooks, and curvilinear boundaries. Status: confirmed.
> Up: [BHaH](index.md)

## Summary

BHaH grid setup starts by registering per-grid numerical parameters, allocating
cell-centered coordinate arrays, and selecting either independent grids or a
multipatch setup. It then optionally populates reference-metric precompute data,
sets curvilinear boundary metadata, chooses a CFL-limited timestep, and registers
coordinate-system wrappers so generated code can dispatch by
`params->CoordSystem_hash`.

Grids with angular coordinates must use each coordinate system's defining
angular bounds. The number of phi gridpoints must be even wherever an inner
boundary identification uses the half-turn `phi -> phi + pi`, at a pole, the
origin, the cylindrical axis, or the `SymTP` focal segment. `SymTP` and
`SinhSymTP` also need an even `xx1` count, and `RingHoleySinhSpherical` needs
`Ntheta >= 2 NGHOSTS`; every other theta count may be odd or even.
[Hard Requirements For Angular Inner Boundaries](#hard-requirements-for-angular-inner-boundaries)
states these requirements before
[Inner Ghost Points In Interpolation Stencils](#inner-ghost-points-in-interpolation-stencils)
explains how interpolation reads the filled ghost points.

## Detail

`register_CFunctions()` in `numerical_grids_and_timestep.py` registers the core
grid parameters: `Nxx*`, `Nxx_plus_2NGHOSTS*`, `xxmin*`, `xxmax*`, `dxx*`,
`invdxx*`, `convergence_factor`, `CoordSystem_hash`, `grid_idx`, and
`gridname`. `numerical_grid_params_Nxx_dxx_xx()` sets default `Nxx` from the
generator's `Nxx_dict`, allows a runtime `Nx[]` override, applies
`commondata->convergence_factor` except along symmetry axes with `Nxx == 2`,
sets `Nxx_plus_2NGHOSTS`, evaluates coordinate extents from reference-metric
parameters, computes `dxx` and `invdxx`, and allocates cell-centered `xx[3]`
arrays. CUDA setup allocates these arrays on device and copies per-grid params
to `d_params`.

`numerical_grids_and_timestep()` supports `"independent grid(s)"` directly and
delegates `"multipatch"` to `multipatch_grids_set_up()`. Independent-grid setup
assigns `gridname`, `CoordSystem_hash`, `grid_physical_size`, calls the
coordinate-grid initializer for each coordinate system, records `NUMGRIDS`, and
sets `grid_idx`. First-time setup initializes `nn`, `nn_0`, `t_0`, and `time`.
When CFL setup is enabled, `cfl_limited_timestep()` computes
`commondata->dt = min(commondata->dt, ds_min * commondata->CFL_FACTOR)` across
active grids. For ordinary coordinate systems, `ds_min_single_pt()` uses
reference-metric scale factors; for `GeneralRFM_fisheyeN*`, it uses the fisheye
metric diagonal. Non-fisheye `GeneralRFM` reports unsupported `ds_min`.

Reference-metric precompute has two paths. `rfm_precompute.py` handles ordinary
non-`GeneralRFM` precompute by discovering coordinate-dependent expressions,
creating `rfm_struct`, and registering `rfm_precompute_malloc`,
`rfm_precompute_defines`, and `rfm_precompute_free`; generated host and CUDA
paths allocate, populate, read, and free lookup arrays. `generalrfm_precompute.py`
handles `GeneralRFM_fisheyeN*` by storing `ghatDD`, `ghatUU`, `detgammahat`,
`ghatDDdD`, and `ghatDDdDD` into `AUXEVOL` gridfunctions. That GeneralRFM
precompute path rejects CUDA and unsupported providers. `generalrfm_cart_to_xx.py`
registers a coordinate-specialized `generalrfm_Cart_to_xx` Newton solve that
inverts a supported GeneralRFM fisheye map and returns nonzero on failure.
`phys_params_to_fisheye.py` registers physical fisheye CodeParameters in
`commondata_struct`, computes internal `fisheye_a*`, `fisheye_R*`,
`fisheye_s*`, and `fisheye_c`, and provides a post-params hook that runs the
conversion after `params_struct_set_to_default()`.

Coordinate wrappers live in two layers. `xx_tofrom_Cart.py` registers the Python
registrars `register_CFunction_xx_to_Cart` and
`register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid`, which emit the
coordinate-specific C functions `xx_to_Cart` and
`Cart_to_xx_and_nearest_i0i1i2_assume_valid`. For an independent grid, the
inverse unshifts Cartesian points by `Cart_origin*` for local-grid calculations
and maps them to `xx`. When the caller supplies `Cart_to_i0i1i2`, the emitted C
function then converts `xx` to nearest indices with `xxmin`, `dxx`, and
`NGHOSTS`; this conversion assumes a valid, in-domain, cell-centered point and
does not perform a bounds check. When `Cart_to_i0i1i2 == NULL`, it returns after
computing `xx`, before any floating-to-integer index conversion. For a
multipatch grid with `params->grid_rotates`, inverse conversion builds the
cumulative rotation matrix and applies `R^T` to the Cartesian vector before
local-origin handling and native-coordinate inversion; forward `xx_to_Cart`
applies `R` after native-to-Cartesian mapping and origin handling. GeneralRFM
fisheye wrappers use a radial Newton solve for inverse mapping and a closed-form
fisheye radius map for forward mapping. Each multipatch registrar also registers
its exact SO(3) dependency closure: inverse registers the matrix builder and
`R^T` vector helper, while forward registers the matrix builder and `R` vector
helper. Registration remains duplicate-safe in either converter order and
independent-grid registration adds no SO(3) helper.

Evidence (descriptive behavior): [xx_tofrom_Cart.py](../../../nrpy/infrastructures/BHaH/xx_tofrom_Cart.py), `register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid` and `register_CFunction_xx_to_Cart` (doctests and generated-source baselines in the deciding module)

Maintenance rule: coordinate-admission checks are request-gated by
[Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md);
the nullable coordinates-only path does not authorize adding them.

Owner full-source baselines cover the default independent-grid inverse and
forward converters plus explicitly named Cartesian multipatch OpenMP/CUDA
variants. The OpenMP GeneralRFM fisheye inverse is intentionally excluded from
full-source golden comparison; its semantic doctest compiles generated OpenMP C
harnesses with `REAL=float` and `REAL=double`, checks finite round-trip fixtures,
and checks recoverable failure for non-finite Cartesian input. No broader
application build, CUDA GeneralRFM runtime, GPU run, restart run, distributed
run, or production numerical-result guarantee follows.

Evidence (generated evidence): [xx_tofrom_Cart.py](../../../nrpy/infrastructures/BHaH/xx_tofrom_Cart.py), `register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid`, `register_CFunction_xx_to_Cart`, and `_run_generalrfm_fisheye_inverse_roundtrip_check`; corroboration: [generalrfm_cart_to_xx.py](../../../nrpy/infrastructures/BHaH/generalrfm_cart_to_xx.py), `register_CFunction_generalrfm_Cart_to_xx`

GeneralRFM CUDA converters are unsupported in these coordinate conversion
registrars.

Evidence (descriptive behavior): [xx_tofrom_Cart.py](../../../nrpy/infrastructures/BHaH/xx_tofrom_Cart.py), `register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid` and `register_CFunction_xx_to_Cart`; corroboration: [generalrfm_cart_to_xx.py](../../../nrpy/infrastructures/BHaH/generalrfm_cart_to_xx.py), `register_CFunction_generalrfm_Cart_to_xx`

`rfm_wrapper_functions.py` creates non-coordinate-specific wrapper functions
that switch on `params->CoordSystem_hash`, calls the matching
coordinate-specific function, and registers uppercase coordinate hash macros in
`BHaH_defines.h`.

Curvilinear boundary registration starts in
`CurviBoundaryConditions/register_all.py`. It registers `outer_bc_type` with
default `"radiation"`, the coordinate-specific `bcstruct_set_up()`, radiation
outer plus inner BCs, inner-only BCs, specific-gridfunction inner-only BCs,
extrapolation outer plus inner BCs, `bc_struct`'s `griddata_struct`
contribution, and BHaH header definitions. CUDA builds also register
`cpyHosttoDevice_bc_struct()`.

The boundary header defines `innerpt_bc_struct`, `outerpt_bc_struct`,
`bc_info_struct`, and `bc_struct`. `bc_struct` owns `inner_bc_array`,
`pure_outer_bc_array[NGHOSTS*3]`, and `bc_info`. Parity arrays are inferred from
registered gridfunction names for evolved gridfunctions, optionally auxiliary
and auxevol gridfunctions; high-rank GeneralRFM auxevol metric fields receive a
placeholder parity only to preserve table indexing.

`bcstruct_set_up()` classifies every ghost point by a double map through the
coordinate system's eigencoordinate: `Spherical` for spherical-like systems,
`Cylindrical`, `SymTP`, or `Cartesian`, which also serves the wedge systems;
`GeneralRFM` uses its own numerical inverse. The ghost point's stored logical
coordinates go to Cartesian coordinates through the eigencoordinate
`xx_to_Cart`, back through the eigencoordinate `Cart_to_xx`, and then to the
nearest stored index. A ghost point whose index maps to itself is a pure outer
boundary point, stored by ghost-zone layer and direction with face signs. Every
other ghost point is an inner boundary point: its dst index is the ghost point,
its src index is the stored point it maps to, and it carries ten parity signs.
The src is an interior point or, where an inner boundary meets the outer
boundary, a pure outer boundary point. Before storing a pair,
`bcstruct_set_up()` maps the dst and src logical coordinates to Cartesian
coordinates with the coordinate system's own map and exits with an error if the
two positions differ by more than a relative tolerance. Each parity sign is the
dot product of the orthonormal unit vectors evaluated at the dst and src logical
coordinates; a dot product farther from +/-1 than its tolerance also exits.
Optional masks skip unset points and identify outer-boundary points; CUDA
masking is rejected by numerical-grid setup.

Evidence (descriptive behavior): [bcstruct_set_up.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/bcstruct_set_up.py), `register_CFunction_bcstruct_set_up`, `Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt`, and `Cfunction__set_parity_for_inner_boundary_single_pt`; corroboration: [apply_bcs_inner_only.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only.py), `register_CFunction_apply_bcs_inner_only` description of inner points that map to interior or pure outer points (generated-source baselines only; no test fills a grid)

Boundary application changes stored gridfunction values only; it does not move a
ghost point, change `xx` or a stored index, or interpolate between positions.
`apply_bcs_inner_only()` loops over the `NUM_EVOL_GFS` gridfunctions of the
array it receives and over every inner boundary point, and sets the dst value to
the src value, with its sign reversed when the parity entry selected by the
gridfunction's parity type in `evol_gf_parity` is -1. CPU-only
`apply_bcs_inner_only_specific_gfs()` does the same for a caller-supplied
gridfunction list and parity array. Parity types 0, 1-3, and 4-9 select the
scalar, vector-component, and symmetric rank-2 tensor-component signs in the
rescaled curvilinear basis. Because a src can be a pure outer boundary point,
the inner copy must follow the outer fill.
`apply_bcs_outerradiation_and_inner()` fills pure outer points first with
radiation BCs using wavespeed and `f_infinity`, then applies inner BCs to RHS
gridfunctions. `apply_bcs_outerextrap_and_inner()` fills pure outer points by
second-order extrapolation, then applies inner BCs; a specific-gridfunction
variant performs the same sequence for selected grid functions. Informal names
such as `apply_inner_bcs` or `apply_inner_boundaries` refer to
`apply_bcs_inner_only()` and `apply_bcs_inner_only_specific_gfs()`.

Evidence (descriptive behavior): [apply_bcs_inner_only.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only.py), `generate_apply_bcs_inner_only__kernel_body`; corroboration: [apply_bcs_inner_only_specific_gfs.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only_specific_gfs.py), `register_CFunction_apply_bcs_inner_only_specific_gfs`

The second-order extrapolation in `apply_bcs_outerextrap_and_inner()` sets
each pure outer point from the three points next to it along the inward normal
of its face. Every direction with an outer face therefore needs at least three
interior gridpoints; with fewer, the extrapolation reads ghost points that do
not yet hold valid values.

Evidence (descriptive behavior): [apply_bcs_outerextrap_and_inner.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerextrap_and_inner.py), `generate_apply_bcs_outerextrap_and_inner_only__kernel_body` (no test fills a grid)

### Hard Requirements For Angular Inner Boundaries

Every angular coordinate must span its coordinate system's defining interval.
Where the table below lists an identification with the half-turn
`phi -> phi + pi`, at a pole, the origin, the cylindrical axis, or the `SymTP`
focal segment, the phi direction must also have an even number of interior
gridpoints; `HoleySinhSpherical` poles use the same half-turn as the
`Spherical` poles. `SymTP` and `SinhSymTP` also need an even `xx1` count, for
a different reason: an odd count puts a row of cell centers at `xx1 = pi/2`,
where the `SymTP` inverse map `Cart_to_xx[1]`, used by both systems through
their `SymTP` eigencoordinate, takes the square root of a difference that
roundoff can make negative, and `bcstruct_set_up()` then converts NaN to an
index.
`RingHoleySinhSpherical` needs `Ntheta >= 2 NGHOSTS`, which keeps every theta
ghost point strictly between the poles; neither of its counts then needs to be
even. A coarser theta grid is an excluded configuration: its theta ghost points
reach a pole, where `bcstruct_set_up()` can exit or can map some inner boundary
points to other inner boundary points, which the single copy pass of
`apply_bcs_inner_only()` can leave unfilled. Every other angular count,
including the wedge counts, may be odd. On a single grid, the 1D setup
`Ntheta = Nphi = 2` on `Spherical`, `SinhSpherical`, `SinhSphericalv2n*`, and
`HoleySinhSpherical`, and the axisymmetric setup with a phi count of 2, meet
these requirements; on `SymTP` and `SinhSymTP`, the axisymmetric setup also
needs an even `xx1` count.
These are required inputs to the inner-boundary mapping, not adjustable
defaults. The bounds are interval endpoints, which are cell faces rather than
the first and last cell centers. `ReferenceMetric` defines them, and
`numerical_grid_params_Nxx_dxx_xx()` copies them into `xxmin*` and `xxmax*`;
`xxmin*`, `xxmax*`, and `PI` are not parameter-file inputs. The counts are the
final `Nxx*` after that function multiplies every count other than 2 by
`convergence_factor` and truncates the product to an integer, so a
`convergence_factor` that makes a count violate these requirements gives an
excluded configuration.

| Coordinate systems | Angular coordinates and defining bounds | Identifications used by inner boundaries |
| --- | --- | --- |
| `Spherical`, `SinhSpherical`, `SinhSphericalv2n*` | theta = `xx1` in `[0, pi]`; phi = `xx2` in `[-pi, pi]` | seam `phi -> phi +/- 2 pi`; poles `theta -> -theta` or `2 pi - theta` with `phi -> phi + pi`; origin `r -> -r`, `theta -> pi - theta`, `phi -> phi + pi` |
| `HoleySinhSpherical` | theta = `xx1` in `[0, pi]`; phi = `xx2` in `[-pi, pi]` | seam and poles; the inner radial face is an outer boundary |
| `RingHoleySinhSpherical` | theta = `xx1` in `[pi/4, 3 pi/4]`; phi = `xx2` in `[-pi, pi]` | seam only; the theta faces and the inner radial face are outer boundaries |
| `Cylindrical`, `SinhCylindrical`, `SinhCylindricalv2n*` | phi = `xx1` in `[-pi, pi]` | seam; axis `rho -> -rho` with `phi -> phi + pi` |
| `SymTP`, `SinhSymTP` | `xx1` in `[0, pi]`; phi = `xx2` in `[-pi, pi]` | seam; poles `xx1 -> -xx1` or `2 pi - xx1` with `phi -> phi + pi`; focal segment `xx0 -> -xx0` with `phi -> phi + pi` |
| `LWedgeHSinhSph`, `UWedgeHSinhSph` | theta = `xx1` in `[pi/4, 3 pi/4]`; phi = `xx2` in `[-pi/4, pi/4]` | none; with eigencoordinate `Cartesian`, every ghost point is a pure outer boundary point |

Cartesian-like and `GeneralRFM` coordinate systems have no angular coordinate.
`SymTP` lists `xx1` in `radial_like_dirns` because near the focal segment it
mainly parametrizes z, but `xx1` keeps the bounds and reflections above.

Stored coordinates are cell centered:
`xx[d][i] = xxmin_d + (i - NGHOSTS + 1/2) dxx_d` with
`dxx_d = (xxmax_d - xxmin_d)/Nxx_d`. On the standard `Spherical` grid with
zero-based interior indices j and k, `theta_j = (j + 1/2) pi/Ntheta` and
`phi_k = -pi + (k + 1/2) 2 pi/Nphi`. The seam identification maps k to
`k +/- Nphi` for any `Nphi`. A shift `phi -> phi + pi` moves k by exactly
`Nphi/2` cells, so it lands on a stored phi value only when `Nphi` is even. The
reflections `theta -> -theta`, `theta -> 2 pi - theta`, and `theta -> pi - theta`
map j to `-1 - j`, `2 Ntheta - 1 - j`, and `Ntheta - 1 - j`, so they map theta
cell centers onto theta cell centers for any `Ntheta`. With the prescribed
bounds and the even counts above, every identification in the table therefore
recovers an existing same-grid point.

Evidence (public/scientific contract): [reference_metric.py](../../../nrpy/reference_metric.py), `spherical_like`, `cylindrical_like`, `prolate_spheroidal_like` (including its `SymTP` `Cart_to_xx`), and `spherical_wedge_like`; corroboration: [numerical_grids_and_timestep.py](../../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py), `register_CFunction_numerical_grid_params_Nxx_dxx_xx`; [apply_bcs_inner_only.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only.py), `generate_apply_bcs_inner_only__kernel_body`; [wave_equation_curvilinear.py](../../../nrpy/examples/wave_equation_curvilinear.py), `Nxx_dict`; [reference_metric_Spherical.py](../../../nrpy/tests/reference_metric_Spherical.py), `trusted_dict` entries `xxmin_1`, `xxmax_1`, `xxmin_2`, and `xxmax_2`; [bcstruct_set_up.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/bcstruct_set_up.py), `Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt` (no test covers excluded angular configurations; none is permitted)

Boundary filling relies on these requirements and cannot repair a grid that
violates them, because it only copies stored values. `bcstruct_set_up()` exits
when a dst position does not match its src position; that exit reports a
misaligned map but does not make it valid. A grid-setup function, exporter,
combiner, or reader that accepts other angular bounds, an odd count where the
requirements above need an even one, or a `RingHoleySinhSpherical` grid with
`Ntheta < 2 NGHOSTS` does not make them a permitted simulation: shifted or
shortened angular intervals, those odd counts, and those coarse ring grids are
excluded configurations. Do not use them to claim a
regression or to justify coordinate bounds checks, interpolation repairs, extra
modes, or coordinate-specific branches;
[Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md)
governs any new coordinate check. Verify another coordinate system against its
own defining bounds in `reference_metric.py`, not against parameter overrides.

### Inner Ghost Points In Interpolation Stencils

Read
[Hard Requirements For Angular Inner Boundaries](#hard-requirements-for-angular-inner-boundaries)
first. On a grid that satisfies them, the seam, poles, origin, cylindrical axis,
and SymTP focal segment are coordinate identifications, and filled inner ghost
points extend the stored logical arrays. A filled ghost value equals the field
at the ghost point's Cartesian position, expressed in unit vectors evaluated at
the ghost point's own logical coordinates, which may lie outside the nominal
interval. A Lagrange stencil may therefore cross the nominal angular interval
and read filled inner ghost points. Selecting a stencil over a ghost-filled
array needs no phi wrapping, theta clipping, pole or origin branch, or
coordinate-dependent stencil shift; coordinate-specific identifications belong
in the coordinate maps and in `bcstruct_set_up()`. A helper whose data contain
no ghost points, such as the geodesic spatial helper, maps its own extended
nodes to stored records instead.

Keep three quantities separate: the logical interpolation target, the selected
stencil center with its surrounding stored nodes, and any nearest stored index
used to classify the target. `Cart_to_xx_and_nearest_i0i1i2_assume_valid`
returns the logical target, with angular coordinates from `acos` in `[0, pi]` or
`atan2` in `(-pi, pi]`, and, when requested, the unclamped nearest index
`(int)((xx - xxmin)/dxx + NGHOSTS)`; it applies no clamp.
`interpolation_3d_general__uniform_src_grid` does not take that index. It
computes its own center from the logical target with the same rounding, and its
Lagrange weights use differences between the target and the stored node
coordinates `xx[d][i]`, ghost nodes included. A target in `[xxmin, xxmax)` gets
a center in `[NGHOSTS, Nxx + NGHOSTS - 1]`, so a stencil of half-width `NGHOSTS`
stays inside the stored array and may read inner ghost points. A target exactly
at an upper endpoint, `theta = pi` on the -z axis or `phi = +pi` from `atan2`
with `y = +0` and `x < 0`, gets center `Nxx + NGHOSTS`, the first upper inner
ghost point. Its stencil then extends past the array, the interpolator returns
`INTERP3D_GENERAL_HORIZON_OUT_OF_BOUNDS`, and that destination value is left
unwritten.

Evidence (descriptive behavior): [xx_tofrom_Cart.py](../../../nrpy/infrastructures/BHaH/xx_tofrom_Cart.py), `register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid`; corroboration: [interpolation_3d_general__uniform_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_3d_general__uniform_src_grid.py), `register_CFunction_interpolation_3d_general__uniform_src_grid` (no test covers endpoint targets)

Lagrange weights depend on differences between the target's logical coordinates
and the stored node coordinates; a physical identification does not make those
differences interchangeable. Near the seam, a target at `phi = +pi` has the raw
nearest index of the first upper inner ghost point, whose stored value is a copy
from near `phi = -pi`. Replacing only the stencil center by that periodic copy
while the target stays at `+pi` would produce coordinate differences near `2 pi`
and turn a local interpolation into distant extrapolation. A stencil centered at
the adjacent interior point `Nxx + NGHOSTS - 1` keeps the target and its nodes
on one coordinate branch and still reads the nearby filled ghost points. The
current selector neither remaps nor clamps the center; it reports the target out
of bounds. A stencil that reaches into inner ghost points differs from a stencil
whose center is an inner ghost point, and the first alone does not show that a
ghost-center path executes. Changing the selector is a coordinate-admission
change governed by
[Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md).

A classification label does not show that a field was filled when it is
interpolated; check the fill calls for the interpolated gridfunctions. Evolved
gridfunctions receive outer and inner fills with initial data and in each Method
of Lines stage ([MoL Time Integration](mol-time-integration.md)). Diagnostic
gridfunctions receive less: `constraints_eval` and `psi4` write interior points,
and `diagnostic_gfs_set` then applies `apply_bcs_inner_only_specific_gfs()` to
`DIAG_PSI4_REGF` and `DIAG_PSI4_IMGF` when Psi4 is enabled and to the constraint
channels `DIAG_HAMILTONIANGF`, `DIAG_MGF`, and `DIAG_LAMBDA_CONSTRAINTGF` when
interpolation diagnostics are enabled. No outer fill reaches these channels, so
their pure outer ghost points, and inner points whose src is a pure outer point,
keep the unset values of the `BHAH_MALLOC` allocation in `diagnostics()`. A
stencil that reaches the outer radial ghost layers reads those values; the
resulting error is in an output sample, not in an evolved field.

Evidence (descriptive behavior): [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostic_gfs_set.py), `register_CFunction_diagnostic_gfs_set`; corroboration: [diagnostics.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics.py), `register_all_diagnostics` allocation of `diagnostic_gfs`

Parity signs act on components in the rescaled curvilinear basis, whose unit
vectors reverse direction across a pole, origin, axis, or focal-segment
identification. Cartesian-basis components need neither a sign nor a rotation
when the same Cartesian position is reached through a different logical
representative: two representatives of one position carry the same Cartesian
tensor. A rotation is correct only between different Cartesian positions related
by a symmetry, such as the final rotation from a stored reference azimuth to the
target azimuth in
[Geodesics And Raytracing Runtime](geodesics-and-raytracing-runtime.md), or for
an explicit change of Cartesian basis.

## Sources

- [coding_style.md](../../../coding_style.md) - `## Coordinate Bounds-Check Prohibition`
- [numerical_grids_and_timestep.py](../../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py) - `register_CFunctions`, `register_CFunction_numerical_grid_params_Nxx_dxx_xx`, `register_CFunction_numerical_grids_and_timestep`, `register_CFunction_cfl_limited_timestep`, `register_CFunction_ds_min_single_pt`, `register_CFunction_ds_min_radial_like_dirns_single_pt`
- [rfm_precompute.py](../../../nrpy/infrastructures/BHaH/rfm_precompute.py) - `ReferenceMetricPrecompute`, `register_CFunctions_rfm_precompute`
- [rfm_wrapper_functions.py](../../../nrpy/infrastructures/BHaH/rfm_wrapper_functions.py) - `get_CoordSystem_hash`, `register_CFunctions_CoordSystem_wrapper_funcs`
- [xx_tofrom_Cart.py](../../../nrpy/infrastructures/BHaH/xx_tofrom_Cart.py) - `register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid`, `register_CFunction_xx_to_Cart`
- [generalrfm_precompute.py](../../../nrpy/infrastructures/BHaH/generalrfm_precompute.py) - `register_CFunction_generalrfm_precompute`, `register_CFunctions_generalrfm_support`
- [generalrfm_cart_to_xx.py](../../../nrpy/infrastructures/BHaH/generalrfm_cart_to_xx.py) - `register_CFunction_generalrfm_Cart_to_xx`
- [phys_params_to_fisheye.py](../../../nrpy/infrastructures/BHaH/fisheye/phys_params_to_fisheye.py) - `register_CFunction_fisheye_params_from_physical_N`, `build_post_params_struct_set_to_default_hook`
- [register_all.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/register_all.py) - `register_C_functions`
- [bcstruct_set_up.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/bcstruct_set_up.py) - `register_CFunction_bcstruct_set_up`
- [BHaH_defines.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/BHaH_defines.py) - `register_BHaH_defines_h`, `BHaH_defines_set_gridfunction_defines_with_parity_types`, `register_griddata_commondata`
- [apply_bcs_outerradiation_and_inner.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerradiation_and_inner.py) - `register_CFunction_apply_bcs_outerradiation_and_inner`, `setup_Cfunction_radiation_bcs`
- [apply_bcs_outerextrap_and_inner.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerextrap_and_inner.py) - `register_CFunction_apply_bcs_outerextrap_and_inner`, `register_CFunction_apply_bcs_outerextrap_and_inner_specific_gfs`, `generate_apply_bcs_outerextrap_and_inner_only__kernel_body`
- [apply_bcs_inner_only.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only.py) - `register_CFunction_apply_bcs_inner_only`, `generate_apply_bcs_inner_only__kernel_body`
- [apply_bcs_inner_only_specific_gfs.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only_specific_gfs.py) - `register_CFunction_apply_bcs_inner_only_specific_gfs`
- [reference_metric.py](../../../nrpy/reference_metric.py) - `ReferenceMetric.__init__`, `spherical_like`, `cylindrical_like`, `prolate_spheroidal_like`, `spherical_wedge_like`
- [reference_metric_Spherical.py](../../../nrpy/tests/reference_metric_Spherical.py) - `trusted_dict`
- [wave_equation_curvilinear.py](../../../nrpy/examples/wave_equation_curvilinear.py) - `Nxx_dict`
- [interpolation_3d_general__uniform_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_3d_general__uniform_src_grid.py) - `register_CFunction_interpolation_3d_general__uniform_src_grid`
- [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostic_gfs_set.py) - `register_CFunction_diagnostic_gfs_set`
- [diagnostics.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics.py) - `register_all_diagnostics`

## See Also

- [BHaH](index.md)
- Depends on: [Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md)
- [Runtime Data, Parameters, Headers, And CLI](runtime-data-parameters-headers-and-cli.md)
- [Reference Metrics](../../core/reference-metrics.md)
- See also: [Geodesics And Raytracing Runtime](geodesics-and-raytracing-runtime.md)
- See also: [BHaHAHA Grid, Interpolation, And Boundaries](bhahaha-grid-interpolation-and-boundaries.md)
