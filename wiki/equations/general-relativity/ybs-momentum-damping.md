# YBS-MOM Timestep-Scaled Momentum Adjustment

> Status: confirmed
> Up: [General Relativity](index.md)

## Summary

YBS-MOM is a default-disabled, direct Yo--Lin--Cao momentum-gradient
adjustment to the BSSN-shaped conformal extrinsic-curvature equation. It adds
no evolved cleaning field. fCCZ4 reuses the same BSSN-shaped term.

For the matter-complete lower conformal momentum residual

\[
\mathcal M_i=\bar\gamma^{jk}\bar D_j\bar A_{ki}
+6\bar A^j{}_i\partial_j\phi
-\frac23\partial_iK-8\pi S_i,
\]

the optional addition is

\[
\left.\partial_t\bar A_{ij}\right|_{\rm YBS\text{-}MOM}
=C_{\rm YBS\_mom}\,\mathrm{CFL\_FACTOR}\,\Delta s_{\min}
\left[
\frac12(\bar D_i\mathcal M_j+\bar D_j\mathcal M_i)
-\frac13\bar\gamma_{ij}\bar\gamma^{kl}\bar D_k\mathcal M_l
\right].
\]

The coefficient follows the CAHD pattern: because
`CFL_FACTOR * DSMINGF` is the local timestep scale, the parabolic coefficient
shrinks with resolution. The explicit timestep limit of the diffusion term is
therefore proportional to the spacing, like the hyperbolic CFL limit, rather
than to its square; the recommended maximum of `C_YBS_mom` is given under
[Strength and timestep bound](#strength-and-timestep-bound). This is
timestep-scaled parabolic damping; it does not change the PDE into a hyperbolic
relaxation system.

Claim evidence:
- Claim: enabled YBS-MOM is the displayed direct, matter-complete, covariant STF momentum-gradient addition with a CAHD-style local timestep prefactor, so the explicit diffusion timestep limit scales with the spacing like the hyperbolic CFL limit; it adds no evolved cleaner field and makes no full-system hyperbolicity or stability guarantee.
- Role: public/scientific contract
- Deciding authority: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch
- Corroboration: [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (56); [Etienne, Phys. Rev. D 110, 064045](https://doi.org/10.1103/PhysRevD.110.064045), Eq. (26)

## Detail

### State and boundaries

No `qD`, `QD`, relaxation equation, cleaner wavespeed, or cleaner-specific
boundary field exists. Initial-data conversion, checkpointing, AMR transfer,
Kreiss--Oliger routing, and radiation/extrapolation boundaries therefore keep
the original evolved-variable set. Enabling YBS-MOM changes existing
`a_rhsDD` expressions only.

### Shared minimum-spacing infrastructure

`DSMINGF` is the one shared CAHD/YBS-MOM auxiliary gridfunction containing raw
local physical spacing. `ds_min_single_pt_exprs` is the common expression
source for both the canonical single-point spacing routine and the all-points
`DSMINGF` fill. Ordinary orthogonal coordinates use
`Abs(scalefactor_orthog[i] * dxx[i])`; fisheye GeneralRFM uses
`sqrt(Abs(ghatDD[i][i])) * Abs(dxx[i])`. Unsupported non-fisheye GeneralRFM
is rejected explicitly. RHS code, not the fill, applies `CFL_FACTOR` and the
consumer coefficient.

Claim evidence:
- Claim: CAHD and YBS-MOM share one `DSMINGF` and one coordinate-aware source of physical spacing expressions; the spacing field stores no consumer-specific coefficient.
- Role: descriptive behavior
- Deciding authority: [numerical_grids_and_timestep.py](../../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py), `ds_min_single_pt_exprs`; [dsmin_gf.py](../../../nrpy/infrastructures/BHaH/general_relativity/dsmin_gf.py), `register_CFunction_dsmin_auxevol_gridfunction`
- Corroboration: [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), CAHD scaling and YBS-MOM registration branches

### Ownership and defaults

`BSSNRHSs` owns \(\mathcal M_i\), its derivative, the covariant STF
projection, and the addition to `Abar_rhsDD`. Its matter branch differentiates
the stress-energy contribution as well. `FCCZ4RHSs` reuses that adjusted BSSN
base. `register_CFunction_rhs_eval` conditionally registers `DSMINGF` and the
runtime `C_YBS_mom` parameter, whose default `0.0` removes the term exactly.
Public APIs default the option to `False`.
The black-hole spectroscopy generator schedules the shared spacing fill when
either CAHD or YBS-MOM needs it.

All affected existing trusted-expression owners enable YBS Gamma and YBS-MOM
together: six BSSN dictionaries, six fCCZ4 dictionaries, and eight BHaH
`rhs_eval` dictionaries. No new test or trusted-dictionary family was added.

### Strength and timestep bound

`C_YBS_mom` is zero by default, which disables the term exactly. A user who
enables the option sets it above zero in the parameter file. The timestep is
`Delta t = CFL_FACTOR * ds_min`, and the coefficient is
`ell_M = C_YBS_mom * CFL_FACTOR * dsmin`, which in BHaH is the same prefactor as
`C_CAHD`. In the cell where `dsmin` equals `ds_min`, the product `ell_M * Delta t / h^2` equals
`C_YBS_mom * CFL_FACTOR^2`; in every coarser cell it is smaller by the ratio
`ds_min / dsmin`.

At principal order the term makes the momentum constraint diffuse, with rates
`-ell_M |k|^2 / 2` (transverse) and `-2 ell_M |k|^2 / 3` (longitudinal). For
the centered finite-difference stencils that NRPy generates, the largest
eigenvalue magnitude of the discrete operator is `lambda_FD * ell_M / h^2`,
with `lambda_FD` equal to 4.00, 5.62, 6.74, and 7.59 for orders 2, 4, 6, and 8.
For the diffusion term alone, an explicit Runge-Kutta step is stable when
`C_YBS_mom * CFL_FACTOR^2 * lambda_FD` does not exceed the real-axis limit of
the method: 2.785 for RK4, 2.513 for RK3, and 2.0 for RK2. The recommended
maximum is half of that limit. For finite-difference orders up to 8 it is
`0.18 / CFL_FACTOR^2` with RK4, and `0.13 / CFL_FACTOR^2` for any method whose
real-axis limit is at least 2. At `CFL_FACTOR = 0.45` these are 0.89 and 0.64.
For orders 4, 6, and 8 the diffusive eigenvalue peaks near 0.75 pi per cell in
each direction; for order 2 it equals its maximum along the whole line with pi
in two directions and any wavenumber in the third. The first-derivative symbol
is nonzero at these points, so propagation at the speed of light does lower the
limit. At `CFL_FACTOR = 0.45`, adding the unit-speed frequency of a
three-dimensional wavevector to every diffusive mode, the RK4 limit on
`C_YBS_mom * CFL_FACTOR^2 * lambda_FD` is 2.77, 2.70, 2.63, and 2.56 for orders
2, 4, 6, and 8 (2.785 without propagation), and the RK3 limit is 2.44, 2.20,
2.06, and 1.96 (2.513 without propagation). The recommended half of the limit
stays below these reduced values.

The analysis is frozen-coefficient and principal-part only. It omits variable
coefficients, larger gauge speeds, and Kreiss--Oliger dissipation, and no
evolution confirmed it; the safety factor of 1/2 is a chosen value, not a
derived bound. The bound also has no lapse factor because the coefficient
carries none.

Claim evidence:
- Claim: the runtime default `C_YBS_mom = 0.0` removes the term exactly, and the recommended maximum is half of the explicit-Runge-Kutta diffusion limit `real-axis limit / (lambda_FD * CFL_FACTOR^2)` at the finest cell; the analysis is frozen-coefficient and principal-part only, and no evolution confirmed it.
- Role: descriptive behavior
- Deciding authority: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch (`ell_M`); [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval` (default value); [numerical_grids_and_timestep.py](../../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py), `register_CFunction_cfl_limited_timestep`
- Corroboration: [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (56), for the operator only; no repository script reproduces the `lambda_FD`, real-axis-limit, and propagation values

### Authority reconciliation

The preserved commissioned source records an earlier auxiliary-relaxation
proposal. It remains frozen historical input, but it does not describe the
implemented contract. The user-directed CAHD-style direct formulation and
living code decide current behavior; [CONTR-0007](../../contradictions.md#contr-0007)
records this resolved mismatch.

### Dendro block spacing

The Dendro examples expose this term with `--ybs-momentum`, independently of
`--ybs-gamma`. Dendro computes `dsmin` once per Cartesian block as
`min(abs(dx), abs(dy), abs(dz))`; it broadcasts the scalar in SIMD kernels.
The generated RHS uses `C_YBS_mom * BSSN_CFL_FACTOR * dsmin`, with the same
CFL factor used for evolution. No spacing gridfunction is needed on a block
with uniform Cartesian spacing. See [BSSN Application Wiring](../../infrastructures/dendro/bssn-application-wiring.md#optional-yo-et-al-adjustments)
for defaults and the shared fCCZ4 interface.

Claim evidence:
- Claim: Dendro evaluates the local momentum coefficient from Cartesian block spacing and the evolution CFL parameter, while reusing the canonical momentum-gradient equation.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`, `register_CFunction_rhs_eval`; `nrpy/infrastructures/Dendro/CodeParameters.py`, `Q1_TOML_PARAMETER_NAMES`.
- Corroboration: `nrpy/examples/dendro_bssn.py` and `nrpy/examples/dendro_fccz4.py`, generation options and RHS registration.

## Sources

- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - direct momentum adjustment
- [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py) - fCCZ4 reuse
- [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py) - option, parameter, and spacing ownership
- [dsmin_gf.py](../../../nrpy/infrastructures/BHaH/general_relativity/dsmin_gf.py) - shared all-points spacing fill
- [numerical_grids_and_timestep.py](../../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py) - shared physical-spacing expressions
- [preserved YBS-MOM specification](../../../raw/source-docs/ybs-momentum-damping-spec.md) - superseded relaxation proposal retained as frozen provenance
- [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2) - direct momentum-gradient adjustment
- [Etienne, Phys. Rev. D 110, 064045](https://doi.org/10.1103/PhysRevD.110.064045) - CAHD timestep-scaling pattern

## See Also

- Parent: [General Relativity](index.md)
- Depends on: [BSSN Family](bssn-family.md)
- Used by: [Fully Covariant Conformal Z4](fccz4.md)
- Implemented by: [BHaH GR Application Wiring](../../infrastructures/bhah/gr-application-wiring.md)
