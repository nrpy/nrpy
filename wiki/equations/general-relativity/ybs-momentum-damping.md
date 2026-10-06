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
=C_{\rm YBS\_mom}\,\mathrm{CFL\_FACTOR}\,\Delta s_{\min}\,W
\left[
\frac12(\bar D_i\mathcal M_j+\bar D_j\mathcal M_i)
-\frac13\bar\gamma_{ij}\bar\gamma^{kl}\bar D_k\mathcal M_l
\right],
\qquad W=e^{-2\phi}.
\]

The coefficient follows the CAHD pattern: because
`CFL_FACTOR * DSMINGF` is the local timestep scale, the parabolic coefficient
shrinks with resolution. The explicit timestep limit of the diffusion term is
therefore proportional to the spacing, like the hyperbolic CFL limit, rather
than to its square. The conformal factor `W` also multiplies the coefficient;
it is 1 far from black holes and vanishes at a puncture (see
[Conformal-factor weight](#conformal-factor-weight)). The recommended value of
`C_YBS_mom` is 1.75, given with its range under
[Strength and timestep bound](#strength-and-timestep-bound). This is
timestep-scaled, `W`-weighted parabolic damping; it does not change the PDE
into a hyperbolic relaxation system.

Claim evidence:
- Claim: enabled YBS-MOM is the displayed direct, matter-complete, covariant STF momentum-gradient addition with a CAHD-style local timestep prefactor and the conformal-factor weight `W`, so the explicit diffusion timestep limit scales with the spacing like the hyperbolic CFL limit; it adds no evolved cleaner field and makes no full-system hyperbolicity or stability guarantee.
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

### Conformal-factor weight

`BSSNRHSs` multiplies the coefficient by the conformal factor
`W = e^{-2 phi}`, built once from the evolved variable in the same way as the
CAKO Kreiss--Oliger weight: `cf` for `W`, `sqrt(cf)` for `chi`, and
`exp(-2 cf)` for `phi`. Far from black holes `W` is 1. At a puncture `W`
vanishes, the fields are under-resolved, and the added term contains second
derivatives of `phi`, which grow like `1/r^2`. The weight multiplies the added
term by `W` and so attenuates it there. It does not by itself make the term
vanish, because `W` times the derivative of the computed residual need not
vanish: in Cartesian coordinates with flat `gammabar_ij`, a constant trace-free
error `delta Abar_ij` gives the residual error `6 delta Abar_ij d_j phi`. For
example, for `delta Abar_ij = epsilon diag(1, -1, 0)` on the +x axis, `W` times
the derivative of that residual error tends to a nonzero constant for
`psi ~ M / (2 r)` (`W ~ r^2`) and grows like `1/r` on the static trumpet
(`W = r / (r + M)`). The weight does not change what the term means
elsewhere: the addition is a multiple of the derivative of the momentum
residual, so it still vanishes on an exact solution for any weight. The
preserved specification
records the paper's coefficient with a lapse window; NRPy uses `W` in its place
and carries no lapse factor. The choice of `W`, like the recommended strength,
is a maintainer decision and not a result derived in the cited sources.

Claim evidence:
- Claim: the YBS-MOM coefficient is multiplied by the conformal factor `W = e^{-2 phi}`, built from the evolved conformal-factor option; the weight is 1 far from black holes and vanishes at a puncture, and the addition still vanishes on an exact solution; the weight attenuates the addition near a puncture but does not by itself make it vanish, because `W` times the derivative of the computed residual need not vanish: for example, for `delta Abar_ij = epsilon diag(1, -1, 0)` with flat `gammabar_ij` in Cartesian coordinates, on the +x axis `W` times the derivative of the residual error `6 delta Abar_ij d_j phi` tends to a nonzero constant for `psi ~ M / (2 r)` and grows like `1/r` on the static trumpet.
- Role: descriptive behavior
- Deciding authority: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch (`W`, `ell_M`)
- Corroboration: [kreiss_oliger_terms.py](../../../nrpy/equations/general_relativity/kreiss_oliger_terms.py), `add_KreissOliger_dissipation_terms` (the same `W` pattern for CAKO); none available for the choice of `W`, which is a maintainer decision; none available for the two puncture scalings, which follow from the momentum-residual expression in `BSSNRHSs` and the stated `W`

### Ownership and defaults

`BSSNRHSs` owns \(\mathcal M_i\), its derivative, the covariant STF
projection, and the addition to `Abar_rhsDD`. Its matter branch differentiates
the stress-energy contribution as well. `FCCZ4RHSs` reuses that adjusted BSSN
base. `register_CFunction_rhs_eval` conditionally registers `DSMINGF` and the
runtime `C_YBS_mom` parameter, whose default `0.0` removes the term exactly and
whose recommended value is 1.75. Public APIs default the option to `False`.
The black-hole spectroscopy generator schedules the shared spacing fill when
either CAHD or YBS-MOM needs it.

All affected existing trusted-expression owners enable YBS Gamma and YBS-MOM
together.

### Strength and timestep bound

`C_YBS_mom` is zero by default, which disables the term exactly. A user who
enables the option sets it above zero in the parameter file; the recommended
value is 1.75. The timestep is `Delta t = CFL_FACTOR * ds_min`, and the
coefficient is `ell_M = C_YBS_mom * CFL_FACTOR * dsmin * W`, whose first three
factors are, in BHaH, the same prefactor as `C_CAHD`. In a cell of spacing
`h = dsmin`, the product `ell_M * Delta t / h^2` equals
`C_YBS_mom * CFL_FACTOR^2 * W * ds_min / dsmin`. Because `ds_min <= dsmin`, it
is at most `C_YBS_mom * CFL_FACTOR^2 * W_max`, where `W_max` is the largest `W`
on the grid at any step. Where `W <= 1` everywhere, as for Brill--Lindquist and
static trumpet initial data, `W_max <= 1` and the product is at most
`C_YBS_mom * CFL_FACTOR^2`, reached only where `W` is near 1 and `dsmin` equals
`ds_min`. Data with `W > 1`, such as Kasner data at physical time `t < 1`
(`W = t^(-1/3)`), can exceed that bound by up to the factor `W_max`.

At principal order the term makes the momentum constraint diffuse, with rates
`-ell_M |k|^2 / 2` (transverse) and `-2 ell_M |k|^2 / 3` (longitudinal). For
the centered finite-difference stencils that NRPy generates, the largest
eigenvalue magnitude of the discrete operator is `lambda_FD * ell_M / h^2`,
with `lambda_FD` equal to 4.00, 5.62, 6.74, and 7.59 for orders 2, 4, 6, and 8.
For the diffusion term alone, an explicit Runge-Kutta step is stable in a cell
when `C_YBS_mom * CFL_FACTOR^2 * lambda_FD * W * ds_min / dsmin` does not
exceed the real-axis limit of the method: 2.785 for RK4, 2.513 for RK3, and 2.0
for RK2.

The condition holds in every cell when
`C_YBS_mom * CFL_FACTOR^2 * lambda_FD * W_max` does not exceed the limit; where
`W <= 1` everywhere it suffices that `C_YBS_mom * CFL_FACTOR^2 * lambda_FD` does
not exceed the limit. At
`CFL_FACTOR = 0.45` with RK3 and `W <= 1` that means `C_YBS_mom` up to about 2.2
for order 4 and 1.6 for order 8 (1.9 and 1.3 once propagation is included,
below); for `W_max > 1`, divide each of these by `W_max`. For `W <= 1`, a larger
value stays stable only while no cell has both `W` and `ds_min / dsmin` near 1.
The finest cells of a black-hole grid lie beside the punctures, where `W` is
small. For `W <= 1` everywhere, the range 0 (off) to about 2 at
`CFL_FACTOR = 0.45` follows from this condition, and its upper end scales as
`1 / CFL_FACTOR^2`. The recommended value 1.75 is a maintainer choice, not a
derived bound.

For orders 4, 6, and 8 the diffusive eigenvalue peaks near 0.75 pi per cell in
each direction; for order 2 it equals its maximum along the whole line with pi
in two directions and any wavenumber in the third. The first-derivative symbol
is nonzero at these points, so propagation at the speed of light does lower the
limit. At `CFL_FACTOR = 0.45`, adding the unit-speed frequency of a
three-dimensional wavevector to every diffusive mode, the RK4 limit on
`C_YBS_mom * CFL_FACTOR^2 * lambda_FD` is 2.77, 2.70, 2.63, and 2.56 for orders
2, 4, 6, and 8 (2.785 without propagation), and the RK3 limit is 2.44, 2.20,
2.06, and 1.96 (2.513 without propagation). These reduced values replace the
real-axis limits in the condition above.

The analysis is frozen-coefficient and principal-part only. It omits variable
coefficients, larger gauge speeds, Kreiss--Oliger dissipation, and the
lower-order terms with derivatives of `phi` and of the connection, which vanish
in flat space and are largest at a puncture; the weight `W` attenuates the term
in that region, and no analysis bounds the remaining term. The coefficient
carries `W` but no lapse factor.

Claim evidence:
- Claim: the runtime default `C_YBS_mom = 0.0` removes the term exactly; the recommended value is 1.75, a maintainer choice and not a derived bound; the explicit-Runge-Kutta diffusion condition is `C_YBS_mom * CFL_FACTOR^2 * lambda_FD * W * ds_min / dsmin` below the real-axis limit in each cell; because `ds_min <= dsmin`, it holds in every cell when `C_YBS_mom * CFL_FACTOR^2 * lambda_FD * W_max` does not exceed the limit, with `W_max` the largest `W` on the grid at any step, and where `W <= 1` everywhere it suffices that `C_YBS_mom * CFL_FACTOR^2 * lambda_FD` does not exceed the limit; the range 0 to about 2 at `CFL_FACTOR = 0.45` is stated for `W <= 1` everywhere; the analysis is frozen-coefficient and principal-part only and omits the strong-field lower-order terms.
- Role: descriptive behavior
- Deciding authority: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch (`ell_M`, `W`); [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval` (default value and parameter description); [numerical_grids_and_timestep.py](../../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py), `register_CFunction_cfl_limited_timestep`
- Corroboration: [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (56), for the operator only; none available for the `lambda_FD`, real-axis-limit, and propagation values, because they follow from NRPy's centered finite-difference stencils and the Runge-Kutta stability limits, which Eq. (56) does not treat; none available for the recommended value, which is a maintainer choice

### Authority reconciliation

The preserved commissioned source records an earlier auxiliary-relaxation
proposal. It remains frozen historical input, but it does not describe the
implemented contract. The user-directed CAHD-style direct formulation and
living code decide current behavior; [CONTR-0007](../../contradictions.md#contr-0007)
records this resolved mismatch. A later maintainer decision multiplies the
coefficient by `W`; it leaves the direct, stateless formulation unchanged.

### Dendro block spacing

The Dendro examples expose this term with `--ybs-momentum`, independently of
`--ybs-gamma`. Dendro computes `dsmin` once per Cartesian block as
`min(abs(dx), abs(dy), abs(dz))`; it broadcasts the scalar in SIMD kernels.
The generated RHS uses `C_YBS_mom * BSSN_CFL_FACTOR * dsmin * W`, with the same
CFL factor used for evolution and `W` built from the evolved conformal factor
(`W` or `chi`). No spacing gridfunction is needed on a block
with uniform Cartesian spacing. See [BSSN Application Wiring](../../infrastructures/dendro/bssn-application-wiring.md#optional-yo-et-al-adjustments)
for defaults and the shared fCCZ4 interface.

Claim evidence:
- Claim: Dendro evaluates the local momentum coefficient from Cartesian block spacing, the evolution CFL parameter, and the conformal-factor weight `W`, while reusing the canonical momentum-gradient equation.
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
