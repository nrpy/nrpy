# YBS-MOM Timestep-Scaled Momentum Adjustment

> Status: confirmed · Last reconciled: 08-30-2026
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
shrinks with resolution and avoids a timestep restriction stronger than the
existing hyperbolic CFL choice. This is timestep-scaled parabolic damping; it
does not change the PDE into a hyperbolic relaxation system.

Claim evidence:
- Claim: enabled YBS-MOM is the displayed direct, matter-complete, covariant STF momentum-gradient addition with a CAHD-style local timestep prefactor; it adds no evolved cleaner field and makes no full-system hyperbolicity or stability guarantee.
- Role: public/scientific contract
- Deciding authority: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch
- Corroboration: [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (56); [Etienne, Phys. Rev. D 110, 064045](https://doi.org/10.1103/PhysRevD.110.064045), Eq. (26)
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction plus BHaH OpenMP/CUDA source generation; precision=exact symbolic source plus 30-significant-digit deterministic trusted sampling; GPU=generation pass, execution not-run; restart=not-applicable because no new state; distributed=not-run; error_path=pass for unsupported non-fisheye GeneralRFM spacing; options=YBS Gamma and YBS-MOM enabled jointly in 6 BSSN, 6 fCCZ4, and 8 BHaH rhs_eval trusted cases; public defaults disabled; date=08-30-2026`

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
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, SymPy 1.14.0; backend=BHaH OpenMP/CUDA source generation; precision=exact symbolic comparison; GPU=generation pass, execution not-run; restart=not-applicable; distributed=not-run; error_path=pass; options=ordinary orthogonal, fisheye GeneralRFM, and rejected non-fisheye GeneralRFM routes; date=08-30-2026`

### Ownership and defaults

`BSSNRHSs` owns \(\mathcal M_i\), its derivative, the covariant STF
projection, and the addition to `Abar_rhsDD`. Its matter branch differentiates
the stress-energy contribution as well. `FCCZ4RHSs` reuses that adjusted BSSN
base. `register_CFunction_rhs_eval` conditionally registers `DSMINGF` and the
runtime `C_YBS_mom=1.0` parameter. Public APIs default the option to `False`.
The black-hole spectroscopy generator schedules the shared spacing fill when
either CAHD or YBS-MOM needs it.

All affected existing trusted-expression owners enable YBS Gamma and YBS-MOM
together: six BSSN dictionaries, six fCCZ4 dictionaries, and eight BHaH
`rhs_eval` dictionaries. No new test or trusted-dictionary family was added.

### Authority reconciliation

The preserved commissioned source records an earlier auxiliary-relaxation
proposal. It remains frozen historical input, but it does not describe the
implemented contract. The user-directed CAHD-style direct formulation and
living code decide current behavior; [CONTR-0007](../../contradictions.md#contr-0007)
records this resolved mismatch.

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
