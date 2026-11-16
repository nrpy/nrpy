# YBS-MOM Timestep-Scaled Momentum Adjustment

> Direct Yo-Lin-Cao symmetric trace-free momentum-gradient damping term with `CFL_FACTOR * dsmin * W` scaling for BSSN and fCCZ4. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

YBS-MOM (`enable_YBS_momentum_constraint_adjustment=True`, default `False`) adds a direct, stateless Yo-Lin-Cao symmetric trace-free gradient of the matter-complete lower conformal momentum residual to `Abar_rhsDD` in BSSN and fCCZ4 without introducing an auxiliary evolved field:

```text
M_i = gammabar^jk Dbar_j Abar_ki + 6 Abar^j_i partial_j phi - (2/3) partial_i K - 8 pi S_i,

(partial_t Abar_ij)|_YBS-MOM = C_YBS_mom * CFL_FACTOR * dsmin * W
                               * [ (Dbar_i M_j + Dbar_j M_i)/2 - (1/3) gammabar_ij gammabar^kl Dbar_k M_l ],
```

where `W = exp(-2 phi)`. Because `CFL_FACTOR * dsmin` scales proportionally to the local grid spacing (like CAHD), the explicit parabolic timestep constraint scales as `O(h)` rather than `O(h^2)` (recommended `C_YBS_mom = 1.75`, runtime default `0.0`).

## Detail

### Stateless Formulation, Shared `DSMINGF`, And Conformal Weight

YBS-MOM adds no `qD`/`QD` relaxation variables, wavespeeds, initial-data steps, checkpoint fields, AMR transfers, KO terms, or boundary conditions; it modifies only `a_rhsDD` in `BSSNRHSs` (including stress-energy derivatives when `enable_T4munu=True`), which `FCCZ4RHSs` reuses directly. All BSSN and fCCZ4 `RHSs` `trusted_dict` modules enable YBS-Gamma and YBS-MOM jointly.

In BHaH, CAHD and YBS-MOM share one auxiliary gridfunction `DSMINGF` populated from `numerical_grids_and_timestep.py::ds_min_single_pt_exprs` (`Abs(scalefactor_orthog[i] * dxx[i])` on orthogonal reference metrics; `sqrt(Abs(ghatDD[i][i])) * Abs(dxx[i])` on fisheye `GeneralRFM`; non-fisheye `GeneralRFM` rejected), while RHS expressions apply `CFL_FACTOR` and `C_YBS_mom`. In Dendro (`--ybs-momentum`), uniform Cartesian blocks compute scalar `dsmin = min(abs(dx), abs(dy), abs(dz))` per block without a gridfunction ([Optional Yo et al. adjustments](../../infrastructures/dendro/bssn-application-wiring.md#optional-yo-et-al-adjustments)).

### Conformal-factor weight

`BSSNRHSs` builds `W = exp(-2 phi)` from the active conformal-factor variable (`cf` for `W`, `sqrt(cf)` for `chi`, `exp(-2 cf)` for `phi`, matching CAKO in `kreiss_oliger_terms.py`, with no lapse factor). The added `a_rhsDD` term is a multiple of the derivative of the momentum residual, so it vanishes on exact constraint solutions for any weight. While `W -> 1` asymptotically and attenuates `1 / r^2` second derivatives of `phi` at punctures (`W -> 0`), `W` times the differentiated discrete residual need not vanish at a puncture: for flat `gammabar_ij` and constant error `delta Abar_ij = epsilon diag(1, -1, 0)` along `+x`, `W` times the derivative of `6 delta Abar_ij partial_j phi` approaches a nonzero constant for `psi ~ M / (2r)` (`W ~ r^2`) and scales as `1 / r` on a static trumpet (`W = r / (r + M)`).

Evidence (scientific/interface requirement): `BSSN_RHSs.py::BSSNRHSs.__init__`; [Yo, Lin, and Cao](https://arxiv.org/pdf/1205.5111v2) Eq. (56); [Etienne](https://doi.org/10.1103/PhysRevD.110.064045) Eq. (26); `BHaH/general_relativity/dsmin_gf.py`.

### Strength and timestep bound

`C_YBS_mom` defaults to `0.0` (exact no-op); `1.75` is a maintainer recommendation rather than a general stability theorem. With global timestep `Delta t = CFL_FACTOR * ds_min`, local spacing `h = dsmin`, and `ell_M = C_YBS_mom * CFL_FACTOR * dsmin * W`,

```text
ell_M * Delta t / h^2 = C_YBS_mom * CFL_FACTOR^2 * W * (ds_min / dsmin)
                     <= C_YBS_mom * CFL_FACTOR^2 * W_max
```

(`W <= 1` on Brill-Lindquist and static-trumpet initial data, whereas Kasner at `t < 1` has `W = t^(-1/3) > 1`).

In a frozen-coefficient Fourier model with `gammabar_ij = delta_ij` and uniform Cartesian spacing `h`, continuum modes decay at `-ell_M |k|^2 / 2` (transverse) and `-2 ell_M |k|^2 / 3` (longitudinal). Because NRPy discretizes pure second derivatives with `-q_i / h^2` and mixed derivatives with `-p_i p_j / h^2` (`q_i != p_i^2`), the discrete 5-component STF operator has maximum eigenvalue magnitude `lambda_FD * ell_M / h^2` with `lambda_FD ~= 4.00, 5.62, 6.74, 7.59` at centered FD orders `2, 4, 6, 8` (peaking near `0.75 pi` per cell for orders `4..8`, and at two `pi` axes for order `2`). Against explicit RK negative-real-axis limits (`~2.785` RK4, `~2.513` RK3, `2.0` RK2), pure-diffusion stability in that model requires:

```text
C_YBS_mom * CFL_FACTOR^2 * lambda_FD * W * (ds_min / dsmin) <= RK_real_limit
```

(giving `C_YBS_mom <= ~2.2` at FD4 and `~1.6` at FD8 for RK3 with `CFL_FACTOR = 0.45`, `W <= 1`, scaling as `1 / (W_max * CFL_FACTOR^2)`; adding unit-speed wave frequencies at `CFL_FACTOR = 0.45` yields `C_YBS_mom * CFL_FACTOR^2 * lambda_FD <= ~2.77, 2.70, 2.63, 2.56` for RK4 and `~2.44, 2.20, 2.06, 1.96` for RK3 across orders `2..8`). Anisotropic conformal metrics exceed the identity-metric bound—for `gammabar_ij = diag(1/16, 4, 4)` (`det = 1`), an `x`-Nyquist `A_xy` mode reaches `-32 ell_M / h^2` vs `-4 ell_M / h^2` at FD2—and the frozen model omits variable coefficients, gauge speeds, boundaries, KO dissipation, and lower-order connection/`phi` couplings.

Evidence (descriptive behavior): `BSSN_RHSs.py::BSSNRHSs.__init__`; `finite_difference.py::compute_fdcoeffs_fdstencl`; `rk_butcher_table_dictionary.py::generate_Butcher_tables`; `BHaH/general_relativity/rhs_eval.py`; `Dendro/general_relativity/rhs_eval.py`.

## Sources

- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs.__init__` (`DbarM_STFDD`, `ell_M`)
- [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py) - `FCCZ4RHSs.__init__` reuse
- [kreiss_oliger_terms.py](../../../nrpy/equations/general_relativity/kreiss_oliger_terms.py) - `add_KreissOliger_dissipation_terms` (`W` construction)
- [finite_difference.py](../../../nrpy/finite_difference.py) - `compute_fdcoeffs_fdstencl`
- [rk_butcher_table_dictionary.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/rk_butcher_table_dictionary.py) - `generate_Butcher_tables`
- [BHaH rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval` (`C_YBS_mom`, `dsmin`)
- [dsmin_gf.py](../../../nrpy/infrastructures/BHaH/general_relativity/dsmin_gf.py) - `register_CFunction_dsmin_auxevol_gridfunction`
- [numerical_grids_and_timestep.py](../../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py) - `ds_min_single_pt_exprs`, `register_CFunction_cfl_limited_timestep`
- [Dendro rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`
- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - `--ybs-momentum`
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - `--ybs-momentum`
- [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2) - Eq. (56)
- [Etienne, Phys. Rev. D 110, 064045](https://doi.org/10.1103/PhysRevD.110.064045) - Eq. (26)

## See Also

- Depends on: [BSSN Family](bssn-family.md)
- Implements: [BHaH GR Evolution, Ricci, And Constraints](../../infrastructures/bhah/gr-evolution-and-constraints.md)
- Implements: [Dendro BSSN Application Wiring](../../infrastructures/dendro/bssn-application-wiring.md)
