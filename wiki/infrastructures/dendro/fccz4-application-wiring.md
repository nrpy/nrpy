# fCCZ4 Application Wiring

> fCCZ4 equations, CAHD and Z4 damping parameters, initialization, and 10-field constraint diagnostics in `Dendro_NRPy_fCCZ4`. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_fccz4` generates `Dendro_NRPy_fCCZ4` with code-generation conformal-factor selection (`W` default or `chi`), appending `Theta_fCCZ4` as the 25th evolved field after the 24 BSSN-ordered components. It uses centered finite differences, slow-start lapse (SSL), constant default `eta = 1`, its own `dt`-scaled fCCZ4 CAHD formula, and separate `Ricci_eval_order_N` and `rhs_eval_order_N` kernels while sharing TwoPunctures, GW/ADM extraction, parameter dumps, and VTU output with BSSN.

## Detail

### fCCZ4 Equations, CAHD, And Runtime Parameters

Each generated directory writes `pars/fccz4.toml` (generated defaults, including `BSSN_CAHD_C = 0.15`) and `pars/q1.par.lowres.toml` (identical to the BSSN copy, overriding `BSSN_CAHD_C` to `0.06`, `ETA_CONST` to `2.0`, and `BSSN_ELE_ORDER` to `6`; see [q1 versus generated defaults](bssn-application-wiring.md#q1-versus-generated-defaults) and [Production Validation And Deferred Checks](production-validation-and-deferred-checks.md)). `--ybs-gamma`, `--ybs-momentum`, and in-script KO follow [Optional Yo et al. adjustments](bssn-application-wiring.md#optional-yo-et-al-adjustments), with the fCCZ4 shift driver receiving the YBS-adjusted fCCZ4 connection RHS rather than Brown's BSSN term.

Unlike BSSN's `dx^2 / dt` scaling, Dendro's fCCZ4 CAHD scales linearly with `dt` using runtime `C_CAHD` (`BSSN_CAHD_C`, default `0.15`; factor `4` vs `2` from `chi = W^2`):

```text
W_rhs   += 2 C_CAHD W H dt     (W evolution),
chi_rhs += 4 C_CAHD chi H dt   (chi evolution).
```

Here `H` is the BSSN-form Hamiltonian expression (also output alongside `H_Z4` by `fCCZ4_constraints_order_N`), not `H_Z4`. Default `OnePlusLog` slicing evolves `d_t alpha = beta^i d_i alpha - 2 alpha (K - 2 Theta)` plus the SSL term from [BSSN Application Wiring](bssn-application-wiring.md) (`BSSN_SSL_H`, `BSSN_SSL_SIGMA`). Z4 damping coefficients `kappa1` (default `0.1`, no extra `alpha` factor) and `kappa2` (default `0`) are read from TOML keys `kappa1` and `kappa2` (omitted in `q1.par.lowres.toml`), while `kappa3 = 1` is fixed.

Evidence (descriptive behavior): `Dendro/general_relativity/rhs_eval.py::register_CFunction_rhs_eval`; `fCCZ4_RHSs.py`; `equations/general_relativity/fCCZ4_constraints.py`; `Dendro/CodeParameters.py::Q1_TOML_PARAMETER_NAMES`.

### Initialization, Block Traversal, And Ten-Field Diagnostics

Initial octrees use the same analytic Dendro-GR puncture seed as BSSN, load the same precomputed TwoPunctures spectral file (`psi = psi_background + u`, `alpha = W = psi^(-2)`), run the same `ADM_to_BSSN -> zip/exchange/unzip -> physical boundary -> initial_data_lambdaU -> zip -> floor/projection` sequence with `Theta_fCCZ4 = 0`, and tag checkpoints distinctly by formulation (`W` vs `chi` and fCCZ4 vs BSSN). After initialization, each RK stage, and AMR transfer, `alpha` (`CHI_FLOOR`) and conformal factor (`sqrt(CHI_FLOOR)` or `CHI_FLOOR`) are floored before `enforce_detgbar_equals_detghat_trAzero`.

Block traversal follows BSSN (`exterior padding -> Ricci -> RHS -> physical-face outgoing radiation BCs`; see [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)) with six `RbarDD` scratch fields. Apparent-horizon searches, `*_ADM.dat`, `*_GW_l<l>_m<m>.dat`, `*_GW_L2.dat`, `*_BHLocations.dat`, and `__PARAM_DUMP__*.toml` follow [Initial Data, Wave Extraction, And Output Files](initial-data-and-wave-diagnostics.md).

`fCCZ4_constraints_order_N` writes ten diagnostic fields to both excised physical-volume and unique-node RMS tables after remeshing: `H_Z4`, `Z4constraintU0..2` (`Lambdatilde^i - DeltaGamma^i = 2 Zbar^i`, reported without enforcement; no `Theta` column is written), `H`, physical lower-index `MD0..MD2`, `M_CONSTRAINT = sqrt(gamma_ij M^i M^j)`, and `LAMBDA_CONSTRAINT = sqrt(gammabar_ij Z4constraintU^i Z4constraintU^j)`. In VTU output, evolved index `24` selects `Theta_fCCZ4` and constraint index `0` selects BSSN-form `H` (`H_Z4` and `Z4constraintU` are not VTU-selectable).

Evidence (descriptive behavior): `Dendro/main_cpp.py`; `Dendro/solver_context.py`; `Dendro/general_relativity/ADM_to_BSSN.py`; `Dendro/general_relativity/fCCZ4_constraints.py::register_CFunction_fCCZ4_constraints`.

## Sources

- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - `parse_args`, `main`
- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - shared TwoPunctures registration
- [param_toml.py](../../../nrpy/infrastructures/Dendro/param_toml.py) - `generate_default_parfile`
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - `Q1_TOML_PARAMETER_NAMES`
- [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py) - `kappa1`, `kappa2` registration and damping terms
- [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py) - `Z4constraintU`, `ZbarU`, `H_Z4`
- [Dendro fCCZ4_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py) - `register_CFunction_fCCZ4_constraints`
- [Ricci_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - conformal Ricci kernel
- [ADM_to_BSSN.py](../../../nrpy/infrastructures/Dendro/general_relativity/ADM_to_BSSN.py) - `Theta_fCCZ4 = 0` initialization
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - block traversal, diagnostics, VTU
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - puncture seed and startup

## See Also

- Depends on: [fCCZ4 Formulation And Gauge Equations](../../equations/general-relativity/fccz4-formulation.md)
- Contrasts with: [BSSN Application Wiring](bssn-application-wiring.md)
- See also: [Initial Data, Wave Extraction, And Output Files](initial-data-and-wave-diagnostics.md)
