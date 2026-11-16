# fCCZ4 Formulation And Gauge Equations

> Three-dimensional reference-metric fCCZ4 continuum equations, algebraic constraints, Z4 Ricci tensor, gauge conditions, and paper-convention resolutions in NRPy. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

NRPy implements the three-dimensional reference-metric fully covariant conformal Z4 (fCCZ4) formulation in spacetime signature `(-,+,+,+)` with spatial dimension three, exact rational coefficients, and full covariance (`kappa3 = 1` fixed). The continuum equations follow Mewes et al. Eqs. (3)-(41) after removing a duplicate connection-vector stretch in their printed shift bracket so the evolved conformal connection `Lambdatilde^i` receives a single Lie-derivative stretch, corroborated in Cartesian coordinates by Alic et al. Eq. (19) and structurally by Sanchis-Gual et al. Eqs. (2.11) and (2.17). Brown supplies the covariant reference-metric connection definitions.

## Detail

### Metric Split, Matter Projections, And Algebraic Constraints

The 3+1 line element, future-pointing unit normal, and stress-energy projections are

```text
ds^2 = -alpha^2 dt^2 + gamma_ij (dx^i + beta^i dt)(dx^j + beta^j dt),
n^mu = alpha^(-1) (1, -beta^i),
rho  = n_mu n_nu T^(mu nu),          S_i  = -gamma_(i mu) n_nu T^(mu nu),
S_ij = gamma_(i mu) gamma_(j nu) T^(mu nu),   S = gamma^ij S_ij,
```

where `rho` corresponds to Mewes et al.'s `E`. Following the displayed evolution system in Mewes et al., `kappa1` damping terms carry no extra factor of `alpha`; `kappa1` has inverse-length dimension (default `0.1`) and dimensionless `kappa2` defaults to `0`.

The conformal decomposition uses

```text
gamma_ij = exp(4 phi) gammabar_ij,
Abar_ij  = exp(-4 phi) (K_ij - gamma_ij K / 3),   K = gamma^ij K_ij,
```

with evolved conformal-factor options `phi`, `W = exp(-2 phi)`, or `chi = exp(-4 phi)` (`BSSNQuantities` raises `ValueError` for other names). The Lagrangian condition `partial_t det(gammabar) = 0` preserves `det(gammabar) = det(gammahat)` analytically from constrained initial data, though equation-module `__main__` checks do not test time preservation directly.

Where symbolic callers apply `BSSN_algebraic_constraints.py::BSSN_algebraic_constraints` (used by backend `enforce_detgbar_equals_detghat_trAzero` kernels; see [GR Evolution, Ricci, And Constraints](../../infrastructures/bhah/gr-evolution-and-constraints.md) and [Generated Backend Comparison](../../generated-backend-comparison.md)), all twelve independent `(hDD, aDD)` components are projected onto `det(gammabar') = det(gammahat)` and `(gammabar')^ij Abar'_ij = 0` without a `trAbar` damping variable or RHS term. On `GeneralRFM*` (`_project_general_rfm`, positive-definite `G = gammabar` and `gammahat`),

```text
q     = (abs(det(gammahat)) / det(G))^(1/3),   G'_ij = q G_ij,
tau'  = (G')^ij A_ij,                          A'_ij = A_ij - G'_ij tau' / 3,
```

whereas on standard orthogonal coordinate systems (`_project_orthogonal_rfm`), diagonal scale factors cancel in `H_ij = delta_ij + h_ij`, unrescaled `H` is inverted once, and `h'_ij = det(H)^(-1/3) H_ij - delta_ij`, `a'_ij = a_ij - H_ij (H^kl a_kl) / 3`. No guard is generated for singular, indefinite, or non-positive-determinant metrics.

With time-independent reference metric `gammahat_ij` (covariant derivative `Dhat_i`, connection `Gammahat^i_jk`) and conformal-metric quantities barred,

```text
DeltaGamma^i_jk = Gammabar^i_jk - Gammahat^i_jk,
DeltaGamma^i    = gammabar^jk DeltaGamma^i_jk,
Lambdatilde^i   = DeltaGamma^i + 2 gammabar^ij Z_j.
```

Evidence (scientific/interface requirement): [Mewes et al.](https://arxiv.org/pdf/2002.06225v2) Eqs. (3)-(10), (38)-(41); `fCCZ4_RHSs.py::FCCZ4RHSs.__init__`; `BSSN_quantities.py::BSSNQuantities.__init__`; `BSSN_algebraic_constraints.py::BSSN_algebraic_constraints`.

### Z4 Vector, Z4 Ricci Tensor, And Constraints

The normal Z4 projection, conformal connection constraint, and spatial Z4 components are

```text
Theta  = -n_mu Z^mu,                 C^i = Lambdatilde^i - DeltaGamma^i,
Zbar^i = C^i / 2,   Z^i = exp(-4 phi) C^i / 2,   Z_i = gammabar_ij C^j / 2.
```

Lowering `Z^i` with physical `gamma_ij` cancels `exp(-4 phi)` and requires the full, potentially off-diagonal `gammabar_ij`. With physical covariant derivative `D_i` and geometric conformal Ricci tensor `RbarGeom_ij`, `RbarZ4_ij = RbarGeom_ij + D_i Z_j + D_j Z_i`. Absorbing `Z_i` into `Lambdatilde^i` and expanding the physical connection `Gamma^k_jl` gives NRPy's derivative-free split:

```text
RbarZ4_ij = RbarNRPy_ij(Lambdatilde) + deltaR_ij,   RbarZ4 = gammabar^ij RbarZ4_ij,
deltaR_ij = -2 (C_i partial_j phi + C_j partial_i phi)
            + 2 gammabar_ij C^l partial_l phi
            + (C^l / 2) (gammabar_ki DeltaGamma^k_jl + gammabar_kj DeltaGamma^k_il),
C_i       = gammabar_ij C^j.
```

The pure-trace term `+2 gammabar_ij C^l partial_l phi` comes from expanding the physical connection and affects both `partial_t K` and `partial_t Theta`.

`FCCZ4Constraints` defines the Z4 Hamiltonian expression without `-2 Theta K` (appended explicitly in `partial_t K` and `partial_t Theta`):

```text
H_Z4 = exp(-4 phi) [RbarZ4 - 8 Dbar^i phi Dbar_i phi - 8 Dbar^i Dbar_i phi]
       + 2 K^2 / 3 - Abar_ij Abar^ij - 16 pi rho.
```

At the symbolic equation level, fCCZ4 corrections vanish and the RHSs reduce to NRPy's BSSN equations when `Theta = C^i = H_Z4 = Abar^i_i = 0` together with vanishing derivatives of those zero fields; full physical constraint satisfaction additionally requires `det(gammabar) = det(gammahat)` and the conformal momentum constraint (`BSSN_constraints.py::BSSNconstraints.MU`):

```text
M^i = exp(-4 phi) [Dbar_j Abar^ij + 6 Abar^ij partial_j phi - 2 gammabar^ij partial_j K / 3]
      - 8 pi S^i = 0.
```

Evidence (scientific/interface requirement): [Mewes et al.](https://arxiv.org/pdf/2002.06225v2) Eqs. (24)-(32); `fCCZ4_constraints.py::FCCZ4Constraints.__init__`; `BSSN_constraints.py::BSSNconstraints.__init__`; [Brown](https://arxiv.org/pdf/0902.3652v2) Eqs. (12a)-(12b), (15).

### Coordinate-Time Evolution Equations

With `D_beta = Dbar_k beta^k`, the evolved system is

```text
partial_t gammabar_ij = beta^k partial_k gammabar_ij + gammabar_ik partial_j beta^k
                        + gammabar_kj partial_i beta^k - 2 gammabar_ij D_beta / 3
                        - 2 alpha Abar_ij,

partial_t phi   = beta^k partial_k phi + (D_beta - alpha K) / 6,
partial_t W     = -2 W partial_t phi,     partial_t chi = -4 chi partial_t phi,

partial_t Abar_ij = beta^k partial_k Abar_ij + Abar_ik partial_j beta^k
                    + Abar_kj partial_i beta^k - 2 Abar_ij D_beta / 3
                    - 2 alpha Abar_ik Abar^k_j + alpha Abar_ij (K - 2 Theta)
                    + exp(-4 phi) [-2 alpha Dbar_i Dbar_j phi + 4 alpha Dbar_i phi Dbar_j phi
                                   + 2(Dbar_i alpha Dbar_j phi + Dbar_j alpha Dbar_i phi)
                                   - Dbar_i Dbar_j alpha + alpha(RbarZ4_ij - 8 pi S_ij)]^TF,

partial_t K     = beta^i partial_i K
                  + exp(-4 phi) [alpha(RbarZ4 - 8 Dbar^i phi Dbar_i phi - 8 Dbar^i Dbar_i phi)
                                 - (2 Dbar^i alpha Dbar_i phi + Dbar^i Dbar_i alpha)]
                  + alpha(K^2 - 2 Theta K) - 3 kappa1(1 + kappa2) Theta
                  + 4 pi alpha(S - 3 rho),

partial_t Theta = beta^i partial_i Theta + alpha(H_Z4 - 2 Theta K) / 2
                  - Z^i partial_i alpha - kappa1(2 + kappa2) Theta.
```

The conformal connection combines `Lambdabar_rhsU_without_Brown_constraint_term` (`LambdaBase^i`) with `deltaLambda^i`:

```text
LambdaBase^i  = beta^k partial_k Lambdatilde^i - Lambdatilde^k partial_k beta^i
                + gammabar^jk Dhat_j Dhat_k beta^i + 2 DeltaGamma^i D_beta / 3
                + gammabar^ij Dbar_j D_beta / 3
                - 2 Abar^ij (partial_j alpha - 6 alpha partial_j phi)
                + 2 alpha Abar^jk DeltaGamma^i_jk - 4 alpha gammabar^ij partial_j K / 3
                - 16 pi alpha gammabar^ij S_j,

deltaLambda^i = 2 gammabar^ij (alpha partial_j Theta - Theta partial_j alpha)
                - 2 alpha K C^i / 3 - kappa1 C^i + 2 C^i Dhat_k beta^k / 3,

partial_t Lambdatilde^i = LambdaBase^i + deltaLambda^i.
```

Because `LambdaBase^i` already stretches the full evolved connection `beta^k Dhat_k Lambdatilde^i - Lambdatilde^k Dhat_k beta^i`, `deltaLambda^i` omits a second `-C^k Dhat_k beta^i` stretch while promoting the shift-divergence coefficient from `DeltaGamma^i` to `Lambdatilde^i`, giving Cartesian `partial_b beta^a` coefficient `-Lambdatilde^b delta^i_a + (2/3) Lambdatilde^i delta^a_b`.

Evidence (scientific/interface requirement): `BSSN_RHSs.py::BSSNRHSs.__init__`; `fCCZ4_RHSs.py::FCCZ4RHSs.__init__`; [Alic et al.](https://arxiv.org/pdf/1106.2254v2) Eq. (19); [Sanchis-Gual et al.](https://arxiv.org/pdf/1403.3653v1) Eqs. (2.11), (2.17).

### Gauge Conditions And Paper-Convention Resolutions

`fCCZ4_gauge_RHSs` accepts the four BSSN lapse conditions (`OnePlusLog` default, `BHSHarmonicSlicing`, `Frozen`, `OnePlusLogAlt`), raising `ValueError` via `BSSN_gauge_RHSs` when `BHSHarmonicSlicing` is paired with `chi`. Only `OnePlusLog` receives an NRPy `Theta` coupling (`OnePlusLogAlt` and `BHSHarmonicSlicing` keep their BSSN forms):

```text
partial_t alpha = beta^i partial_i alpha - 2 alpha (K - 2 Theta).
```

It accepts the seven BSSN shift conditions (`Frozen`, `GammaDriving2ndOrder_NoCovariant`, `GammaDriving2ndOrder_Covariant`, `GammaDriving2ndOrder_Covariant__Hatted` default, `GammaDriving1stOrder_Covariant`, `GammaDriving1stOrder_Covariant__Hatted`, `NonAdvectingGammaDriving`) with default `eta = 2`, passing the complete fCCZ4 `Lambdatilde` RHS so Brown's BSSN-only connection term never enters:

```text
Hatted 2nd-order:  partial_t beta^i = B^i + beta^j Dhat_j beta^i,
                   partial_t B^i    = beta^j Dhat_j B^i
                                      + (3/4)(partial_t Lambdatilde^i - beta^j Dhat_j Lambdatilde^i)
                                      - eta B^i,
Nonadvecting:      partial_t beta^i = B^i,
                   partial_t B^i    = (3/4) partial_t Lambdatilde^i - eta B^i.
```

When `enable_YBS_Gamma_constraint_adjustment=True`, the three second-order drivers and `NonAdvectingGammaDriving` read the YBS-adjusted `Lambdatilde` RHS once; `Frozen` and both first-order drivers are unchanged.

Five explicit paper-convention resolutions govern NRPy's equations:
1. `W = exp(-2 phi)` and `chi = exp(-4 phi)` match NRPy BSSN rather than sign-reversed prose in Mewes et al. Eqs. (33)-(34).
2. Advective `OnePlusLog` retains the lapse factor `-2 alpha (K - 2 Theta)` (Alic et al. Eq. (20); nonadvective Sanchis-Gual et al. Eq. (2.25)), treating the missing `alpha` in Mewes et al. Eqs. (35)-(37) as typographical.
3. `kappa1` terms omit an extra `alpha`, matching Mewes et al.'s displayed evolution system after Eq. (32).
4. `H_Z4` excludes `-2 Theta K`; `partial_t K` and `partial_t Theta` append `-2 Theta K` explicitly.
5. `partial_0 = partial_t - L_beta` in Mewes et al. already stretches `Lambdatilde^i = DeltaGamma^i + C^i`, so NRPy drops their duplicate `-C^k Dhat_k beta^i` bracket term and keeps `+2 C^i Dhat_k beta^k / 3` (corroborated by Alic et al. Eq. (19); Sanchis-Gual et al. Eqs. (2.11)/(2.17) vs spherical Eq. (3.14) differ off `C^i = 0` in the divergence coefficient, while Brown Eqs. (21e)/(22e) supply the BSSN-only term omitted here).

Evidence (descriptive behavior): `fCCZ4_gauge_RHSs.py::fCCZ4_gauge_RHSs`; `BSSN_gauge_RHSs.py::BSSN_gauge_RHSs`; [Mewes et al.](https://arxiv.org/pdf/2002.06225v2) Eqs. (25)-(41); [Alic et al.](https://arxiv.org/pdf/1106.2254v2) Eqs. (19)-(20).

## Sources

- [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2) - Eqs. (3)-(41) and unnumbered evolution system after Eq. (32)
- [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1) - Eqs. (2.1)-(2.29) and spherical Eq. (3.14)
- [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2) - Sec. II, Eqs. (19)-(20), covariance discussion after Eq. (23)
- [Brown, arXiv:0902.3652v2](https://arxiv.org/pdf/0902.3652v2) - Eqs. (12a), (12b), (15), (21e), (22e)
- [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py) - `FCCZ4Constraints`
- [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py) - `FCCZ4RHSs`
- [fCCZ4_gauge_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py) - `fCCZ4_gauge_RHSs`
- [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py) - `BSSNQuantities`
- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - `BSSN_gauge_RHSs`
- [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints`
- [BSSN_algebraic_constraints.py](../../../nrpy/equations/general_relativity/BSSN_algebraic_constraints.py) - `BSSN_algebraic_constraints`

## See Also

- Depends on: [BSSN Family](bssn-family.md)
- Depends on: [Reference Metrics](../../core/reference-metrics.md)
- Implements: [fCCZ4 Implementation, Assembly, And Validation](fccz4-implementation-and-validation.md)
