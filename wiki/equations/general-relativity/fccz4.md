# Fully Covariant Conformal Z4

> Define the three-dimensional reference-metric fCCZ4 equations implemented by NRPy, including algebraic constraints, gauge, and validation boundaries. · Status: confirmed · Last reconciled: 08-28-2026
> Up: [General Relativity](index.md)

## Summary

NRPy's fCCZ4 implementation expresses the fully covariant conformal Z4
system as corrections to option-matched cached reference-metric BSSN
expressions, read without mutation. A canonical constraint module owns the
connection constraint, spatial Z4 vectors, Z4 Ricci tensor, and Hamiltonian
expression; the evolution module consumes that cache and exposes 18 RHS
components while evolving the normal Z4 projection in storage named
`Theta_fCCZ4`. Existing conformal
connection storage represents `LambdatildeU`. A separate gauge module reuses
the BSSN gauge implementation and adds only the fCCZ4 lapse and Gamma-driver
corrections.

The scientific contract follows the three-dimensional formulation in Mewes et
al. Eqs. (3)-(41), except for an apparent convention error in the printed
connection shift sector. Mewes et al. define their time operator using the full
vector Lie derivative, which already stretches the evolved connection, and
then print the constraint-vector stretch a second time. NRPy uses the
first-principles one-stretch composition derived below. Alic et al. corroborate
the full corrected Cartesian coefficient; Sanchis-Gual et al. corroborate the
one-stretch structure, with their differing off-constraint divergence
presentations qualified below. Brown supplies the covariant BSSN
reference-connection construction. Spatial dimension is
three, spacetime signature is `(-,+,+,+)`, exact rational coefficients are
preserved, and full covariance fixes `kappa3=1` rather than exposing it as a
parameter.

## Detail

### Conventions and matter projections

The 3+1 metric and future-pointing normal are

```text
ds^2 = -alpha^2 dt^2
       + gamma_ij (dx^i + beta^i dt)(dx^j + beta^j dt),
n^mu = alpha^(-1) (1, -beta^i).
```

Matter projections are

```text
rho  = n_mu n_nu T^(mu nu),
S_i  = -gamma_(i mu) n_nu T^(mu nu),
S_ij = gamma_(i mu) gamma_(j nu) T^(mu nu),
S    = gamma^ij S_ij.
```

Claim evidence:
- Claim: The adopted three-dimensional fCCZ4 convention uses signature `(-,+,+,+)`, the displayed 3+1 metric and normal, and the displayed matter projections, with `rho` corresponding to Mewes et al.'s `E`.
- Role: public/scientific contract
- Deciding authority: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (3)-(5) and (38)-(41); targeted validation: [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), module `__main__` matter case, and [fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py), `trusted_dict`
- Corroboration: [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), Eqs. (2.1)-(2.29)

Here `rho` is Mewes et al.'s `E`. The implementation uses the damping
convention displayed by Mewes et al.: `kappa1` terms carry no additional lapse
factor. `kappa1` has inverse-length units and defaults to `0.1`; dimensionless
`kappa2` defaults to zero.

Claim evidence:
- Claim: `FCCZ4RHSs` fixes full covariance instead of registering `kappa3`, registers inverse-length `kappa1` with default `0.1` and dimensionless `kappa2` with default zero, and constructs `kappa1` damping terms with no additional lapse factor.
- Role: descriptive behavior
- Deciding authority: [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__`
- Corroboration: [fCCZ4_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_Cartesian.py), `trusted_dict`; [fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py), `trusted_dict`; [fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py), `trusted_dict`
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=30-significant-digit deterministic trusted sampling of final evolution outputs; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=six evolution comparisons across SinhCartesian, SinhSpherical, and representative Cartesian coordinates; date=08-28-2026`

### Conformal and reference-metric variables

The conformal split is

```text
gamma_ij = exp(4 phi) gammabar_ij,
Abar_ij  = exp(-4 phi) (K_ij - gamma_ij K/3),
K        = gamma^ij K_ij.
```

The Lagrangian determinant condition is
`partial_t det(gammabar)=0`. Initializing
`det(gammabar)=det(gammahat)` therefore preserves that equality.

NRPy backend kernels can impose both algebraic conformal constraints at a
caller- or schedule-selected state-repair point. For raw tensors
`G_ij=gammabar_ij` and `A_ij=Abar_ij`, define

```text
q          = (abs(det(gammahat)) / det(G))^(1/3),
G'_ij      = q G_ij,
tau'       = (G')^ij A_ij,
A'_ij      = A_ij - G'_ij tau'/3.
```

For nondegenerate positive-definite spatial `G` and `gammahat`, one loop forms
`G'`, inverts that corrected metric, then forms `A'`; all raw components are
loaded before any corrected component is stored. Since
`abs(det(gammahat))=det(gammahat)` in this domain,
`det(G')=det(gammahat)` and `(G')^ij A'_ij=0` are imposed together, with the
trace computed using the post-enforced metric. This constrained-state contract
applies mathematically to both BSSN and fCCZ4 storage; no `trAbar` damping
variable or RHS correction is retained. Reviewed owner paths and collision
examples use BSSN application wiring and do not establish an end-to-end fCCZ4
projection lifecycle for the new RHS/gauge owners and `Theta_fCCZ4`.
No guard or behavior for singular, indefinite, or non-positive-determinant
metric input is established.

Claim evidence:
- Claim: The conformal variables obey the displayed metric and trace-free extrinsic-curvature split; the source-specified Lagrangian determinant condition preserves `det(gammabar)=det(gammahat)` when imposed initially, but fCCZ4 owner validation does not directly exercise that determinant-preservation statement.
- Role: public/scientific contract
- Deciding authority: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (7)-(10); targeted validation of the exercised conformal variables: [fCCZ4_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_Cartesian.py), `trusted_dict`
- Corroboration: [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), conformal/reference-metric definitions in Eqs. (2.1)-(2.29)

Claim evidence:
- Claim: For nondegenerate positive-definite spatial metrics, BHaH, ETLegacy, and CarpetX provide one algebraic projection that imposes `det(gammabar)=det(gammahat)` and `tr(Abar)=0` using the inverse of the already determinant-corrected metric; the kernel is compatible with shared BSSN/fCCZ4 tensor storage, but reviewed owner paths and collision examples do not establish an end-to-end lifecycle for the new fCCZ4 evolution and `Theta_fCCZ4` storage.
- Role: public/scientific contract
- Deciding authority: backend `register_CFunction_enforce_detgbar_equals_detghat_trAzero` implementations listed in [BSSN Family](bssn-family.md); [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), constrained metric RHS
- Corroboration: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (8)-(10)
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=BHaH, ETLegacy, CarpetX generated C; precision=exact symbolic off-diagonal determinant/trace identities and generated-source structure; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=Cartesian and SinhSpherical reference-metric precompute in all three backends, end-to-end fCCZ4 application wiring not established by inspected paths; date=08-26-2026`

NRPy supports evolving `phi`, `W=exp(-2 phi)`, or `chi=exp(-4 phi)`.

Claim evidence:
- Claim: `FCCZ4RHSs` supports the existing NRPy BSSN conformal-factor options `phi`, `W=exp(-2 phi)`, and `chi=exp(-4 phi)` and rejects other values through the reused BSSN quantities interface.
- Role: descriptive behavior
- Deciding authority: [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__`; [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__`; [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__`
- Corroboration: [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), module `__main__`, exercises `W`, `phi`, and `chi` in both Cartesian and SinhSpherical coordinates
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=30-significant-digit deterministic trusted sampling of final Hamiltonian output; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=W, phi, and chi constraints in SinhCartesian and SinhSpherical coordinates plus representative Cartesian W; date=08-28-2026`

For time-independent reference metric `gammahat_ij`, let `Dhat_i` and
`Gammahat^i_jk` denote its derivative and connection. Barred quantities belong
to `gammabar_ij`. Define

```text
DeltaGamma^i_jk = Gammabar^i_jk - Gammahat^i_jk,
DeltaGamma^i    = gammabar^jk DeltaGamma^i_jk,
Lambdatilde^i  = DeltaGamma^i + 2 gammabar^ij Z_j.
```

### Connection constraint and spatial Z4 vector

The normal projection and conformal connection constraint are

```text
Theta = -n_mu Z^mu,
C^i   = Lambdatilde^i - DeltaGamma^i.
```

The exposed spatial quantities are

```text
Zbar^i = C^i / 2,
Z^i    = exp(-4 phi) C^i / 2,
Z_i    = gammabar_ij C^j / 2.
```

Physical lowering cancels the conformal factor in `Z^i`; lowering must use the
full, potentially off-diagonal conformal metric.

### Z4 Ricci tensor

Here `D_i` is compatible with the physical spatial metric `gamma_ij`, and
`RbarGeom_ij` is the geometric Ricci tensor of `gammabar_ij`. The source
definition is

```text
RbarZ4_ij = RbarGeom_ij + D_i Z_j + D_j Z_i.
```

After absorbing `Z_i` into `Lambdatilde^i`, Mewes et al. give the equivalent
derivative-free form

```text
RbarZ4_ij = RbarNRPy_ij(Lambdatilde)
            - 4 (Z_i partial_j phi + Z_j partial_i phi)
            + [gamma_ki (Gamma^k_jl - Gammahat^k_jl)
               + gamma_kj (Gamma^k_il - Gammahat^k_il)] Z^l.
```

The `gamma` and `Gamma` in the last term are physical. Expanding the physical
connection produces the correction

```text
RbarZ4_ij = RbarNRPy_ij(Lambdatilde) + deltaR_ij,

deltaR_ij = -2 (C_i partial_j phi + C_j partial_i phi)
            + 2 gammabar_ij C^l partial_l phi
            + (C^l / 2)
              (gammabar_ki DeltaGamma^k_jl
               + gammabar_kj DeltaGamma^k_il),
C_i       = gammabar_ij C^j.
```

Set `RbarZ4=gammabar^ij RbarZ4_ij`. The pure-trace term
`2 gammabar_ij C^l partial_l phi` is required by the physical connection;
omitting it changes both the `K` and `Theta` equations.

Claim evidence:
- Claim: In three-dimensional, Lagrangian reference-metric fCCZ4 with `kappa3=1`, the evolved connection, spatial Z4 vector, and absorbed Ricci tensor obey the definitions and `deltaR_ij` expansion displayed above.
- Role: public/scientific contract
- Deciding authority: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (25), (26), (29), (30), and (32), and [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__`
- Corroboration: [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), Eqs. (2.1)-(2.29); [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Sec. II, Eq. (19) and the covariance discussion following Eq. (23); [Brown, arXiv:0902.3652v2](https://arxiv.org/pdf/0902.3652v2), Eqs. (12a), (12b), and (15)

NRPy stores `Lambdatilde^i` in the existing `LambdabarU`/`lambdaU` fields.
`FCCZ4Constraints` constructs every component of `C^i` before lowering with
the full conformal metric. `BSSNQuantities` uses evolved `lambdaU` in the
linear Ricci derivative slot and geometric connection differences in nonlinear
terms.

Claim evidence:
- Claim: NRPy maps `Lambdatilde^i` to existing `LambdabarU`/`lambdaU` storage, completes `C^i` before full-metric lowering, and uses the evolved-linear/geometric-nonlinear Ricci split.
- Role: descriptive behavior
- Deciding authority: [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__`; [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__`
- Corroboration: [fCCZ4_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_Cartesian.py), `trusted_dict`; [fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py), `trusted_dict`
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=source inspection plus 30-significant-digit deterministic sampling of final Hamiltonian output; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=ten constraint comparisons; external producer and independent scientific implementation not verified; date=08-28-2026`

### Constraints

The module stores

```text
H_Z4 = exp(-4 phi)
       [RbarZ4 - 8 Dbar^i phi Dbar_i phi - 8 Dbar^i Dbar_i phi]
       + 2 K^2 / 3 - Abar_ij Abar^ij - 16 pi rho.
```

`H_Z4` does not include `-2 Theta K`; the evolution equations append that term
where needed. The differential connection constraint is `C^i=0`. At the
equation level, fCCZ4 corrections vanish and the implemented RHSs coincide with
their BSSN bases when `Theta=C^i=H_Z4=Abar^i_i=0`, with derivatives of the
identically zero fields also zero. This is not a complete physical constraint
surface: a fully constrained solution additionally satisfies the momentum
constraint `M^i=0` and `det(gammabar)=det(gammahat)`.

Claim evidence:
- Claim: `H_Z4` is the displayed matter-inclusive fCCZ4 Hamiltonian expression and excludes `-2 Theta K`; the implemented fCCZ4 corrections reduce to their BSSN base equations when `Theta=C^i=H_Z4=Abar^i_i=0` and derivatives of the identically zero fields vanish, while a complete physical constraint surface also requires `M^i=0` and `det(gammabar)=det(gammahat)`.
- Role: public/scientific contract
- Deciding authority: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (24)-(32) and the evolution system following Eq. (32); [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__`
- Corroboration: [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), constraints and BSSN reduction in Eqs. (2.1)-(2.29)

`FCCZ4Constraints` owns this expression and the option-matched
`fCCZ4_constraints` cache. It adds the matter source only when
`enable_T4munu=True`; `FCCZ4RHSs` consumes and publicly retains `H_Z4`.

Claim evidence:
- Claim: `FCCZ4Constraints` and `fCCZ4_constraints` own the canonical `H_Z4` expression and option-matched cache, add its matter source only when `enable_T4munu=True`, and provide it for `FCCZ4RHSs` to consume and publicly retain.
- Role: descriptive behavior
- Deciding authority: [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__`, `FCCZ4ConstraintsDict`, and `fCCZ4_constraints`; [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__`
- Corroboration: [fCCZ4_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_Cartesian.py), `trusted_dict`; [fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py), `trusted_dict`
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=30-significant-digit deterministic trusted sampling of final Hamiltonian output; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=ten constraint comparisons across W, phi, chi, external Ricci, precompute, matter, and three coordinate families; date=08-28-2026`

`BSSNconstraints.MU` constructs the standard conformal momentum expression:

```text
M^i = exp(-4 phi)
      [Dbar_j Abar^ij + 6 Abar^ij partial_j phi
       - 2 gammabar^ij partial_j K / 3]
      - 8 pi S^i = 0.
```

`FCCZ4RHSs` does not expose a momentum diagnostic or instantiate
`BSSNconstraints`. Consumers needing the diagnostic can construct that object
and read its unchanged `MU` expression.

Claim evidence:
- Claim: `BSSNconstraints.MU` constructs the displayed standard conformal momentum expression; `FCCZ4RHSs` neither exposes it nor instantiates `BSSNconstraints`, while consumers can obtain the diagnostic by constructing `BSSNconstraints` and reading `MU`.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints.__init__`; [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__`
- Corroboration: [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), constraint system in Eqs. (2.1)-(2.29)
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=exact source inspection; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=vacuum and matter branches inspected; date=08-26-2026`

### Evolution equations

Let `D_beta=Dbar_k beta^k`. The coordinate-time equations are

```text
partial_t gammabar_ij = beta^k partial_k gammabar_ij
                        + gammabar_ik partial_j beta^k
                        + gammabar_kj partial_i beta^k
                        - 2 gammabar_ij D_beta / 3
                        - 2 alpha Abar_ij,

partial_t phi = beta^k partial_k phi + (D_beta - alpha K) / 6,

partial_t Abar_ij = beta^k partial_k Abar_ij
                    + Abar_ik partial_j beta^k
                    + Abar_kj partial_i beta^k
                    - 2 Abar_ij D_beta / 3
                    - 2 alpha Abar_ik Abar^k_j
                    + alpha Abar_ij (K - 2 Theta)
                    + exp(-4 phi)
                      [-2 alpha Dbar_i Dbar_j phi
                       + 4 alpha Dbar_i phi Dbar_j phi
                       + 2(Dbar_i alpha Dbar_j phi
                           + Dbar_j alpha Dbar_i phi)
                       - Dbar_i Dbar_j alpha
                       + alpha(RbarZ4_ij - 8 pi S_ij)]^TF,

partial_t K = beta^i partial_i K
              + exp(-4 phi)
                [alpha(RbarZ4
                       - 8 Dbar^i phi Dbar_i phi
                       - 8 Dbar^i Dbar_i phi)
                 - (2 Dbar^i alpha Dbar_i phi
                    + Dbar^i Dbar_i alpha)]
              + alpha(K^2 - 2 Theta K)
              - 3 kappa1(1 + kappa2) Theta
              + 4 pi alpha(S - 3 rho),

partial_t Theta = beta^i partial_i Theta
                  + alpha(H_Z4 - 2 Theta K) / 2
                  - Z^i partial_i alpha
                  - kappa1(2 + kappa2) Theta.
```

For alternate conformal factors, the implementation applies the exact chain
rule:

```text
partial_t W   = -2 W partial_t phi,
partial_t chi = -4 chi partial_t phi.
```

The evolved conformal connection uses a BSSN-shaped base

```text
LambdaBase^i = beta^k partial_k Lambdatilde^i
      - Lambdatilde^k partial_k beta^i
      + gammabar^jk Dhat_j Dhat_k beta^i
      + 2 DeltaGamma^i D_beta / 3
      + gammabar^ij Dbar_j D_beta / 3
      - 2 Abar^ij(partial_j alpha - 6 alpha partial_j phi)
      + 2 alpha Abar^jk DeltaGamma^i_jk
      - 4 alpha gammabar^ij partial_j K / 3
      - 16 pi alpha gammabar^ij S_j,
```

plus the fCCZ4 correction

```text
deltaLambda^i = 2 gammabar^ij
                (alpha partial_j Theta - Theta partial_j alpha)
                - 2 alpha K C^i / 3 - kappa1 C^i
                + 2 C^i Dhat_k beta^k / 3,

partial_t Lambdatilde^i = LambdaBase^i + deltaLambda^i.
```

This composition has exactly one vector stretch. The base advection-stretch
pair is the full vector Lie pair, which may be written equivalently as

```text
beta^k partial_k Lambdatilde^i
  - Lambdatilde^k partial_k beta^i
= beta^k Dhat_k Lambdatilde^i
  - Lambdatilde^k Dhat_k beta^i.
```

Since `Lambdatilde^i=DeltaGamma^i+C^i`, the covariant form already contains
the constraint-vector stretch `-C^k Dhat_k beta^i`. Adding the same term from
Mewes et al.'s printed `kappa3=1` bracket would count that contribution twice.
The retained `+2 C^i Dhat_k beta^k/3` is required:
the BSSN base divergence uses geometric `DeltaGamma^i`, so this correction
promotes its coefficient to evolved `Lambdatilde^i`. In Cartesian coordinates,
for shift-gradient component `partial_b beta^a`, the aggregate coefficient is

```text
-Lambdatilde^b delta^i_a + (2/3) Lambdatilde^i delta^a_b.
```

Claim evidence:
- Claim: The displayed `gammabar_ij`, `phi`, `Abar_ij`, `K`, `Theta`, and `Lambdatilde^i` equations define the adopted three-dimensional fCCZ4 evolution system; the connection equation contains exactly one evolved-vector stretch and retains the constraint-vector divergence promotion, while `H_Z4` excludes `-2 Theta K`, `kappa1` damping is lapse-unscaled, and the `W` and `chi` chain rules are exact.
- Role: public/scientific contract
- Deciding authority: the definition `Lambdatilde^i=DeltaGamma^i+C^i` and full-vector-Lie time operator in [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (29) and the evolution system after Eq. (32), composed with [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__`, and [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__` and module `__main__`
- Corroboration: [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (19), corroborates the full corrected Cartesian coefficient; [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), Eqs. (2.11) and (2.17), corroborates one stretch but not the promoted-divergence coefficient; [fCCZ4_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_Cartesian.py), `trusted_dict`; [fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py), `trusted_dict`

### Reuse of BSSN expressions

The implementation retrieves parameter-aware cached `BSSNQuantities` and
`BSSNRHSs` objects without mutation. `FCCZ4Constraints` constructs the
fCCZ4-only constraint, Z4-vector, Ricci, and Hamiltonian aggregates;
`FCCZ4ConstraintsDict` and the module-level `fCCZ4_constraints` object cache
them.
`FCCZ4RHSs` retrieves the option-matched constraint object, copies its mutable
tensor containers, and realizes fCCZ4 as these corrections:

```text
(partial_t h_ij)_fCCZ4 = (partial_t h_ij)_BSSN
                         [with tr(Abar)=0 imposed analytically],

(partial_t a_ij)_fCCZ4 = (partial_t a_ij)_BSSN
                         + [-2 alpha Theta Abar_ij
                            + alpha exp(-4 phi) (deltaR_ij)^TF] / ReDD_ij,

(partial_t K)_fCCZ4 = (partial_t K)_BSSN
                      + alpha(H_Z4 - 2 Theta K)
                      - 3 kappa1(1 + kappa2) Theta.
```

`ReDD_ij` is the componentwise reference-metric rescaling factor, not a tensor
contraction. Removing the former BSSN off-constraint trace term and the
compensating fCCZ4 subtraction is a paired change: the fCCZ4 metric RHS remains
the displayed source equation. The conformal-factor RHS is reused unchanged.
The connection RHS adds `deltaLambda^i` before vector rescaling.
When matter is enabled, the base BSSN source terms and the canonical `H_Z4`
matter term use the established `T4munu` source helpers.

Claim evidence:
- Claim: `FCCZ4Constraints` constructs the fCCZ4-only connection, Z4 Ricci, and `H_Z4` aggregates; `FCCZ4ConstraintsDict` and `fCCZ4_constraints` cache them; `FCCZ4RHSs` retrieves option-matched cached BSSN and fCCZ4 constraint objects, copies mutable aggregates without mutation, applies the displayed corrections at reference-rescaled output boundaries, and delegates enabled matter sources to the established `T4munu` helpers.
- Role: descriptive behavior
- Deciding authority: [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__` and `FCCZ4ConstraintsDict`; [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__`
- Corroboration: [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities`; [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs`; [T4munu.py](../../../nrpy/equations/general_relativity/T4munu.py), `BSSN_RHSs_T4UU_source_terms` and `BSSN_constraints_T4UU_source_terms`
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=30-significant-digit deterministic trusted sampling of 18 final evolution outputs; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=six comparisons, prioritizing SinhCartesian and SinhSpherical with one representative Cartesian case; date=08-28-2026`

### Gauge equations

The wrapper accepts the four lapse names supported by BSSN:
`OnePlusLog`, `BHSHarmonicSlicing`, `Frozen`, and `OnePlusLogAlt`. For
`OnePlusLog`, fCCZ4 changes the delegated lapse equation to

```text
partial_t alpha = beta^i partial_i alpha - 2 alpha(K - 2 Theta).
```

It accepts seven shift names: `Frozen`,
`GammaDriving2ndOrder_NoCovariant`, `GammaDriving2ndOrder_Covariant`,
`GammaDriving2ndOrder_Covariant__Hatted`,
`GammaDriving1stOrder_Covariant`,
`GammaDriving1stOrder_Covariant__Hatted`, and
`NonAdvectingGammaDriving`. The hatted covariant second-order Gamma driver is

```text
partial_t beta^i = B^i + beta^j Dhat_j beta^i,
partial_t B^i    = beta^j Dhat_j B^i
                   + 3(partial_t Lambdatilde^i
                       - beta^j Dhat_j Lambdatilde^i) / 4
                   - eta B^i,
```

and the nonadvecting Gamma driver,

```text
partial_t beta^i = B^i,
partial_t B^i    = 3 partial_t Lambdatilde^i / 4 - eta B^i.
```

The default is the hatted covariant second-order driver, and `eta=2` is reused
from BSSN. The wrapper delegates all named choices to `BSSN_gauge_RHSs`. Its
fCCZ4 additions are `4 alpha Theta` for `OnePlusLog` and
`3 deltaLambda^i/4` in the driver RHS for all three second-order drivers and
the nonadvecting driver, before standard vector rescaling. The two first-order
drivers and the other three lapse choices receive no local fCCZ4 addition.

Claim evidence:
- Claim: `fCCZ4_gauge_RHSs` accepts the exact four lapse and seven shift names listed above, defaults to `OnePlusLog` with `GammaDriving2ndOrder_Covariant__Hatted`, reuses the BSSN default `eta=2`, and applies the listed fCCZ4 additions only to `OnePlusLog` and the three second-order plus nonadvecting driver RHSs.
- Role: descriptive behavior
- Deciding authority: [fCCZ4_gauge_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py), `fCCZ4_gauge_RHSs`; [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py), `BSSN_gauge_RHSs`
- Corroboration: [fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_SinhSpherical_rfm_precompute_T4munu.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_SinhSpherical_rfm_precompute_T4munu.py), `trusted_dict`; [fCCZ4_gauge_RHSs_Frozen_Frozen_SinhSpherical.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_Frozen_Frozen_SinhSpherical.py), `trusted_dict`; [fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_Cartesian.py), `trusted_dict`
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=30-significant-digit deterministic trusted sampling of seven final gauge outputs; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=full 4-lapse by 7-shift matrices in SinhCartesian and SinhSpherical plus two representative Cartesian/precompute-matter comparisons; date=08-28-2026`

### NRPy mapping and public outputs

NRPy evolves reference-rescaled components:

```text
gammabar_ij     = gammahat_ij + ReDD_ij hDD_ij,
Abar_ij         = ReDD_ij aDD_ij,
Lambdatilde^i   = ReU_i lambdaU_i,
beta^i          = ReU_i vetU_i,
B^i             = ReU_i betU_i.
```

Claim evidence:
- Claim: NRPy's fCCZ4 public variables use the displayed `ReDD` and `ReU` reference-metric component mappings for the conformal metric, trace-free curvature, evolved connection, shift, and shift driver.
- Role: descriptive behavior
- Deciding authority: [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__`; [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__`; [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__`; [fCCZ4_gauge_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py), `fCCZ4_gauge_RHSs`
- Corroboration: none available; the mappings are established directly by the owner evolution and gauge code
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=30-significant-digit deterministic trusted sampling of final mapped outputs; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=six evolution and 58 gauge comparisons across SinhCartesian, SinhSpherical, and representative Cartesian coordinates; date=08-28-2026`

`FCCZ4RHSs` exposes conceptual attribute `Theta`, `Z4constraintU`, `ZbarU`, `ZU`, `ZD`,
`LambdatildeU`, `RbarZ4DD`, `RbarZ4`, `H_Z4`, rescaled metric and
extrinsic-curvature RHS arrays, and unrescaled and rescaled
conformal-connection RHS arrays. `Theta_dD` is centered; `Theta_dupD` is the upwinded advection
derivative. `fCCZ4_RHSs_varname_to_expr_dict` has 18 sorted entries: six each
for symmetric `hDD` and `aDD`, three for `lambdaU`, and one each for the
conformal factor, `K`, and `Theta_fCCZ4`. `FCCZ4RHSsDict` provides
coordinate-option caching through the module-level `fCCZ4_RHSs` object.
Storage, derivative basenames, and dictionary keys use the formulation-specific
names `Theta_fCCZ4`, `Theta_fCCZ4_dD`, `Theta_fCCZ4_dupD`, and
`Theta_fCCZ4_rhs`. The conceptual Python attributes remain `Theta`,
`Theta_dD`, `Theta_dupD`, and `Theta_rhs`; no generic storage alias is retained.

Claim evidence:
- Claim: `FCCZ4RHSs` uses centered and upwinded formulation-specific Theta derivatives and provides a sorted 18-entry final evolution dictionary containing `Theta_fCCZ4_rhs`, with no generic Theta storage alias.
- Role: descriptive behavior
- Deciding authority: [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs`, `FCCZ4RHSsDict`, and `fCCZ4_RHSs`
- Corroboration: [fCCZ4_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_Cartesian.py), `trusted_dict`; [fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py), `trusted_dict`; ordering and namespace behavior are established by owner-source inspection
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=30-significant-digit deterministic trusted sampling of 18 final outputs; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=SinhCartesian and SinhSpherical baseline/external-Ricci cases, SinhSpherical precompute with matter, and representative Cartesian; date=08-28-2026`

`FCCZ4Constraints` independently exposes the shared connection, spatial-Z4,
Ricci, and Hamiltonian aggregates through the parameter-aware module cache
`fCCZ4_constraints`. Beyond reused base BSSN gridfunctions, enabling matter
incrementally registers the required symmetric `T4UU` AUXEVOL inputs as its
only additional gridfunctions; it does not register `H`, `M`,
`LAMBDA_CONSTRAINT`, `Theta`, or other diagnostic outputs. The reused matter
helper also registers or reuses the `PI` CodeParameter.

Claim evidence:
- Claim: `FCCZ4Constraints` exposes the connection constraint, spatial Z4 vectors, Z4 Ricci tensor/scalar, and `H_Z4` through the parameter-aware `fCCZ4_constraints` cache; beyond reused base BSSN gridfunctions, enabling matter incrementally registers required symmetric `T4UU` AUXEVOL inputs as its only additional gridfunctions and no diagnostic outputs, while the reused matter helper registers or reuses the `PI` CodeParameter.
- Role: descriptive behavior
- Deciding authority: [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints`, `FCCZ4ConstraintsDict`, and `fCCZ4_constraints`; [T4munu.py](../../../nrpy/equations/general_relativity/T4munu.py), `BSSN_constraints_T4UU_source_terms`
- Corroboration: [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__` canonical-cache consumption
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=source inspection plus 30-significant-digit deterministic trusted sampling of final Hamiltonian output; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=nine curved-coordinate comparisons and one representative Cartesian comparison; date=08-28-2026`

When `enable_RbarDD_gridfunctions=True`, `FCCZ4RHSs` interprets supplied
`RbarDD` auxiliaries under the `RbarNRPy_ij(Lambdatilde)` hybrid convention:
evolved `lambdaU` occupies the linear connection-derivative slot, while
nonlinear connection terms remain geometric. It consumes those symbols and
adds `deltaR_ij`. When the six symmetric `RbarDD` components are absent,
`BSSNQuantities_dict` registers them in group `AUXEVOL` and array
`auxevol_gfs`; this creates input storage, not an external value producer.

Claim evidence:
- Claim: When `enable_RbarDD_gridfunctions=True`, `FCCZ4RHSs` interprets supplied `RbarDD` values under the evolved-linear/geometric-nonlinear hybrid Ricci convention, consumes those auxiliary symbols, and adds `deltaR_ij`; `BSSNQuantities_dict` conditionally registers six symmetric `RbarDD` components in `AUXEVOL`/`auxevol_gfs`, but neither that storage registration nor the trusted gridfunction variant constructs or verifies the external producer.
- Role: descriptive behavior
- Deciding authority: [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__` Ricci construction and gridfunction branch, and `BSSNQuantities_dict.__getitem__`/`BSSN_quantities` storage registration; [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__` Z4 Ricci correction; [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__` parameter contract and Ricci consumption
- Corroboration: [fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py), `trusted_dict`; [fCCZ4_RHSs_SinhSpherical_RbarDD_gridfunctions.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhSpherical_RbarDD_gridfunctions.py), `trusted_dict`; neither constructs the external producer
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=30-significant-digit deterministic trusted sampling of final consumer outputs; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=SinhCartesian and SinhSpherical with RbarDD gridfunctions; date=08-28-2026`

### Explicit conventions and resolved source discrepancies

Five source-resolution statements are explicit:

1. NRPy uses `W=exp(-2 phi)` and `chi=exp(-4 phi)`, matching its BSSN
   variables. Sign-reversed conformal-factor prose in Mewes et al. is not
   adopted.
2. Advective 1+log uses `-2 alpha(K-2 Theta)`, matching Alic et al. Eq. (20).
   Sanchis-Gual et al. Eq. (2.25) corroborates the lapse factor only in its
   nonadvective form; the missing lapse in the displayed Mewes et al. gauge
   equation is treated as a typographical omission.
3. Among differing published damping prescriptions, NRPy follows the displayed
   Mewes et al. evolution system, so `kappa1` terms have no extra `alpha`.
4. NRPy's `H_Z4` bookkeeping excludes `-2 Theta K`; evolution equations append
   it explicitly, preventing double counting.
5. NRPy corrects an apparent convention error in Mewes et al.'s printed
   connection shift sector. Their definition
   `partial_0=partial_t-L_beta` uses the full vector Lie derivative, so the
   base advection-stretch pair may equivalently be written with `Dhat` and
   already contains `-C^k Dhat_k beta^i`; the printed `kappa3=1` bracket adds
   that part again. NRPy omits the duplicate and retains
   `+2 C^i Dhat_k beta^k/3`, yielding exactly one evolved-vector stretch. Alic
   et al. corroborate the full corrected Cartesian coefficient. Sanchis-Gual
   et al. corroborate one stretch, but their general and spherical displayed
   divergence terms differ off the connection-constraint surface.

These statements describe the source resolutions made by the current
implementation; they do not claim that the papers use identical conventions.

Claim evidence:
- Claim: The current NRPy implementation applies the five listed fCCZ4 source resolutions where the cited sources differ or contain inconsistent equations or prose; its connection RHS corrects the duplicated Mewes constraint-vector stretch, retains the required divergence promotion, and has the exact aggregate Cartesian coefficient `-Lambdatilde^b delta^i_a+(2/3)Lambdatilde^i delta^a_b`.
- Role: public/scientific contract
- Deciding authority: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__`; [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__`; [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__` and module `__main__`; [fCCZ4_gauge_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py), `fCCZ4_gauge_RHSs`; source crosswalks in [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), and [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1)
- Corroboration: [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (19), corroborates the full corrected Cartesian coefficient; [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), Eqs. (2.11) and (2.17), corroborates one stretch but not the complete promoted-divergence coefficient; no paper is used as an oracle for the exact local coefficient identity
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=30-significant-digit deterministic trusted sampling of final evolution and gauge dictionaries after separate derivation review; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=W/phi/chi constraints and complete 4-by-7 gauge matrices in both curved coordinate families, plus representative Cartesian coverage; date=08-28-2026`

### Source crosswalk

In Mewes et al., `Theta` and `Z_i` are Eqs. (25)-(26),
`Lambdatilde^i` and `Z^i` are Eqs. (29)-(30), and the absorbed Ricci tensor is
Eq. (32). The evolution system immediately after Eq. (32) is unnumbered.
Conformal-factor alternatives are Eqs. (33)-(34), gauge equations are
Eqs. (35)-(37), and matter projections are Eqs. (38)-(41). Sanchis-Gual et
al. Eqs. (2.1)-(2.29) supply variable, constraint, gauge, covariance, and
BSSN-limit cross-checks. Their Eqs. (2.11) and (2.17) corroborate one full-Lie
stretch but leave the general divergence coefficient geometric; spherical
Eq. (3.14) contains the compatible Z4 divergence promotion. The two displayed
divergence sectors therefore are not an unqualified off-constraint cross-check
of each other. Alic et al. supplies the exact advective 1+log form in
Eq. (20) and corroborates the one-stretch Cartesian `kappa3=1` connection
composition through Sec. II, Eq. (19), and the discussion following Eq. (23).
Mewes et al.'s printed off-constraint shift bracket is not authoritative where
it conflicts with the result obtained from its own full-Lie time-operator
definition: it contains the duplicate constraint-vector stretch documented
above. Brown Eqs. (12a),
(12b), and (15) supply the tensorial
reference-connection viewpoint used by both NRPy BSSN and fCCZ4.

### BHaH spectroscopy application boundary

`blackhole_spectroscopy.py` keeps BSSN as its default and accepts `--fccz4`
as an opt-in formulation selector. That Boolean selects the RHS/gauge owners
used by `register_CFunction_rhs_eval` and asks the ADM converter to initialize
`Theta_fCCZ4=0` while constructing fresh data. `initial_data()` attempts
checkpoint loading before that fresh conversion and returns from the checkpoint
branch, so restored Theta storage is not overwritten by fresh-data zeroing.
Both formulations use the same combined determinant/trace projector after
initial-data boundary handling and through the spectroscopy Method of Lines
post-RHS hook.

This wiring establishes source-selected registration and initialization order.
It does not prove checkpoint round trips, restart correctness, long-time
stability, convergence, waveform validity, or scientific accuracy. A bounded
generated-project build and short fCCZ4 startup run, when recorded below, prove
only generation, compilation, and at least one Method of Lines step for that
specific configuration.

Claim evidence:
- Claim: BHaH spectroscopy defaults to BSSN and exposes opt-in `--fccz4`; the selected formulation controls RHS/gauge registration, fresh fCCZ4 data initializes `Theta_fCCZ4` to zero only after the checkpoint branch has declined to return, and both formulations share initial-data and post-RHS determinant/trace projection; no restart or scientific-result guarantee follows.
- Role: public/scientific contract
- Deciding authority: [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py), `parser`, `enable_fCCZ4`, initial-data/RHS/Method of Lines registrations; [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`; [initial_data.py](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py), `register_CFunction_initial_data`; [ADM_Initial_Data_Reader__BSSN_Converter.py](../../../nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py), `register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`
- Corroboration: [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, GCC 13.3.0, GNU Make 4.3; backend=BHaH OpenMP; precision=double; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-run; options=default BSSN generation/build, opt-in fCCZ4 generation/build and t_final=0.1 startup through iteration 1; date=08-28-2026`

### Validation boundary

The three equation owners contain no doctest examples but run the canonical
`doctest.testmod()` gate before trusted comparisons. Their `__main__` paths
configure 74 final-output comparisons:

- six evolution dictionaries, each with 18 mapped RHS outputs;
- ten constraint dictionaries, each with the single final `H_Z4` output; and
- 58 gauge dictionaries, each with seven final lapse/shift outputs.

Evolution covers SinhCartesian baseline and external-Ricci cases,
SinhSpherical baseline, external-Ricci, and precompute-with-matter cases, plus
one lower-priority Cartesian representative. Constraints cover `W`, `phi`, and
`chi` in both curved coordinate families, external Ricci in both, one
SinhSpherical precompute-with-matter case, and one Cartesian representative.
Gauge coverage prioritizes complete 4-by-7 lapse/shift matrices in both
SinhCartesian and SinhSpherical, then adds a Cartesian default representative
and a SinhSpherical precompute-with-matter representative. Each trusted
dictionary contains only the stated final outputs.

These deterministic samples are regression evidence for symbolic construction
under the exercised options. They are not an independent scientific
implementation and do not establish generated C/CUDA compilation, long-time
stability, convergence order, waveform validity, or scientific accuracy.

Claim evidence:
- Claim: The three equation owners configure 74 trusted final-output comparisons—six 18-output evolution dictionaries, ten one-output constraint dictionaries, and 58 seven-output gauge dictionaries—with complete SinhCartesian and SinhSpherical gauge matrices prioritized over representative Cartesian cases; this validation does not establish generated backend builds, numerical stability, convergence, waveform validity, or scientific accuracy.
- Role: descriptive behavior
- Deciding authority: [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py), module `__main__`; [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py), module `__main__`; [fCCZ4_gauge_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py), module `__main__`; [validate_expressions.py](../../../nrpy/validate_expressions/validate_expressions.py), `process_dictionary_of_expressions` and `compare_or_generate_trusted_results`
- Corroboration: [fCCZ4_constraints_SinhSpherical_phi.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_constraints_SinhSpherical_phi.py), `trusted_dict`; [fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py), `trusted_dict`; [fCCZ4_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_Cartesian.py), `trusted_dict`; [fCCZ4_gauge_RHSs_Frozen_Frozen_SinhSpherical.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_Frozen_Frozen_SinhSpherical.py), `trusted_dict`; [fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_Cartesian.py), `trusted_dict`
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=30-significant-digit deterministic trusted sampling of final-output dictionaries; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=6 evolution + 10 constraint + 58 gauge comparisons; date=08-28-2026`

## Sources

- [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2) - Eqs. (3)-(41), including the unnumbered evolution system after Eq. (32); principal reference-metric fCCZ4 source, with the documented apparent connection-sector convention error corrected locally
- [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1) - Eqs. (2.1)-(2.29), reference-metric covariance, constraints, nonadvective lapse-factor corroboration, and BSSN limit
- [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2) - Sec. II, Eqs. (19)-(20) and the covariance discussion following Eq. (23)
- [Brown, arXiv:0902.3652v2](https://arxiv.org/pdf/0902.3652v2) - Eqs. (12a), (12b), and (15), covariant BSSN reference-connection construction
- [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints.MU`, unchanged momentum-constraint diagnostic
- [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py) - `BSSNQuantities`, inherited conformal variables and reference-metric quantities
- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`, inherited BSSN evolution expressions
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - `BSSN_gauge_RHSs`, delegated gauge equations and `eta` default
- [T4munu.py](../../../nrpy/equations/general_relativity/T4munu.py) - `BSSN_RHSs_T4UU_source_terms`, `BSSN_constraints_T4UU_source_terms`
- [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py) - `FCCZ4Constraints`, `FCCZ4ConstraintsDict`, `fCCZ4_constraints`
- [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py) - `FCCZ4RHSs`, `FCCZ4RHSsDict`, `fCCZ4_RHSs`
- [fCCZ4_gauge_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py) - `fCCZ4_gauge_RHSs`
- [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py) - `--fccz4`, formulation-selecting registrations, projection hooks
- [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`
- [initial_data.py](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py) - checkpoint-first `register_CFunction_initial_data`
- [ADM_Initial_Data_Reader__BSSN_Converter.py](../../../nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py) - fresh-data `Theta_fCCZ4` initialization
- [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - shared determinant/trace projector
- [fCCZ4_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_Cartesian.py) - `trusted_dict`
- [fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhCartesian_RbarDD_gridfunctions.py) - `trusted_dict`
- [fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhSpherical_rfm_precompute_T4munu.py) - `trusted_dict`
- [fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_SinhSpherical_rfm_precompute_T4munu.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_SinhSpherical_rfm_precompute_T4munu.py) - `trusted_dict`
- [fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_Cartesian.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted_Cartesian.py) - `trusted_dict`
- [fCCZ4_constraints_SinhSpherical_phi.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_constraints_SinhSpherical_phi.py) - representative curved alternate-conformal-factor `trusted_dict`
- [fCCZ4_RHSs_SinhSpherical_RbarDD_gridfunctions.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_RHSs_SinhSpherical_RbarDD_gridfunctions.py) - curved external-Ricci consumer `trusted_dict`
- [fCCZ4_gauge_RHSs_Frozen_Frozen_SinhSpherical.py](../../../nrpy/equations/general_relativity/tests/fCCZ4_gauge_RHSs_Frozen_Frozen_SinhSpherical.py) - representative curved frozen-gauge `trusted_dict`
- [validate_expressions.py](../../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`

## See Also

- Parent: [General Relativity](index.md)
- Depends on: [BSSN Family](bssn-family.md)
- Depends on: [Reference Metrics](../../core/reference-metrics.md)
- Validated by: [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
- See also: [Metric Conversions And Matter](metric-conversions-and-matter.md)
