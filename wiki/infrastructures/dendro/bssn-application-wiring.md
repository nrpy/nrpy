# BSSN Application Wiring

> Describe BSSN equations, initial data, conversions, and runtime services in `NRPy_BSSN_GR`. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_bssn` emits a complete BSSN application with a
code-generation choice of W (the default) or chi. It uses centered
derivatives, SSL, constant default `eta=1`, Dendro-style CAHD, Brown's
covariant Lambda adjustment, and separate conformal Ricci and RHS kernels.
Runtime parameter input can override the generated eta default, but cannot
change the evolved conformal factor of an already generated binary.

## Detail

The BSSN CAHD contribution depends on the evolved conformal factor:

```text
W_rhs += (C_CAHD/2) W H dx^2 / (dt (1 + 10 dx^2)),  default C_CAHD = 0.06.
chi_rhs += C_CAHD chi H dx^2 / (dt (1 + 10 dx^2)).
```

`C_CAHD` is a runtime parameter; `dx` and `dt` are the current block spacing
and timestep. The factor of two
follows `chi=W^2`; it preserves the same CAHD perturbation to the physical
conformal factor. SSL uses `W=cf` or `W=sqrt(cf)` in its relaxation toward
`alpha=W`. KO has no W-dependent scaling in this profile (`enable_CAKO=False`).
Brown's adjustment appears
exactly once in the covariant `lambdaU` RHS. The generated constraint kernels
are named `BSSN_constraints_order_N`; there is no generic constraint catch-all
source.

Claim evidence:
- Claim: The W and chi BSSN variants receive the corresponding CAHD factor; SSL uses W in either variant, while KO is not scaled by W. The conformal representation is a code-generation choice, not a runtime switch.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`, `register_CFunction_rhs_eval`; `nrpy/examples/dendro_bssn.py`, `parse_args` and `main`.
- Corroboration: `nrpy/infrastructures/Dendro/general_relativity/ADM_to_BSSN.py`, `register_CFunction_ADM_to_BSSN`, uses the selected conformal factor in initial-data conversion.

The initial octree is seeded with Dendro-GR's analytic
`punctureDataPhysicalCoord` prescription, adapted into the generated entry
point with its MIT provenance retained. This is a grid-construction seed,
not the evolved initial data. Single-rank `--tpid` mode solves for the puncture
correction `u`, then writes reusable spectral coefficients. Fresh evolution
loads them and interpolates ADM data using full `psi=psi_background+u` before
initial-grid remeshing. Initialization
sets `alpha=W=psi^(-2)` even when the evolved field is `chi=W^2`, and follows
this sequence on the initial grid and again after an initial-grid remesh:

```text
TwoPunctures -> ADM_to_BSSN -> zip/exchange/unzip -> physical boundary
             -> initial_data_lambdaU -> zip -> owned-node floor/projection
```

Claim evidence:
- Claim: Fresh evolution loads a precomputed TwoPunctures solution, converts its ADM data, computes the initial connection field from exchanged block data, then zips and floors/projects owned evolved nodes.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`; `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::initialize` within `output_solver_context_cpp`.
- Corroboration: `nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py` and `enforce_detgbar_equals_detghat_trAzero.py`, owned-node CFunction signatures.

Apparent-horizon searches interpolate selected evolved BSSN fields and then
convert each search point to ADM variables. `BSSN_to_ADM` reconstructs grid
fields used by ADM surface quantities; waveform extraction evaluates Psi4
directly from the selected BSSN fields. Generated runtime services also provide
physical boundaries, conformal-factor volume-weighted and unique-node constraint diagnostics, checkpoint and
restart, and scheduled output.

After one evolved-state halo exchange and physical-boundary fill, solver
context traverses blocks once. For each block it calls `Ricci_eval_order_N`
immediately followed by `rhs_eval_order_N`. No synchronization or second block
traversal occurs between them. Algebraic projection runs after initialization,
each RK stage before the next exchange, and AMR transfer.
At the same three points, `alpha` is floored at `CHI_FLOOR` and the evolved
conformal factor at `sqrt(CHI_FLOOR)` for W or `CHI_FLOOR` for chi. Thus both
representations use the same physical chi floor. The chi variant has distinct
checkpoint formulation metadata, preventing an incompatible W/chi restart.

The Dendro connector writes covariant momentum components to its fixed `MU`
diagnostic slots, followed by scalar `M_CONSTRAINT` and
`LAMBDA_CONSTRAINT` magnitudes. [Constraints And Diagnostic
Norms](constraints-and-diagnostic-norms.md) gives the contractions and explains
why comparable plots require the same norm, excision, remesh state, and
physical time.

The same parameter file alone does not make generated and native BSSN
evolutions equivalent:

| Choice | Generated `NRPy_BSSN_GR` | Native CPU `BSSN_GR` |
| --- | --- | --- |
| Evolved conformal factor | W by default; `--conformal-factor chi` selects chi at generation | chi |
| TwoPunctures initial lapse | Always `alpha=(psi_background+u)^(-2)=W` | `TPID_REPLACE_LAPSE_WITH_SQRT_CHI=true` is needed to replace the ordinary `INITIAL_LAPSE=2` result with full-psi W |
| Shift-driver damping in the SSL/CAHD path | Spatially constant `eta`; default 1, overridden by `ETA_CONST` in an input file | Radial RIT profile in the CPU RHS; `ETA_CONST` does not select constant damping on this path |
| KO strength | Both generated strengths read `KO_DISS_SIGMA`; the emitted sample sets 0.4 when KO is enabled | Reads `KO_DISS_SIGMA` with CAKO off; when CAKO is enabled, uses `sqrt(chi)` times separate gauge and other CAKO coefficients |

Thus matching the conformal representation, full-psi lapse, and KO parameter
still leaves a gauge difference when native uses its radial eta profile. The
generated BSSN profile does not offer that profile. Native CAKO must also be
off to compare the same base KO parameter; it can be enabled at startup or
after merger. Matching AMR controls and
diagnostic norms also does not prove identical mesh histories or evolved
fields.

Claim evidence:
- Claim: Native CPU SSL/CAHD evolution uses radial RIT eta even when `ETA_CONST` is present, whereas generated BSSN uses constant eta; native requires its lapse-replacement option to match the generated full-psi W initial lapse. Native KO uses `KO_DISS_SIGMA` only with CAKO off; enabling CAKO selects chi-scaled gauge/other coefficients instead.
- Role: public/scientific contract
- Deciding authority: `nrpy/examples/dendro_bssn.py`, `main`; `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`; `BSSN_GR/src/rhs.cpp`, CPU eta, SSL/CAHD include selection, and CAKO branches; `BSSN_GR/src/TwoPunctures.cpp`, lapse replacement.
- Corroboration: `nrpy/infrastructures/Dendro/CodeParameters.py`, q1 parameter mapping; `BSSN_GR/src/eta_RIT.inc.cpp`, radial formula; `BSSN_GR/src/parameters.cpp`, lapse and CAKO settings; `BSSN_GR/src/bssngr_main.cpp`, post-merger CAKO switch.

## Sources

- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - BSSN generation profile.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - BSSN RHS, SSL, CAHD, and Brown adjustment registration.
- [BSSN_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py) - BSSN constraint kernel registration.
- [ADM_to_BSSN.py](../../../nrpy/infrastructures/Dendro/general_relativity/ADM_to_BSSN.py) - ADM conversion.
- [initial_data_lambdaU.py](../../../nrpy/infrastructures/Dendro/general_relativity/initial_data_lambdaU.py) - connection initialization.
- [twopunctures.py](../../../nrpy/infrastructures/Dendro/general_relativity/twopunctures.py) - TwoPunctures data and interpolation.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - traversal and runtime scheduling.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - analytic puncture seed and generated startup sequence.
- [floor_the_lapse_and_conformal_factor.py](../../../nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py) - representation-dependent conformal-factor floor.
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - runtime mapping of eta and KO strengths to q1 parameter names.
- [param_toml.py](../../../nrpy/infrastructures/Dendro/param_toml.py) - emitted sample's eta and KO values.
- [Dendro-GR rhs.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/rhs.cpp) - native CPU eta and SSL/CAHD dispatch.
- [Dendro-GR TwoPunctures.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/TwoPunctures.cpp) - native initial-lapse replacement.
- [Dendro-GR parameters.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/parameters.cpp) - native CAKO and lapse settings.
- [Dendro-GR bssngr_main.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssngr_main.cpp) - post-merger CAKO switch.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- See also: [Constraints And Diagnostic Norms](constraints-and-diagnostic-norms.md)
- Contrasts with: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
