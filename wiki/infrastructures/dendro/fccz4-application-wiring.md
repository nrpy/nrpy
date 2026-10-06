# fCCZ4 Application Wiring

> Describe fCCZ4 equations, initial data, conversions, and runtime services in `Dendro_NRPy_fCCZ4`. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_fccz4` emits a complete fCCZ4 application with W or chi
selected at code-generation time (W by default). Its
24-component BSSN-compatible state is followed by `Theta_fCCZ4`. It uses SSL,
constant default `eta=1`, centered derivatives, its fCCZ4 CAHD contribution,
and separate conformal Ricci and RHS kernels.

## Detail

Each generated directory contains `pars/fccz4.toml` with the generated
runtime defaults and `pars/q1.par.lowres.toml` with the supplied equal-mass
TwoPunctures BBH parameters. The printed single-rank `--tpid` command and the
MPI evolution command both name the q1 file, and its values differ from the
generated defaults, among them `BSSN_CAHD_C` (0.06 in q1, 0.15 in the generated
fCCZ4 file); see [q1 versus generated
defaults](bssn-application-wiring.md#q1-versus-generated-defaults).
[Production Validation And Deferred Checks](production-validation-and-deferred-checks.md)
states which of these files and commands CI exercises. The packaged q1 file
starts a fresh run, sets a large end time of 1000000, omits an explicit
iteration cap, and uses native scaling with base output frequencies of 80. The
copied q1 file is identical in the BSSN and fCCZ4 directories; each executable
constructs its own evolved fields from the same physical initial data.

The independent `--ybs-gamma` and `--ybs-momentum` options, runtime
coefficients, paper references, and in-script KO setting are described in
[BSSN Application Wiring](bssn-application-wiring.md#optional-yo-et-al-adjustments).
The fCCZ4 shift driver receives the adjusted fCCZ4 connection RHS.

The fCCZ4 CAHD contribution is

```text
W_rhs   += 2 C_CAHD W H dt,    for W evolution;
chi_rhs += 4 C_CAHD chi H dt,  for chi evolution.
Default C_CAHD = 0.15.
```

`C_CAHD` is a runtime parameter, read from the key `BSSN_CAHD_C`. The chi term follows from `chi=W^2` and
`chi_rhs=2 W W_rhs`. Here `H` is the BSSN-shaped Hamiltonian expression used
by the RHS adjustment and now also reported separately; it is not the distinct
`H_Z4` constraint. This equation is not the BSSN CAHD equation. fCCZ4 retains
its formulation's Lambda RHS rather than adding Brown's BSSN adjustment. Order-specific
constraint kernels are named `fCCZ4_constraints_order_N`.

The slow-start lapse is the relaxation described in [BSSN Application
Wiring](bssn-application-wiring.md), with the same keys `BSSN_SSL_H` and
`BSSN_SSL_SIGMA`. The Z4 damping coefficients `kappa1` (default 0.1) and `kappa2`
(default 0) are read from keys of the same names; the packaged q1 file does not
set them.

Claim evidence:
- Claim: fCCZ4's W and chi CAHD terms use the BSSN-shaped Hamiltonian expression with a runtime coefficient read from the key `BSSN_CAHD_C` (generated default 0.15; the packaged q1 file sets 0.06), not `H_Z4`; changing the generated conformal factor changes the coefficient from 2 to 4. The Z4 damping coefficients `kappa1` and `kappa2` are read from keys of the same names with defaults 0.1 and 0.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`, `register_CFunction_rhs_eval`; `nrpy/equations/general_relativity/fCCZ4_RHSs.py` (`kappa1`, `kappa2` registration); `nrpy/infrastructures/Dendro/CodeParameters.py`, `Q1_TOML_PARAMETER_NAMES`.
- Corroboration: `nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py`, `register_CFunction_fCCZ4_constraints`, emits both `H` and formulation-specific `H_Z4` diagnostics.

The initial octree uses the same analytic Dendro-GR puncture seed as BSSN.
Fresh evolution loads the same precomputed TwoPunctures coefficients as BSSN.
Evolved initial data uses the resulting ADM data and full
`psi=psi_background+u`, `alpha=W=psi^(-2)`, ADM conversion, halo exchange,
physical-boundary fill, separate `initial_data_lambdaU`, and algebraic
projection sequence as BSSN. `Theta_fCCZ4` receives its formulation-defined
initial value. W and chi builds have distinct checkpoint formulation IDs and
reject cross-formulation restores.

Claim evidence:
- Claim: Fresh fCCZ4 evolution can load the same precomputed TwoPunctures solution as generated BSSN, then converts ADM fields and initializes the fCCZ4 state; checkpoint formulation IDs remain distinct.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`; `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::initialize` within `output_solver_context_cpp`.
- Corroboration: `nrpy/examples/dendro_fccz4.py` and `nrpy/examples/dendro_bssn.py`, use the same entry-point generator and TwoPunctures registration.

After initialization, every RK stage, and AMR transfer, `alpha` is floored at
`CHI_FLOOR` and the evolved conformal factor at `CHI_FLOOR` for chi or
`sqrt(CHI_FLOOR)` for W before algebraic projection.

Solver context takes each block through the same sequence as BSSN: exterior
padding fill, Ricci, RHS, then the outgoing-radiation condition on the physical
faces (see [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)).
Six `RbarDD` values remain scratch storage,
not evolved state. Apparent-horizon searches interpolate the evolved
fCCZ4 fields and convert each search point to ADM variables. ADM surface
quantities use grid fields from `BSSN_to_ADM`; waveform extraction evaluates
Psi4 directly from the evolved fields and writes per-mode files plus the
per-radius real and imaginary L2 file described in [BSSN Application
Wiring](bssn-application-wiring.md). The shared runtime also writes tracked
puncture positions and a per-launch TOML parameter dump, described there.

The diagnostic kernel retains `H_Z4` and the three `Z4constraintU` components
as its first four fields. It then reports `H`, physical lower-index momentum
components `MD0..MD2`, `M_CONSTRAINT=sqrt(gamma_ij M^i M^j)`, and
`LAMBDA_CONSTRAINT=sqrt(gammabar_ij Z4constraintU^i Z4constraintU^j)`.
The connection residual is reported, not enforced. Both excised physical-volume
and unique-node RMS files include all ten fields and are emitted after remeshing.
VTU output uses the shared selection described in [BSSN Application
Wiring](bssn-application-wiring.md): constraint index 0 writes the BSSN-form
`H`, and `H_Z4` and `Z4constraintU` are not selectable for VTU output.

The shared runtime parses puncture-centered AMR controls, tracks centers with
`vetU`, retains puncture history in checkpoints, and computes diagnostics
after any scheduled remesh. This shared path does not make the fCCZ4
constraint or CAHD equations identical to BSSN's.

## Sources

- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - fCCZ4 generation profile.
- [param_toml.py](../../../nrpy/infrastructures/Dendro/param_toml.py) - sample parameter file and TwoPunctures lapse default.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - fCCZ4 RHS and CAHD registration.
- [fCCZ4_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py) - fCCZ4 constraint registration.
- [Ricci_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - separate conformal Ricci kernel.
- [ADM_to_BSSN.py](../../../nrpy/infrastructures/Dendro/general_relativity/ADM_to_BSSN.py) - ADM conversion.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - traversal and runtime scheduling.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - analytic initial-grid seed and startup.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [fCCZ4 System](../../equations/general-relativity/fccz4.md)
- Contrasts with: [BSSN Application Wiring](bssn-application-wiring.md)
- See also: [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)
- See also: [Runtime Parameter Keys](runtime-parameters.md)
