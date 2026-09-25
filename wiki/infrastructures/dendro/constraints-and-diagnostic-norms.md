# Constraints And Diagnostic Norms

> Explain Dendro BSSN and fCCZ4 momentum and connection diagnostics, excision, and two reported constraint norms. · Status: confirmed
> Up: [Dendro](index.md)

## Summary

The Dendro BSSN and fCCZ4 connectors lower NRPy's contravariant momentum constraint
before writing the fixed `MU0..MU2` diagnostic slots. It also reports scalar
momentum and conformal connection residual magnitudes. Separate files contain
conformal-factor volume-weighted and unique-node RMS values; the latter matches the native
Dendro-BSSN node-weighting convention. Diagnostic kernels run on the
post-remesh state, including at the initial grid.

## Detail

The equation module constructs `M^i`. The Dendro connector writes
`M_i = gamma_ij M^j = gammabar_ij M^j / exp(-4 phi)` to the same ABI slots.
Here `exp(-4 phi)=W^2=chi`. This is a tensorial index lowering, not a single
scalar rescaling: off-diagonal conformal-metric terms contribute. The
Hamiltonian expression is unchanged. The conversion occurs immediately before
C code generation, leaving the shared equation module untouched.

Claim evidence:
- Claim: Dendro's BSSN `MU` diagnostic slots contain `M_i = gammabar_ij M^j / exp(-4 phi)`, while the equation module still constructs `M^i`.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py`, `register_CFunction_BSSN_constraints`, `momentum_covariant`.
- Corroboration: `nrpy/equations/general_relativity/BSSN_constraints.py`, `BSSN_constraints.MU`; `nrpy/equations/general_relativity/BSSN_quantities.py`, `BSSN_quantities.exp_m4phi`.

After `H` and `MU0..MU2`, the BSSN diagnostic kernel writes
`M_CONSTRAINT = sqrt(gamma_ij M^i M^j)` and
`LAMBDA_CONSTRAINT = sqrt(gammabar_ij C^i C^j)`, where
`C^i = Lambdabar^i - DeltaGamma^i`. Both contractions include off-diagonal
metric components. The generated state header preserves this field order; the
gridfunction registry's alphabetical order does not determine diagnostic array
indices. fCCZ4 retains `H_Z4` and three `Z4constraintU` components before
these six BSSN-comparable fields; its `LAMBDA_CONSTRAINT` contracts
`Z4constraintU` with `gammabar_ij`. Reporting either residual does not enforce
it or change its evolution equation.

Claim evidence:
- Claim: Dendro BSSN reports the physical-metric scalar momentum magnitude and conformal-metric scalar connection-residual magnitude after the three lower-index momentum components; the generated header and diagnostic kernel use the same indices, and neither residual is enforced by this output.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py`, `register_CFunction_BSSN_constraints`; `nrpy/infrastructures/Dendro/state_h.py`, `BSSN_DIAGNOSTIC_GRIDFUNCTIONS` and `output_state_h`.
- Corroboration: `nrpy/equations/general_relativity/BSSN_constraints.py`, `BSSNconstraints.Msquared`, `LambdaConstraintSquared`, and `LambdaConstraintMagnitude`; `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::diagnostic_output`.

Claim evidence:
- Claim: fCCZ4 reports its four original Z4 diagnostic fields followed by six BSSN-comparable fields, including lower-index momentum and covariant scalar magnitudes; `H` is the BSSN-shaped Hamiltonian, distinct from `H_Z4`.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py`, `register_CFunction_fCCZ4_constraints`; `nrpy/infrastructures/Dendro/state_h.py`, `FCCZ4_DIAGNOSTIC_GRIDFUNCTIONS`.
- Corroboration: `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`, `register_CFunction_rhs_eval`, uses the same BSSN-shaped `H` in CAHD.

`*_Constraints_volweighted.dat` contains a conformal-factor volume-weighted RMS and
pointwise maximum for each diagnostic field, including both new scalars. The
quadrature excludes spherical puncture regions and weights each retained
quadrature point by a tensor-product trapezoidal factor times
`W^-3 dx dy dz` for W evolution or `chi^-3/2 dx dy dz` for chi evolution.
This equals `sqrt(gamma) dx dy dz` when `det(gammabar)=1`, as imposed by the
algebraic projection; interpolation can leave small determinant errors at
diagnostic points. The earlier W-only `W^-6` weight is not the physical
volume element even under that condition. This norm answers a geometric-volume
question; it is not interchangeable with Dendro's node norm.

`*_Constraints.dat`, named as in Dendro-GR `BSSN_GR`, contains one RMS per
diagnostic field over locally owned, unique continuous-Galerkin nodes outside
the same spherical excision regions. It counts a shared node once and weights
each retained node equally, matching the native Dendro-BSSN reporting
convention. The file also records the node count after all diagnostic RMS
columns (six for BSSN and ten for fCCZ4). For BSSN, columns 3-6 (`H`,
`MU0..MU2`) correspond to native `C_HAM` and `C_MOM0..2`; fCCZ4 lists its four
Z4 fields first. Rows in both files begin with the step and physical time. Each
file opens with column labels in the style of BHaHAHA's horizon diagnostics
files: a title line naming the formulation and evolved conformal factor, then
one `# column N = <name>: <meaning>` line per column. `*_Constraints.dat` names
its columns `TimeStep` and `time`, as native `BSSN_GR`'s header does, then the
generated diagnostic names and `unexcised_nodes`; the volume-weighted file names
each field's pair `<name>_rms` and `<name>_max`. The volume-weighted title
states its weight, `W^-3 dx dy dz` or `chi^-3/2 dx dy dz`. Compare against
native results by physical time, not merely step number. Both files print
floating-point values with up to ten significant digits (default notation,
trailing zeros dropped).

Claim evidence:
- Claim: The generated solver reports distinct excised conformal-factor volume-weighted and unique-owned-node RMS norms, and only the latter matches native Dendro-BSSN's node-weighting convention. The unique-node norm goes to `*_Constraints.dat`, native `BSSN_GR`'s file name, labelled `TimeStep` and `time` (native `BSSN_GR`'s names), the generated diagnostic names, and `unexcised_nodes`; for BSSN its columns 3-6 correspond to native `C_HAM` and `C_MOM0..2`. The volume-weighted norm goes to `*_Constraints_volweighted.dat`, which labels an RMS and a maximum absolute value per field. A new or empty file of either kind first receives a title line and one `# column N = <name>: <meaning>` line per column. Both files print floating-point values with up to ten significant digits.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `output_solver_context_cpp` / `Ctx::diagnostic_output`, `open_labeled_output`, and `diagnostic_meanings`; `nrpy/infrastructures/Dendro/general_relativity/diagnostics.py`, `register_CFunction_diagnostics`.
- Corroboration: `BSSN_GR/include/grUtils.tcc`, `bssn::computeConstraintL2Norm(const ot::Mesh*,...)`, native ownership, excision, and RMS reduction, and `bssn::extractConstraints`, native file name and its `TimeStep` and `time` header names; `nrpy/infrastructures/Dendro/state_h.py`, diagnostic component order.

The generated main loop evolves, remeshes and transfers when scheduled, then
advances puncture centers before diagnostic output. This ordering prevents a
pre-remesh diagnostic from being compared with native post-remesh output at
the same reported time. Matching initial data, field representation, mesh,
excision, and norm remain separate prerequisites; agreement of one scalar norm
does not prove equality of the evolved fields or equations.

## Sources

- [BSSN_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py) - Dendro momentum-index lowering at code generation.
- [fCCZ4_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py) - fCCZ4 Z4 and BSSN-comparable diagnostic construction.
- [Equation BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py) - `BSSNconstraints.MU`, the upper-index momentum residual.
- [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py) - conformal metric and `exp_m4phi`.
- [state_h.py](../../../nrpy/infrastructures/Dendro/state_h.py) - canonical diagnostic order and generated header names.
- [diagnostics.py](../../../nrpy/infrastructures/Dendro/general_relativity/diagnostics.py) - conformal-factor volume-weighted accumulation and excision.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - unique-node norm and output.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - post-remesh output order.
- [Dendro-GR grUtils.tcc](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/include/grUtils.tcc) - native unique-owned-node RMS and excision convention, and the constraint file name and header.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [BSSN Application Wiring](bssn-application-wiring.md)
- See also: [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)
