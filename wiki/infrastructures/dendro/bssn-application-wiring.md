# BSSN Application Wiring

> Describe BSSN equations, initial data, conversions, and runtime services in `Dendro_NRPy_BSSN`. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_bssn` emits a complete BSSN application with a
code-generation choice of W (the default) or chi. It uses centered
derivatives, SSL, constant default `eta=1`, Dendro-style CAHD, Brown's
covariant Lambda adjustment, and separate conformal Ricci and RHS kernels.
Runtime parameter input can override the generated eta default, but cannot
change the evolved conformal factor of an already generated binary.

## Detail

Each generated directory contains `pars/bssn.toml` with the generated
runtime defaults and `pars/q1.par.lowres.toml` with the supplied equal-mass
TwoPunctures BBH parameters. Both files work with the single-rank `--tpid`
command and the MPI evolution command. The packaged q1 file starts a fresh
run, sets a large end time of 1000000, omits an explicit iteration cap,
and uses native scaling with base output frequencies of 80. The copied q1 file
is identical in the BSSN and fCCZ4 directories; each executable constructs its
own evolved fields from the same physical initial data.

The BSSN CAHD contribution depends on the evolved conformal factor:

```text
W_rhs += (C_CAHD/2) W H dx^2 / (dt (1 + 10 dx^2)),  default C_CAHD = 0.06.
chi_rhs += C_CAHD chi H dx^2 / (dt (1 + 10 dx^2)).
```

`C_CAHD` is a runtime parameter; `dx` is the smallest of the current block's
three physical axis spacings, and `dt` is the current timestep. This choice
preserves the cubic-grid coefficient and limits directional damping on
rectangular domains. The factor of two
follows `chi=W^2`; it preserves the same CAHD perturbation to the physical
conformal factor. SSL uses `W=cf` or `W=sqrt(cf)` in its relaxation toward
`alpha=W`. KO has no W-dependent scaling in this profile (`enable_CAKO=False`).
Brown's adjustment appears
exactly once in the covariant `lambdaU` RHS. The generated constraint kernels
are named `BSSN_constraints_order_N`; there is no generic constraint catch-all
source.

Claim evidence:
- Claim: The W and chi BSSN variants receive the corresponding CAHD factor using the minimum physical block spacing; SSL uses W in either variant, while KO is not scaled by W. The conformal representation is a code-generation choice, not a runtime switch.
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

`twopunctures` and `ADM_to_BSSN` act on block interiors only: `zip` reads
only interior values, and the following unzip and physical-boundary fill
supply the padding that `initial_data_lambdaU` differentiates.

Claim evidence:
- Claim: Fresh evolution loads a precomputed TwoPunctures solution, converts its ADM data, computes the initial connection field from exchanged block data, then zips and floors/projects owned evolved nodes.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`; `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::initialize` within `output_solver_context_cpp`.
- Corroboration: `nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py` and `enforce_detgbar_equals_detghat_trAzero.py`, owned-node CFunction signatures.

Apparent-horizon searches interpolate selected evolved BSSN fields and then
convert each search point to ADM variables. `BSSN_to_ADM` reconstructs grid
fields used by ADM surface quantities; waveform extraction evaluates Psi4
directly from the selected BSSN fields. Generated runtime services also provide
physical boundaries, conformal-factor volume-weighted and unique-node constraint
diagnostics, checkpoint and restart, and scheduled output. ADM surface
quantities are written to `*_ADM.dat`: at the effective gravitational-wave output cadence,
when `BSSN_GW_RADAII` is non-empty, one row with the step, time, outermost
extraction radius, and the ADM energy, linear momentum, and angular momentum,
printed with up to ten significant digits (default notation, trailing zeros
dropped). The labeled ASCII diagnostic files (`*_ADM.dat`, the
constraint files, the Psi4 mode files, `*_GW_L2.dat`, and
`*_BHLocations.dat`) open with column labels in the style of BHaHAHA's horizon
diagnostics files: when the file is missing or empty, rank 0 first writes a
title line naming the formulation and evolved conformal factor, then one
`# column N = <name>: <meaning>` line per column. Dendro-GR writes
uncommented header lines for its corresponding waveform-norm and puncture
location files. A restart appends below the existing local labels.

Claim evidence:
- Claim: At the effective gravitational-wave output cadence, when `BSSN_GW_RADAII` is non-empty, rank 0 appends one row to `<BSSN_PROFILE_FILE_PREFIX>_ADM.dat` with the step, time, the last listed extraction radius, and the seven ADM surface quantities (energy, three linear-momentum and three angular-momentum components), printed with up to ten significant digits in default notation. This applies to the fCCZ4 application as well.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::adm_output` within `output_solver_context_cpp`; `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`, the `BSSN_GW_EXTRACT_FREQ` and `BSSN_GW_RADAII` reads; `Ctx::update_output_frequencies`.
- Corroboration: `nrpy/infrastructures/Dendro/general_relativity/adm_quantities.py`, `register_CFunction_adm_quantities`, the order of the seven quantities; `nrpy/examples/tests/dendro_application_check.py`, `Leg.check_run_a` column use.

Claim evidence:
- Claim: When `*_ADM.dat`, either constraint file, a `*_GW_l<l>_m<m>.dat` file, `*_GW_L2.dat`, or `*_BHLocations.dat` is missing or empty, rank 0 writes `# <title>`, `#`, and one `# column N = <name>: <meaning>` line per column before the first row; a file that already has content receives no further labels. The step and time columns are named `TimeStep` and `time` (`t` in the Psi4 files); Psi4 radius columns are `r0`, `r1`, ..., as in native `BSSN_GR`'s headers; the constraint and puncture-location columns use descriptive generated labels. Native waveform-norm and puncture-location files instead use uncommented headers. This applies to the fCCZ4 application as well.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `open_labeled_output` and its output paths within `output_solver_context_cpp`, and `diagnostic_meanings`.
- Corroboration: `nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_file_output.py`, the horizon diagnostics header this format follows; `nrpy/examples/tests/dendro_application_check.py`, `parse_table` and `parse_modes`, which read the labels.

Each Psi4 mode goes to its own `*_GW_l<l>_m<m>.dat` file, with Dendro-GR
`BSSN_GR`'s name and row layout: the step, time and one complex `(Re,Im)` pair
per extraction radius. Instead of native's uncommented step-0 header line, the
column labels keep its names (`TimeStep`, `t`, `r0`, `r1`, ...) and give each
radius's value. A reader that takes column names from native's first line, such
as `BSSN_GR/scripts/getstrain.py` with `pandas.read_csv(sep='\t')`, must instead
skip `#` lines and name the columns from the labels. Files are written for every
`l` from 2 to the largest entry of `BSSN_GW_L_MODES`; native `BSSN_GR` writes
only the listed `l` values, which is the same set for the default list.

Claim evidence:
- Claim: Rank 0 appends each (l, m) mode, for l = 2..max(`BSSN_GW_L_MODES`) and m = -l..l, to `<BSSN_PROFILE_FILE_PREFIX>_GW_l<l>_m<m>.dat`: column labels, then one row per extraction step with the step, time and one `(Re,Im)` pair per extraction radius, in scientific notation with 10 digits after the decimal point. Native `BSSN_GR` uses the same name and row layout, with an uncommented step-0 header line instead of the labels, and writes only the listed l values. This applies to the fCCZ4 application as well.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::gravitational_wave_output` within `output_solver_context_cpp`.
- Corroboration: Dendro-GR `BSSN_GR/include/gwExtract.h`, `GW::extractFarFieldPsi4`, per-mode file writer.

At the effective gravitational-wave extraction cadence, rank 0 also appends
`<BSSN_PROFILE_FILE_PREFIX>_GW_L2.dat`. Each row gives the step and time, then
one `(Re,Im)` pair per extraction radius. The two values are
`sqrt(sum(Re(Psi4)^2))` and `sqrt(sum(Im(Psi4)^2))`, with each sum over the
valid Lebedev samples and reduced across MPI ranks. The sums use no Lebedev
weights or `4*pi` factor. Labels name the columns `TimeStep`, `t`, `r0`, `r1`,
and so on, and give each radius. Unlike Dendro-GR's uncommented step-zero
header, the generated file uses comment labels. A fresh run writes step 0. A
restart skips the checkpoint step already written and appends the next
scheduled row.

The same cadence writes
`<BSSN_PROFILE_FILE_PREFIX>_BHLocations.dat` with the step, time, and Cartesian
coordinates of the two tracked puncture centers. These are puncture/excision
center positions, not apparent-horizon centers. The labeled header follows the
other local diagnostic files, while Dendro-GR writes an uncommented header; a
restart appends rows below the existing labels.
The initial fresh-run row uses the initial centers, and evolved rows use the
centers updated before diagnostics.

Claim evidence:
- Claim: At each active gravitational-wave extraction step, rank 0 appends `<BSSN_PROFILE_FILE_PREFIX>_GW_L2.dat` with the step, time and one complex pair per radius. Each pair contains the square roots of the MPI-reduced, unweighted sums of the squared real and imaginary interpolated Psi4 samples over valid Lebedev points. Rank 0 also appends `<BSSN_PROFILE_FILE_PREFIX>_BHLocations.dat` with the step, time and coordinates of both tracked puncture centers. Fresh runs write step 0; restarts do not repeat the checkpoint step. The same behavior is used by the fCCZ4 application.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/gravitational_waves.py`, `register_CFunction_gravitational_waves`; `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::gravitational_wave_output` and `Ctx::black_hole_locations_output`; `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`.
- Corroboration: Dendro-GR `BSSN_GR/include/gwExtract.h`, `GW::extractFarFieldPsi4`; Dendro-GR `BSSN_GR/src/dataUtils.cpp`, `writeBHCoordinates`; `nrpy/examples/tests/dendro_application_check.py`, universal output parsing.

On normal solver startup, rank 0 writes
`<BSSN_PROFILE_FILE_PREFIX>__PARAM_DUMP__YYYY-MM-DD-HH-MM-SS.toml`, using
machine local time in the filename. The TOML retains supplied keys and tables,
and adds consumed host fallback and registered runtime CodeParameter defaults
only when every host fallback and mapped CodeParameter default for a missing
key agrees. A key with differing call-site fallbacks or a host fallback that
differs from its CodeParameter default remains absent so rerunning the file
preserves each call site's current behavior. In particular, omitted
`BSSN_BH1.V_X` and `BSSN_BH1.V_Y` retain their distinct defaults for puncture
tracking and TwoPunctures momentum; omitted `CHI_FLOOR` retains its
`0.1` initial-mesh puncture-seed floor and its `1e-4` evolved-field floor. The
dump writes every inserted floating-point default with enough digits to
recover its original double value, including values inside TOML arrays. A
rerun can therefore reuse the same TwoPunctures solution file when its input
parameters agree. The special `--tpid` utility run does not write an evolution
parameter dump.

Claim evidence:
- Claim: A normal solver launch writes one rank-0 TOML file named `<BSSN_PROFILE_FILE_PREFIX>__PARAM_DUMP__YYYY-MM-DD-HH-MM-SS.toml`. It preserves supplied keys and tables, and adds consumed host fallback and registered runtime CodeParameter defaults only when every host fallback and mapped CodeParameter default for a missing key agrees. Inserted floating-point defaults retain enough digits to recover their double values, so the same TwoPunctures solution file can be reused. Conflicting keys remain absent to preserve their distinct call-site behavior. `--tpid` does not write the file. This applies to the fCCZ4 application as well.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `ParameterFile` and `output_main_cpp`; `nrpy/infrastructures/Dendro/CodeParameters.py`, `output_toml_default_assignments`.
- Corroboration: Dendro-GR `BSSN_GR/src/bssngr_main.cpp`, `writeParamTOMLFile` call; Dendro-GR `BSSN_GR/src/parameters.cpp`, resolved parameter writer.

VTU output follows Dendro-GR `BSSN_GR`'s selection and defaults. With
`BSSN_VTU_Z_SLICE_ONLY` (default true) it writes the elements that touch the
z-normal plane through the domain center (z = 0 for a domain centered on the
origin); otherwise it writes the full volume. The first
`BSSN_NUM_EVOL_VARS_VTU_OUTPUT` entries of `BSSN_VTU_OUTPUT_EVOL_INDICES` select
evolved fields by generated index. The generated BSSN order corresponds to
`BSSN_GR`'s variable order position by position, but the fields are NRPy's
forms: the selected conformal factor in place of chi, `lambdaU`, `hDD` (the
deviation of the conformal metric from flat) and `aDD` in place of `Gt`, `gt`
and `At`, and B^i in NRPy's normalization. fCCZ4 adds `Theta_fCCZ4` at index
24. The first `BSSN_NUM_CONST_VARS_VTU_OUTPUT` entries of
`BSSN_VTU_OUTPUT_CONST_INDICES` select, in `BSSN_GR` numbering, `H`,
`MU0..MU2`, `psi4_real` and `psi4_imag`. Both counts default to 1, so a
parameter file without these keys writes `alpha` and `H`. Requested constraint
or Psi4 fields are computed at the output step. Field names are NRPy's (`H`,
not `C_HAM`).

Claim evidence:
- Claim: With `BSSN_VTU_Z_SLICE_ONLY` true (the default), VTU output writes the z-normal slice through the domain center, and otherwise the full volume. It writes the selected evolved fields followed by the selected constraint and Psi4 fields, and computes constraints or Psi4 only when they are selected. Evolved indices refer to the generated field order; constraint indices use `BSSN_GR`'s order (C_HAM, C_MOM0-2, C_PSI4_REAL, C_PSI4_IMG). Both counts default to 1. This applies to the fCCZ4 application as well.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::write_vtu` within `output_solver_context_cpp`; `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`.
- Corroboration: Dendro-GR `BSSN_GR/src/bssnCtx.cpp`, `BSSNCtx::write_vtu`; `BSSN_GR/src/parameters.cpp`, VTU parameter defaults.

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

| Choice | Generated `Dendro_NRPy_BSSN` | Native CPU `BSSN_GR` |
| --- | --- | --- |
| Evolved conformal factor | W by default; `--conformal-factor chi` selects chi at generation | chi |
| TwoPunctures initial lapse | Always `alpha=(psi_background+u)^(-2)=W`; startup rejects `TPID_REPLACE_LAPSE_WITH_SQRT_CHI = false`, and `INITIAL_LAPSE` and `TPID_INITIAL_LAPSE_PSI_EXPONENT` are reported as having no effect | `TPID_REPLACE_LAPSE_WITH_SQRT_CHI=true` is needed to replace the ordinary `INITIAL_LAPSE=2` result with full-psi W |
| Shift-driver damping in the SSL/CAHD path | Spatially constant `eta`; default 1, overridden by `ETA_CONST` in an input file | Radial RIT profile in the CPU RHS; `ETA_CONST` does not select constant damping on this path |
| Gamma-driver auxiliary `B^i` | `d_t beta^i = B^i + advection`, `d_t B^i = (3/4) d_t Lambdabar^i - eta B^i + advection`; the wavelet refinement test scales `betU` by 4/3 (see [Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)) | `d_t beta^i = (3/4) B^i + advection` with `BSSN_LAMBDA_F = (1, 0)`, `d_t B^i = d_t Gt^i - eta B^i + advection`; with `BSSN_LAMBDA = (1, 1, 1, 1)` and the same `eta`, `B` is 4/3 of the generated `betU`; where the RIT profile differs from the generated constant `eta`, the relation does not hold |
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

### Optional Yo et al. adjustments

Both examples accept `--ybs-gamma` and `--ybs-momentum` independently or
together. Both default to off. `--ybs-gamma` forwards the canonical Gamma
constraint adjustment to the evolution and shift-driver equations. Its runtime
coefficient is `YBS_chi`, defaulting to 0, which disables the term; the
recommended maximum is 2/3 for BSSN and 4/3 for fCCZ4, a customary choice for
BSSN and a matching choice for fCCZ4, not a stability bound (see
[BSSN Family](../../equations/general-relativity/bssn-family.md) and
[fCCZ4](../../equations/general-relativity/fccz4.md)). This is the additional
NRPy coefficient, with Brown's BSSN contribution retained separately. The help cites
[Yo, Baumgarte, and Shapiro, arXiv:gr-qc/0209066](https://arxiv.org/abs/gr-qc/0209066),
Eq. (45), and [Yo, Lin, and Cao, arXiv:1205.5111](https://arxiv.org/abs/1205.5111),
Eq. (47).

`--ybs-momentum` adds the covariant, symmetric trace-free gradient of the
lower conformal momentum residual to the conformal extrinsic-curvature RHS.
The help cites Yo, Lin, and Cao, Eq. (56). NRPy multiplies this term by
`C_YBS_mom * BSSN_CFL_FACTOR * min(abs(dx), abs(dy), abs(dz))`, using the
current Cartesian block spacing and runtime CFL factor; `C_YBS_mom` defaults
to 0, which disables the term, and its recommended maximum with RK4 is
`0.18 / BSSN_CFL_FACTOR^2` (see [YBS-MOM](../../equations/general-relativity/ybs-momentum-damping.md#strength-and-timestep-bound)).
This local coefficient is NRPy's timestep scaling of the paper's term.
It adds no evolved gridfunction. See [YBS-MOM](../../equations/general-relativity/ybs-momentum-damping.md)
for the equation and the limits of a damping or stability claim.

Kreiss--Oliger generation is controlled by
`enable_KreissOliger_dissipation` inside each example, defaulting to `True`.
There are no `--ko` or `--no-ko` arguments. The runtime strength remains
`KO_DISS_SIGMA`.

Claim evidence:
- Claim: The Dendro examples expose independent default-off Gamma and momentum adjustments whose runtime strengths `YBS_chi` and `C_YBS_mom` default to 0, which removes each term; the recommended maxima are advice from the linked equation pages, not a stability guarantee; the momentum coefficient uses the current block's minimum physical spacing and the evolution CFL factor, without adding evolved fields.
- Role: descriptive behavior
- Deciding authority: `nrpy/examples/dendro_bssn.py` and `nrpy/examples/dendro_fccz4.py`, `parse_args` and `main`; `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`, `register_CFunction_rhs_eval`; `nrpy/infrastructures/Dendro/CodeParameters.py`, `Q1_TOML_PARAMETER_NAMES`.
- Corroboration: `nrpy/equations/general_relativity/BSSN_RHSs.py` and `fCCZ4_RHSs.py`, canonical adjustment construction; the cited Yo et al. equations.

## Sources

- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - BSSN generation profile.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - BSSN RHS, SSL, CAHD, and Brown adjustment registration.
- [BSSN_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py) - BSSN constraint kernel registration.
- [ADM_to_BSSN.py](../../../nrpy/infrastructures/Dendro/general_relativity/ADM_to_BSSN.py) - ADM conversion.
- [initial_data_lambdaU.py](../../../nrpy/infrastructures/Dendro/general_relativity/initial_data_lambdaU.py) - connection initialization.
- [twopunctures.py](../../../nrpy/infrastructures/Dendro/general_relativity/twopunctures.py) - TwoPunctures data and interpolation.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - traversal and runtime scheduling.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - puncture seed, effective parameter values, and generated startup sequence.
- [floor_the_lapse_and_conformal_factor.py](../../../nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py) - representation-dependent conformal-factor floor.
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - runtime mapping of eta and KO strengths to q1 parameter names and recorded compatible TOML defaults.
- [gravitational_waves.py](../../../nrpy/infrastructures/Dendro/general_relativity/gravitational_waves.py) - Psi4 interpolation, mode decomposition, and real/imaginary L2 sums.
- [param_toml.py](../../../nrpy/infrastructures/Dendro/param_toml.py) - emitted sample's eta and KO values.
- [Dendro-GR rhs.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/rhs.cpp) - native CPU eta and SSL/CAHD dispatch.
- [Dendro-GR TwoPunctures.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/TwoPunctures.cpp) - native initial-lapse replacement.
- [Dendro-GR parameters.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/parameters.cpp) - native CAKO and lapse settings.
- [Dendro-GR bssngr_main.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssngr_main.cpp) - parameter dump name, startup, and diagnostic cadence.
- [Dendro-GR dataUtils.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/dataUtils.cpp) - native puncture-coordinate output.
- [Dendro-GR bssnCtx.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssnCtx.cpp) - native VTU field selection and slicing.
- [Dendro-GR gwExtract.h](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/include/gwExtract.h) - native Psi4 L2 and per-mode file names and layouts.
- [diagnostics_file_output.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_file_output.py) - BHaHAHA horizon diagnostics header, the column-label format.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- See also: [Constraints And Diagnostic Norms](constraints-and-diagnostic-norms.md)
- Contrasts with: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
