# Production Validation And Deferred Checks

> Define checks for complete generated Dendro applications and state current limits. · Status: provisional
> Up: [Dendro](index.md)

## Summary

Production examples generate complete Dendro applications.

## Detail

Required generator checks, with the layer that enforces each:

- Isolated static analysis of each generator module: the `static-analysis` CI
  job, whose file search skips `*/tests/*` paths among others, so the helper
  itself is not analyzed.
- BSSN and fCCZ4 import and registration, and two byte-identical clean
  generations: the helper generates each W and chi variant twice and compares
  the trees (a CI assertion).
- Presence of all order-specific kernels: a generation-time `ValueError` from
  `output_solver_context_cpp`. The RHS registrar raises at generation when its
  keys differ from the canonical order, and `validate_registered_state` raises
  when the registered names differ from the canonical state.
- Direct Python/CFunction/C++ name correspondence and final-state pointer
  indices: the kernel-presence guard above enforces that the order-specific
  kernels and the five services `BSSN_to_ADM`, the lapse and conformal-factor
  floor, `gravitational_waves`, `adm_quantities`, and `physical_boundary_ghosts`
  are registered; the calls to `twopunctures`,
  `physical_boundary`, `diagnostics`, `apparent_horizon`, and
  `enforce_detgbar_equals_detghat_trAzero` have no generation-time registration
  check. That a Python module's basename equals its registrar suffix, and that
  pointer indices follow the canonical state order beyond the RHS key-order
  guard, are review-time checks.
- Rejection of native `BSSN_GR` sources: a generation-time `ValueError` from
  `output_CFunctions_function_prototypes_and_construct_CMakeLists`. The screen
  applies only to the application-supplied sources; registered CFunction sources
  get duplicate and path checks, and no test triggers the screen.

Claim evidence:
- Claim: The generators raise `ValueError` at generation when an order-specific kernel or one of the named runtime services (`BSSN_to_ADM`, the lapse and conformal-factor floor, `gravitational_waves`, `adm_quantities`, `physical_boundary_ghosts`) is not registered (`output_solver_context_cpp`; the calls to `twopunctures`, `physical_boundary`, `diagnostics`, `apparent_horizon`, and `enforce_detgbar_equals_detghat_trAzero` have no registration check), when the RHS keys differ from the canonical evolved order (`register_CFunction_rhs_eval`), when the registered EVOL or AUXEVOL names differ from the canonical state (`validate_registered_state`), and when an application-supplied CMake source is unsafe, duplicated, or a native `BSSN_GR` source (`output_CFunctions_function_prototypes_and_construct_CMakeLists`); registered CFunction sources get only the duplicate and path checks; the generators contain no check that a Python module's basename equals its registrar suffix and none of pointer indices beyond the canonical-order check of the RHS keys.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/solver_context.py`, `output_solver_context_cpp`; `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`, `register_CFunction_rhs_eval`; `nrpy/infrastructures/Dendro/state_h.py`, `validate_registered_state`; `nrpy/infrastructures/Dendro/CMakeLists.py`, `output_CFunctions_function_prototypes_and_construct_CMakeLists`
- Corroboration: none available; the helper generates every variant, so each guard runs in CI, but no case makes a guard fail

Required application checks configure and build both standalone applications with
Dendrolib. Both BSSN and fCCZ4 checks must generate W and chi variants.
MPI TwoPunctures runs must exercise finite-difference orders 4, 6,
and 8 through initialization and the first completed diagnostic output. Checks
must cover nontrivial finite fields, `alpha=W=sqrt(chi)` for both formulations,
Lambda initialization, algebraic projection, block-boundary continuity,
constraints, wave extraction, apparent horizons, remeshing, checkpoint,
and restart. W/chi checkpoint-formulation mismatches must be rejected.
Temporary direct numerical comparisons must cover interior stencils, block
boundaries, conformal-factor algebra, SSL, CAHD, and separate Ricci/RHS
coupling. A native-control comparison is meaningful only with the same initial data,
post-remesh state, puncture excision, momentum-component convention, and
unique-node RMS; conformal-factor volume-weighted RMS is a distinct diagnostic. Compare by physical time and check
the initial diagnostic before interpreting long evolutions. Record the evolved
conformal-factor choice, native full-psi lapse replacement, eta prescription,
KO strength, native CAKO state, any post-merger CAKO switch, puncture tracker,
time integrator and CFL spacing, and the defaults of omitted keys before
attributing any later difference to the formulation.

Claim evidence:
- Claim: A generated Dendro-BSSN run and a native Dendro-GR run that share one parameter file still differ in the evolved conformal factor, lapse replacement, eta prescription, KO strength and CAKO state (including any post-merger CAKO switch), puncture tracker, time integrator and CFL spacing, omitted-key defaults, constraint-norm weighting, and momentum-component convention, so a comparison that matches only the parameter file does not isolate the formulation. W and chi checkpoints carry distinct formulation identifiers, and a restore across them is rejected.
- Role: descriptive behavior
- Deciding authority: `nrpy/examples/dendro_bssn.py` and `nrpy/examples/dendro_fccz4.py`, `parse_args`; `nrpy/infrastructures/Dendro/checkpoint.py`, `output_checkpoint_cpp` (formulation metadata); `nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py`, momentum lowering; `nrpy/infrastructures/Dendro/solver_context.py`, diagnostic scheduling and node reduction; Dendro-GR `BSSN_GR/src/parameters.cpp`, native lapse, CAKO, eta, and default settings; `BSSN_GR/src/bssngr_main.cpp`, post-merger CAKO switch and time stepping; `BSSN_GR/src/grUtils.cpp`, `computeBHLocations`.
- Corroboration: none available; no CI job builds native Dendro-GR, and the helper's W/chi cross-restore rejection covers only the checkpoint identifiers.

AMR parameter parity is not proof of identical remesh histories. The generated
and native Nyquist paths both use all three components of the puncture
separation, but the two codes track the puncture centers with different
schemes, so the center histories that drive excision, puncture-centered AMR, and
the Nyquist history can differ. Compare the resulting meshes before interpreting
constraint differences as differences between evolution equations.

The `dendro-validation` workflow job runs
`nrpy/examples/tests/dendro_application_check.py` once per formulation and
configures these required checks: W and chi generation (byte-identical repeat
generation), complete builds, MPI TwoPunctures runs at FD4, FD6, and FD8 through
the first diagnostic output with the initial Hamiltonian-constraint norm ordered
by FD order (an ordering check, not a convergence test), finite diagnostics on
every run, constraint output, wave extraction with odd-m modes near zero and the
reflection relation C(l,-m) = (-1)^l conj C(l,m), apparent-horizon irreducible
masses matching the puncture ADM masses, a remesh test scheduled every four steps
(`BSSN_REMESH_TEST_FREQ = 4`) whose evolved state and node counts match a stored
reference, checkpoint and byte-identical diagnostic
restore, point reflection of the checkpointed puncture centers and their motion
along the puncture momenta, and rejection of both W/chi checkpoint-formulation
mismatches. The restore comparison normalizes GridInfo wall time, excludes only
the per-launch `dgr__PARAM_DUMP__*.toml` files under `dat/`, and compares every
other file under `dat/`, `bah/`, and `vtu/` byte for byte.
The stored-reference comparison of run A's evolved constraint, ADM, and horizon
values is a regression check on the evolution, not a correctness proof. Lambda
initialization is asserted through the stored reference, which includes the
step-0 Lambda constraint (and, for fCCZ4, the Z4 diagnostics). A stored value at
round-off level makes that comparison an upper bound: it fails only when the
value exceeds the absolute tolerance. The algebraic projection is asserted
through checkpointing: the checkpoint writer refuses to write, and a restore
refuses to read, a checkpoint whose projection residual is nonfinite or above
its fixed tolerance, so the residual is tested only when a checkpoint is written
or read, and only against that upper bound. Halo exchange between ranks is
covered by 1-, 3-, and 4-rank runs whose diagnostics must agree within a fixed
relative and absolute tolerance, with equal unexcised node counts; block
boundaries within a rank are covered only indirectly, through the numerical
checks and the stored reference. The job does not inspect field data, does not check
`alpha=W=sqrt(chi)` pointwise, and it runs none of the temporary direct
numerical comparisons or native-control comparisons above; those remain
review-time checks. The list above omits several helper checks that [Generated
Project CI](../../validation/generated-project-ci.md) describes: the closed-form
ADM energy and angular momentum, the vanishing ADM momentum and in-plane angular
momentum, the horizon reflection symmetry, the W-and-chi agreement at step 0,
the per-run sanity check, and the full list of rejection and warning cases.

Three limits bound what the comparisons show. The rank comparison of the
wave-mode and GridInfo files needs only one step in common with the compared
steps and compares rows pairwise, so a run that lacks later wave rows passes on
the rows both files have; the stored reference holds no wave values, and the
helper's other checks require the 21 wave files at two radii and step-checked
constraint rows. The "one shared mesh" of the FD-order comparison is shown by
equal element counts after the initial-grid remesh, which does not show equal
octant coordinates or levels, and the W-versus-chi comparison covers the node
count, `E`, `J_z`, and the Psi4 modes but not the constraint or horizon columns.
The hand-run capability test below tests each error as `err > 1e-9`, which is
false for a NaN, and its injected defects use zeroed or shifted data, so a host
that writes NaN would pass the affected axes.

Claim evidence:
- Claim: The `static-analysis` job analyzes the Python modules found by its file search, which skips `__init__.py` files and the paths `./project/`, `./build/`, `*/tests/*`, `*manga*`, and `./nrpy/examples/visualization_scripts/`, so the Dendro generator modules are analyzed and `dendro_application_check.py` is not; the helper generates each W and chi variant twice and compares the trees; the rank comparison of wave-mode and GridInfo files requires only one step in common and compares rows pairwise, the stored reference holds no wave values, the FD-order comparison shows one shared mesh only through equal element counts after the initial-grid remesh, and the W-versus-chi comparison covers the node count, `E`, `J_z`, and the Psi4 modes only.
- Role: CI behavior
- Deciding authority: [main.yml](../../../.github/workflows/main.yml), `static-analysis`; [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py), `Leg.generate_and_build`, `Leg.compare_runs`, `Leg.check_orders`, `Leg.compare_variants`, `Leg.check_reference`
- Corroboration: none available; no case in the helper fails these comparisons by design

Every tolerance in the helper is a literal in `dendro_application_check.py`.
The helper derives none of them; it ties the stored-reference tolerance to the
printed digits of the output columns. A pass shows agreement within those
thresholds, not an error bound for the evolution.

The helper generates only with Kreiss-Oliger dissipation enabled (the example
scripts set it), with `--fd-order 6`, and with neither `--ybs-gamma` nor
`--ybs-momentum`. It evolves the generated `pars/<stem>.toml` with profile
overrides and never reads `pars/q1.par.lowres.toml` or runs the printed run
commands, so the packaged q1 file is outside the CI evidence. Run-time orders 4
and 8 are exercised only as `BSSN_ELE_ORDER` values of a binary generated at FD6;
generation with `--fd-order 4` or `8`, and with either Yo et al. option, is
exercised by no job.

Claim evidence:
- Claim: `dendro-validation` configures the listed required application checks through `dendro_application_check.py`, covers the evolution and the scheduled remesh test through a stored-reference regression check, covers Lambda initialization through the stored step-0 constraint row and the algebraic projection through the checkpoint write-time and restore-time residual checks, covers halo exchange through rank-count agreement and block boundaries within a rank only through the numerical checks and the stored reference, and omits field-data inspection, the pointwise initial-lapse check, the temporary direct numerical comparisons, and native-control comparisons. Its stored-reference and checkpoint-residual comparisons are upper bounds where the stored value is at round-off level, and the helper derives none of its tolerances. The helper generates only with Kreiss-Oliger dissipation enabled, with `--fd-order 6`, and without `--ybs-gamma` or `--ybs-momentum`, and it never reads the packaged q1 parameter file or runs the printed commands.
- Role: CI behavior
- Deciding authority: [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py), `Leg.run_variant`, `Leg.check_run_a`, `Leg.check_reference`, `Leg.check_orders`, `Leg.compare_runs`, `Leg.run_negatives`; [main.yml](../../../.github/workflows/main.yml), `dendro-validation`; [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) and [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py), `main` (Kreiss-Oliger enabled in the script); [checkpoint.py](../../../nrpy/infrastructures/Dendro/checkpoint.py), `output_checkpoint_cpp`, write-time and restore-time projection-residual checks
- Corroboration: [Generated Project CI](../../validation/generated-project-ci.md), Dendro job description

Runtime results belong in active review or CI output, not as KB snapshots.

The hand-run capability test `nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp`
checks, against a chosen Dendrolib build, the host assumptions that the generated
kernels make: the scalar ABI, the padded block dimensions, the padding rule, the
unzip offsets, the variable-major x-fastest layout, the padded origin, and the
validity of the halo. NRPy does not build it, no generated project contains it,
and no CI job runs it; its README in the same directory gives the build and run
commands, to be repeated when the selected Dendrolib changes. `CAPTEST_ORDERS`
selects the element orders (a nonempty comma-separated list of positive even
integers), and `CAPTEST_INJECT` introduces one known defect per axis so that each
checker can be shown to fail; the README tabulates the values. The weekly
`dendrolib-canary.yml` run, not this program, is the CI check against Dendrolib's
`master`.

Claim evidence:
- Claim: The capability test checks the scalar ABI, padded block dimensions, padding rule, unzip offsets, x-fastest layout, padded origin, and halo validity against a chosen Dendrolib build; it is built and run by hand, appears in no workflow, and accepts `CAPTEST_ORDERS` and `CAPTEST_INJECT`; its error comparisons `err > 1e-9` are false for a NaN, and its injected defects use zeroed or shifted data, so a host that writes NaN would pass the affected axes.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp`, the axis checkers and the `CAPTEST_*` reads; `nrpy/infrastructures/Dendro/tests_infra/README.md`, build and run commands; [main.yml](../../../.github/workflows/main.yml) and [dendrolib-canary.yml](../../../.github/workflows/dendrolib-canary.yml), which contain no capability-test step.
- Corroboration: none available; the program is not part of any automated run.

The helper uses fixed output frequencies for its stored-reference profiles by
setting `BSSN_SCALE_VTU_AND_GW_EXTRACTION = false`; the generated production
parameter files enable native scaling. The helper reads the native GridInfo
CSV header and includes grid counts, timestep, and physical time in its
existing restart and rank comparisons. Rank comparisons exclude wall time and active MPI rank count. Exact restart
comparisons normalize only the GridInfo wall-time column and skip only the
per-launch parameter dumps; every other byte of every compared file must match.
These fixed-cadence checks do not establish that native frequency scaling is
correct.

Claim evidence:
- Claim: The helper sets `BSSN_SCALE_VTU_AND_GW_EXTRACTION = false` for its stored-reference profiles and reads the native GridInfo header; rank comparisons exclude wall time and the active rank count; the restart comparison normalizes only the GridInfo wall-time column, skips only the per-launch `dgr__PARAM_DUMP__*.toml` files under `dat/`, and requires every other file under `dat/`, `bah/`, and `vtu/` to match byte for byte.
- Role: CI behavior
- Deciding authority: [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py), `COMMON_OVERRIDES`, `trees_identical`, `Leg.run_variant`, `Leg.compare_runs`
- Corroboration: none available; the configuration shows the comparison rules, not any run outcome

## Sources

- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - complete BSSN application generation.
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - complete fCCZ4 application generation.
- [CMakeLists.py](../../../nrpy/infrastructures/Dendro/CMakeLists.py) - explicit generated source manifest.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - the RHS key-order guard.
- [state_h.py](../../../nrpy/infrastructures/Dendro/state_h.py) - `validate_registered_state`.
- [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py) - configured CI checks for both applications.
- [dendro_application_check_reference.py](../../../nrpy/examples/tests/dendro_application_check_reference.py) - stored run-A reference values.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - production evolution and service scheduling.
- [Constraints And Diagnostic Norms](constraints-and-diagnostic-norms.md) - momentum convention and comparable RMS diagnostics.
- [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md) - remesh path and limits of native parity.
- [main.yml](../../../.github/workflows/main.yml) - currently configured CI jobs.
- [Dendro-GR](https://github.com/paralab/Dendro-GR) - `BSSN_GR/src/parameters.cpp`, `bssngr_main.cpp`, `grUtils.cpp`, and `rhs.cpp`, the native settings that a comparison must record.
- [dendrolib-canary.yml](../../../.github/workflows/dendrolib-canary.yml) - weekly helper run against Dendrolib `master`.
- [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp) and [README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md) - hand-run host-assumption tests.
- [Code Test Policy](../../validation/code-test-policy.md) - permitted test changes and proof limits.
- [BSSN_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py) - momentum lowering.
- [checkpoint.py](../../../nrpy/infrastructures/Dendro/checkpoint.py) - formulation metadata and the residual checks.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Generated Project CI](../../validation/generated-project-ci.md)
- See also: [Project Assembly And Generating Functions](project-assembly-and-emitters.md)
