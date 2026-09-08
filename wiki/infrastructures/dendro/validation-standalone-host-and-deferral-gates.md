# Validation, Standalone Host, And Deferral Gates

> Explain the standalone host vehicle, the generated self-tests and Minkowski lifecycle gates, the CI coverage, and the Dendrolib pin and proven capability axes. · Status: provisional · Last reconciled: 09-07-2026
> Up: [Dendro](index.md)

## Summary

Dendro validation now covers owner checks, the standalone test vehicle, and a
pinned real-host fCCZ4 qualification. Both standalone formulations retain eleven
CTest cases. The generated fCCZ4 solver also builds against pinned Dendro-GR and
Dendrolib, passes distributed transport checks, and completes a checked 100-step
Minkowski evolution on two active MPI ranks. Real-host coverage is limited to
fixed mesh, global RK4, double precision, one CPU thread per rank, and analytic
Minkowski exterior data; the broader gates remain explicit below.

## Detail

### Owner doctests

Every emitter module carries doctests in the production module and runs them
through the standard `__main__` runner. The three emitters that carry a whole
C++ file as a template string cannot be run without a fully registered solver,
so their doctests assert on the template instead: the markers are all present,
and every one of them sits inside a clang-format guard, which is what keeps it
in the shipped artifact.

### Trusted baselines

The small emitted kernels are compared against trusted `.cpp` baselines stored
beside their builders, so a change in the lowered text fails the owning module
rather than surfacing when someone next builds a project. Ten files at
finite-difference order 4, five per formulation: the four initial-data builders
-- Minkowski fill, smooth perturbation, ADM-to-evolved conversion and the
connection initialization -- and the algebraic-constraint enforcement. Each is
captured at the conformal factor its application ships, `chi` for fCCZ4 and `W`
for BSSN, which the file name records; the captured object is the registered
CFunction's `full_function`, as BHaH captures its own, so a baseline carries the
padded-block pointer bindings and the point loop along with the lowered
expressions, and every one is byte-identical to the source the corresponding
generated project ships.

The right-hand side and the constraint diagnostics get no such baseline: they
are hundreds of kilobytes of SymPy-lowered kernel each, which `coding_style.md`
excludes from golden-output files on both sensitivity and size grounds, and no
tracked oracle file in this repository exceeds about eighty kilobytes. They are
pinned the way that rule directs instead -- symbolically, through
`nrpy.validate_expressions.validate_expressions.compare_or_generate_trusted_results`
on the assembled
expression dictionaries, in the same `__main__` sweeps, at the same two shipped
profiles. Four trusted dictionaries, 0.3 to 1.6 kilobytes each, as
`nrpy/infrastructures/CarpetX/general_relativity/rhs_eval.py` pins its own. The
finite-difference order axis is left to `nrpy/finite_difference.py`'s own
oracles, because the assembled expressions carry no stencil; the generated
padding self-test only bounds the emitted padding below by the centered radius
and does not discriminate the reach.

Claim evidence:
- Claim: the Dendro initial-data and algebraic-constraint-enforcement builders capture ten trusted generated-source baselines at finite-difference order 4, five per formulation and each at the conformal factor its application ships, every one byte-identical to the source that project generates; the right-hand side and the constraint diagnostics receive no generated-source baseline and are pinned symbolically instead by four trusted expression dictionaries at the same two shipped profiles.
- Role: generated evidence
- Deciding authority: [initial_data.py](../../../nrpy/infrastructures/Dendro/general_relativity/initial_data.py) and [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/Dendro/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), their `__main__` generated-source sweeps; [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) and [constraints_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py), their `__main__` symbolic sweeps
- Corroboration: [generic.py](../../../nrpy/helpers/generic.py), `validate_strings`, which writes a missing baseline and raises on a mismatch; [validate_expressions.py](../../../nrpy/validate_expressions/validate_expressions.py), `compare_or_generate_trusted_results`, which does the same for the expression dictionaries
- Validation: `inspected=pass; generated=pass; built=not-applicable; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, SymPy 1.14.0; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=--fd-order 4, both formulations at the conformal factor each ships, Kreiss-Oliger off; date=09-06-2026`

### Standalone host vehicle

With `<PREFIX>_STANDALONE_HOST=ON`, the generated solver compiles against
standalone host types in `dendro_standalone_host.h`. That header is not a test fixture: the examples'
inline assembly copies it into every generated project through `copy_files`, so
it is part of the generated solver's build. It and the shared `block_geometry.h` are copied assets rather than formatted
emitter output, as BHaH ships `simd_intrinsics.h`. The generated `tests/` directory registers ten CTest cases
covering the state registry, parameter registry, padding, memory offsets, upwind
selection, the right-hand side, initial data, exact-name selection, det/trace
enforcement, and constraint diagnostics, and the solver's own `CMakeLists.txt`
registers the Minkowski lifecycle as an eleventh: every gate it prints exits
nonzero on failure, so the numerical gates are part of the suite rather than a
demo someone has to remember to run.

With the option `OFF`, the generated context includes real Dendrolib headers
and links the host's `dendro5` and `toml11::toml11` targets. The two build modes
select separate context implementations; merely finding a real target no longer
prevents configuration. Duplicate solver executable names are rejected.

The entry point runs a Minkowski lifecycle whose gates are the det/trace
residual, the maximum constraint violation, the flat-state right-hand side, the
flat-block adapter agreement, the perturbed right-hand-side response, the
observed convergence order, the 100-step drift, and the exact enforcement-pass
count; a run that clears all eight prints its `MINKOWSKI_OK` line, which is a
terminal banner rather than a ninth gate. The lifecycle runs as a registered
CTest case, so a failing gate fails the suite.

Claim evidence:
- Claim: the generated solver builds warning-free against the standalone host and passes its eleven CTest cases, one of which is the Minkowski lifecycle, which also completes on two MPI ranks with observed convergence order 3.977 for BSSN and 3.976 for fCCZ4 at finite-difference order 4; none of this touches the real Dendro-GR host.
- Role: descriptive behavior
- Deciding authority: [self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/self_tests_cpp.py), `output_self_tests_cpp`; [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py), `output_main_cpp`
- Corroboration: [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py), `output_tests_cmake`'s ten `add_test` registrations and `output_solver_cmake`'s `add_test(NAME <stem>_minkowski_lifecycle ...)`
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3, OpenMPI 4.1.6; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=1 and 2 MPI ranks; error_path=not-run; options=--fd-order 4 --no-ko; date=09-06-2026`

### Runtime parameters

The real entry point accepts `-t FILE`. Rank 0 reads the TOML text and broadcasts
it; the host TOML library parses it, and registry-derived bindings populate
`params_struct`. Unknown fields, malformed/nonfinite values, and profile
mismatch fail through a parent-communicator abort. The standalone vehicle still
rejects parameter files. The real runner also accepts `--steps` and `--dt`.
The qualification mesh and initial data are fixed; shared geometry/perturbation
registry defaults do not select another physical problem.

Claim evidence:
- Claim: the real entry point binds registered TOML parameters, forwards kernel parameters through registered signatures, and terminates the MPI job on invalid input; its mesh and initial data remain the fixed qualification problem.
- Role: descriptive behavior
- Deciding authority: [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py), `_REAL_MAIN`; [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py), `output_toml_bindings`; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `_codeparameter_tail` and `_REAL_SOURCE`
- Corroboration: [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), analytic `eta` response check; [README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md), `Generated real-host qualification`
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=GCC 13.3.0, Open MPI 4.1.6, CMake 3.28.3; backend=Dendro real host; precision=double; GPU=not-run; restart=not-run; distributed=2 active MPI ranks; error_path=rank-local NaN and invalid TOML/profile; options=FD4, KO off, fixed mesh; date=09-07-2026`

### CI coverage

Continuous integration exercises this infrastructure through `codegen-ubuntu`,
which generates both projects, configures and builds them with CMake, and runs
their `ctest` suites. The symbolic and emitted-source contracts run as owner
doctests in the static-analysis job.

That job covers the standalone vehicle. The real-host build and numerical
checks below were run locally against pinned checkouts; no real-host CI job is
configured by this change. A prebuilt-host CI route remains separate work.

### Host pin and proven capabilities

Both records that gated real-host work are now closed for the geometry the
generated solver reads. The generated project's `CMakeLists.txt` records the
proven Dendrolib commit, `246043709e806021fcfc011fe657b8bf964cae4c` of
`paralab/Dendro-5.01`, in `<PREFIX>_PROVEN_DENDROLIB_COMMIT`, and warns when
Dendro-GR's own `DENDRO_dendrolib_GIT_TAG` names a different one. A standalone
build defines no tag, so nothing is compared there. Dendro-GR declares that tag
as a cache variable before adding the solver subdirectory, so a generated
`set()` would be discarded, and comparing is the only form that can warn. What a
consumer does with the warning is set Dendro-GR's own
`DENDRO_dendrolib_GIT_TAG` to the proven commit, or re-run
`nrpy/infrastructures/Dendro/tests_infra` against the commit they intend to
build and record the result here before relying on the block-layout claims. The proof
taken against that commit is recorded here rather than in a JSON record beside
the generator: the pin travels with the project a consumer builds, and the
claims travel with the page a reader routes through.

The proof is a standalone MPI harness,
`tests_infra/dendrolib_capability_test.cpp`, which builds an adaptive
octree, unzips three linear fields at `dof` 3, and compares every padded point
against values recomputed from the block record rather than read back from the
call under test. Linear data is reproduced exactly by the unzip interpolation,
so any mismatch is a layout, origin, or halo defect rather than interpolation
error. Seven axes now read `proven`: the scalar ABI, the padded dimensions, the
padding rule, the unzip offsets, the variable-major x-fastest layout, the padded
origin, and halo validity. The octree-to-domain map is not a separate axis: it
is the helper the origin and layout axes recompute their expected values
through, so those two axes fail if it is wrong.

Two results change what the infrastructure may claim. Block padding is the
element order halved, which `ot::Block`'s own constructor sets, so padding is
not independently selectable. And padding 5 is reachable at element order 10
and was proven on this pin, so an eighth-order finite-difference profile is no
longer host-gated; the generated profiles still expose orders 2, 4 and 6, which
is now a generator limit rather than a host limit.

Each of the seven axes was exercised against a known-bad input before its
passing result was recorded: `CAPTEST_INJECT` names one defect per axis, and
each named axis reports `FAILED` under its own injection. A passing axis whose
checker has not been shown to fail is not evidence.

Claim evidence:
- Claim: against Dendrolib commit `246043709e806021fcfc011fe657b8bf964cae4c`, the unzipped block layout is padded variable-major and x-fastest, the padded extent is `eleOrder * 2^(regGridLev - blockLev) + 1 + 2 * padding` per axis, padding is the element order halved, padded index zero sits `padding` cells below the block node's lower corner, and in-domain halo points carry neighbour data; padding 5 is reachable at element order 10.
- Role: public/scientific contract
- Deciding authority: [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp), its `run_order` and per-axis checkers
- Corroboration: `Dendro-GR/BSSN_GR/src/bssnCtx.cpp`, `ptmin[0] = GRIDX_TO_X(...) - PW * dx` and the `getAllocationSzX/Y/Z` block reads, read in an upstream Dendro-GR checkout that this repository does not track
- Validation: `inspected=pass; generated=not-applicable; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=GCC 13.3.0, Open MPI 4.1.6, CMake 3.28.3; backend=Dendrolib; precision=double; GPU=not-applicable; restart=not-run; distributed=1 and 2 MPI ranks; error_path=fault injection on every axis; options=element orders 2, 4, 6, 8 and 10; date=09-06-2026`

### Generated real-host qualification

The generated fCCZ4 context uses actual `ot::Mesh`, `ot::Block`, `ot::DVector`,
and `ts::Ctx` types. `block_geometry` normalizes the padded allocation, per-
component offset, physical padded origin, spacing, and boundary flags. Pointer
arrays retain component bases; the generated block kernel adds the block offset
exactly once. The real RHS callback exchanges halos, evaluates each generated
block kernel, and zips the result. The pinned RK4 host calls `post_timestep` on
four stage states and the accepted state; projection runs there, never on the
RHS vector passed to `post_stage`. A callback failure aborts `MPI_COMM_WORLD`,
including inactive ranks, because the pinned integrator ignores return codes.

The generated solver built against Dendro-GR
`b3261e2a0d3457781b11d63ac5ab38375ffab93b` and Dendrolib
`246043709e806021fcfc011fe657b8bf964cae4c`. On two active ranks, the transport
oracle exercised 18 local blocks on one rank, 27 nonzero-offset blocks overall,
32,796 in-domain halo points, and 3,259 received nodes. Component-distinct affine
fields and nonconstant zip results agreed within `5.92e-12`. Offset and halo
faults on rank 1 fail the oracle; a rank-1 NaN aborts the whole job.

The 100-step Minkowski run reached time `0.10000000000000007`, with maximum RHS
`2.61e-11`, constraints `2.70e-11`, drift `2.76e-13`, and projection residual
`1.00e-15` (rounded upward). Drift must stay below `1e-11`, projection below
`1e-13`, and RHS/constraints below `256*epsilon(double)/h_min^2`, here
`1.31e-10`. The mesh-scaled bound accounts for interpolation roundoff amplified
by second derivatives. Exact step count, elapsed time, finite values, and
`1 + 5*steps` projection passes are checked. This is a Minkowski roundoff check,
not a real-host convergence or general-boundary result.

Claim evidence:
- Claim: the generated fCCZ4 solver builds against the two stated host pins and passes a two-active-rank fixed-mesh RK4 Minkowski run, with independently checked component offsets, physical padded origins, real halo transport, nonconstant zip, parameter response, and rank-local failure termination; broader runtime profiles are not qualified.
- Role: descriptive behavior
- Deciding authority: [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `_REAL_HEADER` and `_REAL_SOURCE`; [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py), `_REAL_MAIN`
- Corroboration: [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), transport oracle and fault modes; [README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md), pinned commands and measured results
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, SymPy 1.14.0, GCC 13.3.0, CMake 3.28.3, Open MPI 4.1.6; backend=Dendro real host; precision=double; GPU=not-run; restart=not-run; distributed=2 active MPI ranks; error_path=rank-local offset, halo, and NaN injection; options=fCCZ4 chi, FD4, KO off, element order 6, fixed mesh, RK4, OMP_NUM_THREADS=1; date=09-07-2026`

### Gates that remain open

General physical boundary conditions and their flag semantics, remeshing and
state transfer, LTS, checkpoint/restart ABI, output selection, GPU execution,
and threaded kernels remain open. The adapter copies boundary flags but this
Minkowski profile prescribes only analytic exterior data by physical position.
The standalone BSSN checks remain valid; this real-host numerical record is
specifically for fCCZ4. Real-host CI is not added here.

## Sources

- [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp) - real context transport, parameter response, and rank-local fault modes
- [block_geometry.h](../../../nrpy/infrastructures/Dendro/block_geometry.h) - shared `BlockGeometry` contract
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - `_REAL_HEADER`, `_REAL_SOURCE`
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - `output_toml_bindings`

- [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp) - `run_order`, `CAPTEST_INJECT` fault injection
- [README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md) - build, run, and checker-exercise instructions
- [dendro_standalone_host.h](../../../nrpy/infrastructures/Dendro/standalone_host/dendro_standalone_host.h) - standalone host types for the generated solver
- [initial_data.py](../../../nrpy/infrastructures/Dendro/general_relativity/initial_data.py) - the `__main__` trusted-baseline sweep
- [self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/self_tests_cpp.py) - `output_self_tests_cpp`
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - `output_main_cpp`
- [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py) - `PROVEN_DENDROLIB_COMMIT`, `output_solver_cmake`, `output_tests_cmake`
- [cmdline_input_and_parfiles.py](../../../nrpy/infrastructures/Dendro/cmdline_input_and_parfiles.py) - `generate_default_parfile`

## See Also

- Parent: [Dendro](index.md)
- See also: [Generated Project CI](../../validation/generated-project-ci.md)
- Depends on: [Project Assembly And Emitters](project-assembly-and-emitters.md)
- Implements: [Code Test Policy](../../validation/code-test-policy.md)
- See also: [Generated Backend Comparison](../../syntheses/generated-backend-comparison.md)
- See also: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
