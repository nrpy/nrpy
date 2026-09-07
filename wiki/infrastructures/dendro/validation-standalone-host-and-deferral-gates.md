# Validation, Standalone Host, And Deferral Gates

> Explain the standalone host vehicle, the generated self-tests and Minkowski lifecycle gates, the CI coverage, and the Dendrolib pin and proven capability axes. · Status: provisional · Last reconciled: 09-07-2026
> Up: [Dendro](index.md)

## Summary

The Dendro infrastructure validates in two places with different reach. Owner
doctests in the emitter modules check emitted text against the registries in
process. A standalone host vehicle lets the emitted C++ compile and run eleven
CTest cases, ten generated self-tests and the Minkowski lifecycle, without any
Dendro-GR checkout. What neither
establishes is behaviour against the real Dendrolib host. The source pin and
the capability proof that gated that work are both recorded on this page; what
remains open is a build against a real Dendro-GR checkout, which is a deliberate
deferral rather than forgotten work.

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

The generated solver compiles against standalone host types in `dendro_standalone_host.h`
rather than the real host. That header is not a test fixture: the examples'
inline assembly copies it into every generated project through `copy_files`, so
it is part of the generated solver's build. It is the one emitted header NRPy
does not run through `clang_format`, exactly as BHaH ships
`simd_intrinsics.h`. The generated `tests/` directory registers ten CTest cases
covering the state registry, parameter registry, padding, memory offsets, upwind
selection, the right-hand side, initial data, exact-name selection, det/trace
enforcement, and constraint diagnostics, and the solver's own `CMakeLists.txt`
registers the Minkowski lifecycle as an eleventh: every gate it prints exits
nonzero on failure, so the numerical gates are part of the suite rather than a
demo someone has to remember to run.

Because the solver is compiled against standalone host types, it refuses to configure
inside a Dendro-GR tree that defines a real `dendro5` target; linking generated
sources built against the standalone-host types to the real host would be a silent type
mismatch.

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

### The sample parameter file is reference only

Parameter-file parsing belongs to the Dendro-GR host, which already owns it;
NRPy emits the parameters and their defaults and does not parse a file. That is
a decision rather than a gap. The entry point therefore refuses a supplied `-t`
file, naming the host as the owner, rather than appearing to apply values it
would ignore, and the emitted parameter table is commented out for the same
reason. The effective values are printed at startup on rank 0: the generated defaults,
except the perturbation wavelength, which is a length and so is derived from the
block extent and spacing the host was given.

### CI coverage

Continuous integration exercises this infrastructure through `codegen-ubuntu`,
which generates both projects, configures and builds them with CMake, and runs
their `ctest` suites. The symbolic and emitted-source contracts run as owner
doctests in the static-analysis job.

What that job establishes is that the emitted C++ compiles and that its own
gates pass against the standalone host: every numerical gate is NRPy's own
kernels checked against NRPy's own host declarations, which says nothing about
the real target. The route that would prove something about the real target is a
container image with Dendro precompiled, generating the solver inside it,
building against the real host, and evolving a small job whose results are
checked. That route is unblocked by the pin and capability records below but has
not been run, so it is recorded as a deferral rather than approximated.

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

### Gates that remain open

This page is the capability record, and it records no proof for the axes the
harness does not exercise. Those stay gates rather than omissions: the boundary flags, the
zip direction, the `ts::Ctx` right-hand-side contracts, remesh, the checkpoint
ABI, output selection, and the thread model. Each names the work that would
close it, and physical boundaries and the checkpoint ABI remain separate
qualified profiles carried on this page.

One deferred item is named rather than merely open: the `BlockGeometry` adapter
proof, meaning a single auditable host function that normalizes
`component_offset` and `pmin_padded` for a block, exercised by a two-block case
and an offset sentinel. Only the `standalone_host/dendro_standalone_host.h`
struct exists so far, so the adapter signatures stay frozen until a build
against a real Dendro-GR checkout runs.

What the pin does not establish is that the generated solver builds against the
real host. The solver still compiles against `dendro_standalone_host.h` and still refuses
to configure inside a Dendro-GR tree defining a real `dendro5` target. The
container route described above is now unblocked by these records, but it has
not been run.

## Sources

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
