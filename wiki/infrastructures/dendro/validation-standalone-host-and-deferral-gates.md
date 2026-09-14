# Validation, Standalone Host, And Deferred Tests

> Explain Dendro module checks, generated standalone validation, real-host tests, CI scope, and tests not yet implemented. · Status: provisional
> Up: [Dendro](index.md)

## Summary

Dendro has reproducible validation at three levels: tests beside the generating
modules, generated standalone CTest cases, and an opt-in real-host test program.
These routes define reproducible checks; this page does
not preserve results, host revisions, dates, environment tuples, or generated
file inventories. A capability claim must be re-established against the
source and host checkout being reviewed.

The standalone route covers both BSSN and fCCZ4 without Dendrolib. The real-host
route exercises the generated fCCZ4 context with Dendro-GR and Dendrolib. General
physical boundaries, remeshing, local time stepping, checkpoint/restart, output
selection, GPU execution, and threaded kernels remain outside that route.

## Detail

### Module Checks And Independent Reference Values

Small doctests stay beside their generating modules. Generated headers, sources,
and CMake interfaces are validated by generating and building complete C++
projects instead of by matching fragments of emitted text.

Equation modules keep trusted symbolic-expression dictionaries under
`nrpy/equations`; VE supplies numerical values for their free symbols. Dendro
RHS and constraint assembly is instead covered by the complete generated C++
project: the nonflat reference checks the large FD kernels numerically, while
the remaining executable sections test their host interfaces.

Claim evidence:
- Claim: equation-level VE stays with the defining equation modules, while complete generated C++ projects test the Dendro generating functions together.
- Role: generated evidence
- Deciding authority: [fCCZ4_constraints.py](../../../nrpy/equations/general_relativity/fCCZ4_constraints.py) and [kreiss_oliger_terms.py](../../../nrpy/equations/general_relativity/kreiss_oliger_terms.py), their VE definitions; [general_relativity/self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py), `output_self_test_artifacts`
- Corroboration: [validate_expressions.py](../../../nrpy/validate_expressions/validate_expressions.py), `compare_or_generate_trusted_results`; [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py), generated CTest registration

### Generated Standalone Test Executable

Every generated solver compiles against `dendro_standalone_host.h` by default.
The fCCZ4 solver exposes `FCCZ4_STANDALONE_HOST=ON` for this selection; BSSN is
standalone-only and exposes no host-selection option. The example generator
copies that header and `block_geometry.h` into the project, so they participate
in the same build as the generated kernels.

Generated CTest cases exercise registry consistency, parameter forwarding,
offsets, derivative selection and reach, RHS and initial-data calls,
algebraic enforcement, constraint diagnostics, and Minkowski evolution. A
nonflat reference evaluates the canonical RHS and constraint-diagnostic
expressions with multiprecision arithmetic on deterministic binary64 samples
and computes expected values independently of generated C++ execution. Its
bound is derived from each expression graph, scales, spacing amplification, and valid
rounding/reassociation effects, not from a measured C++ error.

The addressing test uses unequal spacing, component offsets, sentinel regions,
centered and mixed derivatives, both upwind directions, zero-speed upwinding,
and Kreiss-Oliger response. The GR upwind test perturbs cells at the recorded
reach and one point beyond it, proving the generated RHS both uses the declared
outer point and ignores anything outside it. Fault variants must make the checker
reject the corresponding bad address, halo, or derivative behavior. The evolution tests
check algebraic residuals, constraints, flat-state RHS, adapter agreement,
perturbation response, convergence, drift, and algebraic enforcement. Any
failed check exits nonzero and therefore fails CTest.

Claim evidence:
- Claim: generated Dendro projects contain standalone executable checks for registry consistency, addressing, numerical kernels, and Minkowski evolution for both formulations.
- Role: descriptive behavior
- Deciding authority: [self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/self_tests_cpp.py), `output_self_test_artifacts`; [general_relativity/self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py), `output_self_test_artifacts`; [general_relativity/main_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/main_cpp.py), `standalone_ctest_statements`
- Corroboration: [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py), generated solver/test CMake registration

### Runtime parameters

The real entry point accepts `-t FILE`. Rank zero reads the TOML text and
broadcasts it; the host TOML library parses it, and use-derived bindings
populate the block-RHS members of `params_struct`. Unknown fields, malformed or
nonfinite values, and profile mismatch terminate through the parent
communicator. Parameters belonging only to standalone qualification kernels
are not accepted TOML keys. The standalone test executable rejects parameter files.
The real runner also accepts `--steps` and `--dt`; its mesh and initial-data
profile remain fixed by the real-host test. [Project Assembly And Generating
Functions](project-assembly-and-emitters.md#parameter-selection-and-host-geometry)
explains the parameter selection and geometry.

Claim evidence:
- Claim: the real entry point binds the opted-in parameters used by its block-RHS CFunction, forwards them through that registered signature, rejects other parameter keys, and terminates the MPI job on invalid input.
- Role: descriptive behavior
- Deciding authority: [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py), `_REAL_MAIN`; [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py), `output_toml_bindings`; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `codeparameter_tail` and `_REAL_SOURCE`
- Corroboration: [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), analytic parameter-response check

### Real-host qualification

The opt-in host branch uses actual `ot::Mesh`, `ot::Block`, `ot::DVector`, and
`ts::Ctx` types. `block_geometry` normalizes padded allocation, component
offset, physical padded origin, and spacing. Component pointers retain their
bases; the generated block kernel applies the block offset once. The RHS
callback exchanges halos, evaluates generated block kernels, and zips the
result. Post-timestep projection applies to stage and accepted states according
to the callback sequence required by the selected host integrator.

`dendrolib_capability_test.cpp` derives expected padded geometry and values from
the selected checkout's block records. It checks that Dendrolib padding is half
the element order; the order-ten case is the padding-five probe needed by an
eighth-order finite-difference profile. It also checks scalar ABI, padded
extents, unzip offsets, variable-major x-fastest layout, padded origin, and
in-domain halo values. Component-distinct affine fields separate layout and
transport errors from interpolation error. Each checker has a corresponding
fault-injection mode; a checker is credible only when its fault makes the
qualification fail.

`runtime_integration_test.cpp` exercises generated fCCZ4 callbacks, component
offsets, padded origins, halo exchange, nonconstant zip, parameter response,
finite-value handling, and rank-local failure termination. The Minkowski route
checks its configured evolution invariants using mesh-scaled numerical bounds.
These checks qualify only the selected source and host checkout; the KB stores
neither a revision fingerprint nor the run outcome.

Claim evidence:
- Claim: Dendro supplies reproducible capability and generated-runtime tests that can qualify the selected real host without embedding host snapshots in the KB.
- Role: tested host behavior
- Deciding authority: [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp), its geometry/value checkers and fault modes; [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), its transport, parameter, evolution, and failure checks
- Corroboration: [README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md), reproduction procedure; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), real-host adapter and callback implementation

### CI And Tests Not Yet Implemented

Module doctests run through static analysis. The configured
`dendro-validation` GitHub job generates default fourth-order, KO-enabled
fCCZ4. It runs the generated nonflat RHS-and-diagnostics multiprecision reference check, then links
the same solver into a checksum-verified real host and runs one Minkowski step
on exactly two MPI ranks. The real state is initialized through Dendro-GR's
Minkowski routine. The job is capped at 30 minutes; the real run is capped at
five minutes and uses one OpenMP thread per rank. These are configured tests,
not a stored execution result.

General application boundary semantics, remeshing and state transfer, local time
stepping, checkpoint/restart ABI, output selection, GPU execution, and threaded
kernels remain open. The real-host test uses Dendrolib block-boundary flags
to prescribe constant analytic exterior data; it does not expose a general
application boundary-condition interface. Its application qualification is
fCCZ4-specific; BSSN retains the standalone route until its generating module adds
an equivalent qualification.

The numerical checks use analytic or property oracles rather than a frozen
runtime output: multiprecision evaluation of the canonical RHS expressions,
resolvable KO contributions, and roundoff-scaled Minkowski bounds. The
configured generation does not establish nondefault finite-difference orders,
convergence, distributed transport, long-time or nonlinear evolution, or broad
physics validation.

## Sources

- [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp) - real-context transport, parameter, evolution, and fault checks
- [main.yml](../../../.github/workflows/main.yml) - configured `dendro-validation` job
- [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp) - host capability checks and fault injection
- [README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md) - reproduction procedure
- [block_geometry.h](../../../nrpy/infrastructures/Dendro/block_geometry.h) - shared `block_geometry_struct` field definitions
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - standalone and real-host contexts
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - TOML bindings
- [dendro_standalone_host.h](../../../nrpy/infrastructures/Dendro/standalone_host/dendro_standalone_host.h) - standalone host types
- [self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/self_tests_cpp.py) - generated self-test source files
- [general_relativity/self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py) - nonflat reference test
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - generic process shell
- [general_relativity/main_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/main_cpp.py) - GR evolution and CTest registration
- [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py) - solver/test CMake emission
- [parfile.py](../../../nrpy/infrastructures/Dendro/parfile.py) - default parameter file

## See Also

- Parent: [Dendro](index.md)
- See also: [Generated Project CI](../../validation/generated-project-ci.md)
- Depends on: [Project Assembly And Generating Functions](project-assembly-and-emitters.md)
- Implements: [Code Test Policy](../../validation/code-test-policy.md)
- See also: [Generated Backend Comparison](../../syntheses/generated-backend-comparison.md)
- See also: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
