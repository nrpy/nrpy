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
route supports either generated context with Dendro-GR and Dendrolib. General
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

Every generated solver compiles against `dendro_standalone_host.h` when its
directory is configured as the top-level CMake project. The example generator
copies that header and `block_geometry.h` into the project, so they participate
in the same build as the generated kernels. Embedded mode instead builds the
production library against the host's `dendro5` target.

Generated CTest cases exercise registry consistency, parameter forwarding,
offsets, derivative selection and reach, RHS and initial-data calls,
algebraic enforcement, constraint diagnostics, and Minkowski evolution. A
nonflat reference evaluates the canonical RHS and constraint-diagnostic
expressions with multiprecision arithmetic on deterministic binary64 samples
and computes expected values independently of generated C++ execution. Its
bound is derived from each expression graph, scales, spacing amplification, and valid
rounding/reassociation effects, not from a measured C++ error.

The addressing test uses unequal spacing, component offsets, sentinel regions,
centered first and second derivatives, mixed derivatives, and Kreiss-Oliger
response. Its coefficient values come from an independent finite-difference
weight construction rather than the code-generation helper. It perturbs cells
at the recorded reach and one point beyond it, proving the generated RHS both
uses the declared outer point and ignores anything outside it. Fault variants must make the checker
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
callbacks divide ownership as follows:

- `rhs` owns the halo exchange, exterior values, traversal of every local
  block, and zip back to a packed vector.
- `rhs_blkwise` accepts already-unzipped component arrays and a list of local
  block identifiers. It performs no exchange or zip and writes only the
  selected block interiors, using each Dendrolib component offset.
- `rhs_blk` accepts one component-major block-local slab. It rebases that
  block's component offset to zero before calling the same numerical kernel.
- `pre_stage_blk`, `post_stage_blk`, and `pre_timestep_blk` are byte-preserving
  no-ops. `post_timestep_blk` applies the algebraic BSSN projection to one
  block-local slab with the same zero-offset rule as `rhs_blk`.

The direct runtime check calls these functions on real Dendrolib storage and
compares the block and whole-vector results with `B^0=0.01`, `eta=1`, and an
explicit nonzero RHS-magnitude requirement. The nonzero shift-driver RHS makes
the comparison sensitive to component routing and input selection. This proves
their numerical and storage behavior, but does not prove that a particular
Berger-Oliger or local time-stepping scheduler invokes them in the required
sequence. Scheduler-driven
stage projection remains unqualified until that sequence is observed in the
selected Dendrolib time integrator.

`dendrolib_capability_test.cpp` derives expected padded geometry and values from
the selected checkout's block records. It checks that Dendrolib padding is half
the element order for generated element orders 4, 6, and 8, corresponding to 2,
3, and 4 points per side. It also checks scalar ABI, padded
extents, unzip offsets, variable-major x-fastest layout, padded origin, and
in-domain halo values. Component-distinct affine fields separate layout and
transport errors from interpolation error. Each checker has a corresponding
fault-injection mode; a checker is credible only when its fault makes the
qualification fail.

`runtime_integration_test.cpp` exercises a selected generated formulation's
whole-vector and block callbacks, component offsets, block-local zero offsets,
selected-block writes, padded origins, halo exchange, nonconstant zip,
parameter response, block/whole projection equivalence, byte-preserving block
hooks, argument rejection, finite-value handling, and rank-local failure
termination. It runs on one rank to isolate local block addressing and on two
ranks to include distributed transport. The Minkowski route checks its
configured evolution invariants using mesh-scaled numerical bounds. These
checks qualify only the selected source and host checkout; the KB stores
neither a revision fingerprint nor the run outcome.

Claim evidence:
- Claim: Dendro supplies reproducible capability and generated-runtime tests that can qualify the selected real host without embedding host snapshots in the KB.
- Role: tested host behavior
- Deciding authority: [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp), its geometry/value checkers and fault modes; [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), its transport, parameter, evolution, and failure checks
- Corroboration: [README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md), reproduction procedure; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), real-host adapter and callback implementation

### CI Coverage And Remaining Tests

Module doctests run through static analysis. The checked-in
`dendro-validation` GitHub job generates BSSN and fCCZ4 projects for
finite-difference orders 4, 6, and 8 with KO dissipation enabled and disabled,
then runs every generated standalone CTest check. It adds the default profiles
to one Dendro-GR build with `NRPY_DENDRO_BUILD_DRIVERS=ON`, builds the current
production libraries and `nrpy_<formulation>_dendro_qualify` executables, and
runs `runtime_integration_test.cpp` on one and two MPI ranks for each
formulation. The job also requires each injected transport and callback defect
to fail with its expected diagnostic and runs one two-rank Minkowski step with
each qualification executable. Workflow configuration proves this check
sequence, not a latest successful run.

Claim evidence:
- Claim: `dendro-validation` configures the complete standalone formulation/order/KO matrix and real-host default-profile checks described above, without recording a run result in the KB.
- Role: CI behavior
- Deciding authority: [main.yml](../../../.github/workflows/main.yml), `dendro-validation`
- Corroboration: [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), real-host checks and injected defects; [general_relativity/self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py), generated standalone numerical checks

General application boundary semantics, remeshing and state transfer, local time
stepping, checkpoint/restart ABI, output selection, GPU execution, and threaded
kernels remain open. The real-host qualification uses Dendrolib block-boundary flags
to prescribe constant analytic exterior data; it does not expose a general
application boundary-condition interface.

The numerical checks use analytic or property oracles rather than a frozen
runtime output: multiprecision evaluation of the canonical RHS expressions,
resolvable KO contributions, and roundoff-scaled Minkowski bounds. The
generated matrix covers finite-difference orders 4, 6, and 8 with KO both
disabled and enabled. It does not establish remeshing, long-time or nonlinear
evolution, or broad physics validation.

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
