# Validation, Standalone Host, And Deferral Gates

> Explain Dendro owner checks, generated standalone validation, real-host qualification, CI scope, and explicit runtime deferrals. · Status: provisional
> Up: [Dendro](index.md)

## Summary

Dendro has durable validation routes at three boundaries: production-module
owner checks, generated standalone CTest fixtures, and an opt-in real-host
qualification harness. These routes define reproducible checks; this page does
not preserve results, host revisions, dates, environment tuples, or generated
artifact inventories. A capability claim must be re-established against the
source and host checkout being reviewed.

The standalone route covers both BSSN and fCCZ4 without Dendrolib. The real-host
route exercises the generated fCCZ4 context with Dendro-GR and Dendrolib. General
physical boundaries, remeshing, local time stepping, checkpoint/restart, output
selection, GPU execution, and threaded kernels remain outside that route.

## Detail

### Owner checks and durable oracles

Emitter modules keep doctests beside production code and use the normal module
`__main__` runner. Emitters containing whole-file C++ templates inspect required
markers and clang-format guards without requiring a registered solver.

The initial-data and algebraic-constraint emitters compare registered
`CFunction.full_function` text with owner-local trusted generated-source
fixtures. RHS and constraint diagnostics use trusted symbolic-expression
dictionaries instead: their generated kernels are large products of common
lowering, while the expressions are the formulation-owned contract. Finite-
difference stencil behavior remains owned by the finite-difference validation
route.

Claim evidence:
- Claim: Dendro emitter owners validate small stable artifacts as generated source and validate RHS/constraint systems at the symbolic-expression boundary.
- Role: generated evidence
- Deciding authority: [initial_data.py](../../../nrpy/infrastructures/Dendro/general_relativity/initial_data.py), [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/Dendro/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), and [constraints_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py), their `__main__` validation routes
- Corroboration: [generic.py](../../../nrpy/helpers/generic.py), `validate_strings`; [validate_expressions.py](../../../nrpy/validate_expressions/validate_expressions.py), `compare_or_generate_trusted_results`

### Generated standalone vehicle

With `<PREFIX>_STANDALONE_HOST=ON`, a generated solver compiles against
`dendro_standalone_host.h`. The example generator copies that header and
`block_geometry.h` into the project, so they participate in the same build as
the generated kernels.

Generated CTest fixtures exercise registry consistency, parameter forwarding,
padding and offsets, derivative selection, RHS and initial-data calls,
algebraic enforcement, constraint diagnostics, and a Minkowski lifecycle. An
independent nonflat reference evaluates the canonical symbolic expressions on
deterministic binary64 samples and computes expected values without calling the
generated kernel or its flat adapter. Its bound is derived from the expression
graph, scales, spacing amplification, and valid rounding/reassociation effects,
not from a measured C++ error.

The address fixture uses unequal spacing, component offsets, sentinel regions,
centered and mixed derivatives, both upwind directions, zero-speed upwinding,
and Kreiss-Oliger response. Fault variants must make the checker reject the
corresponding bad address, halo, or derivative behavior. The lifecycle gates
check algebraic residuals, constraints, flat-state RHS, adapter agreement,
perturbation response, convergence, drift, and enforcement invocation. Any
failed gate exits nonzero and therefore fails CTest.

Claim evidence:
- Claim: generated Dendro projects contain standalone executable checks for registry, addressing, numerical-kernel, and lifecycle contracts for both formulations.
- Role: descriptive behavior
- Deciding authority: [self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/self_tests_cpp.py), `output_self_test_artifacts`; [general_relativity/self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py), `_nonflat_reference_cpp`; [general_relativity/main_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/main_cpp.py), `standalone_ctest_statements`
- Corroboration: [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py), generated solver/test CMake registration

### Runtime parameters

The real entry point accepts `-t FILE`. Rank zero reads the TOML text and
broadcasts it; the host TOML library parses it, and registry-derived bindings
populate `params_struct`. Unknown fields, malformed or nonfinite values, and
profile mismatch terminate through the parent communicator. The standalone
vehicle rejects parameter files. The real runner also accepts `--steps` and
`--dt`; its mesh and initial-data profile remain fixed by the qualification
vehicle.

Claim evidence:
- Claim: the real entry point binds registered TOML parameters, forwards kernel parameters through registered signatures, and terminates the MPI job on invalid input.
- Role: descriptive behavior
- Deciding authority: [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py), `_REAL_MAIN`; [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py), `output_toml_bindings`; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `_codeparameter_tail` and `_REAL_SOURCE`
- Corroboration: [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), analytic parameter-response check

### Real-host qualification

The opt-in host branch uses actual `ot::Mesh`, `ot::Block`, `ot::DVector`, and
`ts::Ctx` types. `block_geometry` normalizes padded allocation, component
offset, physical padded origin, and spacing. Component pointers retain their
bases; the generated block kernel applies the block offset once. The RHS
callback exchanges halos, evaluates generated block kernels, and zips the
result. Post-timestep projection applies to stage and accepted states according
to the selected host integrator's callback contract.

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
- Claim: Dendro supplies reproducible capability and generated-runtime harnesses that can qualify the selected real host without embedding host snapshots in the KB.
- Role: public interface contract
- Deciding authority: [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp), its geometry/value checkers and fault modes; [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), its transport, parameter, lifecycle, and failure checks
- Corroboration: [README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md), reproduction procedure; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), real-host adapter and callback implementation

### CI and open gates

Owner doctests route through static analysis. Generated standalone and real-host
build/runtime commands are local qualification routes; the configured GitHub
workflow does not run them. Adding either route requires explicit authorization
to change the protected workflow.

Physical boundary semantics, remeshing and state transfer, local time stepping,
checkpoint/restart ABI, output selection, GPU execution, and threaded kernels
remain open. The real-host vehicle prescribes analytic exterior data by physical
position and does not expose a general boundary-flag interface. Its application
qualification is fCCZ4-specific; BSSN retains the standalone route until a
real-host owner adds an equivalent qualification.

## Sources

- [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp) - real-context transport, parameter, lifecycle, and fault checks
- [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp) - host capability checks and fault injection
- [README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md) - reproduction procedure
- [block_geometry.h](../../../nrpy/infrastructures/Dendro/block_geometry.h) - shared `block_geometry_struct` contract
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - standalone and real-host contexts
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - TOML bindings
- [dendro_standalone_host.h](../../../nrpy/infrastructures/Dendro/standalone_host/dendro_standalone_host.h) - standalone host types
- [self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/self_tests_cpp.py) - generated self-test artifacts
- [general_relativity/self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py) - nonflat reference fixture
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - generic process shell
- [general_relativity/main_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/main_cpp.py) - GR lifecycle and CTest registration
- [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py) - solver/test CMake emission
- [cmdline_input_and_parfiles.py](../../../nrpy/infrastructures/Dendro/cmdline_input_and_parfiles.py) - default parameter file

## See Also

- Parent: [Dendro](index.md)
- See also: [Generated Project CI](../../validation/generated-project-ci.md)
- Depends on: [Project Assembly And Emitters](project-assembly-and-emitters.md)
- Implements: [Code Test Policy](../../validation/code-test-policy.md)
- See also: [Generated Backend Comparison](../../syntheses/generated-backend-comparison.md)
- See also: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
