# Validation, Host Mock, And Deferral Gates

> Explain the mock host vehicle, the generated self-tests and Minkowski lifecycle gates, the absent CI coverage, and the unproven Dendrolib pin and capability gates. · Status: provisional · Last reconciled: 09-05-2026
> Up: [Dendro](index.md)

## Summary

The Dendro infrastructure validates in two places with different reach. Owner
doctests in the emitter modules check emitted text against the registries in
process. A mock host vehicle lets the emitted C++ compile, run ten CTest cases,
and execute a Minkowski lifecycle without any Dendro-GR checkout. What neither
establishes is behaviour against the real Dendrolib host: that is gated behind
a source pin and a capability proof that are both currently unrecorded, and the
gate is deliberate rather than outstanding work.

## Detail

### Owner doctests

Each emitter that can be exercised cheaply carries doctests in the production
module and runs them through the standard `__main__` runner: the naming
decorations, the loop helpers, the parameter-header and parameter-file
emitters, the CMake emitters, the state header and its enum ordering, the role
sidecar, the generation-parameter validator, and all four fCCZ4 builders. An
emitter that would need a full registered fCCZ4 environment to exercise carries
no runner, because an empty runner is not coverage.

### Mock host vehicle

The generated solver compiles against mock host types in `dendro_mock.h`
rather than the real host. That header is not a test fixture: `output_project`
copies it into every generated project through `copy_files`, so it is part of
the generated solver's build. The generated `tests/` directory registers ten
CTest cases covering the state registry, parameter registry, padding, memory
offsets, upwind selection, the right-hand side, initial data, exact-name
selection, algebraic projection, and constraint diagnostics.

Because the solver is compiled against mock host types, it refuses to configure
inside a Dendro-GR tree that defines a real `dendro5` target; linking generated
sources built against the mock types to the real host would be a silent type
mismatch.

The entry point runs a Minkowski lifecycle whose gates are the projection
residual, the maximum constraint violation, the flat-state right-hand side, the
flat-block adapter agreement, the perturbed right-hand-side response, the
observed convergence order, and the 100-step drift. It is run by hand.

Claim evidence:
- Claim: the generated solver builds warning-free against the mock host, passes its ten CTest cases, and completes its Minkowski lifecycle on one and two MPI ranks with observed convergence order 3.977 at finite-difference order 4; none of this touches the real Dendro-GR host.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/Dendro_self_tests_cpp.py`, `output_Dendro_self_tests_cpp`, and `nrpy/infrastructures/Dendro/Dendro_main_cpp.py`, `output_Dendro_main_cpp`
- Corroboration: `nrpy/infrastructures/Dendro/cmake_helpers.py`, `output_tests_cmake` `add_test` registrations
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3, OpenMPI 4.1.6; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=1 and 2 MPI ranks; error_path=not-run; options=--fd-order 4 --no-ko; date=09-05-2026`

### The sample parameter file is reference only

The generated profile has no parameter-file binding: the entry point refuses a
supplied `-t` file rather than appear to apply values it would ignore, and the
emitted parameter table is commented out for the same reason. The effective
values are the generated defaults, printed at startup.

### No CI coverage

Nothing in continuous integration exercises this infrastructure. The symbolic
and emitted-source contracts run as owner doctests in the static-analysis job,
but the generated self-tests and the Minkowski lifecycle are run by hand.

A job that built the solver against the mock host would establish only that the
emitted C++ compiles: every numerical gate would be NRPy's own kernels checked
against NRPy's own stub, which says nothing about the real target. The route
that would prove something is a container image with Dendro precompiled,
generating the solver inside it, building against the real host, and evolving a
small job whose results are checked. That waits on the pin and capability
records below, so it is recorded as a deferral rather than approximated. The
accepted cost is that a codegen change producing invalid C++ now surfaces only
when someone builds by hand.

### Deferral gates

Two records gate the real host integration, and both are open:

- the Dendrolib source pin is `UNPINNED`. No Dendrolib source is vendored here,
  and the record requires a full upstream commit to replace a moving branch
  default before the adapter signatures and the padding and element-order
  claims may be treated as settled; and
- the capability record is `UNPROVEN`. The block layout convention, scalar ABI,
  offsets, dimensions, ghost validity, and origin are all recorded as unproven,
  and the maximum proven padding is four, so an eighth-order profile that needs
  padding five is capability-gated and currently rejected.

Physical boundaries and the checkpoint ABI are separate qualified profiles, and
the generated solver's README carries an explicit deferral record for them
rather than a silent omission. Treating these as recorded gates rather than as
missing features is the point: the generated solver states what it has not
proven, so a reader is not left inferring coverage from a green build.

## Sources

- [dendrolib_pin.json](../../../nrpy/infrastructures/Dendro/dendrolib_pin.json) - `status`, `required`, `gate`
- [dendrolib_capabilities.json](../../../nrpy/infrastructures/Dendro/dendrolib_capabilities.json) - `status`, `layout`, `padding`
- [dendro_mock.h](../../../nrpy/infrastructures/Dendro/host_mock/dendro_mock.h) - mock host types for the generated solver
- [Dendro_self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/Dendro_self_tests_cpp.py) - `output_Dendro_self_tests_cpp`
- [Dendro_main_cpp.py](../../../nrpy/infrastructures/Dendro/Dendro_main_cpp.py) - `output_Dendro_main_cpp`
- [cmake_helpers.py](../../../nrpy/infrastructures/Dendro/cmake_helpers.py) - `output_tests_cmake`, `output_solver_cmake`
- [Dendro_README_md.py](../../../nrpy/infrastructures/Dendro/Dendro_README_md.py) - `output_solver_README_md`
- [Dendro_parameter_file.py](../../../nrpy/infrastructures/Dendro/Dendro_parameter_file.py) - `output_Dendro_parameter_file`

## See Also

- Parent: [Dendro](index.md)
- See also: [Generated Project CI](../../validation/generated-project-ci.md)
- Depends on: [Project Assembly And Emitters](project-assembly-and-emitters.md)
- Implements: [Code Test Policy](../../validation/code-test-policy.md)
- See also: [Generated Backend Comparison](../../syntheses/generated-backend-comparison.md)
- See also: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
