# Validation, Host Mock, And Deferral Gates

> Explain generated-project verification, the scaffold scan, the mock host vehicle and its self-tests, the Dendro CI job, and the unproven Dendrolib pin and capability gates. · Status: provisional · Last reconciled: 09-05-2026
> Up: [Dendro](index.md)

## Summary

The Dendro backend validates in three places with different reach. A rendered
verifier checks a generated tree against its own manifests and ownership rules.
A mock host vehicle lets the emitted C++ compile, run its self-tests, and
execute a Minkowski lifecycle without any Dendro-GR checkout, and a scoped CI
job runs exactly that. What none of them establishes is behaviour against the
real Dendrolib host: that is gated behind a source pin and a capability proof
that are both currently unrecorded, and the gate is deliberate rather than
outstanding work.

## Detail

### Generated-project verification

`validation.py` is the single source of the generated `tools/verify_generated_project.py`;
the exporter renders that file from this source with a generated-file banner.
All checks must pass for exit zero:

- every file listed in the generated file-hash artifacts exists and matches, and
  every file on disk is listed, noting that the three hash artifacts cannot list
  themselves;
- every generated module file carries the ownership banner and the module ABI,
  apart from a short exemption set recorded on
  [Frozen Snapshot And Pure Translators](frozen-snapshot-and-pure-translators.md);
- every manifest parses as JSON;
- the CMake generated source list equals, as a set, the source paths derived
  from every registered CFunction record, which is the operator-versus-inventory
  cross-check that keeps one source per registered CFunction; and
- no concrete EVOL name, bare or role-prefixed, and no NRPy-emitted numerical
  loop appears in a scaffold file.

The scaffold scan has a stated limit worth carrying: the EVOL-name scan is the
mandated one, the loop check is an additional heuristic keyed on the NRPy loop
helper's index names, and three prohibitions — a physics parameter default, a
finite-difference coefficient, and a loop written with other index names —
are not mechanically enforced and remain a review obligation.

### Mock host vehicle

The generated module compiles against mock host types in `dendro_mock.hpp`
rather than the real host. That header is not a test fixture: the exporter
copies it into every generated project, so it is part of the generated module's
build. The generated `tests/` directory registers ten CTest cases covering the
state registry, parameter registry, padding, memory offsets, upwind selection,
the right-hand side, initial data, exact-name selection, algebraic projection,
and constraint diagnostics.

Because the module is compiled against mock host types, it refuses to configure
inside a Dendro-GR tree that defines a real `dendro5` target; linking generated
sources built against the mock types to the real host would be a silent type
mismatch.

### No CI coverage yet

Nothing in continuous integration exercises this module. The symbolic and
emitted-source contracts run as owner doctests in the static-analysis job, but
the generated self-tests and the Minkowski lifecycle are run by hand.

A job that built the module against the mock host would establish only that the
emitted C++ compiles: every numerical gate would be NRPy's own kernels checked
against NRPy's own stub, which says nothing about the real target. The route
that would prove something is a container image with Dendro precompiled,
generating the module inside it, building against the real host, and evolving a
small job whose results are checked. That waits on the pin and capability
records below, so it is recorded as a deferral rather than approximated.

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
the generated project carries an explicit deferral record for boundaries rather
than a silent omission. Treating these as recorded gates rather than as missing
features is the point: the generated module states what it has not proven, so a
reader is not left inferring coverage from a green build.

### What the current evidence supports

A passing CI run establishes that the generated module regenerates
deterministically, compiles warning-free under an optimized build, satisfies its
ten self-tests, and completes its Minkowski lifecycle across two ranks. None of
that touches the real Dendro-GR host, so it establishes nothing about behavior
on the real target.

Claim evidence:
- Claim: the generated project ships ten CTest cases and a Minkowski lifecycle whose gates cover constraint violation, perturbed-RHS response, observed convergence order and 100-step drift; they are run by hand against the mock host and no continuous-integration job executes them.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/templates/FCCZ4_GR_tests_CMakeLists.in`, the registered `add_test` set
- Corroboration: `nrpy/infrastructures/Dendro/templates/fccz4_main_cpp.in`, the gate list in its header comment
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=2 MPI ranks; error_path=not-run; options=--fd-order 4 --no-ko; date=09-05-2026`

## Sources

- [validation.py](../../../nrpy/infrastructures/Dendro/validation.py) - module docstring check list, `verify`
- [copy_tool.py](../../../nrpy/infrastructures/Dendro/copy_tool.py) - `--check`, `--verify-only`, `preconditions`
- [dendrolib_pin.json](../../../nrpy/infrastructures/Dendro/dendrolib_pin.json) - `status`, `required`, `gate`
- [dendrolib_capabilities.json](../../../nrpy/infrastructures/Dendro/dendrolib_capabilities.json) - `status`, `layout`, `padding`
- [dendro_mock.hpp](../../../nrpy/infrastructures/Dendro/host_mock/dendro_mock.hpp) - mock host types for the generated module
- [FCCZ4_GR_tests_CMakeLists.in](../../../nrpy/infrastructures/Dendro/templates/FCCZ4_GR_tests_CMakeLists.in) - the registered `add_test` set
- [FCCZ4_GR_README_md.in](../../../nrpy/infrastructures/Dendro/templates/FCCZ4_GR_README_md.in) - `Deferred (recorded, not dropped)`

## See Also

- Parent: [Dendro](index.md)
- See also: [Generated Project CI](../../validation/generated-project-ci.md)
- Depends on: [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
- Implements: [Code Test Policy](../../validation/code-test-policy.md)
- See also: [Generated Backend Comparison](../../syntheses/generated-backend-comparison.md)
- See also: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
