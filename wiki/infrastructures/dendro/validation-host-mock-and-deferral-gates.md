# Validation, Host Mock, And Deferral Gates

> Explain generated-project verification, the scaffold scan, the mock host vehicle and its self-tests, the Dendro CI job, and the unproven Dendrolib pin and capability gates. · Status: provisional · Last reconciled: 09-04-2026
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

### Scoped CI

The `dendro-fccz4` job generates the module, builds it under an optimized
configuration with `-Wall -Wextra -Werror`, runs the ten self-tests, and runs
the Minkowski lifecycle gates under two MPI ranks. It covers only what needs a
compiler and an MPI runtime; the symbolic and emitted-source contracts run as
owner doctests in the static-analysis job. See
[Generated Project CI](../../validation/generated-project-ci.md) for the
configured job map and what that job does not establish.

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
ten self-tests, and completes its Minkowski lifecycle across two ranks. What it
does not establish is enumerated once, on
[Generated Project CI](../../validation/generated-project-ci.md).

Claim evidence:
- Claim: a passing `dendro-fccz4` run establishes deterministic regeneration, a warning-free optimized build under `-Wall -Wextra -Werror`, ten passing generated self-tests, and a completed Minkowski lifecycle across two MPI ranks; it establishes neither pointwise agreement with another NRPy backend, nor Kreiss-Oliger emitted-source compilation, nor any behavior of the real Dendro-GR host.
- Role: CI behavior
- Deciding authority: `.github/workflows/dendro-fccz4-validation.yml`, job `dendro-fccz4`
- Corroboration: `nrpy/infrastructures/Dendro/templates/FCCZ4_GR_tests_CMakeLists.in`, the registered `add_test` set; the tuple below records a local reproduction of that configured recipe, not a hosted-runner result
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=2 MPI ranks; error_path=not-run; options=--fd-order 4 --no-ko; date=09-04-2026`

## Sources

- [validation.py](../../../nrpy/infrastructures/Dendro/validation.py) - module docstring check list, `verify`
- [copy_tool.py](../../../nrpy/infrastructures/Dendro/copy_tool.py) - `--check`, `--verify-only`, `preconditions`
- [dendrolib_pin.json](../../../nrpy/infrastructures/Dendro/dendrolib_pin.json) - `status`, `required`, `gate`
- [dendrolib_capabilities.json](../../../nrpy/infrastructures/Dendro/dendrolib_capabilities.json) - `status`, `layout`, `padding`
- [dendro_mock.hpp](../../../nrpy/infrastructures/Dendro/host_mock/dendro_mock.hpp) - mock host types for the generated module
- [FCCZ4_GR_tests_CMakeLists.in](../../../nrpy/infrastructures/Dendro/templates/FCCZ4_GR_tests_CMakeLists.in) - the registered `add_test` set
- [dendro-fccz4-validation.yml](../../../.github/workflows/dendro-fccz4-validation.yml) - `dendro-fccz4`
- [FCCZ4_GR_README_md.in](../../../nrpy/infrastructures/Dendro/templates/FCCZ4_GR_README_md.in) - `Deferred (recorded, not dropped)`

## See Also

- Parent: [Dendro](index.md)
- Validated by: [Generated Project CI](../../validation/generated-project-ci.md)
- Depends on: [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
- Implements: [Code Test Policy](../../validation/code-test-policy.md)
- See also: [Generated Backend Comparison](../../syntheses/generated-backend-comparison.md)
- See also: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
