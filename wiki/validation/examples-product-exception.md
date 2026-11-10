# Narrow Examples Product Exception

> Conditions under which a new command-line helper under `nrpy/examples/**` may test an NRPy-generated standalone product, and how the checked-in helpers compare with them. · Status: provisional
> Up: [Validation](index.md)

## Summary

A new command-line helper under `nrpy/examples/**` may test an NRPy-generated
standalone product only when all twelve conditions below hold. They include
that the tested behavior is observable only at the generated program's process
boundary, that no cheaper module-level check proves the same regression,
same-change CI wiring, and `shell=False` subprocesses with timeouts and bounded
output. The page also records how the checked-in SEOB and Dendro application
helpers compare with these conditions.

Like [Code Test Policy](code-test-policy.md), this page sets prospective
policy. Its provisional status reflects incomplete durable support for the
commissioned norms; the checked-in SEOB helpers do not meet every condition.

## Detail

Outside this exception, examples are not ordinary policy evidence. A new
command-line helper under `nrpy/examples/**` is allowed only when all twelve
conditions hold:

1. The target is an NRPy-generated standalone product. Invoking an example
   generator may be setup only; generator argument behavior belongs in an API
   check or doctest in the module that defines it unless it is inseparable from
   the behavior required of the generated program when it runs as its own
   process.
2. An external compiler, build/runtime toolchain, or launcher is used only when
   the generated program's required behavior needs it.
3. The required behavior is observable only at the generated program's process
   boundary: its command-line arguments, parameter file, runtime state, or
   parsed output.
4. No cheaper API check in the module that defines the regressed behavior, and
   no cheaper symbolic, source-oracle, or semantic check, proves the same
   regression.
5. The same change wires an explicit CI invocation and target build/run.
6. Python helper subprocesses use `shell=False` argument vectors. Any
   unavoidable fixed shell orchestration stays in audited CI YAML and never
   interpolates helper or user data.
7. Every helper subprocess has an explicit timeout, small fixed resource or job
   limits, checked failure, and bounded or capped captured output. The CI job
   has `timeout-minutes` or equivalent.
8. Cwd and environment restoration and output cleanup are unconditional.
9. Inputs and randomness are deterministic; generated output is CI-owned or
   temporary.
10. Network, installation, and external access are explicitly authorized.
    Controllable dependencies use an exact pin, tag, or commit; otherwise
    record resolved version and reason and label evidence non-reproducible
    beyond the named run.
11. The helper asserts concrete output, state, semantic property, or numerical
    criterion beyond status, file existence, file count, or import success. A
    status-only build or crash assertion on the generated program must also
    satisfy all six conditions of the [Status-Only Build Or Crash
    Gate](code-test-policy.md#status-only-build-or-crash-gate) and keep the
    narrow compatibility/crash claim.
12. The assertion names the exact validation layer proved.

If any condition fails, the helper is prohibited. This exception never
authorizes general pytest or unittest adoption in examples.

The checked-in SEOB helpers are finite-input numerical-and-argument precedents:
they use deterministic samples, argument-vector subprocesses with checked
failure, parsed stdout, and a numerical failure criterion, and the workflow
invokes them explicitly. They lack explicit subprocess timeouts, bounded or
capped output capture, unconditional cwd restoration and cleanup, and parts of
the current resource-limit policy; their CI jobs also lack `timeout-minutes`.
None of those gaps may be copied.

Evidence (descriptive behavior): [sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py), `run_sebob`, `process_input_set`; [sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py), `run_sebobv2`, `process_input_set`; [main.yml](../../.github/workflows/main.yml), `sebob-consistency-test` and `sebobv2-consistency-test`

The Dendro application helper passes argument vectors with `shell=False`,
gives every generation, configure, build, run, and rejection subprocess an
explicit timeout, sets build jobs and MPI ranks through options with defaults of
4, writes subprocess output to log files and prints only a bounded tail on
failure, never changes its own cwd, passes the generator a copied environment
whose `XDG_CACHE_HOME` is a new directory inside its work directory, removes its
work directory in `finally`, records the resolved tool and package versions,
and prints each assertion with its validation layer and numerical bound. Its
workflow job sets `timeout-minutes`.

Evidence (descriptive behavior): [dendro_application_check.py](../../nrpy/examples/tests/dendro_application_check.py), `run_logged`, `run_checked`, `tail`, `Report.check`, `Leg.record_versions`, `Leg.generate_and_build`, `Leg.run`, `main`; [main.yml](../../.github/workflows/main.yml), `dendro-validation`

## Sources

- [sebob_consistency_check.py](../../nrpy/examples/tests/sebob_consistency_check.py) - `run_sebob`, `process_input_set`
- [sebobv2_consistency_check.py](../../nrpy/examples/tests/sebobv2_consistency_check.py) - `run_sebobv2`, `process_input_set`
- [dendro_application_check.py](../../nrpy/examples/tests/dendro_application_check.py) - `run_logged`, `run_checked`, `tail`, `Report.check`, `Leg.record_versions`, `Leg.generate_and_build`, `Leg.run`, `main`
- [main.yml](../../.github/workflows/main.yml) - `sebob-consistency-test`, `sebobv2-consistency-test`, `dendro-validation`

## See Also

- Parent: [Validation](index.md)
- Depends on: [Code Test Policy](code-test-policy.md)
- See also: [Generated Project CI](generated-project-ci.md)
