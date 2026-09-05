# Lifecycle And Project Assembly

> Explain the Dendro generation transaction, the staged pure export of the FCCZ4_GR module, its CMake and manifests, and the safe installer into a Dendro-GR checkout. · Status: provisional · Last reconciled: 09-04-2026
> Up: [Dendro](index.md)

## Summary

A Dendro run has four ordered stages: a generation transaction snapshots the
process-global NRPy registries and restores them on failure; builders register
their kernels; the environment is frozen; and a pure exporter turns that frozen
snapshot into a complete `FCCZ4_GR` module. The exporter accepts no formulation
configuration, writes through a staged temporary sibling directory, and refuses
a target it did not itself produce. Installing that module into a real Dendro-GR
checkout is a separate, explicitly invoked tool that leaves a receipt and can
remove exactly what it installed.

## Detail

### Generation transaction

NRPy's registries are process-global. `dendro_generation_transaction` adds
rollback and test isolation on top of the production pattern of one profile per
fresh process: it verifies a safe clean starting state, snapshots the relevant
registries, and restores them if generation raises. `registry_digest` gives the
transaction a comparison value for that restore. The transaction is a
generation-time safety mechanism, not part of the generated product.

### Staged pure export

`output_project` is read-only over the registries and verifies them before and
after the staged commit. Before committing it renders the whole artifact set a
second time in memory and requires the two renders to be byte-identical, which
is what makes regeneration of an unchanged environment deterministic rather than
merely usually equal. It writes artifacts to a temporary sibling directory,
checks them against the snapshot, and then atomically replaces the target
project directory. An existing non-generated target is refused; a target
produced by a previous run of the exporter is detected by its
`GENERATION_RECEIPT.json` and may be replaced. `verify_generated_project`
re-checks a finished tree.

The emitted module is laid out as `Dendro-GR/FCCZ4_GR` inside the project
directory, with `generated/` holding the pure snapshot renders,
`include/` and `src/` holding the runtime context and entry point, `pars/`
holding the generated parameter profile, `tests/` holding the generated
self-tests, and `host_mock/` holding the mock host header the module compiles
against. The rendered tools — the verifier, the installer, and their wrappers —
are written at the *project* root rather than inside the module, so a reader
looking for `verify_generated_project.py` or `copy_into_dendro.py` should look
beside `Dendro-GR/`, not under `FCCZ4_GR/`.

### Template policy

Fixed templates under `templates/` may contain only structures the NRPy core
registries do not model naturally: context scaffolding, Dendro vector lifecycle
calls, CMake target declarations, `main()` structure, and host adapter
functions. A template must not contain a concrete state-field name, a physics
parameter default, an fCCZ4 expression, a finite-difference coefficient, a
numerical loop, or a source inventory; `$`-style placeholders are the only
profile-dependent content. The generated-project self-test template is the one
designated exception, so it may name evolved fields and carry its own loops
over a test block. The scaffold scan covers a fixed set — both READMEs, the module and test
CMake files, the context header and source, the entry point, the generated
preamble header, and the parameter profile — rather than every rendered
template; see
[Validation, Host Mock, And Deferral Gates](validation-host-mock-and-deferral-gates.md).

### Module CMake

`render_module_cmake` emits a module CMake that contains no field list, no
physics parameter, no finite-difference coefficient, and no numerical loop. The
generated source list is derived from the frozen CFunction set through
`render_nrpy_generated_sources_cmake`, so the build never hardcodes a source
inventory. The module compiles with `-Wall -Wextra -Werror` and requires only
MPI and a C++17 compiler. It also refuses to configure inside a tree that
defines a real `dendro5` target; the reason is recorded once, on
[Validation, Host Mock, And Deferral Gates](validation-host-mock-and-deferral-gates.md).

### Manifests

`all_manifests` renders the manifest set as pure JSON from the frozen snapshot:
sorted keys, stable arrays, exact rational coefficient strings, UTF-8 with LF
endings, no absolute source paths, and no timestamps inside hashed content.
Separate renderers cover the module record, generation parameters, stencils,
CFunctions, provenance, projection, initial data, diagnostics, and the file
inventory. Stencil operator records come from the builder sidecar, whose
coefficients were computed once from NRPy's `compute_fdcoeffs_fdstencl`, so
there is one coefficient source rather than a second stencil generator.

### Safe installer and receipt

`copy_tool` is the single source of the generated `tools/copy_into_dendro.py`.
It defaults to a dry run. Under `--execute` it verifies the staged module,
optionally backs up an existing installed module, stages the module in a
temporary sibling directory, writes the installation receipt into the staged tree,
atomically renames that tree into place, updates the root `CMakeLists.txt`
through a temporary file when `--apply-root-cmake-patch` is also passed, runs
post-install verification, and prints the exact changed paths. `--execute`
alone therefore installs the module without touching the root `CMakeLists.txt`. The
root CMake marker block is owned solely by the installer, which is the one
place that applies it.

`--remove` deletes only files the receipt identifies and removes only the exact
marker block. The complete installed file set is validated before any unlink,
so a refusal never leaves a half-removed module. Receipt paths are
module-relative, which makes an installed copy self-verifying independently of
the generated project it came from.

Claim evidence:
- Claim: `--execute` alone installs the module without editing the root `CMakeLists.txt`, which is applied only when `--apply-root-cmake-patch` is also passed; `--remove` deletes only files the installation receipt identifies and removes only the exact root CMake marker block, and the complete installed file set is validated before any unlink, so a refusal leaves the installed module intact rather than half-removed.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/copy_tool.py`, `execute`, `remove`, and `receipt_files`
- Corroboration: `nrpy/infrastructures/Dendro/copy_tool.py` module docstring, `Removal (--remove)` paragraph; note that the docstring's numbered step order is stale where it disagrees with `execute`
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3; backend=Dendro; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=--execute, --apply-root-cmake-patch, --remove; date=09-04-2026`

The installer and the manifests compute SHA-256 over generated artifacts to
detect local modification and to bind a removal to what was installed. These are
integrity checks over generated products, not source tracking; the distinction
is recorded once, on
[Frozen Snapshot And Pure Translators](frozen-snapshot-and-pure-translators.md).

### Artifact boundary

The exported module, its binaries, its manifests, and its receipt are generated
products, not source evidence. Cite the Python generators and the registry
symbols instead; see
[Generated Output Boundaries](../../architecture/generated-output-boundaries.md).

Claim evidence:
- Claim: the exported module, its binaries, its manifests, and its installation receipt are generated products rather than source evidence, so citations name the Python generators and registry symbols instead of emitted files.
- Role: normative rule
- Deciding authority: `wiki/architecture/generated-output-boundaries.md`, `Summary` and `Detail`
- Corroboration: `nrpy/infrastructures/Dendro/project.py`, `output_project` staged-write and banner behavior

## Sources

- [project.py](../../../nrpy/infrastructures/Dendro/project.py) - module docstring template policy, `output_project`, `verify_generated_project`
- [transaction.py](../../../nrpy/infrastructures/Dendro/transaction.py) - `GenerationTransaction`, `dendro_generation_transaction`, `registry_digest`
- [cmake.py](../../../nrpy/infrastructures/Dendro/cmake.py) - `render_module_cmake`, `render_nrpy_generated_sources_cmake`
- [manifests.py](../../../nrpy/infrastructures/Dendro/manifests.py) - `all_manifests`, `render_stencils_json`, `render_receipt_json`
- [copy_tool.py](../../../nrpy/infrastructures/Dendro/copy_tool.py) - transaction steps, `marker_lines`, `remove`, `--force-generated-replace`
- [FCCZ4_GR_README_md.in](../../../nrpy/infrastructures/Dendro/templates/FCCZ4_GR_README_md.in) - `What this module contains`, `Deferred (recorded, not dropped)`
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - generation entry point and command-line profile

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Frozen Snapshot And Pure Translators](frozen-snapshot-and-pure-translators.md)
- Validated by: [Validation, Host Mock, And Deferral Gates](validation-host-mock-and-deferral-gates.md)
- Contrasts with: [superB Lifecycle And Project Assembly](../superb/lifecycle-and-project-assembly.md)
- See also: [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
- See also: [Generated Backend Comparison](../../syntheses/generated-backend-comparison.md)
