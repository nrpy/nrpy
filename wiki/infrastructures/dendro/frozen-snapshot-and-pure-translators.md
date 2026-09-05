# Frozen Snapshot And Pure Translators

> Explain the Dendro freeze boundary, the immutable snapshot records and semantic hashes, and the translators that render state, parameter, and CFunction artifacts from the snapshot alone. · Status: provisional · Last reconciled: 09-04-2026
> Up: [Dendro](index.md)

## Summary

The Dendro backend splits generation into two halves separated by one
boundary. Before the boundary, builders register gridfunctions, CodeParameters,
and CFunctions into NRPy's mutable registries. At the boundary,
`freeze_nrpy_dendro_environment` copies those registered objects into an
immutable snapshot, validates the frozen environment's invariants, computes
independent semantic hashes, and seals Dendro role registration so a late
Dendro registration raises. Any other late registry mutation is caught rather
than prevented: the exporter re-checks the live registries against the snapshot
before and after its staged commit. After the boundary, every target translator is pure: it reads only
the snapshot, imports no equation module, and authors no field name, parameter
default, or expression of its own.

## Detail

### What the snapshot holds

`FrozenNRPyDendroSnapshot` aggregates byte-stable copies of the authoritative
registered objects as `FrozenGridFunction`, `FrozenCodeParameter`, and
`FrozenCFunction` records, together with the Dendro role metadata, the recorded
access captures, and a builder sidecar in `extras`. `FrozenHashes` carries the
semantic hashes described below. `assert_mutable_registries_match` re-checks the
live registries against the snapshot, which is what lets the exporter verify the
registries both before and after its staged commit.

### Semantic hashes and the module ABI

Freeze computes seven values over canonicalized record lines, each covering one
independent contract:

- `state_schema_hash` covers the ordered gridfunction records and the scalar
  contract.
- `parameter_schema_hash` covers the used `CodeParameter` records.
- `equation_hash` covers canonical right-hand-side expression digests and
  formulation choices.
- `stencil_hash` covers finite-difference operators, exact coefficients, and
  access offsets.
- `cfunction_api_hash` covers registered names, prototypes, roles, and the call
  graph.
- `cfunction_source_hash` covers normalized registered `full_function` bytes.
- `module_abi_hash` combines all preceding values.

`canonical_expression_digest` produces the per-expression value the equation
hash is built from. The module ABI appears in the banner of every generated
module file the verifier scans, so a generated tree states which frozen
environment produced it, and regenerating an unchanged environment reproduces
the same value. The verifier walks the whole module and exempts only JSON,
Markdown, and the verbatim-copied mock host header, so the parameter TOML
profile carries the banner alongside the C++ and CMake artifacts. Among the files under
`manifest/`, `module.json` and `provenance.json` record the ABI value and the
others do not.

Claim evidence:
- Claim: the ownership banner carrying the module ABI is required on every generated module file except JSON, Markdown, and the verbatim-copied mock host header, so the parameter TOML profile is covered; separately, `module.json` and `provenance.json` are the manifests that record the ABI value.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/validation.py`, `_check_banners`
- Corroboration: `nrpy/infrastructures/Dendro/manifests.py`, `render_provenance_json` and `render_module_json`
- Validation: `inspected=pass; generated=pass; built=pass; run=not-run; result_checked=not-run`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=--fd-order 4 --no-ko; date=09-04-2026`

Separating the hashes means an intentional rename or reorder moves `state_schema_hash` without being confused with a change
of expressions or of emitted source bytes; that hash is the recorded basis for
rejecting an incompatible checkpoint once the checkpoint ABI lands, which is a
deferred profile rather than something the generated module does today.

These are semantic identities computed over registered NRPy objects and emitted
artifacts. They are not source-tracking fingerprints for KB sources, and the KB
stores no digest values for them.

### Purity rule for translators

Each translator turns one part of the snapshot into generated files and is
constrained to that input:

- `gridfunction_output` renders the state header, component bindings, the state
  JSON, the name index, and the provenance header. It contains no list of
  fCCZ4 field names and imports no equation module: every name, order, count,
  and metadata value comes from the frozen gridfunction records.
- `CodeParameters_output` renders the parameter header, the parameter JSON, and
  a sample TOML profile from the frozen `CodeParameter` records. No physics
  parameter table is authored there.
- `CFunction_output` writes one source file per registered CFunction. The
  registered `full_function` becomes the file content plus a deterministic
  banner, declarations are emitted from the registered prototypes, and the CMake
  source list is derived from the same objects. No signature is reconstructed,
  no loop is wrapped, and no registered name is changed.

The practical consequence is that a formulation change is impossible to make in
the translator layer. To change what the generated module computes, a builder
must change what it registers before the freeze.

Claim evidence:
- Claim: the three target translators read only the frozen snapshot; none imports an equation module or contains an fCCZ4 field name, so a formulation change cannot be made in the translator layer.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/gridfunction_output.py`, `CodeParameters_output.py`, and `CFunction_output.py`, their module docstrings and import lists
- Corroboration: `nrpy/infrastructures/Dendro/freeze.py`, `FrozenNRPyDendroSnapshot`
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3; backend=Dendro; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=--fd-order 4 --no-ko; date=09-04-2026`

### Why the boundary exists

NRPy's registries are mutable process-global state, and a Dendro module's
correctness depends on the exact set, order, and content of what was registered.
Freezing converts that mutable state into a value that can be validated, hashed,
compared, and rendered repeatedly without risk that a later import or
registration silently changes the result.

## Sources

- [freeze.py](../../../nrpy/infrastructures/Dendro/freeze.py) - `freeze_nrpy_dendro_environment`, `FrozenNRPyDendroSnapshot`, `FrozenHashes`, `canonical_expression_digest`, `assert_mutable_registries_match`
- [gridfunction_output.py](../../../nrpy/infrastructures/Dendro/gridfunction_output.py) - `render_fccz4_state_hpp`, `render_state_json`, `render_fccz4_provenance_hpp`, `state_name_index`
- [CodeParameters_output.py](../../../nrpy/infrastructures/Dendro/CodeParameters_output.py) - `render_fccz4_parameters_hpp`, `render_parameters_json`, `render_parameter_toml_sample`
- [CFunction_output.py](../../../nrpy/infrastructures/Dendro/CFunction_output.py) - `render_CFunction_source`, `render_CFunction_declarations`, `CFunction_cmake_source_list`
- [ADR_dendro_names.md](../../../nrpy/infrastructures/Dendro/ADR_dendro_names.md) - `Decision`, `Consequences`

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions, Naming, Access Capture, And Loops](gridfunctions-naming-access-capture-and-loops.md)
- See also: [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
- See also: [C Function Registry](../../core/c-function-registry.md)
- See also: [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
