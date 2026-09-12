# New Infrastructure Conformance

> Rules a new NRPy infrastructure must satisfy to match the established ones, each with a right example, a wrong example, and a mechanical test. · Status: confirmed
> Up: [Infrastructures](index.md)

## Summary

Patterns match existing code. Before writing a mechanism into a new
infrastructure, find its counterpart in
`nrpy/infrastructures/{BHaH,ETLegacy,CarpetX,superB}` and determine whether its
requirements apply. Repeated compatible implementations strengthen a
convention, but no invented occurrence threshold decides authority. A mechanism
without precedent needs a concrete host or project requirement and belongs in
the design record.

Every rule below carries the rule, a right example from real NRPy code, a
wrong example drawn from a mistake actually made in the Dendro effort, and a
mechanical test where one exists. The wrong examples came from the Dendro
effort and document failure modes the rules prevent.

## Detail

### Conformance is one-way

A new infrastructure conforms to the established ones. Never propose changing
BHaH, ETLegacy, CarpetX, or superB to match a newcomer, and never offer that as
an option in a review brief: it invites a reviewer to justify the code under
review instead of measuring it. Divergence is permitted only where the new host
genuinely requires it, stated concretely. "Cleaner" is not a requirement.

### The generic layer carries no formulation name

**Rule.** An infrastructure's top-level modules are named for the artifact they
emit and contain no formulation name. Physics lives under
`<Infrastructure>/general_relativity/`.

**Right** — BHaH's generic layer: `BHaH_defines_h.py`, `main_c.py`,
`Makefile_helpers.py`, `CodeParameters.py`.

**Wrong** — a generic layer with formulation-named templates and
`Dendro-GR/FCCZ4_GR/` hardcoded into path construction. A second formulation
could not be lowered through that layer without editing it, which means the
abstraction did not exist.

**Test.**

```bash
grep -ri "<formulation>" nrpy/infrastructures/<Infrastructure>/*.py
```

Any match outside `general_relativity/` needs a demonstrated generic-layer
reason.

### Read the registries directly

**Rule.** An infrastructure reads `gri.glb_gridfcs_dict`,
`par.glb_code_params_dict`, and `cfc.CFunction_dict` at the point of use. It
does not copy them into a parallel record set.

**Right** — `nrpy/infrastructures/BHaH/BHaH_defines_h.py`:

```python
for cp_name, code_param in par.glb_code_params_dict.items():
```

The emitter names the registry it reads. `Makefile_helpers.py` reads
`cfc.CFunction_dict` the same way, and ETLegacy, CarpetX, and superB follow the
same pattern.

**Wrong** — an emitter that takes a snapshot of the registries as a parameter:

```python
def render_state_header(snapshot: FrozenNRPyDendroSnapshot) -> str:
    for fg in snapshot.gridfunctions:
```

This looks disciplined — immutable input, no global reads, easy to test — but it
is still wrong: it duplicates the authoritative registries into `Frozen*`
records and threads a `snapshot=` parameter through emitters only to hand back
values the registries already hold.

Pitfalls worth naming, because each one is what made the invention feel like an
improvement:

- *The parallel structure arrives with a virtue attached.* Immutability,
  purity, and determinism are real virtues, and they are the reason an invented
  layer reads as an upgrade rather than a divergence. The question is not
  whether the mechanism is good; it is whether NRPy already does this.
- *The invention hides behind a plausible seam.* Freezing before emitting reads
  as a lifecycle stage, not as a second registry, so nobody asks where the data
  came from.
- *Nobody checks, because no brief asks.* A review brief that never poses "does
  NRPy do it this way?" cannot surface the answer, however many reviewers read
  it.

**Test.**

```bash
grep -rn "glb_gridfcs_dict\|glb_code_params_dict\|CFunction_dict" \
  nrpy/infrastructures/<Infrastructure>/
grep -rn "snapshot\|Frozen[A-Z]" nrpy/infrastructures/<Infrastructure>/
```

The first search must identify each registry the emitter consumes. Every match
from the second search needs a demonstrated purpose other than duplicating an
authoritative registry.

### Names for the generated unit are function arguments

**Rule.** The name of the generated project, thorn, or solver is threaded as a
function argument. NRPy registers none of these as a `CodeParameter`.

**Right** — BHaH threads `project_name` and `exec_or_library_name`; ETLegacy
and CarpetX thread `thorn_name`.

**Wrong** — registering `Dendro_module_name` as a `CodeParameter` and then
hardcoding the directory anyway, so the registered parameter had no effect.

**Test.**

```bash
grep -rn "register_param.*_name" nrpy/infrastructures/<Infrastructure>/  # expect no unit-name parameters
```

### The host's vocabulary governs emitted identifiers

**Rule.** Namespaces, target names, directory names, and file prefixes follow
the host's own conventions, read from the host's source.

**Right** — Cactus says thorn, so ETLegacy says `thorn_name`. Dendro's
`BSSN_GR/CMakeLists.txt` header says "BSSN SOLVER", so Dendro says
`solver_name`; Dendro namespaces solvers by lowercase formulation
(`namespace bssn` in `BSSN_GR`, alongside `fluid`, `ode`, `solver`, `timer`),
so a generated Dendro solver does the same.

Claim evidence:
- Claim: Dendro-GR uses the solver name BSSN and lowercase `bssn` namespace; these host conventions govern the generated solver name and namespace, without prescribing an occurrence count.
- Role: normative rule
- Deciding authority: this page, `The host's vocabulary governs emitted identifiers`, Rule
- Corroboration: registered Dendro-GR source, `BSSN_GR/CMakeLists.txt` header and `BSSN_GR` namespace declarations establish the host vocabulary

**Wrong** — replacing the C++ namespace `fccz4::generated` with
`Dendro::generated` on the reasoning that "fccz4" is a formulation name and
must go. No `Dendro` namespace exists anywhere in Dendro. The original name was
correct by the host's convention; the actual defect was that `fccz4` was
hardcoded in the generic layer rather than threaded from the caller.

Removing a formulation name is not the same as making a layer generic. The test
is whether the name is *threaded* or *hardcoded*, not whether it appears in the
source.

### Generated filenames are infrastructure-prefixed

**Rule.** A generated file is named for the infrastructure that emits it, never
for the project instance.

**Right** — BHaH emits `BHaH_defines.h`, never `<project_name>_defines.h`.

**Precedence.** Where the host itself names solver files for the formulation,
the host-vocabulary rule above governs and this one yields: Dendro-GR ships
`bssnCtx.cpp` and `bssn_constraints.h`, so a generated Dendro solver emits
`<solver_stem>Ctx.cpp`. This rule still governs artifacts that are
infrastructure-generic rather than host-named. The host's vocabulary reaches
only host-side identifiers — the solver directory, namespace, executable,
context class and CMake variables. An NRPy-emitted kernel takes the name BHaH
and ETLegacy already use for that operation (`rhs_eval`, `constraints_eval`,
`enforce_detgbar_equals_detghat_trAzero`), prefixed with the stem the way
ETLegacy prefixes its thorn name.

### Module naming

**Rule.** Modules are named for what they emit or do.

**Right** — `BHaH_defines_h.py`, `main_c.py`, `write_checkpoint.py`,
`Makefile_helpers.py`.

**Wrong** — `project.py`, `freeze.py`, `manifests.py`, `validation.py`. The
current `coding_style.md` wording, "snake_case naming that directly describes
their purpose", was too weak to prevent any of them.

## Sources

- [BHaH_defines_h.py](../../nrpy/infrastructures/BHaH/BHaH_defines_h.py) - `par.glb_code_params_dict` iteration
- [Makefile_helpers.py](../../nrpy/infrastructures/BHaH/Makefile_helpers.py) - `cfc.CFunction_dict` iteration
- [main_c.py](../../nrpy/infrastructures/BHaH/main_c.py) - `register_CFunction_main_c`
- [dendro_fccz4.py](../../nrpy/examples/dendro_fccz4.py) - `main`, the inline project assembly
- [coding_style.md](../../coding_style.md) - `## Python Coding Style`, module naming

## See Also

- Parent: [Infrastructures](index.md)
- Depends on: [Infrastructure Code Style](infrastructure-code-style.md)
- Depends on: [Python Coding Style](../architecture/python-coding-style.md)
- See also: [Generated Output Boundaries](../architecture/generated-output-boundaries.md)
- Example: [Dendro Project Assembly And Emitters](dendro/project-assembly-and-emitters.md)
