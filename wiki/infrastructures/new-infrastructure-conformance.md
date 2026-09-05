# New Infrastructure Conformance

> Rules a new NRPy infrastructure must satisfy to match the established ones, each with a right example, a wrong example, and a mechanical test. · Status: confirmed · Last reconciled: 09-05-2026
> Up: [Infrastructures](index.md)

## Summary

Patterns match existing code. Before writing a mechanism into a new
infrastructure, find its counterpart in
`nrpy/infrastructures/{BHaH,ETLegacy,CarpetX,superB}` and count the instances.
Three or more is settled convention, and a new infrastructure ignoring it is
wrong by default. One instance is weak evidence and should be called out as
such. A mechanism with **no** instance anywhere is an invention: that absence
is evidence against it and belongs in the design record, not in a reviewer's
third-round discovery.

Every rule below carries the rule, a right example from real NRPy code, a
wrong example drawn from a mistake actually made in the Dendro effort, and a
mechanical test where one exists. The wrong examples are real and each passed
at least one review, which is the argument that the rule is needed.

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
`Makefile_helpers.py`, `CodeParameters.py`, with two incidental formulation
mentions across the whole layer.

**Wrong** — a generic layer with 152 `fccz4` occurrences, nine of twelve
templates fCCZ4-named, and `Dendro-GR/FCCZ4_GR/` hardcoded into path
construction. A second formulation could not be lowered through that layer
without editing it, which means the abstraction did not exist.

**Test.**

```bash
grep -rio "<formulation>" nrpy/infrastructures/<Infrastructure>/*.py | wc -l
```

Expect zero outside `general_relativity/`.

### Read the registries directly

**Rule.** An infrastructure reads `gri.glb_gridfcs_dict`,
`par.glb_code_params_dict`, and `cfc.CFunction_dict` at the point of use. It
does not copy them into a parallel record set.

**Right** — `nrpy/infrastructures/BHaH/BHaH_defines_h.py`:

```python
for cp_name, code_param in par.glb_code_params_dict.items():
```

The emitter names the registry it reads. `Makefile_helpers.py` reads
`cfc.CFunction_dict` the same way, and ETLegacy, CarpetX, and superB all follow
suit: four infrastructures, one pattern.

**Wrong** — an emitter that takes a snapshot of the registries as a parameter:

```python
def render_state_header(snapshot: FrozenNRPyDendroSnapshot) -> str:
    for fg in snapshot.gridfunctions:
```

This looks disciplined — immutable input, no global reads, easy to test — and
that is exactly why it survived six independent reviews. It is still wrong: it
duplicated three registries into `Frozen*` records, threaded a `snapshot=`
parameter through seven modules, and added 916 lines whose entire job was to
hand back values `gri.glb_gridfcs_dict` already held.

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
  nrpy/infrastructures/<Infrastructure>/ | wc -l    # expect non-zero
grep -rn "snapshot\|Frozen[A-Z]" nrpy/infrastructures/<Infrastructure>/ | wc -l  # expect zero
```

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
(`namespace bssn`, 14 occurrences in `BSSN_GR`, alongside `fluid`, `ode`,
`solver`, `timer`), so a generated Dendro solver does the same.

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
infrastructure-generic rather than host-named.

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
- [output_project.py](../../nrpy/infrastructures/Dendro/output_project.py) - `output_project` argument list
- [coding_style.md](../../coding_style.md) - `## Python Coding Style`, module naming

## See Also

- Parent: [Infrastructures](index.md)
- Depends on: [Infrastructure Code Style](infrastructure-code-style.md)
- Depends on: [Python Coding Style](../architecture/python-coding-style.md)
- See also: [Generated Output Boundaries](../architecture/generated-output-boundaries.md)
- Example: [Dendro Project Assembly And Emitters](dendro/project-assembly-and-emitters.md)
