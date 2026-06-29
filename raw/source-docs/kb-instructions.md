# Instructions: Build a Self-Maintaining Markdown Knowledge Base From Scratch

These instructions describe, in full, how to recreate this repository's
plain-markdown, agent-oriented knowledge base (KB) in **another repo** —
including replacing that repo's existing `AGENTS.md`. They are a synthesis of
four independent instruction sets, keeping the strongest treatment of each topic.

The approach is **domain-agnostic**. Only the *content* (branch names, leaf
topics, glossary terms, cited source files) changes per repo. The *tree shape,
file roles, page contracts, citation rules, and lint checks stay identical*.

**Fixed scaffolding — reuse these names verbatim; do not rename per repo.** The
directories `wiki/` and `raw/` and the metadata files below are constants of the
pattern. Recreate them with these exact names in every repo:

| Fixed path | Role |
| --- | --- |
| `AGENTS.md` | root router (the only root-level KB doc) |
| `wiki/` | the compiled, model-owned page tree |
| `wiki/SCHEMA.md` | the rules/governance hub |
| `wiki/workflows.md` | ingest/query/lint procedures |
| `wiki/glossary.md` | canonical-terms hub |
| `wiki/lint/CHECKS.md` | the lint checklist |
| `raw/` | source manifest + relocated source docs |
| `raw/SOURCES.md` | the source manifest |
| `raw/source-docs/` | source notes moved out of the repo root |

Only the **branch directories and leaf filenames inside `wiki/`** (and the cited
paths in `raw/SOURCES.md`) are adapted to the target repo's subject matter.

The end state: an agent answers questions by starting at `AGENTS.md`, following
a few router links to one or two compiled leaf pages, and synthesizing from
them — instead of grepping the whole tree. If normal work requires whole-tree
search, the routing has drifted and must be repaired.

---

## 1. Two kinds of "layers" — do not conflate them

The word "layer" means two unrelated things here. Keep them separate.

**Ownership layers — exactly three, by design.** A *who-owns-what* split, not a
tree depth:

| Layer | Owner | Location | Editable |
| --- | --- | --- | --- |
| Immutable/living sources | human/build/system | code tree, config/param files, `raw/source-docs/`, selected logs | never edited just to make the wiki cleaner |
| Compiled wiki | model/maintainer | `wiki/` | yes — synthesized pages with citations and cross-links |
| Schema/governance | shared | `wiki/SCHEMA.md`, `wiki/workflows.md`, `wiki/lint/CHECKS.md` | rarely — defines contracts and checks |

**Navigation layers — as many as the subject needs.** The compiled wiki is
itself a deep, recursively nested tree of router directories above the content
leaves. There is **no maximum depth**; a large repo may need five or more router
hops before a leaf. The only fixed rule is: *context accumulates at the leaves,
while `AGENTS.md` and every `index.md` act as directory routers.*

`SCHEMA.md` must itself state this distinction explicitly, so future agents do
not read "three layers" as a depth cap.

Operating principles:
1. **Compile once, navigate forever.** Ingest each source once into structured
   pages; do not re-answer a query by chunk retrieval.
2. **Routers thin, leaves substantive.** All real facts live in leaves; each
   leaf cites its sources and links to neighbors.
3. **Depth on demand.** Keep splitting routers downward and promoting overgrown
   leaves into directories so every router stays scannable and every leaf stays
   focused (§3).

---

## 2. The KB is a navigation tree (the centerpiece)

`AGENTS.md` is the root directory of the KB — a file only because agents
auto-load it, conceptually a table of contents. It points to near-trunk branch
routers. Those point to deeper routers or to leaves. Interior routers may point
to still deeper routers. **Leaves are the only normal nodes that carry full
contextual detail.**

Intended route shape (arbitrary depth):

```text
AGENTS.md
  -> wiki/<near-trunk>/index.md            (which part of the repo am I in?)
       -> wiki/<near-trunk>/<area>/index.md        (which component/workflow?)
            -> wiki/.../<topic>/index.md            (which topic inside it?)
                 -> wiki/.../<topic>/<leaf>.md      (the actual answer)
  -> wiki/glossary.md          (global hub)
  -> wiki/lint/CHECKS.md       (global hub)
  -> raw/SOURCES.md            (source manifest)
  -> wiki/SCHEMA.md            (rules)
```

### 2.1 Node roles

| Node type | File pattern | Purpose | Carries facts? |
| --- | --- | --- | --- |
| Root router | `AGENTS.md` | Entry point and task router | No |
| Near-trunk router | `wiki/<branch>/index.md` | Route broad question classes | No |
| Interior router | `wiki/**/index.md` (any depth) | Route narrower topic groups | No |
| Leaf | `wiki/**/<topic>.md` (not `index.md`) | Synthesized contextual detail | Yes |
| Hub | `wiki/SCHEMA.md`, `wiki/workflows.md`, `wiki/glossary.md`, `wiki/lint/CHECKS.md` | Governance / canonical terms / checks | Yes, for its governance role only |
| Source manifest | `raw/SOURCES.md` | Provenance and source state | Source metadata only |

### 2.2 Navigation invariants

- `AGENTS.md` links primarily to near-trunk `index.md` nodes and global hubs —
  not to every leaf.
- Every directory under `wiki/` that contains pages has an `index.md` router
  (the `wiki/lint/` utility directory is the one allowed exception).
- Near-trunk and interior routers are broad directory routers, not fact pages.
- A router lists only the **next useful layer** (immediate children), never the
  whole subtree.
- Routers may link to child routers and child leaves.
- Unique technical facts belong in leaves, not in `AGENTS.md` or any `index.md`.
- A `Where Do I Start?` table may shortcut to important leaves, but the primary
  hierarchy still flows through routers; shortcuts enter the graph, not bypass it.
- Every leaf is reachable from `AGENTS.md` by following ordinary relative links
  through routers.
- A leaf gains child pages only by being **promoted into a router**, with its
  content moved down into child leaves.

### 2.3 The two structural invariants

- **Invariant 1 — every node is either a router or a leaf, never both.** A node
  *with children* is a router (pointers only, no `Detail`, no citable facts). A
  node *with content* is a leaf (full page contract, §7.3).
- **Invariant 2 — recursive split keeps the tree deep and the nodes small.**
  Depth is created on demand from both ends:
  - *From the top:* when a router's direct child list exceeds **~7 entries**,
    group its children into sub-routers (each a new `index.md` directory) and
    add a layer. Routers grow *downward*, never into long flat lists.
  - *From the bottom:* when a leaf outgrows one focused topic, **promote it** —
    turn `topic.md` into `topic/index.md` (router) and split its content into
    child leaves.

### 2.4 Placement rule (where does a given topic go?)

1. A major mental model for the project → a near-trunk directory.
2. A topic with several recurring subtasks/concepts/phases/failure modes →
   an interior directory with its own `index.md`.
3. Repeat step 2 until each route decision is small and obvious.
4. A single coherent body of context → a leaf.
5. A router that starts accruing unique facts → move those facts into a leaf.
6. A leaf that becomes a list of unrelated contexts → split into child leaves
   and promote the old file to an `index.md` router.

### 2.5 Recursive split test (router vs leaf, at a glance)

- If a page would mainly tell the reader *where to go next* → it is a router.
- If a page would *explain a coherent topic* with sources, status, constraints,
  and neighbor links → it is a leaf.
- If a leaf needs many unrelated headings to stay readable → it is too broad;
  split into a deeper directory.
- If a router needs paragraphs to be useful → it is carrying leaf content; move
  that prose into one or more leaves.

### 2.6 Depth is real, not padding

- **No padding layers.** Depth must come from genuine branching (a router with
  several children), never a chain of single-child directories. Be deep where
  the subject fans out; never artificially deep.
- **Do not flatten a domain because it "only" has a few pages today.** If it has
  multiple natural question classes or is likely to grow, create the directory
  node now and keep it sparse — but never create empty leaves.
- **Target shape:** each router scannable at a glance (≲7 links), each leaf one
  topic, and several `index.md` hops between root and a typical leaf.

A generic deep route, to calibrate (no domain names):

```text
AGENTS.md
  -> wiki/subsystems/index.md
       -> wiki/subsystems/<component-family>/index.md
            -> wiki/subsystems/<component-family>/<component>/index.md
                 -> wiki/subsystems/<component-family>/<component>/failure-modes/index.md
                      -> .../failure-modes/startup.md      (leaf)
                      -> .../failure-modes/steady-state.md (leaf)
                      -> .../failure-modes/shutdown.md     (leaf)
```

---

## 3. Directory layout

`AGENTS.md` is the **only** markdown KB document at the repo root. Existing
project files (`README.md`, `LICENSE`, etc.) stay as project files — do not turn
them into routers. The directory names `wiki/` and `raw/` and the metadata
filenames (`AGENTS.md`, `wiki/SCHEMA.md`, `wiki/workflows.md`, `wiki/glossary.md`,
`wiki/lint/CHECKS.md`, `raw/SOURCES.md`, `raw/source-docs/`) are reused verbatim;
only the branch/leaf names inside `wiki/` are repo-specific.

Minimal skeleton (create before large-scale ingest):

```text
AGENTS.md
raw/
  SOURCES.md
  source-docs/
wiki/
  SCHEMA.md
  workflows.md
  glossary.md
  lint/
    CHECKS.md
  <branch>/
    index.md
    <leaf>.md
```

A mature tree expands the minimum into a layered structure (names match the
target repo; the *shape* is what is fixed):

```text
AGENTS.md
raw/
  SOURCES.md
  source-docs/
    original-agents.md            # if replacing an existing AGENTS.md (§4)
    <imported-design-note>.md
    <selected-log-or-evidence>.md
wiki/
  SCHEMA.md
  workflows.md
  glossary.md
  lint/
    CHECKS.md
  architecture/
    index.md
    overview.md
    build-and-run.md
    data-structures.md
    generated-vs-handwritten.md
  subsystems/
    index.md
    <small-subsystem>.md          # flat leaf when the area is genuinely small
    <component-family>/
      index.md                    # interior router
      <component>/
        index.md                  # deeper router
        concepts.md               # leaf
        algorithms/
          index.md                # still deeper router
          initialization.md       # leaf
          runtime-loop.md         # leaf
          error-handling.md       # leaf
        known-issues.md           # leaf
  experiments/
    index.md
    validation/
      index.md
      current-status.md
      known-failures.md
    campaigns/
      index.md
      <campaign-name>.md
  operations/                     # optional near-trunk branch
    index.md
    deployment.md
    observability.md
```

This is a *scaffold*, not a mandate to create empty pages. Create a directory
only when it owns at least one real router or leaf; create a leaf only when you
can fill it with compiled, sourced context.

Bootstrap the skeleton mechanically:

```bash
cd <target-repo-root>
mkdir -p wiki/lint raw/source-docs
for b in architecture subsystems experiments <other-branches>; do mkdir -p "wiki/$b"; done
mkdir -p wiki/subsystems/<component-family>   # add interior dirs as the design needs
# preserve and relocate any existing root KB doc (see §4)
git mv AGENTS.md raw/source-docs/original-agents.md 2>/dev/null \
  || mv AGENTS.md raw/source-docs/original-agents.md
```

Choosing branches: pick near-trunk routes that match how maintainers ask
questions, then subdivide below them. Common branches: `architecture`
(lifecycle, build/run, generated-vs-handwritten boundaries, core data
structures), `subsystems` (component behavior and failure modes), `experiments`
(validation status, campaigns, known outcomes), `operations` (deploy, runbooks,
incidents), `api` (public interfaces and contracts), plus any near-independent
subproject. Always include the meta files: `SCHEMA.md`, `workflows.md`,
`glossary.md`, `lint/CHECKS.md`.

---

## 4. Replacing the existing `AGENTS.md`

The root file becomes a short router, but its durable content must not be lost.

1. Read the existing `AGENTS.md` completely.
2. Copy it verbatim to `raw/source-docs/original-agents.md` (or another clearly
   named source file).
3. Register that file in `raw/SOURCES.md` (status `frozen`, ingest `registered`).
4. Compile its durable, repo-specific instructions (build/test/debug steps, etc.)
   into the right wiki leaves — e.g. `wiki/architecture/build-and-run.md` or
   `wiki/operations/debugging.md` — citing the source and linking the leaves
   from the relevant routers. Mark the original-agents source `ingested` once done.
5. Replace root `AGENTS.md` with the one-screen router (§7.1).

Do not delete useful instructions just because the root file must become short —
the root becomes the entry point; durable facts move into compiled pages.

---

## 5. Status and ingest vocabularies (use these exact words)

**Page status** (every leaf header):

| Status | Meaning |
| --- | --- |
| `confirmed` | Directly supported by current registered sources and reconciled. |
| `provisional` | Plausible but not fully sourced or not fully verified. |
| `contested` | Sources disagree, or current behavior conflicts with an older log. |
| `stale` | Known to need re-ingest because a living source has moved. |

**Source status** (`raw/SOURCES.md`):

| Status | Meaning |
| --- | --- |
| `frozen` | A snapshot/imported note/immutable evidence; not meant to change. |
| `living` | May change; drift must trigger re-ingest of dependent leaves. |

**Ingest state** (`raw/SOURCES.md`):

| State | Meaning |
| --- | --- |
| `registered` | Tracked as an authority but not yet compiled. |
| `partial` | Some of it compiled (e.g. an aggregate code-tree row). |
| `ingested` | Its relevant facts are compiled into one or more leaves. |

Dates are ISO `YYYY-MM-DD`. A leaf's `Last reconciled` is the date its content
was last checked against its sources.

---

## 6. Citation and confidence rules

- Cite code by path **plus a stable symbol** — function, struct, macro, or
  header name. Cite docs by path plus a stable heading or document role.
- **Do not rely on generated line numbers.** If a line number appears in a local
  note, treat it as advisory only.
- A non-trivial claim with no source makes the page `provisional`.
- If sources disagree, mark the page `contested` and state which supersedes.
- If a living source changed since reconciliation, mark/update the page `stale`.
- Use ordinary relative markdown links. **Do not** use Obsidian-style `[[wikilinks]]`.
- A page that only says "see the design doc" or "see the whitepaper" fails the
  compiled-page contract.

Good source lines:

```markdown
- [main.c](../../main.c) - `main`
- [src/config.h](../../src/config.h) - `DEFAULT_TIMEOUT_MS`
- [design-notes.md](../../raw/source-docs/design-notes.md) - retry policy section
```

Avoid:

```markdown
- See the code.
- main.c:123
- [[Design Notes]]
```

---

## 7. File-by-file specifications (templates)

### 7.1 Root router — `AGENTS.md`

One screen. Routing only — a `Router` table to each branch and global hub, plus
a task-oriented `Where Do I Start?` table. No unique facts, no `Detail` section.

```markdown
# <Project> Knowledge Base

This repository carries a plain-markdown knowledge base for <project purpose>.
Start here, follow a few links, and synthesize from the compiled pages instead
of grepping the whole tree first.

## Router

| Go to | Use it for |
| --- | --- |
| [Architecture](wiki/architecture/index.md) | End-to-end flow, build/run shape, boundaries, and core data structures. |
| [Subsystems](wiki/subsystems/index.md) | Major functional areas and component behavior. |
| [Experiments](wiki/experiments/index.md) | Validation status, logs, campaigns, and known outcomes. |
| [Glossary](wiki/glossary.md) | Canonical spellings and meanings for recurring entities. |
| [Lint Checks](wiki/lint/CHECKS.md) | Mechanical and review checks for the KB. |
| [Sources](raw/SOURCES.md) | Immutable/living source manifest with mtimes and hashes. |
| [Schema](wiki/SCHEMA.md) | Ownership, page contracts, and maintenance rules. |

## Where Do I Start?

| Task | Read first |
| --- | --- |
| Understand how to build or run | [Build And Run](wiki/architecture/build-and-run.md) |
| Debug <common failure class> | [<Failure Page>](wiki/subsystems/<area>/known-issues.md) |
| Change <common subsystem> | [<Subsystem>](wiki/subsystems/<subsystem>/index.md) |
| Inspect validation status | [Experiments](wiki/experiments/index.md) |

Rules for maintaining this KB live in [wiki/SCHEMA.md](wiki/SCHEMA.md).
```

### 7.2 Router pages — `wiki/**/index.md`

A title, a one-line blockquote naming what the branch routes, a `Page | Go here
when...` table, and navigation back to parent/root. No `Detail`, no unique facts,
no `Sources` (unless the router itself carries an unusual durable claim — avoid).
Each row must state *when to choose that route*; never a bare link list. Link
only to pages that exist; update a router only when a route actually changes.

```markdown
# <Branch>

> Router for <branch purpose>.

| Page | Go here when... |
| --- | --- |
| [Overview](overview.md) | You need the high-level shape or lifecycle. |
| [Specific Topic](specific-topic.md) | You need a focused implementation area. |
| [Sub-area](sub-area/index.md) | You need everything under <sub-area>. |

Back to [AGENTS.md](../../AGENTS.md).
```

Nested router (add an `Up` link; fix the `../` depth so the back-link reaches
root):

```markdown
Up to [<Parent>](../index.md). Back to [AGENTS.md](../../../AGENTS.md).
```

### 7.3 Leaf pages — the page contract

Every non-router content page uses exactly this shape and section order:

```markdown
# <Title>

> One-sentence purpose. · Status: <confirmed|provisional|contested|stale> · Last reconciled: <YYYY-MM-DD>
> Up: [<Parent Router>](index.md)

## Summary

Short synthesized answer — useful even when read alone. A few sentences.

## Detail

Compiled prose: behavior, invariants, decisions, edge cases, current status.
Not a pointer dump, not pasted raw notes.

## Sources

- [path/to/source.ext](../../path/to/source.ext) - `stable_symbol_or_header`
- [source-doc.md](../../raw/source-docs/source-doc.md) - stable section/role

## See Also

- [Neighbor Page](neighbor.md)
- [Other Branch Page](../other-branch/page.md)
```

Leaf rules:
- Header is **two blockquote lines**: line 1 = purpose + `· Status: … · Last
  reconciled: …`; line 2 = `Up: <parent router link>`.
- `Summary`, `Detail`, `Sources`, `See Also` appear in that order.
- `Summary` is synthesized, not a list of links.
- `Detail` is compiled prose, not retrieval notes.
- `Sources` is **never empty**; one source per line.
- `See Also` links to the parent and **≥1 neighboring** wiki page (the `Up:`
  link does not count toward the neighbor minimum).
- The hubs `SCHEMA.md` and `glossary.md` are exempt from the neighbor minimum.

### 7.4 Source manifest — `raw/SOURCES.md`

A pointer/provenance ledger, not content. Distinguish aggregate sources, moved
source documents, cited code/config files, selected evidence directories, and
exclusions.

| Field | Meaning |
| --- | --- |
| Source | Relative path or aggregate name. |
| Provenance | What the source represents and where it came from. |
| Status | `frozen` or `living`. |
| Mtime | Last-modified timestamp at audit time. |
| Hash | Stable content hash (`sha256:<hash>`). |
| Ingest | `registered`, `partial`, or `ingested`. |

```markdown
# Source Manifest

> Pointer manifest for source material. Root-level documentation sources live
> under `raw/source-docs/` so `AGENTS.md` is the only root KB document. Code,
> config, fixtures, selected logs, and build inputs stay in place. Status is
> `frozen` when a source is meant not to change and `living` when drift must
> trigger re-ingest. Last audited: <YYYY-MM-DD>.

## Aggregate Sources

| Source | Provenance | Status | Mtime | Hash | Ingest |
| --- | --- | --- | --- | --- | --- |
| `<source tree>` | in-tree implementation files, primary executable record | living | <YYYY-MM-DD HH:MM:SS> | sha256:<hash> | partial |

## Source Documents Moved Below Root

| Source | Provenance | Status | Mtime | Hash | Ingest |
| --- | --- | --- | --- | --- | --- |
| `raw/source-docs/original-agents.md` | previous root agent instructions | frozen | <YYYY-MM-DD HH:MM:SS> | sha256:<hash> | ingested |

## Cited Code And Config Sources

Exact cited files are registered below; they may also be covered by an aggregate
row above.

| Source | Status | Mtime | Hash |
| --- | --- | --- | --- |
| `src/example.ext` | living | <YYYY-MM-DD HH:MM:SS> | sha256:<hash> |

## Exclusions

Build artifacts, generated binaries, object files, generated PDFs used only as
build output, transient plans, bulk scratch logs, and latest-snapshot reports
are excluded unless a selected evidence directory is registered explicitly.
```

Manifest rules: every cited authority appears here directly or through a clearly
named aggregate row; `living` sources trigger re-ingest when mtime/hash changes;
fix a wrong `frozen` source by adding a new source, not editing history.

Collect metadata mechanically (never use placeholder hashes):

```bash
date +%F                                   # audit date
date -r path/to/source "+%Y-%m-%d %H:%M:%S"   # (or: stat -c '%y' path/to/source)
sha256sum path/to/source                   # take field 1, prefix "sha256:"
# Deterministic aggregate hash for a whole code tree:
find src include -type f \( -name '*.c' -o -name '*.h' \) -print0 \
  | sort -z | xargs -0 sha256sum | sha256sum
```

### 7.5 `wiki/SCHEMA.md` — the rules hub

Holds the full maintenance contract so `AGENTS.md` can stay a router. Required
sections:

```markdown
# Knowledge Base Schema

> Rules for the self-maintaining markdown KB. · Status: confirmed · Last reconciled: <YYYY-MM-DD>

## Summary
## Commitments
## Ownership
## Page Contract
## Citation Rules
## Router Rules
## Canonical Terms
## Not In The KB
## Sources
## See Also
```

The `Summary` must state that the ownership "layers" are **not** navigation
layers — navigation is a router tree of unbounded depth. Use these commitments:

1. Sources are immutable for wiki purposes; wrong sources are corrected by
   adding/updating a source, not by rewriting history through the wiki.
2. Sources, compiled pages, and schema are separate layers with separate owners.
3. The model owns the compiled wiki: summarize, cross-reference, file, update
   neighbors, and lint.
4. Compile once into structured pages; do not re-answer every query by chunk
   retrieval.
5. Ingest one source at a time and trace its implications across the linked graph.
6. Link everything with ordinary relative markdown links.
7. Navigate by indexes from `AGENTS.md` through as many router layers as needed;
   if a query requires reading nearly everything, the indexes have drifted.
8. Lint for broken links, orphans, dead code references, stale claims, low
   confidence, and entity drift.
9. Start with a small, high-value source set and grow conventions deliberately.

### 7.6 `wiki/workflows.md` — procedures

**Ingest** (one source at a time):
1. Register the source in `raw/SOURCES.md` (provenance, `frozen`/`living`,
   mtime, hash, ingest state).
2. Decide which branch owns the compiled facts.
3. Update the owning leaf in synthesized prose; cite exact files + stable symbols.
4. From glossary entities the source touches, follow existing links one hop and
   update affected neighbors.
5. Update a router only when a new leaf or route exists.
6. Run `lint/CHECKS.md` on the touched subgraph. Keep one-off reports out of the
   committed KB unless they add a durable checklist rule.

**Query:** start at `AGENTS.md`; follow the narrow route through routers to the
owning leaf; read one or two leaves; answer from compiled pages. If a normal
question needs broad tree search, record or fix index drift.

**Living-source change:** recompute mtime/hash; re-ingest the affected leaf;
update one-hop neighbors sharing changed glossary entities; mark `stale` if it
cannot be reconciled immediately. Never edit source files to make the wiki cleaner.

**Debugging (evidence-first):** do not propose a root cause or patch for a
runtime failure until it has been reproduced or explicitly marked unreproduced.
If a fix depends on a cached/optimized/approximate/generated helper matching
runtime behavior, instrument the reproduced case and compare the exact
runtime-visible quantities before changing code. Keep durable evidence as a
frozen source under `raw/source-docs/`; compile only the durable lesson.

### 7.7 `wiki/glossary.md` — canonical-terms hub

```markdown
# Glossary

> Canonical spellings for recurring code, domain, and experiment entities. · Status: confirmed · Last reconciled: <YYYY-MM-DD>

## Summary

Use these spellings in compiled pages to reduce entity drift. When a source uses
a looser spelling, compile it back to the canonical term here.

## Terms

| Term | Meaning |
| --- | --- |
| `<Term>` | Canonical meaning and stable context. |

## Sources

- [<source>](../<source>) - stable symbol or header

## See Also

- [Schema](SCHEMA.md)
```

Add or update an entry when: a term recurs across branches; a symbol has multiple
spellings; a generated name encodes domain meaning; or a user-facing name differs
from an internal function/file name. As a hub it is exempt from the neighbor
minimum.

### 7.8 `wiki/lint/CHECKS.md` — the lint checklist

Treat the KB like code. **Mechanical checks** (must pass for touched pages):
- *Broken links*: every relative `.md` link target exists.
- *Dead code references*: every cited file exists; a cited symbol appears in it.
- *Orphans*: every `wiki/**/*.md` is reachable from `AGENTS.md` through routers.
- *Directory routers*: every content-bearing directory under `wiki/` has an
  `index.md` (except exempt utility directories such as `wiki/lint/`).
- *Router discipline*: `AGENTS.md` and every `index.md` are routing lists only,
  with no `Detail` section and no duplicated leaf facts.
- *Root discipline*: `AGENTS.md` points to near-trunk routers and global hubs;
  task shortcuts enter the directory graph rather than bypass it.
- *Leaf edges*: every leaf links to its parent and ≥1 neighbor; `SCHEMA.md` and
  `glossary.md` hubs are exempt from the neighbor minimum.
- *Provenance*: every leaf has a non-empty `Sources` section.
- *Source registration*: every cited authority appears in `raw/SOURCES.md`
  directly or through an aggregate row.
- *Staleness*: no listed `living` source has a newer mtime/hash than the leaf's
  `Last reconciled` date.
- *Transient size checks*: generate token counts on request, but do not commit
  token-count reports (counts drift with every documentation edit).

**Review checks** (flag, do not block):
- *Naming drift*: reconcile split spellings through `glossary.md`.
- *Entity edges*: if a term appears in five or more leaves but no page links its
  relevant neighbor, add the missing edge.
- *Contradictions*: preserve both sources, link them, state which supersedes.
- *Low confidence*: mark `provisional`/`contested`/`stale` rather than hide it.

Carry an `Exclusions` list (§8) and a short `Current Baseline` note recording
that no build artifacts or scratch dumps are committed.

---

## 8. Constraints — what must NOT become a KB page

Do not compile pages for, and exclude from lint/size checks:

- Candidate/ranking/plan artifacts and transient working notes.
- Build products and LaTeX by-products (`.aux`, `.log`, `.toc`, `.out`, `.fls`,
  `.fdb_latexmk`), and generated PDFs used only as build output.
- Object files, executables, and other binaries.
- Bulk data dumps and bulk scratch logs — **unless** a specific evidence
  directory is registered in `raw/SOURCES.md` because a compiled page depends on it.
- "Latest"/token-count snapshot reports (they drift with every edit).

These are sources at most (often not even that); never compiled pages.

`raw/source-docs/` is for source notes that should not live at repo root: design
references, imported investigation summaries, frozen debug evidence, durable
experiment logs, and documents moved out of root. Bad candidates: one-off
prompts, accepted plan transcripts, generated by-products, rendered PDFs when the
`.tex` exists, uncited bulk logs.

---

## 9. Build procedure from scratch (ordered passes)

Build in passes; do not write every leaf before the router tree exists.

1. **Inventory the repo.** Identify durable sources: the existing `AGENTS.md` and
   other root docs; build/run docs (`README`, build files, CI config, manifests,
   deploy files); major source directories and near-independent components;
   generated-code boundaries; config/schema/migrations/fixtures; important issue,
   experiment, validation, or incident logs. Register only high-value sources first.
2. **Design the router tree.** Choose near-trunk routes from the repo's
   conceptual shape and how maintainers ask questions. Sketch at least one
   interior layer per branch; add more wherever a branch would otherwise point to
   many unrelated leaves (§2.4).
3. **Create the skeleton** (§3): `wiki/`, `wiki/lint/`, `raw/`,
   `raw/source-docs/`, every planned router directory, and only the first leaves
   you can fill. No empty leaves.
4. **Preserve and replace `AGENTS.md`** (§4): copy the old one to
   `raw/source-docs/original-agents.md`, register it, compile durable facts into
   leaves, then install the root router.
5. **Write the meta files** (they define the contract): `wiki/SCHEMA.md` and
   `wiki/lint/CHECKS.md` first, then seed `wiki/glossary.md` and `wiki/workflows.md`.
6. **Build `raw/SOURCES.md`**: register the code tree as an aggregate row, then
   the specific docs/files you intend to cite, with real mtimes/hashes (§7.4).
7. **Ingest one source at a time** (§7.6): write/extend the owning leaf at the
   correct depth, add citations and `See Also` edges, update glossary terms,
   update routers only when a route changes.
8. **Write routers bottom-up**: each interior `index.md`, then each branch
   `index.md`, then root `AGENTS.md` — each pointing only at its direct children.
9. **Test routing with real maintainer questions.** Each must route from
   `AGENTS.md` through routers to a leaf without whole-tree search, e.g.:
   - How do I build and run this? Where does execution start?
   - Where is configuration parsed?
   - What are the major subsystems and the most important data structures?
   - What generated files or external systems are involved?
   - What are the known failures, caveats, or current validation results?
   - Where do I start for the most common debug task?
   If a question cannot be routed, add router entries and leaves.
10. **Lint** (§7.8) the touched subgraph and fix every mechanical failure. Fix
    router drift immediately if navigation feels broad or vague.

Grow one source at a time, adding router layers whenever a source reveals durable
internal structure. Start small; do not stop at a shallow outline if the project
has deeper structure.

---

## 10. Maintenance rules

**When adding or changing a page:** update `Last reconciled`; update `Sources`;
update `See Also`; update the parent router only if the route shape changed;
update `raw/SOURCES.md` if sources were added/changed; run touched-subgraph lint.

**When a living source changes:** recompute mtime/hash; re-ingest affected
leaves; follow glossary/entity links one hop and update neighbors; mark
unresolved conflicts `contested` rather than silently choosing a convenient answer.

**When deleting or moving pages:** update routers; update inbound links; update
`raw/SOURCES.md` if source paths changed; confirm no orphan pages remain. Before
overwriting or deleting a page, confirm it is what its description claims; if the
content contradicts that, surface it rather than proceeding.

---

## 11. Anti-patterns to avoid

- Root-level KB docs other than `AGENTS.md`.
- Routers that carry unique facts or a `Detail` section.
- Leaves without sources, or pages that only point to a source document.
- Unregistered evidence logs; citing unregistered authorities.
- Generated-artifact churn committed under `raw/` or `wiki/`.
- Claims based on old/generated line numbers.
- "Latest" snapshots committed to the KB.
- Fix recommendations based only on plausible theories (not reproduced).
- Broad grep-first answering when a router route exists.
- Single-child padding layers, or flattening a deep domain into one giant router.

---

## 12. Writing style

- Write what is true, why it matters, and where it comes from.
- Prefer compiled prose synthesis over lists of facts; keep summaries short and
  directly useful.
- Prefer stable names over line numbers; keep exact code identifiers in backticks.
- Preserve uncertainty explicitly behind page/source status.
- Use ASCII punctuation for tokens with symbolic forms unless the repo already
  uses another convention.
- Use tables for routers and manifests.
- Avoid long quotations, raw grep dumps, marketing copy, and tutorials inside
  routers.
- Do not overfit to one investigation: if a debug finding yields a general rule,
  put the rule in `workflows.md`, `SCHEMA.md`, or the relevant subsystem leaf,
  and keep detailed evidence in `raw/source-docs/`.

---

## 13. Acceptance checklist

The new KB matches this approach when:

- [ ] `AGENTS.md` is the only root KB doc, a one-screen router, and links to
      `wiki/SCHEMA.md`, `wiki/lint/CHECKS.md`, and `raw/SOURCES.md`.
- [ ] Every content-bearing `wiki/` directory has a thin `index.md`; every node
      is either a router or a leaf, never both.
- [ ] Routing works at arbitrary depth: `AGENTS.md -> near-trunk -> child
      router(s) -> leaf`, and each router lists only the next useful layer.
- [ ] The tree is deep where the subject fans out: no router exceeds ~7 child
      links, no leaf carries more than one focused topic, a typical leaf sits
      several `index.md` hops below root, and there are no single-child padding layers.
- [ ] No router contains unique facts or a `Detail` section.
- [ ] Every leaf follows the page contract (two-line header; `Summary`,
      `Detail`, non-empty `Sources`, ≥1 `See Also` neighbor).
- [ ] Every cited file/symbol exists; every relative markdown link resolves; no
      orphans (every page reachable from `AGENTS.md`).
- [ ] Every cited authority is registered in `raw/SOURCES.md` with provenance,
      status, real mtime + `sha256:` hash, and ingest state.
- [ ] Durable terms are centralized in `wiki/glossary.md`; every non-trivial
      claim is sourced or marked `provisional`/`contested`/`stale`.
- [ ] No excluded artifacts (build products, binaries, dumps, token reports) were
      compiled into pages.
- [ ] The first source set is small, useful, and actually compiled into pages;
      future queries start from the router graph, not whole-tree search.
```
