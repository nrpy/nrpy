# Dendro infrastructure normalization plan

## Current state

The active worktree already contains the first normalization pass:

- `CFunction_roles.py`, `block_kernel_helpers.py`,
  `gridfunction_name_decorations.py`, and `bssn_host_adapter.py` are deleted;
- `Ricci_eval.py` exists as a separate generator;
- BSSN and fCCZ4 examples use one parallel-codegen wave;
- `DendroGridFunction` maps NRPy conformal-factor storage to
  `cf_W_or_chi` at the Dendro connector;
- RHS, Ricci, constraint, initial-data, solver-context, state, parameter, and
  CMake generators still contain abstractions or duplicate layers that must be
  removed;
- `NRPy_fCCZ4_GR/` is incomplete generated output, not accepted source.

The current generated fCCZ4 output also exposes a correctness defect that must
be fixed before any structural cleanup is accepted: parallel workers bind
`RbarDD00`--`RbarDD22` to worker-local AUXEVOL slots 0--5, while the merged
state places them in slots 18--23.  The generated Ricci and RHS kernels
therefore agree with each other but access the wrong final fields.  Freeze the
complete field layout before parallel lowering and validate every generated
binding against the final state enumeration.

Preserve valid current work. Replace only remaining nonconforming structure.

## Goal

Refactor `nrpy/infrastructures/Dendro` to follow established ETLegacy and BHaH
generator structure. Generated `NRPy_BSSN_GR` and `NRPy_fCCZ4_GR` remain complete
Dendro-GR applications beside `BSSN_GR`; generator code remains in NRPy.

Running either example emits the complete application, not a kernel library,
standalone mock host, or adapter around `BSSN_GR`. With no output argument, the
example locates the enclosing Dendro-GR checkout by finding the nearest ancestor
containing both root `CMakeLists.txt` and `BSSN_GR/`, then writes
`NRPy_BSSN_GR/` or `NRPy_fCCZ4_GR/` there. An explicit output path remains
available for isolated deterministic-generation checks. Generation fails if no
Dendro-GR root is found; it never silently writes a nested application under
the NRPy repository.

## Required numerical behavior

- Generate finite-difference orders 4, 6, and 8 in one parallel-codegen wave.
- Use centered derivatives. Do not introduce upwinding.
- Use regular/KO pairs `4/2`, `6/4`, and `8/6`; stencil radii remain 2, 3,
  and 4 points.
- Preserve W evolution, SSL, Dendro-style CAHD, Brown's covariant Lambda
  adjustment, constant `eta = 1`, TwoPunctures `alpha = W`, boundary
  treatment, field ordering, and Dendro-GR solver behavior.
- Keep Ricci evaluation separate from RHS evaluation for BSSN and fCCZ4.
  Call Ricci immediately before RHS evaluation without another synchronization
  or block traversal.
- After initialization, every Runge--Kutta stage, and AMR transfer, floor
  `alpha` at `CHI_FLOOR` and W at `sqrt(CHI_FLOOR)` before algebraic
  projection. This is the W-form equivalent of Dendro-BSSN's chi floor and
  prevents interpolation undershoots from entering the next RHS evaluation.
- Do not change `nrpy/equations/general_relativity/**` or `/work/BSSN_GR/**`.

Pin the formulation-specific CAHD expressions explicitly; these are different
formulations, not two spellings of one host coefficient:

```text
BSSN:   W_rhs += (C_CAHD/2) W H dx^2 / (dt (1 + 10 dx^2)), C_CAHD=0.06
fCCZ4:  W_rhs += 2 C_CAHD W H dt,                         C_CAHD=0.15
```

Here `dx` is the current refinement-level spacing and `dt` is the current
time-step size. Register exactly one generated `C_CAHD` parameter per
formulation. Do not assign an absent host member or multiply by a second
host-side coefficient. Validate each formula numerically at every generated
finite-difference order and at two refinement levels.

Runtime dispatch uses the configured Dendro element order, restricted to
`4`, `6`, or `8`. It selects the matching generated kernel, requires padding
of at least `order/2`, and requires KO order `order-2`; any mismatch is a hard
configuration error. Each profile receives an application smoke run, not only
a source-presence check.

## One Python generator per generated numerical kernel

When one Python module owns one generated C++ kernel, these four names match
except for required language prefixes and the `.cpp` suffix:

1. Python file basename;
2. `register_CFunction_<name>()` suffix;
3. registered CFunction name;
4. emitted C++ source basename.

Examples:

```text
rhs_eval.py       register_CFunction_rhs_eval()       rhs_eval       rhs_eval.cpp
Ricci_eval.py     register_CFunction_Ricci_eval()     Ricci_eval     Ricci_eval.cpp
ADM_to_BSSN.py    register_CFunction_ADM_to_BSSN()    ADM_to_BSSN    ADM_to_BSSN.cpp
```

An order-specific kernel uses the same order suffix in CFunction and source
name, for example `rhs_eval_order_6` and `rhs_eval_order_6.cpp`. One module may
register orders 4, 6, and 8 because they are variants of one operation.

Dendro does not require `initial_data.cpp`. Delete catch-all `initial_data.py`.
Connection initialization needs neighboring converted fields, so preserve its
separate pass and exact names:

```text
initial_data_lambdaU.py
register_CFunction_initial_data_lambdaU()
initial_data_lambdaU
initial_data_lambdaU.cpp
```

## Core general-relativity modules

Use this layout:

```text
general_relativity/
  ADM_to_BSSN.py
  BSSN_constraints.py
  BSSN_to_ADM.py
  Ricci_eval.py
  fCCZ4_constraints.py
  rhs_eval.py
  initial_data_lambdaU.py
  enforce_detgbar_equals_detghat_trAzero.py
  floor_the_lapse_and_conformal_factor.py
```

Each file contains one public registration routine, apart from finite-difference
order variants of one operation. Each routine reads linearly:

1. register owned `CodeParameter`s and gridfunctions;
2. construct symbolic expressions;
3. lower expressions to C++;
4. define `desc`, `cfunc_type`, `name`, `params`, and `body`;
5. call `cfc.register_CFunction()`.

Inline one-use expression constructors and `_impl` registration functions.
Do not add expression records, build records, role records, or host adapters.

### Ricci and RHS execution

`Ricci_eval.py` registers only a per-block Ricci kernel. `rhs_eval.py` registers
only the matching per-block RHS kernel. Solver context owns block traversal and
calls both kernels within that traversal:

```text
Ricci_eval(block)
rhs_eval(block)
```

Do not register all-block or flat-block wrappers. Do not synchronize between
these calls.

Use the actual Dendro ABI directly. Do not retain `block_geometry_struct` or
introduce another transport structure. Each order-specific kernel receives the
current `ot::Block`, component-pointer tables, generated parameter object,
stage time where needed, and physical-domain bounds. The kernel derives its
offset, padded dimensions, spacings, and interior loop bounds directly from the
block. The generated solver context owns the unzipped EVOL input/output and
six-component Ricci scratch arrays and passes the correct block offsets.

This fused traversal is valid only after one EVOL halo exchange and physical
boundary fill.  Add generation-time checks that Ricci produces exactly the six
`RbarDD` components and that the RHS takes no derivative of an `RbarDD`
component.  The solver then calls `Ricci_eval(block)` immediately followed by
`rhs_eval(block)` for each block, with no exchange, synchronization, or second
block traversal between them.

`RbarDD00`, `RbarDD01`, `RbarDD02`, `RbarDD11`, `RbarDD12`, and `RbarDD22`
form a dedicated six-component scratch layout, not a suffix of a mixed AUXEVOL
registry. Ricci writes that layout and RHS reads the same layout. No other field
may enter it. This removes the current worker-local/final-registry index
ambiguity.

### Initial data

`ADM_to_BSSN.py` performs pointwise conversion from ADM fields to evolved
fields. `initial_data_lambdaU.py` computes derivative-dependent conformal
connection values after all converted fields are available. Solver context owns
required synchronization between these passes.

The executable initialization sequence is fixed:

```text
TwoPunctures ADM data
  -> pointwise ADM_to_BSSN
  -> zip / halo exchange / unzip
  -> physical-boundary fill
  -> initial_data_lambdaU
  -> algebraic projection
```

`twopunctures.py` uses NRPy's existing BHaH TwoPunctures registrations to emit
application-local TwoPunctures sources and a module-owned Dendro interpolation
driver. The driver exposes `gammaDD`, `KDD`, `betaU`, `BU`, and full
`psi=psi_background+u`; these feed generated `ADM_to_BSSN`. It sets
`alpha=W=(psi_background+u)^(-2)` through TwoPunctures' `W` lapse branch. The
generated CMake manifest includes these local sources.
It explicitly excludes native `/work/BSSN_GR/src/TwoPunctures.cpp`,
`TPID.cpp`, `parameters.cpp`, and `adm2bssn.h`: those write chi and `GT`
through `bssn::VAR` and cannot serve a W/Lambda application. Do not copy or
edit `/work/BSSN_GR`. Validate the connection constraint on at least two MPI
ranks so a block/refinement boundary is exercised.

After every completed time step, interpolate the evolved shift at both
puncture centers and integrate `dx^i/dt=-beta^i`, matching native
Dendro-BSSN. Use the moving centers for constraint excision, black-hole AMR,
and apparent-horizon fallback locations. Checkpoint and restore both centers.

Algebraic projection runs after initial-data construction, after every RK
stage before the next ghost exchange/RHS evaluation, and after every AMR
transfer. Checkpoints contain already-projected state. Restore verifies the
stored algebraic residual and preserves the exact stored state; it neither
silently reprojects nor changes Lambda.

Minkowski and smooth-perturbation qualification builders do not belong in
production initial-data generation. TwoPunctures remains production initial
data. `BSSN_to_ADM.py` reconstructs ADM fields for wave extraction and other
geometric calculations. BHaHAHA instead interpolates the evolved W-BSSN fields
and applies the W-BSSN-to-ADM conversion at each search point. Nonlinear
conversion before interpolation changes the low-resolution horizon equation.

## Complete generated-application ownership

Each generated sibling owns its executable, solver context, namespaces, state
enumerations, parameters, parameter reader, evolution scheduling, AMR transfer
hooks, physical-boundary dispatch, checkpoint/restart metadata, output and
diagnostic scheduling, constraint calls, wave extraction calls, apparent-
horizon calls, and formulation conversions. No generated source includes or
compiles `bssnCtx.cpp`, `bssngr_main.cpp`, `rhs.cpp`, `bssn_constraints.cpp`,
`grUtils.cpp`, or another chi-specific BSSN application implementation.

The generated CMake manifest contains no source from `/work/BSSN_GR` and never
links `bssn_common`. Dendrolib is the shared runtime dependency. All
formulation/application sources are generated module-owned; there is no
preprocessor renaming, adapter, or implicit source glob.

| Runtime service | Generated sibling responsibility |
| --- | --- |
| TwoPunctures | Configure and solve common puncture data; expose ADM fields and full `psi_background+u`; convert through generated `ADM_to_BSSN` |
| Time integration | Own Dendro RK context and call generated Ricci/RHS/projection kernels |
| Ghosts and boundaries | Own zip/exchange/unzip order and generated state-aware physical-boundary dispatch |
| AMR/remesh | Own refinement indicators and transfer using the generated canonical field enumeration |
| Checkpoint/restart | Write and validate formulation name, field names/order/count, parameters, iteration, and time before restore |
| Constraints | Call generated BSSN or fCCZ4 constraint kernels; never Dendro-BSSN's chi-specific constraint kernel |
| Geometry consumers | Convert generated W/formulation state through `BSSN_to_ADM` for waveform extraction and ADM diagnostics; BHaHAHA converts after interpolation |
| Output | Own VTU, scalar diagnostic, timing, and horizon output schedules and generated field-name tables |

Emitter-to-output ownership is explicit:

| Python emitter | Generated source responsibility |
| --- | --- |
| `general_relativity/twopunctures.py` | `twopunctures.cpp` plus application-local BHaH TwoPunctures sources and ADM interpolation |
| `general_relativity/physical_boundary.py` | `physical_boundary.cpp` |
| `general_relativity/diagnostics.py` | `diagnostics.cpp` |
| `general_relativity/apparent_horizon.py` | `apparent_horizon.cpp`, including the module-owned Dendro/BHaHAHA bridge without `bssnCtx` |
| `general_relativity/psi4_eval.py` | `psi4_eval_order_4.cpp`, `psi4_eval_order_6.cpp`, and `psi4_eval_order_8.cpp` |
| `general_relativity/gravitational_waves.py` | `gravitational_waves.cpp` |
| `general_relativity/adm_quantities_surface_data.py` | `adm_quantities_surface_data_order_4.cpp`, `adm_quantities_surface_data_order_6.cpp`, and `adm_quantities_surface_data_order_8.cpp` |
| `general_relativity/adm_quantities.py` | `adm_quantities.cpp` |
| `checkpoint.py` | `checkpoint.cpp` |
| `solver_context.py` | generated context, RK staging, AMR/remesh/transfer, and service scheduling |
| `main_cpp.py` | generated executable entry point |

`CMakeLists.py` writes an explicit source manifest containing the executable,
context, parameter/state support, local TwoPunctures implementation, all
runtime services above, projection/conversion kernels, and order-4/6/8 Ricci,
RHS, and constraint kernels. No source glob and no BSSN_GR source is permitted.

The generated parameter file contains every parameter needed by these paths.
The existing low-resolution TwoPunctures parfile must run without a translation
script or hand edit. Checkpoint restore must reject the wrong formulation or
field layout rather than reinterpret storage.

### Canonical component order

Never derive runtime storage order from alphabetical registry sorting. Define
one explicit, formulation-specific sequence before parallel work. For BSSN,
preserve Dendro's established semantic order while replacing its chi slot by W
and its non-tensor connection slot by the covariant NRPy variable:

```text
alpha, cf_W_or_chi, trK,
lambdaU0, lambdaU1, lambdaU2,
vetU0, vetU1, vetU2,
betU0, betU1, betU2,
hDD00, hDD01, hDD02, hDD11, hDD12, hDD22,
aDD00, aDD01, aDD02, aDD11, aDD12, aDD22
```

The fCCZ4 sequence is exactly the same 24-component BSSN-compatible base,
followed by `Theta_fCCZ4` at component 24. Define explicit ordered lists for
EVOL, six-component Ricci scratch,
constraints/diagnostics, and any auxiliary state. Generate enum values, name
tables, pointer bindings, zip/unzip, transfer, output selection, checkpoint
metadata, and diagnostic maps from those same lists. Generation fails if a
registered component is absent, duplicated, reordered, or silently added.

## Remove unsupported abstraction layers

Delete these modules or keep existing deletion:

- `CFunction_roles.py`;
- `block_kernel_helpers.py`;
- `gridfunction_name_decorations.py`;
- `gridfunction_bindings.py`;
- `generated_file_banner.py`;
- `header_guards.py`;
- `general_relativity/generation_parameters.py`;
- `general_relativity/bssn_host_adapter.py`;
- `standalone_host/`.

Do not replace them with renamed helpers.

### Gridfunction access

Keep `DendroGridFunction`; finite-difference code generation needs a Dendro
storage expression. Keep only substantial access methods. Remove pointer-name
helpers and short formatting helpers.

Use `DendroGridFunction.dendro_name` directly. At the Dendro connector,
canonical NRPy gridfunction `cf` maps to `cf_W_or_chi`. Do not change equation
module gridfunction names.

Each numerical registrar emits its short pointer-binding loop directly. Bind
registered fields once in registry order. Do not scan expressions to infer a
read set. Construct RHS and diagnostic maps with registered gridfunction names
as keys. Never recover names through string replacement or regular expressions.

### Loop construction

Keep `simple_loop.py`; Dendro needs padded block bounds, x-fastest indexing, and
no inner OpenMP region. Inline the serial-parallelization check into
`simple_loop()`. Delete `block_loop()` after wrapper kernels are removed.

## Generated application files

Keep one substantial emitter for each generated file. Small banner,
header-guard, path, and naming helpers are not separate modules.

```text
Dendro/
  CMakeLists.py
  CodeParameters.py
  Dendro_defines_h.py
  constants_h.py
  main_cpp.py
  param_toml.py
  simple_loop.py
  solver_context.py
  state_h.py
  types_h.py
```

Retain `block_geometry.h` only until the direct-`ot::Block` kernel ABI is
emitted and all consumers compile; then delete it. It is not part of the final
layout.

- `CMakeLists.py` corresponds to ETLegacy `make_code_defn.py`.
- `param_toml.py` corresponds to ETLegacy `param_ccl.py`.
- `state_h.py` corresponds to ETLegacy `interface_ccl.py`.
- `solver_context.py` replaces Cactus scheduling and MoL registration.
- `types_h.py`, `constants_h.py`, and `Dendro_defines_h.py` provide definitions
  supplied by Cactus in ETLegacy.

Merge duplicate root and `general_relativity` `main_cpp.py` and
`solver_context.py` emitters. Emit exact identifiers directly. Remove placeholder
replacement tables and blanket `str.replace()` passes.

`main_cpp.py` emits the complete Dendro entry point. `solver_context.py` emits
the complete runtime context described above, including initialization,
evolution, remesh, checkpoint, diagnostics, wave extraction, and horizon
dispatch. These are substantial application emitters, not host adapters.

Exclude generated standalone self-test programs from production generation.
Retain both existing `self_tests_cpp.py` files byte-for-byte unchanged. They
are historical test emitters, not imported by either production example; do
not execute them after their private legacy dependencies are removed. Use
actual generated Dendro applications for compilation and runtime checks.
Preserve numerical coverage with temporary, uncommitted old/new kernel
comparisons for BSSN and fCCZ4 at orders 4, 6, and 8. Compare interior stencil
values, conformal-factor algebra, SSL/CAHD terms, separate-Ricci coupling, and
block-boundary values. Do not add a replacement test file or test case without
the user's express permission.

## Parallel code generation

Examples remain thin orchestration modules:

1. set formulation and finite-difference parameters;
2. enqueue expensive registration calls for orders 4, 6, and 8;
3. call `pcg.do_parallel_codegen()` once;
4. perform inexpensive file emission.

All symbolic expression construction and lowering occur in workers. Parent code
must not rebuild RHS, Ricci, or constraint expressions after worker completion.
Use the deterministic CSE sorting algorithm already used by Psi4 when normal
sorting is too slow.

Parallel registry merging must be deterministic. Sort task and registry keys
before merging. Build all merged registries in temporary dictionaries, reject
unequal duplicate definitions, and commit all registries only after every
merge succeeds. Every worker saves and restores global finite-difference order
with `try/finally`.

Before workers start, register and freeze the complete EVOL, AUXEVOL, DIAG, and
AUX component order used by every task.  Workers lower array accesses against
that shared order.  After merging, compare every pointer binding in each kernel
with the final generated state enumeration; generation fails on any mismatch.

The parent snapshots every global registry before starting workers. Each worker
starts from that snapshot and returns definitions without mutating the parent.
The parent sorts task results and registry keys, merges into temporary
dictionaries, verifies that duplicate keys have identical complete
definitions, validates all state bindings against the frozen canonical lists,
then commits all registries in one step. Any worker failure, unequal duplicate,
or layout mismatch leaves every parent registry byte-for-byte unchanged.

## Exact file boundary

This table supersedes the older boundary in `plan_nrpy_refactor.md` where later
user decisions required inlining, direct Dendro names, complete sibling
applications, and retention of existing tests.

Modify only:

- `/work/CMakeLists.txt`;
- `/work/nrpy/.github/single_file_static_analysis.sh`;
- `/work/nrpy/nrpy/examples/dendro_bssn.py`;
- `/work/nrpy/nrpy/examples/dendro_fccz4.py`;
- `/work/nrpy/nrpy/grid.py`;
- `/work/nrpy/nrpy/helpers/parallel_codegen.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/CodeParameters.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/Dendro_defines_h.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/__init__.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/constants_h.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/main_cpp.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/simple_loop.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/solver_context.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/state_h.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/types_h.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/__init__.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/enforce_detgbar_equals_detghat_trAzero.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`;
- `/work/nrpy/wiki/architecture/python-coding-style.md`;
- `/work/nrpy/wiki/catalog.md`;
- `/work/nrpy/wiki/contradictions.md`;
- `/work/nrpy/wiki/source-map.md`;
- `/work/nrpy/wiki/infrastructures/dendro/index.md`;
- `/work/nrpy/wiki/infrastructures/dendro/project-assembly-and-emitters.md`;
- `/work/nrpy/wiki/infrastructures/dendro/bssn-application-wiring.md`;
- `/work/nrpy/wiki/infrastructures/dendro/fccz4-application-wiring.md`;
- `/work/nrpy/wiki/infrastructures/dendro/finite-difference-profiles-and-dendro-conformance.md`;
- `/work/nrpy/wiki/infrastructures/dendro/gridfunctions-naming-and-loops.md`.
- `/work/nrpy/wiki/infrastructures/infrastructure-code-style.md`;

Create only:

- `/work/nrpy/nrpy/infrastructures/Dendro/CMakeLists.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/param_toml.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/checkpoint.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/ADM_to_BSSN.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/BSSN_to_ADM.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/initial_data_lambdaU.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/twopunctures.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/adm_quantities.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/adm_quantities_surface_data.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/physical_boundary.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/diagnostics.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/apparent_horizon.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/gravitational_waves.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/psi4_eval.py`;
- `/work/nrpy/wiki/infrastructures/dendro/grid-amr-and-time-stepping.md`;
- `/work/nrpy/wiki/infrastructures/dendro/validation-standalone-host-and-deferral-gates.md`;
- generated `/work/NRPy_BSSN_GR/`;
- generated `/work/NRPy_fCCZ4_GR/`.

Delete after all consumers are migrated:

- `/work/nrpy/nrpy/infrastructures/Dendro/CFunction_roles.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/block_kernel_helpers.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/block_geometry.h`;
- `/work/nrpy/nrpy/infrastructures/Dendro/cmake_helpers.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/generated_file_banner.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/gridfunction_bindings.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/gridfunction_name_decorations.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/header_guards.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/parfile.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/standalone_host/__init__.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/standalone_host/dendro_standalone_host.h`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/bssn_host_adapter.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/generation_parameters.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/initial_data.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/main_cpp.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/solver_context.py`;
- incidental generated `/work/nrpy/NRPy_fCCZ4_GR/`.

Preserve unchanged:

- `/work/nrpy/nrpy/equations/general_relativity/**`;
- `/work/BSSN_GR/**`;
- `/work/nrpy/.github/workflows/main.yml`;
- `/work/nrpy/nrpy/infrastructures/Dendro/self_tests_cpp.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py`;
- `/work/nrpy/nrpy/infrastructures/Dendro/tests_infra/**`;
- every other test file and test case;
- unrelated user files and run output.

Do not stage or commit generated applications into the nested NRPy repository.
The examples write them directly as siblings of `/work/BSSN_GR`.  The root
Dendro-GR `CMakeLists.txt` owns only conditional `add_subdirectory()` entries
for those two siblings; it does not absorb generated source inventories.
Remove the incidental nested `/work/nrpy/NRPy_fCCZ4_GR/` only after confirming
it is untracked generated output and the two sibling applications have been
generated successfully.

## Tests and validation

Do not add or extend any test file, test case, doctest prompt, trusted output,
or generated oracle without express user permission.  Existing prompts may
receive only meaning-preserving identifier/path updates required by the
refactor.  Do not add assert-only doctests.

Create disposable full-tree `HEAD` and candidate checkouts and place all
caches/output outside both trees. First run `black .` from the isolated
candidate checkout, inspect its diff, reject any formatting change outside the
exact boundary above, and transfer only accepted in-boundary formatting to the
working tree/candidate checkout. Then, for every changed handwritten Python
file present in the candidate checkout, run NRPy's
`.github/single_file_static_analysis.sh <file>` once there. New Python files
must score exactly `10.00/10.00`. For
every changed tracked handwritten Python file that exists in both trees, run
the same wrapper, configuration, arguments, and environment once in each
`HEAD` and candidate checkout; compare exit status and diagnostics per file,
with no candidate regression. Check new files only in the candidate; deleted
files require no wrapper run. Before every wrapper call, resolve and validate
its argument as a
repository-relative path beneath that checkout: reject an absolute path,
`..`, symlink escape, missing file, or path outside the exact boundary above.
Inspect `git status` after all checks. Then perform:

1. an isolated, uncommitted merge-failure experiment that snapshots every
   registry, injects unequal duplicate definitions, requires `ValueError`, and
   proves no registry changed;
2. import and registration checks for BSSN and fCCZ4;
3. two clean generations into separate temporary directories with identical
   arguments and environment, followed by byte-for-byte tree comparison;
4. confirm all 4th-, 6th-, and 8th-order kernels exist;
5. confirm Python module, registrar, CFunction, and `.cpp` names correspond;
6. confirm generated pointer indices equal the final state enumerations;
7. inspect each generated CMake source list: no `bssn_common`, BSSN application
   context, BSSN RHS, BSSN constraints, or chi-specific utility is present;
8. configure from `/work` into `/work/build-nrpy`, build both generated targets
   with `cmake --build /work/build-nrpy --parallel 96`, and also configure each
   generated application independently;
9. run the exact BSSN production gate
   `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 mpirun -n 2 -x OMP_NUM_THREADS -x OPENBLAS_NUM_THREADS -x MKL_NUM_THREADS --mca pml ob1 --mca btl self,vader,tcp /work/build-nrpy/NRPy_BSSN_GR/nrpyBssnSolver /work/q1.par.lowres.toml`
   through
   initialization and the first completed diagnostic dump; require finite
   nontrivial data, `alpha=W`, successful Lambda initialization/projection,
   and no connection-constraint discontinuity at a block boundary;
10. force one remesh, write a checkpoint, restart it, and compare field names,
    values, iteration, time, constraints, and diagnostic scheduling across the
    restart;
11. exercise constraints, wave extraction, and apparent-horizon scheduling in
    the BSSN smoke run; require their output files and no chi/Lambda/DGamma
    reinterpretation;
12. run the mandatory fCCZ4 production gate
   `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 mpirun -n 2 -x OMP_NUM_THREADS -x OPENBLAS_NUM_THREADS -x MKL_NUM_THREADS --mca pml ob1 --mca btl self,vader,tcp /work/build-nrpy/NRPy_fCCZ4_GR/nrpyFccz4Solver /work/q1.par.lowres.toml`
   through initialization and the first completed diagnostic dump; require
   finite nontrivial fields, `alpha=W`, successful Lambda initialization and
   algebraic projection, and no connection-constraint discontinuity across a
   block/refinement boundary; no Minkowski or unsupported-path substitution is
   acceptable;
13. run `python tools/kb_lint.py`, review changed prose for direct
    computational-physics language, reconcile the catalog/source-map entries,
    and run `git diff --check` in both repositories.

For application gates 9 and 12, make temporary parameter copies outside the
repository selecting element orders 4, 6, and 8. Each must complete its first
diagnostic dump with the corresponding FD/KO profile. Delete those temporary
copies after recording results.

Equation validation does not treat the current generated fCCZ4 tree as an
oracle. Compare against protected equation-level expressions through an
explicit canonical field mapping and an independent temporary direct numerical
evaluation; use existing trusted fixtures where available. No new test file,
test case, prompt, or oracle is added.

Do not modify `.github/workflows/main.yml`.

## Trialectic delta from the current state

The implementation begins from the partially normalized tree, not from the
pre-refactor baseline.  The already-deleted role, block-helper, naming-helper,
and BSSN-host-adapter layers stay deleted.  Preserve the existing direct
``DendroGridFunction.dendro_name`` connector mapping, the separate
``Ricci_eval.py``, W-only SSL/CAHD expressions, deterministic
``cse_sorting="none"``, and the atomic-registry work already present.  Replace
the remaining wrapper kernels, catch-all modules, duplicate GR application
emitters, standalone host, and incomplete generated fCCZ4 tree.

The Trialectic science seat fixed the numerical contracts: the six Ricci
components are non-evolved scratch, no RHS differentiates them, BSSN receives
Brown's term exactly once, fCCZ4 retains its own pre-Brown Lambda RHS, and the
two distinct CAHD formulas above use the live block spacing and time step.  The
integration seat fixed the application boundary: generated modules may use
``dendro5`` and its public PR-13 BHaHAHA interface, but no source, include, or
target from ``BSSN_GR``.  The standards seat found that deleting renamed
modules also requires meaning-preserving KB-link updates in
``wiki/glossary.md``, ``wiki/syntheses/generated-backend-comparison.md``, and
``wiki/validation/generated-project-ci.md``; those three documentation files
are added to the modification boundary.

``ADM_to_BSSN.py`` emits ``ADM_to_BSSN_order_N.cpp`` directly.  Dendro imposes
no ``initial_data.cpp`` name, so no compatibility wrapper or catch-all
``initial_data.py`` remains.  The derivative-dependent connection pass stays
separate as ``initial_data_lambdaU.py``.

The existing q1 file's explicit eta value remains a runtime override; the
generated default is eta=1.  Thus an unchanged parameter file retains normal
TOML precedence without weakening the requested default.

## Implementation order from current state

1. Freeze the exact file boundary above. Preserve unrelated user changes.
2. Pre-register the canonical component lists, separate the six Ricci scratch
   components, and make parallel registry merging atomic; prove final component
   indices before further experiments.
3. Normalize core kernel modules and enforce one-to-one names.
4. Remove gridfunction-binding and expression-builder layers.
5. Merge duplicate solver-context and main emitters; replace the standalone
   mock host and BSSN adapter with the complete module-owned Dendro application;
   implement the explicit TwoPunctures-to-Lambda initialization sequence and
   all runtime-service ownership above.
6. Normalize CMake, parameter, state, and header emitters.
7. Update both examples to emit sibling `NRPy_BSSN_GR` and
   `NRPy_fCCZ4_GR` applications.
8. Remove obsolete imports, files, adapters, wrappers, and stale documentation.
9. Run single-file static analysis after each Python edit.
10. Generate twice, compare outputs, build with 96 jobs, and run focused smoke
   tests.
11. Submit final diff and validation evidence to independent scientific,
    integration/simplification, and standards/release reviewers. Resolve every
    substantive finding before delivery.
