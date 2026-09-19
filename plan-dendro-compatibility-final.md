# Final plan: Dendro-compatible NRPy BSSN and fCCZ4 modules

## Purpose

Modify NRPy so its Dendro examples emit host-linkable BSSN and fCCZ4 modules with the numerical profiles and data ownership expected by Dendro-GR:

- centered finite differences with shift upwinding disabled;
- finite-difference orders 4, 6, and 8, with order 6 as the default;
- two, three, and four ghost points, respectively;
- Kreiss-Oliger (KO) finite-difference order two below the regular finite-difference order;
- no fixed block dimensions or physical-cell counts;
- Dendro-owned mesh, blocks, padding, offsets, coordinates, spacing, storage, communication, refinement, and time integration; and
- the same generated-module structure for BSSN and fCCZ4.

The public [Dendro-GR repository](https://github.com/paralab/Dendro-GR) is the read-only conformance reference. Generated qualification projects must be written outside its checkout.

## Status and qualification coverage

Sections 1 through 8 record the implemented baseline. The configured
qualification now:

- exercises real-host BSSN and fCCZ4 FD4, FD6, and FD8 modules with KO enabled,
  plus both formulations at FD6 with KO disabled;
- runs sanitizer-enabled block-offset and flat-storage checks;
- injects host element-order and TOML profile mismatches and checks their
  diagnostics;
- prints transient FD6 block-RHS time, allocation count, generated object and
  source size, and compile time outside the KB; and
- compares independent KO coefficients and signs with the public Dendro-GR
  `ko_deriv21`, `ko_deriv42`, and radius-four `ko_deriv64` functions.

These are configured checks, not a stored CI run result. Real-host runtime
measurements remain diagnostic values without permanent pass thresholds.

The twelve-profile standalone matrix, independent nonflat expression check,
manufactured centered-derivative convergence check, centered-stencil reach
check, and structural absence of directional advection are already present.

## Resolved design decisions

### Generate one numerical profile per project

Each invocation emits one formulation and one finite-difference order. The generator supports all three profiles, but does not place three complete RHS kernels and a runtime dispatcher in one module.

This is the minimal Dendro design. Dendro-GR selects one derivative family for a build and runs a mesh with one element order and padding profile. Its block callback does not provide the Einstein Toolkit runtime `fd_order` selector used by BaikalVacuum. A Baikal-style combined project would add three large kernels, multiple metadata sets, runtime dispatch, and mesh/profile mismatch paths without enabling mixed-order Dendro blocks.

Generate separate projects or libraries for orders 4, 6, and 8 when one enclosing application needs all three build choices. Reconsider a combined module only after a concrete Dendro runtime interface requires order switching without module selection.

### Use an unambiguous KO-order convention

Use `ko_fd_order` for the order supplied to NRPy's KO construction, before the existing `dKOD` increment by two. The Dendro setting is always:

```text
ko_fd_order = fd_order - 2
```

| Centered `fd_order` | `ko_fd_order` | Effective KO difference order | KO reach | Required padding |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 2 | 4 | 2 | 2 |
| 6 | 4 | 6 | 3 | 3 |
| 8 | 6 | 8 | 4 | 4 |

The effective KO difference order is `ko_fd_order + 2`, matching the current NRPy `dKOD` convention. This naming prevents “KO order” from ambiguously meaning either the caller's accuracy setting or the actual even difference. It also makes the generated radius equal the centered-derivative radius for every supported profile.

Dendro-GR currently selects `ko_deriv21`, `ko_deriv42`, and `ko_pw4_deriv42` for its common 4/6/8 paths, but also contains the radius-four `ko_deriv64` family. The NRPy FD8 profile follows the required `ko_fd_order=6` rule and must verify its coefficients and signs against the radius-four family rather than copying the current FD8 selector.

### Center the shift-advection terms inside the Dendro builder

Do not change the continuum BSSN or fCCZ4 equation modules. After a formulation assembles its RHS and optional KO terms, the Dendro builder copies the expression dictionaries and changes derivative symbols in shift-advection terms from `_dupD<i>` or `_ddnD<i>` to `_dD<i>`. It then lowers the copied centered expressions without an upwind control vector.

Use parsed derivative symbols and SymPy substitution. Do not use unrestricted string replacement or mutate cached shared equations. Assert after normalization and after derivative discovery that no upwind, downwind, full-upwind, or full-downwind operator remains. The physical shift and Gamma-driver equations remain present; only their spatial derivative stencil changes.

## Implemented baseline requirements

### 1. Add a narrowly scoped KO finite-difference input

Modify:

- `nrpy/finite_difference.py`
- `nrpy/c_codegen.py`

Add optional keyword-only `ko_fd_order` inputs to `compute_fdcoeffs_fdstencl`, `stencil_reach_per_axis`, `CCodeGen`, and `c_codegen`.

Required behavior:

- `None` preserves every existing caller's behavior: `dKOD` starts from the regular `fd_order` and internally adds two.
- An explicit value replaces the starting order only for `dKOD`; the existing increment by two still produces the effective KO difference.
- Require a positive even value.
- Non-KO centered, mixed, upwind, and full-upwind derivatives are unchanged.
- Thread the value through all coefficient construction and reconstruction routes, including stored derivatives, prototype expressions, finite-difference helper functions, and helper naming.
- Include the KO setting in any helper identity or cache key. Do not reuse an upwind/KO helper merely because the regular `fd_order` matches.
- Disable the existing upwind/downwind-to-KO algebraic reuse whenever an explicit `ko_fd_order` differs from regular `fd_order`. That identity assumes the default `dKOD(fd_order + 2)` stencil. If the general coefficient path cannot safely emit an explicit KO order together with an upwind control vector, reject that combination; never apply the identity with a mismatched order.
- Make `stencil_reach_per_axis` use the same value as emitted coefficient generation.

If correct helper naming would require a broad registry redesign, the minimal acceptable first implementation may reject explicit `ko_fd_order` with `enable_fd_functions=True`, because the Dendro examples use inline finite differences. It must fail explicitly rather than silently reuse a helper with different coefficients.

Tests in the defining modules must prove the unchanged legacy default, exact effective orders and radii for Dendro's 2/4/6 KO settings, invalid-input rejection, inline/reconstructed coefficient agreement, equality of emitted and reported reach, and correct rejection or non-reuse for an explicit KO order combined with an upwind control vector.

### 2. Define the Dendro numerical profiles once

Modify:

- `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`
- `nrpy/infrastructures/Dendro/block_kernel_helpers.py`
- `nrpy/infrastructures/Dendro/CFunction_roles.py` only if existing calling code demonstrably needs another scalar

Use a small immutable mapping, not a general profile framework:

```text
4 -> ko_fd_order 2, padding 2
6 -> ko_fd_order 4, padding 3
8 -> ko_fd_order 6, padding 4
```

Set the builder default to 6 and reject every other regular order. Pass the mapped KO setting to code generation and stencil-reach calculation.

Continue computing padding from the derivative operators actually present. Check the computed result against the profile rather than replacing reach analysis with a hard-coded return. Require:

- widest emitted reach equals the profile padding;
- `dKOD` exists exactly when KO is enabled;
- KO reach never exceeds centered reach;
- centered/advection normalization leaves no directional derivative family; and
- recorded order, KO setting, KO switch, and padding agree.

### 3. Normalize both formulations to centered advection

Modify `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`.

Normalize the BSSN and fCCZ4 expression dictionaries inline in their shared Dendro RHS build flow before free-symbol discovery, derivative discovery, C common-subexpression elimination, code generation, and storage in the RHS build result. This is one symbolic setup step with one call site, so it does not introduce a private single-use helper. Record no Dendro upwind metadata and pass the no-upwinding sentinel expected by `c_codegen`.

In the same slice, remove the Dendro role APIs and fields for `set_upwind_control_fields` and `upwind_control_fields`, the generated state-header control arrays, `UPWIND_ALG`, and their dependent doctests. Convert all internal oracle callers before removing the fields; do not preserve an empty tuple as a public compatibility interface. Generic NRPy upwind support remains unchanged.

Tests must cover representative scalar, vector, and tensor shift-advection terms; preservation of target/evolved-field mappings; absence of directional operators; and presence of the centered `beta^i partial_i u` contribution. Shared equation files remain unchanged.

### 4. Update the two Dendro example entry points

Modify:

- `nrpy/examples/dendro_bssn.py`
- `nrpy/examples/dendro_fccz4.py`

For both:

- change `--fd-order` choices to `(4, 6, 8)`;
- set the default to 6;
- derive `ko_fd_order=fd_order-2` and padding from the checked profile;
- pass `ko_fd_order` only to RHS coefficient generation and reach calculation, pass `fd_order` to RHS and any derivative-emitting constraint or diagnostic code, pass explicit scalar profile values to constants, parameter files, host checks, and tests, and pass neither order to pointwise algebraic projection;
- report the regular order, KO finite-difference order, effective KO difference order when enabled, and required padding; and
- remove order-2, order-10, and five-ghost-zone guidance.

Preserve the current switch defaults unless a separate scientific decision changes them: BSSN KO off, fCCZ4 KO on. Qualify both KO states for every order.

Each invocation still writes one project. Regeneration replaces the selected generated project. Callers retaining several profiles use distinct output directories; NRPy does not add profile-versioning machinery.

Pass `ko_fd_order` only to RHS coefficient generation and derivative-reach calculation. Pass ordinary `fd_order` only to code that actually emits derivatives. Give constants, parameter files, host checks, and corresponding tests explicit scalar profile values. Do not thread a finite-difference profile through pointwise algebraic projection or constraint code that emits no derivatives, and do not extend `CFunction_roles` without existing calling code.

### 5. Emit and validate the compiled profile

Modify:

- `nrpy/infrastructures/Dendro/constants_h.py`
- `nrpy/infrastructures/Dendro/parfile.py`
- `nrpy/infrastructures/Dendro/main_cpp.py`
- corresponding doctests

Emit unambiguous constants for:

- `FD_ORDER`;
- `KO_FD_ORDER`;
- `KO_EFFECTIVE_DIFFERENCE_ORDER`;
- `KO_ENABLED`; and
- `REQUIRED_PADDING`.

Retain profile values unconditionally: `KO_FD_ORDER=FD_ORDER-2` and `KO_EFFECTIVE_DIFFERENCE_ORDER=FD_ORDER`, whether KO is enabled or disabled. `KO_ENABLED` alone states whether the operator is present. Use this convention consistently in headers, parameter files, startup checks, and diagnostics.

The parameter file describes and validates the compiled kernel; it does not select another order at runtime. A mismatch in any profile value must fail before RHS evaluation and name the expected and observed values.

### 6. Keep all production geometry under Dendro ownership

Modify only where current code violates this rule:

- `nrpy/infrastructures/Dendro/solver_context.py`
- `nrpy/infrastructures/Dendro/general_relativity/solver_context.py`
- `nrpy/infrastructures/Dendro/main_cpp.py`
- affected runtime qualification tests

Retain `block_geometry()` as the single conversion from each `ot::Block`. Read every block's allocation sizes, one-dimensional padding, component offset, spacing, node coordinates, and physical padded origin from the host. Keep boundary flags in the host adapter for physical-exterior treatment; add them to `block_geometry_struct` only if a generated numerical kernel demonstrably consumes them.

At context/startup validation:

- require host element order to equal compiled `FD_ORDER`;
- require each block's padding to equal `REQUIRED_PADDING`;
- validate `nx`, `ny`, and `nz` independently and require each to exceed twice its padding;
- validate component intervals against the host allocation; and
- require finite positive spacing.

Never infer element order as `2 * REQUIRED_PADDING` in production code, choose an element order for the host, derive block dimensions from padding, assume equal axes, or encode seven physical cells, thirteen padded points, or any other fixed extent.

Keep the host-owned flow:

```text
Dendro vector
  -> Dendro exchange and unzip
  -> physical-exterior boundary fill
  -> host local-block traversal
  -> one generated block RHS call with host geometry and offset
  -> Dendro zip
```

Generated kernels do not allocate a mesh, exchange MPI data, partition/refine the octree, own block storage, or choose timesteps. Demonstration programs and numerical tests may use explicit test dimensions, but must include noncubic dimensions without treating them as production block sizes.

### 7. Complete the proper Dendro-GR module boundary

Modify:

- `nrpy/infrastructures/Dendro/solver_context.py`
- `nrpy/infrastructures/Dendro/general_relativity/solver_context.py`
- `nrpy/infrastructures/Dendro/cmake_helpers.py`
- both examples

Expose only the whole-vector and Berger-Oliger block callback compatibility signatures declared by the configured `ts::Ctx` interface in `BSSN_GR/include/bssnCtx.h`; do not invent overloads. The corresponding reference implementations are currently disabled in `BSSN_GR/src/bssnCtx.cpp`, so test the generated callbacks directly and claim scheduler qualification only if the selected Dendrolib actually invokes them. Reuse one numerical block kernel and the existing flat-block adapter.

- Whole-vector `rhs` performs one host unzip, evaluates the host's local blocks, and performs one zip.
- `rhs_blkwise` evaluates only host-selected block identifiers and neither unzips nor zips.
- `rhs_blk` validates the block identifier and field count, derives geometry from the referenced host block, treats its supplied component-major slab as block-local, and passes `component_offset=0` exactly once.
- `pre_stage_blk`, `post_stage_blk`, and `pre_timestep_blk` are explicit no-ops. `post_timestep_blk` applies the same generated algebraic projection as whole-vector `post_timestep`, using only block-local storage with component offset zero. Physical-exterior filling remains host-owned and is not duplicated in these hooks.
- Host-supplied block times and schedules remain authoritative.

Add whole-vector versus per-block `post_timestep` projection equivalence tests for BSSN and fCCZ4. The block tests must prove selected-block-only writes and preserve sentinels outside the selected block.

Emit production targets `nrpy_bssn_dendro` and `nrpy_fccz4_dendro` without aliases. Their normal production link interface is `PUBLIC dendro5`, which carries Dendrolib's public MPI, OpenMP, and configuration dependencies. Fail configuration with a named diagnostic if `dendro5` is unavailable. Neither production target links `bssn_common` or `toml11::toml11`. Qualification drivers may link either only when their initial-data or parameter-parsing code actually uses it; otherwise use generated NRPy initial data and parameters.

Make the generated CMake module safe under `add_subdirectory`:

- use `CMAKE_CURRENT_LIST_DIR` instead of a parent `${SRC}` variable;
- keep include paths, compile definitions, and linked libraries target-scoped;
- avoid unconditional top-level `project()` and `enable_testing()` effects;
- place `nrpy_bssn_dendro_qualify` and `nrpy_fccz4_dendro_qualify` executables behind `NRPY_DENDRO_BUILD_DRIVERS`;
- place generated tests behind `NRPY_DENDRO_BUILD_TESTS`;
- default both options on only when the generated project is top-level and off when included as a subproject; and
- verify both generated formulation libraries can coexist in one parent build.

Enable real-host generation for both formulations after the BSSN collision is removed. The production result is a linkable library/context, not an executable that creates a mesh. Keep mesh construction and initial-data adapters in qualification drivers. Preserve CPU-only scope. Use a qualification wrapper outside read-only `BSSN_GR`.

### 8. Replace obsolete upwind tests and extend scientific checks

Modify:

- `nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py`
- `nrpy/infrastructures/Dendro/self_tests_cpp.py`
- `nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp`
- `nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp`
- `nrpy/infrastructures/Dendro/tests_infra/README.md`

Remove tests that require three upwind controls or choose stencils from shift signs. Centered-advection checks perturb both sides, compare with independently computed centered coefficients, and prove that points beyond declared reach do not contribute. Absence of directional derivative operators and upwind-control state establishes shift-sign independence structurally; do not add a duplicate numerical sign-selection test.

Generalize the independent nonflat RHS and constraint oracle beyond order 4. It must not call `compute_fdcoeffs_fdstencl` or reuse generated coefficient tables. Construct radius-2/3/4 rational KO coefficients from the closed binomial formula or fixed rational tables, then compare them with both NRPy output and Dendro-GR's `ko_deriv21`, `ko_deriv42`, and `ko_deriv64` families. Use fields that make centered and KO contributions nonzero. Derive test extents and checked points from padding, but vary their interior sizes and axes so test dimensions cannot become production assumptions.

Add equivalence and bounds checks for:

- whole-vector RHS versus selected-block evaluation;
- pointer-array block kernels versus flat `rhs_blk`;
- nonzero whole-mesh offsets versus zero block-local offsets;
- a noncubic standalone block and multiple actual host blocks without assuming variable host allocation dimensions;
- selected-block-only writes and sentinel regions; and
- early rejection of order, padding, field-count, offset, and block-ID mismatches.

## Qualification matrix

Generate all twelve standalone profiles in fresh disposable directories:

| Formulation | Regular order | KO state |
| --- | --- | --- |
| BSSN | 4, 6, 8 | off and on |
| fCCZ4 | 4, 6, 8 | off and on |

For every profile:

1. Verify generated metadata and coefficient-derived padding.
2. Verify absence of upwind/downwind operators and controls.
3. For KO on, verify exact effective difference order, radius, coefficients, and a nonzero KO response.
4. For KO off, verify absence of `dKOD`.
5. Build with warnings enabled and run generated tests.
6. Run the existing independent nonflat RHS/constraint check and manufactured centered-derivative convergence check.
7. Run sanitizer-enabled block-offset and flat-storage checks where supported.

For real-host qualification, build both formulations against the read-only Dendro-GR source as separate, collision-free subprojects. Exercise host element orders 4, 6, and 8; KO on for all six formulation/order pairs; and KO off at least for both default-order modules, with the full KO-off matrix remaining mandatory in standalone tests. Run whole-vector and block callbacks on one and two MPI ranks where available. Query actual block geometry; do not assert fixed dimensions.

For the default FD6 release profile, record transient qualification measurements for block-RHS time, allocation count, generated object/code size, and compile time. Require `rhs_blk` to allocate no mesh-sized storage and perform no unzip, zip, or MPI operation. These measurements detect gross regressions; they do not establish permanent thresholds and must not be stored as KB snapshots.

If the available Dendrolib exposes callback signatures but does not actively schedule local Berger-Oliger stepping, directly compile and call the exact callbacks and report local-time-stepping scheduling as unqualified. Do not simulate a successful local-time-stepping result.

## Qualification limits

Current tests do not qualify general physical boundaries, AMR remeshing and
state transfer, Dendrolib NUTS scheduling, checkpoint/restart, output selection,
GPU execution, threaded kernels, long-time or nonlinear evolution, or broad
physics results. They also do not establish exact equality with every
Dendro-GR derivative implementation. Keep these limits until separate tests
exercise each named behavior against the intended Dendro application.

## Documentation and delivery checks

Update these directly affected pages:

- `wiki/core/finite-difference.md`;
- `wiki/infrastructures/dendro/bssn-application-wiring.md`;
- `wiki/infrastructures/dendro/fccz4-application-wiring.md`;
- `wiki/infrastructures/dendro/gridfunctions-naming-and-loops.md`;
- `wiki/infrastructures/dendro/project-assembly-and-emitters.md`;
- `wiki/infrastructures/dendro/validation-standalone-host-and-deferral-gates.md`; and
- `wiki/infrastructures/dendro/index.md` only if its routing text changes.

Revise them for:

- supported profiles and the order-6 default;
- the `fd_order`/`ko_fd_order` distinction;
- centered shift advection;
- Dendro ownership of mesh geometry, storage, communication, and time;
- whole-vector versus block-local offset rules;
- one profile per generated project;
- production libraries versus qualification drivers; and
- any remaining Berger-Oliger qualification limit.

The centered-advection slice removes obsolete Dendro-only upwind state from `nrpy/infrastructures/Dendro/CFunction_roles.py`, `state_h.py`, and `Dendro_defines_h.py`. Generated-source checks must confirm that no Dendro function references the deleted interfaces.

Use direct computational-physics language. Store no source revisions, timestamps, run-result snapshots, counts, or maintenance logs. Modify `.github/workflows/main.yml` only with explicit user permission for that file; this Dendro validation work has that permission.

Run focused doctests for every modified Python module, existing finite-difference/upwind/KO regression tests, the generated standalone matrix, available real-host tests, focused static analysis, formatting checks, `git diff --check`, and `python tools/kb_lint.py` if KB pages change. Compare with the pre-change baseline; add no suppressions, skipped tests, or removed assertions.

## Acceptance criteria

The change is complete when:

- BSSN and fCCZ4 examples default to centered FD6 and accept only FD4/6/8.
- `ko_fd_order` is exactly 2/4/6 and effective KO reach is exactly 2/3/4.
- Required padding is exactly 2/3/4 with KO either off or on.
- No generated BSSN or fCCZ4 RHS contains shift-upwind/downwind operators.
- Existing NRPy callers retain old KO behavior when `ko_fd_order` is omitted.
- Every production kernel consumes per-block runtime geometry supplied by Dendro and assumes no grid size.
- Whole-mesh offsets and block-local zero offsets are each applied exactly once.
- Host-selected block callbacks touch only selected host storage.
- Host/generated profile mismatches fail before numerical evaluation.
- Both formulations build as collision-free libraries against the available Dendro interface.
- All twelve standalone profiles pass their numerical checks.
- The public Dendro-GR reference checkout remains unchanged, and authorized `.github/workflows/main.yml` edits remain limited to Dendro validation.

## Deferred optimizations

After correctness and Dendro compatibility are established, profile before considering stored/shared derivative workspaces, altered common-subexpression elimination, SIMD, block tasking, GPU kernels, or combined multi-order projects. Measure sixth-order time per interior point, memory traffic, generated code size, compile time, and register pressure. Any optimization must preserve one numerical block body, Dendro-owned geometry/storage, the callback interface, and numerical equivalence for all qualified profiles.

## Trialectic reconciliation basis

All three independent seats selected one finite-difference profile per generated project and agreed on centered Dendro-local lowering, host-owned geometry, collision-free libraries, and a 4/6/8 standalone matrix. The proposals differed on FD8 KO reach. The later numerical requirement resolves that difference: `ko_fd_order=fd_order-2`, which NRPy's existing `dKOD` convention turns into effective differences 4/6/8 and reaches 2/3/4. This final plan is authoritative where an earlier independent proposal used another KO convention.
