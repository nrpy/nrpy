# Geodesics audit synthesis — Revision 4

Status: Final — unanimously ratified by all four agents. This is a commissioned
coordination artifact, not an implementation change. It synthesizes
`/tmp/analysis1.md` through `/tmp/analysis4.md`; it does not concatenate them.

## Scope and decision rules

The four reviews covered the equation core, massive and photon generated
runtimes, numerical interpolation and its binary producer, examples, Python
consumers, and the directly relevant NRPy/BHaH infrastructure routed from
`AGENTS.md`.

Use these reachability labels throughout:

- **Default-live:** shipped default generation/runtime reaches the defect or
  misleading behavior.
- **Valid edge:** supported input or geometry reaches it, but defaults usually
  do not.
- **Invalid input:** reject once only when the state is noncoordinate or the user
  explicitly authorizes the named coordinate target. Otherwise document the
  unsupported coordinate precondition and request target-specific authorization;
  do not add hot-loop guards.
- **Latent/unwired:** an advertised/exported surface exists, but no in-repo
  product currently assembles it end to end.
- **Conditional:** act only after the stated measurement or product gate.

Central engineering decision: remove causes, dead state, and unsupported
surfaces. Do not add compatibility shims, generalized frameworks, epsilon
branches, or per-ray defensive checks. Every behavior-changing slice requires
an exact or generated-runtime acceptance test before source-golden updates.

## Executive decision table

| ID | Priority | Confidence | Reach / practical likelihood | Impact | Effort | Must precede / depends on |
| --- | --- | --- | --- | --- | --- | --- |
| GEO-001 | P1, first safety slice | Very high | Default-live read; staged value is unused, so numerical corruption is unobserved; sanitizer failure likely | One-past-end read / C UB | S | Before GEO-012 |
| GEO-002 | P1 | Very high | Default-live byte nondeterminism/data exposure; visibly harmful on failure/no-window valid paths | Indeterminate records and invalid-field consumption | S-M | Before GEO-007 |
| GEO-003 | P1 | Very high | Valid generation/composition edge; shipped aggregator order is deterministic | Supported module composition can silently change integration/termination defaults | M | Before numerical before/after comparisons |
| GEO-004 | P1 | High | Default-live semantic/config debt; `t_integration_max` has no effect | Users configure/observe the wrong quantity | M | Same slice as GEO-003 |
| GEO-005 | P1 | Very high | Valid exact-axis edge; default pixels unlikely to land exactly there | Carter `Q=NaN` at a regular point | M | Before GEO-012 and oracle curation |
| GEO-006 | P1 / partly authorization-gated | Very high | Valid near-collinear or scaled-normal edge; degenerate geometry unsupported | Inconsistent rays/hit coordinates and source classification | M | Before photon orchestration cleanup |
| GEO-007 | P1 mandatory / P3 portability | High | Same-build default format; partial-write failure is environmental; cross-machine persistence is unproven | Corrupt/ambiguous artifacts; conditional portability debt | S-M mandatory; L portable | After GEO-002 invariant |
| GEO-008 | P1 behavior correction | Very high for last-column loss; high for ordinal choice | Default-live renderer defect | Confirmed last-column loss; dense ordinal ownership is the evidenced correction | M | Atomic with GEO-007 Phase A |
| GEO-009 | P2 decision | Very high | Latent/unwired in repository; external consumers uncertain | Large unproved subsystem and NaN-only failure channel | S inventory; M delete; L retain | After live P1 work |
| GEO-010 | P1 if retained | High | Latent until GEO-009 retains numerical mode | Silent scientific misinterpretation of a wrong binary layout | M | GEO-009 retain decision; before GEO-011 |
| GEO-011 | P2 if retained / mapping refactor gated | High | Numerical hot path only | Duplicate boundary logic and potentially dominant rotation cost | M-L | GEO-010 fixture; product and authorization gates |
| GEO-012 | P2, separate slices | High | `Lx/Ly` are default-live dead outputs; unsupported `Numerical` is latent | Misleading/nonconserved outputs and invalid generated API | S-M each | GEO-001; external API/owner inventory; after GEO-005 |
| GEO-013 | P2 | High | Dead/self-validation-only equation APIs; one active BSSN caller | 1,257-line mixed owner, duplicated transforms and goldens | M-L | Delete dead surfaces, reassess; move is conditional |
| GEO-014 | P2 | Very high | Massive default-live | Repeated unused metric work and speculative batch initializer | M | Can follow GEO-003; coordinate equation aliases |
| GEO-015 | P2 local / conditional refactor | High | Photon CPU default-live; CUDA benefit unmeasured | Dead state/no-op/avoidable atomics; conditional pipeline complexity | S local, L refactor | GEO-002/GEO-006A; benchmark only for refactor |
| GEO-016 | P2 | High | Registration/generation paths | Duplicated registrar mechanics and misleading particle distinction; geometry mismatch remains theoretical | M | Coordinate with GEO-013 |
| GEO-017 | P1 cross-cutting | Very high | Existing checks miss confirmed defects; several runners execute zero tests | False confidence and unstable huge oracles | M-L | Exact tests accompany every preceding slice |
| GEO-018 | P2 | High | Valid `M_scale != 1` configuration | Spin parameter documentation changes physical geometry intent | S for docs; L for semantic migration | Decide before final oracle regeneration |
| GEO-019 | P3 conditional | Medium-high | Platform/data dependent | Renderer peak memory/complexity; small visualization failures | S-M | Benchmark or direct fixture first |

Effort: S = focused change, M = several owners/callers, L = generated-product or
format migration.

## Detailed issues and implementation instructions

### GEO-001 — Massive diagnostics stage nine components from an eight-component state

**Evidence and why it matters.**
`nrpy/examples/mass_geodesic_integrator.py:134` declares `double y[8]` and exposes
it through `PhotonStateSoA.f`. Both CPU and CUDA paths in
`nrpy/infrastructures/BHaH/general_relativity/geodesics/conserved_quantities.py:240-264`
copy components `m=0..8`. The shipped massive example calls this before and
after integration. Diagnostic expressions consume only coordinates and four
contravariant tangent components, state indices `0..7`; the ninth photon path
length is never needed. The extra value is currently unused, explaining the
low visible-failure rate, but reading `y[8]` is still C undefined behavior.

**Do not:** pad massive state to nine, add a bounds branch inside the copy loop,
or preserve photon-named universal types. Those hide the ownership error.

**Implementation plan.**

1. Make conserved diagnostics' generated state width exactly eight, since that
   is its invariant mathematical input for both particle types. If another
   shared consumer genuinely needs particle-specific width, generate one
   compile-time width from a single owner; do not use a runtime hierarchy.
2. Generate staging allocation, loop bound, copy size, and comments from that
   one width. Rename the shared view/arguments to `GeodesicStateSoA` or another
   trajectory-neutral name.
3. Leave photon master storage at nine components in the photon integrator; only
   its diagnostic projection narrows to eight.
4. Update generated prototypes, docs, and full-source evidence atomically.

**Reuse.** Keep existing BHaH allocation/copy macros and particle-kernel
generation utilities. No new allocator or gridfunction abstraction fits.

**Tests / acceptance.** Generate the massive example in an isolated tree and
run CPU under ASan/UBSan; require no invalid read and unchanged diagnostic
values/layout. Inspect generated massive diagnostic bounds/allocation as eight.
Run photon diagnostics and prove its integrator still preserves state component
eight. Land this smallest safety slice before GEO-012 removes result fields.

**Priority ruling.** Reachable UB is always a release defect, but P0 here would
imply immediate user-visible or consumed-state corruption that evidence does
not show: staged component eight is unused. Classify P1/high and land first among
safety changes. Promote to P0 if generated code consumes that slot or a default
crash/corrupted used value is reproduced.

### GEO-002 — Blueprint records serialize indeterminate fields

**Evidence and why it matters.**
`photon/batch_integrator_numerical.py:333-335` allocates `d_results_buffer`
without initialization. Window, source, and escape handlers each write only
their own field group; finalization copies the entire prior struct. Photon
`main.py:225-241` shifts every `y_w/z_w` and serializes each full record.
Therefore irrelevant fields are indeterminate under defaults, and valid paths
that terminate before a window event can feed indeterminate coordinates to the
analysis path. Ordinary default rays normally hit the window, so false rendered
pixels are an edge, not an established default outcome.

**Do not:** scatter field checks through render/event hot loops, use sparse
records, or rely on `memset` for enum/floating representation assumptions.

**Implementation plan.**

1. Define one invariant next to `blueprint_data_t`: termination is always valid;
   window/source/sphere groups are valid only for their documented event; all
   other fields are deterministic zero. Preserve `y_w/z_w` as actual window
   intersection coordinates; do not overload them with launch-pixel identity.
2. Use zero-initialized `L_w` as the minimal existing window-valid sentinel.
   Trace confirms `L_w` is written only after successful window intersection;
   affine starts at zero and accepted steps add positive `h`. Enforce positive
   finite ordered noncoordinate step policy
   `0<h_min<=initial_h<=h_max`. Document the existing distinct camera/window
   plane precondition; do not add a coincidence/coordinate guard without the
   request-gated authorization required by KB policy. Add a serialized
   `window_valid` field only if a supported counterexample to `L_w>0` is produced.
3. Initialize every complete record once after allocation using the existing
   generated CPU/CUDA kernel utility.
4. Have each event write its complete group. Apply tile offsets exactly as
   `if (L_w>0) { y_w += offset_x; z_w += offset_y; }`, so an absent window group
   stays deterministic zero. GEO-008 derives image placement independently from
   tile metadata plus record ordinal.
5. Keep one fixed-size record per ray. Document which event owns each group.

**Reuse.** Use
`nrpy.helpers.parallelization.utilities.generate_kernel_and_launch_code`; do not
invent a second launcher.

**Tests / acceptance.** Prefill storage with two hostile byte patterns, then run
source, window, escape, early RK failure both before and after a window event,
time-limit, and generic-failure cases
on CPU and generated CUDA branches. Deserialize and require byte-identical
records across patterns/runs and status-specific field invariants. Require the
owner test to prove no terminal `ACTIVE/REJECTED`, and make consumers reject
such records strictly; do not add a production finalization scan,
unless a live terminal path proves such a check necessary. CPU MSan/ASan and
compute-sanitizer where available. Prove all successful window crossings have
`L_w>0` and all absent crossings retain zero. GEO-007
then formalizes this invariant on disk.

### GEO-003 — Shared integration CodeParameters have multiple owners

**Evidence and why it matters.**
`batch_integrator_numerical.py:44-56` registers `p_t_max=1e3` and
`numerical_initial_h=0.1`; `event_detection_manager_kernel.py:24-33` registers
`p_t_max=1e5`; `rkf45_finalize_and_control_kernel.py:98-118` registers
`numerical_initial_h=1.0`. `nrpy/params.py:491-496` overwrites existing entries,
so alternative supported module composition/order can change public defaults
and runtime behavior. The shipped batch aggregator currently registers event,
finalizer, then batch deterministically, so ordinary default nondeterminism is
not established; the one-owner rule and composition hazard remain real.

**Do not:** merely standardize the duplicate literal values or depend on one
current aggregator order. Leaf kernels must not own shared orchestration policy.

**Implementation plan.**

1. Choose existing setup owners: event policy owns the renamed momentum limit;
   RKF45 configuration owns `numerical_initial_h`; a package/example aggregator
   establishes dependencies and order.
2. Register every shared name once with reviewed type, default, units, and
   semantics. Preserve the shipped batch winners unless a separate behavior
   change is approved: `pU0_max=1e3` and `numerical_initial_h=0.1`. Leaf kernels
   consume symbols without re-registration.
3. Replace example writes to `par.glb_code_params_dict` with
   `par.adjust_CodeParam_default`.
4. Add registry assertions for owner module, metadata, `commondata`, and default
   after each supported CPU/CUDA registration sequence.

**Tests / acceptance.** First enumerate actually supported registration paths.
For each, metadata/defaults must remain identical or missing-owner use must fail
clearly; do not imply arbitrary leaf permutations are supported. Verify the
shipped aggregator's before/after generated commondata/default parfiles and the
two exact preserved values. Complete this
before trajectory/performance comparisons so alternate generated assemblies do
not compare different policies.

### GEO-004 — Photon configuration misnames state `p^0` and advertises a dead time limit

**Evidence and why it matters.**
Equation and reversal code evolve/write contravariant `pU[0]`, but
`set_initial_conditions_kernel.py`, event detection, batch integration, enums,
and dashboards call it covariant `p_t`. The threshold is coordinate-component
blow-up policy, not a conserved energy. Separately, `t_integration_max` is
registered and printed but no control flow reads it; slot-manager bounds govern
the interval. Single photon and analogous massive examples also register a
canonical parfile parser but do not call it, then hard-code runtime values.

**Do not:** retain old aliases or both time endpoints, or mix a behavior change
to the threshold with the semantic rename.

**Implementation plan.**

1. Rename the active contract atomically to `pU0_max` and
   `FAILURE_PU0_TOO_BIG` (final spelling may follow project convention). Because
   this is a source-spelling correction, keep the serialized integer value
   unchanged. Atomically rename C/Python constants, parfiles, dashboard, docs,
   and tests with no shim. Bundle with GEO-007 only if maintainers intentionally
   renumber the code; no current evidence supports renumbering.
2. Delete `t_integration_max` and print the effective slot interval. If product
   requirements instead demand a public end time, derive slot bounds from that
   one value and delete the independent endpoint.
3. In single examples, call `cmdline_input_and_parfile_parser` immediately after
   commondata defaults; move project defaults to generation time and remove
   hard-coded runtime assignments. Remove unused `blueprint_path`.
4. Keep threshold numerical behavior unchanged during the rename.

**Reuse.** Directly use `cmdline_input_and_parfile_parser` and
`par.adjust_CodeParam_default`; wrappers add no value.

**Tests / acceptance.** Generated source must map state indices `4..7` to
`pU0..pU3`; schema/parfile snapshots use only the new name. Two endpoint values
must alter terminal time/status if an endpoint survives, and the dashboard must
show that same value. CLI/parfile overrides must change generated runtime state.

### GEO-005 — Carter `Q` has a removable rotation-axis singularity

**Evidence and why it matters.**
`geodesic_diagnostics/conserved_quantities.py:257-299` constructs
`rho=sqrt(x^2+y^2)`, `p_theta/rho`, and `1/sin(theta)^2`. At the regular point
`a=0,(x,y,z)=(0,0,2),pU=(1,2,3,4)`, the current result is `NaN`; the required
Schwarzschild angular result is `Lx^2+Ly^2=52`. The physical Kerr ring is not
involved (`r^2=4`). A custom axis-aligned trajectory is valid, though default
pixel centers probably miss the exact axis. Near-axis material loss is not yet
demonstrated and must remain a measured question.

With spatial covector `pD`, `r2` the Kerr radius squared, `R2=r2+a^2`,
`A=x*p_x+y*p_y`, and `mu2=1` massive or `0` photon, use:

```text
Q = z^2*R2/r2*(p_x^2+p_y^2)
    - 2*z*A*p_z
    + r2*(x^2+y^2)/R2*p_z^2
    + z^2*a^2/r2*(mu2-E^2)
```

The reviews independently verified algebraic identity off-axis. It removes
only the removable `rho` denominator where `r2>0`. In this Kerr-Schild branch,
`r2=0` across the excluded disk `z=0, x^2+y^2<=a^2`; the curvature singularity
is its boundary ring. The existing analytic metric is already undefined on that
branch disk, so this issue adds no new guard or domain admission.

**Do not:** add epsilon clamps, `Piecewise`, an axis branch, or claim a proven
near-axis production failure.

**Implementation plan.**

1. Lower `pU` to `pD` once and reuse it for `E`, `Lz`, and `Q`.
2. Replace angle-based construction by the Cartesian formula; select `mu2`
   explicitly by particle type.
3. Reuse the analytic owner's Kerr-radius and spin symbols; retain/document the
   existing `r2=0` Kerr-Schild branch-disk exclusion and distinguish its
   physically singular boundary ring. Add no new guard.
4. Apply GEO-012 output contraction only after this identity is accepted.

**Tests / acceptance.** Require exact symbolic old/new equality away from the
old declared denominators on `r2>0`; exact `a=0:Q=Lx^2+Ly^2`, including regular
axis points with `r2>0`, for
both particles; both time roots and signs of `z`; random high-precision off-axis
and multiple-azimuth near-axis samples; generated-C comparison; no generated
`1/(x*x+y*y)` equivalent. Measure SymPy operation/CSE count, codegen time,
generated C size, and runtime so the algebraic simplification does not regress
the pipeline.

### GEO-006 — Camera and source frames have inconsistent ownership

**Evidence and why it matters.**
Camera basis construction appears in `photon/main.py:153-185`,
`set_initial_conditions_kernel.py:390-422`, and
`handle_window_plane_intersection.py:89-125` with different collinearity
thresholds, fallback vectors, and component tests. The interval
`1e-10 <= |up x forward| < 1e-9` deterministically selects different branches.
`handle_source_plane_intersection.py:99-147` does not normalize its plane normal,
so `(0,0,2)` and `(0,0,1)` describe the same plane but change projected
coordinates and annulus membership.

**Do not:** normalize/check inside each ray, or force symbolic SO(3) helpers into
a host/device runtime-frame contract they do not implement.

**GEO-006A — ungated valid-input correction.** Normalize a valid, nonzero source
normal once and reuse that normalized value for plane crossing and projection.
Add no zero/nonfinite check or new fallback. Preserve/document the existing
nondegenerate-normal precondition and current camera/source fallback semantics.
Test invariance under positive normal scaling and characterize existing
threshold behavior without changing it.

**GEO-006B — authorization-gated consolidation.** Only after the user expressly
authorizes the camera/source-frame coordinate target:

1. Add one private, dependency-free
   `orthonormal_basis3_t` containing only three axes and one setup builder used
   at exactly two call sites. Build one camera basis and one independent source
   basis:
   normalize normal/forward, project configured up orthogonally, choose the
   least-aligned Cartesian axis deterministically if collinear, and complete a
   right-handed orthonormal basis.
2. Compute camera axes once from the original master center and camera. Keep the
   per-tile `window_center` origin in its existing owner and pair it with the
   invariant axes at projection; do not rebuild axes per tile. Source origin
   likewise remains separate from its fixed axes.
3. Keep both immutable bases as ordinary setup values, not CodeParameters or
   parfile state. Pass `const orthonormal_basis3_t *` through host orchestration
   and a device-visible/thread-local pointer/value to kernels. Prototype before
   considering 18 scalar arguments, which likely worsen the interface/register
   surface. Remove unused `window_center_out` and `n_*_out` outputs.
4. Add any zero/nonfinite/coincident rejection only if that exact behavior is
   included in the authorization.

**Authorization fence.** Consolidating or altering collinearity fallback and
adding zero/nonfinite/coincident geometry rejection changes coordinate-derived
admission/classification and is request-gated by KB policy. The current user
request does not name that target. Document unsupported degenerate geometry and
seek exact authorization before implementing those parts. Even reusing one of
the existing camera-axis implementations selects behavior in their threshold
disagreement region, so it belongs to GEO-006B. GEO-006A does not authorize new
guards.

**Tests / acceptance.** GEO-006A requires invariant source results under nonzero
normal scaling and no fallback behavior change. GEO-006B, after authorization,
adds generic/principal-axis, near/exact-collinearity, and threshold regression
tests; require unit norms, zero dot products, positive determinant, and identical
axes at former call sites. The private two-instance value is the abstraction
ceiling; do not generalize it into a rotation/frame framework.

### GEO-007 — Blueprint writing/ABI/consumers are unchecked; portability is a separate decision

**Evidence and why it matters.**
Photon `main.py:230-249` uses `sprintf`, unchecked `fwrite`, and
`system("zip ...")`. The packed C record contains an implementation-defined enum
width; Python assumes native `int32` plus ten native-endian doubles. There is no
C offset assertion. Renderer
and analysis duplicate discovery; missing, empty, multiple, truncated, or
partial tile artifacts are accepted inconsistently.

**Do not:** reuse checkpoint/camera payload structs wholesale, silently allow
partial images, or invoke shell compression from simulation C.

**Phase A — mandatory same-build contract.**

1. Add a small local checked-write function, `snprintf`, and checked
   open/write/close. The checkpoint `BHAH_safe_write_impl` is evidence for the
   contract but currently a local generated prefunc, not a callable shared API.
   Factor shared BHaH binary-I/O support only as a separately reviewed follow-up
   if a second live caller makes it simpler overall.
2. Serialize termination as fixed `int32_t`; keep doubles native-endian and use
   explicit Python native dtypes (`=i4`, `=f8`). Define the packed native header
   in this exact order: `char magic[8]`; `uint32_t native_schema_version`;
   `uint32_t header_size`; `uint32_t record_size`; `uint32_t tx`;
   `uint32_t ty`; `uint32_t tiles_w`; `uint32_t tiles_h`;
   `uint32_t scan_density`; `uint64_t record_count`. No reserved or implicit
   padding bytes are allowed. This is same-build structural metadata, not Phase B.
3. Add C `_Static_assert`s for header size and every header field offset, plus
   record size/offset assertions. Match them in Python with native-byte-order,
   standard-size, no-alignment `struct` format `=8sIIIIIIIIQ` (or an exactly
   equivalent dtype) and explicit size/offset assertions. Validate
   `record_count==scan_density*scan_density`, expected tile identity/set, and
   exact header-plus-payload length using checked 64-bit arithmetic. Validate
   header magic/version/sizes, positive discrete dimensions, and tile index
   range as storage/ownership metadata.
4. In one atomic change, migrate the producer, renderer, analysis consumer,
   generated copied scripts, and KB to raw `.bin` plus the shared header reader.
   Delete runtime ZIP/system calls, ZIP discovery, empty/no-result fallbacks,
   and dual raw/ZIP compatibility. If later archival is needed, make it a
   separate post-generation product decision.
5. Create one Python artifact reader shared by renderer and analysis. Require
   exactly one payload per expected tile and reject missing/empty/multiple/
   truncated/misaligned artifacts before rendering.

**Phase B — conditional portable persistence.** First obtain a product
requirement for cross-machine/cross-build archives. If required, extend the
native header contract with endian/real markers and explicit
little-endian encoding or byte swapping in generated C. Merely declaring NumPy
`<f8` does not make raw C writes little-endian. Migrate producer and both
consumers atomically. Do not build a schema generator for one record format.

**Reuse.** Reuse the checkpoint helper's checked-write behavior, locally at
first. The versioned `output_raytracing_data` header is a Phase B pattern, not a
drop-in payload.

**Tests / acceptance.** Phase A: same-build cross-language known-record round
trip, fault-injected checked-write failure, wrong tile/version/header size/
record size/count and truncated record,
complete/missing/empty/multiple tiles, identical consumer acceptance, and no
shell command in generated C. Phase B, only if authorized: wrong endian/version/
record size plus an actual cross-endian encoding fixture.

### GEO-008 — Renderer should map dense ray identity, not floating crossing coordinates

**Evidence and why it matters.**
`render_lensed_image.py:195-204` admits half-open
`y_w_min <= y_w < y_w_max`, then computes
`int(frac*(output_pixel_width-1))`. Because `frac<1`, the last column is not
selected except accidental rounding. Replacing it with `int(frac*W)` is still
unsafe: `nextafter(upper,lower)` can pass admission while normalization rounds
to `1.0`, producing `W` or row `-1`. A coordinate clamp would fix storage but is
request-gated and unauthorized here.

The confirmed defect is last-column loss under current endpoint mapping. The
chosen correction, strongly evidenced but still a product-semantics inference,
uses the generator's stronger structural identity:
`set_initial_conditions_kernel.py:198-253` maps record `i` to dense
`row=i/scan_density,col=i%scan_density` and aims the initial ray at that sample.
Records preserve master-ray order. `blueprint_analysis.py` describes the window
plot as validating this uniform initial distribution. Image pixels therefore
represent initial camera directions; later curved `y_w/z_w` crossings remain
diagnostics, not bin ownership. Existing renderer prose must be migrated because
it currently calls the crossing an endpoint-to-pixel mapping.

**Do not:** add or plan an unauthorized coordinate clamp/admission change, or
overwrite diagnostic `y_w/z_w` with launch coordinates.

**Implementation plan.** Consume GEO-007's native header and record ordinal. For
local `(row,col)` and global tile/sample indices, use checked 64-bit products and
center-preserving integer resampling:

```text
Nw = tiles_w*scan; gx = tx*scan + col
Nh = tiles_h*scan; gy = ty*scan + row
px = ((2*gx + 1)*W) // (2*Nw)
py = H - 1 - ((2*gy + 1)*H) // (2*Nh)
```

Track record offset while streaming chunks. Stop filtering image placement on
`y_w/z_w`; use `L_w>0` only for window diagnostics. This is a documented
behavior-contract correction. If product owners instead require actual curved
crossing coordinates, implementation is blocked pending authorization naming
renderer pixel mapping; then use a proven terminal clamp.

**Tests / acceptance.** Identity resolution, up/downsampling, every tile and
seam, chunk offsets, `1x1`/`2x2`, count mismatch, and checked-product overflow.
Ordinal mapping structurally covers every row/column when
`W<=tiles_w*scan_density` and `H<=tiles_h*scan_density`; assert actual color
population only with an all-renderable-ray fixture because failed rays may leave
pixels blank. All indices remain structurally in bounds. Preserve crossing
coordinates exactly for analysis.

### GEO-009 — Decide whether numerical-spacetime ray integration is a product

**Evidence and why it matters.**
No in-repo product calls `register_CFunction_numerical_interpolation`, initializes
`NumericalTimeWindowManager`, maps slots, and routes results through RK stages.
The analytic photon kernel always calls analytic metric/connection workers.
The exported numerical wrapper returns `void`, collapses stage failures to NaNs,
and has no immediate per-ray status channel. KB prose currently describes an
operational path. Downstream users outside the repository remain possible.

**Do not:** add helpers or optimize the latent path before product ownership is
decided; do not treat NaN as the sole control channel in a retained product.

**Implementation plan.**

1. Inventory downstream/deprecation commitments and name a product owner.
2. If unsupported now, delete the interpolation registrar family and optional
   numerical-window RK cap in one reviewed migration; correct KB claims. Keep
   exporter/combiner only if independently supported.
3. If an owner/current consumer exists, write and approve a product contract
   covering lifecycle, file compatibility, failure signaling, and numerical
   accuracy before implementation investment.
4. If retained, assemble one explicit numerical photon mode: validate/init
   manager, map active slots before launch, prevent remap while readers run,
   clean up on every exit, and return per-ray interpolation status consumed
   immediately by RK orchestration. NaNs may remain diagnostic payload only.

**Tests / acceptance if retained.** A minimal container analytic in time and
space with exactly known interpolation; first/last supported slots, spatial
boundaries, malformed input, deterministic terminal reason, cleanup, ASan, and
tensor/geodesic-step comparison against analytic values. If deleted, no public
choice/docs/registration remains and independent combiner tests still pass.
This P2 decision follows live P1 defects; only release commitment or a runnable
supported assembly promotes it.

### GEO-010 — Numerical container reader trusts unstated binary assumptions

**Evidence and why it matters.**
`time_window_manager_numerical.py` hard-codes offsets and validates selected
sizes/counts, while interpolation assumes exactly 53 little-endian doubles,
3 coordinate values, 10 Cartesian metric values, 40 Cartesian Christoffels,
Spherical/native lookup, axisymmetry, and two phi planes. It does not validate
magic/version/endian/header/real size, record layout, basis, component counts,
geometry version, or symmetry metadata. The context's caller-supplied
`stored_phi_samples` can disagree with header grid metadata.

**Do not:** create a generalized serialization framework or keep independent
writer/reader constants.

**Implementation plan.**

1. Add one compact combiner-owned
   `raytracing_container_schema.py` containing constants and layout enums. The
   Python writer imports them; the manager registrar imports the same values at
   code-generation time and emits literal C checks. Do not create a runtime
   schema abstraction or multi-language generator without a second format.
2. Before `mmap`, validate allowed structural/categorical fields: magic/version,
   endian and f64 size, header/record sizes, 53/10/40 counts, Cartesian
   basis/layout, no payload ghosts, one Spherical grid tag, native-lookup tag,
   z-axisymmetry tag, canonical two-plane IDs, and geometry versions.
3. Numeric checks/canonicalization of `xxmin/dxx/phi`, logical-coordinate ranges,
   or coordinate-derived stencil semantics are request-gated. Deriving
   cell-centered phi values and removing caller-supplied `stored_phi_samples`
   may proceed only after authorization names this numerical interpolation
   target; otherwise document and characterize the existing coordinate contract.
4. Use checked arithmetic for expected byte counts and prove file size before
   `mmap`/offset access.
5. Return specific status classes for I/O, schema incompatibility, time range,
   and mmap failure.

**Tests / acceptance.** Valid minimal container plus negative cases for every
allowed structural/categorical class, including byte-swapped, wrong
basis/layout/count/version and canonical plane IDs. Coordinate-value rejection/
migration tests require the same target authorization. Writer fields and
generated C reader must agree from the same source.
Complete before GEO-011 optimizes transformations.

### GEO-011 — Boundary completion and axisymmetry rotation need a measured ownership decision

**Evidence and why it matters.**
The numerical interpolator maps negative-radius/polar extended nodes and toggles
phi per stencil sample, then expands dense 4D metric/connection arrays and
performs a 4D rank-3 rotation. Defaults imply roughly 130 rotations/ray and
532,480 mostly-zero connection summands. Canonical `bcstruct` already owns
curvilinear ghost mappings, but assumes ghosted registered gridfunctions and
standard parity metadata; current payload is read-only interior mmap AoS with
serialized Cartesian 4D metric/connections. A per-ray `bcstruct` call is both
wrong and slower.

The two stored planes are intended to differ by exactly pi, but reader tolerance
currently also admits genuinely near-pi data. Diagonal signs
`s=(+1,-1,-1,+1)` are exact only if the writer guarantees exact opposite planes
or the reader explicitly canonicalizes accepted metadata. Then
`g'_{mn}=s_m s_n g_mn` and
`Gamma'^a_mn=s_a s_m s_n Gamma^a_mn`; constant Cartesian rotation has no
inhomogeneous term.

**Do not:** call boundary setup/fill per ray, treat Christoffels as ordinary
tensors under general coordinate transformations, use 3D SO(3) rank-2 helpers
as a fake 4D connection replacement, or silently round real near-pi samples.

**Implementation plan.**

1. First establish GEO-010's trustworthy fixture and a categorical canonical
   two-plane axisymmetry contract. Writer and reader derive opposite planes from
   the same schema/plane IDs; do not require bitwise
   `phi1-phi0==M_PI`, which cell-centered floating construction may not satisfy.
2. Make retained interior mmap plus the current mapper the correctness baseline.
   Cross-validate mapped indices against canonical Spherical `bcstruct` only on
   the overlap domain `n<=NGHOSTS`.
3. Apply serialized 10/40 exact-pi sign tables directly during weighted
   accumulation. This removes the 125 node-alignment transforms out of roughly
   130 default rotations. Profile again before touching the five final
   temporal-slice rotations.
4. Retain the current mapper as preferred baseline. `bcstruct` covers a fixed
   `NGHOSTS` shell, while geodesics accepts runtime interpolation half-width
   `n` under `2*n+1<=Nxx`; no current contract proves `n<=NGHOSTS`. Cross-check
   canonical mappings only on the overlap domain. After profiling, a compact
   lookup may precompute either local-mapper results for every accepted raw
   index, or `bcstruct` results only where `n<=NGHOSTS`; it must prove key size,
   full accepted-order coverage, lower LOC/branches, and no public regression.
5. Full 50-field ghost staging is deferred, not an equal first prototype. Try
   it only if the retained/lookup design misses a named performance target and a
   memory budget exists; it still needs explicit 4D connection rotation and may
   add transposition/copy complexity.
6. Keep the residual arbitrary-z dense loop unless post-sign profiling shows it
   material. Only then prototype sparse symbolic `ixp` expressions with
   `ccg.c_codegen(enable_cse=True)` and require smaller/clearer/faster generated
   output.

**Reuse.** Keep existing `Cart_to_xx...` coordinate query and
`interpolation_lagrange_uniform.h`; both are exact fits. `bcstruct` is a mapping
source with an adapter gap. Existing SO(3) helpers cover 3D vector/rank-2 only.
Refactoring/precomputing the coordinate-derived mapper itself is request-gated;
seek authorization naming this interpolation target before implementation.

**Tests / acceptance.** Old/new `g4DD` and `Gamma4UDD` on interior,
negative-radius, both polar crossings, both planes, axis/stencil edge, exact
categorical opposite-plane sign derivation, `+pi/-pi`, and identity. Reject or
migrate noncanonical legacy metadata only after coordinate-target authorization.
Check metric and
lower-connection symmetry and a known axisymmetric analytic metric. Benchmark
generated LOC/op count, codegen and compile time, stack/register use, setup
memory/time, and per-ray runtime at default/max widths. A production lookup is
considered only after product retention, profiling, explicit authorization, and
complete supported-order coverage; ghost staging remains further deferred.
Retain the current residual loop/mapper if replacement lacks a clear win.

### GEO-012 — Separate unsupported diagnostics from result contraction

**Evidence and why it matters.**
Equation/runtime diagnostics compute `E,Lx,Ly,Lz,Q`, but repository consumers
read only `E,Lz,Q`; `Lx/Ly` are not conserved in Kerr. Momentum is lowered in
three separate routines. The API advertises `spacetime="Numerical"`, creates free
metric symbols, but the BHaH generator never hydrates them, so generated C is
not self-contained. Only Kerr-Schild has an in-repo caller. Emitted free symbols
also iterate a Python set.

**Do not:** preserve empty optional fields, compatibility aliases, or a generic
numerical Killing-constant claim without explicit symmetry and metric inputs.

**Implementation plan — do not bundle these migrations.**

1. **Local deterministic cleanup:** sort emitted free symbols by stable name.
2. **Latent-mode decision (P2):** inventory an owner/current consumer, then
   remove `Numerical` diagnostic mode unless one is identified.
   If required, accept one explicit local metric bundle plus declared spacetime
   symmetries; never leave free metric symbols.
3. **Result contraction (P2/public ABI):** after GEO-001/GEO-005 and external
   consumer inventory, reduce result/equation state atomically to nonoptional
   `E_expr,Lz_expr,Q_expr`; keep `Lx/Ly` local only in the Schwarzschild test.
4. Restrict current analytic registration to supported Kerr-Schild and both
   particle types only in the latent-mode slice.

**Tests / acceptance.** Exact result key/field set; old/new `E,Lz`, Carter tests,
fixed massive/photon SymPy values, and early errors for unsupported keys.
Generated C contains no unresolved metric names. Perform public deprecation
inventory before field removal; no shim after migration.

### GEO-013 — Equation core mixes dead generic recipes with the active BSSN boundary

**Evidence and why it matters.**
`geodesics.py` is 1,257 lines spanning analytic equations, numerical/BSSN basis
recipes, and a 315-line validation program. No production caller uses
`symbolic_numerical_christoffel_recipe`,
`symbolic_christoffel_recipe_from_grid_basis`,
`_GAMMA4UDD_METRIC_RECIPE`, or per-instance
`Gamma4UDD_from_metric_recipe`. The one external recipe consumer is
`BHaH/diagnostics/output_raytracing_data.py:143-159`, using two BSSN-specialized
metric/Christoffel functions. Local metric transformation duplicates canonical
`BasisTransforms`; connection transformation has a genuine inhomogeneous term
and is not a tensor transform.

**Do not:** move dead APIs into a new module, leave forwarding aliases, or route
Christoffels through tensorDD helpers.

**Implementation plan.**

1. Re-run the call inventory, then delete the four no-caller surfaces, aliases,
   attributes, validation, trusted keys, and KB claims.
2. Reassess cohesion after deletion. Move the active BSSN functions and
   connection transform to a `numerical_spacetime` owner only if independently
   changing analytic/BSSN responsibilities still remain inseparable or the move
   deletes additional duplication. One caller alone does not justify churn.
3. Whether retained or moved, compose the four-metric transform from canonical `BasisTransforms`: retain
   `g00`, transform `g0i` as covector and `gij` as tensorDD, mirror `g0i`.
4. Keep the connection transform local, including its inverse-Jacobian
   derivative term. Add no generic transform framework without a second caller.
5. Keep numerical transform validation with its actual owner after the cohesion
   decision; keep analytic geodesic invariants in the core module.

**Tests / acceptance.** Exact old/new metric transform components for Cartesian
and Spherical maps; BSSN identity; nonidentity connection fixture; flat
Spherical-to-Cartesian metric/Christoffel; caller import/generation. Any new
module passes repository static checks. Delete dead code before deciding a move.

### GEO-014 — Massive runtime performs unused metric work and overgeneralizes scalar initialization

**Evidence and why it matters.**
Massive equation RHS depends on state and `Gamma4UDD`, not `g4DD`, yet every GSL
RHS call evaluates a metric and passes an unused buffer. `u0_massive` has one
scalar caller but exposes batch/ray indices, two outputs, offload pragmas,
registration-time prints, example-defined indexing macros, and string
replacement over generated C.

**Do not:** retain an unused metric argument for symmetry with photon code, or
wrap a scalar operation in a speculative batch API.

**Implementation plan.**

1. Remove metric hydration/parameter from `calculate_ode_rhs_massive` and metric
   evaluation/buffer from the GSL wrapper. Retain metric use for initial root,
   normalization, and diagnostics.
2. Replace `u0_massive` by a host scalar interface over `metric_local[10]` and
   `state[8]` (or a clean scalar return). One owner writes state index four.
3. Emit local aliases consumed directly by `ccg.c_codegen`; delete text
   `.replace()`, custom batch index macros, offload pragmas, prints, and duplicate
   output assignment.
4. In equation construction, validate particle type first; share acceleration
   and quadratic-root helpers only across their two demonstrated massive/photon
   call sites while preserving symbol basenames and root choices.

**Tests / acceptance.** Generated RHS has state/connection inputs only; exact
RHS comparison and fixed full trajectory/conservation results; RHS evaluation
benchmark. For initialization require flat
`u0=sqrt(1+|u|^2)`, Kerr `g_uu=-1`, forward-root sign, and documented invalid
discriminant behavior.

### GEO-015 — Photon orchestrator carries unmeasured CUDA complexity in the CPU path

**Evidence and why it matters.**
The 1,152-line batch generator has two copies of host/device workspaces and
duplicated current/next pack-transfer-integrate-event-unpack logic. OpenMP
`memcpy_async` is synchronous, so two streams cannot overlap there. Slot-list
mutation occurs in sequential host loops but uses GCC CAS/atomic counters.
Stage six launches an update whose body intentionally does nothing. Four event
history pointers in `PhotonStateSoA` are neither allocated nor consumed;
initial metric scans run even when only telemetry needs them.

**Do not:** delete the CUDA pipeline without measurement, replace particle loops
with structured-grid `simple_loop`, or keep dead fields for hypothetical future
history.

**GEO-015A — proven local cleanup.** After GEO-002/GEO-006A:

1. Remove the stage-six no-op launch while preserving exact stage/status order.
2. Delete the four unallocated/unconsumed event pointers. Gate initial metric/
   state telemetry behind its diagnostic/conservation setting.
3. Re-run repository and generated-call inventory, then independently remove
   host CAS from list operations only if all supported callers remain the five
   observed sequential batch sites. The helpers are static inline definitions;
   keep the public/external claim narrow.

**GEO-015B — benchmark-gated pipeline refactor.**

1. Capture CPU/CUDA wall time and peak memory for representative `1x1` and
   `2x2` tiles and bundle remainders.
2. Generate one chunk-stage sequence and one numerical oracle shared by CPU and
   CUDA. A simple OpenMP path may use one workspace and synchronous copies;
   CUDA may use two workspaces/streams only for measured overlap. Do not hand
   maintain two orchestration truths.
3. Introduce `PhotonChunkWorkspace` only if it reduces generated LOC/argument
   surface without performance loss; it is not mandatory architecture.

**Reuse.** Keep BHaH allocators and the existing particle-kernel generator.
`simple_loop` owns 3D grids and is not a fit.

**Tests / acceptance.** Deterministic slot-list fixture; before/after state,
status, blueprint, and conservation arrays over multiple bundle remainders;
OpenMP TSAN; CUDA race/sanitizer where available. CPU ends with one workspace
only if GEO-015B's benchmark supports that rewrite; GEO-015A removes local dead
work independently. Retained CUDA complexity must have measured benefit.

### GEO-016 — Metric and connection registrars duplicate and can mismatch geometry truth

**Evidence and why it matters.**
`g4DD_metric.py` and `connections.py` independently flatten tensors, scan
symbols, obtain coordinates, generate preambles, and register C. `PARTICLE`
changes only a comment/validation while function names omit it and both state
layouts put coordinates in the same indices. The API could pair expressions
from one geometry with another name, but no live mismatch was found and NRPy
generators commonly accept expressions. Both owners have zero-example runners.

**Do not:** make one giant public registrar or force a geometry-bundle/key-only
API to prevent a theoretical caller lie.

**Implementation plan.**

1. Remove `PARTICLE` from both APIs and descriptions.
2. Validate tensor ranks, dimensions, component counts, and supported coordinate
   hydration at registration.
3. Share only the demonstrated private mechanics: symmetric flattening and
   used-coordinate preamble. Keep separate public C registrations.
4. Add a geometry bundle/key-only boundary only after a reproduced mismatch or
   second consumer contract demonstrates need.

**Tests / acceptance.** Exact flatten order and component counts, deterministic
symbol order, generated compile for supported
geometry. Remove empty doctest runners; do not replace them with placeholders.

### GEO-017 — Validation overfits full generated text and misses decisive contracts

**Evidence and why it matters.**
Equation owners serialize entire `object.__dict__`; massive/photon trusted files
duplicate hundreds of metric/Christoffel/internal keys. Several symmetry checks
prove values explicitly mirrored during construction. Exact root substitution,
axis Carter behavior, runtime state bounds, record initialization, mmap schema,
and interpolation math are absent. Many `__main__` doctest runners execute zero
examples. Source goldens cannot catch GEO-001, GEO-002, GEO-010, or hot-path
behavior.

**Do not:** add ceremonial doctests, accept missing oracle files as success, or
put all runtime cases into one expensive end-to-end test.

**Implementation plan.**

1. Add exact owner contracts first: Kerr radial identity/null vector, root
   substitution and flat signs, Carter identity/axis, transform identity and
   nonidentity, state/result layout, schema fields, and pixel mapping.
2. Curate stable expression dictionaries: analytic metric outputs; shared
   Christoffels only if sampled value remains useful; particle RHS/root;
   diagnostic `E,Lz,Q`. Delete coordinates, duplicated metric, normalization,
   dead recipes, and internal attributes.
3. Keep full-source oracles only for stable handwritten wrappers. Add scoped
   generated compile/runtime tests: massive sanitizer, blueprint determinism,
   and—only if retained—minimal numerical interpolation.
4. Remove zero-attempt runners. Every behavior/API slice co-lands its narrow
   regression; only final relocation/curation of numerical-transform validation
   waits for GEO-013's ownership decision.
5. Regenerate trusted expressions using the required isolated two-process
   procedure: fresh candidate, inspect full diff/math, second fresh comparison.

**Tests / acceptance.** Every test declares what it proves; exact key sets are
asserted; CPU/OpenMP and generated CUDA source branches are covered where they
differ. Runtime CUDA claims require hardware evidence. Do not edit protected
`.github/workflows/main.yml`; request explicit exact-file authorization for any
future CI hook.

### GEO-018 — `a_spin` equations are dimensional but examples call it dimensionless

**Evidence and why it matters.**
The metric adds `a_spin^2` directly to coordinate-length squared and uses it
with `M_scale`; current equations implement Kerr `a=J/M` with length units.
Examples call `a_spin=0.9` dimensionless. This is numerically hidden at
`M_scale=1` but changes intended geometry when users vary mass scale.

**Do not:** silently reinterpret the existing symbol or add a second alias.

**Implementation plan.**

1. Prefer the behavior-preserving decision: document `a_spin` as dimensional
   Kerr `a`, add parameter descriptions/units, and correct examples, telemetry,
   KB, and comments.
2. If maintainers require dimensionless `chi`, make it a separate API/behavior
   migration: compute `a=M*chi` exactly once in the analytic owner, pass that
   canonical value to diagnostics, and atomically update callers/oracles/docs.
3. Ensure diagnostics reuse the analytic owner's parameter symbol rather than
   recreating it by spelling.

**Tests / acceptance.** State the units/affine convention, then test the full
Kerr scaling map: scale `(M,a,t,x,y,z)` and every other dimensional input by one
factor, transform vector/covector momentum according to that convention, and
compare dimensionless metric/invariant observables. Varying `M` alone is not a
covariance test. Require metric and `Q` to consume identical physical `a` and
generated parfile descriptions to match equations. Decide before final GEO-017
oracle regeneration.

### GEO-019 — Conditional Python/visualization simplifications

These are real but lower leverage and should not delay correctness work.

1. Renderer ProcessPool “zero-copy” is fork-dependent; spawn copies large
   textures to up to eight workers. Benchmark sequential/thread/process modes at
   1, 4, and 16 tiles with wall time and peak RSS. Because Numba releases the
   GIL, use a bounded ThreadPool or sequential path only if measured simpler and
   competitive. Otherwise document fork/shared-memory behavior accurately.
2. `visualize_trajectory.py` assumes `np.loadtxt` is 2-D; use
   `np.atleast_2d`, validate columns, and test one-row/empty/malformed files.
3. Correct the horizon note: `r=2M` is exact for Schwarzschild and approximate
   for spinning Kerr. Correct producer azimuth documentation to `[-pi,pi]`.

**Acceptance.** Concurrency choice has measured benefit, bounded memory, and
bitwise-equal image/count accumulators. One-row trajectory input succeeds or
fails with a contract message; mathematical labels match producers.

## Required landing order

For every numbered slice below, its narrow GEO-017 exact/runtime regression
lands in the same change. Only broad oracle curation waits until interfaces
settle.

1. Land independent quick safety defect GEO-001.
2. GEO-003/GEO-004: stabilize parameter ownership and baseline behavior.
3. GEO-002: deterministic records and `L_w` validity proof; then land GEO-007
   Phase A plus ordinal-owned GEO-008 as one atomic raw artifact contract across
   producer, both consumers, copied scripts, and KB.
4. GEO-005 and ungated GEO-006A. GEO-006B remains blocked until exact
   camera/source-frame coordinate-target authorization.
5. GEO-009 P2 product gate after live work. If retained, GEO-010 before GEO-011.
6. GEO-012's deterministic symbol ordering may land locally; its latent
   `Numerical` decision follows owner inventory, while `Lx/Ly` contraction waits
   for GEO-001/GEO-005 and external ABI inventory.
7. GEO-015A may land after GEO-002/GEO-006A; it does not wait for equation work.
8. GEO-013/GEO-014/GEO-016. Run GEO-015B whenever its dedicated CPU/CUDA
   baseline exists; it has no dependency on those equation/registrar changes.
9. GEO-018 decision, followed only by final broad GEO-017 oracle curation.
10. GEO-019 only after higher-priority work or when its benchmark is requested.

Prefer small independently reviewable changes; do not combine algebra changes,
public API removal, and generated-performance refactors into one commit.

## Explicit non-findings and fences

- No sign/index defect was established in geodesic acceleration or the RKF45
  Fehlberg tableau. Massive `(-B-sqrt(D))/(2A)` is the forward-time root for
  `g00<0`; photon `(-B+sqrt(D))/(2A)` is the documented reverse-ray root. Add
  exact root tests, but do not change signs without new evidence.
- Photon Eulerian path length `alpha*p^0` matches its documented signed,
  slicing-dependent diagnostic.
- Temporal Lagrange normalization matches the uniform integer-node convention
  under its stated strictly increasing/uniform precondition.
- The Christoffel transform's inhomogeneous derivative term and sign are
  genuine. Never replace it by an ordinary tensor transform.
- Metric and connection upper-triangle serialization order matched inspected
  consumers.
- Symbolic construction caches are justified; narrow keys/surfaces rather than
  deleting caches for style.
- `acos(z/r)` needs no escape-handler zero guard: supported positive escape
  radius and escape status imply `r>0`. Test that invariant.
- Noncoordinate time-slot invariants use setup validation. Frame normals and
  coordinate guards remain documented preconditions/request-gated targets, not
  new validation in this plan.
- `.github/workflows/main.yml` is protected and must remain unchanged without
  explicit authorization naming that exact file.

## Uncertainty ledger

- Generated projects were not compiled/run by every reviewer. Source dataflow
  proves GEO-001/GEO-002 causes; sanitizer and deterministic binary fixtures are
  still acceptance evidence, not completed work.
- Numerical interpolation has no in-repo assembly, but external consumers may
  exist. GEO-009 requires deprecation inventory before deletion.
- Exact-pi sign transforms are invalid for genuinely near-pi sampled data. They
  require the categorical writer/reader opposite-plane contract; no numeric
  angle canonicalization is authorized by this plan.
- Ghost staging is deferred because it may be worse than read-only mmap AoS.
  GEO-011 first tests categorical signs. Any lookup requires retained-product
  profiling, explicit coordinate-target authorization, and complete supported-
  order coverage; the current mapper remains baseline.
- Removing `Lx/Ly`, API names, or generic recipes requires a final external
  consumer inventory. Repository consumers inspected do not use them.
- CUDA runtime validation depends on available hardware; generated-source-only
  evidence must be labeled as such.

## Dialectic log

### Round 0 — Scribe synthesis choices

- Deduplicated four reports into 19 implementation workstreams and one common
  priority/dependency model.
- Classified default-live, valid-edge, invalid-input, latent, and conditional
  claims so low-likelihood edges do not outrank ordinary behavior.
- Initially treated full ghost completion and retained mmap mapping as equal
  candidates; Round 1 rejected that weighting.
- Made opposite-plane sign optimization conditional on exact/canonical pi
  metadata and a full generated-code/runtime benchmark.
- Preserved negative findings and the protected-workflow fence.

### Round 1 — Hostile review and Revision 1 rulings

All three reviewers rejected Revision 0. The matrix records each material claim;
duplicate attacks remain separate so agreement/dissent is visible.

| Claimant | Core objection | Opponents / support | Evidence | Scribe ruling | Material Revision 1 change |
| --- | --- | --- | --- | --- | --- |
| Photon P1-A | “GEO-003 Default-live overstates reach; shipped batch order is deterministic.” | Root/runtime support | Example registers event, finalizer, batch; params overwrite entries | Accept | P0/default-live -> P1 valid-composition edge; test supported assemblies only |
| Photon P1-B | “GEO-002 lacks an actionable window-valid carrier; launch identity changes semantics.” | Root/runtime support; scribe's launch interpretation rejected | `window_event_found` is runtime-only; `y_w/z_w` documented as intersections | Accept objection; reject new bit unless needed | Preserve intersection semantics; zero-init and use proven `L_w>0` sentinel |
| Photon P1-C | “GEO-007 mixes mandatory safety with unproved portable little-endian v2.” | Root/runtime support | Producer/consumer are local; NumPy dtype cannot endian-convert raw C writes | Accept | Split mandatory native same-build Phase A from conditional portable Phase B |
| Photon P1-D | “GEO-006 carrier is ambiguous and could pollute CodeParameters.” | Root supports abstraction ceiling | Frames are derived setup values, not user parameters | Accept with ceiling | Private `plane_frame_t`, exactly two instances, const pass; scalar fallback only for CUDA |
| Photon P1-E | “Separate CPU/CUDA integrators create numerical drift.” | Runtime/root support | Same RK/event sequence; only copy/synchronization differs | Accept | One generated sequence/oracle; workspace representation is conditional |
| Photon P1-F | “GEO-008 must state the vertical formula; naive inversion can be OOB.” | No dissent | Half-open proof gives `H-1-int(frac_z*H)` | Accept | Exact row/column formulas and positive-dimension setup check |
| Correctness C1 | “Reject GEO-001=P0: copied component is provably unused.” | Runtime supports; photon explicitly prefers P0 | Default read is certain; no consumed value/crash is shown | Accept by impact criterion | P1/high, first safety slice; explicit falsifier for P0 promotion |
| Correctness C2 | “Reject launch-pixel values in event `y_w/z_w`.” | Runtime and later photon support | Curvature means launch sample and crossing can differ | Accept | Removed launch overload; event ownership retained |
| Correctness C3 | “Split checked ABI from version/header/shared extraction.” | Photon/runtime support | Safe writer is local prefunc; portability requirement absent | Accept | GEO-007 two phases; local checked writer first |
| Correctness C4 | “GEO-008 must not depend on artifact-reader work.” | No dissent | Pixel formula is exact, local, default-live | Accept | Independent early landing |
| Correctness C5 | “Add bcstruct-derived compact mapping without copying 50 fields.” | Photon supports; runtime prefers retained baseline/lookup | `bcstruct` can supply mappings while mmap AoS remains interior-only | Accept | Compact setup lookup is preferred candidate; ghost staging deferred |
| Correctness C6 | “GEO-003 default reach unproved.” | Photon/runtime support | Only deterministic shipped order identified | Accept | Same as P1-A |
| Correctness C7 | “No generalized frame abstraction or per-ray checks.” | Photon supports | Only camera/source frames need the operation | Accept | Private two-instance value is abstraction ceiling; setup preconditions only |
| Runtime R1 | “GEO-001 P0/first language conflicts with unused slot and landing order.” | Correctness supports; photon dissents on priority | Exact OOB, unused ninth staged component | Accept priority; retain quick first landing | P1/high, ASan + unchanged diagnostics acceptance; removed unrelated finiteness claim |
| Runtime R2 | “GEO-009 P1 overinvests in an unwired subsystem.” | No dissent | No in-repo registrar/manager/RK assembly | Accept | P2 owner/deprecation gate after live fixes; contract before retention work |
| Runtime R3 | “GEO-010 schema generator overbuilds one format.” | No dissent | One Python writer and one generation-time C reader | Accept | Compact combiner constants/layout module; emit literals; checked file-size arithmetic |
| Runtime R4 | “Full ghost staging does not delete hard 4D transforms and adds copies.” | Photon partially; correctness asks it remain third/deferred design | Read-only mmap AoS, no ghosts, 50 Cartesian fields | Accept | Retained mmap baseline; compact lookup first; ghost staging only after missed target/budget |
| Runtime R5 | “Exact-pi signs remove 125/130 transforms; arbitrary-z rewrite is premature.” | No dissent | Only five final temporal rotations remain after node signs | Accept | Sign table first, profile again, symbolic CSE only if residual material |
| Runtime R6 | “Launch-pixel initialization overloads intersection semantics.” | Correctness/photon support | Field comments and handler define actual crossing | Accept | Same as P1-B/C2; `L_w` sentinel chosen after trace |
| Runtime R7 | “GEO-003 P0 unsupported until two shipped orders differ.” | Photon/correctness support | One shipped deterministic order found | Accept | Same as P1-A/C6 |

**Window-validity adjudication evidence.** `L_w` is written only inside a
successful window handler. Affine state initializes to zero; accepted RK steps
add positive `h`. With `camera != window_center`, the camera is not on the plane
whose normal is `center-camera`, so a real crossing cannot occur at affine zero.
Revision 1 therefore uses zero as absence and `L_w>0` as validity, with positive
step and distinct camera/window setup preconditions. An explicit fixed-width bit
is the fallback only if a counterexample is demonstrated.

**Round 1 status.** No approval claimed. Revision 1 incorporates the accepted
attacks while preserving the photon review's minority P0 preference for
GEO-001 as recorded above. Round 2 must test these rulings with evidence and may
reopen them.

### Round 2 — Rebuttal, new falsifiers, and Revision 2 rulings

Round 2 conceded the main Revision 1 priority/architecture rulings, then found
new contract and floating-point defects. “Accept” below means the document was
materially changed, not that the full synthesis is approved.

| Claimant | Claim / concession | Support / dissent | Evidence | Scribe ruling | Revision 2 change |
| --- | --- | --- | --- | --- | --- |
| Correctness A | Terminal `ACTIVE/REJECTED` must be a test invariant, not a production scan | Photon supports | No live terminal path proved | Accept | Owner test + strict consumer rejection only |
| Correctness B | Camera axes are master-constant while tile origin changes | Photon/runtime support | `main.py` shifts center per tile; handlers project from local center | Accept | Axes-only basis; origins remain in existing owners |
| Correctness C | Headerless Phase A cannot self-validate expected count | Photon/runtime support | Raw payload has no count/config authority | Accept, superseded structurally | Minimal native header now carries tile/scan/count |
| Correctness D | CAS deletion requires final caller inventory | Photon provides current inventory | Five observed calls are sequential, but generated/external claims need care | Accept | Re-run repo/generated inventory before deletion |
| Correctness E | Scaling `M` alone is not Kerr covariance | Photon/runtime support | All dimensional quantities must scale consistently | Accept | Full stated scaling map and momentum convention required |
| Correctness F | Exact tests must land with each slice | No dissent | Deferred validation permits unproved intermediate changes | Accept | Landing order now makes scoped GEO-017 tests atomic |
| Photon R2-1 | Concede `L_w>0`; remove production final-status guard | Correctness/runtime support | Root lies in nondecreasing sign-bracketed affine interval; first crossing positive | Accept | Positive step policy, sentinel tests, no production scan |
| Photon R2-2 | Rev1 float pixel formula fails at `nextafter(upper,lower)` | Runtime supports | Normalized fraction can round to 1.0 | Accept defect; reject unauthorized clamp | Replace coordinate mapping with ordinal integer ownership |
| Photon R2-3 | `plane_frame_t(origin+axes)` is wrong/overwide | Correctness/runtime support | Camera basis fixed, tile center varies | Accept | `orthonormal_basis3_t`; origin separate; pointer/value prototype |
| Photon R2-4 | Phase A count claim needs header or caller metadata | Correctness/runtime support | Visualization lacks authoritative scan count | Accept header alternative | Native tile/scan/count header in Phase A |
| Photon R2-5 | Current CAS callers are exactly five sequential batch sites | Correctness qualifies public scope | Repository search and line inventory | Accept with recheck | Narrow conditional CAS removal |
| Photon R2-6 | Concede GEO-001 P1 action | Correctness/runtime already support | P1 label no longer delays first landing | Record concession | No change |
| Photon A | GEO-012 improperly bundles latent mode, public ABI contraction, Carter work | No dissent | Different reach, owners, and dependencies | Accept | Three independent P2/local slices; pD math remains GEO-005 |
| Photon B | GEO-016 geometry bundle prevents only a theoretical lie | No dissent | No mismatch found; expression-taking APIs are common | Accept | Remove dead `PARTICLE`; validate shape; defer bundle/key-only boundary |
| Photon C | GEO-018 covariance must transform every dimensional input | Correctness/runtime support | Dimensional Kerr scaling | Accept | Exact full-scale acceptance wording |
| Runtime W1 | `L_w` sentinel works, but camera/window rejection is unauthorized | Correctness agrees on no new guard | KB request-gates coordinate admission | Accept | Document distinct-plane precondition; validate only noncoordinate step policy |
| Runtime W2 | `bcstruct` lookup cannot cover accepted `n>NGHOSTS` | Root and focused peers support retained mapper | Geodesics checks `2*n+1<=Nxx`; bcstruct maps fixed ghost shell | Accept | Current mapper baseline; overlap cross-check only; mapping refactor authorization-gated |
| Runtime W3 | Opposite-plane contract must be mathematical, not bitwise pi | No dissent | Cell-centered floating construction can round differently | Accept | Categorical plane IDs/shared derivation |
| Runtime W4 | Headerless count is not intrinsic | Correctness/photon support | Exact count needs external metadata | Accept native header | Same as R2-4 |
| Runtime N1 | Rev1 pixel formula is not IEEE-safe | Photon supports | `nextafter` produces W/-1 | Accept defect; structural alternative chosen | Ordinal integer formula; no coordinate clamp/check |
| Runtime N2 | Frame carrier conflates basis and varying origin | Correctness/photon support | Tile center mutates | Accept | Same as R2-3 |
| Runtime N3 | GEO-013 precommits a module move after deletion | No dissent | Dead-code deletion may restore cohesion by itself | Accept | Reassess after deletion; move only with demonstrated remaining boundary |
| Runtime N4 | GEO-018 scale test underspecified | Correctness/photon support | Scale all units/inputs consistently | Accept | Same as Correctness E |

**Ordinal pixel-ownership adjudication.** Both runtime and photon agents traced
the same invariant: ray `i` is launched at dense sample
`(i/scan_density,i%scan_density)`, event writes retain `master_idx`, and output
serialization preserves order. Image pixels physically label initial camera
directions. No adaptive, sparse, compacted, or reordered result path was found.
Revision 2 therefore adopts a minimal native header and center-preserving
integer ordinal mapping. `y_w/z_w/L_w/t_w` remain curved window-event
diagnostics. This corrects current renderer/documentation behavior; it is not
presented as a no-op. Falsifiers are an explicit product requirement that later
curved crossing owns the detector pixel, adaptive/multiple samples per output
pixel, result reordering/compaction, or downstream dependence on crossing-based
accumulation. None was found.

**Authorization adjudication.** The KB expressly gates changes that validate,
clamp, classify, or alter control flow for coordinates or coordinate-derived
indices. Revision 2 removes the proposed renderer clamp, new coincident-geometry
check, and unconditional mapper refactor. Integer storage arithmetic over
serialized ordinals/counts remains allowed. Frame fallback consolidation and
interpolation-mapper refactoring are marked blocked until the user names those
exact coordinate targets.

**Round 2 status.** GEO-001 priority, GEO-003 reach, GEO-007 phase split,
GEO-009 priority, GEO-010 authority, and sign-first rotation ordering now have
explicit concessions from prior opponents. Surviving alternatives for Round 3
are chiefly product-policy choices: ordinal versus crossing-owned pixels if a
contrary requirement appears; delete versus retain numerical interpolation;
and whether post-deletion equation cohesion warrants a module move.

### Round 3 — Consistency adjudication and Revision 3 rulings

All reviewers conditionally accepted the settled mathematics/priority decisions
but withheld approval for contract contradictions. Revision 3 applies every
blocker without adding architecture.

| Claimant(s) | Blocker | Evidence / criterion | Ruling and exact Revision 3 change |
| --- | --- | --- | --- |
| Runtime B1 / correctness 2 | Global invalid-input rule and frame-normal nonfinding contradicted request-gated coordinate policy | KB expressly gates coordinate admission/classification | Invalid-input rule now rejects only noncoordinate/authorized states; frame guards are documented/gated |
| Runtime B2 | Absent window group would lose zero invariant when main shifts every record | Current main offsets all `y_w/z_w` | Tile shift is now exactly conditional on `L_w>0` |
| Runtime B3 / photon 2 | GEO-004 implied numeric enum migration and conflicted with landing order | Source rename does not change serialized integer | Preserve integer value; atomically rename spelling/docs; bundle only if an intentional future renumber occurs |
| Runtime B4 / photon 1 / correctness 1 | Native header fields and raw transport migration were ambiguous | Reader needs parseable sizes/count authority | Explicit magic/version/header size/record size/tile grid/scan/count fields; raw producer+both consumers+copied scripts+KB migrate atomically; no ZIP compatibility |
| Runtime B5 / photon 9 | GEO-008 overclaimed ordinal semantics and color population | Last-column loss is proven; ordinal ownership is inferred from dense launch/master index | Separate confirmed defect from chosen correction; retain falsifiers; structural coverage and all-renderable color fixture wording |
| Runtime B6 | GEO-010 mixed allowed storage checks with gated coordinate-value checks and false GEO-012 dependency | KB distinguishes discrete storage from coordinate classification | Structural/categorical validation allowed; numeric coordinate/stencil canonicalization gated; dependency only on GEO-011 |
| Runtime B7 / correctness 3 / photon 5 | GEO-011 retained stale near-pi/lookup claims | Canonical plane IDs replace tolerant angle semantics; `n` may exceed `NGHOSTS` | Categorical-ID tests; current mapper baseline; lookup only after product/profile/auth/domain gates; ledger aligned |
| Runtime B8 / photon 4 / correctness 4 | GEO-006 blurred valid normal scaling with gated fallback/admission consolidation | Scaling fix needs no classification; threshold selection does | Split GEO-006A ungated normalization from GEO-006B authorization-gated axes/fallback work; only A lands now |
| Runtime B9 / photon 7 / correctness 6 | GEO-017 said tests wait for GEO-013 | Every behavior/API change needs its own proof | Narrow regression co-lands with every slice; only final ownership relocation/curation waits |
| Runtime B10 / photon 6 | GEO-007+008 and GEO-015 dependencies were inconsistent | Header/ordinal reader is one contract; local photon deletions do not depend on equations | GEO-007+008 atomic; GEO-015A early local cleanup, GEO-015B independent benchmark refactor |
| Runtime B11 / photon 3 | GEO-003 did not choose behavior-preserving defaults | Shipped deterministic aggregator resolves to `1e3` and `0.1` | Values stated explicitly with before/after generated baseline acceptance |
| Runtime B12 | GEO-015 mixed proved deletions with speculative workspace rewrite | Stage-six/dead fields/CAS inventory are local; streams need measurements | Split A/B and removed unrelated equation dependencies |
| Runtime B13 | GEO-008 title treated chosen semantics as previously proven defect | Current renderer intentionally used crossings | Title/evidence now call ordinal mapping an evidenced product correction, not historical fact |
| Runtime N3 / correctness 5 | GEO-013 still risked precommitting a move | Deletion may restore cohesion | Table and plan make move conditional after post-deletion reassessment; only one GEO-013 heading exists |
| All | Final plan must not invent consensus on open product choices | No numerical product owner/consumer commitment found | GEO-009 remains delete/retain gate; crossing-owned renderer alternative and conditional module move retain explicit falsifiers |

**Round 3 status.** Revision 3 is pending ratification, not approved by the
scribe. Settled claims remain: Carter Cartesian formula, GEO-001 P1-first,
`L_w` sentinel, behavior-preserving parameter defaults, ordinal mapping as the
best evidenced renderer correction, numerical interpolation as a P2 product
gate, retained interpolation mapper, categorical sign optimization, and
conditional equation-module movement. Open choices are intentionally actionable
decision gates rather than fabricated consensus.

### Ratification audit — photon blockers and Revision 4 rulings

| Blocker | Evidence | Ruling / Revision 4 change |
| --- | --- | --- |
| GEO-007 header types/padding remained implicit | Cross-language same-build parsing still requires exact signedness, order, size, and offsets | Exact packed `magic[8]`, eight ordered `uint32_t` fields, and `uint64_t record_count`; no reserved/padding bytes; C header+record size/offset assertions and matching Python `=8sIIIIIIIIQ` assertions |
| GEO-005 called `r2=0` only the physical ring | At `z=0`, the chosen Kerr-Schild root gives `r2=0` throughout `x^2+y^2<=a^2`; curvature singularity is the boundary ring | Carter rewrite domain is now explicitly `r2>0`; document the already-excluded branch disk and physical boundary ring; add no guard |

**Revision 4 status.** Unanimously ratified by all four agents. No priority,
formula, implementation ordering, authorization fence, or open product gate was
otherwise changed.

### Final ratification

- Equation-core simplification/scribe agent: approved Revision 4 after the
  renewed header, Kerr-domain, consistency, and status audit.
- Runtime simplification agent: approved Revision 4 after full active-plan and
  dialectic-tail verification.
- Photon simplification agent: approved Revision 4 after exact header and
  Kerr-domain corrections.
- Correctness agent: approved Revision 4 after independent header-offset,
  branch-disk, heading, workflow-fence, and active-plan verification.

Unanimous approval preserves the explicit open gates: numerical-product
delete/retain ownership, target-specific coordinate authorization, optional
portable format, crossing-owned renderer falsifier, benchmark-only runtime and
rotation work, external API inventories, and conditional post-deletion module
movement.
