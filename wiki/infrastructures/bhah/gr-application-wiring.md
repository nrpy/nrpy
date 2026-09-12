# GR Application Wiring

> Map how BHaH registers generated CFunctions that connect GR equations, initial data, diagnostics, and basis transforms. Status: confirmed. Last reconciled: 09-12-2026
> Up: [BHaH](index.md)

## Summary

BHaH GR wiring is a code-generation layer over the symbolic GR modules. It does
not rederive the equations. It chooses coordinate systems and feature flags,
pulls symbolic expressions from the BSSN or fCCZ4, ADM, Psi4, and initial-data modules,
wraps them in BHaH loop/kernel infrastructure, and registers concrete
CFunctions such as `rhs_eval`, `Ricci_eval`, `constraints_eval`,
`initial_data`, `diagnostic_gfs_set`, and `psi4`.

The black-hole and GRHD examples show the dataflow shape: register initial-data
import/conversion, grid and diagnostic helpers, reference-metric precompute,
Ricci/RHS/algebraic constraint projection/constraints, Method of Lines step glue,
coordinate and basis transforms, wrapper dispatchers, headers, parser, `main`,
and cleanup.

## Detail

`register_CFunction_rhs_eval(..., enable_fCCZ4=False,
enable_YBS_Gamma_constraint_adjustment=False,
enable_YBS_momentum_constraint_adjustment=False)` generates the shared RHS
CFunction and preserves BSSN as the public default. The fCCZ4 opt-in branch pulls
non-gauge expressions from `fCCZ4_RHSs.get_rhs(...)` and gauge expressions from
`fCCZ4_gauge_RHSs`; the default branch keeps `BSSN_RHSs[...]` and
`BSSN_gauge_RHSs`. Both branches copy the cached owner dictionary before adding
gauge, dissipation, CAHD, or slow-start-lapse terms. It builds a sorted local name-to-expression dictionary,
maps each RHS name to the matching `rhs_gfs` gridfunction with
`BHaHGridFunction.access_gf`, and emits an interior `simple_loop` with finite
difference codegen, optional SIMD/CUDA intrinsics, optional finite-difference
helper functions, and upwinding controlled by beta. The generated function
signature changes with `enable_rfm_precompute`: it receives either
`rfmstruct` or coordinate arrays. It can include `RbarDD` gridfunctions,
`T4munu`, Kreiss-Oliger dissipation, curvature-aware KO, Hamiltonian-constraint
damping, and slow-start lapse through code-generation flags and commondata
parameters.

Claim evidence:
- Claim: `register_CFunction_rhs_eval` defaults to BSSN and accepts `enable_fCCZ4=True` to select the fCCZ4 non-gauge and gauge owners, while preserving the same generated CFunction boundary and applying optional local terms to a copied expression dictionary.
- Role: public/scientific contract
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`
- Corroboration: [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py), formulation-selecting call
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, GCC 13.3.0, GNU Make 4.3; backend=BHaH OpenMP; precision=double; GPU=not-run; restart=not-applicable; distributed=not-run; error_path=not-run; options=default BSSN generation/build, enable_fCCZ4=True generation/build and one-step startup; date=08-28-2026`

The historical YBS Gamma option applies to either formulation. When enabled,
the registrar leaves the coordinate/options cache string unchanged, registers
the runtime `YBS_chi` parameter after the parallel-registration guard, and
forwards the Boolean to the selected nongauge and gauge owners. BSSN defines
the symbolic addition once; fCCZ4 reuses that adjusted base, with the
connection slot interpreted as `LambdatildeU`. Each formulation keeps enabled
expressions in a separate internal cache without encoding the Boolean in the
coordinate name. The equation layer uses a plain real SymPy symbol, while this
BHaH registrar owns the conditional runtime parameter. The standalone and
superB black-hole generators expose no command-line switch for this option:
each declares an explicit false source constant and forwards it to this shared
registrar.

Claim evidence:
- Claim: the shared registrar conditionally owns the runtime YBS parameter and forwards the opt-in flag to the selected BSSN or fCCZ4 equation and gauge owners; the four black-hole generators keep it explicitly false by default and share this registrar.
- Role: descriptive behavior
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`
- Corroboration: [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py), [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py), [superB_two_blackholes_collide.py](../../../nrpy/examples/superB_two_blackholes_collide.py), and [superB_blackhole_spectroscopy.py](../../../nrpy/examples/superB_blackhole_spectroscopy.py), forwarded `enable_YBS_Gamma_constraint_adjustment` constants
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=BHaH and superB source wiring inspected only; precision=not-applicable; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=default-disabled callers plus enabled BSSN and fCCZ4 source branches; date=08-28-2026`

The separate `enable_YBS_momentum_constraint_adjustment` option controls the
default-disabled timestep-scaled momentum adjustment for either formulation.
When enabled, the registrar adds the shared raw-spacing `DSMINGF` auxiliary
gridfunction and runtime `C_YBS_mom` parameter. The selected equation owner
changes existing `a_rhsDD` outputs; no evolved cleaner state, cleaner RHS,
initial-data path, boundary path, or KO route is added.
`blackhole_spectroscopy.py` keeps a false source constant, forwards it only to
RHS registration, and registers and schedules one generic local-spacing helper
when either CAHD or YBS-MOM is enabled. CAHD consumes that same raw spacing in
its RHS coefficient. The continuum equations and validation limits remain
owned by [YBS-MOM Timestep-Scaled Momentum Adjustment](../../equations/general-relativity/ybs-momentum-damping.md).

Claim evidence:
- Claim: `register_CFunction_rhs_eval` exposes an independent default-false YBS-MOM option for BSSN or fCCZ4, conditionally owns shared raw `DSMINGF` and `C_YBS_mom`, changes existing `a_rhsDD` expressions without new evolved state, and shares one local-spacing helper with CAHD in the black-hole spectroscopy example.
- Role: descriptive behavior
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`
- Corroboration: [dsmin_gf.py](../../../nrpy/infrastructures/BHaH/general_relativity/dsmin_gf.py), `register_CFunction_dsmin_auxevol_gridfunction`; [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py), shared CAHD/YBS-MOM registration and scheduling gate; [representative BHaH rhs_eval trusted output](../../../nrpy/infrastructures/BHaH/general_relativity/tests/rhs_eval_OnePlusLog_GammaDriving2ndOrder_Covariant_SinhSpherical_RbargfsFalse_T4munuFalse_ImprovementsFalse.py), `trusted_dict`
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, SymPy 1.14.0; backend=BHaH symbolic RHS registration plus OpenMP/CUDA source generation; precision=30-significant-digit deterministic trusted sampling and exact symbolic source; GPU=generation pass, execution not-run; restart=not-applicable because no new state; distributed=not-run; error_path=pass for unsupported non-fisheye GeneralRFM spacing; options=both YBS Gamma and YBS-MOM enabled jointly in all eight existing BHaH rhs_eval trusted cases; public defaults remain disabled; date=08-30-2026`

`register_CFunction_Ricci_eval` emits `Ricci_eval` from
`BSSN_quantities[CoordSystem + "_rfm_precompute"].Ricci_exprs` by default. It always uses
the rfm-precompute expression family, stores either ordinary `RBARDD*GF`
auxiliary gridfunctions or `DIAG_RBARDD*GF` channels for the host-only
diagnostics version, and wraps the interior loop with BHaH kernel/launch code.
CUDA generation is rejected for `GeneralRFM`; `host_only_version=True`
temporarily forces OpenMP generation so CUDA applications can still compute
host-side diagnostic Ricci data as `Ricci_eval_host`.

`register_CFunction_hDDdD_eval` (`hDDdD_eval.py`) stores the first derivatives of
`hDD` as the `SCRATCH` gridfunctions `hDDdD`, each direction over the interior grown by
`fd_order/2` points in the directions transverse to its own stencil, and
`register_CFunction_Ricci_eval(..., enable_hDDdD_gridfunctions=True)` then reads them through
an extra `scratch_gfs` argument, rebuilding every mixed second derivative of `hDD` as a single
first derivative of the stored gridfunction. The registration helper returns the eligible
gridfunction names, which Ricci passes as `stored_first_derivatives` to `c_codegen()`.
Finite-difference lowering selects their reads and stencils while retaining the original
mathematical derivative temporaries. Equation construction and cache keys are unchanged. No additional array is allocated during evolution for `hDDdD`: `blackhole_spectroscopy.py` calls
`hDDdD_eval(params, RK_INPUT_GFS, RK_OUTPUT_GFS)` and
`Ricci_eval(params, rfmstruct, RK_INPUT_GFS, RK_OUTPUT_GFS, auxevol_gfs)` before `rhs_eval`
overwrites `RK_OUTPUT_GFS` in the same substep, and enables this only for CUDA double
precision, where it was measured to help.

Claim evidence:
- Claim: Stored-derivative selection belongs to finite-difference lowering. Ricci and RHS consumers explicitly supply the registered fields available to each kernel; centered first derivatives read storage, canonical mixed second derivatives apply a first-derivative stencil to it, and diagonal/upwind/KO derivatives retain their operators. Equations and equation caches do not select storage.
- Role: descriptive behavior
- Deciding authority: [finite_difference.py](../../../nrpy/finite_difference.py), `select_stored_first_derivatives`; [c_codegen.py](../../../nrpy/c_codegen.py), `gridfunction_management_and_FD_codegen`
- Corroboration: [Ricci_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/Ricci_eval.py) and [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), explicit consumer selections; owner doctests verify selection and polynomial evaluation.
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3; backend=BHaH CUDA; precision=double; GPU=RTX 4060 Ti; restart=not-run; distributed=not-applicable; error_path=unregistered selection doctest; options=SinhCylindrical, no symmetry axes, FD8, 64x64x128, both stored options, nearest checkpoint to t=0.4 is t=0.40063545122319866, 18136512 finite values compared with committed option-on build, maximum absolute difference 9.12e-14 within atol=rtol=1e-12; date=09-11-2026`


Both `register_CFunction_Ricci_eval` and `register_CFunction_rhs_eval` expose a
default-false `enable_cpu_tiling` choice. Ordinary registrations therefore emit only
the full-grid functions. In `blackhole_spectroscopy.py`, one eligibility predicate
enables both tile producers, registers `rhs_eval_with_Ricci`, and replaces the
separate full-grid calls with that coordinator; the predicate requires OpenMP,
separate Ricci/RHS evaluation, reference-metric precompute, BSSN, and neither
stored-derivative option. The full-grid functions remain registered for diagnostics
and other callers.

If explicitly enabled for OpenMP,
`register_CFunction_diagnostic_gfs_set(..., enable_hDDdD_gridfunctions=True)` computes
fresh derivatives from `y_n_gfs` in temporary scratch storage, passes that storage to
Ricci, and frees it before evaluating constraints. CUDA diagnostics retain the
unstored `Ricci_eval_host` path. The example registers the CPU tiled scheduler only
when neither stored-derivative option is enabled, matching its call-site gate.

Claim evidence:
- Claim: Ricci and RHS tile registration is default-disabled; eligible OpenMP spectroscopy uses one predicate for both tile producers, the `rhs_eval_with_Ricci` coordinator, and call replacement while retaining full-grid functions. With stored hDD derivatives explicitly enabled for OpenMP, GR diagnostics computes fresh derivatives from the current solution in temporary scratch, passes scratch to Ricci, and frees it before constraints. CUDA diagnostics retain the unstored host Ricci path. Compilation does not establish runtime numerical correctness.
- Role: descriptive behavior
- Deciding authority: [Ricci_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/Ricci_eval.py), `register_CFunction_Ricci_eval`; [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval` and `register_CFunction_rhs_eval_with_Ricci`; [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostic_gfs_set.py), `register_CFunction_diagnostic_gfs_set`
- Corroboration: [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py), shared tile/coordinator predicate and call replacement; [hDDdD_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/hDDdD_eval.py), `register_CFunction_hDDdD_eval`, supplies fresh derivatives and stencil halos
- Validation: `inspected=pass; generated=pass; built=pass; run=not-run; result_checked=not-run`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, GCC 13.3.0, GNU Make 4.3; backend=BHaH and superB OpenMP generation, BHaH OpenMP build; precision=double; GPU=not-run; restart=not-run; distributed=not-applicable; error_path=allocation failure not-run; options=default blackhole spectroscopy tile coordinator, default two-black-hole collision, Ricci-only TOV, default superB collision, plus isolated stored-derivative configuration; date=09-12-2026`

`register_CFunction_cfdD_alphadD_vetUdD_eval` (`cfdD_alphadD_vetUdD_eval.py`) applies the
same tensor-product identity to the right-hand sides. It stores the first derivatives of `cf`,
`alpha` and `vetU` that their fifteen mixed second derivatives are built from as the `AUXEVOL`
gridfunctions `cfdD`, `alphadD` and `vetUdD`, and
`register_CFunction_rhs_eval(..., enable_cfdD_alphadD_vetUdD_gridfunctions=True)` selects
the registered storage through the same `c_codegen()` option. BSSN/fCCZ4 RHSs,
their gauge equations, and CAHD constraints retain their mathematical derivative
symbols and ordinary cache entries. The producer gets its output expressions from
`cfdD_alphadD_vetUdD_gridfunction_expressions()`, without rewriting completed equations.
Unmixed and upwind derivatives keep their original stencils. Only the
directions a mixed second derivative differentiates are stored, which is `partial_0` and
`partial_1` with the index pair canonicalized to `j < k`, so ten gridfunctions are registered
and each is produced over the interior grown by `fd_order/2` points in the directions that
differentiate it. These cannot be `SCRATCH` gridfunctions like `hDDdD`: `rhs_eval` reads them
with a stencil while writing the Method of Lines buffer, so a pointwise store would overwrite a
neighbor's stencil point. They therefore cost memory: `NUM_AUXEVOL_GFS` goes 6 to 16, about 61 MB
and 10% of the run's footprint on the standard grid, in device memory and again in the host
mirror. The BHaH BSSN examples that build for CUDA (`blackhole_spectroscopy.py`,
`two_blackholes_collide.py`, `spinning_blackhole.py`, `hydro_without_hydro.py`,
`kasner_exact_evolution.py`) expose it as `enable_cfdD_alphadD_vetUdD_gridfunctions_for_GPU`
next to their other options; it defaults to off because of the memory cost, requires `--cuda`,
and when set calls `cfdD_alphadD_vetUdD_eval(params, RK_INPUT_GFS, auxevol_gfs)` before
`rhs_eval` within each substep.

`register_CFunction_constraints_eval` emits the diagnostics-side Hamiltonian,
momentum, and conformal connection-constraint evaluator. It temporarily forces
OpenMP, reads
`BSSN_constraints[CoordSystem + "_rfm_precompute_RbarDD_gridfunctions" +
optional "_T4munu"]`, writes `H`, the physical momentum-constraint magnitude
`sqrt(gamma_ij M^i M^j)`, and the conformal connection-constraint magnitude to
`DIAG_HAMILTONIANGF`, `DIAG_MGF`, and `DIAG_LAMBDA_CONSTRAINTGF`, respectively,
and places the function in the `diagnostics/` subdirectory. When the original
parallelization is CUDA, generated references to `RBARDD` and optional `T4UU`
auxiliary gridfunctions are rewritten to diagnostic channels so host-side
constraint evaluation consumes the diagnostic buffer filled for output.

Claim evidence:
- Claim: BHaH `register_CFunction_constraints_eval` writes Hamiltonian, `sqrt(BSSNconstraints.Msquared)`, and `BSSNconstraints.LambdaConstraintMagnitude` to `DIAG_HAMILTONIANGF`, `DIAG_MGF`, and `DIAG_LAMBDA_CONSTRAINTGF`; CUDA host-side generation retains the diagnostic-buffer input rewrites.
- Role: descriptive behavior
- Deciding authority: [constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py), `register_CFunction_constraints_eval`
- Corroboration: [BSSN_constraints.py](../../../nrpy/equations/general_relativity/BSSN_constraints.py), `BSSNconstraints.__init__`; [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostic_gfs_set.py), `register_CFunction_diagnostic_gfs_set`
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=BHaH OpenMP registration with CUDA rewrite path inspected; precision=symbolic code generation; GPU=not-run; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=Cartesian registration contract plus source inspection of optional T4munu and CUDA rewrites; date=08-28-2026`

`register_CFunction_enforce_detgbar_equals_detghat_trAzero` emits the combined
algebraic-constraint CFunction. In one all-points loop it loads all independent
`hDD` and `aDD` components, rescales reconstructed `gammabarDD` so
`det(gammabar)=det(gammahat)`, inverts that corrected metric, and then projects
`AbarDD` trace-free using the same metric. All twelve corrected components are
formed before stores. The generated kernel accepts either `rfmstruct` or
coordinate arrays according to `enable_rfm_precompute`, plus read-only
`auxevol_gfs`.

`register_CFunction_initial_data(..., enable_conformal_projection=True)` adds
projection after checkpoint boundary/interpatch repair and after fresh-data
boundary handling; the option defaults to `False`. Its independent
`enable_fCCZ4=False` option also preserves existing callers. When enabled, the
ADM converter writes `Theta_fCCZ4=0` during fresh conversion. The checkpoint
branch runs before that converter and returns after repair/projection, so loaded
Theta storage is not replaced by fresh-data initialization.

The spectroscopy generator enables projection for both its default BSSN and
opt-in fCCZ4 paths and places the same combined projector in its caller-supplied
Method of Lines `post_rhs_string`. Separately, the collision example enables
the same two projection hooks but remains explicitly BSSN-specific. Neither
initial-data registration nor Method of Lines enables projection for every
caller.

Claim evidence:
- Claim: BHaH preserves default-disabled `enable_conformal_projection` and `enable_fCCZ4` initial-data options; fresh fCCZ4 conversion initializes `Theta_fCCZ4` only after the checkpoint branch has declined to return, while checkpoint data are repaired/projected without fresh Theta overwrite; spectroscopy shares initial and post-RHS projection across BSSN and opt-in fCCZ4, while the separate collision example keeps its BSSN-specific projection wiring.
- Role: public/scientific contract
- Deciding authority: [initial_data.py](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py), `register_CFunction_initial_data`; [ADM_Initial_Data_Reader__BSSN_Converter.py](../../../nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py), `register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`; [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py), formulation and projection registrations; [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py), BSSN initial-data and Method of Lines registrations
- Corroboration: [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, GCC 13.3.0, GNU Make 4.3; backend=BHaH OpenMP source registration; precision=double; GPU=not-run; restart=source-path-inspected only, no restart run; distributed=not-run; error_path=not-run; options=default BSSN generation/build, opt-in fCCZ4 generation/build and one-step startup, spectroscopy shared projection, collision path inspection only; date=08-28-2026`

`register_CFunction_initial_data` is the application-level initial-data
assembler. For built-in exact data it instantiates `InitialData_Cartesian` or
`InitialData_Spherical`, registers the exact ADM provider by
`register_CFunction_exact_ADM_ID_function`, then registers per-coordinate
ADM-to-BSSN readers through
`register_CFunctions_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`.
The generated `initial_data()` optionally attempts `read_checkpoint()` first;
on restart it applies inner boundary conditions and optional interpatch
interpolation before returning. Without restart it creates an `ID_persist`
struct, optionally populates it, loops over grids, optionally runs
`generalrfm_precompute`, calls
the selected initial-data reader conversion function, applies
outer-extrapolation plus inner boundary conditions, and optionally frees
persistent initial-data storage.

The ADM reader converter is a generated prefunc chain. `register_BHaH_defines_h`
adds `initial_data_struct` and `ID_persist_struct` to `BHaH_defines.h`.
`Cfunction_ADM_SphorCart_to_Cart` transforms ADM variables from the input
spherical, Cartesian, or GeneralRFM basis to Cartesian. `Cfunction_ADM_Cart_to_BSSN_Cart`
converts Cartesian ADM data to Cartesian BSSN fields. `Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm`
transforms those BSSN tensors/vectors to the destination reference-metric basis
and applies BSSN rescalings. `build_initial_data_conversion_loop` writes
`alpha`, `cf`, `trK`, `hDD`, `aDD`, `vetU`, `betU`, optional `T4UU`, and
opt-in fresh-data `Theta_fCCZ4` into MoL gridfunction arrays. YBS-MOM adds no
initial-data state.
`build_lambdaU_zeroing_block` initializes `lambdaU`,
`build_apply_inner_bcs_block` applies parity-sensitive inner boundary
conditions, and `Cfunction_initial_data_lambdaU_grid_interior` computes
`lambdaU` by finite differencing the initialized conformal metric.

`register_CFunction_diagnostic_gfs_set` bridges evolved GR state to diagnostics.
It registers diagnostic gridfunctions, builds parity metadata for them, and
generates `diagnostic_gfs_set(commondata, griddata, diagnostic_gfs)`. Runtime
flow is per grid: compute Ricci into `auxevol_gfs` or diagnostics, evaluate
constraints into diagnostics, optionally compute Psi4 and apply inner boundary
conditions to `DIAG_PSI4_RE/IM`, optionally apply inner boundary conditions to
constraint diagnostics before interpolation, then copy lapse, conformal factor,
and grid index into diagnostic channels. `register_CFunction_diagnostics_nearest`
and `register_CFunction_diagnostics_volume_integration` consume these
diagnostic buffers without owning their memory.

The nearest and volume diagnostic wiring is deliberately generic after
`diagnostic_gfs_set`. `diagnostics_nearest()` exposes user-editable `which_gfs`
arrays and dispatches to 0D, 1D, and 2D nearest samplers. Volume diagnostics
build recipes from diagnostic enum tokens and call
`diags_integration_execute_recipes`. The volume-element helper is
coordinate-specialized by `register_CFunction_sqrt_detgammahat_d3xx_volume_element`
when diagnostics are registered.

Basis transforms are registered through
`basis_transforms.register_all.register_CFunctions`. The two production modules
emit private per-coordinate single-point kernels with coordinate-system suffixes.
The public unsuffixed runtime dispatchers are emitted later by
`rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs`. The
rfm-to-Cartesian path reads rescaled BSSN storage, reconstructs barred
quantities, transforms `gammabarDD`, `AbarDD`, beta, B, and Lambdabar to
Cartesian, and writes Cartesian-basis storage. The Cartesian-to-rfm path reads
Cartesian storage, transforms tensors/vectors into the destination basis,
subtracts the destination reference metric where needed, rescales by `ReDD` and
`ReU`, and writes native BSSN storage.

TwoPunctures is wired as an external compact-object initial-data family.
`TwoPunctures_lib.register_C_functions_explicit` registers
`initialize_ID_persist_struct`, `TP_CoordTransf`, `TP_Equations`,
`TP_FuncAndJacobian`, `TP_Newton`, `TP_Interp`, `TP_solve`, and
`TP_utilities`. `ID_persist_str` registers commondata controls for binary
description, spectral resolution, and optional bare masses, then contributes
the persistent spectral-solve fields. `initialize_ID_persist_struct` sets
defaults, copies or derives binary masses, separation, momenta, spin, center
offset, spectral grid sizes, and orientation. The explicit orientation choices
are `native_cartesian_xy_plane` and `legacy_swap_xz`; the older
`register_C_functions(enable_xy_plane=...)` shim maps the boolean interface to
those names. The shared initial-lapse default remains `psi^n`. `TP_Interp`
also recognizes `W` and assigns `1 / (psi1 / static_psi)^2`, so the spectral
correction `U` in `psi1` participates in that initial lapse. The top-level
standalone and superB spectroscopy generators explicitly select `W` for SSL;
callers without an override, including paper-reproduction examples, retain
`psi^n`.

Claim evidence:
- Claim: TwoPunctures retains its shared `psi^n` initial-lapse default, supports an optional `W` selector evaluated from the corrected total conformal factor, and the top-level standalone and superB spectroscopy generators explicitly select `W` for SSL while paper-reproduction examples continue to use the shared initializer without a lapse override; generation and compilation alone do not establish evolution behavior or scientific accuracy.
- Role: descriptive behavior
- Deciding authority: [ID_persist_struct.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/ID_persist_struct.py), `register_CFunction_initialize_ID_persist_struct`; [TP_interp.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/TP_interp.py), `register_CFunction_TP_Interp`; [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py) and [superB_blackhole_spectroscopy.py](../../../nrpy/examples/superB_blackhole_spectroscopy.py), `populate_ID_persist_struct_str`
- Corroboration: [superB_blackhole_spectroscopy_8Mseparation.py](../../../nrpy/examples/superb_paper2/superB_blackhole_spectroscopy_8Mseparation.py) and [superB_blackhole_spectroscopy_last_orbit.py](../../../nrpy/examples/superb_paper2/superB_blackhole_spectroscopy_last_orbit.py), unchanged shared-initializer calls without a lapse override
- Validation: `inspected=pass; generated=pass; built=pass; run=not-run; result_checked=not-run`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, GCC 13.3.0, GNU Make 4.3; backend=BHaH OpenMP; precision=double; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=top-level standalone W generation/build, top-level superB W generation, paper examples inspected unchanged with shared psi^n default; date=09-04-2026`

If TwoPunctures/TOVola setup detail grows beyond routing and
dataflow, split it into a future `compact-object-initial-data.md` leaf.

TOVola is the corresponding single-star initial-data path. `TOVola.ID_persist_str`
registers central density, polytropic EOS constants, ODE controls, and
interpolation stencil sizes, then contributes radial-table pointers and count
fields to `ID_persist_struct`. `register_CFunction_TOVola_solve` registers the
GSL ODE integration driver, stores the solved radial data in the persistent
arrays, and frees temporary solve storage. `register_CFunction_TOVola_interp`
registers the pointwise interpolation provider that maps Cartesian points to
isotropic radius, interpolates the radial table, and fills ADM-like
`initial_data_struct` fields including lapse, spherical spatial metric,
zero extrinsic curvature, and stress-energy components. The GRHD TOV example
uses this by registering `TOVola_interp`, `TOVola_solve`, and a spherical
ADM-to-BSSN reader whose `ID_persist_struct_str` comes from TOVola.

Psi4 wiring has two layers. `psi4.register_CFunction_psi4` temporarily forces
OpenMP generation, requires the grid origin to be zero, calls
`generate_CFunction_psi4_tetrad` and
`generate_CFunction_psi4_metric_deriv_quantities` to create local helper
kernels, then loops over interior points to write `DIAG_PSI4_REGF` and
`DIAG_PSI4_IMGF`. The tetrad helper uses `Psi4Tetrads` to compute
`mre4U`, `mim4U`, and `n4U` from metric gridfunction values; the derivative
helper computes `gammaDDdDD`, `GammaUDD`, and `KDDdD` arrays used by the
symbolic `Psi4` expression. `psi4_spinweightm2_decomposition` then interpolates
diagnostic Psi4 data from grid 0 onto spherical extraction shells, calls
`spin_weight_minus2_sph_harmonics` for each `(l,m)`, numerically integrates
over the shell, and appends `R_ext * psi4_{l,m}` time series to radius- and
mode-tagged text files.

`spin_weight_minus2_sph_harmonics` is registered in the BHaH special-functions
branch. It registers `swm2sh_maximum_l_mode_to_compute` in commondata and emits
a switch over every `m` in the inclusive source range `[-l, +l]` for every
generated `l`. Out-of-range requests print an error and exit, so the
decomposition generator must register enough modes for the requested
extraction. This inclusive infrastructure generator is distinct from the
equation helper's narrower current test sweep; the equation-side validation gap
does not remove `m=+l` cases from emitted BHaH C.

## Sources

- [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`
- [dsmin_gf.py](../../../nrpy/infrastructures/BHaH/general_relativity/dsmin_gf.py) - `register_CFunction_dsmin_auxevol_gridfunction`
- [Ricci_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/Ricci_eval.py) - `register_CFunction_Ricci_eval`
- [BSSN_quantities.py](../../../nrpy/equations/general_relativity/BSSN_quantities.py) - `BSSNQuantities`, mathematical derivative construction
- [hDDdD_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/hDDdD_eval.py) - `register_CFunction_hDDdD_eval`, `register_hDDdD_gridfunctions`
- [cfdD_alphadD_vetUdD_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/cfdD_alphadD_vetUdD_eval.py) - `register_CFunction_cfdD_alphadD_vetUdD_eval`, `cfdD_alphadD_vetUdD_gridfunction_expressions`, `register_cfdD_alphadD_vetUdD_gridfunctions`
- [constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py) - `register_CFunction_constraints_eval`
- [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- [initial_data.py](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py) - `register_CFunction_initial_data`
- [ADM_Initial_Data_Reader__BSSN_Converter.py](../../../nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py) - `register_CFunction_exact_ADM_ID_function`, `register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`, `Cfunction_ADM_SphorCart_to_Cart`, `Cfunction_ADM_Cart_to_BSSN_Cart`, `Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm`, `Cfunction_initial_data_lambdaU_grid_interior`
- [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostic_gfs_set.py) - `register_CFunction_diagnostic_gfs_set`
- [diagnostics_nearest.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostics_nearest.py) - `register_CFunction_diagnostics_nearest`
- [diagnostics_volume_integration.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostics_volume_integration.py) - `register_CFunction_diagnostics_volume_integration`
- [register_all.py](../../../nrpy/infrastructures/BHaH/general_relativity/basis_transforms/register_all.py) - `register_CFunctions`
- [basis_transform_BSSN_rfm_to_Cartesian_single_point.py](../../../nrpy/infrastructures/BHaH/general_relativity/basis_transforms/basis_transform_BSSN_rfm_to_Cartesian_single_point.py) - `register_CFunction_basis_transform_BSSN_rfm_to_Cartesian_single_point`
- [basis_transform_BSSN_Cartesian_to_rfm_single_point.py](../../../nrpy/infrastructures/BHaH/general_relativity/basis_transforms/basis_transform_BSSN_Cartesian_to_rfm_single_point.py) - `register_CFunction_basis_transform_BSSN_Cartesian_to_rfm_single_point`
- [TwoPunctures_lib.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/TwoPunctures_lib.py) - `register_C_functions_explicit`, `register_C_functions`
- [ID_persist_struct.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/ID_persist_struct.py) - `ID_persist_str`, `register_CFunction_initialize_ID_persist_struct`
- [TP_interp.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/TP_interp.py) - `register_CFunction_TP_Interp`
- [TOVola/ID_persist_struct.py](../../../nrpy/infrastructures/BHaH/general_relativity/TOVola/ID_persist_struct.py) - `ID_persist_str`
- [TOVola_solve.py](../../../nrpy/infrastructures/BHaH/general_relativity/TOVola/TOVola_solve.py) - `register_CFunction_TOVola_solve`
- [TOVola_interp.py](../../../nrpy/infrastructures/BHaH/general_relativity/TOVola/TOVola_interp.py) - `register_CFunction_TOVola_interp`
- [psi4.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4/psi4.py) - `register_CFunction_psi4`
- [compute_psi4_metric_deriv.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4/compute_psi4_metric_deriv.py) - `generate_CFunction_psi4_metric_deriv_quantities`
- [compute_psi4_tetrad.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4/compute_psi4_tetrad.py) - `generate_CFunction_psi4_tetrad`
- [psi4_spinweightm2_decomposition.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4_spinweightm2_decomposition.py) - `register_CFunction_psi4_spinweightm2_decomposition`, `lowlevel_decompose_psi4_into_swm2_modes`
- [spin_weight_minus2_spherical_harmonics.py](../../../nrpy/infrastructures/BHaH/special_functions/spin_weight_minus2_spherical_harmonics.py) - `register_CFunction_spin_weight_minus2_sph_harmonics`
- [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py) - `BHaH.general_relativity.rhs_eval.register_CFunction_rhs_eval`, `BHaH.general_relativity.basis_transforms.register_all.register_CFunctions`
- [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py) - `--fccz4`, formulation-selecting RHS/initial-data registration, shared projection hooks
- [superB_two_blackholes_collide.py](../../../nrpy/examples/superB_two_blackholes_collide.py) - explicit default-disabled YBS flag forwarded to shared RHS registration
- [superB_blackhole_spectroscopy.py](../../../nrpy/examples/superB_blackhole_spectroscopy.py) - explicit default-disabled YBS flag forwarded to shared RHS registration
- [groovy_TOV_BSSN.py](../../../nrpy/examples/groovy_TOV_BSSN.py) - `BHaH.general_relativity.TOVola.TOVola_interp.register_CFunction_TOVola_interp`, `BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter.register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`
- [SOURCES.md](../../../raw/SOURCES.md) - `infrastructure-modules-and-embedded-headers`

## See Also

- Parent: [BHaH](index.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- Depends on: [YBS-MOM Timestep-Scaled Momentum Adjustment](../../equations/general-relativity/ybs-momentum-damping.md)
- Depends on: [Initial Data](../../equations/general-relativity/initial-data.md)
- Depends on: [Metric Conversions And Matter](../../equations/general-relativity/metric-conversions-and-matter.md)
- Depends on: [Psi4 And Tetrads](../../equations/general-relativity/psi4-and-tetrads.md)
- Depends on: [Geometry And Special-Function Support](../../equations/geometry-and-special-function-support.md)
- See also: [Diagnostics Output And Checkpointing](diagnostics-output-and-checkpointing.md)
