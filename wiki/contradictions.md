# Contradictions

> Register for contested and stale KB claims. Plain Markdown only. · Status: confirmed

Known contested/stale claims are tracked below. A row records
source-side truth and containment; it does not imply that NRPy source was fixed.

## Rules

- Follow the contradiction contract in [Schema](SCHEMA.md#contradiction-contract).
- Open a row when registered sources disagree, when a living source supersedes
  a KB claim, or when a page is marked `contested` or `stale`.
- Update affected-page links and exact inline markers together.
- Resolve a row only after its resolution test passes and all affected pages,
  reverse dependents, aliases, typed neighbors, and targeted wiki hits are
  reconciled.
- Look up one ID with `rg -n -F 'CONTR-0014' wiki/contradictions.md`.
- Resolved rows are removed once the affected pages record the outcome and no
  page links to the row; git history keeps them. Retired IDs: CONTR-0004,
  CONTR-0006, CONTR-0008 to CONTR-0011. New IDs continue after the highest ID
  used.

## Register

| ID | Claim | Claim status | Source A | Source B | Authority decision | Affected pages | Page-status rationale | Owner/trigger | Resolution test | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CONTR-0001 | `manga_bhah_lib` tells users to run `./bhah_lib` and says `bhah_lib.par` exists, but its assembly path requests a library target and does not call the parser/default-parfile writer. | stale | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) final two `print()` calls | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) commented parser registration and `output_CFunctions_function_prototypes_and_construct_Makefile(create_lib=True)` call; [`Makefile_helpers.py`](../nrpy/infrastructures/BHaH/Makefile_helpers.py) `output_CFunctions_function_prototypes_and_construct_Makefile` | For descriptive current behavior, precise registrations and writer calls decide. They emit `libbhah_lib.so` on Linux or `libbhah_lib.dylib` on Darwin, not an executable `bhah_lib`; no invoked writer emits `bhah_lib.par`. Both generic final prints are stale. | [BHaH Lifecycle And Project Assembly](infrastructures/bhah/lifecycle-and-project-assembly.md); [Matter TOV Workflows](examples/matter-tov-workflows.md) | Both pages remain `confirmed`: these are bounded final-message claims, corrected inline with explicit non-guarantees; their principal lifecycle/workflow answers remain supported. | BHaH/example owners; trigger when `manga_bhah_lib.py` final messaging, parser/parfile registration, or Makefile arguments change. | In an isolated output root, run the generator, capture stdout, inspect emitted files and the Makefile `all` target, and build. Pass only if stdout gives library-consumer guidance, never instructs executing `./bhah_lib`, mentions a parfile only if one is emitted, and `make` produces the platform library named by the Makefile. | Configured CI does not cover this resolution test. |
| CONTR-0002 | `sebobv1_jax` emits `Commondata(..., a_f=a_f)`, but its batch field registration supplies fourteen names/descriptions and thirteen dtypes/defaults; `zip()` therefore omits `a_f` from the generated dataclass. | resolved | [`SEOBNRv5_aligned_spin_coefficients.py`](../nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py) `register_PyFunction_SEOBNRv5_aligned_spin_coefficients` emitted return | [`sebobv1_jax.py`](../nrpy/examples/sebobv1_jax.py) `register_commondata_params` call and [`commondata.py`](../nrpy/infrastructures/JAX/commondata.py) `register_commondata_params` / `generate_commondata_dataclass` | Descriptive behavior is decided by all three code paths together. `register_commondata_params` now raises `ValueError` when its four lists differ in length, `sebobv1_jax` supplies fourteen entries in each list, and the generated `Commondata` therefore has every keyword, including `a_f`, that the emitted function passes. | [SEBOBv1 JAX Workflow](infrastructures/jax/sebobv1-jax-workflow.md); [Commondata And PyFunction Registry](infrastructures/jax/commondata-and-pyfunction-registry.md); [Waveform JAX PN Generators](examples/waveform-jax-pn-generators.md); [Example Generator Catalog](examples/example-generator-catalog.md) | The inline markers are removed. `SEBOBv1 JAX Workflow` returns to `confirmed` because its coefficient interface now matches the registered fields; the other pages already were `confirmed` and now state the corrected registration. | JAX/example owners; trigger when batch registration loses its length check, the four lists change, or the emitted return signature changes. | In an isolated output root, generate the project; assert batch metadata passes equal-length validation; import the generated coefficient function; call it with a valid finite input set; assert every returned keyword is a generated `Commondata` field and every returned numerical field is finite; if `a_f` remains in the intended return contract, assert that `a_f` is present and finite; then run the generated smoke test. All applicable checks must pass. | Configured CI installs the generated package and calls the coefficient function once, but asserts no returned field or value, so this resolution test is not run in CI. Generated code became runnable only together with the `py_codegen()` change that prints right-hand sides without `assign_to` and `NRPyJaxPrinter._print_Max`; see [Python Codegen](core/python-codegen.md) and [CSE And Printer Support](core/helpers/cse-and-printer-support.md). |
| CONTR-0003 | Pylint policy distinguishes grandfathered existing tracked handwritten files from newly added handwritten files, but both configs use `fail-under=10`, the local wrapper uses a blanket `9.91` floor, and the workflow uses a blanket `9.5` floor; none enforces the tiered rule. | stale | [`coding_style.md`](../coding_style.md) `## Static Analysis Configuration` | [`.pylintrc`](../.pylintrc) and [`.pylintrc_python36`](../.pylintrc_python36) `[MASTER]`; [`single_file_static_analysis.sh`](../.github/single_file_static_analysis.sh) `run_test_step`; [`main.yml`](../.github/workflows/main.yml) `static-analysis` | Existing tracked handwritten files retain their pre-change score without regression, including legacy scores at or below `9.5`; newly added handwritten files require `10.00/10.00`. Enforcement must classify the path from the base revision, compare legacy base and proposed ratings under the same tool/config, require exact `10.00` for new paths, and reject missing or unparseable ratings. | [Static Analysis](validation/static-analysis.md) | The page remains `confirmed`: the enforcement mismatch is one bounded, directly sourced Detail claim with an inline marker and explicit non-guarantee; the principal command and configuration answer remains reliable. | Validation and contributor-guidance owners; trigger when the tiered Pylint rule, either config, local wrapper classification/rating handling, or workflow classification/rating handling changes. | Deterministically inspect contributor guidance, both configs, local wrapper, and workflow. Pass only if documentation states the tiered rule and every enforcement path classifies paths from a base revision, rejects a lower legacy rating, requires exactly `10.00` for a new handwritten file, and rejects missing or unparseable ratings; generated trusted data remains exempt. | Configs and wrappers do not implement the required classification. |
| CONTR-0005 | Mewes et al.'s prose around its conformal-factor alternatives reverses the signs implied by its displayed evolution equations and used by NRPy's established BSSN variables. | resolved | [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (33)-(34) and surrounding prose | [`BSSN_quantities.py`](../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__` | NRPy retains its established `W=exp(-2 phi)` and `chi=exp(-4 phi)` definitions, whose chain rules agree with the displayed evolution equations; sign-reversed prose is not adopted. | [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md) | Page remains `confirmed`: definitions and authority choice are explicit and owner validation covers all three conformal-factor options. | BSSN/fCCZ4 equation owners; trigger when conformal-factor definitions, source crosswalk, or option validation changes. | Run fCCZ4 constraint/RHS owner validation across `W`, `phi`, and `chi`; verify exact chain-rule definitions and trusted expressions. | This resolves a source-internal prose/equation mismatch without changing NRPy's established variables. |
| CONTR-0007 | The frozen commissioned YBS-MOM source proposes an auxiliary hyperbolic-relaxation field, but the corrected implementation contract requires the Yo--Lin--Cao term directly with CAHD-style timestep scaling and no new evolved state. | resolved | [preserved YBS-MOM specification](../raw/source-docs/ybs-momentum-damping-spec.md), `Proposed system` and auxiliary-field analysis | [`BSSN_RHSs.py`](../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch; commissioned user clarification | The user clarification controls the commissioned design: apply the direct covariant STF momentum gradient with `C_YBS_mom*CFL_FACTOR*DSMINGF`, exactly analogous to CAHD's timestep-scaled parabolic coefficient. Retain the frozen source as historical provenance; do not add `qD` or any other evolved cleaner field. A later maintainer decision multiplies the coefficient by the conformal factor `W`; the direct, stateless formulation is unchanged. | [YBS-MOM Timestep-Scaled Momentum Adjustment](equations/general-relativity/ybs-momentum-damping.md); [BSSN Family](equations/general-relativity/bssn-family.md); [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md); [BHaH GR Application Wiring](infrastructures/bhah/gr-application-wiring.md) | Pages remain `confirmed`: code and compiled pages now agree on the direct no-state contract; the superseded proposal is explicitly bounded. | BSSN/fCCZ4/BHaH owners; trigger when YBS-MOM variables, coefficient, momentum residual, STF projection, shared spacing, defaults, or source authority changes. | Verify no `qD`/`q_rhsD` gridfunction or output exists; compare the lower residual with canonical `BSSNconstraints.MU`; run all existing jointly YBS-enabled trusted dictionaries; generate ordinary and fisheye OpenMP/CUDA spacing paths; run static analysis and KB lint. | Resolution rejects auxiliary-state hyperbolization. The term remains parabolic in principal character; CAHD-style local timestep scaling satisfies the intended CFL compatibility. |
| CONTR-0012 | Python codegen and ordinary scalar C codegen expose substitution options, but a surviving postprocessed CSE temporary makes a nonempty substitution dictionary raise `TypeError`. | contested | [py_codegen.py](../nrpy/py_codegen.py), `PyCodeGen` and `py_codegen` option and returned-code descriptions; [c_codegen.py](../nrpy/c_codegen.py), `CCodeGen` and `c_codegen` substitution option | [py_codegen.py](../nrpy/py_codegen.py) and [c_codegen.py](../nrpy/c_codegen.py), assignment to `common_subexpression[1]`; [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `cse_postprocess` tuple entries | The implementation decides current behavior; tuple assignment prevents code return for these combinations. The KB documents the failure without treating it as successful substitution. | [Python Codegen](core/python-codegen.md); [C Codegen](core/c-codegen.md) | Both pages are `contested` because a public option combination fails. | C and Python codegen sources; recheck when substitution handling or CSE replacement representation changes. | For Python and ordinary non-FD, non-SIMD C generation, retain a common subexpression and supply a nonempty substitution dictionary with matching and nonmatching keys; require code return and agreement with the intended substituted expressions. Retain no-CSE and no-temporary cases. | This claim does not extend the scalar C failure to SIMD or FD paths. |
| CONTR-0013 | The public Python-codegen description maps input expressions to output names, but an output name could overwrite a CSE temporary or input symbol required by a later output. | resolved | [py_codegen.py](../nrpy/py_codegen.py), `py_codegen` input/output description | [py_codegen.py](../nrpy/py_codegen.py), numbered temporary construction and sequential output assignment | The implementation decides. `py_codegen()` now parses each printed output right-hand side and raises `ValueError` when an earlier output name is a Python name that a later one reads, with CSE enabled or disabled. This covers CSE temporaries, input symbols, identifiers renamed by `NRPyJaxPrinter` such as `lambda_`, and printer-emitted names such as `jnp`, so it never returns assignments that read an overwritten name. Names that no later right-hand side reads remain accepted. | [Python Codegen](core/python-codegen.md) | The inline marker is removed. The page stays `contested` only because of CONTR-0012. | Python-codegen source; recheck when temporary naming, output assignment order, `_check_no_read_after_overwrite`, or `NRPyJaxPrinter` identifier renaming changes. | Generate `[x**3 + sin(x**3), x**3]` with output names `["tmp0", "b"]`: with CSE enabled, require `ValueError`, because the second output reads the temporary `tmp0`; with CSE disabled, no temporary exists, so require returned assignments that evaluate to `x**3 + sin(x**3)` and `x**3` (at `x = 2`, `8 + sin(8)` and `8`). Generate `[x + 1, 2*x]` with `["x", "y"]` with CSE enabled and disabled; require `ValueError` in each case. Generate `[1, 2*lambda]` with `["lambda_", "b"]` with CSE enabled and disabled; require `ValueError`. Generate `[2*x, x + 1]` with `["y", "x"]`; require the unchanged assignments `y = 2*x` and `x = x + 1`. | Resolved by rejection, not by renaming temporaries; generated variable names are unchanged. The `py_codegen` doctests cover the resolution test. The check reads names from the printed code, so printer renaming of reserved words is covered. |
| CONTR-0014 | CSE sorting can preserve input-list order for tied structural keys, and unchecked generated constant or temporary names can capture input symbols. | contested | [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), documented helper behavior | [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `cse_preprocess, sort_cse_output_deterministically` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [CSE And Printer Support](core/helpers/cse-and-printer-support.md); [C Codegen](core/c-codegen.md) | Pages are `contested`: ordering and expression preservation require qualifications. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Permute replacements for `x**y` and `y**x` and require equal renamed output. Restore preprocessed constants in `2*x + _Integer_2`, and expand renamed output for `([(a0,x**3)], [a0+tmp0])`; require the original expressions or explicit collision rejection. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0015 | SIMD symbolic rewrites can omit the reciprocal numerator constant or emit Pow/Cbrt/Log/Sign names missing from the header. | contested | [simd.py](../nrpy/helpers/simd.py), `expr_convert_to_simd_intrins` | [simd_intrinsics.h](../nrpy/helpers/simd_intrinsics.h), instruction-set definitions | The implementation decides current behavior; missing numerator registration and definitions prevent a general supported-expression guarantee. | [SIMD And Intrinsic Support](core/helpers/simd-and-intrinsic-support.md) | Page is `contested`: missing constants and definitions affect generated calculations. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Verify declared constants for negative half-powers and compile supported Pow/Cbrt/Log/Sign expressions. | Including simd_intrinsics.h alone does not supply the missing SIMD names. |
| CONTR-0016 | The symbol-name helper claims SymPy free-symbol equivalence but removes free occurrences of bound-symbol names in limits or substitution points. | contested | [expression_utils.py](../nrpy/helpers/expression_utils.py), documented helper behavior | [expression_utils.py](../nrpy/helpers/expression_utils.py), `_get_free_symbol_names_cached, get_params_commondata_symbols_from_expr_list` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [Symbolic Expression Utilities](core/helpers/symbolic-expression-utilities.md) | Page is `contested`: symbol discovery, rational conversion, and product behavior each have source discrepancies. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Compare names with SymPy free_symbols for Integral, Lambda, and Subs, including `Integral(x,(x,0,x))`; verify parameter grouping retains registered free names. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0017 | f2r claims high-precision float conversion but appends zeros after scientific-notation exponents, changing values or creating invalid input; padding cannot recover lost float precision. | contested | [float_to_rational.py](../nrpy/helpers/float_to_rational.py), documented helper behavior | [float_to_rational.py](../nrpy/helpers/float_to_rational.py), `f2r` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [Symbolic Expression Utilities](core/helpers/symbolic-expression-utilities.md) | Page is `contested`: multiple helper guarantees require qualifications. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Require scientific and fixed-point strings for the same decimal value to produce equal rationals; cover `1e-7`, `1.2e-7`, and explicit padding values. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0018 | product describes Cartesian tuples but flattens iterable-valued elements, consumes one-shot inner iterators, and mishandles single scalar-valued inputs and repeat=0. | contested | [functional.py](../nrpy/helpers/functional.py), documented helper behavior | [functional.py](../nrpy/helpers/functional.py), `product` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [Symbolic Expression Utilities](core/helpers/symbolic-expression-utilities.md) | Page is `contested`: multiple helper guarantees require qualifications. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Compare with itertools.product for reusable and one-shot inputs, iterable-valued elements, single integer-valued inputs, and repeat=0. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0019 | Kernel call generation corrupts empty-argument calls; the CUDA wrapper ignores cfunc_type and may leave rewritten parameter accesses without a declared streamid. | contested | [utilities.py](../nrpy/helpers/parallelization/utilities.py), documented helper behavior | [gpu_kernel.py](../nrpy/helpers/parallelization/gpu_kernel.py), `GPUKernel.c_function_call and GPUKernel.__init__` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [Loop Kernel And Device Helpers](core/helpers/loop-kernel-and-device-helpers.md) | Page is `contested`: public wrapper inputs can produce invalid calls or a different return type. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Generate and compile an empty-argument host call, a CUDA helper with a requested non-void return type, and custom CUDA launch settings with parameter reads; require a valid call, the requested return type, and a declared streamid. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0020 | Coordinate-bound helper docstrings promise Boolean comparisons, but shifted quotients yield fractional values, shifted thresholds, and singular denominators; the SEOBNR ringdown frequency inherits those weights. | contested | [Min_Max_and_Piecewise_Expressions.py](../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py), coordinate-helper docstrings | [Min_Max_and_Piecewise_Expressions.py](../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py), coordinate-helper return expressions | Return expressions decide current behavior; exact comparison semantics are not guaranteed. | [GRHD](equations/grhd.md); [Precessing Rotations And Ringdown](equations/seobnr/precessing-rotations-and-ringdown.md) | GRHD remains `confirmed`: one bounded helper claim is qualified inline. The SEOBNR page is `contested` because the frequency, Euler angles, and QNM shifts inherit the discrepancy. | GRHD comparison helpers and SEOBNR ringdown constructor; recheck on changes to offsets, formulas, or documented comparison semantics. | Interpret nrpyAbs as absolute value and compare with true predicates on both sides of zero and at plus/minus TINYDOUBLE; require finite 0-or-1 results matching the predicates or an explicitly restricted interface. Verify the corresponding ringdown frequency weights and downstream angle/QNM expressions. | Trusted samples replacing nrpyAbs with sin, cos, or exp do not decide these semantics. |
| CONTR-0021 | Wave constructors can leave required gridfunctions or CodeParameters unregistered when existing names trigger incomplete or incorrect guards. | contested | [WaveEquation_RHSs.py](../nrpy/equations/wave_equation/WaveEquation_RHSs.py), constructor registration comments; [WaveEquation_Solutions_InitialData.py](../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py), constructor registration comments | [WaveEquation_RHSs.py](../nrpy/equations/wave_equation/WaveEquation_RHSs.py), [WaveEquationCurvilinear_RHSs.py](../nrpy/equations/wave_equation/WaveEquationCurvilinear_RHSs.py), and [WaveEquation_Solutions_InitialData.py](../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py), constructors; [params.py](../nrpy/params.py), CodeParameter registration | Actual registry checks decide current behavior; callers must provide compatible entries when reusing state. | [Wave Equation](equations/wave-equation.md) | Page is `contested`: independent gridfunction and parameter-registration conditions invalidate unconditional setup claims. | Wave equation constructors; recheck when registry guards or registration validation changes. | In separate isolated processes, pre-register uu without vv, or with incompatible metadata, and pre-register NRPy parameters named wavespeed or time; require compatible gridfunction and CodeParameter entries or explicit rejection. Keep a clean-registry case. | Curvilinear RHS wavespeed registration is unconditional; its uu-only gridfunction guard still applies. |
| CONTR-0022 | A zero Lorentz boost should be the identity, but direct zero-velocity construction produces undefined spatial entries. | contested | [LorentzBoost.py](../nrpy/equations/general_relativity/LorentzBoost.py), documented Lorentz-transformation purpose | [LorentzBoost.py](../nrpy/equations/general_relativity/LorentzBoost.py), `_compute_LorentzMatrix` spatial block | The implemented division by squared speed decides current behavior; zero input is not handled by its identity limit. | [Metric Conversions And Matter](equations/general-relativity/metric-conversions-and-matter.md) | Page remains `confirmed`: this is one bounded helper defect, qualified inline; other metric and matter conversions remain valid. | LorentzBoost constructor; recheck when speed-domain handling or the spatial matrix formula changes. | Construct the direct zero boost and its inverse; require finite identity matrices and unchanged transformed tensors, retaining nonzero subluminal cases. | This KB correction documents the source defect without changing Python. |
| CONTR-0023 | The geodesic azimuthal-symmetry spatial interpolation rotates each reflected stencil record about z by +/-pi before accumulation, although the record already holds the Cartesian-basis tensor at the reflected node's own position. | contested | [azimuthal_symmetry_spatial_lagrange_interpolation.py](../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/azimuthal_symmetry_spatial_lagrange_interpolation.py), helper description of rotating reflected nodes to the reference plane | [output_raytracing_data.py](../nrpy/infrastructures/BHaH/diagnostics/output_raytracing_data.py), module description of stored Cartesian tensor components; [azimuthal_symmetry_spatial_lagrange_interpolation.py](../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/azimuthal_symmetry_spatial_lagrange_interpolation.py), `azimuthal_symmetry_spatial_lagrange_map_extended_node` reflection rules | The generated helper decides current behavior; because a reflected record is the Cartesian tensor at the raw node's position, the per-node rotation is not adopted as correct. | [Geodesics And Raytracing Runtime](infrastructures/bhah/geodesics-and-raytracing-runtime.md) | Page remains `confirmed`: one bounded defect in an optional interpolation helper, qualified inline; the other geodesic and export claims remain valid. | Geodesic interpolation helper; recheck when the reflection map, the per-node rotation, or the record basis changes. | On `Spherical` data with theta in `[0, pi]`, phi in `[-pi, pi]`, and `Nphi = 2`, interpolate a nonconstant z-axisymmetric Cartesian-basis `g4DD` and `Gamma4UDD` field at targets whose (r, theta) stencil crosses theta = 0, theta = pi, or r = 0; require the convergence rate of targets whose stencil stays in the stored interior. | Output-sample defect; evolved fields are unaffected. |

### CONTR-0001

The resolution test deliberately requires generation, artifact inspection, and
a build because source inspection alone cannot establish the user-facing
library contract.

### CONTR-0002

The resolution test deliberately asserts the returned fields and their finite
values because the configured CI call checks only that the generated function
returns without error.

The mismatch no longer occurs: `sebobv1_jax` now supplies a dtype and default
for `a_f`, and `register_commondata_params` rejects lists of unequal length
instead of letting `zip()` drop trailing fields. Running the generated function
also required `py_codegen()` to print each right-hand side without passing the
output name as `assign_to`, which had dropped negative powers such as `m2/m1`,
and `NRPyJaxPrinter` to print the `Max` in the final mass without `functools`;
see [Python Codegen](core/python-codegen.md) and [CSE And Printer
Support](core/helpers/cse-and-printer-support.md). The large integers in this
function appear only as numerators and denominators of rational coefficients,
which Python divides to floats, so this function does not depend on the printer
rules for large integers and rationals.

### CONTR-0003

Existing tracked handwritten files retain their pre-change Pylint score without
regression, including legacy scores at or below `9.5`; newly added handwritten
files require `10.00/10.00`. Enforcement must classify paths from the base
revision and reject missing or unparseable ratings.

Evidence (normative rule): registered [`coding_style.md`](../coding_style.md), `## Static Analysis Configuration`, as updated by the commissioned policy decision; corroboration: [Static Analysis](validation/static-analysis.md) and [Code Test Policy](validation/code-test-policy.md) apply the rule; registered [`single_file_static_analysis.sh`](../.github/single_file_static_analysis.sh), `run_test_step`, and [`main.yml`](../.github/workflows/main.yml), `static-analysis`, identify enforcement surfaces but do not implement the distinction

The configs and wrappers use blanket floors. None classifies paths from a base
revision or checks legacy score regression.

### CONTR-0005

NRPy uses `W=exp(-2 phi)` and `chi=exp(-4 phi)`. These definitions agree with
Mewes et al.'s displayed conformal-factor equations; the neighboring
sign-reversed prose is not part of NRPy's scientific contract.

Evidence (public/scientific contract): [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (33)-(34); [`BSSN_quantities.py`](../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__`; corroboration: [`fCCZ4_constraints.py`](../nrpy/equations/general_relativity/fCCZ4_constraints.py), module `__main__`, exercises `W`, `phi`, and `chi` across Cartesian and SinhSpherical owner cases

### CONTR-0007

The frozen commissioned source interpreted “hyperbolize” as introducing an
auxiliary relaxation covector. The later explicit user clarification instead
defines the intended analogy to CAHD: keep the direct parabolic
momentum-gradient adjustment and multiply its coefficient by the local
timestep scale. Current YBS-MOM therefore applies the direct covariant
symmetric trace-free gradient of the matter-complete lower momentum residual,
multiplied by `C_YBS_mom*CFL_FACTOR*DSMINGF*W`; the weight `W` is a later
maintainer decision. The frozen auxiliary-relaxation proposal is superseded,
and no cleaner state is permitted. Current code therefore adds no evolved field
and leaves initial data, restart, boundary, AMR, and dissipation state routes
unchanged.

Evidence (public/scientific contract): commissioned user clarification composed with [`BSSN_RHSs.py`](../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch; corroboration: [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (56), and [Etienne, Phys. Rev. D 110, 064045](https://doi.org/10.1103/PhysRevD.110.064045), Eq. (26)

### CONTR-0012

The postprocessed CSE replacements are tuples, but a nonempty
`postproc_substitution_dict` makes `py_codegen()` assign a new right-hand side
into `common_subexpression[1]`. The ordinary non-FD, non-SIMD CSE loop in
`c_codegen()` performs the same assignment. A surviving temporary therefore causes
`TypeError` even when no substitution key matches. The no-CSE path and the
final-expression path without surviving temporaries do not perform that tuple
assignment. This is an implementation limitation, not evidence of successful
support for every combination of options.

Evidence (descriptive behavior): [py_codegen.py](../nrpy/py_codegen.py), `py_codegen`, and [c_codegen.py](../nrpy/c_codegen.py), `c_codegen`, scalar common-subexpression substitution loops; corroboration: [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `cse_postprocess`, constructs the returned replacement entries as tuples

### CONTR-0013

Resolved. `py_codegen()` emits output assignments sequentially into one Python
namespace. Previously, for `[x**3 + sin(x**3), x**3]` and targets
`["tmp0", "b"]`, the first output assignment replaced the temporary that the
second still read, and an output could likewise overwrite an input symbol read
later. `_check_no_read_after_overwrite` now raises `ValueError` whenever a
later printed right-hand side reads an earlier output name, on the CSE and
no-CSE paths, before any code is returned. The test is whether a later
right-hand side reads the name, not how the output is spelled: with CSE
disabled, the `tmp0` example creates no temporary, so its second output reads
only `x` and the request is accepted with correct values. The input-symbol
example `[x + 1, 2*x]` with targets `["x", "y"]` is rejected on both paths.
Requests in which no later right-hand side reads an output name emit the same
assignments as before; the generated SEBOBv1 coefficient module is unchanged. The fix rejects
colliding names rather than renaming temporaries, so generated variable names
do not change.

The check parses each printed output right-hand side with Python's `ast`
module, so it compares output names with the names the generated code reads,
not with SymPy symbol names. A symbol named `lambda` prints as `lambda_`, so
`[sp.Integer(1), 2*sp.Symbol("lambda")]` with output names `["lambda_", "b"]`
is rejected with CSE enabled or disabled. Printer-emitted names such as `jnp`
are covered the same way.

Evidence (descriptive behavior): [py_codegen.py](../nrpy/py_codegen.py), `_check_no_read_after_overwrite` and its calls in `py_codegen`, and the module-level `printer` (the `py_codegen` doctests in the deciding module exercise the rejected and accepted cases)

### CONTR-0014

Evidence (descriptive behavior): [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `cse_preprocess, sort_cse_output_deterministically`

### CONTR-0015

Evidence (descriptive behavior): [simd.py](../nrpy/helpers/simd.py), `expr_convert_to_simd_intrins`; [simd_intrinsics.h](../nrpy/helpers/simd_intrinsics.h), instruction-set definitions

### CONTR-0016

Evidence (descriptive behavior): [expression_utils.py](../nrpy/helpers/expression_utils.py), `_get_free_symbol_names_cached, get_params_commondata_symbols_from_expr_list`

### CONTR-0017

Evidence (descriptive behavior): [float_to_rational.py](../nrpy/helpers/float_to_rational.py), `f2r`

### CONTR-0018

Evidence (descriptive behavior): [functional.py](../nrpy/helpers/functional.py), `product`

### CONTR-0019

Evidence (descriptive behavior): [gpu_kernel.py](../nrpy/helpers/parallelization/gpu_kernel.py), `GPUKernel.c_function_call` and `GPUKernel.__init__`; [utilities.py](../nrpy/helpers/parallelization/utilities.py), `generate_kernel_and_launch_code`

### CONTR-0020

With positive `TINYDOUBLE`, strict coordinate comparisons can return fractional
values, non-strict comparisons shift their threshold, and all four helpers have
a shifted zero-denominator point. The SEOBNR ringdown constructor directly uses
the strict weights in its precession frequency.

Evidence (descriptive behavior): [Min_Max_and_Piecewise_Expressions.py](../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py), coordinate-helper return expressions; corroboration: [SEOBNRv5_merger_ringdown.py](../nrpy/equations/seobnr/SEOBNRv5_merger_ringdown.py), `SEOBNRv5_MergerRingdown_Quantities.__init__`, substitutes `Abs` and combines the strict helper values for `omega_prec`

### CONTR-0021

Wave constructors check only `uu` before paired gridfunction registration; the
Cartesian RHS and initial-data constructors check `glb_params_dict` before
CodeParameter registration, so existing names can suppress required entries.

Evidence (descriptive behavior): [WaveEquation_RHSs.py](../nrpy/equations/wave_equation/WaveEquation_RHSs.py), `WaveEquation_RHSs.__init__`; [WaveEquationCurvilinear_RHSs.py](../nrpy/equations/wave_equation/WaveEquationCurvilinear_RHSs.py), `WaveEquationCurvilinear_RHSs.__init__`; [WaveEquation_Solutions_InitialData.py](../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py), `WaveEquation_solution_Cartesian.__init__`; corroboration: [params.py](../nrpy/params.py), `CodeParameter._register_globally`, stores CodeParameters in the separate `glb_code_params_dict`

### CONTR-0022

Direct zero-velocity `LorentzBoost` construction yields undefined spatial
entries instead of the identity because the implementation divides by the
squared speed without a zero branch; direct zero substitution exposes `0/0`.

Evidence (descriptive behavior): [LorentzBoost.py](../nrpy/equations/general_relativity/LorentzBoost.py), `LorentzBoost._compute_LorentzMatrix` (the trusted tensor samples do not assert the identity limit)

### CONTR-0023

A reflected stencil node reads the record at its own Cartesian position, so the
per-node rotation by `phi_ref - node_phi_ref` reverses components with an odd
number of x and y indices; the final rotation by `target_phi - phi_ref` is
correct.

Evidence (descriptive behavior): [azimuthal_symmetry_spatial_lagrange_interpolation.py](../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/azimuthal_symmetry_spatial_lagrange_interpolation.py), `register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation`; corroboration: [output_raytracing_data.py](../nrpy/infrastructures/BHaH/diagnostics/output_raytracing_data.py), module description of stored Cartesian tensor components (no test covers reflected stencil nodes)
