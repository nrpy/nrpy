# Contradictions

> Register for contested and stale KB claims. Plain Markdown only. · Status: confirmed

Known contested/stale claims are tracked below. A row records
source-side truth and containment; it does not imply that NRPy source was fixed.

## Register

| ID | Claim | Claim status | Source A | Source B | Authority decision | Affected pages | Page-status rationale | Owner/trigger | Resolution test | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CONTR-0001 | `manga_bhah_lib` tells users to run `./bhah_lib` and says `bhah_lib.par` exists, but its assembly path requests a library target and does not call the parser/default-parfile writer. | stale | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) final two `print()` calls | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) commented parser registration and `output_CFunctions_function_prototypes_and_construct_Makefile(create_lib=True)` call; [`Makefile_helpers.py`](../nrpy/infrastructures/BHaH/Makefile_helpers.py) `output_CFunctions_function_prototypes_and_construct_Makefile` | For descriptive current behavior, precise registrations and writer calls decide. They emit `libbhah_lib.so` on Linux or `libbhah_lib.dylib` on Darwin, not an executable `bhah_lib`; no invoked writer emits `bhah_lib.par`. Both generic final prints are stale. | [BHaH Lifecycle And Project Assembly](infrastructures/bhah/lifecycle-and-project-assembly.md); [Matter TOV Workflows](examples/matter-tov-workflows.md) | Both pages remain `confirmed`: these are bounded final-message claims, corrected inline with explicit non-guarantees; their principal lifecycle/workflow answers remain supported. | BHaH/example owners; trigger when `manga_bhah_lib.py` final messaging, parser/parfile registration, or Makefile arguments change. | In an isolated output root, run the generator, capture stdout, inspect emitted files and the Makefile `all` target, and build. Pass only if stdout gives library-consumer guidance, never instructs executing `./bhah_lib`, mentions a parfile only if one is emitted, and `make` produces the platform library named by the Makefile. | Configured CI does not cover this resolution test. |
| CONTR-0002 | `sebobv1_jax` emits `Commondata(..., a_f=a_f)`, but its batch field registration supplies fourteen names/descriptions and thirteen dtypes/defaults; `zip()` therefore omits `a_f` from the generated dataclass. | resolved | [`SEOBNRv5_aligned_spin_coefficients.py`](../nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py) `register_PyFunction_SEOBNRv5_aligned_spin_coefficients` emitted return | [`sebobv1_jax.py`](../nrpy/examples/sebobv1_jax.py) `register_commondata_params` call and [`commondata.py`](../nrpy/infrastructures/JAX/commondata.py) `register_commondata_params` / `generate_commondata_dataclass` | Descriptive behavior is decided by all three code paths together. `register_commondata_params` now raises `ValueError` when its four lists differ in length, `sebobv1_jax` supplies fourteen entries in each list, and the generated `Commondata` therefore has every keyword, including `a_f`, that the emitted function passes. | [SEBOBv1 JAX Workflow](infrastructures/jax/sebobv1-jax-workflow.md); [Commondata And PyFunction Registry](infrastructures/jax/commondata-and-pyfunction-registry.md); [Waveform JAX PN Generators](examples/waveform-jax-pn-generators.md); [Example Generator Catalog](examples/example-generator-catalog.md) | The inline markers are removed. `SEBOBv1 JAX Workflow` returns to `confirmed` because its coefficient interface now matches the registered fields; the other pages already were `confirmed` and now state the corrected registration. | JAX/example owners; trigger when batch registration loses its length check, the four lists change, or the emitted return signature changes. | In an isolated output root, generate the project; assert batch metadata passes equal-length validation; import the generated coefficient function; call it with a valid finite input set; assert every returned keyword is a generated `Commondata` field and every returned numerical field is finite; if `a_f` remains in the intended return contract, assert that `a_f` is present and finite; then run the generated smoke test. All applicable checks must pass. | Configured CI installs the generated package and calls the coefficient function once, but asserts no returned field or value, so this resolution test is not run in CI. Generated code became runnable only together with the `py_codegen()` change that prints right-hand sides without `assign_to` and `NRPyJaxPrinter._print_Max`; see [Python Codegen](core/python-codegen.md) and [CSE And Printer Support](core/helpers/cse-and-printer-support.md). |
| CONTR-0003 | Pylint policy distinguishes grandfathered existing tracked handwritten files from newly added handwritten files, but both configs use `fail-under=10`, the local wrapper uses a blanket `9.91` floor, and the workflow uses a blanket `9.5` floor; none enforces the tiered rule. | stale | [`coding_style.md`](../coding_style.md) `## Static Analysis Configuration` | [`.pylintrc`](../.pylintrc) and [`.pylintrc_python36`](../.pylintrc_python36) `[MASTER]`; [`single_file_static_analysis.sh`](../.github/single_file_static_analysis.sh) `run_test_step`; [`main.yml`](../.github/workflows/main.yml) `static-analysis` | Existing tracked handwritten files retain their pre-change score without regression, including legacy scores at or below `9.5`; newly added handwritten files require `10.00/10.00`. Enforcement must classify the path from the base revision, compare legacy base and proposed ratings under the same tool/config, require exact `10.00` for new paths, and reject missing or unparseable ratings. | [Static Analysis](validation/static-analysis.md) | The page remains `confirmed`: the enforcement mismatch is one bounded, directly sourced Detail claim with an inline marker and explicit non-guarantee; the principal command and configuration answer remains reliable. | Validation and contributor-guidance owners; trigger when the tiered Pylint rule, either config, local wrapper classification/rating handling, or workflow classification/rating handling changes. | Deterministically inspect contributor guidance, both configs, local wrapper, and workflow. Pass only if documentation states the tiered rule and every enforcement path classifies paths from a base revision, rejects a lower legacy rating, requires exactly `10.00` for a new handwritten file, and rejects missing or unparseable ratings; generated trusted data remains exempt. | Configs and wrappers do not implement the required classification. |
| CONTR-0004 | Under its own full-Lie time-operator definition, Mewes et al.'s printed connection shift bracket adds the constraint-vector stretch a second time. | resolved | [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), `partial_0=partial_t-L_beta`, Eq. (29), and the evolution system after Eq. (32) | [`BSSN_RHSs.py`](../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__`; [`fCCZ4_constraints.py`](../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__`; [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (19) | fCCZ4 reuses `Lambdabar_rhsU_without_Brown_constraint_term`; its full Lie pair already contains `-C^k Dhat_k beta^i`. NRPy omits the duplicate printed term and retains `+2 C^i Dhat_k beta^k/3`. | [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md) | Page remains `confirmed`: the corrected equation, derivation, source discrepancy, and exact validation boundary are documented with its central scientific equations. | fCCZ4 equation modules; trigger when the evolved-connection definition, BSSN connection shift sector, fCCZ4 correction, or source crosswalk changes. | Inspect the full-Lie derivation and aggregate coefficient `-Lambdatilde^b delta^i_a+(2/3)Lambdatilde^i delta^a_b`; run the fCCZ4 RHS and gauge trusted-expression comparisons. For a BSSN-only Brown change, require every fCCZ4 value to remain unchanged. | The printed Mewes-minus-corrected residual is exactly the duplicated `-C^k Dhat_k beta^i`; Alic et al. corroborate the full corrected Cartesian coefficient, while Sanchis-Gual et al. corroborate one stretch and their spherical display contains the compatible divergence promotion. |
| CONTR-0005 | Mewes et al.'s prose around its conformal-factor alternatives reverses the signs implied by its displayed evolution equations and used by NRPy's established BSSN variables. | resolved | [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (33)-(34) and surrounding prose | [`BSSN_quantities.py`](../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__` | NRPy retains its established `W=exp(-2 phi)` and `chi=exp(-4 phi)` definitions, whose chain rules agree with the displayed evolution equations; sign-reversed prose is not adopted. | [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md) | Page remains `confirmed`: definitions and authority choice are explicit and owner validation covers all three conformal-factor options. | BSSN/fCCZ4 equation owners; trigger when conformal-factor definitions, source crosswalk, or option validation changes. | Run fCCZ4 constraint/RHS owner validation across `W`, `phi`, and `chi`; verify exact chain-rule definitions and trusted expressions. | This resolves a source-internal prose/equation mismatch without changing NRPy's established variables. |
| CONTR-0006 | Mewes et al.'s displayed advective 1+log equation omits the lapse factor present in Alic et al.'s exact advective form; Sanchis-Gual et al. corroborates the lapse factor only in a nonadvective form. | resolved | [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eq. (35) | [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (20) | NRPy adopts `partial_0 alpha=-2 alpha(K-2 Theta)` by reusing BSSN 1+log and adding `4 alpha Theta`; the missing lapse in Mewes Eq. (35) is treated as a typographical omission. | [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md) | Page remains `confirmed`: the adopted gauge formula and source decision are explicit and validated across Cartesian and curved reference metrics. | fCCZ4 gauge modules; trigger when lapse options, 1+log correction, source authority, or gauge validation changes. | Run the fCCZ4 gauge trusted-expression comparisons and verify the 1+log correction is exactly `4 alpha Theta`. | Sanchis-Gual et al. Eq. (2.25) corroborates the lapse factor but not the advective operator; frozen lapse is unaffected. |
| CONTR-0007 | The frozen commissioned YBS-MOM source proposes an auxiliary hyperbolic-relaxation field, but the corrected implementation contract requires the Yo--Lin--Cao term directly with CAHD-style timestep scaling and no new evolved state. | resolved | [preserved YBS-MOM specification](../raw/source-docs/ybs-momentum-damping-spec.md), `Proposed system` and auxiliary-field analysis | [`BSSN_RHSs.py`](../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch; commissioned user clarification | The user clarification controls the commissioned design: apply the direct covariant STF momentum gradient with `C_YBS_mom*CFL_FACTOR*DSMINGF`, exactly analogous to CAHD's timestep-scaled parabolic coefficient. Retain the frozen source as historical provenance; do not add `qD` or any other evolved cleaner field. A later maintainer decision multiplies the coefficient by the conformal factor `W`; the direct, stateless formulation is unchanged. | [YBS-MOM Timestep-Scaled Momentum Adjustment](equations/general-relativity/ybs-momentum-damping.md); [BSSN Family](equations/general-relativity/bssn-family.md); [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md); [BHaH GR Application Wiring](infrastructures/bhah/gr-application-wiring.md) | Pages remain `confirmed`: code and compiled pages now agree on the direct no-state contract; the superseded proposal is explicitly bounded. | BSSN/fCCZ4/BHaH owners; trigger when YBS-MOM variables, coefficient, momentum residual, STF projection, shared spacing, defaults, or source authority changes. | Verify no `qD`/`q_rhsD` gridfunction or output exists; compare the lower residual with canonical `BSSNconstraints.MU`; run all existing jointly YBS-enabled trusted dictionaries; generate ordinary and fisheye OpenMP/CUDA spacing paths; run static analysis and KB lint. | Resolution rejects auxiliary-state hyperbolization. The term remains parabolic in principal character; CAHD-style local timestep scaling satisfies the intended CFL compatibility. |
| CONTR-0008 | The Dendro BSSN leaf claimed the emitted `bssn_constraints*` kernels matched Dendro-GR's `bssn_constraints.cpp`. Upstream, that file's `enforce_bssn_constraints` is the det(gammabar)/tr(Abar) *enforcement*, and the Hamiltonian/momentum diagnostics live in `physcon.cpp` as `physical_constraints`; the claimed precedent named the wrong operation. | resolved | [`bssn-application-wiring.md`](infrastructures/dendro/bssn-application-wiring.md), "Dendro's own vocabulary" (superseded wording) | Dendro-GR `BSSN_GR/include/bssn_constraints.h`, `enforce_bssn_constraints`; `BSSN_GR/src/physcon.cpp`, `physical_constraints` | Host vocabulary governs host-side identifiers only. NRPy emits order-specific `BSSN_constraints` or `fCCZ4_constraints` diagnostic kernels and a separate `enforce_detgbar_equals_detghat_trAzero` algebraic-projection kernel. | [BSSN Application Wiring](infrastructures/dendro/bssn-application-wiring.md); [fCCZ4 Application Wiring](infrastructures/dendro/fccz4-application-wiring.md); [New Infrastructure Conformance](infrastructures/new-infrastructure-conformance.md) | Pages stay `provisional`/`confirmed`: the false sentence is rewritten and the conformance page now scopes host vocabulary explicitly. | Dendro owner; trigger when an emitted Dendro kernel is named for an upstream source file. | Generate and build both applications, then inspect their named constraint and algebraic-projection kernels. | The host and NRPy operation vocabularies are intentionally distinct. |
| CONTR-0009 | The Dendro BSSN leaf claimed the BSSN constraint builder "removes"/"deletes" the newly added AUX gridfunction names it does not write. The builder deletes nothing; it suppresses the two registrations at the source by setting `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` to `False` across the factory construction and restoring it in a `finally` block. | resolved | [`bssn-application-wiring.md`](infrastructures/dendro/bssn-application-wiring.md), "The DIAG-before-factory ordering, and why" and its claim-evidence block (superseded wording) | [`BSSN_constraints.py`](../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py), the `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` save/set/restore around the equation construction | For descriptive current behavior the code decides. No Dendro builder deletes a registered name from `glb_gridfcs_dict`. The flag is a construction-parameter memo key, so setting it forces a rebuild that never registers the two names. The `NUM_AUX_GFS = 0` outcome the leaf reports is unaffected. | [BSSN Application Wiring](infrastructures/dendro/bssn-application-wiring.md) | The page stays `provisional`: one bounded mechanism sentence and its claim-evidence block are rewritten, and the page's routed answer about what the BSSN builders do is unchanged. | Dendro owner; trigger when the AUX suppression mechanism, the `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` flag, or the emitted `NUM_AUX_GFS` changes. | Inspect `glb_gridfcs_dict` mutations in Dendro owners; generate the complete BSSN generated application and verify `NUM_AUX_GFS = 0`; build it and run authorized MPI diagnostics. | The superseded wording attested only the outcome, not the mechanism. |
| CONTR-0010 | The Dendro BSSN leaf called the private-builder-per-formulation layout "BHaH's own arrangement". BHaH's module holding both formulations' RHS systems holds no private builder and no private registrar. | resolved | [`bssn-application-wiring.md`](infrastructures/dendro/bssn-application-wiring.md), the Summary's layout paragraph (superseded wording) | [BHaH rhs_eval.py](../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval` | The code decides. BHaH's module holding both formulations' RHS systems uses public `register_CFunction_rhs_eval`, branches inline on `enable_fCCZ4`, and holds no private builder or registrar. Dendro's module layout does follow BHaH; the split inside the module does not. | [BSSN Application Wiring](infrastructures/dendro/bssn-application-wiring.md) | The page stayed `provisional` for this row: the layout paragraph is rewritten to separate the conforming module layout from the Dendro-only intra-module split, and the rest of the page is unaffected. CONTR-0011 subsequently opened and was resolved by consolidation. | Dendro owner; trigger when either Dendro leaf describes the intra-module layout, or when BHaH's `rhs_eval.py` gains a private per-formulation builder. | Require `rg '^def ' nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py` to identify only `register_CFunction_rhs_eval`, and require the Dendro BSSN leaf to attribute the per-formulation split to no peer. | The superseded Dendro wording conflated module layout with intra-module structure. |
| CONTR-0011 | The Dendro RHS and constraint modules duplicated private builders and registrars per formulation without a host requirement. | resolved | [New Infrastructure Conformance](infrastructures/new-infrastructure-conformance.md), `### Conformance is one-way` and `## Summary` | [rhs_eval.py](../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`; [BSSN_constraints.py](../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`; [fCCZ4_constraints.py](../nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py), `register_CFunction_fCCZ4_constraints` | Each generated operation has one public registrar. Formulation-specific constraint equations remain in their named modules. | [BSSN Application Wiring](infrastructures/dendro/bssn-application-wiring.md) | The affected page returns to provisional: this conflict is closed, while broader runtime qualification remains explicitly scoped. | Dendro owner; trigger when either module's per-formulation expression assembly or registrar structure changes, when a host requirement for a split is identified, or when the shared registrar changes. | Inspect for private per-formulation builder or registrar pairs. Generate and build both complete generated applications, then run authorized MPI diagnostics. | Resolved by one public registrar per operation. |
| CONTR-0012 | Python codegen and ordinary scalar C codegen expose substitution options, but a surviving postprocessed CSE temporary makes a nonempty substitution dictionary raise `TypeError`. | contested | [py_codegen.py](../nrpy/py_codegen.py), `PyCodeGen` and `py_codegen` option and returned-code descriptions; [c_codegen.py](../nrpy/c_codegen.py), `CCodeGen` and `c_codegen` substitution option | [py_codegen.py](../nrpy/py_codegen.py) and [c_codegen.py](../nrpy/c_codegen.py), assignment to `common_subexpression[1]`; [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `cse_postprocess` tuple entries | The implementation decides current behavior; tuple assignment prevents code return for these combinations. The KB documents the failure without treating it as successful substitution. | [Python Codegen](core/python-codegen.md); [C Codegen](core/c-codegen.md) | Both pages are `contested` because a public option combination fails. | C and Python codegen sources; recheck when substitution handling or CSE replacement representation changes. | For Python and ordinary non-FD, non-SIMD C generation, retain a common subexpression and supply a nonempty substitution dictionary with matching and nonmatching keys; require code return and agreement with the intended substituted expressions. Retain no-CSE and no-temporary cases. | This claim does not extend the scalar C failure to SIMD or FD paths. |
| CONTR-0013 | The public Python-codegen description maps input expressions to output names, but an output name could overwrite a CSE temporary or input symbol required by a later output. | resolved | [py_codegen.py](../nrpy/py_codegen.py), `py_codegen` input/output description | [py_codegen.py](../nrpy/py_codegen.py), numbered temporary construction and sequential output assignment | The implementation decides. `py_codegen()` now raises `ValueError` when an output name is a free symbol of a later printed output right-hand side, with CSE enabled or disabled, so it never returns assignments that read an overwritten temporary or input. Names that no later right-hand side reads remain accepted. | [Python Codegen](core/python-codegen.md) | The inline marker is removed. The page stays `contested` only because of CONTR-0012. | Python-codegen source; recheck when temporary naming, output assignment order, or `_check_no_read_after_overwrite` changes. | Generate `[x**3 + sin(x**3), x**3]` with output names `["tmp0", "b"]`, and `[x + 1, 2*x]` with `["x", "y"]` with CSE enabled and disabled; require `ValueError` in each case. Generate `[2*x, x + 1]` with `["y", "x"]`; require the unchanged assignments `y = 2*x` and `x = x + 1`. | Resolved by rejection, not by renaming temporaries; generated variable names are unchanged. The `py_codegen` doctests cover the resolution test. |
| CONTR-0014 | CSE sorting can preserve input-list order for tied structural keys, and unchecked generated constant or temporary names can capture input symbols. | contested | [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), documented helper behavior | [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `cse_preprocess, sort_cse_output_deterministically` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [CSE And Printer Support](core/helpers/cse-and-printer-support.md); [C Codegen](core/c-codegen.md) | Pages are `contested`: ordering and expression preservation require qualifications. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Permute replacements for `x**y` and `y**x` and require equal renamed output. Restore preprocessed constants in `2*x + _Integer_2`, and expand renamed output for `([(a0,x**3)], [a0+tmp0])`; require the original expressions or explicit collision rejection. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0015 | SIMD symbolic rewrites can omit the reciprocal numerator constant or emit Pow/Cbrt/Log/Sign names missing from the header. | contested | [simd.py](../nrpy/helpers/simd.py), `expr_convert_to_simd_intrins` | [simd_intrinsics.h](../nrpy/helpers/simd_intrinsics.h), instruction-set definitions | The implementation decides current behavior; missing numerator registration and definitions prevent a general supported-expression guarantee. | [SIMD And Intrinsic Support](core/helpers/simd-and-intrinsic-support.md) | Page is `contested`: missing constants and definitions affect generated calculations. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Verify declared constants for negative half-powers and compile supported Pow/Cbrt/Log/Sign expressions. | Including simd_intrinsics.h alone does not supply the missing SIMD names. |
| CONTR-0016 | The symbol-name helper claims SymPy free-symbol equivalence but removes free occurrences of bound-symbol names in limits or substitution points. | contested | [expression_utils.py](../nrpy/helpers/expression_utils.py), documented helper behavior | [expression_utils.py](../nrpy/helpers/expression_utils.py), `_get_free_symbol_names_cached, get_params_commondata_symbols_from_expr_list` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [Symbolic Expression Utilities](core/helpers/symbolic-expression-utilities.md) | Page is `contested`: symbol discovery, rational conversion, and product behavior each have source discrepancies. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Compare names with SymPy free_symbols for Integral, Lambda, and Subs, including `Integral(x,(x,0,x))`; verify parameter grouping retains registered free names. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0017 | f2r claims high-precision float conversion but appends zeros after scientific-notation exponents, changing values or creating invalid input; padding cannot recover lost float precision. | contested | [float_to_rational.py](../nrpy/helpers/float_to_rational.py), documented helper behavior | [float_to_rational.py](../nrpy/helpers/float_to_rational.py), `f2r` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [Symbolic Expression Utilities](core/helpers/symbolic-expression-utilities.md) | Page is `contested`: multiple helper guarantees require qualifications. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Require scientific and fixed-point strings for the same decimal value to produce equal rationals; cover `1e-7`, `1.2e-7`, and explicit padding values. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0018 | product describes Cartesian tuples but flattens iterable-valued elements, consumes one-shot inner iterators, and mishandles single scalar-valued inputs and repeat=0. | contested | [functional.py](../nrpy/helpers/functional.py), documented helper behavior | [functional.py](../nrpy/helpers/functional.py), `product` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [Symbolic Expression Utilities](core/helpers/symbolic-expression-utilities.md) | Page is `contested`: multiple helper guarantees require qualifications. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Compare with itertools.product for reusable and one-shot inputs, iterable-valued elements, single integer-valued inputs, and repeat=0. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0019 | Kernel call generation corrupts empty-argument calls; the CUDA wrapper ignores cfunc_type and may leave rewritten parameter accesses without a declared streamid. | contested | [utilities.py](../nrpy/helpers/parallelization/utilities.py), documented helper behavior | [gpu_kernel.py](../nrpy/helpers/parallelization/gpu_kernel.py), `GPUKernel.c_function_call and GPUKernel.__init__` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [Loop Kernel And Device Helpers](core/helpers/loop-kernel-and-device-helpers.md) | Page is `contested`: public wrapper inputs can produce invalid calls or a different return type. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Generate and compile an empty-argument host call, a CUDA helper with a requested non-void return type, and custom CUDA launch settings with parameter reads; require a valid call, the requested return type, and a declared streamid. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0020 | Coordinate-bound helper docstrings promise Boolean comparisons, but shifted quotients yield fractional values, shifted thresholds, and singular denominators; the SEOBNR ringdown frequency inherits those weights. | contested | [Min_Max_and_Piecewise_Expressions.py](../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py), coordinate-helper docstrings | [Min_Max_and_Piecewise_Expressions.py](../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py), coordinate-helper return expressions | Return expressions decide current behavior; exact comparison semantics are not guaranteed. | [GRHD](equations/grhd.md); [Precessing Rotations And Ringdown](equations/seobnr/precessing-rotations-and-ringdown.md) | GRHD remains `confirmed`: one bounded helper claim is qualified inline. The SEOBNR page is `contested` because the frequency, Euler angles, and QNM shifts inherit the discrepancy. | GRHD comparison helpers and SEOBNR ringdown constructor; recheck on changes to offsets, formulas, or documented comparison semantics. | Interpret nrpyAbs as absolute value and compare with true predicates on both sides of zero and at plus/minus TINYDOUBLE; require finite 0-or-1 results matching the predicates or an explicitly restricted interface. Verify the corresponding ringdown frequency weights and downstream angle/QNM expressions. | Trusted samples replacing nrpyAbs with sin, cos, or exp do not decide these semantics. |
| CONTR-0021 | Wave constructors can leave required gridfunctions or CodeParameters unregistered when existing names trigger incomplete or incorrect guards. | contested | [WaveEquation_RHSs.py](../nrpy/equations/wave_equation/WaveEquation_RHSs.py), constructor registration comments; [WaveEquation_Solutions_InitialData.py](../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py), constructor registration comments | [WaveEquation_RHSs.py](../nrpy/equations/wave_equation/WaveEquation_RHSs.py), [WaveEquationCurvilinear_RHSs.py](../nrpy/equations/wave_equation/WaveEquationCurvilinear_RHSs.py), and [WaveEquation_Solutions_InitialData.py](../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py), constructors; [params.py](../nrpy/params.py), CodeParameter registration | Actual registry checks decide current behavior; callers must provide compatible entries when reusing state. | [Wave Equation](equations/wave-equation.md) | Page is `contested`: independent gridfunction and parameter-registration conditions invalidate unconditional setup claims. | Wave equation constructors; recheck when registry guards or registration validation changes. | In separate isolated processes, pre-register uu without vv, or with incompatible metadata, and pre-register NRPy parameters named wavespeed or time; require compatible gridfunction and CodeParameter entries or explicit rejection. Keep a clean-registry case. | Curvilinear RHS wavespeed registration is unconditional; its uu-only gridfunction guard still applies. |
| CONTR-0022 | A zero Lorentz boost should be the identity, but direct zero-velocity construction produces undefined spatial entries. | contested | [LorentzBoost.py](../nrpy/equations/general_relativity/LorentzBoost.py), documented Lorentz-transformation purpose | [LorentzBoost.py](../nrpy/equations/general_relativity/LorentzBoost.py), `_compute_LorentzMatrix` spatial block | The implemented division by squared speed decides current behavior; zero input is not handled by its identity limit. | [Metric Conversions And Matter](equations/general-relativity/metric-conversions-and-matter.md) | Page remains `confirmed`: this is one bounded helper defect, qualified inline; other metric and matter conversions remain valid. | LorentzBoost constructor; recheck when speed-domain handling or the spatial matrix formula changes. | Construct the direct zero boost and its inverse; require finite identity matrices and unchanged transformed tensors, retaining nonzero subluminal cases. | This KB correction documents the source defect without changing Python. |

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

Claim evidence:
- Claim: Existing tracked handwritten files retain their pre-change Pylint score without regression, including legacy scores at or below `9.5`; newly added handwritten files require `10.00/10.00`.
- Role: normative rule
- Deciding authority: registered [`coding_style.md`](../coding_style.md), `## Static Analysis Configuration`, as updated by the commissioned policy decision
- Corroboration: [Static Analysis](validation/static-analysis.md) and [Code Test Policy](validation/code-test-policy.md) apply the rule; registered [`single_file_static_analysis.sh`](../.github/single_file_static_analysis.sh), `run_test_step`, and [`main.yml`](../.github/workflows/main.yml), `static-analysis`, identify enforcement surfaces but do not implement the distinction

The configs and wrappers use blanket floors. None classifies paths from a base
revision or checks legacy score regression.

### CONTR-0004

Mewes et al. define
`Lambdatilde^i=DeltaGamma^i+C^i` and
`partial_0=partial_t-L_beta`. The full contravariant-vector Lie pair obeys

```text
beta^k partial_k Lambdatilde^i - Lambdatilde^k partial_k beta^i
= beta^k Dhat_k Lambdatilde^i - Lambdatilde^k Dhat_k beta^i.
```

fCCZ4 reuses the BSSN base before Brown's BSSN-only adjustment, retaining the
full `-Lambdatilde^k Dhat_k beta^i` stretch. The companion
`+2 C^i Dhat_k beta^k/3` promotes the reused geometric divergence coefficient.

NRPy corrects the apparent paper error. Its aggregate Cartesian coefficient
for `partial_b beta^a` is exactly
`-Lambdatilde^b delta^i_a+(2/3)Lambdatilde^i delta^a_b`. That property was
settled by the derivation and review recorded here. The equation owner now
locks the complete reviewed expression set through trusted dictionaries; it
does not add a Mewes-specific reconstruction, coefficient, or shear test.

Claim evidence:
- Claim: fCCZ4 reuses the pre-Brown BSSN base, retains `+2 C^i Dhat_k beta^k/3`, and produces exactly one evolved-vector stretch.
- Role: public/scientific contract
- Deciding authority: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), its `partial_0` definition, Eq. (29), and the evolution system following Eq. (32), composed with [`BSSN_RHSs.py`](../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__`, [`fCCZ4_constraints.py`](../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__`, and [`fCCZ4_RHSs.py`](../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__` and module `__main__`
- Corroboration: [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (19), corroborates the full corrected Cartesian coefficient; [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), Eqs. (2.11) and (2.17), corroborates one stretch, while spherical Eq. (3.14) contains the compatible divergence promotion but is not identical to the general display off constraint

### CONTR-0005

NRPy uses `W=exp(-2 phi)` and `chi=exp(-4 phi)`. These definitions agree with
Mewes et al.'s displayed conformal-factor equations; the neighboring
sign-reversed prose is not part of NRPy's scientific contract.

Claim evidence:
- Claim: NRPy resolves the Mewes conformal-factor prose/equation mismatch by retaining `W=exp(-2 phi)` and `chi=exp(-4 phi)`, matching the displayed equations and established BSSN variables.
- Role: public/scientific contract
- Deciding authority: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (33)-(34); [`BSSN_quantities.py`](../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__`
- Corroboration: [`fCCZ4_constraints.py`](../nrpy/equations/general_relativity/fCCZ4_constraints.py), module `__main__`, exercises `W`, `phi`, and `chi` across Cartesian and SinhSpherical owner cases

### CONTR-0006

NRPy adopts advective 1+log as
`partial_0 alpha=-2 alpha(K-2 Theta)`. Its fCCZ4 gauge owner obtains this by
reusing the BSSN branch and adding `4 alpha Theta`; it does not adopt the
lapse-free right-hand side printed in Mewes et al. Eq. (35).

Claim evidence:
- Claim: NRPy treats the missing lapse in Mewes et al. Eq. (35) as a typographical omission and adopts `partial_0 alpha=-2 alpha(K-2 Theta)` for its supported fCCZ4 1+log branch.
- Role: public/scientific contract
- Deciding authority: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eq. (35); [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (20); [`fCCZ4_gauge_RHSs.py`](../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py), `fCCZ4_gauge_RHSs`
- Corroboration: [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), Eq. (2.25), corroborates the lapse factor in a nonadvective form

### CONTR-0007

The frozen commissioned source interpreted “hyperbolize” as introducing an
auxiliary relaxation covector. The later explicit user clarification instead
defines the intended analogy to CAHD: keep the direct parabolic
momentum-gradient adjustment and multiply its coefficient by the local
timestep scale. Current code therefore adds no evolved field and leaves
initial data, restart, boundary, AMR, and dissipation state routes unchanged.

Claim evidence:
- Claim: Current YBS-MOM uses the direct covariant STF gradient of the matter-complete lower momentum residual, multiplied by `C_YBS_mom*CFL_FACTOR*DSMINGF*W`, where the weight `W` is a later maintainer decision; the frozen auxiliary-relaxation proposal is superseded and no cleaner state is permitted.
- Role: public/scientific contract
- Deciding authority: commissioned user clarification composed with [`BSSN_RHSs.py`](../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch
- Corroboration: [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (56), and [Etienne, Phys. Rev. D 110, 064045](https://doi.org/10.1103/PhysRevD.110.064045), Eq. (26)

## Rules

- Follow the contradiction contract in [Schema](SCHEMA.md#contradiction-contract).
- Open a row when registered sources disagree, when a living source supersedes
  a KB claim, or when a page is marked `contested` or `stale`.
- Update affected-page links and exact inline markers together.
- Resolve a row only after its resolution test passes and all affected pages,
  reverse dependents, aliases, typed neighbors, and targeted wiki hits are
  reconciled.

### CONTR-0008

The Dendro BSSN leaf recorded that the emitted `bssn_constraints*` CFunctions
matched `bssn_constraints.cpp` in Dendro-GR's own BSSN solver. Reading that
source shows the match was in name only: upstream's `enforce_bssn_constraints`
rescales the conformal metric and removes the trace of the conformal extrinsic
curvature — the operation BHaH and ETLegacy call
`enforce_detgbar_equals_detghat_trAzero` — while the Hamiltonian and momentum
constraint diagnostics are computed in `physcon.cpp`. The kernels are now named
for their operation with the NRPy vocabulary, and the host vocabulary is scoped
to host-side identifiers.

Claim evidence:
- Claim: Dendro-GR's `bssn_constraints.cpp` implements the algebraic det/trace enforcement, not the Hamiltonian/momentum constraint diagnostics; NRPy emits order-specific `BSSN_constraints` or `fCCZ4_constraints` kernels and a separate algebraic-projection kernel.
- Role: descriptive behavior
- Deciding authority: Dendro-GR `BSSN_GR/include/bssn_constraints.h`, `enforce_bssn_constraints`; `BSSN_GR/src/physcon.cpp`, `physical_constraints`; [BHaH enforce_detgbar_equals_detghat_trAzero.py](../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `name`
- Corroboration: [BSSN_constraints.py](../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`; [fCCZ4_constraints.py](../nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py), `register_CFunction_fCCZ4_constraints`; [enforce_detgbar_equals_detghat_trAzero.py](../nrpy/infrastructures/Dendro/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`

### CONTR-0009

The Dendro BSSN leaf described the constraint builder as removing, and in its
claim-evidence block as deleting, the newly added AUX gridfunction names the
kernel does not write. The builder does no such thing. It reads
`register_M_and_LAMBDA_CONSTRAINT_gridfunctions`, sets it to `False`, constructs
`BSSN_constraints` inside a `try`, and restores the previous value in a
`finally` block, so `M` and `LAMBDA_CONSTRAINT` are never registered in the
first place. The code comment at that site says exactly this. Suppressing the
registration is also what makes a deletion pass unnecessary: the same
construction pulls in the evolved state, so a pass that deleted newly added
names would have to tell the two constraint names apart from it.

Claim evidence:
- Claim: the Dendro BSSN constraint builder suppresses the `M` and `LAMBDA_CONSTRAINT` registrations by setting `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` to `False` across the `BSSN_constraints` construction and restoring the previous value afterwards; no Dendro production builder deletes registered fields.
- Role: descriptive behavior
- Deciding authority: [BSSN_constraints.py](../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py), the `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` save/set/restore around the equation construction
- Corroboration: [BSSN_constraints.py](../nrpy/equations/general_relativity/BSSN_constraints.py), the `group="AUX"` registrations gated on that CodeParameter

### CONTR-0010

The Dendro BSSN leaf attributed the private-builder-per-formulation layout to
BHaH. BHaH's `general_relativity/rhs_eval.py` defines the public
`register_CFunction_rhs_eval` and branches inline on `enable_fCCZ4`.
What Dendro does take from BHaH is the module layout -- one
module per artifact, both formulations behind one boolean -- and that half of
the sentence stands; the private builder and private registrar per formulation
inside the module were Dendro's own. CONTR-0011 subsequently removed that
duplication; this entry preserves the attribution correction.

Claim evidence:
- Claim: `BHaH/general_relativity/rhs_eval.py` holds one public registration function and no private per-formulation builder or registrar; the former private-builder-per-formulation split inside a Dendro module was Dendro's own, not BHaH's arrangement.
- Role: descriptive behavior
- Deciding authority: [BHaH rhs_eval.py](../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval` and its `enable_fCCZ4` branches
- Corroboration: none available; BHaH's own module is the whole of the evidence, and no Dendro leaf records the per-formulation split as a peer practice

### CONTR-0011

Resolved. Both modules now contain one public registrar. In `rhs_eval.py`,
`enable_fCCZ4` selects the expression set, the CAHD parameter registration, the
extra `Theta_fCCZ4` right-hand-side entry, the evolved-field list and its
validation, and the formulation label in the description; lowering, wrappers,
and registration share one path. The previous duplication had no host
requirement and was removed. Other top-level Dendro modules use
`enable_fCCZ4` to select state, checkpoint, kernel, and label text, as
[New Infrastructure Conformance](infrastructures/new-infrastructure-conformance.md)
lists.

Claim evidence:
- Claim: the unsupported per-formulation builder and registrar duplication is removed in both Dendro modules; each operation now has one public registrar, resolving CONTR-0011.
- Role: normative rule
- Deciding authority: [new-infrastructure-conformance.md](infrastructures/new-infrastructure-conformance.md), `Conformance is one-way`; [rhs_eval.py](../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`; [BSSN_constraints.py](../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py), `register_CFunction_BSSN_constraints`; [fCCZ4_constraints.py](../nrpy/infrastructures/Dendro/general_relativity/fCCZ4_constraints.py), `register_CFunction_fCCZ4_constraints`
- Corroboration: [BHaH rhs_eval.py](../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), inline formulation branching; complete generated BSSN and fCCZ4 project tests exercise both consolidated paths


### CONTR-0012

The postprocessed CSE replacements are tuples, but a nonempty
`postproc_substitution_dict` makes `py_codegen()` assign a new right-hand side
into `common_subexpression[1]`. The ordinary non-FD, non-SIMD CSE loop in
`c_codegen()` performs the same assignment. A surviving temporary therefore causes
`TypeError` even when no substitution key matches. The no-CSE path and the
final-expression path without surviving temporaries do not perform that tuple
assignment. This is an implementation limitation, not evidence of successful
support for every combination of options.

Claim evidence:
- Claim: In Python codegen and ordinary non-FD, non-SIMD C codegen, a nonempty substitution dictionary and a surviving postprocessed CSE temporary cause tuple-item assignment to raise `TypeError`; the no-CSE and no-surviving-temporary cases avoid that assignment.
- Role: descriptive behavior
- Deciding authority: [py_codegen.py](../nrpy/py_codegen.py), `py_codegen`, and [c_codegen.py](../nrpy/c_codegen.py), `c_codegen`, scalar common-subexpression substitution loops
- Corroboration: [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `cse_postprocess`, constructs the returned replacement entries as tuples.

### CONTR-0013

Resolved. `py_codegen()` emits output assignments sequentially into one Python
namespace. Previously, for `[x**3 + sin(x**3), x**3]` and targets
`["tmp0", "b"]`, the first output assignment replaced the temporary that the
second still read, and an output could likewise overwrite an input symbol read
later. `_check_no_read_after_overwrite` now raises `ValueError` in both cases,
on the CSE and no-CSE paths, before any code is returned. Requests in which no
later right-hand side reads an output name emit the same assignments as
before; the generated SEBOBv1 coefficient module is unchanged. The fix rejects
colliding names rather than renaming temporaries, so generated variable names
do not change.

Claim evidence:
- Claim: Python codegen raises `ValueError` when an earlier output name is a free symbol of a later printed output right-hand side, covering both CSE temporaries and input symbols, and otherwise emits unchanged assignments.
- Role: descriptive behavior
- Deciding authority: [py_codegen.py](../nrpy/py_codegen.py), `_check_no_read_after_overwrite` and its calls in `py_codegen`
- Corroboration: none available; the `py_codegen` doctests in the deciding module exercise the rejected and accepted cases but are not a separate source.

### CONTR-0014

Claim evidence:
- Claim: CSE sorting can preserve input-list order for tied structural keys, and unchecked generated constant or temporary names can capture input symbols.
- Role: descriptive behavior
- Deciding authority: [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `cse_preprocess, sort_cse_output_deterministically`
- Corroboration: `none available`; the stated counterexamples follow from the implementation, while its docstring or paired symbolic definition supplies the intended behavior.

### CONTR-0015

Claim evidence:
- Claim: SIMD symbolic rewrites can omit the reciprocal numerator constant or emit Pow/Cbrt/Log/Sign names missing from the header.
- Role: descriptive behavior
- Deciding authority: [simd.py](../nrpy/helpers/simd.py), `expr_convert_to_simd_intrins`; [simd_intrinsics.h](../nrpy/helpers/simd_intrinsics.h), instruction-set definitions
- Corroboration: `none available`; the stated counterexamples follow from the implementation, while its docstring or paired symbolic definition supplies the intended behavior.

### CONTR-0016

Claim evidence:
- Claim: The symbol-name helper claims SymPy free-symbol equivalence but removes free occurrences of bound-symbol names in limits or substitution points.
- Role: descriptive behavior
- Deciding authority: [expression_utils.py](../nrpy/helpers/expression_utils.py), `_get_free_symbol_names_cached, get_params_commondata_symbols_from_expr_list`
- Corroboration: `none available`; the stated counterexamples follow from the implementation, while its docstring or paired symbolic definition supplies the intended behavior.

### CONTR-0017

Claim evidence:
- Claim: f2r claims high-precision float conversion but appends zeros after scientific-notation exponents, changing values or creating invalid input; padding cannot recover lost float precision.
- Role: descriptive behavior
- Deciding authority: [float_to_rational.py](../nrpy/helpers/float_to_rational.py), `f2r`
- Corroboration: `none available`; the stated counterexamples follow from the implementation, while its docstring or paired symbolic definition supplies the intended behavior.

### CONTR-0018

Claim evidence:
- Claim: product describes Cartesian tuples but flattens iterable-valued elements, consumes one-shot inner iterators, and mishandles single scalar-valued inputs and repeat=0.
- Role: descriptive behavior
- Deciding authority: [functional.py](../nrpy/helpers/functional.py), `product`
- Corroboration: `none available`; the stated counterexamples follow from the implementation, while its docstring or paired symbolic definition supplies the intended behavior.

### CONTR-0019

Claim evidence:
- Claim: Kernel call generation corrupts empty-argument calls; the CUDA wrapper ignores cfunc_type and may leave rewritten parameter accesses without a declared streamid.
- Role: descriptive behavior
- Deciding authority: [gpu_kernel.py](../nrpy/helpers/parallelization/gpu_kernel.py), `GPUKernel.c_function_call` and `GPUKernel.__init__`; [utilities.py](../nrpy/helpers/parallelization/utilities.py), `generate_kernel_and_launch_code`
- Corroboration: `none available`; the stated counterexamples follow from the implementation, while its docstring or paired symbolic definition supplies the intended behavior.

### CONTR-0020

Claim evidence:
- Claim: With positive `TINYDOUBLE`, strict coordinate comparisons can return fractional values, non-strict comparisons shift their threshold, and all four helpers have a shifted zero-denominator point. The SEOBNR ringdown constructor directly uses the strict weights in its precession frequency.
- Role: descriptive behavior
- Deciding authority: [Min_Max_and_Piecewise_Expressions.py](../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py), coordinate-helper return expressions
- Corroboration: [SEOBNRv5_merger_ringdown.py](../nrpy/equations/seobnr/SEOBNRv5_merger_ringdown.py), `SEOBNRv5_MergerRingdown_Quantities.__init__`, substitutes `Abs` and combines the strict helper values for `omega_prec`.

### CONTR-0021

Claim evidence:
- Claim: Wave constructors check only `uu` before paired gridfunction registration; the Cartesian RHS and initial-data constructors check `glb_params_dict` before CodeParameter registration, so existing names can suppress required entries.
- Role: descriptive behavior
- Deciding authority: [WaveEquation_RHSs.py](../nrpy/equations/wave_equation/WaveEquation_RHSs.py), `WaveEquation_RHSs.__init__`; [WaveEquationCurvilinear_RHSs.py](../nrpy/equations/wave_equation/WaveEquationCurvilinear_RHSs.py), `WaveEquationCurvilinear_RHSs.__init__`; [WaveEquation_Solutions_InitialData.py](../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py), `WaveEquation_solution_Cartesian.__init__`
- Corroboration: [params.py](../nrpy/params.py), `CodeParameter._register_globally`, stores CodeParameters in the separate `glb_code_params_dict`.

### CONTR-0022

Claim evidence:
- Claim: Direct zero-velocity LorentzBoost construction yields undefined spatial entries instead of the identity because the implementation divides by the squared speed without a zero branch.
- Role: descriptive behavior
- Deciding authority: [LorentzBoost.py](../nrpy/equations/general_relativity/LorentzBoost.py), `LorentzBoost._compute_LorentzMatrix`
- Corroboration: none available for zero-velocity handling; direct zero substitution exposes `0/0`, and the current trusted tensor samples do not assert the identity limit.
