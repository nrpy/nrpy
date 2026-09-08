# Contradictions

> Register for contested and stale KB claims. Plain Markdown only. · Status: confirmed · Last reconciled: 09-07-2026

Known contested/stale claims as of 09-07-2026 are tracked below. A row records
source-side truth and containment; it does not imply that NRPy source was fixed.

## Register

| ID | Claim | Claim status | Source A | Source B | Authority decision | Affected pages | Page-status rationale | Owner/trigger | Resolution test | Opened | Resolved | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CONTR-0001 | `manga_bhah_lib` tells users to run `./bhah_lib` and says `bhah_lib.par` exists, but its assembly path requests a library target and does not call the parser/default-parfile writer. | stale | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) final two `print()` calls | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) commented parser registration and `output_CFunctions_function_prototypes_and_construct_Makefile(create_lib=True)` call; [`Makefile_helpers.py`](../nrpy/infrastructures/BHaH/Makefile_helpers.py) `output_CFunctions_function_prototypes_and_construct_Makefile` | For descriptive current behavior, precise registrations and writer calls decide. They emit `libbhah_lib.so` on Linux or `libbhah_lib.dylib` on Darwin, not an executable `bhah_lib`; no invoked writer emits `bhah_lib.par`. Both generic final prints are stale. | [BHaH Lifecycle And Project Assembly](infrastructures/bhah/lifecycle-and-project-assembly.md); [Matter TOV Workflows](examples/matter-tov-workflows.md) | Both pages remain `confirmed`: these are bounded final-message claims, corrected inline with explicit non-guarantees; their principal lifecycle/workflow answers remain supported. | BHaH/example owners; trigger when `manga_bhah_lib.py` final messaging, parser/parfile registration, or Makefile arguments change. | In an isolated output root, run the generator, capture stdout, inspect emitted files and the Makefile `all` target, and build. Pass only if stdout gives library-consumer guidance, never instructs executing `./bhah_lib`, mentions a parfile only if one is emitted, and `make` produces the platform library named by the Makefile. | 06-30-2026 | - | Current CI commands for this generator are commented out, so generation/build results are `not-run` here. |
| CONTR-0002 | `sebobv1_jax` emits `Commondata(..., a_f=a_f)`, but its batch field registration supplies fourteen names/descriptions and thirteen dtypes/defaults; `zip()` therefore omits `a_f` from the generated dataclass. | contested | [`SEOBNRv5_aligned_spin_coefficients.py`](../nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py) `register_PyFunction_SEOBNRv5_aligned_spin_coefficients` emitted return | [`sebobv1_jax.py`](../nrpy/examples/sebobv1_jax.py) `register_commondata_params` call and [`commondata.py`](../nrpy/infrastructures/JAX/commondata.py) `register_commondata_params` / `generate_commondata_dataclass` | Descriptive behavior is decided by all three code paths together: Python `zip()` truncation registers thirteen fields, dataclass emission follows that registry, and the emitted function passes unsupported keyword `a_f`. Project generation alone does not establish callable coefficient initialization. | [SEBOBv1 JAX Workflow](infrastructures/jax/sebobv1-jax-workflow.md); [Commondata And PyFunction Registry](infrastructures/jax/commondata-and-pyfunction-registry.md); [Waveform JAX PN Generators](examples/waveform-jax-pn-generators.md); [Example Generator Catalog](examples/example-generator-catalog.md) | `SEBOBv1 JAX Workflow` is `contested` because the mismatch invalidates its sole generated coefficient interface and changes its Summary-level routed answer. Other affected pages remain `confirmed` only where they describe this as a bounded noncentral detail and make no callable-function guarantee. | JAX/example owners; trigger when batch registration gains length validation, the four lists become equal, or the emitted return signature changes. | In an isolated output root, generate the project; assert batch metadata passes equal-length validation; import the generated coefficient function; call it with a valid finite input set; assert every returned keyword is a generated `Commondata` field and every returned numerical field is finite; if `a_f` remains in the intended return contract, assert that `a_f` is present and finite; then run the generated smoke test. All applicable checks must pass. | 07-06-2026 | - | Workflow config invokes generation only; generated-package install, import, coefficient call, and result checks are `not-run` here. |
| CONTR-0003 | Pylint policy distinguishes grandfathered existing tracked handwritten files from newly added handwritten files, but both configs use `fail-under=10`, the local wrapper uses a blanket `9.91` floor, and the workflow uses a blanket `9.5` floor; none enforces the tiered rule. | stale | [`coding_style.md`](../coding_style.md) `## Static Analysis Configuration` | [`.pylintrc`](../.pylintrc) and [`.pylintrc_python36`](../.pylintrc_python36) `[MASTER]`; [`single_file_static_analysis.sh`](../.github/single_file_static_analysis.sh) `run_test_step`; [`main.yml`](../.github/workflows/main.yml) `static-analysis` | Existing tracked handwritten files retain their pre-change score without regression, including legacy scores at or below `9.5`; newly added handwritten files require `10.00/10.00`. Enforcement must classify the path from the base revision, compare legacy base and proposed ratings under the same tool/config, require exact `10.00` for new paths, and reject missing or unparseable ratings. | [Static Analysis](validation/static-analysis.md) | The page remains `confirmed`: the enforcement mismatch is one bounded, directly sourced Detail claim with an inline marker and explicit non-guarantee; the principal command and configuration answer remains reliable. | Validation and contributor-guidance owners; trigger when the tiered Pylint rule, either config, local wrapper classification/rating handling, or workflow classification/rating handling changes. | Deterministically inspect contributor guidance, both configs, local wrapper, and workflow. Pass only if documentation states the tiered rule and every enforcement path classifies paths from a base revision, rejects a lower legacy rating, requires exactly `10.00` for a new handwritten file, and rejects missing or unparseable ratings; generated trusted data remains exempt. | 07-13-2026 | - | Pylint 4.0.6 gave both base and proposed versions of the two touched legacy files `10.00/10.00`; wrapper and workflow classification/enforcement remain `not-run`. |
| CONTR-0004 | Under its own full-Lie time-operator definition, Mewes et al.'s printed connection shift bracket adds the constraint-vector stretch a second time. | resolved | [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), `partial_0=partial_t-L_beta`, Eq. (29), and the evolution system after Eq. (32) | [`BSSN_RHSs.py`](../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__`; [`fCCZ4_constraints.py`](../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__`; [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (19) | Compose the variable definition and full Lie pair consistently: the reused BSSN base pair may be written with `Dhat` and already contains `-C^k Dhat_k beta^i`, so NRPy omits the duplicate printed term and retains `+2 C^i Dhat_k beta^k/3` to promote the geometric divergence coefficient. This corrects an apparent paper error, not a discretionary formulation departure. | [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md) | Page remains `confirmed`: the corrected equation, derivation, source discrepancy, and exact validation boundary are documented in its central scientific contract. | fCCZ4 equation owners; trigger when the evolved-connection definition, BSSN connection shift sector, fCCZ4 correction, or source crosswalk changes. | Inspect the full-Lie derivation and aggregate coefficient `-Lambdatilde^b delta^i_a+(2/3)Lambdatilde^i delta^a_b`; process and compare the complete five RHS and twelve gauge expression dictionaries, with affected moving-driver values updated and frozen-shift dictionaries unchanged. | 08-26-2026 | 08-27-2026 | The printed Mewes-minus-corrected residual is exactly the duplicated `-C^k Dhat_k beta^i`; Alic et al. corroborate the full corrected Cartesian coefficient, while Sanchis-Gual et al. corroborate one stretch and their spherical display contains the compatible divergence promotion. |
| CONTR-0005 | Mewes et al.'s prose around its conformal-factor alternatives reverses the signs implied by its displayed evolution equations and used by NRPy's established BSSN variables. | resolved | [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (33)-(34) and surrounding prose | [`BSSN_quantities.py`](../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__` | NRPy retains its established `W=exp(-2 phi)` and `chi=exp(-4 phi)` definitions, whose chain rules agree with the displayed evolution equations; sign-reversed prose is not adopted. | [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md) | Page remains `confirmed`: definitions and authority choice are explicit and owner validation covers all three conformal-factor options. | BSSN/fCCZ4 equation owners; trigger when conformal-factor definitions, source crosswalk, or option validation changes. | Run fCCZ4 constraint/RHS owner validation across `W`, `phi`, and `chi`; verify exact chain-rule definitions and trusted expressions. | 08-26-2026 | 08-26-2026 | This resolves a source-internal prose/equation mismatch without changing NRPy's established variables. |
| CONTR-0006 | Mewes et al.'s displayed advective 1+log equation omits the lapse factor present in Alic et al.'s exact advective form; Sanchis-Gual et al. corroborates the lapse factor only in a nonadvective form. | resolved | [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eq. (35) | [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (20) | NRPy adopts `partial_0 alpha=-2 alpha(K-2 Theta)` by reusing BSSN 1+log and adding `4 alpha Theta`; the missing lapse in Mewes Eq. (35) is treated as a typographical omission. | [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md) | Page remains `confirmed`: the adopted gauge formula and source decision are explicit and validated across Cartesian and curved reference metrics. | fCCZ4 gauge owners; trigger when lapse options, 1+log correction, source authority, or gauge validation changes. | Run all twelve fCCZ4 gauge trusted variants and verify the 1+log correction is exactly `4 alpha Theta`. | 08-26-2026 | 08-26-2026 | Sanchis-Gual et al. Eq. (2.25) corroborates the lapse factor but not the advective operator; frozen lapse is unaffected. |
| CONTR-0007 | The frozen commissioned YBS-MOM source proposes an auxiliary hyperbolic-relaxation field, but the corrected implementation contract requires the Yo--Lin--Cao term directly with CAHD-style timestep scaling and no new evolved state. | resolved | [preserved YBS-MOM specification](../raw/source-docs/ybs-momentum-damping-spec.md), `Proposed system` and auxiliary-field analysis | [`BSSN_RHSs.py`](../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch; commissioned user clarification | The user clarification controls the commissioned design: apply the direct covariant STF momentum gradient with `C_YBS_mom*CFL_FACTOR*DSMINGF`, exactly analogous to CAHD's timestep-scaled parabolic coefficient. Retain the frozen source as historical provenance; do not add `qD` or any other evolved cleaner field. | [YBS-MOM Timestep-Scaled Momentum Adjustment](equations/general-relativity/ybs-momentum-damping.md); [BSSN Family](equations/general-relativity/bssn-family.md); [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md); [BHaH GR Application Wiring](infrastructures/bhah/gr-application-wiring.md) | Pages remain `confirmed`: code and compiled pages now agree on the direct no-state contract; the superseded proposal is explicitly bounded. | BSSN/fCCZ4/BHaH owners; trigger when YBS-MOM variables, coefficient, momentum residual, STF projection, shared spacing, defaults, or source authority changes. | Verify no `qD`/`q_rhsD` gridfunction or output exists; compare the lower residual with canonical `BSSNconstraints.MU`; run all existing jointly YBS-enabled trusted dictionaries; generate ordinary and fisheye OpenMP/CUDA spacing paths; run static analysis and KB lint. | 08-30-2026 | 08-30-2026 | Resolution rejects auxiliary-state hyperbolization. The term remains parabolic in principal character; CAHD-style local timestep scaling satisfies the intended CFL compatibility. |
| CONTR-0008 | The Dendro BSSN leaf claimed the emitted `bssn_constraints*` kernels matched Dendro-GR's `bssn_constraints.cpp`. Upstream, that file's `enforce_bssn_constraints` is the det(gammabar)/tr(Abar) *enforcement*, and the Hamiltonian/momentum diagnostics live in `physcon.cpp` as `physical_constraints`; the claimed precedent named the wrong operation. | resolved | [`bssn-application-wiring.md`](infrastructures/dendro/bssn-application-wiring.md), "Dendro's own vocabulary" (superseded wording) | Dendro-GR `BSSN_GR/include/bssn_constraints.h`, `enforce_bssn_constraints`; `BSSN_GR/src/physcon.cpp`, `physical_constraints` | Host vocabulary governs host-side identifiers only; NRPy-emitted kernels take the BHaH/ETLegacy operation names with the stem prefix. The constraint-diagnostics kernel is `<stem>_constraints_eval`, the enforcement kernel `<stem>_enforce_detgbar_equals_detghat_trAzero`. | [BSSN Application Wiring](infrastructures/dendro/bssn-application-wiring.md); [fCCZ4 Application Wiring](infrastructures/dendro/fccz4-application-wiring.md); [New Infrastructure Conformance](infrastructures/new-infrastructure-conformance.md) | Pages stay `provisional`/`confirmed`: the false sentence is rewritten and the conformance page now scopes host vocabulary explicitly. | Dendro owner; trigger when an emitted Dendro kernel is named for an upstream source file. | `grep -rn bssn_constraints nrpy/infrastructures/Dendro` returns no CFunction name; the two generated projects build and pass ctest under the renamed kernels. | 09-06-2026 | 09-06-2026 | Found by the DI01 naming audit; the untracked Dendro-GR checkout was consulted read-only. |
| CONTR-0009 | The Dendro BSSN leaf claimed the BSSN constraint builder "removes"/"deletes" the newly added AUX gridfunction names it does not write. The builder deletes nothing; it suppresses the two registrations at the source by setting `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` to `False` across the factory construction and restoring it in a `finally` block. | resolved | [`bssn-application-wiring.md`](infrastructures/dendro/bssn-application-wiring.md), "The DIAG-before-factory ordering, and why" and its claim-evidence block (superseded wording) | [`constraints_eval.py`](../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py), the gate save/set/restore around the `BSSN_constraints` construction | For descriptive current behavior the code decides. No Dendro builder deletes a registered name from `glb_gridfcs_dict`; the subpackage's only production-code clear of that registry is the whole-registry sweep reset in `trusted_capture.reset_generation_state`, which no production path reaches; that registry's other `clear()` mentions in the subpackage are doctest fixtures. The gate is a construction-parameter memo key, so setting it forces a rebuild that never registers the two names. The `NUM_AUX_GFS = 0` outcome the leaf reports is unaffected. | [BSSN Application Wiring](infrastructures/dendro/bssn-application-wiring.md) | The page stays `provisional`: one bounded mechanism sentence and its claim-evidence block are rewritten, and the page's routed answer about what the BSSN builders do is unchanged. | Dendro owner; trigger when the AUX suppression mechanism, the `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` gate, or the emitted `NUM_AUX_GFS` changes. | `grep -rn "glb_gridfcs_dict" nrpy/infrastructures/Dendro` shows no `del` and no `pop`, and its only non-doctest `clear()` is the whole-registry sweep reset in `trusted_capture.py`; the BSSN generator emits `NUM_AUX_GFS = 0`; the two generated projects build and pass ctest. | 09-07-2026 | 09-07-2026 | Found by the inconsistencies survey as K1. The superseded wording carried a full `inspected=pass` validation tuple, so the tuple attested the outcome and not the mechanism. |
| CONTR-0010 | The Dendro BSSN leaf called the private-builder-per-formulation layout "BHaH's own arrangement". BHaH's only module holding both formulations' RHS systems holds no private builder and no private registrar. | resolved | [`bssn-application-wiring.md`](infrastructures/dendro/bssn-application-wiring.md), the Summary's layout paragraph (superseded wording) | [BHaH rhs_eval.py](../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval` | The code decides. BHaH's only module holding both formulations' RHS systems is one public `register_CFunction_rhs_eval` that branches inline on `enable_fCCZ4`; it holds no private builder and no private registrar. Dendro's module layout does follow BHaH; the split inside the module does not. | [BSSN Application Wiring](infrastructures/dendro/bssn-application-wiring.md) | The page stayed `provisional` for this row: the layout paragraph is rewritten to separate the conforming module layout from the Dendro-only intra-module split, and the rest of the page is unaffected. CONTR-0011 subsequently opened and was resolved by consolidation. | Dendro owner; trigger when either Dendro leaf describes the intra-module layout, or when BHaH's `rhs_eval.py` gains a private per-formulation builder. | Deterministic inspection: `grep -n "^def " nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py` returns exactly `register_CFunction_rhs_eval`, and the Dendro BSSN leaf attributes the per-formulation split to no peer. | 09-07-2026 | 09-07-2026 | Found by the inconsistencies survey as K2. That survey also reported the Dendro fCCZ4 leaf as stating the opposite correctly; it does not. Its then-current "That split is Dendro's own" sentence (superseded wording) was about the pure-`build_*`-plus-`register_CFunctions_*` pairing, a different split, so this was a leaf-versus-code contradiction rather than a disagreement between two leaves. |
| CONTR-0011 | The Dendro RHS and constraint modules duplicated private builders and registrars per formulation without a host requirement. | resolved | [New Infrastructure Conformance](infrastructures/new-infrastructure-conformance.md), `### Conformance is one-way` and `## Summary` | [rhs_eval.py](../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `build_rhs_eval` and `register_CFunctions_rhs_eval`; [constraints_eval.py](../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py), `build_constraints_eval` and `register_CFunctions_constraints_eval` | Consolidated each pair into one public builder and registrar. Formulation-specific expressions remain separate; shared lowering and registration follow established conventions. | [BSSN Application Wiring](infrastructures/dendro/bssn-application-wiring.md) | The affected page returns to provisional: this conflict is closed, while broader runtime qualification remains explicitly scoped. | Dendro owner; trigger when either module's per-formulation builder or registrar structure changes, when a host requirement for the split is identified, or when the two pairs are consolidated. | Inspection finds no private per-formulation builder or registrar pairs. Both trusted-expression sweeps pass, and registered full functions and Dendro metadata compare exactly before and after at both shipped profiles. | 09-07-2026 | 09-07-2026 | Resolved by the bounded consolidation; pure builder-plus-registrar pairing remains established convention.

### CONTR-0001

Evidence tuple checks: `inspected=pass` (registered source paths and target
naming inspected); `generated=not-run`; `built=not-run`; `run=not-run`;
`result_checked=not-run`. Dimensions: `platform=not-run`;
`tool version=not-run`; `backend=pass` (BHaH library-generation source path);
`precision=not-run`; `GPU=not-run`; `restart=not-run`;
`distributed=not-run`; `error path=not-run` (stale stdout behavior predicted
from source only); `options=pass` (`create_lib=True` inspected);
`date=pass` (`07-12-2026` inspection). No compiler, generated output, MANGA
consumer, or error path was exercised.

### CONTR-0002

Evidence tuple checks: `inspected=pass` (names `14`, descriptions `14`, dtypes
`13`, defaults `13`, and `zip()` consumer inspected); `generated=not-run`;
`built=not-run`; `run=not-run`; `result_checked=not-run`. Dimensions:
`platform=not-run`; `tool version=not-run`; `backend=pass` (JAX/Commondata
source path); `precision=not-run`; `GPU=not-run`; `restart=not-run`;
`distributed=not-run`; `error path=not-run` (generated call not exercised);
`options=pass` (batch-registration lists inspected); `date=pass`
(`07-12-2026` inspection). Workflow configuration proves generation job shape,
not generated-package execution or a latest successful run.

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

Inspection found `fail-under=10` in both configs, `9.91` in the local wrapper,
and `9.5` in the workflow. None classifies paths from a base revision or checks
legacy score regression. Pylint 4.0.6 gave both base and proposed versions of
the two touched legacy files `10.00/10.00`; the wrapper, workflow, and
classification behavior were not run.

### CONTR-0004

Mewes et al. define
`Lambdatilde^i=DeltaGamma^i+C^i` and
`partial_0=partial_t-L_beta`. The full contravariant-vector Lie pair obeys

```text
beta^k partial_k Lambdatilde^i - Lambdatilde^k partial_k beta^i
= beta^k Dhat_k Lambdatilde^i - Lambdatilde^k Dhat_k beta^i.
```

Decomposing the covariant form of this pair shows that its
`-C^k Dhat_k beta^i` contribution is already present. The same term in the
printed `kappa3=1` bracket therefore duplicates the constraint-vector stretch.
The companion divergence
term is different: the reused BSSN base contains `DeltaGamma^i` in that slot,
so `+2 C^i Dhat_k beta^k/3` is required to obtain the evolved
`Lambdatilde^i` coefficient.

NRPy corrects the apparent paper error. Its aggregate Cartesian coefficient
for `partial_b beta^a` is exactly
`-Lambdatilde^b delta^i_a+(2/3)Lambdatilde^i delta^a_b`. That property was
settled by the derivation and review recorded here. The equation owner now
locks the complete reviewed expression set through trusted dictionaries; it
does not add a Mewes-specific reconstruction, coefficient, or shear test.

Claim evidence:
- Claim: Under Mewes et al.'s own full-Lie time-operator and evolved-connection definitions, its printed `kappa3=1` bracket duplicates `-C^k Dhat_k beta^i`; NRPy corrects that apparent paper error by omitting the duplicate while retaining `+2 C^i Dhat_k beta^k/3`, producing exactly one evolved-vector stretch.
- Role: public/scientific contract
- Deciding authority: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), its `partial_0` definition, Eq. (29), and the evolution system following Eq. (32), composed with [`BSSN_RHSs.py`](../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__`, [`fCCZ4_constraints.py`](../nrpy/equations/general_relativity/fCCZ4_constraints.py), `FCCZ4Constraints.__init__`, and [`fCCZ4_RHSs.py`](../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__` and module `__main__`
- Corroboration: [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (19), corroborates the full corrected Cartesian coefficient; [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), Eqs. (2.11) and (2.17), corroborates one stretch, while spherical Eq. (3.14) contains the compatible divergence promotion but is not identical to the general display off constraint
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=first-principles derivation review followed by 30-significant-digit deterministic trusted sampling of complete expression dictionaries; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=five fCCZ4 evolution variants and twelve gauge variants, including four unchanged frozen-shift cases; date=08-27-2026`

### CONTR-0005

NRPy uses `W=exp(-2 phi)` and `chi=exp(-4 phi)`. These definitions agree with
Mewes et al.'s displayed conformal-factor equations; the neighboring
sign-reversed prose is not part of NRPy's scientific contract.

Claim evidence:
- Claim: NRPy resolves the Mewes conformal-factor prose/equation mismatch by retaining `W=exp(-2 phi)` and `chi=exp(-4 phi)`, matching the displayed equations and established BSSN variables.
- Role: public/scientific contract
- Deciding authority: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (33)-(34); [`BSSN_quantities.py`](../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__`
- Corroboration: [`fCCZ4_constraints.py`](../nrpy/equations/general_relativity/fCCZ4_constraints.py), module `__main__`, exercises `W`, `phi`, and `chi` across Cartesian and SinhSpherical owner cases
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=exact definitions plus deterministic trusted sampling; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=W, phi, and chi across Cartesian and SinhSpherical owner cases; date=08-26-2026`

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
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=exact 4 alpha Theta correction plus deterministic trusted sampling; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=1+log with all three supported shift branches in Cartesian and SinhSpherical coordinates; date=08-26-2026`

### CONTR-0007

The frozen commissioned source interpreted “hyperbolize” as introducing an
auxiliary relaxation covector. The later explicit user clarification instead
defines the intended analogy to CAHD: keep the direct parabolic
momentum-gradient adjustment and multiply its coefficient by the local
timestep scale. Current code therefore adds no evolved field and leaves
initial data, restart, boundary, AMR, and dissipation state routes unchanged.

Claim evidence:
- Claim: Current YBS-MOM uses the direct covariant STF gradient of the matter-complete lower momentum residual, multiplied by `C_YBS_mom*CFL_FACTOR*DSMINGF`; the frozen auxiliary-relaxation proposal is superseded and no cleaner state is permitted.
- Role: public/scientific contract
- Deciding authority: commissioned user clarification composed with [`BSSN_RHSs.py`](../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS momentum branch
- Corroboration: [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (56), and [Etienne, Phys. Rev. D 110, 064045](https://doi.org/10.1103/PhysRevD.110.064045), Eq. (26)
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04 x86_64; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction plus BHaH OpenMP/CUDA source generation; precision=exact symbolic checks plus 30-significant-digit deterministic trusted sampling; GPU=generation pass, execution not-run; restart=not-applicable because no new state; distributed=not-run; error_path=pass for unsupported non-fisheye GeneralRFM spacing; options=both YBS flags enabled jointly in 6 BSSN, 6 fCCZ4, and 8 BHaH rhs_eval existing trusted cases; date=08-30-2026`

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
- Claim: Dendro-GR's `bssn_constraints.cpp` implements the algebraic det/trace enforcement, not the Hamiltonian/momentum constraint diagnostics; NRPy-emitted Dendro kernels are named `<stem>_constraints_eval` and `<stem>_enforce_detgbar_equals_detghat_trAzero`.
- Role: descriptive behavior
- Deciding authority: Dendro-GR `BSSN_GR/include/bssn_constraints.h`, `enforce_bssn_constraints`; `BSSN_GR/src/physcon.cpp`, `physical_constraints`; [BHaH enforce_detgbar_equals_detghat_trAzero.py](../nrpy/infrastructures/BHaH/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `name`
- Corroboration: [constraints_eval.py](../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py), `CONSTRAINTS_EVAL_ALL_BLOCKS_SUFFIX`; [enforce_detgbar_equals_detghat_trAzero.py](../nrpy/infrastructures/Dendro/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `ENFORCE_DETGBAR_EQUALS_DETGHAT_TRAZERO_ALL_BLOCKS_SUFFIX`
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3, GCC 13.3.0, CMake 3.28.3, OpenMPI 4.1.6; backend=Dendro; precision=double; GPU=not-applicable; restart=not-applicable; distributed=1 and 2 MPI ranks; error_path=not-run; options=--fd-order 4 --no-ko; date=09-06-2026`

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
- Claim: the Dendro BSSN constraint builder suppresses the `M` and `LAMBDA_CONSTRAINT` registrations by setting `register_M_and_LAMBDA_CONSTRAINT_gridfunctions` to `False` across the `BSSN_constraints` construction and restoring the previous value afterwards; no Dendro builder deletes an entry from `gri.glb_gridfcs_dict`, and the subpackage's only production-code clear of that registry is the whole-registry sweep reset in `trusted_capture.reset_generation_state`, which no production path reaches; that registry's remaining `clear()` mentions in the subpackage are doctest fixtures.
- Role: descriptive behavior
- Deciding authority: [constraints_eval.py](../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py), the gate save/set/restore around the `BSSN_constraints` construction in the BSSN builder
- Corroboration: [BSSN_constraints.py](../nrpy/equations/general_relativity/BSSN_constraints.py), the `group="AUX"` registrations gated on that CodeParameter
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3; backend=Dendro; precision=not-run; GPU=not-applicable; restart=not-applicable; distributed=not-run; error_path=not-run; options=not-run; date=09-07-2026`

### CONTR-0010

The Dendro BSSN leaf attributed the private-builder-per-formulation layout to
BHaH. BHaH's `general_relativity/rhs_eval.py` defines exactly one function,
the public `register_CFunction_rhs_eval`, and branches inline on `enable_fCCZ4`
at two sites. What Dendro does take from BHaH is the module layout -- one
module per artifact, both formulations behind one boolean -- and that half of
the sentence stands; the private builder and private registrar per formulation
inside the module were Dendro's own. CONTR-0011 subsequently removed that
duplication; this entry preserves the attribution correction.

Claim evidence:
- Claim: `BHaH/general_relativity/rhs_eval.py` holds one public registration function and no private per-formulation builder or registrar; the former private-builder-per-formulation split inside a Dendro module was Dendro's own, not BHaH's arrangement.
- Role: descriptive behavior
- Deciding authority: [BHaH rhs_eval.py](../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval` and its two `enable_fCCZ4` branches
- Corroboration: none available; BHaH's own module is the whole of the evidence, and no Dendro leaf records the per-formulation split as a peer practice
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=Ubuntu 24.04; tool_version=Python 3.12.3; backend=Dendro; precision=not-run; GPU=not-applicable; restart=not-applicable; distributed=not-run; error_path=not-run; options=not-run; date=09-07-2026`

### CONTR-0011

Resolved on 09-07-2026. Both modules now contain one public builder and one
public registrar. Only expression assembly and formulation-specific validation
branch on `enable_fCCZ4`; lowering, wrappers, metadata, and registration share
one path. The previous duplication had no host requirement and was removed.

The pure `build_*` plus `register_CFunctions_*` pairing is unchanged: established
infrastructures already use both forms. This row concerned duplication within
each pair, not the separation between building and registration.

Claim evidence:
- Claim: the unsupported per-formulation builder and registrar duplication is removed in both Dendro modules, resolving CONTR-0011.
- Role: normative rule
- Deciding authority: [new-infrastructure-conformance.md](infrastructures/new-infrastructure-conformance.md), `Conformance is one-way`; [rhs_eval.py](../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `build_rhs_eval` and `register_CFunctions_rhs_eval`; [constraints_eval.py](../nrpy/infrastructures/Dendro/general_relativity/constraints_eval.py), `build_constraints_eval` and `register_CFunctions_constraints_eval`
- Corroboration: [BHaH rhs_eval.py](../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), inline formulation branching; both Dendro trusted-expression sweeps and direct before/after full-function and metadata comparisons
