# Contradictions

> Register for contested and stale KB claims. Plain Markdown only. · Status: confirmed · Last reconciled: 08-26-2026

Known contested/stale claims as of 08-26-2026 are tracked below. A row records
source-side truth and containment; it does not imply that NRPy source was fixed.

## Register

| ID | Claim | Claim status | Source A | Source B | Authority decision | Affected pages | Page-status rationale | Owner/trigger | Resolution test | Opened | Resolved | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CONTR-0001 | `manga_bhah_lib` tells users to run `./bhah_lib` and says `bhah_lib.par` exists, but its assembly path requests a library target and does not call the parser/default-parfile writer. | stale | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) final two `print()` calls | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) commented parser registration and `output_CFunctions_function_prototypes_and_construct_Makefile(create_lib=True)` call; [`Makefile_helpers.py`](../nrpy/infrastructures/BHaH/Makefile_helpers.py) `output_CFunctions_function_prototypes_and_construct_Makefile` | For descriptive current behavior, precise registrations and writer calls decide. They emit `libbhah_lib.so` on Linux or `libbhah_lib.dylib` on Darwin, not an executable `bhah_lib`; no invoked writer emits `bhah_lib.par`. Both generic final prints are stale. | [BHaH Lifecycle And Project Assembly](infrastructures/bhah/lifecycle-and-project-assembly.md); [Matter TOV Workflows](examples/matter-tov-workflows.md) | Both pages remain `confirmed`: these are bounded final-message claims, corrected inline with explicit non-guarantees; their principal lifecycle/workflow answers remain supported. | BHaH/example owners; trigger when `manga_bhah_lib.py` final messaging, parser/parfile registration, or Makefile arguments change. | In an isolated output root, run the generator, capture stdout, inspect emitted files and the Makefile `all` target, and build. Pass only if stdout gives library-consumer guidance, never instructs executing `./bhah_lib`, mentions a parfile only if one is emitted, and `make` produces the platform library named by the Makefile. | 06-30-2026 | - | Current CI commands for this generator are commented out, so generation/build results are `not-run` here. |
| CONTR-0002 | `sebobv1_jax` emits `Commondata(..., a_f=a_f)`, but its batch field registration supplies fourteen names/descriptions and thirteen dtypes/defaults; `zip()` therefore omits `a_f` from the generated dataclass. | contested | [`SEOBNRv5_aligned_spin_coefficients.py`](../nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py) `register_PyFunction_SEOBNRv5_aligned_spin_coefficients` emitted return | [`sebobv1_jax.py`](../nrpy/examples/sebobv1_jax.py) `register_commondata_params` call and [`commondata.py`](../nrpy/infrastructures/JAX/commondata.py) `register_commondata_params` / `generate_commondata_dataclass` | Descriptive behavior is decided by all three code paths together: Python `zip()` truncation registers thirteen fields, dataclass emission follows that registry, and the emitted function passes unsupported keyword `a_f`. Project generation alone does not establish callable coefficient initialization. | [SEBOBv1 JAX Workflow](infrastructures/jax/sebobv1-jax-workflow.md); [Commondata And PyFunction Registry](infrastructures/jax/commondata-and-pyfunction-registry.md); [Waveform JAX PN Generators](examples/waveform-jax-pn-generators.md); [Example Generator Catalog](examples/example-generator-catalog.md) | `SEBOBv1 JAX Workflow` is `contested` because the mismatch invalidates its sole generated coefficient interface and changes its Summary-level routed answer. Other affected pages remain `confirmed` only where they describe this as a bounded noncentral detail and make no callable-function guarantee. | JAX/example owners; trigger when batch registration gains length validation, the four lists become equal, or the emitted return signature changes. | In an isolated output root, generate the project; assert batch metadata passes equal-length validation; import the generated coefficient function; call it with a valid finite input set; assert every returned keyword is a generated `Commondata` field and every returned numerical field is finite; if `a_f` remains in the intended return contract, assert that `a_f` is present and finite; then run the generated smoke test. All applicable checks must pass. | 07-06-2026 | - | Workflow config invokes generation only; generated-package install, import, coefficient call, and result checks are `not-run` here. |
| CONTR-0003 | Pylint policy distinguishes grandfathered existing tracked handwritten files from newly added handwritten files, but both configs use `fail-under=10`, the local wrapper uses a blanket `9.91` floor, and the workflow uses a blanket `9.5` floor; none enforces the tiered rule. | stale | [`coding_style.md`](../coding_style.md) `## Static Analysis Configuration` | [`.pylintrc`](../.pylintrc) and [`.pylintrc_python36`](../.pylintrc_python36) `[MASTER]`; [`single_file_static_analysis.sh`](../.github/single_file_static_analysis.sh) `run_test_step`; [`main.yml`](../.github/workflows/main.yml) `static-analysis` | Existing tracked handwritten files retain their pre-change score without regression, including legacy scores at or below `9.5`; newly added handwritten files require `10.00/10.00`. Enforcement must classify the path from the base revision, compare legacy base and proposed ratings under the same tool/config, require exact `10.00` for new paths, and reject missing or unparseable ratings. | [Static Analysis](validation/static-analysis.md) | The page remains `confirmed`: the enforcement mismatch is one bounded, directly sourced Detail claim with an inline marker and explicit non-guarantee; the principal command and configuration answer remains reliable. | Validation and contributor-guidance owners; trigger when the tiered Pylint rule, either config, local wrapper classification/rating handling, or workflow classification/rating handling changes. | Deterministically inspect contributor guidance, both configs, local wrapper, and workflow. Pass only if documentation states the tiered rule and every enforcement path classifies paths from a base revision, rejects a lower legacy rating, requires exactly `10.00` for a new handwritten file, and rejects missing or unparseable ratings; generated trusted data remains exempt. | 07-13-2026 | - | Pylint 4.0.6 gave both base and proposed versions of the two touched legacy files `10.00/10.00`; wrapper and workflow classification/enforcement remain `not-run`. |
| CONTR-0004 | Mewes et al. and Alic et al. display connection shift sectors that agree on the connection-constraint surface but differ off it. | resolved | [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), evolution system after Eq. (32) | [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (19) | Mewes et al. decide NRPy's reference-metric fCCZ4 shift sector. At equal `kappa3=1`, exact Cartesian subtraction gives `S_Mewes^i-S_Alic^i=-C^k partial_k beta^i`; Alic is an on-constraint covariance corroboration, not an off-constraint replacement. | [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md) | Page remains `confirmed`: authority and exact source-level residual are documented, while Cartesian owner validation independently rebuilds both shift deltas relative to the same BSSN base and checks the residual coefficient and selected Mewes coefficients. | fCCZ4 equation owners; trigger when connection-variable conventions or shift-sector terms change. | Document the authority decision and exact source-level residual; verify the Cartesian equal-`kappa3=1` Mewes/Alic crosswalk coefficient and the implemented Mewes shift-gradient coefficients `(2/3) C^i delta^b_a-C^b delta^i_a`; owner trusted-expression validation must pass. | 08-26-2026 | 08-26-2026 | Sanchis-Gual et al.'s general and spherical displayed sectors also differ off constraint and are treated as qualified corroboration. |
| CONTR-0005 | Mewes et al.'s prose around its conformal-factor alternatives reverses the signs implied by its displayed evolution equations and used by NRPy's established BSSN variables. | resolved | [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eqs. (33)-(34) and surrounding prose | [`BSSN_quantities.py`](../nrpy/equations/general_relativity/BSSN_quantities.py), `BSSNQuantities.__init__` | NRPy retains its established `W=exp(-2 phi)` and `chi=exp(-4 phi)` definitions, whose chain rules agree with the displayed evolution equations; sign-reversed prose is not adopted. | [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md) | Page remains `confirmed`: definitions and authority choice are explicit and owner validation covers all three conformal-factor options. | BSSN/fCCZ4 equation owners; trigger when conformal-factor definitions, source crosswalk, or option validation changes. | Run fCCZ4 constraint/RHS owner validation across `W`, `phi`, and `chi`; verify exact chain-rule definitions and trusted expressions. | 08-26-2026 | 08-26-2026 | This resolves a source-internal prose/equation mismatch without changing NRPy's established variables. |
| CONTR-0006 | Mewes et al.'s displayed advective 1+log equation omits the lapse factor present in Alic et al.'s exact advective form; Sanchis-Gual et al. corroborates the lapse factor only in a nonadvective form. | resolved | [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), Eq. (35) | [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (20) | NRPy adopts `partial_0 alpha=-2 alpha(K-2 Theta)` by reusing BSSN 1+log and adding `4 alpha Theta`; the missing lapse in Mewes Eq. (35) is treated as a typographical omission. | [Fully Covariant Conformal Z4](equations/general-relativity/fccz4.md) | Page remains `confirmed`: the adopted gauge formula and source decision are explicit and validated across Cartesian and curved reference metrics. | fCCZ4 gauge owners; trigger when lapse options, 1+log correction, source authority, or gauge validation changes. | Run all twelve fCCZ4 gauge trusted variants and verify the 1+log correction is exactly `4 alpha Theta`. | 08-26-2026 | 08-26-2026 | Sanchis-Gual et al. Eq. (2.25) corroborates the lapse factor but not the advective operator; frozen lapse is unaffected. |

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

Mewes et al. are authoritative for NRPy's fully covariant reference-metric
fCCZ4 evolution system. Direct Cartesian expansion of the published Mewes and
Alic connection shift sectors at equal `kappa3=1` gives
`S_Mewes^i-S_Alic^i=-C^k partial_k beta^i`. This vanishes on `C^i=0` but is not
identically zero off that surface. Cartesian fCCZ4 owner validation rebuilds
both shift deltas relative to the same BSSN base, samples them with unconstrained
symbolic `C^i`, and checks the exact residual and selected Mewes
shift-gradient coefficients. This is not an independent scientific
implementation, explicit field-data fixture, or evolution sample.

Claim evidence:
- Claim: Mewes et al. decide the implemented connection shift sector; at equal `kappa3=1`, direct Cartesian expansion of the two published sectors gives `S_Mewes^i-S_Alic^i=-C^k partial_k beta^i`, so Alic is corroboration only on the connection-constraint surface for this sector.
- Role: public/scientific contract
- Deciding authority: [Mewes et al., arXiv:2002.06225v2](https://arxiv.org/pdf/2002.06225v2), evolution system following Eq. (32); [fCCZ4_RHSs.py](../nrpy/equations/general_relativity/fCCZ4_RHSs.py), `FCCZ4RHSs.__init__`
- Corroboration: [Alic et al., arXiv:1106.2254v2](https://arxiv.org/pdf/1106.2254v2), Eq. (19), supplies the comparison sector; [Sanchis-Gual et al., arXiv:1403.3653v1](https://arxiv.org/pdf/1403.3653v1), Eqs. (2.17) and (3.14), is qualified on-constraint corroboration
- Validation: `inspected=pass; generated=pass; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=SymPy expression construction; precision=Cartesian equal-kappa3 construction of both shift deltas relative to the same BSSN base, deterministic sampling with unconstrained symbolic C^i, and exact differentiation of their residual and the implemented correction with respect to shift-gradient symbols; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=five fCCZ4 evolution owner variants, literal Alic crosswalk limited to Cartesian, explicit off-constraint field-data fixture and evolution sample not-run; date=08-26-2026`

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

## Rules

- Follow the contradiction contract in [Schema](SCHEMA.md#contradiction-contract).
- Open a row when registered sources disagree, when a living source supersedes
  a KB claim, or when a page is marked `contested` or `stale`.
- Update affected-page links and exact inline markers together.
- Resolve a row only after its resolution test passes and all affected pages,
  reverse dependents, aliases, typed neighbors, and targeted wiki hits are
  reconciled.
