# Contradictions

> Register for contested and stale KB claims. Plain Markdown only. · Status: confirmed · Last reconciled: 07-13-2026

Known contested/stale claims as of 07-13-2026 are tracked below. A row records
source-side truth and containment; it does not imply that NRPy source was fixed.

## Register

| ID | Claim | Claim status | Source A | Source B | Authority decision | Affected pages | Page-status rationale | Owner/trigger | Resolution test | Opened | Resolved | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CONTR-0001 | `manga_bhah_lib` tells users to run `./bhah_lib` and says `bhah_lib.par` exists, but its assembly path requests a library target and does not call the parser/default-parfile writer. | stale | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) final two `print()` calls | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) commented parser registration and `output_CFunctions_function_prototypes_and_construct_Makefile(create_lib=True)` call; [`Makefile_helpers.py`](../nrpy/infrastructures/BHaH/Makefile_helpers.py) `output_CFunctions_function_prototypes_and_construct_Makefile` | For descriptive current behavior, precise registrations and writer calls decide. They emit `libbhah_lib.so` on Linux or `libbhah_lib.dylib` on Darwin, not an executable `bhah_lib`; no invoked writer emits `bhah_lib.par`. Both generic final prints are stale. | [BHaH Lifecycle And Project Assembly](infrastructures/bhah/lifecycle-and-project-assembly.md); [Matter TOV Workflows](examples/matter-tov-workflows.md) | Both pages remain `confirmed`: these are bounded final-message claims, corrected inline with explicit non-guarantees; their principal lifecycle/workflow answers remain supported. | BHaH/example owners; trigger when `manga_bhah_lib.py` final messaging, parser/parfile registration, or Makefile arguments change. | In an isolated output root, run the generator, capture stdout, inspect emitted files and the Makefile `all` target, and build. Pass only if stdout gives library-consumer guidance, never instructs executing `./bhah_lib`, mentions a parfile only if one is emitted, and `make` produces the platform library named by the Makefile. | 06-30-2026 | - | Current CI commands for this generator are commented out, so generation/build results are `not-run` here. |
| CONTR-0002 | `sebobv1_jax` emits `Commondata(..., a_f=a_f)`, but its batch field registration supplies fourteen names/descriptions and thirteen dtypes/defaults; `zip()` therefore omits `a_f` from the generated dataclass. | contested | [`SEOBNRv5_aligned_spin_coefficients.py`](../nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py) `register_PyFunction_SEOBNRv5_aligned_spin_coefficients` emitted return | [`sebobv1_jax.py`](../nrpy/examples/sebobv1_jax.py) `register_commondata_params` call and [`commondata.py`](../nrpy/infrastructures/JAX/commondata.py) `register_commondata_params` / `generate_commondata_dataclass` | Descriptive behavior is decided by all three code paths together: Python `zip()` truncation registers thirteen fields, dataclass emission follows that registry, and the emitted function passes unsupported keyword `a_f`. Project generation alone does not establish callable coefficient initialization. | [SEBOBv1 JAX Workflow](infrastructures/jax/sebobv1-jax-workflow.md); [Commondata And PyFunction Registry](infrastructures/jax/commondata-and-pyfunction-registry.md); [Waveform JAX PN Generators](examples/waveform-jax-pn-generators.md); [Example Generator Catalog](examples/example-generator-catalog.md) | `SEBOBv1 JAX Workflow` is `contested` because the mismatch invalidates its sole generated coefficient interface and changes its Summary-level routed answer. Other affected pages remain `confirmed` only where they describe this as a bounded noncentral detail and make no callable-function guarantee. | JAX/example owners; trigger when batch registration gains length validation, the four lists become equal, or the emitted return signature changes. | In an isolated output root, generate the project; assert batch metadata passes equal-length validation; import the generated coefficient function; call it with a valid finite input set; assert every returned keyword is a generated `Commondata` field and every returned numerical field is finite; if `a_f` remains in the intended return contract, assert that `a_f` is present and finite; then run the generated smoke test. All applicable checks must pass. | 07-06-2026 | - | Workflow config invokes generation only; generated-package install, import, coefficient call, and result checks are `not-run` here. |
| CONTR-0003 | Pylint policy distinguishes grandfathered existing tracked handwritten files from newly added handwritten files, but both configs use `fail-under=10`, the local wrapper uses a blanket `9.91` floor, and the workflow uses a blanket `9.5` floor; none enforces the tiered rule. | stale | [`coding_style.md`](../coding_style.md) `## Static Analysis Configuration` | [`.pylintrc`](../.pylintrc) and [`.pylintrc_python36`](../.pylintrc_python36) `[MASTER]`; [`single_file_static_analysis.sh`](../.github/single_file_static_analysis.sh) `run_test_step`; [`main.yml`](../.github/workflows/main.yml) `static-analysis` | Existing tracked handwritten files retain their pre-change score without regression, including legacy scores at or below `9.5`; newly added handwritten files require `10.00/10.00`. Enforcement must classify the path from the base revision, compare legacy base and proposed ratings under the same tool/config, require exact `10.00` for new paths, and reject missing or unparseable ratings. | [Static Analysis](validation/static-analysis.md) | The page remains `confirmed`: the enforcement mismatch is one bounded, directly sourced Detail claim with an inline marker and explicit non-guarantee; the principal command and configuration answer remains reliable. | Validation and contributor-guidance owners; trigger when the tiered Pylint rule, either config, local wrapper classification/rating handling, or workflow classification/rating handling changes. | Deterministically inspect contributor guidance, both configs, local wrapper, and workflow. Pass only if documentation states the tiered rule and every enforcement path classifies paths from a base revision, rejects a lower legacy rating, requires exactly `10.00` for a new handwritten file, and rejects missing or unparseable ratings; generated trusted data remains exempt. | 07-13-2026 | - | Pylint 4.0.6 gave both base and proposed versions of the two touched legacy files `10.00/10.00`; wrapper and workflow classification/enforcement remain `not-run`. |

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

## Rules

- Follow the contradiction contract in [Schema](SCHEMA.md#contradiction-contract).
- Open a row when registered sources disagree, when a living source supersedes
  a KB claim, or when a page is marked `contested` or `stale`.
- Update affected-page links and exact inline markers together.
- Resolve a row only after its resolution test passes and all affected pages,
  reverse dependents, aliases, typed neighbors, and targeted wiki hits are
  reconciled.
