# Contradictions

> Register for contested and stale KB claims. Plain Markdown only. · Status: confirmed

Known contested/stale claims are tracked below. A row records
source-side truth and containment; it does not imply that NRPy source was fixed.

## Register

| ID | Claim | Claim status | Source A | Source B | Authority decision | Affected pages | Page-status rationale | Owner/trigger | Resolution test | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CONTR-0001 | `manga_bhah_lib` tells users to run `./bhah_lib` and says `bhah_lib.par` exists, but its assembly path requests a library target and does not call the parser/default-parfile writer. | stale | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) final two `print()` calls | [`manga_bhah_lib.py`](../nrpy/examples/manga_bhah_lib.py) commented parser registration and `output_CFunctions_function_prototypes_and_construct_Makefile(create_lib=True)` call; [`Makefile_helpers.py`](../nrpy/infrastructures/BHaH/Makefile_helpers.py) `output_CFunctions_function_prototypes_and_construct_Makefile` | For descriptive current behavior, precise registrations and writer calls decide. They emit `libbhah_lib.so` on Linux or `libbhah_lib.dylib` on Darwin, not an executable `bhah_lib`; no invoked writer emits `bhah_lib.par`. Both generic final prints are stale. | [BHaH Lifecycle And Project Assembly](infrastructures/bhah/lifecycle-and-project-assembly.md); [Matter TOV Workflows](examples/matter-tov-workflows.md) | Both pages remain `confirmed`: these are bounded final-message claims, corrected inline with explicit non-guarantees; their principal lifecycle/workflow answers remain supported. | BHaH/example owners; trigger when `manga_bhah_lib.py` final messaging, parser/parfile registration, or Makefile arguments change. | In an isolated output root, run the generator, capture stdout, inspect emitted files and the Makefile `all` target, and build. Pass only if stdout gives library-consumer guidance, never instructs executing `./bhah_lib`, mentions a parfile only if one is emitted, and `make` produces the platform library named by the Makefile. | Configured CI does not cover this resolution test. |
| CONTR-0003 | Pylint policy distinguishes grandfathered existing tracked handwritten files from newly added handwritten files, but both configs use `fail-under=10`, the local wrapper uses a blanket `9.91` floor, and the workflow uses a blanket `9.5` floor; none enforces the tiered rule. | stale | [`coding_style.md`](../coding_style.md) `## Static Analysis Configuration` | [`.pylintrc`](../.pylintrc) and [`.pylintrc_python36`](../.pylintrc_python36) `[MASTER]`; [`single_file_static_analysis.sh`](../.github/single_file_static_analysis.sh) `run_test_step`; [`main.yml`](../.github/workflows/main.yml) `static-analysis` | Existing tracked handwritten files retain their pre-change score without regression, including legacy scores at or below `9.5`; newly added handwritten files require `10.00/10.00`. Enforcement must classify the path from the base revision, compare legacy base and proposed ratings under the same tool/config, require exact `10.00` for new paths, and reject missing or unparseable ratings. | [Static Analysis](validation/static-analysis.md) | The page remains `confirmed`: the enforcement mismatch is one bounded, directly sourced Detail claim with an inline marker and explicit non-guarantee; the principal command and configuration answer remains reliable. | Validation and contributor-guidance owners; trigger when the tiered Pylint rule, either config, local wrapper classification/rating handling, or workflow classification/rating handling changes. | Deterministically inspect contributor guidance, both configs, local wrapper, and workflow. Pass only if documentation states the tiered rule and every enforcement path classifies paths from a base revision, rejects a lower legacy rating, requires exactly `10.00` for a new handwritten file, and rejects missing or unparseable ratings; generated trusted data remains exempt. | Configs and wrappers do not implement the required classification. |
| CONTR-0012 | Python codegen and ordinary scalar C codegen expose substitution options, but a surviving postprocessed CSE temporary makes a nonempty substitution dictionary raise `TypeError`. | contested | [py_codegen.py](../nrpy/py_codegen.py), `PyCodeGen` and `py_codegen` option and returned-code descriptions; [c_codegen.py](../nrpy/c_codegen.py), `CCodeGen` and `c_codegen` substitution option | [py_codegen.py](../nrpy/py_codegen.py) and [c_codegen.py](../nrpy/c_codegen.py), assignment to `common_subexpression[1]`; [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `cse_postprocess` tuple entries | The implementation decides current behavior; tuple assignment prevents code return for these combinations. The KB documents the failure without treating it as successful substitution. | [Python Codegen](core/python-codegen.md); [C Codegen](core/c-codegen.md) | Both pages are `contested` because a public option combination fails. | C and Python codegen sources; recheck when substitution handling or CSE replacement representation changes. | For Python and ordinary non-FD, non-SIMD C generation, retain a common subexpression and supply a nonempty substitution dictionary with matching and nonmatching keys; require code return and agreement with the intended substituted expressions. Retain no-CSE and no-temporary cases. | This claim does not extend the scalar C failure to SIMD or FD paths. |
| CONTR-0014 | CSE sorting can preserve input-list order for tied structural keys, because the structural key sorts the arguments of commutative-valued expressions such as powers. | contested | [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), documented helper behavior | [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `sort_cse_output_deterministically` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [CSE And Printer Support](core/helpers/cse-and-printer-support.md); [C Codegen](core/c-codegen.md) | Pages are `contested`: output ordering requires a qualification. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Permute replacements for `x**y` and `y**x` and require equal renamed output. | Source limitation remains; this KB correction does not change the implementation. |
| CONTR-0015 | SIMD symbolic rewrites can emit Pow/Cbrt/Log/Sign names missing from the header. | contested | [simd.py](../nrpy/helpers/simd.py), `expr_convert_to_simd_intrins` | [simd_intrinsics.h](../nrpy/helpers/simd_intrinsics.h), instruction-set definitions | The implementation decides current behavior; the header defines no `PowSIMD`, `CbrtSIMD`, `LogSIMD`, or `SignSIMD`, so a general supported-expression guarantee is not established. | [SIMD And Intrinsic Support](core/helpers/simd-and-intrinsic-support.md) | Page is `contested`: missing definitions affect generated calculations. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Compile supported Pow/Cbrt/Log/Sign expressions. | Including simd_intrinsics.h alone does not supply the missing SIMD names. |
| CONTR-0019 | The CUDA wrapper ignores cfunc_type and may leave rewritten parameter accesses without a declared streamid. | contested | [utilities.py](../nrpy/helpers/parallelization/utilities.py), documented helper behavior | [utilities.py](../nrpy/helpers/parallelization/utilities.py), `generate_kernel_and_launch_code`; [gpu_kernel.py](../nrpy/helpers/parallelization/gpu_kernel.py), `GPUKernel.__init__` | The implementation decides current behavior; the documented general guarantee is not established for these inputs. | [Loop Kernel And Device Helpers](core/helpers/loop-kernel-and-device-helpers.md) | Page is `contested`: public wrapper inputs can produce a different return type or an undeclared `streamid`. | Helper sources; recheck when the named implementation or its documented guarantee changes. | Generate and compile a CUDA helper with a requested non-void return type and custom CUDA launch settings with parameter reads; require the requested return type and a declared streamid. | Source limitations remain; this KB correction does not change the implementation. |
| CONTR-0020 | Coordinate-bound helper docstrings promise Boolean comparisons, but shifted quotients yield fractional values, shifted thresholds, and singular denominators; the SEOBNR ringdown frequency inherits those weights. | contested | [Min_Max_and_Piecewise_Expressions.py](../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py), coordinate-helper docstrings | [Min_Max_and_Piecewise_Expressions.py](../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py), coordinate-helper return expressions | Return expressions decide current behavior; exact comparison semantics are not guaranteed. | [GRHD](equations/grhd.md); [Precessing Rotations And Ringdown](equations/seobnr/precessing-rotations-and-ringdown.md) | GRHD remains `confirmed`: one bounded helper claim is qualified inline. The SEOBNR page is `contested` because the frequency, Euler angles, and QNM shifts inherit the discrepancy. | GRHD comparison helpers and SEOBNR ringdown constructor; recheck on changes to offsets, formulas, or documented comparison semantics. | Interpret nrpyAbs as absolute value and compare with true predicates on both sides of zero and at plus/minus TINYDOUBLE; require finite 0-or-1 results matching the predicates or an explicitly restricted interface. Verify the corresponding ringdown frequency weights and downstream angle/QNM expressions. | Trusted samples replacing nrpyAbs with sin, cos, or exp do not decide these semantics. |
| CONTR-0021 | Wave constructors can leave required gridfunctions or CodeParameters unregistered when existing names trigger incomplete or incorrect guards. | contested | [WaveEquation_RHSs.py](../nrpy/equations/wave_equation/WaveEquation_RHSs.py), constructor registration comments; [WaveEquation_Solutions_InitialData.py](../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py), constructor registration comments | [WaveEquation_RHSs.py](../nrpy/equations/wave_equation/WaveEquation_RHSs.py), [WaveEquationCurvilinear_RHSs.py](../nrpy/equations/wave_equation/WaveEquationCurvilinear_RHSs.py), and [WaveEquation_Solutions_InitialData.py](../nrpy/equations/wave_equation/WaveEquation_Solutions_InitialData.py), constructors; [params.py](../nrpy/params.py), CodeParameter registration | Actual registry checks decide current behavior; callers must provide compatible entries when reusing state. | [Wave Equation](equations/wave-equation.md) | Page is `contested`: independent gridfunction and parameter-registration conditions invalidate unconditional setup claims. | Wave equation constructors; recheck when registry guards or registration validation changes. | In separate isolated processes, pre-register uu without vv, or with incompatible metadata, and pre-register NRPy parameters named wavespeed or time; require compatible gridfunction and CodeParameter entries or explicit rejection. Keep a clean-registry case. | Curvilinear RHS wavespeed registration is unconditional; its uu-only gridfunction guard still applies. |
| CONTR-0022 | A zero Lorentz boost should be the identity, but direct zero-velocity construction produces undefined spatial entries. | contested | [LorentzBoost.py](../nrpy/equations/general_relativity/LorentzBoost.py), documented Lorentz-transformation purpose | [LorentzBoost.py](../nrpy/equations/general_relativity/LorentzBoost.py), `_compute_LorentzMatrix` spatial block | The implemented division by squared speed decides current behavior; zero input is not handled by its identity limit. | [Metric Conversions And Matter](equations/general-relativity/metric-conversions-and-matter.md) | Page remains `confirmed`: this is one bounded helper defect, qualified inline; other metric and matter conversions remain valid. | LorentzBoost constructor; recheck when speed-domain handling or the spatial matrix formula changes. | Construct the direct zero boost and its inverse; require finite identity matrices and unchanged transformed tensors, retaining nonzero subluminal cases. | This KB correction documents the source defect without changing Python. |

### CONTR-0001

The resolution test deliberately requires generation, artifact inspection, and
a build because source inspection alone cannot establish the user-facing
library contract.

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

### CONTR-0014

The structural key sorts the arguments of any commutative-valued expression, so
`x**y` and `y**x` can have equal keys, and ties retain input order.

Claim evidence:
- Claim: CSE sorting can preserve input-list order for tied structural keys, because the structural key sorts the arguments of commutative-valued expressions such as powers.
- Role: descriptive behavior
- Deciding authority: [cse_preprocess_postprocess.py](../nrpy/helpers/cse_preprocess_postprocess.py), `sort_cse_output_deterministically`
- Corroboration: `none available`; the stated counterexamples follow from the implementation, while its docstring or paired symbolic definition supplies the intended behavior.

### CONTR-0015

`expr_convert_to_simd_intrins()` can emit `PowSIMD`, `CbrtSIMD`, `LogSIMD`, and
`SignSIMD`, but `simd_intrinsics.h` defines none of them.

Claim evidence:
- Claim: Emitted Pow/Cbrt/Log/Sign SIMD names lack definitions in the SIMD header.
- Role: descriptive behavior
- Deciding authority: [simd.py](../nrpy/helpers/simd.py), `expr_convert_to_simd_intrins`; [simd_intrinsics.h](../nrpy/helpers/simd_intrinsics.h), instruction-set definitions
- Corroboration: `none available`; the stated counterexamples follow from the implementation, while its docstring or paired symbolic definition supplies the intended behavior.

### CONTR-0019

`generate_kernel_and_launch_code()` ignores `cfunc_type` on the CUDA path, and
custom launch settings without a `stream` entry leave rewritten
`d_params[streamid].` accesses without a declared `streamid`.

Claim evidence:
- Claim: The CUDA wrapper ignores cfunc_type and may leave rewritten parameter accesses without a declared streamid.
- Role: descriptive behavior
- Deciding authority: [utilities.py](../nrpy/helpers/parallelization/utilities.py), `generate_kernel_and_launch_code`; [gpu_kernel.py](../nrpy/helpers/parallelization/gpu_kernel.py), `GPUKernel.__init__`
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

## Rules

- Follow the contradiction contract in [Schema](SCHEMA.md#contradiction-contract).
- Open a row when registered sources disagree, when a living source supersedes
  a KB claim, or when a page is marked `contested` or `stale`.
- Update affected-page links and exact inline markers together.
- Delete a row and its subsection once its resolution test passes and all
  affected pages, reverse dependents, aliases, typed neighbors, and targeted
  wiki hits are reconciled, and remove every pointer to its ID.
