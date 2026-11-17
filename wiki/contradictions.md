# Contradictions

> Register of active contested and stale NRPy KB claims. · Status: confirmed

Active source discrepancies and unreconciled drift are tracked below. Recording an entry documents current code versus specification or docstring behavior without implying a source patch.

## Rules

- Follow [KB Rules](SCHEMA.md#contradiction-rules): keep only `contested` and `stale` claims in `## Register`, paired with `### CONTR-XXXX` subsections and inline markers (`Claim status: contested; contradiction: CONTR-0001.` or `stale`) on every affected leaf.
- Look up one entry with `rg -n -F '| CONTR-0014 |' wiki/contradictions.md` and read its `### CONTR-XXXX` subsection.
- Retire resolved entries in the resolving commit once affected leaves are reconciled and all `#contr-*` links are removed. Retired IDs: CONTR-0002, CONTR-0004, CONTR-0005, CONTR-0006, CONTR-0007, CONTR-0008 to CONTR-0011, CONTR-0013.

## Register

| ID | Claim status | Affected pages | Summary |
| --- | --- | --- | --- |
| CONTR-0001 | stale | [Assembly](infrastructures/bhah/lifecycle-and-project-assembly.md); [TOV](examples/matter-tov-workflows.md) | `manga_bhah_lib.py` prints `./bhah_lib`/`.par` while building only `libbhah_lib`. |
| CONTR-0003 | stale | [Static Analysis](validation/static-analysis.md) | `coding_style.md` states tiered Pylint rules; configs/wrappers enforce blanket floors. |
| CONTR-0012 | contested | [Python Codegen](core/python-codegen.md); [C Codegen](core/c-codegen.md) | Nonempty `postproc_substitution_dict` raises `TypeError` on tuple CSE entries. |
| CONTR-0014 | contested | [CSE/Printer](core/helpers/cse-and-printer-support.md); [C Codegen](core/c-codegen.md) | CSE sort keeps input order on ties, and `_Integer_*`/`tmp*` names can capture inputs. |
| CONTR-0015 | contested | [SIMD](core/helpers/simd-and-intrinsic-support.md) | `expr_convert_to_simd_intrins` omits `_Integer_1` or emits missing `Pow`/`Cbrt`/`Log`/`Sign`. |
| CONTR-0016 | contested | [Expression Utils](core/helpers/symbolic-expression-utilities.md) | `_get_free_symbol_names_cached` drops free occurrences matching bound variable names. |
| CONTR-0017 | contested | [Expression Utils](core/helpers/symbolic-expression-utilities.md) | `f2r` appends zeros after scientific-notation exponents, corrupting values. |
| CONTR-0018 | contested | [Expression Utils](core/helpers/symbolic-expression-utilities.md) | `product` flattens iterables, exhausts iterators, and mishandles scalars/`repeat=0`. |
| CONTR-0019 | contested | [Loop/Device](core/helpers/loop-kernel-and-device-helpers.md) | `GPUKernel` corrupts empty-arg calls, ignores `cfunc_type`, and can omit `streamid`. |
| CONTR-0020 | contested | [GRHD](equations/grhd.md); [Ringdown](equations/seobnr/precessing-rotations-and-ringdown.md) | Shifted coordinate-bound quotients have singular denominators affecting `omega_prec`. |
| CONTR-0021 | contested | [Wave Equation](equations/wave-equation.md) | Wave constructors guard only `uu` and check `glb_params_dict` for `CodeParameter`s. |
| CONTR-0022 | contested | [Metric/Matter](equations/general-relativity/metric-conversions-and-matter.md) | Zero-velocity `LorentzBoost` divides by `0/0` instead of returning identity. |
| CONTR-0023 | contested | [Geodesics](infrastructures/bhah/geodesics-and-raytracing-runtime.md) | Azimuthal geodesic interpolation rotates reflected Cartesian records by `+/-pi`. |

### CONTR-0001

- Competing sources & authority: `nrpy/examples/manga_bhah_lib.py` final `print()` lines versus its `output_CFunctions_function_prototypes_and_construct_Makefile(create_lib=True)` call and `Makefile_helpers.py` (invoked writers decide: shared library only, no `bhah_lib.par`). Both affected leaves stay `confirmed`.
- Recheck trigger: edits to `manga_bhah_lib.py` final messages, parameter-file registration, or `Makefile_helpers.py` `create_lib` handling.
- Resolution check: generate in an isolated directory, verify stdout gives library instructions without `./bhah_lib` or ungenerated `bhah_lib.par`, and verify `make` builds `libbhah_lib.{so,dylib}`.

Evidence (descriptive behavior): `nrpy/examples/manga_bhah_lib.py` final `print()` calls vs `create_lib=True` in `BHaH/Makefile_helpers.py::output_CFunctions_function_prototypes_and_construct_Makefile`.

### CONTR-0003

- Competing sources & authority: `coding_style.md` (`## Static Analysis Configuration`: legacy handwritten files keep pre-change score without regression; new handwritten files require `10.00/10.00`) versus `.pylintrc`/`.pylintrc_python36` (`fail-under=10`), `.github/single_file_static_analysis.sh` (`9.91`), and `.github/workflows/main.yml` (`9.5`). `coding_style.md` decides contributor policy; [Static Analysis](validation/static-analysis.md) stays `confirmed`.
- Recheck trigger: edits to `coding_style.md`, `.pylintrc*`, `single_file_static_analysis.sh`, or `main.yml` `static-analysis`.
- Resolution check: inspect all four files and verify base-revision new-vs-legacy classification (`10.00` for new handwritten files, non-regression for legacy files, rejection of unparseable scores).

Evidence (normative rule): `coding_style.md::Static Analysis Configuration`; corroboration: `.pylintrc`, `.github/single_file_static_analysis.sh::run_test_step`, and `.github/workflows/main.yml::static-analysis`.

### CONTR-0012

- Competing sources & authority: `nrpy/py_codegen.py` and `nrpy/c_codegen.py` substitution option docstrings versus their scalar CSE loops assigning into tuple `common_subexpression[1]` returned by `cse_postprocess`. Implementation decides current behavior; both codegen leaves are `contested` (SIMD and FD C paths are unaffected).
- Recheck trigger: edits to `postproc_substitution_dict` handling or `cse_postprocess` entry types.
- Resolution check: run `py_codegen` and non-FD, non-SIMD `c_codegen` with a surviving CSE temporary and nonempty `postproc_substitution_dict` (matching and nonmatching keys), requiring valid substituted code return.

Evidence (descriptive behavior): `nrpy/py_codegen.py::py_codegen` and `nrpy/c_codegen.py::c_codegen` scalar CSE substitution loops; corroboration: `nrpy/helpers/cse_preprocess_postprocess.py::cse_postprocess`.

### CONTR-0014

- Competing sources & authority: `nrpy/helpers/cse_preprocess_postprocess.py` determinism/equivalence docstrings versus `cse_preprocess` and `sort_cse_output_deterministically`. Implementation decides; both affected leaves are `contested`.
- Recheck trigger: edits to `cse_preprocess` or `sort_cse_output_deterministically`.
- Resolution check: permute replacements for `x**y` and `y**x` requiring identical output, and test `2*x + _Integer_2` and `([(a0, x**3)], [a0 + tmp0])` requiring preserved expressions or explicit collision rejection.

Evidence (descriptive behavior): `nrpy/helpers/cse_preprocess_postprocess.py::cse_preprocess` and `sort_cse_output_deterministically`.

### CONTR-0015

- Competing sources & authority: `nrpy/helpers/simd.py` (`expr_convert_to_simd_intrins`) versus `nrpy/helpers/simd_intrinsics.h` macro definitions. Implementation decides; the SIMD leaf is `contested`.
- Recheck trigger: edits to `expr_convert_to_simd_intrins` or `simd_intrinsics.h`.
- Resolution check: verify `_Integer_1` declaration for negative half-powers and compile emitted `Pow`/`Cbrt`/`Log`/`Sign` intrinsics with `simd_intrinsics.h`.

Evidence (descriptive behavior): `nrpy/helpers/simd.py::expr_convert_to_simd_intrins`; corroboration: `nrpy/helpers/simd_intrinsics.h`.

### CONTR-0016

- Competing sources & authority: `nrpy/helpers/expression_utils.py` SymPy `free_symbols` equivalence docstring versus `_get_free_symbol_names_cached` and `get_params_commondata_symbols_from_expr_list`. Implementation decides; leaf is `contested`.
- Recheck trigger: edits to `_get_free_symbol_names_cached` or parameter extraction helpers.
- Resolution check: compare extracted names against `expr.free_symbols` on `Integral(x, (x, 0, x))`, `Lambda`, and `Subs`, verifying registered free parameter names are retained.

Evidence (descriptive behavior): `nrpy/helpers/expression_utils.py::_get_free_symbol_names_cached` and `get_params_commondata_symbols_from_expr_list`.

### CONTR-0017

- Competing sources & authority: `nrpy/helpers/float_to_rational.py` high-precision conversion docstring versus `f2r` string zero-padding on scientific notation. Implementation decides; leaf is `contested`.
- Recheck trigger: edits to `f2r`.
- Resolution check: require equal rationals for scientific (`1e-7`, `1.2e-7`) and fixed-point decimal strings across padding widths.

Evidence (descriptive behavior): `nrpy/helpers/float_to_rational.py::f2r`.

### CONTR-0018

- Competing sources & authority: `nrpy/helpers/functional.py` Cartesian-product docstring versus `product` recursive flattening/iteration. Implementation decides; leaf is `contested`.
- Recheck trigger: edits to `functional.py::product`.
- Resolution check: match `itertools.product` behavior (or reject invalid scalar inputs) across one-shot iterators, tuple-valued elements, single integers, and `repeat=0`.

Evidence (descriptive behavior): `nrpy/helpers/functional.py::product`.

### CONTR-0019

- Competing sources & authority: `nrpy/helpers/parallelization/utilities.py` wrapper docs versus `nrpy/helpers/parallelization/gpu_kernel.py` (`GPUKernel.__init__`, `GPUKernel.c_function_call`). Implementation decides; leaf is `contested`.
- Recheck trigger: edits to `GPUKernel` or `generate_kernel_and_launch_code`.
- Resolution check: generate and compile an empty-argument host call, a non-`void` CUDA kernel, and custom CUDA launch settings reading parameters, verifying valid call syntax, `cfunc_type`, and declared `streamid`.

Evidence (descriptive behavior): `nrpy/helpers/parallelization/gpu_kernel.py::GPUKernel.__init__` and `GPUKernel.c_function_call`; corroboration: `utilities.py::generate_kernel_and_launch_code`.

### CONTR-0020

- Competing sources & authority: `nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py` Boolean coordinate-helper docstrings versus their `TINYDOUBLE`-shifted quotient formulas used in `SEOBNRv5_merger_ringdown.py` (`omega_prec`). Expressions decide; `grhd.md` stays `confirmed` while `precessing-rotations-and-ringdown.md` is `contested`.
- Recheck trigger: edits to `coord_less_bound`/`coord_greater_bound` offsets/formulas or `SEOBNRv5_MergerRingdown_Quantities.__init__`.
- Resolution check: evaluate with `nrpyAbs(x) = |x|` across `0` and `+/-TINYDOUBLE`, requiring finite `{0, 1}` step values (or restricted documented semantics) and consistent `omega_prec` weights.

Evidence (descriptive behavior): `nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py` coordinate helpers; corroboration: `nrpy/equations/seobnr/SEOBNRv5_merger_ringdown.py::SEOBNRv5_MergerRingdown_Quantities.__init__`.

### CONTR-0021

- Competing sources & authority: `nrpy/equations/wave_equation/` constructor comments versus `uu`-only gridfunction guards and `glb_params_dict` checks (while `CodeParameter` entries live in `glb_code_params_dict`). Implementation decides; `wave-equation.md` is `contested`.
- Recheck trigger: edits to wave-equation constructor guards or `CodeParameter` registration.
- Resolution check: in isolated processes, pre-register `uu` without `vv` (or incompatible metadata) and pre-register `NRPyParameter` `wavespeed`/`time`, requiring valid registrations or explicit `ValueError`.

Evidence (descriptive behavior): `WaveEquation_RHSs.py::WaveEquation_RHSs.__init__`, `WaveEquationCurvilinear_RHSs.py::__init__`, and `WaveEquation_Solutions_InitialData.py::__init__`; corroboration: `nrpy/params.py::CodeParameter._register_globally`.

### CONTR-0022

- Competing sources & authority: `nrpy/equations/general_relativity/LorentzBoost.py` Lorentz-transformation docstring versus `_compute_LorentzMatrix` dividing by squared boost speed without a `v = 0` branch. Implementation decides; `metric-conversions-and-matter.md` stays `confirmed`.
- Recheck trigger: edits to `LorentzBoost._compute_LorentzMatrix`.
- Resolution check: evaluate `LorentzBoost` and its inverse at `vU = [0, 0, 0]`, requiring finite `4x4` identity matrices while preserving subluminal nonzero cases.

Evidence (descriptive behavior): `nrpy/equations/general_relativity/LorentzBoost.py::LorentzBoost._compute_LorentzMatrix` (trusted samples test nonzero boost only).

### CONTR-0023

- Competing sources & authority: `azimuthal_symmetry_spatial_lagrange_interpolation.py` docstring claiming reflected nodes need rotation to the reference plane versus `output_raytracing_data.py` storing Cartesian-basis tensor components at each node's own position. Emitted C code decides; `geodesics-and-raytracing-runtime.md` stays `confirmed`.
- Recheck trigger: edits to `azimuthal_symmetry_spatial_lagrange_interpolation.py` reflection/rotation logic or `output_raytracing_data.py` basis conventions.
- Resolution check: on `Spherical` data (`Nphi = 2`), interpolate a nonconstant `z`-axisymmetric Cartesian `g4DD` and `Gamma4UDD` across `theta in {0, pi}` or `r = 0` and verify interior convergence order.

Evidence (descriptive behavior): `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/azimuthal_symmetry_spatial_lagrange_interpolation.py::register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation`; corroboration: `output_raytracing_data.py`.
