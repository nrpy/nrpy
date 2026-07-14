# Equation Setup Style

> Symbolic equation construction, naming, validation, and dependency rules for NRPy equation modules. · Status: provisional · Last reconciled: 07-13-2026
> Up: [Equations](index.md)

## Summary

Equation modules build symbolic expressions explicitly with SymPy and
`nrpy.indexedexp`, avoid broad simplification and pattern substitution, use
standard tensor and derivative suffixes, and validate outputs through trusted
numerical dictionaries only when sampled regression is meaningful. The style
contract favors explicit loops, exact fractions, stable result-key names, and
exact or semantic invariants before sampled golden data.

## Detail

### SymPy Usage

Avoid `sp.simplify()` in equation-building code. It is acceptable in
test/validation code or explicit identity checks, and some modules document
"No SymPy .subs() or simplify() calls" as a local enforcement rule.

`nrpy.helpers.cached_functions.cached_simplify()` is a narrow exception for
small local subexpressions when simplification materially improves downstream
symbolic tractability, code-generation stability, or expression readability,
especially in general-relativity equations. Do not use it as blanket
permission to simplify whole tensors or large assembled RHS expressions. Add a
short comment when the need is not obvious.

Do not use `sp.subs()` or `sp.replace()` for pattern-based expression
transformation in core equation building. Allowed `.subs()` cases are
coordinate substitution in elliptic source terms, face-value substitution in
GRHD fluxes, systematic `nrpyABS` to `sp.Abs` conversion, surface-radius
substitution in horizon modules, and evaluating an expression at a specific
parameter value, such as waveform attachment-point calculations.

Use `sp.sympify(0)` or `sp.sympify(1)` for accumulator initialization,
`sp.Rational()` for exact fractions, and `sp.symbols(..., real=True)` for
scalar quantities.

### Indexed Expressions And Explicit Accumulation

Use `nrpy.indexedexp` imported as `ixp` for tensor-like symbolic arrays.
Standard helpers include `ixp.zerorank1()` through `ixp.zerorank4()`,
`ixp.declarerank1()` through `ixp.declarerank4()`, symmetric matrix inversions
for 2x2 through 4x4 matrices, and 3D Levi-Civita symbol/tensor helpers.

Use supported symmetry strings such as `sym01`, `sym12`, `sym01_sym23`, and
`nosym`. Build tensor expressions with explicit nested loops and component
accumulation. Do not use Einstein summation notation or matrix multiplication
operators for equation construction.

### Naming Conventions

Tensor and derivative suffixes encode meaning:

- `U` means contravariant, as in `betaU`, `sU`, and `LambdabarU`.
- `D` means covariant or first derivative context, as in `alpha_dD`, `sD`, and
  `SD`.
- `DD` and `UU` mean rank-2 covariant or contravariant tensors, such as
  `gammaDD` and `AbarUU`.
- `dD` and `dDD` mean first and second partial derivatives, such as `cf_dD`
  and `cf_dDD`.
- `dupD` means upwinded derivative, as in `trK_dupD`.
- `dBarD` and `dHatD` mean conformal or reference-metric covariant derivative.
- `rhs` means evolution-equation right-hand side, as in `alpha_rhs` and
  `gammabar_rhsDD`.

Derivative-like symbols intended to represent finite-difference-like
derivatives must include the `dD` or `dDD` components. Do not invent derivative
names that omit those markers. Use `nrpy.indexedexp` declaration helpers so
SymPy symbols carry expected suffix and component encoding.

### Trusted Expression Validation

Prefer an exact analytic, symbolic, or semantic invariant. When sampled
trusted-result regression is appropriate, equation modules build a dictionary
that maps stable descriptive names to SymPy expressions, process it with
`ve.process_dictionary_of_expressions(...)` using fixed substitutions for free
symbols, and pass the processed results to
`ve.compare_or_generate_trusted_results(...)`. A direct dictionary
`ve.assert_equal(...)` call rejects non-identical raw key sets, different list
nesting, and flattened-name collisions among entries retained for numerical
processing; `funcform` keys are omitted from numerical processing and collision
checks. Positional dictionary comparison is unsupported.

Claim evidence:
- Claim: Direct dictionary `ve.assert_equal(...)` calls reject non-identical raw key sets, different list nesting, and flattened-name collisions among numerically processed entries; `funcform` keys are omitted from numerical processing and collision checks, and positional dictionary comparison is unsupported.
- Role: descriptive behavior
- Deciding authority: [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py), `assert_equal`
- Corroboration: [test_parse_BSSN.py](../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py), `test_example_BSSN`, exercises direct dictionary comparison across scalar, vector, and matrix expression values; error-path tests remain colocated with the deciding helper
- Validation: `inspected=pass; generated=not-run; built=not-run; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=not-applicable; precision=30 decimal digits; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=pass; options=68 validator doctests plus explicit dictionary-structure and funcform-collision probes; date=07-13-2026`

Expression keys must match the corresponding `trusted_dict` keys in
`tests/<module>*.py`. Use established names such as `*_rhs`, `*_expr_list_*`,
or the owning module family's existing convention. Do not invent new result-key
conventions without aligning the trusted files.

Trusted-value files under `*/tests/` are generated artifacts. They contain only
the needed `mpf` or `mpc` import with `# type: ignore` and a `trusted_dict`
mapping keys to high-precision values. They have no module docstrings,
functions, or classes. Do not hand-edit trusted values. For legitimate
algorithm changes, follow the isolated, independently reviewed two-process
procedure in [Test Oracles And Safe
Updates](../validation/test-oracles-and-safe-updates.md). That page owns store
admission, format, and update safety.

Do not add golden data merely to satisfy a test count. A meaningful doctest-only
equation module, including quaternion tensor rotation, is valid when no sampled
oracle is needed. [Code Test Policy](../validation/code-test-policy.md) owns
owner-runner placement and the meaningful-contract gate.

Key validation APIs are `assert_equal(vardict_1, vardict_2)`, `check_zero`,
`process_dictionary_of_expressions`, `compare_against_trusted`,
`output_trusted`, and `compare_or_generate_trusted_results`.

### Dependencies

Do not import `re` for simple string manipulation that direct string methods
can handle. Regex is allowed only for genuine pattern matching, such as
coordinate-system variant detection or variable text layouts that cannot be
handled robustly with plain replacement. If a regex is needed, follow existing
legitimate examples and add a comment explaining why `.replace()` is
insufficient.

Core NRPy code must not depend on numpy. The core workflow is symbolic
expression to generated C, with computation done through SymPy and code
generation. Visualization and post-processing scripts may use numpy for image
handling, plotting, binary parsing, or similar non-core analysis.

## Sources

- [original-agents.md](../../raw/source-docs/original-agents.md) - `## Equation Setup Rules`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### SymPy`, `### Indexed Expressions`, `### Expression Construction`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### Expression Validation`, `### Prohibited / Restricted Dependencies`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `assert_equal`, `compare_or_generate_trusted_results`
- [tensor_rotation.py](../../nrpy/equations/quaternion_rotations/tensor_rotation.py) - `rotate`, module `__main__` path

## See Also

- Parent: [Equations](index.md)
- Depends on: [Indexed Expressions](../core/indexed-expressions.md)
- Validated by: [Trusted Expression Pipeline](trusted-expression-pipeline.md)
- Validated by: [Expression Validation Helpers](../validation/expression-validation-helpers.md)
- Depends on: [Test Oracles And Safe Updates](../validation/test-oracles-and-safe-updates.md)
- Depends on: [Code Test Policy](../validation/code-test-policy.md)
- See also: [Contribution Style And Static Analysis](../architecture/contribution-style-and-static-analysis.md)
