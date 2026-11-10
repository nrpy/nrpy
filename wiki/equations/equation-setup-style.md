# Equation Setup Style

> Symbolic equation construction, naming, validation, and dependency rules for NRPy equation modules. · Status: provisional
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

Use `dD` and `dDD` for centered first and second partial derivatives, and
`dupD` for upwinded derivatives as explicitly allowed by `coding_style.md`.
The guide's specific upwind convention qualifies its broader `dD`/`dDD`
wording; valid `dupD` names must not be renamed to satisfy that broader
sentence. Consult [Finite Difference](../core/finite-difference.md) for the
recognized derivative families and which can generate C code. Use
`nrpy.indexedexp` declaration helpers for the suffix and component encoding.

Claim evidence:
- Claim: The explicit `dupD` upwind naming rule qualifies the style guide's blanket `dD`/`dDD` sentence.
- Role: normative rule
- Deciding authority: [coding_style.md](../../coding_style.md), `### Derivative Variable Naming Conventions`, the explicit upwind bullet
- Corroboration: [finite_difference.py](../../nrpy/finite_difference.py), `C_CODEGEN_DERIVATIVE_FAMILIES`, recognizes `_dupD` for generated derivatives.

### Trusted Expression Validation

Prefer an exact analytic, symbolic, or semantic invariant. When sampled
trusted-result regression is appropriate, equation modules build a dictionary
that maps stable descriptive names to SymPy expressions, process it with
`ve.process_dictionary_of_expressions(...)` using fixed substitutions for free
symbols, and pass the processed results to
`ve.compare_or_generate_trusted_results(...)`. Direct dictionary comparisons
must also use stable descriptive keys. [Expression Validation
Helpers](../validation/expression-validation-helpers.md) owns helper behavior
and limitations.

Processed result keys must match the corresponding `trusted_dict` keys in
`tests/<module>*.py`. List entries become `<name>_<flat_index>`; `funcform`
keys and unsupported values are skipped, and colliding processed names silently
overwrite earlier entries. Choose names so every intended numerical output
survives processing with a distinct key. Use established names such as `*_rhs`,
`*_expr_list_*`, or the owning family's convention, and align any intentional
name changes with the trusted files.

Claim evidence:
- Claim: Trusted dictionaries compare processed names; list flattening, skipped values, and collisions can remove intended input outputs before comparison.
- Role: descriptive behavior
- Deciding authority: [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py), `process_dictionary_of_expressions`
- Corroboration: `none available`; the named implementation decides the stated behavior.

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

### Dependencies

Do not add `import re` anywhere under `nrpy/`; use plain string methods instead.
The legacy uses identified in `coding_style.md` are grandfathered rather than
precedent. Remove their regular-expression use when each file is next materially
edited; do not attach that removal to unrelated changes elsewhere. Repository
tooling outside `nrpy/`, such as `tools/kb_lint.py`, is outside this restriction.

Claim evidence:
- Claim: New code under `nrpy/` does not import `re`; the legacy sites named in `coding_style.md` are grandfathered until their next material edit and do not authorize new uses. Repository tooling outside `nrpy/` is excluded.
- Role: normative rule
- Deciding authority: [coding_style.md](../../coding_style.md), `### Prohibited Dependencies`, the `import re` bullet
- Corroboration: none available; the current style guide supersedes the frozen historical dependency guidance.

Core NRPy code must not depend on numpy. The core workflow is symbolic
expression to generated C, with computation done through SymPy and code
generation. Visualization and post-processing scripts may use numpy for image
handling, plotting, binary parsing, or similar non-core analysis.

## Sources

- [coding_style.md](../../coding_style.md) - `### Prohibited Dependencies`, `### Derivative Variable Naming Conventions`
- [finite_difference.py](../../nrpy/finite_difference.py) - `C_CODEGEN_DERIVATIVE_FAMILIES`

- [original-agents.md](../../raw/source-docs/original-agents.md) - `## Equation Setup Rules`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### SymPy`, `### Indexed Expressions`, `### Expression Construction`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### Expression Validation`, `### Prohibited / Restricted Dependencies`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `assert_equal`, `compare_or_generate_trusted_results`
- [tensor_rotation.py](../../nrpy/equations/quaternion_rotations/tensor_rotation.py) - `rotate`, module `__main__` path

## See Also

- Parent: [Equations](index.md)
- Depends on: [Indexed Expressions](../core/indexed-expressions.md)
- Validated by: [Trusted Expression Pipeline](trusted-expression-pipeline.md)
- Validated by: [Expression Validation Helpers](../validation/expression-validation-helpers.md)
- Depends on: [Test Oracles And Safe Updates](../validation/test-oracles-and-safe-updates.md)
- Depends on: [Code Test Policy](../validation/code-test-policy.md)
- See also: [Contribution Style And Static Analysis](../architecture/contribution-style-and-static-analysis.md)
