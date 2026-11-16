# Equation Setup Style

> SymPy construction rules, `nrpy.indexedexp` loops, tensor/derivative suffix conventions, `trusted_dict` key stability, and prohibited dependencies in `nrpy/equations/`. · Status: provisional
> Up: [Equations](index.md)

## Summary

NRPy equation modules construct symbolic expressions with explicit nested loops over `nrpy.indexedexp` (`ixp`) arrays, exact `sp.Rational` fractions, `real=True` scalars, and standard tensor/derivative suffixes (`U`, `D`, `DD`, `dD`, `dDD`, `dupD`, `rhs`), avoiding `sp.simplify()`, general `.subs()`, `import re`, and `numpy`. Modules validate exact algebraic identities first and use `validate_expressions.py` `trusted_dict` checks when sampled high-precision regression adds signal.

## Detail

### SymPy Rules, `nrpy.indexedexp` Loops, And Suffix Conventions

Avoid `sp.simplify()` in equation builders (permitted in tests/identity checks; `nrpy.helpers.cached_functions.cached_simplify()` is reserved for small local subexpressions that materially improve GR symbolic tractability, never whole tensors or assembled RHSs). Avoid pattern-based `sp.subs()` / `sp.replace()`; permitted `.subs()` uses are elliptic coordinate substitution, GRHD face-value substitution, `nrpyABS -> sp.Abs`, horizon radius substitution, and evaluating waveform expressions at an attachment point. Initialize accumulators with `sp.sympify(0)`/`sp.sympify(1)`, fractions with `sp.Rational()`, and scalars with `sp.symbols(..., real=True)`.

Import `nrpy.indexedexp as ixp` (`zerorank1..4`, `declarerank1..4`, `2x2..4x4` symmetric inversions, 3D Levi-Civita helpers; symmetries `sym01`, `sym12`, `sym01_sym23`, `nosym`) and accumulate contractions with explicit Python `for` loops—never Einstein summation strings or matrix `@` operators.

Standard suffixes: `U` (contravariant), `D` (covariant), `DD`/`UU` (rank-2), `dD`/`dDD` (centered 1st/2nd partial derivatives), `dupD` (upwinded 1st derivative, explicitly permitted alongside `dD`/`dDD` by `coding_style.md` and `finite_difference.py::C_CODEGEN_DERIVATIVE_FAMILIES`), `dBarD`/`dHatD` (conformal/reference covariant derivatives), and `rhs` (time-evolution RHS); see [Finite Difference](../core/finite-difference.md).

Evidence (normative rule): `coding_style.md` (`## Equation Setup Patterns`, `### Derivative Variable Naming Conventions`); `finite_difference.py::C_CODEGEN_DERIVATIVE_FAMILIES`.

### `trusted_dict` Validation Keys And Prohibited Dependencies

Prefer exact identities first; when sampled regression is appropriate, map stable keys to SymPy expressions, call `ve.process_dictionary_of_expressions(..., fixed_mpfs_for_free_symbols=True)`, and compare via `ve.compare_or_generate_trusted_results(...)` ([Expression Validation Helpers](../validation/expression-validation-helpers.md)). Because list entries flatten to `<name>_<idx>`, `funcform` and unsupported types are skipped, and duplicate processed keys overwrite silently, choose distinct names matching `tests/<module>*.py` (`trusted_dict` files hold only `mpf`/`mpc` `# type: ignore` imports + `trusted_dict` and are updated only via the two-process procedure in [Test Oracles And Safe Updates](../validation/test-oracles-and-safe-updates.md); doctest-only modules such as `quaternion_rotations/tensor_rotation.py` need no artificial `trusted_dict`).

Core equation modules must not `import numpy` (reserved for post-processing/visualization scripts) or `import re` (grandfathered legacy sites in `coding_style.md` are removed only when substantively edited; repo tools outside `nrpy/` are exempt).

Evidence (normative rule): `coding_style.md` (`### Expression Validation via Trusted Dictionaries`, `### Prohibited Dependencies`); `validate_expressions.py::process_dictionary_of_expressions`.

## Sources

- [coding_style.md](../../coding_style.md) - `## Equation Setup Patterns`, `### Derivative Variable Naming Conventions`, `### Expression Validation via Trusted Dictionaries`, `### Prohibited Dependencies`
- [finite_difference.py](../../nrpy/finite_difference.py) - `C_CODEGEN_DERIVATIVE_FAMILIES`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `assert_equal`, `compare_or_generate_trusted_results`
- [tensor_rotation.py](../../nrpy/equations/quaternion_rotations/tensor_rotation.py) - `rotate`, doctest-only `__main__`

## See Also

- Depends on: [Indexed Expressions](../core/indexed-expressions.md)
- Depends on: [Finite Difference](../core/finite-difference.md)
- Depends on: [Test Oracles And Safe Updates](../validation/test-oracles-and-safe-updates.md)
