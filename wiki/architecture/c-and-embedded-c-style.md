# C And Embedded C Style

> C/H formatting, Doxygen, embedded-C string, and generated-C body rules. · Status: provisional · Last reconciled: 07-20-2026
> Up: [Architecture](index.md)

## Summary

New or modified handwritten C/CUDA/H follows the project indentation, naming,
brace, header-guard, and documentation conventions. For C/CUDA/H strings stored
or assembled in Python generators, existing `clang_format`/`clang-format`
normalization exclusively owns final layout. Raw generator-string cosmetics
for indentation, tabs, spacing, alignment, wrapping, and brace layout must never
be reviewed, enforced, or repaired. Mandatory semantic C/C++ `//` line-comment
`END` marker presence, correct construct keyword, colon separator, and accurate
meaningful high-signal description of at most five words remain enforced for
generated code. No
generated-identifier naming style is reviewed or enforced.
Generated semantics, interfaces, documentation content and syntax,
compiler behavior, and runtime behavior remain in scope.

This leaf records contributor rules, not a claim that every handwritten or
generated legacy C fragment already conforms. Existing exceptions are evidence
about current source, not automatic precedent for new code.

## Detail

### Baseline Handwritten C/CUDA/H Style

These baseline layout and naming rules apply only to handwritten C/CUDA/H. Use
2 spaces and no tabs. Keep function, variable, and struct names in
snake_case; structs use a `_struct` suffix. Macros and enum values use
UPPER_SNAKE_CASE. Place opening braces on the same line for functions, control
structures, and struct definitions.

Header guards use traditional `#ifndef` and `#define` pairs with UPPER_SNAKE_CASE
names. Both simple names such as `BHAHAHA_HEADER_H` and legacy double-underscore
forms such as `__SIMD_INTRINSICS_H__` appear in the codebase. Put standard
library includes first, then project headers in quotes, with platform-specific
headers behind conditional compilation such as `#if defined(__AVX512F__)`.

Function-like macros parenthesize arguments and macro parameter uses.
Conditional compilation should use `#if defined(...)`, `#elif defined(...)`,
and `#endif`. Header-defined helper functions use `static inline`; keep return
type and function name on one line, and use `const` and `restrict` qualifiers
where they describe actual use.

Declare variables at first use, prefer `const` for immutable values, use
`restrict` for performance-sensitive pointers, and group declarations when that
improves scanability. C line length is about 100 characters, while Python line
wrapping is left to Black.

### Embedded C In Python

Use `r"""..."""` for embedded C bodies so backslashes are literal. Use
`rf"""...{var}..."""` only when interpolation is needed, double literal C
braces as `{{` and `}}`, and keep interpolation Python-3.7-compatible by moving
complex expressions to preceding locals.

For C, CUDA, and header text stored or assembled in Python generators, raw
tabulation or tabs, indentation,
dedentation or reindentation, spacing or alignment, wrapping, brace placement or
layout cosmetics must never be reviewed, enforced, or repaired. Final generated
layout belongs exclusively to the
existing `clang_format`/`clang-format` normalization. If an emission path lacks
final normalization, do not compensate with raw-string cosmetics or helpers.
This exemption does not cover `// END ...` closing-brace comments: handwritten
C/CUDA/H must follow the full rules below, and Python-generated C/CUDA/H must
follow the semantic marker subset below.

When Python options conditionally include C body sections, build the body in
execution order with `body += ...`. Prefer top-to-bottom construction that
matches generated C execution. Do not add a Python function, method, lambda,
wrapper, utility, transform, postprocessor, or other special routine whose
purpose is to format generated C/CUDA/H strings cosmetically by hand.

Do not review or enforce any casing, prefix, suffix, `snake_case`, or other
naming style for generated identifiers. Exact current API, registry, prototype,
compiler, collision-avoidance, and semantic requirements remain enforceable.
This exemption is formatting-only: raw/rf string and doubled-brace validity,
interpolation correctness, semantic execution order, Doxygen and interface
content accuracy, syntax, compiler behavior, and runtime behavior remain in
scope. Generated `#include` behavior remains an interface/compiler concern, not
a layout exception.

Claim evidence:
- Claim: Existing clang-format normalization exclusively owns indentation, tabs, spacing, alignment, wrapping, and brace layout for C/CUDA/H strings stored or assembled in Python generators; those raw cosmetics must never be reviewed, enforced, or repaired, Python-only cosmetic formatting routines are prohibited, and no generated-identifier naming style is reviewed or enforced, while mandatory generated semantic C/C++ `//` line-comment `END` marker presence, correct construct keyword, colon separator, accurate meaningful high-signal description of at most five words, and exact API, registry, prototype, semantic, interface, documentation, syntax, compiler, collision, and runtime requirements remain enforceable.
- Role: normative rule
- Deciding authority: `coding_style.md` - `#### Embedded C Code String Conventions`, `## C/H Coding Style`
- Corroboration: none available; the frozen historical style source conflicts where it applies handwritten layout rules to generated strings.

### End-Curly-Brace Comments

For Python-generated C/CUDA/H, review and enforce only this semantic marker
subset for every closing brace that ends a non-trivial block: an `END` marker
must be carried by a C/C++ `//` line comment, use the correct
syntactic-construct keyword, include a colon separator, and carry an accurate,
meaningful, high-signal description of at most five words. Do not review or enforce its exact
whitespace, alignment, wrapping, placement, or brace shape.

For new handwritten C/CUDA/H, every closing brace that ends a non-trivial block
must carry an informative `// END ...` comment. The keyword is all caps,
followed by a colon and a concise description. Use the syntactic construct as
the keyword, keep the description to at most five words, and preserve its
highest-signal semantic context:

- `} // END FUNCTION: compute_spin`
- `} // END LOOP: for h over all horizons`
- `} // END WHILE: refining spin until convergence`
- `} // END IF: destination-point allocation failed`
- `} // END ELSE IF: num_resolutions_multigrid > 0`
- `} // END ELSE: not enable_BBH_mode`
- `} // END SWITCH: select inner-boundary gridfunctions`
- `} // END OMP PARALLEL: reduce per-thread diagnostics`
- `} // END OMP CRITICAL: update shared centroid diagnostics`
- `} // END BLOCK: Cartesian-to-grid conversion`

For handwritten `do...while`, put the comment after the semicolon:
`} while (condition); // END DO-WHILE: brief description`.

`END BLOCK` is for anonymous scoping blocks only. Do not use generic labels
when the real purpose can be stated. Prefer domain descriptions such as `theta
points on the horizon surface` over low-signal text such as `grid index`.
For handwritten C/CUDA/H only, do not add braces around one-statement
control-flow bodies in new code. Write the single statement directly under the
control-flow header instead.

### Doxygen And `desc=` Text

Documentation-content requirements apply to handwritten and generated
functions. Every C function whose body exceeds roughly 10 lines requires
accurate Doxygen documentation; shorter functions may omit it if the name and
signature are self-documenting. The raw placement, blank-line, alignment, and
wrapping rules below apply only to handwritten C/CUDA/H. If a handwritten
function is declared in a header and defined in a `.c` file, document the
declaration; otherwise document the definition.

Doxygen comments start with `/**` on its own line and close with ` */` on its
own line. The first content line is the brief description; do not use `@brief`.
Insert a blank ` *` line after the brief, and after any extended description,
before tags. For complex functions, prefer a numbered step list over dense
prose; avoid extended descriptions near the 10-line threshold unless the logic
is genuinely non-obvious.

Tag order is `@param`, then `@return`, then a blank ` *` line before
`@note`, `@warning`, or `@pre`. Keep `@param` descriptions to one line and do
not put a dash after the parameter name. `const` pointer parameters are
`@param[in]`; non-`const` pointer parameters use `@param[out]` or
`@param[in,out]`; pass-by-value parameters use plain `@param`. Omit `@return`
for `void` functions. For integer error-code returns, enumerate outcomes or
name the enum values. Use standalone `@warning` tags for correctness hazards,
not inline warnings inside parameter descriptions.

`desc=` strings passed to `register_CFunction()` must provide accurate interface
content and valid Doxygen tag syntax, but their raw blank lines, indentation,
alignment, and wrapping are not enforced. The framework emits the `/**` and
` */` delimiters, so the string itself must not include them. Multiline `desc`
strings use triple double quotes whenever possible, with raw f-triple-quoted
strings only when interpolation is required.

### Common C Patterns

In handwritten C/CUDA/H, struct fields are grouped by purpose with visible
`//==========================` section separators, and inline or block comments
inside functions use `//` rather than `/* */`, except for Doxygen comments.

Common generated-code patterns include:

- `#ifndef REAL` / `#define REAL double` type flexibility.
- `#ifndef M_PI` portability guards.
- `#ifdef STANDALONE` self-contained test programs.
- `#pragma GCC optimize("unroll-loops")` before performance-critical functions
  and `#pragma GCC reset_options` after them.
- `MAYBE_UNUSED` for variables unused on some code paths.
- `BHAH_MALLOC(ptr, size)` and `BHAH_FREE(ptr)` for tracked allocation, with
  immediate null checks and error returns.
- Standard `malloc()` and `free()` for simple cases.
- Single large heap allocations preferred over many small allocations.
- VLAs for stack-allocated temporary buffers where used by existing C.
- OpenMP patterns such as `#pragma omp parallel for`, reductions,
  `collapse(2)`, `critical`, and generated `LOOP_OMP(...)` loops.
- SIMD patterns such as `sum_lagrange_x0_simd()`, `#pragma omp simd`, and
  `simd_intrinsics.h` copied to projects when needed.

## Sources

- [coding_style.md](../../coding_style.md) - `#### Embedded C Code String Conventions`, `## C/H Coding Style`, `## Style Comparison Summary`; current authority for formatter ownership and generated naming scope
- [original-agents.md](../../raw/source-docs/original-agents.md) - historical `## C/H Style`; current `coding_style.md` decides conflicts
- [original-agents.md](../../raw/source-docs/original-agents.md) - historical `### Embedded C in Python Strings`, `### C Function Registration from Python`; current `coding_style.md` decides conflicts
- [original-agents.md](../../raw/source-docs/original-agents.md) - historical `### Preprocessor / Comment Patterns`, `## Additional Project Rules`, `## Quick Reference`

## See Also

- Parent: [Architecture](index.md)
- Depends on: [Python Coding Style](python-coding-style.md)
- Depends on: [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
- See also: [Infrastructure Code Style](../infrastructures/infrastructure-code-style.md)
