# C And Embedded C Style

> C/H formatting, Doxygen, embedded-C string, and generated-C body rules. · Status: confirmed · Last reconciled: 07-05-2026
> Up: [Architecture](index.md)

## Summary

C and embedded C use 2-space indentation, K&R braces, snake_case functions and
variables, UPPER_SNAKE_CASE macros and enum values, traditional header guards,
and Doxygen comments for non-trivial functions. Generated C emitted from Python
registration functions follows the same conventions, with raw triple-quoted
strings for embedded bodies and informative `// END ...` comments on closing
braces.

## Detail

### Baseline C/H Style

Use 2 spaces and no tabs. Keep function, variable, and struct names in
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

Inside raw strings, C code should use 2-space indentation when practical, but
do not churn valid historical indentation for its own sake. `#include`
directives inside generated bodies use `#include "set_CodeParameters.h"` only
when the generator intentionally uses `include_CodeParameters_h=True`.

When Python options conditionally include C body sections, build the body in
execution order with `body += ...`. Prefer top-to-bottom construction that
matches generated C execution. Avoid temporary fragments or helper functions
whose only purpose is indentation, dedentation, reindentation, or cosmetic
formatting that `clang_format` will normalize later.

### End-Curly-Brace Comments

Every closing brace that ends a non-trivial block in new code should carry an
informative `// END ...` comment. The keyword is all caps, followed by a colon
and a concise description. Use the syntactic construct as the keyword and keep
high-signal semantic context in the description:

- `} // END FUNCTION: compute_spin`
- `} // END LOOP: for h over all horizons`
- `} // END WHILE: refining spin until convergence`
- `} // END IF: destination-point allocation failed`
- `} // END ELSE IF: num_resolutions_multigrid > 0`
- `} // END ELSE: not enable_BBH_mode`
- `} // END SWITCH: select gridfunctions requiring inner boundary conditions`
- `} // END OMP PARALLEL: reduce per-thread diagnostics`
- `} // END OMP CRITICAL: update shared centroid diagnostics`
- `} // END BLOCK: theta/phi stencil bounds and center-index sanity checks`

For `do...while`, put the comment after the semicolon:
`} while (condition); // END DO-WHILE: brief description`.

`END BLOCK` is for anonymous scoping blocks only. Do not use generic labels
when the real purpose can be stated. Prefer domain descriptions such as `theta
points on the horizon surface` over low-signal text such as `grid index`.
Omit end comments only when a multi-statement block body is fewer than five
lines and the opening brace remains visible without scrolling.

Do not add braces around one-statement control-flow bodies in new code. Write
the single statement directly under the control-flow header instead.

### Doxygen And `desc=` Text

Every C function whose body exceeds roughly 10 lines requires a Doxygen comment
immediately above it. Shorter functions may omit the comment if the name and
signature are self-documenting. If a function is declared in a header and
defined in a `.c` file, document the declaration; otherwise document the
definition.

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

`desc=` strings passed to `register_CFunction()` follow the same tag rules, but
the framework emits the `/**` and ` */` delimiters, so the string itself must
not include them. Multiline `desc` strings use triple double quotes whenever
possible, with raw f-triple-quoted strings only when interpolation is required.

### Common C Patterns

Struct fields are grouped by purpose with visible `//==========================`
section separators. Use `//` for inline and block comments inside functions,
not `/* */`, except for Doxygen comments.

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

- [original-agents.md](../../raw/source-docs/original-agents.md) - `## C/H Style`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### Embedded C in Python Strings`, `### C Function Registration from Python`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `### Preprocessor / Comment Patterns`, `## Additional Project Rules`, `## Quick Reference`

## See Also

- Parent: [Architecture](index.md)
- Depends on: [Python Coding Style](python-coding-style.md)
- See also: [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
- See also: [Infrastructure Code Style](../infrastructures/infrastructure-code-style.md)
