# C And Embedded C Style

> Handwritten C/CUDA/H layout, `clang-format` control of embedded C strings, `// END ...` closing-brace comments, and Doxygen rules. · Status: provisional
> Up: [Architecture](index.md)

## Summary

Handwritten C/CUDA/H follows 2-space indentation, `snake_case` functions/structs (`_struct` suffix), `UPPER_SNAKE_CASE` macros/enums, `#ifndef` header guards, and Doxygen comments. For C/CUDA/H strings assembled in Python generators, `clang_format()`/`clang-format` controls whitespace and brace layout exclusively—raw generator-string cosmetics and generated-identifier naming styles are never reviewed or repaired—while semantic `// END <CONSTRUCT>: <desc>` comments (`<= 5` words except `nrpy/helpers/loop.py` footers), Doxygen accuracy, and compiler/runtime behavior remain enforced.

## Detail

### Handwritten C/CUDA/H Layout, Includes, And Header Guards

Handwritten C/CUDA/H uses 2 spaces (no tabs), `~100`-character lines, K&R same-line opening braces, no braces around single-statement control bodies, `//` inline comments (`//==========================` struct section dividers), declarations at first use with `const`/`restrict`, `static inline` header helpers with return type and name on one line, parenthesized macro parameters, and `#if defined(...)` / `#elif defined(...)` / `#endif` conditionals. Standard-library `#include` directives precede quoted project headers, with ISA headers behind `#if defined(__AVX512F__)`.

Handwritten headers use `#ifndef`/`#define` `UPPER_SNAKE_CASE` guards (e.g., `BHAHAHA_HEADER_H`; legacy `__SIMD_INTRINSICS_H__` also appears). By contrast, BHaH's generated `set_CodeParameters*.h` fragments are intentionally unguarded local declarations included inside multiple function scopes (`BHaH/CodeParameters.py::write_CodeParameters_h_files`).

Evidence (descriptive behavior): `coding_style.md` (`## C/H Coding Style`); `BHaH/CodeParameters.py::write_CodeParameters_h_files`; `c_function.py::CFunction.generate_full_function`.

### Embedded C Strings In Python And `// END ...` Comments

Use `r"""..."""` for embedded C (`rf"""...{var}..."""` with doubled `{{`/`}}` and Python-3.7-compatible local variables only when interpolating), appending conditional blocks in execution order via `body += ...`. Never review, enforce, or write Python helpers for raw embedded-C indentation, wrapping, or brace cosmetics, and never enforce naming styles on generated C identifiers.

Every non-trivial closing brace requires a `// END <KEYWORD>: <description>` comment (`FUNCTION`, `LOOP`, `WHILE`, `DO-WHILE` after `;`, `IF`, `ELSE IF`, `ELSE`, `SWITCH`, `OMP PARALLEL`, `OMP CRITICAL`, or anonymous-scope `BLOCK`) with a domain-specific description of **at most 5 words**:
- **Loop-footer exception**: `nrpy/helpers/loop.py` emits `} // END LOOP: for <idx> over [<low>, <high>)` with literal bound expressions exceeding 5 words across all backends; hand-assembled markers around `c_codegen` output still obey the 5-word cap.
- **C++ namespace guard**: C++ namespace closers (emitted in Dendro) wrap their `// END NAMESPACE: ...` comment inside `// clang-format off` / `// clang-format on` so `FixNamespaceComments` does not strip `END` and `:`.
- **Short handwritten block exception**: handwritten multi-statement blocks `< 5` lines with visible opening braces may omit `// END ...` (one-statement braced blocks remain disallowed).

Evidence (normative rule): `coding_style.md` (`#### Embedded C Code String Conventions`, `### 10. End-Curly-Brace Comments`); `nrpy/helpers/loop.py`.

### Doxygen Comments, `desc=`, And Common Generated C Idioms

Functions exceeding `~10` lines require Doxygen comments (placed on header declarations when split across `.h`/`.c`): open with `/**` and close with ` */` on separate lines, lead with the brief line (no `@brief`), separate brief/extended numbered steps from tags with blank ` *` lines, and order tags `@param` / `@param[in]` (`const *`) / `@param[out]` / `@param[in,out]` (1 line each, no dash after name) -> `@return` (omitted for `void`; enumerate error codes) -> blank ` *` -> `@note` / `@warning` / `@pre`. `register_CFunction(desc=...)` strings omit `/** ... */` delimiters (added by `CFunction`) and use `"""..."""` unless interpolation requires `rf"""..."""`.

Recurring C idioms across generated projects: `#ifndef REAL` / `#ifndef M_PI` / `#ifdef STANDALONE`, `#pragma GCC optimize("unroll-loops")` + `reset_options`, `MAYBE_UNUSED`, `BHAH_MALLOC`/`BHAH_FREE` with immediate null checks (preferring single large allocations, or `malloc`/`free` and stack VLAs where established), OpenMP `LOOP_OMP` / `parallel for` / `collapse(2)` / `critical`, and `#pragma omp simd` / `simd_intrinsics.h`.

Evidence (normative rule): `coding_style.md` (`## C/H Coding Style`, `#### Embedded C Code String Conventions`).

## Sources

- [coding_style.md](../../coding_style.md) - `## C/H Coding Style`, `#### Embedded C Code String Conventions`, `### 10. End-Curly-Brace Comments`
- [CodeParameters.py](../../nrpy/infrastructures/BHaH/CodeParameters.py) - `write_CodeParameters_h_files`
- [c_function.py](../../nrpy/c_function.py) - `CFunction.generate_full_function`
- [loop.py](../../nrpy/helpers/loop.py) - `loop1D` closing-brace comment format

## See Also

- Depends on: [Python Coding Style](python-coding-style.md)
- Depends on: [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
- See also: [Infrastructure Code Style](../infrastructures/infrastructure-code-style.md)
