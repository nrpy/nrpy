# Maintenance And Validation Helpers

> Utility helpers for formatting, trusted string files, caches, conditional writes, colored output, and runtime annotation checks. · Status: confirmed
> Up: [Helper APIs](index.md)

## Summary

These helpers support generated-text maintenance and developer feedback paths. `generic.py` provides ordered uniqueness, `clang-format` integration, trusted string validation, string diffs, and package-data copying. `cached_functions.py` stores hashed pickle caches under the user cache directory and provides cached SymPy simplification. `ConditionalFileUpdater` avoids unnecessary file writes by comparing generated content against existing content. `colorize_text.py` and `type_annotation_utilities.py` provide terminal color selection and lightweight runtime validation utilities.

## Detail

### Formatting, `validate_strings`, And Package File Copying

`superfast_uniq()` removes duplicate list elements while preserving first occurrence order. `copy_files()` reads named package-data files with `pkgutil.get_data()`, creates the requested project subdirectory, and writes those files there as binary output. Copied project files are generated or packaged output files, not separate compiled wiki sources unless they are deliberately registered as evidence.

`clang_format()` reads the `clang_format_options` NRPy parameter, constructs a `clang-format` command, and caches successful formatted output using the exact input string plus option string as the cache key. When the option string contains `file`, so that the style comes from a `.clang-format` file found relative to the working directory, the output depends on more than the cache key and the helper neither reads nor writes the cache. It raises `RuntimeError` when `clang-format` exits unsuccessfully and raises `OSError` with platform-aware installation guidance when the executable is not found. The helper therefore does not imply `clang-format` is always installed. It formats emitted source strings, though its `Popen.communicate()` call sets no timeout; see [Test Oracles And Safe Updates](../../validation/test-oracles-and-safe-updates.md) for subprocess timeout rules.

`validate_strings()` validates raw string output, usually generated C or similar text. It derives the caller's directory and function name, creates a sibling `tests/` directory when needed, and uses `tests/<caller>_<string_desc>.<file_ext>` as the trusted file path. If the trusted file exists, the function reads UTF-8 text with universal-newline translation and compares that string with `to_check` and raises a diff-producing `ValueError` on mismatch. If the file is missing, it writes the provided string as the new trusted file. That branch only captures candidate data; it does not compare, validate truth, compile, or run output. Full-output selection, normalization, independent review, and two-process acceptance belong to Test Oracles And Safe Updates.

Evidence (descriptive behavior): `generic.py::validate_strings`.

Trusted string files from `validate_strings()` are separate from trusted expression dictionaries. `validate_strings()` compares one raw string against one caller-derived text file and does not process symbolic expression dictionaries; the equation trusted-expression workflow is documented separately in [Expression Validation Helpers](../../validation/expression-validation-helpers.md).

`diff_strings()` uses `difflib.ndiff()` and reports only added or removed lines, omitting unchanged lines and intraline marker lines. This gives `validate_strings()` a compact mismatch report.

### `.nrpycache` Pickles, Conditional Writes, And Type Checks

`cached_functions.py` maps cache identifiers to SHA-256 names ending in `.nrpycache` below `appdirs.user_cache_dir("nrpy")`. `cache_file()` creates the directory if needed; `is_cached()`, `read_cached()`, and `write_cached()` manage pickle-backed entries. `cached_simplify()` hashes a pickled SymPy expression, returns a cached simplification when present, writes `sp.simplify()` results when absent, returns zero directly for the zero expression, and falls back to uncached `sp.simplify()` only for `pickle.PicklingError`. Other serialization exceptions, corrupt-cache reads, and filesystem errors propagate. These cache files are excluded from the KB. Tests must use isolated temporary cache directories without clearing or overwriting shared user caches; see [Test Oracles And Safe Updates](../../validation/test-oracles-and-safe-updates.md).

Evidence (descriptive behavior): `cached_functions.py::cached_simplify`.

`ConditionalFileUpdater` is a context manager that captures generated text in a `StringIO`, optionally runs `clang_format()`, reads the existing file if present, and compares stripped old and new content. It writes only when content changes and the module-level `nochange` flag is false. When module-level `verbose` is true, it prints a context diff before writing. Status text is colorized through `colorize_text`. `__exit__()` returns without comparing or writing when the context body raised an exception, so partial generated output never replaces the old file and a file that did not exist is not created; the exception then propagates unchanged.

Evidence (descriptive behavior): `conditional_file_updater.py::ConditionalFileUpdater.__exit__`.

`colorize_text.py` defines `ColorNames` as the allowed literal color names and maps them to ANSI escape sequences. `apply_colorization()` validates the color name and wraps stringified output; `leave_text_alone()` returns plain text. The module selects the exported `colorize` function at import time based on whether stdout looks like a terminal or a Jupyter output stream.

`type_annotation_utilities.py` provides `is_type_literal()` for recognizing `Literal` annotations across supported Python versions. `validate_literal_arguments()` inspects the calling function, checks recognized, evaluated `Literal[...]` annotations, and raises when a value is outside the allowed set. It locates the executing function by comparing code objects along the method resolution order of `self`, so a method inherited from a base class is checked against its own annotations, and it resolves quoted or postponed annotations with `typing.get_type_hints()`. If that resolution fails for the function, string annotations are left unresolved and skipped. `generate_class_representation()` builds a simple `ClassName(field=value, ...)` representation from public string and integer attributes on the calling object.

Evidence (descriptive behavior): `type_annotation_utilities.py::validate_literal_arguments`.

## Sources

- [nrpy/helpers/generic.py](../../../nrpy/helpers/generic.py) - `superfast_uniq`, `clang_format`, `validate_strings`, `diff_strings`, `copy_files`
- [nrpy/helpers/cached_functions.py](../../../nrpy/helpers/cached_functions.py) - `get_hash`, `cache_file`, `is_cached`, `read_cached`, `write_cached`, `cached_simplify`
- [nrpy/helpers/conditional_file_updater.py](../../../nrpy/helpers/conditional_file_updater.py) - `ConditionalFileUpdater`, `verbose`, `nochange`
- [nrpy/helpers/colorize_text.py](../../../nrpy/helpers/colorize_text.py) - `ColorNames`, `leave_text_alone`, `apply_colorization`, `colorize`
- [nrpy/helpers/type_annotation_utilities.py](../../../nrpy/helpers/type_annotation_utilities.py) - `is_type_literal`, `validate_literal_arguments`, `generate_class_representation`

## See Also

- Depends on: [Test Oracles And Safe Updates](../../validation/test-oracles-and-safe-updates.md)
- See also: [Symbolic Expression Utilities](symbolic-expression-utilities.md)
- See also: [Workflows](../../workflows.md)
