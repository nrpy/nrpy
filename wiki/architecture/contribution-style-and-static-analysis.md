# Contribution Style And Static Analysis

> Preserve the practical style, artifact, trusted-value, and static-analysis rules for code changes. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Architecture](index.md)

## Summary

The KB files project style rules into focused leaves for Python style, C/embedded-C style, equation setup style, infrastructure code style, and static analysis. Python changes require Black formatting before commit and single-file static analysis for every modified Python file, except generated trusted-value files under `*/tests/*.py`. Documentation-only Markdown changes do not require Python static analysis unless source files are accidentally edited.

## Detail

For Python source changes, run `black .` before commit. Then run `./.github/single_file_static_analysis.sh <path-to-file.py>` on each modified Python file. The local script checks Black, isort, mypy in strict mode with untyped calls allowed, pylint with a local threshold of 9.91, pydocstyle, darglint, and doctests by executing the target file.

The broader CI static-analysis job scans Python files while excluding generated projects, build directories, trusted-value test directories, visualization scripts, and selected special cases. CI runs across multiple Ubuntu and Python-version combinations, adjusts some tools by Python version, and also performs generated-project validation in separate jobs.

Trusted numerical value files are a special source class. They live under `*/tests/`, contain generated reference values, skip single-file static analysis, and are regenerated from their owning module rather than hand-edited. Legitimate trusted-output changes should be explained in the commit message.

Style expectations include Black formatting, isort import grouping, canonical NRPy aliases such as `import nrpy.c_codegen as ccg` and `import nrpy.c_function as cfc`, Sphinx/reST docstrings, explicit return annotations, and the established doctest runner pattern for runnable equation modules. C and embedded C follow the project C/H style, including 2-space indentation, Doxygen-style function docs where required, and informative `// END ...` comments on non-trivial closing braces. Infrastructure generators add extra rules for `CodeParameter` registration scope, meaningful doctests, parallel codegen discovery guards, and BHaH symbolic-codegen helper use.

Artifact rules apply to documentation and code work: do not add binary files, images, archives, compiled outputs, generated projects, scratch logs, or other non-text assets unless maintainers approve. Generated trusted-value files are evidence, not prose pages, and generated runtime projects are outputs unless explicitly selected and registered as frozen evidence.

## Sources

- [raw/source-docs/original-agents.md](../../raw/source-docs/original-agents.md) - `## Required Checks`, `## Equation Setup Rules`, `## Quick Reference`
- [.github/single_file_static_analysis.sh](../../.github/single_file_static_analysis.sh) - `run_test_step`, static-analysis steps
- [.github/workflows/main.yml](../../.github/workflows/main.yml) - `static-analysis`, `codegen-ubuntu`, `codegen-mac`

## See Also

- [Architecture](index.md)
- [Build And Run](build-and-run.md)
- [Python Coding Style](python-coding-style.md)
- [C And Embedded C Style](c-and-embedded-c-style.md)
- [Equation Setup Style](../equations/equation-setup-style.md)
- [Infrastructure Code Style](../infrastructures/infrastructure-code-style.md)
- [Generated Output Boundaries](generated-output-boundaries.md)
- [Symbolic Codegen Lifecycle](symbolic-codegen-lifecycle.md)
