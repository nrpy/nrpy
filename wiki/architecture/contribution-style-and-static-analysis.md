# Contribution Style And Static Analysis

> Preserve the practical style, artifact, trusted-value, and static-analysis rules for code changes. · Status: confirmed · Last reconciled: 07-13-2026
> Up: [Architecture](index.md)

## Summary

The KB files contributor rules for new or modified code into focused leaves for Python style, C/embedded-C style, equation setup style, infrastructure code style, and static analysis. These rules are not claims that every legacy file already conforms. Handwritten Python changes require Black formatting in an isolated user-owned intended-change tree, plus single-file static analysis and Pylint **10.00/10.00** for every modified file. Documentation-only Markdown changes do not require Python static analysis unless source files are accidentally edited.

## Detail

For Python source changes, run `black .` only in an isolated, user-owned
intended-change worktree or copy with no unrelated modifications, then inspect
its diff. Run
`./.github/single_file_static_analysis.sh <path-to-file.py>` on each modified
handwritten Python file. The local script checks Black, isort, mypy in strict mode with
untyped calls allowed, Pylint, pydocstyle, and darglint, then executes the target
file. Each handwritten file must report Pylint **10.00/10.00**. Inspect direct-
execution effects before running the script; [Static
Analysis](../validation/static-analysis.md) owns exact mechanics and current
enforcement gaps.

The broader CI static-analysis job scans Python files while excluding generated
projects, build directories, trusted-value test directories, the obsolete
`./nrpy/examples/visualization_scripts/*` path, and selected special cases. That
exclusion does not match current `nrpy/examples/geodesic_visualizations/**`, so
current geodesic visualization companions remain in the scan. Configured CI
checks differ from the local script; Static Analysis owns those details.
Generated-project validation runs in separate jobs.

Trusted numerical value files are a special source class. They live under
`*/tests/`, contain generated reference values, and are regenerated from their
owning module rather than hand-edited. Generated trusted data is exempt from the
single-file script, but its handwritten owner is not. [Test Oracles And Safe
Updates](../validation/test-oracles-and-safe-updates.md) owns store admission
and update safety.

Style expectations include Black formatting, isort import grouping, canonical NRPy aliases such as `import nrpy.c_codegen as ccg` and `import nrpy.c_function as cfc`, Sphinx/reST docstrings, explicit return annotations, and the established doctest runner pattern for runnable equation modules. C and embedded C follow the project C/H style, including 2-space indentation, Doxygen-style function docs where required, and informative `// END ...` comments on non-trivial closing braces. Infrastructure generators add extra rules for `CodeParameter` registration scope, meaningful doctests, parallel codegen discovery guards, and BHaH symbolic-codegen helper use.

Artifact rules apply to documentation and code work: do not add binary files, images, archives, compiled outputs, generated projects, scratch logs, or other non-text assets unless maintainers approve. Generated trusted-value files are evidence, not prose pages, and generated runtime projects are outputs unless explicitly selected and registered as frozen evidence.

## Sources

- [coding_style.md](../../coding_style.md) - `## Python Coding Style`, `### Formatting`
- [raw/source-docs/original-agents.md](../../raw/source-docs/original-agents.md) - historical `## Required Checks`, plus `## Equation Setup Rules` and `## Quick Reference`; current `coding_style.md` decides conflicts
- [.github/single_file_static_analysis.sh](../../.github/single_file_static_analysis.sh) - `run_test_step`, static-analysis steps
- [.github/workflows/main.yml](../../.github/workflows/main.yml) - `static-analysis`, `codegen-ubuntu`, `codegen-mac`

## See Also

- Parent: [Architecture](index.md)
- See also: [Build And Run](build-and-run.md)
- See also: [Python Coding Style](python-coding-style.md)
- See also: [C And Embedded C Style](c-and-embedded-c-style.md)
- See also: [Equation Setup Style](../equations/equation-setup-style.md)
- See also: [Infrastructure Code Style](../infrastructures/infrastructure-code-style.md)
- Depends on: [Code Test Policy](../validation/code-test-policy.md)
- Depends on: [Test Oracles And Safe Updates](../validation/test-oracles-and-safe-updates.md)
- Validated by: [Static Analysis](../validation/static-analysis.md)
- Depends on: [Generated Output Boundaries](generated-output-boundaries.md)
- See also: [Symbolic Codegen Lifecycle](symbolic-codegen-lifecycle.md)
