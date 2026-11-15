# Language For Computational Physicists

> Translate vague software-management wording into the direct language used by computational physicists. · Status: confirmed
> Up: [NRPy Knowledge Base](../AGENTS.md)

## Summary

Write prose by naming the equation, quantity, operation, numerical method, file, or result under discussion. This rule governs every explanation read by people or agents: KB pages, documentation, comments, docstrings, CLI help, diagnostics, and text emitted by generators. Follow the NRPy tutorial style: state the mathematical problem, identify inputs and outputs, specify numerical assumptions, and name the generated or evaluated code directly.

## Detail

### Start From The Computation

Before writing a sentence, identify its concrete subject in NRPy's technical vocabulary: symbolic expression, scalar, tensor, indexed expression, gridfunction, parameter, coordinate system, finite-difference derivative, stencil, right-hand side, initial data, boundary condition, generated C code, source file, executable, numerical grid, exact solution, convergence test, and validation result.

State directly what is differentiated, contracted, substituted, registered, generated, compiled, evolved, compared, or written, giving tensor indices, coordinate basis, units, finite-difference order, inputs, outputs, and failure conditions where relevant.

### Translation Table

Choose the specific physical, numerical, or file-level noun that matches the computation:

| Avoid in explanatory prose | Prefer |
| --- | --- |
| artifact | generated file, output file, executable, plot, data table, test result, build output, or named source file |
| contract | required behavior, interface definition, mathematical assumption, numerical requirement, or expected result |
| schema | file format, data layout, field definitions, parameter list, table columns, or page layout |
| surface | public functions, command-line options, configuration parameters, or host interface |
| payload | values, input data, output data, message fields, or array contents |
| plumbing | call sequence, data path, registration path, build steps, or input/output path |
| pipeline | ordered calculation stages, sequence of kernels, or named numerical method |
| bridge | host-transfer array, conversion function, coupling term, or named interface |
| bundle | photon chunk, component array, group of fields, or another named collection |
| gate | required check, acceptance condition, comparison, or test |
| harness | test program, test driver, comparison script, or generated example |
| fixture | test input, reference data, initial data, or expected output |
| consumer | calling function, generated solver, host code, reader, or downstream calculation |
| owner | defining module, generating function, equation module, or responsible test |
| lifecycle | generation, initialization, evolution, output, and cleanup steps, naming only the steps that apply |

Do not replace one vague noun with another. Rewrite “the generation artifact satisfies the schema contract” as a checkable statement such as “the generated TOML file contains every runtime parameter read by `rhs_eval_block`.”

### Exact And Mathematical Names

Preserve exact code identifiers and established mathematical or physical terminology: “contract the upper index with the first lower index,” apparent-horizon or ADM two-surfaces, normal or tangent bundles, named algorithms such as Split-Pipeline, and explicit memory, POSIX file, or distributed subdomain ownership when the sentence names what is owned. Exact filenames (`schema.json`), GitHub Action names (`actions/upload-artifact`), public symbol names, and direct quotations stay unchanged inside code spans or quotation marks. Outside exact names and quotations, reserve “artifact” for its archaeological meaning.

### Required Review And Conciseness

For every new or revised passage:

1. Replace left-column software metaphors (including inflections such as `owning`, `owned`, `owns`, `ownership`) with the concrete mathematical, numerical, or file-level subject.
2. Eliminate prose paragraphs that restate adjacent Markdown table cells; put per-row evidence directly in table cells with compact ``(Ev: `module.py::symbol`)`` tags when needed.
3. Omit repeated generic symbolic-versus-runtime disclaimers already covered by [KB Rules](SCHEMA.md#volatile-information-policy); state only module-specific untested paths or failure modes.
4. Compare recurring terminology with [Glossary](glossary.md) and review meaning rather than performing blind global word replacements.

Evidence (normative rule): `AGENTS.md::Human-Facing Language`; corroboration: NRPy+ scalar-wave and 10-minute overview tutorials.

## Sources

- [AGENTS.md](../AGENTS.md) - `## Human-Facing Language`
- [NRPy+ 10-Minute Overview](https://github.com/zachetienne/nrpytutorial/blob/master/Tutorial-NRPyPlus_10_Minute_Overview.ipynb) - problem statements, symbolic expressions, tensors, gridfunctions, finite differences, and generated C code
- [Start-to-Finish Scalar Wave](https://github.com/zachetienne/nrpytutorial/blob/master/Tutorial-Start_to_Finish-ScalarWave.ipynb) - equations, initial data, numerical methods, exact-solution comparison, and convergence validation

## See Also

- Depends on: [Glossary](glossary.md)
- See also: [KB Rules](SCHEMA.md)
- Validated by: [Workflows](workflows.md#deterministic-and-review-checks)
