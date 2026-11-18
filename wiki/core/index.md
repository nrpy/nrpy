# Core APIs

> Router for NRPy's core symbolic-to-generated-code modules.

| Page | Go here when... |
| --- | --- |
| [C Codegen](c-codegen.md) | You need SymPy-to-C entry points, codegen options, or finite-difference-aware C code emission. |
| [C Function Registry](c-function-registry.md) | You need `CFunction`, `CFunction_dict`, or how infrastructure modules register generated C functions. |
| [Python Codegen](python-codegen.md) | You need JAX-compatible Python assignment emission from SymPy expressions, `PyCodeGen`, or `py_codegen`. |
| [Python Function Registry](python-function-registry.md) | You need generated Python function text assembly, `PyFunction`, `PyFunction_dict`, or `register_PyFunction`. |
| [Gridfunctions And Parameters](gridfunctions-and-parameters.md) | You need gridfunction registration, code parameters, or generated grid/common data structs. |
| [Indexed Expressions](indexed-expressions.md) | You need tensor declaration helpers, `create_tensor_symbolic`, `symmetry_axes`, symmetry strings, matrix inversion, or Levi-Civita helpers. |
| [Reference Metrics](reference-metrics.md) | You need coordinate-system reference metrics, supported `CoordSystem` names, `rfm_dict`, reference-metric validation, or precompute behavior. |
| [Finite Difference](finite-difference.md) | You need finite-difference coefficients, derivative extraction, gridfunction memory reads, or FD helper functions. |
| [Helper APIs](helpers/index.md) | You need shared helper modules for expressions, CSE/printers, SIMD/intrinsics, loops, device kernels, parallel codegen, or maintenance utilities. |

Back to [AGENTS.md](../../AGENTS.md).
