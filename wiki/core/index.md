# Core APIs

> Router for NRPy's core symbolic-to-generated-C API surfaces.

| Page | Go here when... |
| --- | --- |
| [C Codegen](c-codegen.md) | You need the SymPy-to-C entry points, codegen options, or finite-difference-aware code emission. |
| [C Function Registry](c-function-registry.md) | You need to understand `CFunction`, `CFunction_dict`, or how infrastructure modules register generated C functions. |
| [Python Codegen](python-codegen.md) | You need JAX-compatible Python assignment emission from SymPy expressions, `PyCodeGen`, or `py_codegen`. |
| [Python Function Registry](python-function-registry.md) | You need generated Python function text assembly, `PyFunction`, `PyFunction_dict`, or `register_PyFunction`. |
| [Gridfunctions And Parameters](gridfunctions-and-parameters.md) | You need gridfunction registration, code parameters, or generated grid/common data structs. |
| [Indexed Expressions](indexed-expressions.md) | You need tensor declaration helpers, `create_tensor_symbolic`, `indexedexp::symmetry_axes`, symmetry strings, zero tensors, matrix inversion, or Levi-Civita helpers. |
| [Reference Metrics](reference-metrics.md) | You need coordinate-system reference metrics, supported CoordSystems, the `rfm_dict` cache, module-level reference-metric validation, or reference-metric precompute behavior. |
| [Finite Difference](finite-difference.md) | You need finite-difference coefficients, derivative extraction, gridfunction memory reads, or FD helper functions. |
| [Helper APIs](helpers/index.md) | You need shared helper modules for expressions, CSE/printers, SIMD/intrinsics, loops, device utilities, parallel codegen, or maintenance support. |

Back to [AGENTS.md](../../AGENTS.md).
