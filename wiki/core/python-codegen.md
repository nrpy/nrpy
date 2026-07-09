# Python Codegen

> Core route for turning SymPy expressions into JAX-compatible Python assignment text. · Status: confirmed · Last reconciled: 06-30-2026
> Up: [Core APIs](index.md)

## Summary

`py_codegen()` emits JAX-compatible Python assignment strings from SymPy expressions. `PyCodeGen` is the option object behind that public function: a stripped-down Python counterpart to `CCodeGen` for JAX-oriented code, without the finite-difference memory-read path used by C kernels.

## Detail

`PyCodeGen` accepts `prestring`, `poststring`, `verbose`, `enable_cse`, `cse_sorting`, `cse_varprefix`, and `postproc_substitution_dict`. The final returned string is the optional verbose comment block followed by `prestring`, generated assignment text, and `poststring`.

The constructor reads `par.parval_from_str("Infrastructure")` and requires the value to be exactly `"JAX"`. If the value differs, it raises `ValueError("Infrastructure must be 'jax' for py_codegen")`. This guard matches the implementation boundary: the module constructs a module-level `NRPyJaxPrinter`, and this Python codegen path targets JAX-compatible assignment output rather than general infrastructure-specific generated packages.

`py_codegen()` rejects tuple inputs for both the SymPy expression argument and the output-name argument. Non-list scalar inputs are normalized into one-element lists, list inputs are copied so later processing does not mutate caller-owned lists, and the expression list length must match the output-name list length.

When `verbose` is enabled, `py_codegen()` emits a Python comment block before generated code. The block records the original SymPy expression or expressions paired with their output assignment names; plural output uses bracketed `"[name = expression]"` lines.

With `enable_cse=False`, each expression is emitted independently. If `postproc_substitution_dict` is nonempty, `apply_substitution_dict()` first rewrites matching free-symbol names by appending the configured suffix. The selected expression is then passed to `printer.doprint(expr, output_name)`, and the printed assignment is appended to the output string.

With `enable_cse=True`, `py_codegen()` collects the expressions and output names, then calls SymPy CSE with numbered temporaries from `cse_varprefix + "tmp"` and the configured `cse_sorting` order. For SymPy versions before 1.3, the implementation prints a warning and uses the raw `sp.cse()` result. For SymPy 1.3 and newer, it passes the `sp.cse()` result through `cse_postprocess()`. CSE temporaries and final reduced expressions both receive optional `apply_substitution_dict()` processing, are expanded with `sp.expand()`, and are emitted through `printer.doprint()`.

This page owns the public `py_codegen()` and `PyCodeGen` contract. Helper pages own CSE helper and printer internals; JAX infrastructure pages own generated package assembly.

## Sources

- [nrpy/py_codegen.py](../../nrpy/py_codegen.py) - `PyCodeGen`, `py_codegen`, `printer`
- [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `apply_substitution_dict`

## See Also

- Parent: [Core APIs](index.md)
- Depends on: [CSE And Printer Support](helpers/cse-and-printer-support.md)
- Contrasts with: [C Codegen](c-codegen.md)
- See also: [JAX Project Generation Lifecycle](../infrastructures/jax/project-generation-lifecycle.md)
