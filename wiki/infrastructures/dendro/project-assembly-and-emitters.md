# Project Assembly And Generating Functions

> Explain how NRPy emits complete sibling Dendro-GR applications. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_bssn` and `nrpy.examples.dendro_fccz4` generate
`NRPy_BSSN_GR/` and `NRPy_fCCZ4_GR/` beside `BSSN_GR/`. Each directory is a
complete Dendro application, not a kernel library, mock host, or adapter around
Dendro-GR BSSN sources. With no explicit output directory, each example finds
the nearest ancestor containing `CMakeLists.txt` and `BSSN_GR/`.

## Detail

Each example registers the complete canonical state before parallel work. It
then queues every expensive Ricci, RHS, constraint, conversion, initial-data,
and runtime-service registrar for finite-difference orders 4, 6, and 8. One
`do_parallel_codegen()` call constructs and lowers all queued expressions.
The parent process does not reconstruct those expressions.

Worker results merge in sorted order into temporary registries. Unequal
duplicate definitions or final-state index mismatches abort before any parent
registry changes. Inexpensive emitters then write headers, parameter input,
context, executable entry point, checkpoint support, prototypes, and CMake.

One Python module owns each generated numerical operation. Python basename,
registrar suffix, CFunction name, and C++ basename correspond directly. Thus
`Ricci_eval.py` registers `Ricci_eval_order_6` and emits
`Ricci_eval_order_6.cpp`; `ADM_to_BSSN.py` emits order-specific
`ADM_to_BSSN` sources; `initial_data_lambdaU.py` owns only the separate
connection-initialization pass.

`CMakeLists.py` writes an explicit source list. It includes the generated
context, entry point, checkpoint support, local TwoPunctures implementation,
runtime services, conversions, projection, and every order-specific numerical
kernel. It uses no source glob, `bssn_common`, or source from `BSSN_GR/`.
Dendrolib is the shared runtime dependency.

## Sources

- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - BSSN registration wave and file emission.
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - fCCZ4 registration wave and file emission.
- [parallel_codegen.py](../../../nrpy/helpers/parallel_codegen.py) - worker execution and deterministic registry merge.
- [CMakeLists.py](../../../nrpy/infrastructures/Dendro/CMakeLists.py) - prototype and explicit CMake source emission.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - application entry point.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - generated Dendro runtime context.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Validated by: [Production Validation And Deferred Checks](validation-standalone-host-and-deferral-gates.md)
