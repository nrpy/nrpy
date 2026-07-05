# Parallel Codegen Orchestration

> Helper leaf for generation-time multiprocessing registration, execution, and global-registry merge behavior. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Helper APIs](index.md)

## Summary

`nrpy.helpers.parallel_codegen` records code-generation function calls, executes them through Python multiprocessing when enabled, and merges each worker's NRPy environment back into the parent process. This is generation-time orchestration for Python code generation, not runtime OpenMP or CUDA execution in generated C/CUDA code.

## Detail

At import time, `parallel_codegen.py` registers two NRPy parameters: `enable_parallel_codegen`, defaulting to `False`, and `parallel_codegen_stage`, defaulting to `"register"`. These parameters control whether calls are collected and when the collected calls are executed.

`ParallelCodeGen` stores one deferred call. It splits a dotted path into `module_path` and `function_name`, then stores the argument dictionary that will later be passed as keyword arguments. `ParallelCodeGen_dict` is the module-level registry for these deferred calls. `register_func_call()` keys entries by `name + str(args)` and raises if the same key has already been registered.

`NRPyEnv()` snapshots the global dictionaries that generated code usually mutates: `par.glb_params_dict`, `par.glb_code_params_dict`, `cfc.CFunction_dict`, `pyfc.PyFunction_dict`, `gri.glb_gridfcs_dict`, and `par.glb_extras_dict`. The explicit return type records the same shape, including concrete gridfunction variants.

`deep_update()` recursively merges nested dictionaries. `unpack_NRPy_environment_dict()` walks the worker result dictionary and merges each returned environment into the parent process: parameters, code parameters, C functions, Python functions, and gridfunctions are updated directly, while extras are merged recursively through `deep_update()`.

`pcg_registration_phase()` is true only when `enable_parallel_codegen` is enabled and `parallel_codegen_stage` is `"register"`. Call sites can use that predicate to choose registration instead of immediate generation. `register_func_call()` only performs duplicate detection and storage; it does not itself check the phase predicate, so callers remain responsible for using it in the intended phase.

`get_nested_function()` dynamically imports the stored module path, walks dot-separated attributes inside the function name, verifies the result is callable, and raises import, attribute, or type errors with context when lookup fails. `parallel_function_call()` uses that lookup, calls the target with the stored keyword arguments, and expects the target to return an `NRPyEnv`-shaped tuple.

`wrapper_func()` is the process-pool worker wrapper. It receives a shared dictionary, key, and `ParallelCodeGen` value; runs `parallel_function_call()`, stores the returned environment under the key, logs elapsed time, and wraps worker failures in `RuntimeError` with the failing key.

`do_parallel_codegen()` is the execution boundary. If `enable_parallel_codegen` is false, it returns immediately. Otherwise it switches `parallel_codegen_stage` to `"codegen"`, forces the multiprocessing start method to `"fork"`, creates a `Manager().dict()` to collect worker environments, maps all registered calls through a `Pool`, then unpacks the collected environments back into the parent globals.

The forced `"fork"` start method is an implementation fact. The source comments explain that NRPy relies on inherited global state during parallel code generation and explicitly chooses `fork` for that behavior. This page does not treat that choice as a portability guarantee across Python platforms or libraries.

Because this helper mutates global registries, it is closely related to the C-function, gridfunction, parameter, and Python-function registries. The ownership remains separate: this page owns the multiprocessing orchestration and merge contract, while the cited registry modules define the meaning and validation rules of the stored objects.

## Sources

- [nrpy/helpers/parallel_codegen.py](../../../nrpy/helpers/parallel_codegen.py) - `ParallelCodeGen`, `ParallelCodeGen_dict`, `NRPyEnv`, `deep_update`, `unpack_NRPy_environment_dict`, `pcg_registration_phase`, `register_func_call`, `get_nested_function`, `parallel_function_call`, `wrapper_func`, `do_parallel_codegen`
- [nrpy/params.py](../../../nrpy/params.py) - `glb_params_dict`, `glb_code_params_dict`, `register_param`, `parval_from_str`, `set_parval_from_str`, `glb_extras_dict`
- [nrpy/c_function.py](../../../nrpy/c_function.py) - `CFunction_dict`
- [nrpy/py_function.py](../../../nrpy/py_function.py) - `PyFunction`, `PyFunction_dict`
- [nrpy/grid.py](../../../nrpy/grid.py) - `glb_gridfcs_dict`, `GridFunction`, `BHaHGridFunction`, `ETLegacyGridFunction`, `CarpetXGridFunction`

## See Also

- [Helper APIs](index.md)
- [Loop Kernel And Device Helpers](loop-kernel-and-device-helpers.md)
- [C Function Registry](../c-function-registry.md)
- [Python Function Registry](../python-function-registry.md)
- [Gridfunctions And Parameters](../gridfunctions-and-parameters.md)
- [BHaH Lifecycle](../../infrastructures/bhah-lifecycle.md)
