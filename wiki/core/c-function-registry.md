# C Function Registry

> Core route for generated C function objects and registration. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Core APIs](index.md)

## Summary

`nrpy.c_function` provides the in-memory registry used by NRPy infrastructure code to collect generated C functions before a project writer emits them. A `CFunction` stores the components of one generated C function, `CFunction_dict` maps registered names to `CFunction` instances, and `register_CFunction()` is the standard entry point that creates and stores them. BHaH infrastructure modules use this registry for executable `main()` generation and for library-style entry points.

## Detail

`CFunction` stores the output subdirectory, includes, optional `prefunc` and `postfunc` text, Doxygen-style description, return type, function name, parameter string, body, optional code-parameter include behavior, and infrastructure-specific metadata such as coordinate-system wrapper naming and Einstein Toolkit scheduling fields. Construction requires `name`, `desc`, and `body`; it validates include shape and forbids passing `griddata_struct` into coordinate-specific wrapper functions. It then builds a prototype, raw function text, and clang-formatted `full_function`.

The generated function text is assembled from includes, optional `prefunc`, a Doxygen comment derived from `desc`, the function signature, an optional `set_CodeParameters*.h` include, the body, an `END FUNCTION` comment, and optional `postfunc`. This makes a `CFunction` both a structured registry item and the source of the final emitted C function text.

`register_CFunction()` computes the actual name and subdirectory, including `__rfm__<CoordSystem>` suffixes when a BHaH coordinate-system wrapper is requested. It rejects duplicate names in `CFunction_dict` and stores the new `CFunction` under the actual generated name.

BHaH shows the registry in normal use. `register_CFunction_main_c()` checks that required functions are already present in `CFunction_dict`, then registers the generated `main` function. `register_CFunctions_bhah_lib()` registers library-facing initialization, evolution, diagnostics, and finalization functions through the same registration mechanism.

`GPUKernel` is a helper-side producer of `CFunction` text and launch calls. Use [Loop Kernel And Device Helpers](helpers/loop-kernel-and-device-helpers.md) for CUDA launch dictionaries, host/device wrapper choices, and GPU-kernel call-string behavior; this page remains the owner of the registry contract itself.

## Sources

- [nrpy/c_function.py](../../nrpy/c_function.py) - `CFunction`, `CFunction_dict`, `register_CFunction`, `function_name_and_subdir_with_CoordSystem`
- [nrpy/infrastructures/BHaH/main_c.py](../../nrpy/infrastructures/BHaH/main_c.py) - `register_CFunction_main_c`
- [nrpy/infrastructures/BHaH/bhah_lib.py](../../nrpy/infrastructures/BHaH/bhah_lib.py) - `register_CFunctions_bhah_lib`
- [coding_style.md](../../coding_style.md) - C function registration pattern

## See Also

- [Core APIs](index.md)
- [C Codegen](c-codegen.md)
- [Gridfunctions And Parameters](gridfunctions-and-parameters.md)
- [Loop Kernel And Device Helpers](helpers/loop-kernel-and-device-helpers.md)
- [Parallel Codegen Orchestration](helpers/parallel-codegen-orchestration.md)
