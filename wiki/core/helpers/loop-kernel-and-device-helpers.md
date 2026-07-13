# Loop Kernel And Device Helpers

> Helper leaf for generic loop emitters, GPU kernel wrappers, and host/device code-generation utilities. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Helper APIs](index.md)

## Summary

`nrpy.helpers.loop` emits generic C loop strings. `GPUKernel` builds CUDA-oriented kernel `CFunction` objects and their launch calls. `nrpy.helpers.parallelization.utilities` chooses host versus CUDA access, allocation, synchronization, wrapper, and loop-parameter snippets for generated kernels. BHaH-specific loop-region policy remains owned by the BHaH infrastructure page.

## Detail

`loop1D()` validates that its index, bounds, increment, and pragma inputs are strings, then returns a C `for`-loop header and matching footer comment. The pragma is optional, and an increment of `"1"` is emitted as `++`; other increments become `+= <increment>`.

`loop()` normalizes scalar loop inputs into lists, requires the index, lower-bound, upper-bound, increment, and pragma lists to have the same length, and emits nested loops. Without a body it returns `(header, footer)`. With `loop_body`, it returns one full loop string. When `tile_size` is supplied, it adds an outer block loop with a `B` suffix and uses `NRPYMIN(<upper>, <block> + <tile>)` for the inner-loop upper bound; the helper emits that macro use but does not define the macro.

`GPUKernel` wraps a generated kernel body in a `CFunction`. Its constructor records the body, decorators, parameter dictionary, generated name, launch dictionary, optional stream parameter, CUDA error-check policy, and BHaH thread-tiling suffix. A launch dictionary is required only for the exact default decorator string `decorators == "__global__"`. If `streamid_param` is enabled and the decorator string does not contain `__host__`, `streamid` is prepended to the generated argument dictionary. This is a substring test, not semantic classification of a host-only function.

`generate_launch_block()` derives CUDA launch setup only when a launch dictionary is present and the decorator string contains `__global__`. `threads_per_block` is padded to three dimensions or defaults to `32,1,1`. The launch dictionary must include `blocks_per_grid`: a nonempty list is padded to three dimensions and used directly, while an empty list requests computed grid dimensions from `params->Nxx_plus_2NGHOSTS*` and the thread counts. Optional `stream` and shared-memory (`sm`) entries add launch arguments, and the final launch settings string becomes the CUDA triple-chevron suffix used by `c_function_call()`.

When the active infrastructure is `BHaH`, `GPUKernel` also participates in BHaH device-thread macro bookkeeping. Non-default `thread_tiling_macro_suffix` values are checked against `par.glb_extras_dict["DEVICE_THREAD_MACROS"]`, mismatches raise, and the thread counts are replaced by `BHAH_THREADS_IN_<DIR>_DIR_<SUFFIX>` macro names. This page records that helper integration point but does not own the BHaH lifecycle.

`c_function_call()` emits either an ordinary function call or a CUDA kernel launch call from the stored parameters and launch settings. Default `__global__` kernels with launch settings get CUDA triple-chevron calls; CUDA helper or device functions whose decorators do not include `__global__` can produce direct calls without triple-chevron launch settings. CUDA error checking is appended only when `cuda_check_error` remains enabled and the decorators include `__global__`. `GPU_Kernel` is a backward-compatible alias for `GPUKernel`.

`parallelization.utilities` centralizes small host-versus-CUDA string choices. `get_params_access()` returns `d_params[streamid].` for CUDA and `params->` otherwise; `get_commondata_access()` returns `d_commondata.` for CUDA and `commondata->` otherwise. Allocation, free, device synchronization, and error-check helpers similarly choose `cudaMalloc`/`malloc`, `cudaFree`/`free`, `cudaDeviceSynchronize()`/empty string, and CUDA error-check strings only for the CUDA path.

`generate_kernel_and_launch_code()` builds a generated kernel definition plus its launch body. For CUDA it rewrites `params->` accesses in the kernel body to the CUDA params access form, uses the BHaH CUDA default launch dictionary when no launch dictionary is provided, creates a `<kernel_name>_gpu` `GPUKernel`, and returns that kernel's full function text as `prefunc` with its launch block and call as `body`. For non-CUDA paths it strips CUDA decorators, creates a `<kernel_name>_host` function without launch settings or CUDA error checks, and returns a direct host call body.

`get_loop_parameters()` emits common per-kernel declarations for generated loop bodies: `Nxx_plus_2NGHOSTS*`, inverse grid spacings, optional SIMD/CUDA intrinsic constants, and, on CUDA, `tid*` and `stride*` values based on block, grid, and thread indices. CUDA output replaces the `SIMD` token family with `CUDA`.

`nrpy.infrastructures.BHaH.simple_loop` is a one-hop consumer of the generic loop helper. It selects BHaH loop regions, OpenMP pragmas, CUDA thread starts and strides, optional SIMD increments, coordinate reads, and reference-metric reads before calling `nrpy.helpers.loop.loop()`. Those BHaH loop-region choices stay in the infrastructure lifecycle page, not in this helper leaf.

## Sources

- [nrpy/helpers/loop.py](../../../nrpy/helpers/loop.py) - `loop1D`, `loop`
- [nrpy/helpers/parallelization/gpu_kernel.py](../../../nrpy/helpers/parallelization/gpu_kernel.py) - `GPUKernel`, `generate_launch_block`, `c_function_call`, `GPU_Kernel`
- [nrpy/helpers/parallelization/utilities.py](../../../nrpy/helpers/parallelization/utilities.py) - `get_params_access`, `get_commondata_access`, `get_memory_malloc_function`, `get_memory_free_function`, `get_check_errors_str`, `get_device_sync_function`, `generate_kernel_and_launch_code`, `get_loop_parameters`
- [nrpy/infrastructures/BHaH/simple_loop.py](../../../nrpy/infrastructures/BHaH/simple_loop.py) - `simple_loop`
- [nrpy/c_function.py](../../../nrpy/c_function.py) - `CFunction`
- [nrpy/params.py](../../../nrpy/params.py) - `parval_from_str`, `glb_extras_dict`

## See Also

- [Helper APIs](index.md)
- [SIMD And Intrinsic Support](simd-and-intrinsic-support.md)
- [Parallel Codegen Orchestration](parallel-codegen-orchestration.md)
- [C Function Registry](../c-function-registry.md)
- [Lifecycle And Project Assembly](../../infrastructures/bhah/lifecycle-and-project-assembly.md)
- [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
