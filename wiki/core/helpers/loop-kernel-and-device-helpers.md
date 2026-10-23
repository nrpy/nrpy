# Loop Kernel And Device Helpers

> Helper leaf for generic loop emitters, GPU kernel wrappers, and host/device code-generation utilities. · Status: contested
> Up: [Helper APIs](index.md)

## Summary

`nrpy.helpers.loop` emits generic C loop strings. `GPUKernel` builds CUDA-oriented kernel `CFunction` objects and their launch calls. `nrpy.helpers.parallelization.utilities` chooses host versus CUDA access, allocation, synchronization, wrapper, and loop-parameter snippets for generated kernels. BHaH-specific loop-region policy remains owned by the BHaH infrastructure page.

## Detail

`loop1D()` validates that its index, bounds, increment, and pragma inputs are strings, then returns a C `for`-loop header and matching footer comment. The pragma is optional, and an increment of `"1"` is emitted as `++`; other increments become `+= <increment>`.

`loop()` normalizes scalar inputs into lists and checks equal lengths for indices, bounds, increments, and pragmas. It returns `(header, footer)` without a body and one C string with a body. Tiling adds outer block loops and inner upper bounds `NRPYMIN(<upper>, <block> + <tile>)`; callers must define `NRPYMIN`. The length check excludes `tile_size`, so supply one tile entry per dimension; a scalar tile is not broadcast across nested loops. For positive constant increments, positive tile widths divisible by the increment preserve each dimension's untiled index set. Otherwise restarting at block boundaries can change the visited indices: bounds `[0,8)`, increment `2`, and tile width `3` visit `0,2,3,5,6`.

Claim evidence:
- Claim: Tiling requires an entry per dimension and compatible widths and increments; the helper neither checks tile-list length nor preserves arbitrary untiled index sets.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/loop.py](../../../nrpy/helpers/loop.py), `loop`
- Corroboration: `none available`; the named implementation decides this behavior.

`GPUKernel` wraps a generated kernel body in a `CFunction`. Its constructor records the body, decorators, parameter dictionary, generated name, launch dictionary, optional stream parameter, CUDA error-check policy, and BHaH thread-tiling suffix. A launch dictionary is required only for the exact default decorator string `decorators == "__global__"`. If `streamid_param` is enabled and the decorator string does not contain `__host__`, `streamid` is prepended to the generated argument dictionary. This is a substring test, not semantic classification of a host-only function.

`generate_launch_block()` derives CUDA launch setup only when a launch dictionary is present and the decorator string contains `__global__`. `threads_per_block` is padded to three dimensions or defaults to `32,1,1`. The launch dictionary must include `blocks_per_grid`: a nonempty list is padded to three dimensions and used directly, while an empty list requests computed grid dimensions from `params->Nxx_plus_2NGHOSTS*` and the thread counts. Optional `stream` and shared-memory (`sm`) entries add launch arguments, and the final launch settings string becomes the CUDA triple-chevron suffix used by `c_function_call()`.

When the active infrastructure is `BHaH`, `GPUKernel` replaces launch thread counts with `BHAH_THREADS_IN_<DIR>_DIR_<SUFFIX>` macros for every suffix, including `DEFAULT`. Non-default suffixes register values in `par.glb_extras_dict["DEVICE_THREAD_MACROS"]` and reject mismatches. `DEFAULT` skips registration and checking, so supplied numeric counts do not configure those macros. Use a non-default suffix for custom counts or configure the default macros separately.

Claim evidence:
- Claim: BHaH replaces numeric thread counts with macro names even for `DEFAULT`, which skips value registration and checking.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/parallelization/gpu_kernel.py](../../../nrpy/helpers/parallelization/gpu_kernel.py), `GPUKernel.__init__`
- Corroboration: `none available`; the named implementation decides this behavior.

`c_function_call()` constructs an ordinary function call or a CUDA launch from stored parameters and launch settings. Global kernels with launch settings use triple chevrons; other device/helper functions use direct calls. CUDA error checking is appended when enabled and the decorators contain `__global__`. Empty argument dictionaries without an injected stream argument are unsupported: unconditional removal of the final two characters can turn `foo(` into `fo);`. This also affects an empty-argument host wrapper. `GPU_Kernel` remains an alias for `GPUKernel`.

Claim status: contested; contradiction: CONTR-0019. See [CONTR-0019](../../contradictions.md#contr-0019).

Claim evidence:
- Claim: Empty-argument calls are malformed because `c_function_call()` unconditionally trims two characters before closing the call.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/parallelization/gpu_kernel.py](../../../nrpy/helpers/parallelization/gpu_kernel.py), `GPUKernel.c_function_call`
- Corroboration: `none available`; the named implementation decides this behavior.

`parallelization.utilities` centralizes small host-versus-CUDA string choices. `get_params_access()` returns `d_params[streamid].` for CUDA and `params->` otherwise; `get_commondata_access()` returns `d_commondata.` for CUDA and `commondata->` otherwise. Allocation, free, device synchronization, and error-check helpers similarly choose `cudaMalloc`/`malloc`, `cudaFree`/`free`, `cudaDeviceSynchronize()`/empty string, and CUDA error-check strings only for the CUDA path.

`generate_kernel_and_launch_code()` returns function-definition text and a call body. The CUDA path rewrites `params->` to `d_params[streamid].`, chooses the BHaH default launch dictionary when none is supplied, and constructs `<kernel_name>_gpu`. With custom launch settings, automatic stream-argument insertion depends on a `stream` entry and the decorator check in `GPUKernel`. If rewritten parameter accesses require `streamid`, ensure it is declared in the generated function. The CUDA path also ignores the supplied `cfunc_type`, using `GPUKernel`'s default `static void`; this argument cannot select a non-void device helper. The host path strips CUDA decorators, forwards `cfunc_type`, and returns a direct `<kernel_name>_host` call without CUDA error checks.

Claim evidence:
- Claim: CUDA parameter rewriting requires a declared `streamid`; CUDA ignores `cfunc_type`, whereas the host path forwards it.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/parallelization/utilities.py](../../../nrpy/helpers/parallelization/utilities.py), `generate_kernel_and_launch_code`
- Corroboration: `none available`; the named implementation decides this behavior.

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

- Parent: [Helper APIs](index.md)
- See also: [SIMD And Intrinsic Support](simd-and-intrinsic-support.md)
- See also: [Parallel Codegen Orchestration](parallel-codegen-orchestration.md)
- See also: [C Function Registry](../c-function-registry.md)
- See also: [Lifecycle And Project Assembly](../../infrastructures/bhah/lifecycle-and-project-assembly.md)
- See also: [Generated Output Boundaries](../../architecture/generated-output-boundaries.md)
