# C Function Registry

> Core route for generated C function objects and registration. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Core APIs](index.md)

## Summary

`nrpy.c_function` provides the in-memory registry used by NRPy infrastructure code to collect generated C functions before a project writer emits them. A `CFunction` stores the components of one generated C function, `CFunction_dict` maps registered names to `CFunction` instances, and `register_CFunction()` is the standard entry point that creates and stores them. BHaH infrastructure modules use this registry for executable `main()` generation and for library-style entry points.

## Detail

`CFunction` stores the output subdirectory, SIMD flag, includes, optional `prefunc` and `postfunc` text, Doxygen-style description, optional decorators, return type, function name, parameter string, body, optional code-parameter include behavior, coordinate-system wrapper metadata, and Einstein Toolkit metadata. Construction requires `name`, `desc`, and `body`; missing any of those fields raises a `ValueError`. It also validates that `includes` is a list and forbids passing `griddata_struct` into coordinate-specific wrapper functions.

Include handling preserves the caller's intent. During `generate_full_function()`, include entries containing `<` are emitted as angle-bracket includes such as `#include <math.h>`, while all other include strings become quoted includes such as `#include "set_CodeParameters.h"`. Non-string include entries raise a `TypeError`.

Comment formatting is centralized in `prefix_with_star()`. It strips only enclosing newlines, trims trailing whitespace per line, normalizes existing leading ` * ` markers, preserves meaningful indentation, and prefixes each description line with ` *`. `generate_full_function()` wraps that result in a Doxygen block before the function signature.

When `include_CodeParameters_h` is true, the generated body begins with a code-parameter include. It uses `set_CodeParameters-simd.h` if `enable_simd` is true or if the function body contains `SIMD_WIDTH`; otherwise it uses `set_CodeParameters.h`. This include is inserted inside the function body before the stripped body text.

`generate_full_function()` creates three outputs and stores them on the object: `function_prototype`, `raw_function`, and clang-formatted `full_function`. The prototype includes decorators, return type, generated name, and parameters. The raw function assembles includes, optional prefunc text, Doxygen comment, signature, optional code-parameter include, stripped body, an `// END FUNCTION: <name>` marker, and optional postfunc text. `full_function` is `clang_format(raw_function)`.

Coordinate-system wrappers use `function_name_and_subdir_with_CoordSystem()`. When `CoordSystem_for_wrapper_func` is set, the generated function goes under `subdirectory/<CoordSystem>` and its actual registered name becomes `<name>__rfm__<CoordSystem>`; otherwise the original subdirectory and name are used. Constructor validation rejects coordinate-specific wrapper parameters that contain `griddata_struct`.

Einstein Toolkit metadata is stored directly on the `CFunction`: `ET_thorn_name`, `ET_schedule_bins_entries`, `ET_current_thorn_CodeParams_used`, and `ET_other_thorn_CodeParams_used`. These fields let ET infrastructure code later know which thorn owns the function, which schedule.ccl entries it needs, which CodeParameters belong to the current thorn's param.ccl, and which referenced CodeParameters belong to other thorns.

`register_CFunction()` computes the actual name and subdirectory, including any `__rfm__<CoordSystem>` suffix, rejects duplicate actual names already present in `CFunction_dict`, constructs the `CFunction`, and stores it under the actual generated name. BHaH shows the registry in normal use: `register_CFunction_main_c()` checks that required functions are already present in `CFunction_dict`, then registers the generated `main` function; `register_CFunctions_bhah_lib()` registers library-facing initialization, evolution, diagnostics, and finalization functions through the same registration mechanism.

`GPUKernel` is a helper-side producer of `CFunction` text and launch calls. Use [Loop Kernel And Device Helpers](helpers/loop-kernel-and-device-helpers.md) for CUDA launch dictionaries, host/device wrapper choices, and GPU-kernel call-string behavior; this page remains the owner of the registry contract itself.

## Sources

- [nrpy/c_function.py](../../nrpy/c_function.py) - `CFunction`, `CFunction_dict`, `prefix_with_star`, `generate_full_function`, `register_CFunction`, `function_name_and_subdir_with_CoordSystem`
- [nrpy/infrastructures/BHaH/main_c.py](../../nrpy/infrastructures/BHaH/main_c.py) - `register_CFunction_main_c`
- [nrpy/infrastructures/BHaH/bhah_lib.py](../../nrpy/infrastructures/BHaH/bhah_lib.py) - `register_CFunctions_bhah_lib`

## See Also

- [Core APIs](index.md)
- [C Codegen](c-codegen.md)
- [Gridfunctions And Parameters](gridfunctions-and-parameters.md)
- [Loop Kernel And Device Helpers](helpers/loop-kernel-and-device-helpers.md)
- [Parallel Codegen Orchestration](helpers/parallel-codegen-orchestration.md)
