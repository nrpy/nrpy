"""
Register diagnostics_integration_weights C function.

Determines the appropriate integration weights based on grid size inputs Nxx1 and Nxx2.
Since the grid is cell-centered, only midpoint methods are applicable, and the function selects
the highest possible order of integration based on divisibility of Nxx1 and Nxx2:
- If both Nxx1 and Nxx2 are divisible by 8, 8th-order weights are used.
- If both Nxx1 and Nxx2 are divisible by 4, 4th-order weights are used.
- Otherwise, defaults to 2nd-order weights.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc


def register_CFunction_diagnostics_integration_weights() -> None:
    """
    Register the C function to determine integration weights and stencil size based on Nxx1 and Nxx2.

    The integration weights are selected based on the fact that the grid is cell-centered,
    making only midpoint methods applicable. The function chooses the highest possible order:
    - If both Nxx1 and Nxx2 are divisible by 8, 8th-order weights are used.
    - If both Nxx1 and Nxx2 are divisible by 4, 4th-order weights are used.
    - Otherwise, 2nd-order weights are used.

    DocTests:
    >>> register_CFunction_diagnostics_integration_weights()
    """
    includes = ["BHaH_defines.h"]
    desc = """
Determines the appropriate integration weights based on grid size inputs Nxx1 and Nxx2.
Since the grid is cell-centered, only midpoint methods are applicable, and the function selects
the highest possible order of integration based on divisibility of Nxx1 and Nxx2:
- If both Nxx1 and Nxx2 are divisible by 8, 8th-order weights are used.
- If both Nxx1 and Nxx2 are divisible by 4, 4th-order weights are used.
- Otherwise, defaults to 2nd-order weights.
"""
    cfunc_type = "void"
    name = "diagnostics_integration_weights"
    params = "int Nxx1, int Nxx2, const REAL *restrict *weights, int *restrict weight_stencil_size"
    body = r"""
// Define weights for different orders
static const REAL weights_2nd_order[1] = {1.0};
static const REAL weights_4th_order[4] = {13.0 / 12.0, 11.0 / 12.0, 11.0 / 12.0, 13.0 / 12.0};
static const REAL weights_8th_order[8] = {295627.0 / 241920.0, 71329.0 / 241920.0, 17473.0 / 8960.0,   128953.0 / 241920.0, //
                                          128953.0 / 241920.0, 17473.0 / 8960.0,   71329.0 / 241920.0, 295627.0 / 241920.0};

// Default to 2nd-order weights
*weights = weights_2nd_order;
*weight_stencil_size = 1;

// Prefer higher order if possible by checking divisibility
if (Nxx1 % 8 == 0 && Nxx2 % 8 == 0) {
    // Use 8th-order weights
    *weights = weights_8th_order;
    *weight_stencil_size = 8;
} else if (Nxx1 % 4 == 0 && Nxx2 % 4 == 0) {
    // Use 4th-order weights
    *weights = weights_4th_order;
    *weight_stencil_size = 4;
}
"""

    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        raise RuntimeError(
            f"Doctest failed: {results.failed} of {results.attempted} test(s)"
        )
    print(f"Doctest passed: All {results.attempted} test(s) passed")
