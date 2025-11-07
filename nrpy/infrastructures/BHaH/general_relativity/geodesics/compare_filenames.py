"""
Generates a C helper function for sorting snapshot filenames.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def compare_filenames() -> None:
    """
    Generate and register a C helper function for qsort to sort filenames.

    This function generates the C code for `compare_filenames`, a comparison
    function required by the standard C library's `qsort`. It extracts the
    numerical timestamp from two filenames and compares them, enabling
    chronological sorting of snapshot files.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["stdio.h", "<stdlib.h>"]
    desc = r"""@brief Comparison function for qsort to sort snapshot filenames numerically.
    @details Extracts integer timestamps from filenames of the format
             "mass_blueprint_t_%%d.kdtree.bin" and compares them.
    @param a Void pointer to the first filename string.
    @param b Void pointer to the second filename string.
    @return -1 if a < b, 0 if a == b, 1 if a > b.
    """
    name = "compare_filenames"
    cfunc_type = "int"
    params = "const void *a, const void *b"

    # The body is algorithmic, not symbolic, so it is defined as a raw C string.
    body = r"""
    const char *str_a = *(const char **)a;
    const char *str_b = *(const char **)b;
    int num_a, num_b;

    // Parse the integer timestamp from each filename.
    sscanf(str_a, "mass_blueprint_t_%d.kdtree.bin", &num_a);
    sscanf(str_b, "mass_blueprint_t_%d.kdtree.bin", &num_b);

    // Return -1, 0, or 1 based on the comparison.
    return (num_a > num_b) - (num_a < num_b);
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        cfunc_type=cfunc_type,
        params=params,
        body=body,
    )
