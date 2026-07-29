"""
Register the C entry point for a single massive-particle integrator.

The helper emits a standard ``main(argc, argv)`` that forwards command-line
arguments to the registered integrator.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc


def main_single(integrator_name: str) -> None:
    """
    Register a C ``main`` that forwards to a single-particle integrator.

    :param integrator_name: C function taking ``argc, argv`` and returning ``int``.
    :raises ValueError: If ``integrator_name`` is not a valid C identifier.
    """
    if not integrator_name.isidentifier():
        raise ValueError(
            f"Invalid single-integrator C function name: {integrator_name}"
        )

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc="Forward the executable entry point to the selected massive-particle integrator.",
        cfunc_type="int",
        name="main",
        params="int argc, const char *argv[]",
        body=f"return {integrator_name}(argc, argv);",
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
