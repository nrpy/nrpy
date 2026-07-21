"""
Register the C entry point for single-photon integrator projects.

Single-ray generators keep descriptive C function names. This module supplies
the required C ``main`` symbol as a minimal forwarding entry point.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc


def main_single(integrator_name: str) -> None:
    """
    Register a C ``main`` that forwards to a single-ray integrator.

    :param integrator_name: C function receiving ``argc`` and ``argv``.
    :raises ValueError: If ``integrator_name`` is not a valid C identifier.
    """
    if not integrator_name.isidentifier():
        raise ValueError(
            f"Invalid single-integrator C function name: {integrator_name}"
        )

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc="Forward the executable entry point to the selected single-photon integrator.",
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
