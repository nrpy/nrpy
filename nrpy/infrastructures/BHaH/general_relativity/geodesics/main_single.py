"""
Register the C entry point for a single-particle geodesic integrator.

The generated ``main(argc, argv)`` forwards command-line arguments to the
selected massive-particle or photon integrator.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc


def main_single(integrator_name: str) -> None:
    """
    Register a C ``main`` that forwards to a single-particle integrator.

    :param integrator_name: C function taking ``argc, argv`` and returning ``int``.
    :raises ValueError: If ``integrator_name`` is not a valid C identifier.

    Doctests:
    >>> import os
    >>> import tempfile
    >>> import nrpy.c_function as cfc
    >>> cache_dir = tempfile.TemporaryDirectory(dir=os.getcwd())
    >>> old_cache_home = os.environ.get("XDG_CACHE_HOME")
    >>> _ = os.environ.__setitem__("XDG_CACHE_HOME", cache_dir.name)
    >>> cfc.CFunction_dict.clear()
    >>> main_single("single_integrator_analytical")
    >>> generated = cfc.CFunction_dict["main"].full_function
    >>> "return single_integrator_analytical(argc, argv);" in generated
    True
    >>> main_single("not-a-C-identifier")
    Traceback (most recent call last):
        ...
    ValueError: Invalid single-integrator C function name: not-a-C-identifier
    >>> if old_cache_home is None:
    ...     _ = os.environ.pop("XDG_CACHE_HOME", None)
    ... else:
    ...     _ = os.environ.__setitem__("XDG_CACHE_HOME", old_cache_home)
    >>> cache_dir.cleanup()
    """
    if not integrator_name.isidentifier():
        raise ValueError(
            f"Invalid single-integrator C function name: {integrator_name}"
        )

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc="Forward the executable entry point to the selected single-particle geodesic integrator.",
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
