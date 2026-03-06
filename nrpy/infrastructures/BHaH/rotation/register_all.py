"""Register all rotation C functions in this directory."""

from __future__ import annotations

from nrpy.infrastructures import BHaH


def register_CFunctions() -> None:
    r"""
    Register rotation helpers and their associated code parameters.

    The called registration functions register required cumulative-hat
    CodeParameters as part of their normal setup. This top-level aggregator also
    pulls in the GR rotation helpers so generated-C validation artifacts exist in
    both rotation test directories.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunctions()
    >>> for func_name in (
    ...     "unrotate_xCart_to_fixed_frame",
    ...     "unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame",
    ... ):
    ...     generated_str = cfc.CFunction_dict[func_name].full_function
    ...     validation_desc = f"{func_name}__openmp"
    ...     validate_strings(generated_str, validation_desc, file_ext="c")
    """
    BHaH.rotation.unrotate_xCart_to_fixed_frame.register_CFunction_unrotate_xCart_to_fixed_frame()
    BHaH.rotation.unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame.register_CFunction_unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
