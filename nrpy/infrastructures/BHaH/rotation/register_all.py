"""Register all rotation C functions in this directory."""

from nrpy.infrastructures import BHaH


def register_CFunctions() -> None:
    """
    Register rotation helpers and their associated code parameters.

    The called registration functions register required cumulative-hat
    CodeParameters as part of their normal setup.
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
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
