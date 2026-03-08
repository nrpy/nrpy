"""Register all GR basis transform C functions in this directory."""

from __future__ import annotations

from typing import Set

from nrpy.infrastructures import BHaH


def register_CFunctions(set_of_CoordSystems: Set[str]) -> None:
    r"""
    Register all general-relativity basis_transforms helper C functions.

    :param set_of_CoordSystems: Coordinate systems for which basis_transforms helpers
        should be registered.

    - each production module in ``general_relativity/basis_transform`` emits one public
      ``register_CFunction_*`` entry point,
    - the two single-point basis transforms register private
      ``__rfm__{CoordSystem}`` kernels,
    - ``rfm_wrapper_functions`` then emits the unsuffixed public runtime
      dispatchers used by pointwise and bulk routines.
    """
    BHaH.general_relativity.basis_transforms.basis_transform_BSSN_rfm_to_Cartesian_single_point.register_CFunction_basis_transform_BSSN_rfm_to_Cartesian_single_point(
        set_of_CoordSystems=set_of_CoordSystems
    )
    BHaH.general_relativity.basis_transforms.basis_transform_BSSN_Cartesian_to_rfm_single_point.register_CFunction_basis_transform_BSSN_Cartesian_to_rfm_single_point(
        set_of_CoordSystems=set_of_CoordSystems
    )
    BHaH.general_relativity.basis_transforms.BSSN_to_ADM_Cartesian.register_CFunction_BSSN_to_ADM_Cartesian()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
