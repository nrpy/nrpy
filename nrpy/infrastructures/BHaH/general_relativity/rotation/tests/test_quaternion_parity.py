"""Parity and convention gates for BSSN Cartesian-basis rotation helpers."""

from __future__ import annotations

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.BHaH.general_relativity.rotation.rotate_BSSN_Cartesian_basis import (
    assert_SO3_convention_in_text,
    register_CFunction_rotate_BSSN_Cartesian_basis,
    register_CFunction_rotate_BSSN_Cartesian_basis_by_R,
    verify_quaternion_interface_parity,
    verify_quaternion_interface_parity_randomized,
)


def run_quaternion_parity_and_convention_gates() -> bool:
    r"""
    Run deterministic parity checks and convention-text consistency gates.

    Doctests:
    >>> run_quaternion_parity_and_convention_gates()
    True
    """
    verify_quaternion_interface_parity()
    verify_quaternion_interface_parity_randomized()

    par.set_parval_from_str("parallelization", "openmp")
    cfc.CFunction_dict.clear()
    register_CFunction_rotate_BSSN_Cartesian_basis_by_R()
    register_CFunction_rotate_BSSN_Cartesian_basis()

    assert_SO3_convention_in_text(
        cfc.CFunction_dict["rotate_BSSN_Cartesian_basis_by_R"].full_function,
        "rotate_BSSN_Cartesian_basis_by_R",
    )
    assert_SO3_convention_in_text(
        cfc.CFunction_dict["rotate_BSSN_Cartesian_basis"].full_function,
        "rotate_BSSN_Cartesian_basis",
    )
    return True


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
