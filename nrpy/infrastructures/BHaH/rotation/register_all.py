"""Register all public rotation C functions in this directory."""

from __future__ import annotations

from nrpy.infrastructures import BHaH


def register_CFunctions() -> None:
    """Register SO(3) rotation helpers for BHaH rotation call sites."""
    BHaH.rotation.so3_matrix_ops.register_CFunctions()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
