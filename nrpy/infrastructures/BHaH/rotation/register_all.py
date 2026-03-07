"""
Register all SO(3) rotation helper C functions in this directory.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from pathlib import Path

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures import BHaH


def register_CFunctions() -> None:
    """Register all public SO(3) helper C functions used by the BHaH rotation infrastructure."""
    BHaH.rotation.so3_build_R_from_hats.register_CFunction_so3_build_R_from_hats()
    BHaH.rotation.so3_axis_angle_to_R.register_CFunction_so3_axis_angle_to_R()
    BHaH.rotation.so3_apply_R_to_vector.register_CFunction_so3_apply_R_to_vector()
    BHaH.rotation.so3_apply_RT_to_vector.register_CFunction_so3_apply_RT_to_vector()
    BHaH.rotation.so3_relative_R_dst_from_src.register_CFunction_so3_relative_R_dst_from_src()
    BHaH.rotation.so3_left_multiply_hats_with_R.register_CFunction_so3_left_multiply_hats_with_R()


if __name__ == "__main__":
    import doctest
    import sys

    cfc.CFunction_dict.clear()
    par.set_parval_from_str("parallelization", "openmp")
    register_CFunctions()
    tests_dir = Path(__file__).resolve().parent / "tests"
    for func_name in (
        "so3_build_R_from_hats",
        "so3_axis_angle_to_R",
        "so3_apply_R_to_vector",
        "so3_apply_RT_to_vector",
        "so3_relative_R_dst_from_src",
        "so3_left_multiply_hats_with_R",
    ):
        generated_str = cfc.CFunction_dict[func_name].full_function
        trusted_path = tests_dir / f"{func_name}_{func_name}__openmp.c"
        trusted_str = trusted_path.read_text(encoding="utf-8")
        if trusted_str != generated_str:
            raise ValueError(f"Mismatch in trusted file: {trusted_path.name}")

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
