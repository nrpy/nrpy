"""
Register the BHaH SO(3) helper that computes a relative rotation matrix.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import List, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.equations.rotation.SO3_rotations import SO3Expressions
from nrpy.helpers.generic import validate_strings


def register_CFunction_so3_relative_R_dst_from_src() -> None:
    """
    Register the generated C function ``so3_relative_R_dst_from_src``.

    High-level behavior:
    - Interprets ``R_dst`` and ``R_src`` as cumulative rotating-to-fixed
      matrices for the destination and source frames.
    - Computes the relative rotation that maps source rotating-basis components
      to destination rotating-basis components.

    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Compute the relative rotation from a source rotating basis to a destination rotating basis.

Convention:
- ``R_dst`` and ``R_src`` map rotating-frame components to fixed-frame components.
- ``DeltaR_dst_from_src = R_dst^T R_src``.
- ``DeltaR_dst_from_src`` maps source rotating-basis components to destination
  rotating-basis components.
- C layout statement: ``R[i][j]`` and ``DeltaR[i][j]`` are row-major by index.

@param[in] R_dst Destination cumulative rotation matrix.
@param[in] R_src Source cumulative rotation matrix.
@param[out] DeltaR Relative rotation matrix ``R_dst^T R_src``.
"""
    cfunc_type = "void"
    name = "so3_relative_R_dst_from_src"
    params = "const REAL R_dst[3][3], const REAL R_src[3][3], REAL DeltaR[3][3]"

    R_dst_sym: List[List[sp.Expr]] = ixp.zerorank2()
    R_src_sym: List[List[sp.Expr]] = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            R_dst_sym[i][j] = cast(sp.Expr, sp.Symbol(f"R_dst[{i}][{j}]"))
            R_src_sym[i][j] = cast(sp.Expr, sp.Symbol(f"R_src[{i}][{j}]"))
    delta_expr = SO3Expressions.relative_rotation_dst_from_src(R_dst_sym, R_src_sym)
    exprs = [delta_expr[i][j] for i in range(3) for j in range(3)]
    lhses = [f"DeltaR[{i}][{j}]" for i in range(3) for j in range(3)]

    body = r"""
  // Form the destination-from-source relative rotation matrix.
""" + ccg.c_codegen(
        exprs,
        lhses,
        include_braces=False,
        verbose=False,
    )

    cfc.register_CFunction(
        subdirectory="rotation",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=False,
    )


if __name__ == "__main__":
    import doctest
    import sys

    cfc.CFunction_dict.clear()
    par.set_parval_from_str("parallelization", "openmp")
    register_CFunction_so3_relative_R_dst_from_src()
    generated_str = cfc.CFunction_dict["so3_relative_R_dst_from_src"].full_function
    validate_strings(generated_str, "so3_relative_R_dst_from_src__openmp", file_ext="c")

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
