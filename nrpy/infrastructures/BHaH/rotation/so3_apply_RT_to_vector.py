"""
Register the BHaH SO(3) helper that applies a transposed rotation matrix.

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


def register_CFunction_so3_apply_RT_to_vector() -> None:
    """
    Register the generated C function ``so3_apply_RT_to_vector``.

    High-level behavior:
    - Applies ``R^T`` to ``vU`` in place.
    - Uses an alias-safe temporary copy of the input vector before overwriting
      ``vU``.

    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Apply the transpose of a rotation matrix to a vector in place.

Convention:
- ``R`` maps rotating-frame components to fixed-frame components.
- ``R^T`` maps fixed-frame components to rotating-frame components.
- Vector update in Einstein notation: ``v^i_out = (R^T)^i{}_j v^j_in``.
- C layout statement: ``R[i][j]`` is row ``i``, column ``j``.

This routine is alias-safe: it snapshots the input vector before overwriting
``vU``.

@param[in] R Rotation matrix whose transpose is applied.
@param[in,out] vU Vector updated in place with ``R^T * vU``.
"""
    cfunc_type = "void"
    name = "so3_apply_RT_to_vector"
    params = "const REAL R[3][3], REAL vU[3]"

    R_sym: List[List[sp.Expr]] = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            R_sym[i][j] = cast(sp.Expr, sp.Symbol(f"R[{i}][{j}]"))
    v_in_sym = [cast(sp.Expr, sp.Symbol(f"v_in[{i}]")) for i in range(3)]
    v_expr = SO3Expressions.apply_RT_to_vector(R_sym, v_in_sym)
    exprs = [v_expr[i] for i in range(3)]
    lhses = [f"vU[{i}]" for i in range(3)]

    body = r"""
  // Snapshot the input vector so the in-place update is alias-safe.
  const REAL v_in[3] = {vU[0], vU[1], vU[2]};
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
    register_CFunction_so3_apply_RT_to_vector()
    generated_str = cfc.CFunction_dict["so3_apply_RT_to_vector"].full_function
    validate_strings(generated_str, "so3_apply_RT_to_vector__openmp", file_ext="c")

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
