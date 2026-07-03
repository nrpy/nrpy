"""
Register the BHaH SO(3) helper that applies a rotation matrix to a rank-2 tensor.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import List, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp
from nrpy.equations.rotation.SO3_rotations import SO3Expressions


def register_CFunction_so3_apply_R_to_tensorDD() -> None:
    """
    Register the generated C function ``so3_apply_R_to_tensorDD``.

    High-level behavior:
    - Applies ``R`` to both indices of ``tDD`` in place.
    - Uses an alias-safe temporary copy of all input tensor entries before
      overwriting ``tDD``.
    - Writes all 9 tensor entries.

    """
    includes = ["BHaH_defines.h"]
    desc = r"""
Apply a rotation matrix to a rank-2 covariant tensor in place.

Convention:
- ``R`` maps rotating-frame components to fixed-frame components.
- Tensor update in Einstein notation: ``T_out_ij = R_i^k R_j^l T_in_kl``.
- C layout statement: ``R[i][j]`` is row ``i``, column ``j``.

This routine is alias-safe: it snapshots the input tensor before overwriting
``tDD``. It writes all 9 components; callers may repack symmetric storage after
the call.

@param[in] R Rotation matrix.
@param[in,out] tDD Tensor updated in place with ``R * tDD * R^T``.
"""
    cfunc_type = "void"
    name = "so3_apply_R_to_tensorDD"
    params = "const REAL R[3][3], REAL tDD[3][3]"

    R_sym: List[List[sp.Expr]] = ixp.zerorank2()
    t_in_sym: List[List[sp.Expr]] = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            R_sym[i][j] = cast(sp.Expr, sp.Symbol(f"R[{i}][{j}]"))
            t_in_sym[i][j] = cast(sp.Expr, sp.Symbol(f"t_in[{i}][{j}]"))
    t_expr = SO3Expressions.apply_R_to_tensorDD(R_sym, t_in_sym)
    exprs = [t_expr[i][j] for i in range(3) for j in range(3)]
    lhses = [f"tDD[{i}][{j}]" for i in range(3) for j in range(3)]

    body = r"""
  // Snapshot the input tensor so the in-place update is alias-safe.
  const REAL t_in[3][3] = {
      {tDD[0][0], tDD[0][1], tDD[0][2]},
      {tDD[1][0], tDD[1][1], tDD[1][2]},
      {tDD[2][0], tDD[2][1], tDD[2][2]},
  };
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

    import nrpy.params as par
    from nrpy.helpers.generic import clang_format, validate_strings

    cfc.CFunction_dict.clear()
    par.set_parval_from_str("parallelization", "openmp")
    register_CFunction_so3_apply_R_to_tensorDD()
    generated_str = cfc.CFunction_dict["so3_apply_R_to_tensorDD"].full_function
    validate_strings(
        clang_format(generated_str), "so3_apply_R_to_tensorDD__openmp", file_ext="c"
    )

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
