"""
Register the BHaH SO(3) helper that updates cumulative hats in place.

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


def register_CFunction_so3_left_multiply_hats_with_R() -> None:
    """
    Register the generated C function ``so3_left_multiply_hats_with_R``.

    High-level behavior:
    - Treats ``xhatU``, ``yhatU``, and ``zhatU`` as the columns of the current
      cumulative rotation matrix.
    - Applies ``DeltaR`` on the left so the updated hats represent
      ``R_new = DeltaR * R_old``.
    - Uses alias-safe temporary copies of all three input hat vectors.

    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Update cumulative hat vectors by left-multiplying with a rotation matrix.

Convention:
- ``R_old`` is the cumulative rotating-to-fixed matrix whose columns are the
  input hat vectors.
- ``R_new = DeltaR * R_old``.
- Each updated hat vector is one column of ``R_new``.
- C layout statement: ``DeltaR[i][j]`` is row ``i``, column ``j``.

This routine is alias-safe: it snapshots all input hat vectors before
overwriting them in place.

@param[in,out] xhatU Cumulative rotating-frame x-basis vector, updated in place.
@param[in,out] yhatU Cumulative rotating-frame y-basis vector, updated in place.
@param[in,out] zhatU Cumulative rotating-frame z-basis vector, updated in place.
@param[in] DeltaR Rotation matrix that left-multiplies the cumulative basis.
"""
    cfunc_type = "void"
    name = "so3_left_multiply_hats_with_R"
    params = "REAL xhatU[3], REAL yhatU[3], REAL zhatU[3], const REAL DeltaR[3][3]"

    DeltaR_sym: List[List[sp.Expr]] = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            DeltaR_sym[i][j] = cast(sp.Expr, sp.Symbol(f"DeltaR[{i}][{j}]"))
    xhat_in_sym = [cast(sp.Expr, sp.Symbol(f"xhat_in[{i}]")) for i in range(3)]
    yhat_in_sym = [cast(sp.Expr, sp.Symbol(f"yhat_in[{i}]")) for i in range(3)]
    zhat_in_sym = [cast(sp.Expr, sp.Symbol(f"zhat_in[{i}]")) for i in range(3)]
    xhat_expr = SO3Expressions.apply_R_to_vector(DeltaR_sym, xhat_in_sym)
    yhat_expr = SO3Expressions.apply_R_to_vector(DeltaR_sym, yhat_in_sym)
    zhat_expr = SO3Expressions.apply_R_to_vector(DeltaR_sym, zhat_in_sym)
    exprs = [*xhat_expr, *yhat_expr, *zhat_expr]
    lhses = [f"xhatU[{i}]" for i in range(3)]
    lhses += [f"yhatU[{i}]" for i in range(3)]
    lhses += [f"zhatU[{i}]" for i in range(3)]

    body = r"""
  // Snapshot the input hat vectors so the in-place update is alias-safe.
  const REAL xhat_in[3] = {xhatU[0], xhatU[1], xhatU[2]};
  const REAL yhat_in[3] = {yhatU[0], yhatU[1], yhatU[2]};
  const REAL zhat_in[3] = {zhatU[0], zhatU[1], zhatU[2]};
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
    register_CFunction_so3_left_multiply_hats_with_R()
    generated_str = cfc.CFunction_dict["so3_left_multiply_hats_with_R"].full_function
    validate_strings(
        generated_str, "so3_left_multiply_hats_with_R__openmp", file_ext="c"
    )

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
