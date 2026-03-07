"""
Register the BHaH SO(3) helper that builds a rotation matrix from hat vectors.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.equations.rotation.SO3_rotations import SO3Expressions
from nrpy.helpers.generic import validate_strings


def register_CFunction_so3_build_R_from_hats() -> None:
    """
    Register the generated C function ``so3_build_R_from_hats``.

    High-level behavior:
    - Interprets ``xhatU``, ``yhatU``, and ``zhatU`` as the cumulative rotating
      basis vectors expressed in the fixed basis.
    - Packs those hat vectors as the three columns of ``R``.

    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Build a cumulative rotation matrix from cumulative hat vectors.

Convention:
- ``R`` maps rotating-frame components to fixed-frame components.
- ``R[:,0]=xhat``, ``R[:,1]=yhat``, ``R[:,2]=zhat``.
- C layout statement: ``R[i][j]`` is row ``i``, column ``j``.

@param[in] xhatU Cumulative rotating-frame x-basis vector in fixed components.
@param[in] yhatU Cumulative rotating-frame y-basis vector in fixed components.
@param[in] zhatU Cumulative rotating-frame z-basis vector in fixed components.
@param[out] R Rotation matrix assembled from the supplied hat vectors.
"""
    cfunc_type = "void"
    name = "so3_build_R_from_hats"
    params = (
        "const REAL xhatU[3], const REAL yhatU[3], const REAL zhatU[3], REAL R[3][3]"
    )

    xhat_sym = [cast(sp.Expr, sp.Symbol(f"xhatU[{i}]")) for i in range(3)]
    yhat_sym = [cast(sp.Expr, sp.Symbol(f"yhatU[{i}]")) for i in range(3)]
    zhat_sym = [cast(sp.Expr, sp.Symbol(f"zhatU[{i}]")) for i in range(3)]
    R_expr = SO3Expressions.build_rotation_matrix_from_hats(
        xhat_sym, yhat_sym, zhat_sym
    )
    exprs = [R_expr[i][j] for i in range(3) for j in range(3)]
    lhses = [f"R[{i}][{j}]" for i in range(3) for j in range(3)]

    body = r"""
  // Assemble the matrix column-by-column from cumulative hat vectors.
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
    register_CFunction_so3_build_R_from_hats()
    generated_str = cfc.CFunction_dict["so3_build_R_from_hats"].full_function
    validate_strings(generated_str, "so3_build_R_from_hats__openmp", file_ext="c")

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
