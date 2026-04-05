"""
Register the BHaH SO(3) helper that maps one unit vector to another.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import List, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.equations.rotation.SO3_rotations import SO3Expressions
from nrpy.helpers.generic import validate_strings


def register_CFunction_so3_find_nU_and_dphi_from_unit_vectors() -> None:
    """
    Register the generated C function ``so3_find_nU_and_dphi_from_unit_vectors``.

    High-level behavior:
    - Treats ``aU`` and ``bU`` as unit vectors.
    - Recovers the axis-angle rotation that maps ``aU`` to ``bU``.
    - Handles parallel and antiparallel limits deterministically.

    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Recover axis-angle data that maps one unit vector to another.

Convention:
- ``aU`` and ``bU`` are interpreted as unit vectors in a common basis.
- ``nU`` and ``dphi`` satisfy ``bU = R(nU,dphi) aU``.
- Parallel input returns the deterministic identity pair
  ``nU=(1,0,0)``, ``dphi=0``.
- Antiparallel input returns ``dphi=pi`` with a deterministic axis
  orthogonal to ``aU``.

@param[in] aU Source unit vector.
@param[in] bU Destination unit vector.
@param[out] nU Unit rotation axis.
@param[out] dphi Rotation angle in radians.
"""
    cfunc_type = "void"
    name = "so3_find_nU_and_dphi_from_unit_vectors"
    params = "const REAL aU[3], const REAL bU[3], REAL nU[3], REAL *restrict dphi"

    aU_sym: List[sp.Expr] = [cast(sp.Expr, sp.Symbol(f"aU[{i}]")) for i in range(3)]
    bU_sym: List[sp.Expr] = [cast(sp.Expr, sp.Symbol(f"bU[{i}]")) for i in range(3)]
    cross_expr = SO3Expressions.cross_product3(aU_sym, bU_sym)
    dot_expr = SO3Expressions.dot_product3(aU_sym, bU_sym)
    cross_dot_codegen = ccg.c_codegen(
        [cross_expr[0], cross_expr[1], cross_expr[2], dot_expr],
        ["crossU[0]", "crossU[1]", "crossU[2]", "dot"],
        include_braces=False,
        verbose=False,
    )

    body = f"""
  REAL crossU[3];
  REAL dot;
{cross_dot_codegen}
  const REAL cross_norm = sqrt(crossU[0] * crossU[0] + crossU[1] * crossU[1] + crossU[2] * crossU[2]);

  if (dot > 1.0)
    dot = 1.0;
  if (dot < -1.0)
    dot = -1.0;

  if (cross_norm < 1e-14) {{
    if (dot > 0.0) {{
      nU[0] = 1.0;
      nU[1] = 0.0;
      nU[2] = 0.0;
      *dphi = 0.0;
      return;
    }} // END IF vectors are parallel.

    int basis_dir = 0;
    if (fabs(aU[1]) < fabs(aU[basis_dir]))
      basis_dir = 1;
    if (fabs(aU[2]) < fabs(aU[basis_dir]))
      basis_dir = 2;
    const REAL eU[3] = {{basis_dir == 0 ? 1.0 : 0.0, basis_dir == 1 ? 1.0 : 0.0, basis_dir == 2 ? 1.0 : 0.0}};
    const REAL axisU[3] = {{
        aU[1] * eU[2] - aU[2] * eU[1],
        aU[2] * eU[0] - aU[0] * eU[2],
        aU[0] * eU[1] - aU[1] * eU[0],
    }};
    const REAL axis_norm = sqrt(axisU[0] * axisU[0] + axisU[1] * axisU[1] + axisU[2] * axisU[2]);
    if (axis_norm < 1e-14) {{
      fprintf(stderr, "ERROR: %s failed to build a deterministic antiparallel rotation axis.\\n", __func__);
      exit(1);
    }} // END IF deterministic antiparallel axis is degenerate.
    nU[0] = axisU[0] / axis_norm;
    nU[1] = axisU[1] / axis_norm;
    nU[2] = axisU[2] / axis_norm;
    *dphi = M_PI;
    return;
  }} // END IF parallel/antiparallel branch.

  nU[0] = crossU[0] / cross_norm;
  nU[1] = crossU[1] / cross_norm;
  nU[2] = crossU[2] / cross_norm;
  *dphi = acos(dot);
"""

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
    register_CFunction_so3_find_nU_and_dphi_from_unit_vectors()
    generated_str = cfc.CFunction_dict[
        "so3_find_nU_and_dphi_from_unit_vectors"
    ].full_function
    validate_strings(
        generated_str,
        "so3_find_nU_and_dphi_from_unit_vectors__openmp",
        file_ext="c",
    )

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
