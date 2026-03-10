"""
Register the BHaH SO(3) helper that converts axis-angle data to a matrix.

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


def register_CFunction_so3_axis_angle_to_R() -> None:
    """
    Register the generated C function ``so3_axis_angle_to_R``.

    High-level behavior:
    - Interprets ``nU`` and ``dphi`` as an axis-angle rotation.
    - Normalizes ``nU`` when possible and builds ``R`` via Rodrigues' formula.
    - Falls back to the identity matrix when the axis norm or angle is
      effectively zero.

    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Build a rotation matrix from axis-angle input.

Convention:
- ``R`` maps rotating-frame components to fixed-frame components.
- ``v_fixed = R v_rot`` and ``v_rot = R^T v_fixed``.
- C layout statement: ``R[i][j]`` is row ``i``, column ``j``.

If ``|dphi|`` is effectively zero or ``|nU|`` is negligible, this routine
returns the identity matrix as a deterministic no-op rotation.

@param[in] nU Rotation axis, which need not be normalized.
@param[in] dphi Rotation angle in radians.
@param[out] R Rotation matrix computed from the axis-angle pair.
"""
    cfunc_type = "void"
    name = "so3_axis_angle_to_R"
    params = "const REAL nU[3], const REAL dphi, REAL R[3][3]"

    nU_unit_sym: List[sp.Expr] = [
        cast(sp.Expr, sp.Symbol("nx")),
        cast(sp.Expr, sp.Symbol("ny")),
        cast(sp.Expr, sp.Symbol("nz")),
    ]
    R_expr = SO3Expressions.rodrigues_matrix_from_unit_axis(
        nU_unit_sym, cast(sp.Expr, sp.Symbol("dphi"))
    )
    exprs = [R_expr[i][j] for i in range(3) for j in range(3)]
    lhses = [f"R[{i}][{j}]" for i in range(3) for j in range(3)]
    R_codegen = ccg.c_codegen(
        exprs,
        lhses,
        include_braces=False,
        verbose=False,
    )

    body = (
        r"""
  const REAL n0 = nU[0];
  const REAL n1 = nU[1];
  const REAL n2 = nU[2];
  const REAL nnorm = sqrt(n0 * n0 + n1 * n1 + n2 * n2);
  const REAL safe_nnorm = (nnorm < 1e-300) ? 1.0 : nnorm;

  // Normalize the input axis before evaluating Rodrigues' formula.
  const REAL nx = n0 / safe_nnorm;
  const REAL ny = n1 / safe_nnorm;
  const REAL nz = n2 / safe_nnorm;
"""
        + R_codegen
        + r"""
  if (fabs(dphi) < 1e-15 || nnorm < 1e-300) {
    // Deterministic identity fallback for degenerate axis-angle input.
    R[0][0] = 1.0;
    R[0][1] = 0.0;
    R[0][2] = 0.0;
    R[1][0] = 0.0;
    R[1][1] = 1.0;
    R[1][2] = 0.0;
    R[2][0] = 0.0;
    R[2][1] = 0.0;
    R[2][2] = 1.0;
  } // END IF: degenerate axis-angle input
"""
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
    register_CFunction_so3_axis_angle_to_R()
    generated_str = cfc.CFunction_dict["so3_axis_angle_to_R"].full_function
    validate_strings(generated_str, "so3_axis_angle_to_R__openmp", file_ext="c")

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
