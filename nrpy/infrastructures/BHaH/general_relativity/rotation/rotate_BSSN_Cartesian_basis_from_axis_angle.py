"""
Generate C code for rotating BSSN Cartesian-basis fields from axis-angle input.

This module registers a wrapper C routine that:
1. takes an axis-angle pair ``(nU, dphi)``,
2. builds a 3x3 rotation matrix using Rodrigues' formula, and
3. delegates the actual field rotation to the matrix-based helper
   ``rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src``.

If the provided axis has near-zero norm, the wrapper uses an identity matrix
so the update is a deterministic no-op.

SO(3) convention used here:
- ``R`` maps rotating-frame components to fixed-frame components.
- ``R[:,0]=xhat``, ``R[:,1]=yhat``, ``R[:,2]=zhat`` (all in fixed basis).
- ``v_fixed = R v_rot`` and ``v_rot = R^T v_fixed``.
- ``T_fixed = R T_rot R^T`` and ``T_rot = R^T T_fixed R``.
- ``DeltaR_dst_from_src = R_dst^T R_src``.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
from nrpy.equations.rotation.SO3_rotations import SO3Expressions


def register_CFunction_rotate_BSSN_Cartesian_basis_from_axis_angle() -> None:
    """
    Register the generated C function ``rotate_BSSN_Cartesian_basis_from_axis_angle``.

    High-level behavior:
    - Interprets ``nU`` and ``dphi`` as an axis-angle rotation.
    - Normalizes ``nU`` when possible and builds ``R`` via Rodrigues.
    - Falls back to identity when ``|nU|`` is effectively zero.
    - Applies the rotation to ``vetU``, ``betU``, ``lambdaU``, ``hDD``, and
      ``aDD`` by calling the DeltaR-based helper.

    In-place contract:
    - All vector and tensor arguments are read and overwritten in place.
    - Symmetric tensor storage behavior is inherited from the callee:
      only upper-triangular entries are authoritative and updated.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_rotate_BSSN_Cartesian_basis_from_axis_angle()
    >>> generated_str = cfc.CFunction_dict["rotate_BSSN_Cartesian_basis_from_axis_angle"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "rotate_BSSN_Cartesian_basis_from_axis_angle_openmp", file_ext="c")
    """
    # Step 1: Set C function metadata and interface.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""
@brief Rotate BSSN Cartesian-basis fields from axis-angle input.

This routine converts (nU,dphi) to DeltaR_dst_from_src and then calls
rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src().

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout statement for helper calls: R[i][j] is row i, column j.
- Einstein notation for the callee update:
  v^i_dst = (DeltaR)^i{}_j v^j_src,
  T^dst_{ij} = (DeltaR)_i{}^k (DeltaR)_j{}^l T^src_{kl}.
- Symmetric tensor storage contract inherited from
  rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src(): only upper-triangular components
  (i <= j) are updated.

@param[in,out] vetU BSSN rescaled shift vector, rotated in place.
@param[in,out] betU BSSN rescaled driver vector, rotated in place.
@param[in,out] lambdaU BSSN conformal connection vector, rotated in place.
@param[in,out] hDD BSSN conformal metric perturbation tensor, rotated in place.
@param[in,out] aDD BSSN trace-free extrinsic curvature tensor, rotated in place.
@param[in] nU Rotation axis (need not be normalized).
@param[in] dphi Rotation angle in radians.
"""
    cfunc_type = "void"
    name = "rotate_BSSN_Cartesian_basis_from_axis_angle"
    params = (
        "REAL vetU[3], REAL betU[3], REAL lambdaU[3], "
        "REAL hDD[3][3], REAL aDD[3][3], const REAL nU[3], const REAL dphi"
    )

    # Step 2: Build symbolic Rodrigues assignments from equations module.
    nU_unit_sym: List[sp.Expr] = [sp.Symbol("nx"), sp.Symbol("ny"), sp.Symbol("nz")]
    R_expr = SO3Expressions.rodrigues_matrix_from_unit_axis(
        nU_unit_sym, sp.Symbol("dphi")
    )
    R_exprs = [R_expr[i][j] for i in range(3) for j in range(3)]
    R_lhses = [f"R[{i}][{j}]" for i in range(3) for j in range(3)]
    R_codegen = ccg.c_codegen(R_exprs, R_lhses, include_braces=False, verbose=False)

    # Step 3: Build Rodrigues matrix from (nU,dphi), then delegate to
    #         rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src().
    body = (
        r"""
      const REAL n0 = nU[0];
      const REAL n1 = nU[1];
      const REAL n2 = nU[2];
      const REAL nnorm = sqrt(n0 * n0 + n1 * n1 + n2 * n2);
    
      REAL R[3][3];
      if (nnorm < 1e-300) {
        // Deterministic fallback for degenerate axis input.
        R[0][0] = 1.0;
        R[0][1] = 0.0;
        R[0][2] = 0.0;
        R[1][0] = 0.0;
        R[1][1] = 1.0;
        R[1][2] = 0.0;
        R[2][0] = 0.0;
        R[2][1] = 0.0;
        R[2][2] = 1.0;
      } else {
        const REAL nx = n0 / nnorm;
        const REAL ny = n1 / nnorm;
        const REAL nz = n2 / nnorm;
    """
        + R_codegen
        + r"""
  }

  rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src(vetU, betU, lambdaU, hDD, aDD, R);
"""
    )
    cfc.register_CFunction(
        subdirectory="general_relativity/rotation",
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

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
