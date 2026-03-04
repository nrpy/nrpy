"""
Register C helper that rotates Cartesian-basis BSSN vectors and tensors.

Global SO(3) convention used here:
- R maps rotating-frame components -> fixed-frame components.
- Columns: R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (all in fixed basis).
- v_fixed = R v_rot
- v_rot = R^T v_fixed
- T_fixed = R T_rot R^T
- T_rot = R^T T_fixed R
- DeltaR_dst_from_src = R_dst^T R_src

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
from nrpy.equations.rotation.SO3_rotations import SO3Expressions


def _rank2_exprs_lhses(
    rank2: List[List[sp.Expr]], varname: str
) -> tuple[List[sp.Expr], List[str]]:
    """
    Flatten a rank-2 tensor into expression and LHS-name lists.

    :param rank2: Rank-2 tensor.
    :param varname: Output C variable name.
    :return: Tuple of flattened expression and LHS-name lists.
    """
    exprs: List[sp.Expr] = []
    lhses: List[str] = []
    for i in range(3):
        for j in range(3):
            exprs.append(rank2[i][j])
            lhses.append(f"{varname}[{i}][{j}]")
    return exprs, lhses


def register_CFunction_rotate_BSSN_Cartesian_basis() -> None:
    """
    Register C function ``rotate_BSSN_Cartesian_basis``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_rotate_BSSN_Cartesian_basis()
    >>> generated_str = cfc.CFunction_dict["rotate_BSSN_Cartesian_basis"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "rotate_BSSN_Cartesian_basis_openmp", file_ext="c")
    """
    # Step 1: Ensure matrix-based helper is registered first.
    from nrpy.infrastructures.BHaH.general_relativity.rotation.rotate_BSSN_Cartesian_basis_by_R import (
        register_CFunction_rotate_BSSN_Cartesian_basis_by_R,
    )

    if "rotate_BSSN_Cartesian_basis_by_R" not in cfc.CFunction_dict:
        register_CFunction_rotate_BSSN_Cartesian_basis_by_R()

    # Step 2: Set C function metadata and interface.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""
@brief Rotate BSSN Cartesian-basis fields from axis-angle input.

This routine converts (nU,dphi) to DeltaR_dst_from_src and then calls
rotate_BSSN_Cartesian_basis_by_R().

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
  rotate_BSSN_Cartesian_basis_by_R(): only upper-triangular components
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
    name = "rotate_BSSN_Cartesian_basis"
    params = (
        "REAL vetU[3], REAL betU[3], REAL lambdaU[3], "
        "REAL hDD[3][3], REAL aDD[3][3], const REAL nU[3], const REAL dphi"
    )

    # Step 3: Build symbolic Rodrigues assignments from equations module.
    nU_unit_sym: List[sp.Expr] = [sp.Symbol("nx"), sp.Symbol("ny"), sp.Symbol("nz")]
    R_expr = SO3Expressions.rodrigues_matrix_from_unit_axis(
        nU_unit_sym, sp.Symbol("dphi")
    )
    R_exprs, R_lhses = _rank2_exprs_lhses(R_expr, "R")
    R_codegen = ccg.c_codegen(R_exprs, R_lhses, include_braces=False, verbose=False)

    # Step 4: Build Rodrigues matrix from (nU,dphi), then delegate to
    #         rotate_BSSN_Cartesian_basis_by_R().
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

  rotate_BSSN_Cartesian_basis_by_R(vetU, betU, lambdaU, hDD, aDD, R);
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
