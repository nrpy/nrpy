"""
Register C helpers that rotate Cartesian-basis BSSN vectors and tensors.

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

from typing import Dict, List

import sympy as sp

# Step 1: Import core NRPy modules needed for symbolic expression construction
#         and C code generation.
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp
from nrpy.equations.rotation.SO3_rotations import SO3Expressions

# Step 2: Import SO(3) infrastructure helper registrations consumed by this module.
from nrpy.infrastructures.BHaH.rotation.SO3_matrix_ops import (
    register_CFunction_SO3_apply_R_to_tensorDD,
    register_CFunction_SO3_apply_R_to_vector,
)


def _codegen_no_braces(exprs: list[sp.Expr], lhses: list[str]) -> str:
    """
    Generate C code from expressions without wrapping braces.

    :param exprs: Expression list.
    :param lhses: Output variable-name list.
    :return: C code assignment string.
    """
    return ccg.c_codegen(exprs, lhses, include_braces=False, verbose=False)


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


def _evaluate_rank1_expression(
    rank1_expr: List[sp.Expr], substitutions: Dict[sp.Symbol, sp.Expr]
) -> List[sp.Expr]:
    """
    Evaluate a rank-1 symbolic expression under substitutions.

    :param rank1_expr: Symbolic rank-1 expression.
    :param substitutions: Mapping from symbols to numeric values.
    :return: Numeric rank-1 value list.
    """
    out = ixp.zerorank1()
    for i in range(3):
        out[i] = sp.N(rank1_expr[i].subs(substitutions), 50)
    return out


def _evaluate_rank2_expression(
    rank2_expr: List[List[sp.Expr]], substitutions: Dict[sp.Symbol, sp.Expr]
) -> List[List[sp.Expr]]:
    """
    Evaluate a rank-2 symbolic expression under substitutions.

    :param rank2_expr: Symbolic rank-2 expression.
    :param substitutions: Mapping from symbols to numeric values.
    :return: Numeric rank-2 value list.
    """
    out = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            out[i][j] = sp.N(rank2_expr[i][j].subs(substitutions), 50)
    return out


def _rodrigues_matrix_eval(
    nU_input: List[sp.Expr], dphi_input: sp.Expr
) -> List[List[sp.Expr]]:
    """
    Evaluate Rodrigues matrix expression numerically from axis-angle input.

    :param nU_input: Rotation-axis components.
    :param dphi_input: Rotation angle in radians.
    :return: Numeric 3x3 rotation matrix.
    """
    axis_sq = (
        nU_input[0] * nU_input[0]
        + nU_input[1] * nU_input[1]
        + nU_input[2] * nU_input[2]
    )
    if axis_sq < sp.sympify("1e-600"):
        return [
            [sp.Integer(1), sp.Integer(0), sp.Integer(0)],
            [sp.Integer(0), sp.Integer(1), sp.Integer(0)],
            [sp.Integer(0), sp.Integer(0), sp.Integer(1)],
        ]

    nU_sym = ixp.declarerank1("nU")
    dphi_sym = sp.Symbol("dphi")
    R_expr = SO3Expressions.rodrigues_matrix_from_axis_angle(nU_sym, dphi_sym)

    substitutions: Dict[sp.Symbol, sp.Expr] = {dphi_sym: sp.sympify(dphi_input)}
    for i in range(3):
        substitutions[nU_sym[i]] = sp.sympify(nU_input[i])
    return _evaluate_rank2_expression(R_expr, substitutions)


def _apply_rotation_to_vector_eval(
    R_input: List[List[sp.Expr]], v_input: List[sp.Expr]
) -> List[sp.Expr]:
    """
    Evaluate v_dst = R v_src from symbolic expressions.

    :param R_input: Numeric rotation matrix.
    :param v_input: Numeric input vector.
    :return: Numeric rotated vector.
    """
    R_sym = ixp.declarerank2("R", symmetry="nosym")
    v_sym = ixp.declarerank1("v")
    v_expr = SO3Expressions.apply_R_to_vector(R_sym, v_sym)

    substitutions: Dict[sp.Symbol, sp.Expr] = {}
    for i in range(3):
        substitutions[v_sym[i]] = sp.sympify(v_input[i])
        for j in range(3):
            substitutions[R_sym[i][j]] = sp.sympify(R_input[i][j])
    return _evaluate_rank1_expression(v_expr, substitutions)


def _apply_rotation_to_tensor_eval(
    R_input: List[List[sp.Expr]], t_input: List[List[sp.Expr]]
) -> List[List[sp.Expr]]:
    """
    Evaluate T_dst = R T_src R^T from symbolic expressions.

    :param R_input: Numeric rotation matrix.
    :param t_input: Numeric covariant rank-2 tensor.
    :return: Numeric rotated tensor.
    """
    R_sym = ixp.declarerank2("R", symmetry="nosym")
    t_sym = ixp.declarerank2("t", symmetry="nosym")
    t_expr = SO3Expressions.apply_R_to_tensorDD(R_sym, t_sym)

    substitutions: Dict[sp.Symbol, sp.Expr] = {}
    for i in range(3):
        for j in range(3):
            substitutions[R_sym[i][j]] = sp.sympify(R_input[i][j])
            substitutions[t_sym[i][j]] = sp.sympify(t_input[i][j])
    return _evaluate_rank2_expression(t_expr, substitutions)


def verify_quaternion_interface_parity() -> None:
    r"""
    Verify axis-angle matrix convention against quaternion tensor rotation interface.

    Doctests:
    >>> from nrpy.equations.quaternion_rotations.tensor_rotation import rotate
    >>> axis = [sp.Rational(3, 10), -sp.Rational(2, 5), sp.Rational(1, 2)]
    >>> dphi = sp.Rational(6, 5)
    >>> vec = [1, -2, sp.Rational(1, 2)]
    >>> ten = [[2, sp.Rational(1, 10), -sp.Rational(3, 10)], [sp.Rational(1, 10), sp.Rational(3, 2), sp.Rational(7, 10)], [-sp.Rational(3, 10), sp.Rational(7, 10), sp.Rational(4, 5)]]
    >>> R = _rodrigues_matrix_eval(axis, dphi)
    >>> vec_by_R = _apply_rotation_to_vector_eval(R, vec)
    >>> vec_by_q = [complex(v).real for v in rotate(vec, axis, dphi)]
    >>> max(abs(a - b) for a, b in zip(vec_by_R, vec_by_q)) < 1e-13
    True
    >>> ten_by_R = _apply_rotation_to_tensor_eval(R, ten)
    >>> ten_by_q = [[complex(v).real for v in row] for row in rotate(ten, axis, dphi)]
    >>> max(abs(ten_by_R[i][j] - ten_by_q[i][j]) for i in range(3) for j in range(3)) < 1e-13
    True
    >>> # Near-pi case
    >>> dphi = sp.pi - sp.sympify("1e-10")
    >>> R = _rodrigues_matrix_eval(axis, dphi)
    >>> vec_by_R = _apply_rotation_to_vector_eval(R, vec)
    >>> vec_by_q = [complex(v).real for v in rotate(vec, axis, dphi)]
    >>> max(abs(a - b) for a, b in zip(vec_by_R, vec_by_q)) < 1e-11
    True
    """


def verify_quaternion_interface_parity_randomized() -> None:
    r"""
    Stress-test matrix/quaternion agreement across deterministic random cases.

    Doctests:
    >>> from nrpy.equations.quaternion_rotations.tensor_rotation import rotate
    >>> import random
    >>> rng = random.Random(20260303)
    >>> def rand_axis():
    ...     while True:
    ...         axis = [rng.uniform(-1.0, 1.0) for _ in range(3)]
    ...         if axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2] > 1e-24:
    ...             return axis
    >>> angles = [0.0, 1e-14, -1e-14, sp.pi - sp.sympify("1e-12"), -sp.pi + sp.sympify("1e-12"), 0.7, -1.1]
    >>> max_vec_err = 0.0
    >>> max_ten_err = 0.0
    >>> for idx in range(64):
    ...     axis = rand_axis()
    ...     dphi = angles[idx % len(angles)] if idx < len(angles) else rng.uniform(-3.141592653589793, 3.141592653589793)
    ...     vec = [rng.uniform(-2.0, 2.0) for _ in range(3)]
    ...     ten = [[rng.uniform(-2.0, 2.0) for _ in range(3)] for _ in range(3)]
    ...     # Symmetrize to match BSSN tensor usage.
    ...     ten = [[0.5 * (ten[i][j] + ten[j][i]) for j in range(3)] for i in range(3)]
    ...     R = _rodrigues_matrix_eval(axis, dphi)
    ...     vec_by_R = _apply_rotation_to_vector_eval(R, vec)
    ...     vec_by_q = [complex(v).real for v in rotate(vec, axis, dphi)]
    ...     ten_by_R = _apply_rotation_to_tensor_eval(R, ten)
    ...     ten_by_q = [[complex(v).real for v in row] for row in rotate(ten, axis, dphi)]
    ...     max_vec_err = max(max_vec_err, max(abs(a - b) for a, b in zip(vec_by_R, vec_by_q)))
    ...     max_ten_err = max(max_ten_err, max(abs(ten_by_R[i][j] - ten_by_q[i][j]) for i in range(3) for j in range(3)))
    >>> max_vec_err < 2e-11
    True
    >>> max_ten_err < 5e-11
    True
    """


def register_CFunction_rotate_BSSN_Cartesian_basis_by_R() -> None:
    """
    Register C function ``rotate_BSSN_Cartesian_basis_by_R``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> from nrpy.infrastructures.BHaH.rotation.SO3_matrix_ops import assert_SO3_convention_in_text
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_rotate_BSSN_Cartesian_basis_by_R()
    >>> generated_str = cfc.CFunction_dict["rotate_BSSN_Cartesian_basis_by_R"].full_function
    >>> assert_SO3_convention_in_text(generated_str, "rotate_BSSN_Cartesian_basis_by_R")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "rotate_BSSN_Cartesian_basis_by_R_openmp", file_ext="c")
    """
    if "SO3_apply_R_to_vector" not in cfc.CFunction_dict:
        register_CFunction_SO3_apply_R_to_vector()
    if "SO3_apply_R_to_tensorDD" not in cfc.CFunction_dict:
        register_CFunction_SO3_apply_R_to_tensorDD()

    # Step 1: Basic C function metadata.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""
@brief Rotate Cartesian-basis BSSN vectors and symmetric tensors using a matrix.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout statement for helper calls: R[i][j] is row i, column j.
- DeltaR_dst_from_src maps source rotating-basis components to destination
  rotating-basis components.
- Vector update in Einstein notation:
  v^i_dst = (DeltaR)^i{}_j v^j_src.
- Covariant rank-2 update in Einstein notation:
  T^dst_{ij} = (DeltaR)_i{}^k (DeltaR)_j{}^l T^src_{kl}.
- Symmetric tensor storage contract:
  only upper-triangular components (i <= j) are authoritative on input
  and are overwritten on output; lower-triangular entries are untouched.
- Internally, full symmetric local tensors are reconstructed from upper
  components before calling SO3_apply_R_to_tensorDD().

@param[in,out] vetU BSSN rescaled shift vector, rotated in place.
@param[in,out] betU BSSN rescaled driver vector, rotated in place.
@param[in,out] lambdaU BSSN conformal connection vector, rotated in place.
@param[in,out] hDD BSSN conformal metric perturbation tensor, rotated in place.
@param[in,out] aDD BSSN trace-free extrinsic curvature tensor, rotated in place.
@param[in] DeltaR_dst_from_src Relative basis rotation matrix.
"""
    cfunc_type = "void"
    name = "rotate_BSSN_Cartesian_basis_by_R"
    params = (
        "REAL vetU[3], REAL betU[3], REAL lambdaU[3], "
        "REAL hDD[3][3], REAL aDD[3][3], const REAL DeltaR_dst_from_src[3][3]"
    )
    body = r"""
  REAL vetU_out[3], betU_out[3], lambdaU_out[3];
  REAL hDD_sym[3][3], aDD_sym[3][3];
  REAL hDD_out[3][3], aDD_out[3][3];

  // Reconstruct full symmetric tensors from authoritative upper-triangular input.
  hDD_sym[0][0] = hDD[0][0];
  hDD_sym[0][1] = hDD[0][1];
  hDD_sym[0][2] = hDD[0][2];
  hDD_sym[1][0] = hDD[0][1];
  hDD_sym[1][1] = hDD[1][1];
  hDD_sym[1][2] = hDD[1][2];
  hDD_sym[2][0] = hDD[0][2];
  hDD_sym[2][1] = hDD[1][2];
  hDD_sym[2][2] = hDD[2][2];

  aDD_sym[0][0] = aDD[0][0];
  aDD_sym[0][1] = aDD[0][1];
  aDD_sym[0][2] = aDD[0][2];
  aDD_sym[1][0] = aDD[0][1];
  aDD_sym[1][1] = aDD[1][1];
  aDD_sym[1][2] = aDD[1][2];
  aDD_sym[2][0] = aDD[0][2];
  aDD_sym[2][1] = aDD[1][2];
  aDD_sym[2][2] = aDD[2][2];

  SO3_apply_R_to_vector(DeltaR_dst_from_src, vetU, vetU_out);
  SO3_apply_R_to_vector(DeltaR_dst_from_src, betU, betU_out);
  SO3_apply_R_to_vector(DeltaR_dst_from_src, lambdaU, lambdaU_out);
  SO3_apply_R_to_tensorDD(DeltaR_dst_from_src, hDD_sym, hDD_out);
  SO3_apply_R_to_tensorDD(DeltaR_dst_from_src, aDD_sym, aDD_out);

  // Write back vectors completely.
  vetU[0] = vetU_out[0];
  vetU[1] = vetU_out[1];
  vetU[2] = vetU_out[2];
  betU[0] = betU_out[0];
  betU[1] = betU_out[1];
  betU[2] = betU_out[2];
  lambdaU[0] = lambdaU_out[0];
  lambdaU[1] = lambdaU_out[1];
  lambdaU[2] = lambdaU_out[2];

  // Write back only upper-triangular tensor components.
  hDD[0][0] = hDD_out[0][0];
  hDD[0][1] = hDD_out[0][1];
  hDD[0][2] = hDD_out[0][2];
  hDD[1][1] = hDD_out[1][1];
  hDD[1][2] = hDD_out[1][2];
  hDD[2][2] = hDD_out[2][2];

  aDD[0][0] = aDD_out[0][0];
  aDD[0][1] = aDD_out[0][1];
  aDD[0][2] = aDD_out[0][2];
  aDD[1][1] = aDD_out[1][1];
  aDD[1][2] = aDD_out[1][2];
  aDD[2][2] = aDD_out[2][2];
"""
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


def register_CFunction_rotate_BSSN_Cartesian_basis() -> None:
    """
    Register C function ``rotate_BSSN_Cartesian_basis``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> from nrpy.infrastructures.BHaH.rotation.SO3_matrix_ops import assert_SO3_convention_in_text
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_rotate_BSSN_Cartesian_basis()
    >>> generated_str = cfc.CFunction_dict["rotate_BSSN_Cartesian_basis"].full_function
    >>> assert_SO3_convention_in_text(generated_str, "rotate_BSSN_Cartesian_basis")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "rotate_BSSN_Cartesian_basis_openmp", file_ext="c")
    """
    # Step 1: Ensure matrix-based helper is registered first.
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
    nU_unit_sym = [sp.Symbol("nx"), sp.Symbol("ny"), sp.Symbol("nz")]
    R_expr = SO3Expressions.rodrigues_matrix_from_unit_axis(
        nU_unit_sym, sp.Symbol("dphi")
    )
    R_exprs, R_lhses = _rank2_exprs_lhses(R_expr, "R")
    R_codegen = _codegen_no_braces(R_exprs, R_lhses)

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
