"""
Register C helper for matrix-first unrotation of Cartesian vectors.

This hot path uses cumulative hats -> R directly and avoids
axis-angle two-stage solve machinery.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
# Step 1: Import core NRPy modules needed for C code generation.
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.equations.rotation.SO3_rotations import SO3Expressions

SO3_CONVENTION_REQUIRED_SNIPPETS = (
    "R maps rotating-frame components -> fixed-frame components.",
    "R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat",
    "v_fixed = R v_rot, v_rot = R^T v_fixed.",
    "T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.",
    "DeltaR_dst_from_src = R_dst^T R_src.",
    "R[i][j]",
    "row i, column j.",
)


def assert_SO3_convention_in_text(text: str, context: str) -> None:
    """
    Assert that a text block includes all locked SO(3) convention snippets.

    :param text: Text block that must contain the locked convention snippets.
    :param context: Context label used in the assertion message.
    :raises AssertionError: If one or more required snippets are missing.
    """
    missing = [
        snippet for snippet in SO3_CONVENTION_REQUIRED_SNIPPETS if snippet not in text
    ]
    if missing:
        raise AssertionError(
            f"{context} is missing locked SO(3) convention snippets:\n- "
            + "\n- ".join(missing)
        )


def _rank1_exprs_lhses(rank1: List[sp.Expr], varname: str) -> tuple[List[sp.Expr], List[str]]:
    """
    Flatten a rank-1 vector into expression and LHS-name lists.

    :param rank1: Rank-1 vector.
    :param varname: Output C variable name.
    :return: Tuple of flattened expression and LHS-name lists.
    """
    exprs: List[sp.Expr] = []
    lhses: List[str] = []
    for i in range(3):
        exprs.append(rank1[i])
        lhses.append(f"{varname}[{i}]")
    return exprs, lhses


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


def verify_unrotate_hot_path_antiparallel_regression() -> None:
    r"""
    Check the direct-matrix hot path on a deterministic anti-parallel basis case.

    This verifies numerical behavior for a 180-degree basis flip without any
    axis-angle decomposition in the hot path. With the locked SO(3) convention,
    the operation tested here is
    ``x^i_fixed = R^i{}_j x^j_rot``.

    Doctests:
    >>> from nrpy.equations.quaternion_rotations.tensor_rotation import rotate
    >>> xhat = [-1.0, 0.0, 0.0]
    >>> yhat = [0.0, 1.0, 0.0]
    >>> zhat = [0.0, 0.0, -1.0]
    >>> R = [[xhat[0], yhat[0], zhat[0]], [xhat[1], yhat[1], zhat[1]], [xhat[2], yhat[2], zhat[2]]]
    >>> v_rot = [1.25, -0.5, 2.0]
    >>> v_fixed = [sum(R[i][j] * v_rot[j] for j in range(3)) for i in range(3)]
    >>> v_fixed
    [-1.25, -0.5, -2.0]
    >>> v_expected = [complex(v).real for v in rotate(v_rot, [0.0, 1.0, 0.0], 3.141592653589793)]
    >>> max(abs(a - b) for a, b in zip(v_fixed, v_expected)) < 1e-14
    True
    """


def register_CFunction_unrotate_xCart_to_fixed_frame() -> None:
    r"""
    Register C function ``unrotate_xCart_to_fixed_frame``.

    This wrapper consumes cumulative hats from ``commondata``, constructs
    ``R^i{}_j``, and applies
    ``x^i_fixed = R^i{}_j x^j_rot`` in place.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> from nrpy.infrastructures.BHaH.rotation.register_all import register_CFunctions
    >>> from nrpy.infrastructures.BHaH.rotation.unrotate_xCart_to_fixed_frame import assert_SO3_convention_in_text
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunctions()
    >>> generated_str = cfc.CFunction_dict["unrotate_xCart_to_fixed_frame"].full_function
    >>> assert_SO3_convention_in_text(generated_str, "unrotate_xCart_to_fixed_frame")
    >>> "unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame(" in generated_str
    False
    >>> "build_R_from_cumulative_hats(" in generated_str or "SO3_apply_R_to_vector(" in generated_str
    False
    >>> "const REAL x_in[3]" in generated_str and "xCart[0]" in generated_str and "detR" in generated_str
    True
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "openmp", file_ext="c")
    """
    # Step 1: Register cumulative-hat commondata CodeParameters.
    _ = par.register_CodeParameter(
        "REAL[3]",
        __name__,
        "cumulative_regrid_xhatU",
        1e300,
        commondata=True,
        add_to_set_CodeParameters_h=False,
        description="Cumulative rotating-frame xhat expressed in fixed frame.",
    )
    _ = par.register_CodeParameter(
        "REAL[3]",
        __name__,
        "cumulative_regrid_yhatU",
        1e300,
        commondata=True,
        add_to_set_CodeParameters_h=False,
        description="Cumulative rotating-frame yhat expressed in fixed frame.",
    )
    _ = par.register_CodeParameter(
        "REAL[3]",
        __name__,
        "cumulative_regrid_zhatU",
        1e300,
        commondata=True,
        add_to_set_CodeParameters_h=False,
        description="Cumulative rotating-frame zhat expressed in fixed frame.",
    )

    # Step 2: Set C function metadata.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""
@brief Map a Cartesian vector from rotating-frame coordinates to fixed-frame coordinates.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout statement for helper calls: R[i][j] is row i, column j.
- Einstein notation for this routine:
  x^i_fixed = R^i{}_j x^j_rot.
- This routine applies x_fixed = R * x_rot directly, i.e. it does not
  pass through an axis-angle decomposition.

@param[in] commondata Commondata structure with cumulative regrid basis vectors.
@param[in,out] xCart Cartesian point/vector in rotating basis on input,
                    overwritten with fixed-basis components on output.
"""
    cfunc_type = "void"
    name = "unrotate_xCart_to_fixed_frame"
    params = "const commondata_struct commondata, REAL xCart[3]"

    # Step 3: Build equation-driven snippets for invariants, R construction, and vector map.
    xhat_sym = ixp.declarerank1("xhatU")
    yhat_sym = ixp.declarerank1("yhatU")
    zhat_sym = ixp.declarerank1("zhatU")
    hat_invariants = SO3Expressions.hat_validation_invariants(xhat_sym, yhat_sym, zhat_sym)
    R_from_hats = SO3Expressions.build_rotation_matrix_from_hats(xhat_sym, yhat_sym, zhat_sym)

    R_sym = ixp.declarerank2("R", symmetry="nosym")
    x_in_sym = ixp.declarerank1("x_in")
    x_fixed_expr = SO3Expressions.apply_R_to_vector(R_sym, x_in_sym)

    sub_hat: dict[sp.Symbol, sp.Expr] = {}
    for i in range(3):
        sub_hat[xhat_sym[i]] = sp.Symbol(f"xhatU[{i}]")
        sub_hat[yhat_sym[i]] = sp.Symbol(f"yhatU[{i}]")
        sub_hat[zhat_sym[i]] = sp.Symbol(f"zhatU[{i}]")
        sub_hat[x_in_sym[i]] = sp.Symbol(f"x_in[{i}]")
        for j in range(3):
            sub_hat[R_sym[i][j]] = sp.Symbol(f"R[{i}][{j}]")

    invariant_exprs = [
        hat_invariants["x_dot_y"],
        hat_invariants["x_dot_z"],
        hat_invariants["y_dot_z"],
        hat_invariants["x_norm"],
        hat_invariants["y_norm"],
        hat_invariants["z_norm"],
        hat_invariants["detR"],
    ]
    invariant_lhses = [
        "x_dot_y",
        "x_dot_z",
        "y_dot_z",
        "x_norm",
        "y_norm",
        "z_norm",
        "detR",
    ]
    invariant_codegen = ccg.c_codegen(
        [expr.subs(sub_hat) for expr in invariant_exprs],
        invariant_lhses,
        include_braces=False,
        verbose=False,
    )

    R_exprs, R_lhses = _rank2_exprs_lhses(R_from_hats, "R")
    R_codegen = ccg.c_codegen(
        [expr.subs(sub_hat) for expr in R_exprs],
        R_lhses,
        include_braces=False,
        verbose=False,
    )

    x_exprs, x_lhses = _rank1_exprs_lhses(x_fixed_expr, "xCart")
    x_codegen = ccg.c_codegen(
        [expr.subs(sub_hat) for expr in x_exprs],
        x_lhses,
        include_braces=False,
        verbose=False,
    )

    # Step 3: Validate hats, build R from hats, and apply x_fixed = R x_rot.
    body = (
        r"""
  const REAL ortho_tol = 1e-12;
  const REAL norm_tol = 1e-12;
  const REAL det_tol = 0.999999999999;

  REAL xhatU[3] = {commondata.cumulative_regrid_xhatU[0], commondata.cumulative_regrid_xhatU[1], commondata.cumulative_regrid_xhatU[2]};
  REAL yhatU[3] = {commondata.cumulative_regrid_yhatU[0], commondata.cumulative_regrid_yhatU[1], commondata.cumulative_regrid_yhatU[2]};
  REAL zhatU[3] = {commondata.cumulative_regrid_zhatU[0], commondata.cumulative_regrid_zhatU[1], commondata.cumulative_regrid_zhatU[2]};
  REAL x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR;
"""
        + invariant_codegen
        + r"""

  const int hats_invalid =
      (fabs(x_dot_y) > ortho_tol || fabs(x_dot_z) > ortho_tol || fabs(y_dot_z) > ortho_tol ||
       fabs(x_norm - 1.0) > norm_tol || fabs(y_norm - 1.0) > norm_tol || fabs(z_norm - 1.0) > norm_tol || detR < det_tol);

  if (hats_invalid) {
    fprintf(stderr,
            "ERROR in %s: invalid hats detected at consuming boundary. "
            "dots=(%.17e, %.17e, %.17e) norms=(%.17e, %.17e, %.17e) det=%.17e\n",
            __func__, x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR);
    exit(1);
  }

  REAL R[3][3];
"""
        + R_codegen
        + r"""

  // Alias-safe copy: permits xCart output to overwrite xCart input.
  const REAL x_in[3] = {xCart[0], xCart[1], xCart[2]};
"""
        + x_codegen
    )
    # Step 4: Register the C function.
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

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
