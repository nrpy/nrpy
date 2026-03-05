"""
Generate ``unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame`` C code.

This module emits a helper that converts the cumulative rotating basis stored
in ``commondata`` into axis-angle outputs for the two-slot unrotation API:
- Slot 1 stores a recovered axis-angle pair whose Rodrigues matrix equals the
  cumulative rotation matrix built from ``xhat``, ``yhat``, ``zhat``.
- Slot 2 is deterministically set to identity ``((1,0,0), 0)``.

Rotation convention used throughout:
- ``R`` maps rotating-frame components to fixed-frame components.
- ``R[:,0]=xhat``, ``R[:,1]=yhat``, ``R[:,2]=zhat`` in the fixed basis.
- ``v_fixed = R v_rot`` and ``v_rot = R^T v_fixed``.
- ``T_fixed = R T_rot R^T`` and ``T_rot = R^T T_fixed R``.
- ``DeltaR_dst_from_src = R_dst^T R_src``.
- In generated C arrays, ``R[i][j]`` means row ``i``, column ``j``.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg

# Step 1: Import core NRPy modules needed for C code generation.
import nrpy.c_function as cfc
import nrpy.params as par
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


def register_CFunction_unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame() -> (
    None
):
    r"""
    Register the two-slot axis-angle recovery helper as generated C source.

    Generated routine contract:
    - Rebuild ``R`` from cumulative hats in ``commondata``.
    - Validate SO(3) invariants on incoming hats before conversion.
    - Recover a robust axis-angle pair for slot 1, including explicit
      near-zero and near-pi handling with deterministic sign choices.
    - Reconstruct ``R_check(n, phi)`` and fail fast if it does not match ``R``.
    - Set slot 2 to deterministic identity output.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> from nrpy.infrastructures.BHaH.rotation.register_all import register_CFunctions
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunctions()
    >>> generated_str = cfc.CFunction_dict["unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame"].full_function
    >>> "Rcheck[3][3]" in generated_str and "max_abs_err" in generated_str
    True
    >>> "nU_part2[0] = 1.0;" in generated_str and "*dphi_part2 = 0.0;" in generated_str
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
@brief Compute axis-angle outputs from cumulative hats.

This routine computes a single unrotation from the cumulative basis matrix and converts it to
one axis-angle pair in (@p nU_part1, @p dphi_part1). The second pair is set to
identity: (@p nU_part2=(1,0,0), @p dphi_part2=0).

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout statement for helper calls: R[i][j] is row i, column j.
- Axis-angle pair is recovered directly from this R with robust near-pi handling.
- Einstein notation for the matrix map represented by part 1:
  v^i_fixed = R^i{}_j v^j_rot.

@param[in] commondata Commondata structure containing cumulative regrid basis vectors.
@param[out] nU_part1 Rotation axis for first rotation.
@param[out] dphi_part1 Rotation angle for first rotation.
@param[out] nU_part2 Rotation axis for second rotation.
@param[out] dphi_part2 Rotation angle for second rotation.
"""
    cfunc_type = "void"
    name = "unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame"
    params = (
        "const commondata_struct commondata, REAL nU_part1[3], "
        "REAL *restrict dphi_part1, REAL nU_part2[3], REAL *restrict dphi_part2"
    )

    # Step 3: Build equation-driven snippets.
    xhat_sym: List[sp.Expr] = [sp.Symbol(f"xhatU[{i}]") for i in range(3)]
    yhat_sym: List[sp.Expr] = [sp.Symbol(f"yhatU[{i}]") for i in range(3)]
    zhat_sym: List[sp.Expr] = [sp.Symbol(f"zhatU[{i}]") for i in range(3)]
    hat_invariants = SO3Expressions.hat_validation_invariants(
        xhat_sym, yhat_sym, zhat_sym
    )
    R_from_hats = SO3Expressions.build_rotation_matrix_from_hats(
        xhat_sym, yhat_sym, zhat_sym
    )

    R_sym: List[List[sp.Expr]] = [
        [sp.Symbol(f"R[{i}][{j}]") for j in range(3)] for i in range(3)
    ]
    trace_expr = SO3Expressions.matrix_trace3(R_sym)
    denom_sym = sp.Symbol("denom")
    axis_general_expr = SO3Expressions.axis_angle_general_branch_axis(R_sym, denom_sym)
    pi_terms_expr = SO3Expressions.axis_angle_pi_branch_diagonal_terms(R_sym)

    n_sym: List[sp.Expr] = [sp.Symbol(f"nU_part1[{i}]") for i in range(3)]
    dphi_local_sym = sp.Symbol("dphi_local")
    Rcheck_expr = SO3Expressions.rodrigues_matrix_from_unit_axis(n_sym, dphi_local_sym)

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
        "detR_hats",
    ]
    invariant_codegen = ccg.c_codegen(
        invariant_exprs,
        invariant_lhses,
        include_braces=False,
        verbose=False,
    )

    R_exprs, R_lhses = _rank2_exprs_lhses(R_from_hats, "R")
    R_codegen = ccg.c_codegen(
        R_exprs,
        R_lhses,
        include_braces=False,
        verbose=False,
    )

    trace_codegen = ccg.c_codegen(
        [trace_expr],
        ["traceR"],
        include_braces=False,
        verbose=False,
    )

    axis_general_codegen = ccg.c_codegen(
        axis_general_expr,
        ["nU_part1[0]", "nU_part1[1]", "nU_part1[2]"],
        include_braces=False,
        verbose=False,
    )

    pi_terms_codegen = ccg.c_codegen(
        pi_terms_expr,
        ["n0", "n1", "n2"],
        include_braces=False,
        verbose=False,
    )

    Rcheck_exprs, Rcheck_lhses = _rank2_exprs_lhses(Rcheck_expr, "Rcheck")
    Rcheck_codegen = ccg.c_codegen(
        Rcheck_exprs,
        Rcheck_lhses,
        include_braces=False,
        verbose=False,
    )

    # Step 3: Validate hats, build R, convert to axis-angle, then set identity slot.
    body = (
        r"""
  const REAL ortho_tol = 1e-12;
  const REAL norm_tol = 1e-12;
  const REAL det_tol = 0.999999999999;

  REAL xhatU[3] = {commondata.cumulative_regrid_xhatU[0], commondata.cumulative_regrid_xhatU[1], commondata.cumulative_regrid_xhatU[2]};
  REAL yhatU[3] = {commondata.cumulative_regrid_yhatU[0], commondata.cumulative_regrid_yhatU[1], commondata.cumulative_regrid_yhatU[2]};
  REAL zhatU[3] = {commondata.cumulative_regrid_zhatU[0], commondata.cumulative_regrid_zhatU[1], commondata.cumulative_regrid_zhatU[2]};
  REAL x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR_hats;
"""
        + invariant_codegen
        + r"""

  const int hats_invalid =
      (fabs(x_dot_y) > ortho_tol || fabs(x_dot_z) > ortho_tol || fabs(y_dot_z) > ortho_tol ||
       fabs(x_norm - 1.0) > norm_tol || fabs(y_norm - 1.0) > norm_tol || fabs(z_norm - 1.0) > norm_tol || detR_hats < det_tol);

  if (hats_invalid) {
    fprintf(stderr,
            "ERROR in %s: invalid hats detected at consuming boundary. "
            "dots=(%.17e, %.17e, %.17e) norms=(%.17e, %.17e, %.17e) det=%.17e\n",
            __func__, x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR_hats);
    exit(1);
  }

  REAL R[3][3];
"""
        + R_codegen
        + r"""

  // Robust matrix-to-axis-angle conversion.
  const REAL pi = acos(-1.0);
  REAL traceR;
"""
        + trace_codegen
        + r"""
  REAL c = 0.5 * (traceR - 1.0);
  if (c > 1.0)
    c = 1.0;
  if (c < -1.0)
    c = -1.0;

  REAL phi = acos(c);
  *dphi_part1 = phi;

  if (phi < 1e-12) {
    nU_part1[0] = 1.0;
    nU_part1[1] = 0.0;
    nU_part1[2] = 0.0;
    *dphi_part1 = 0.0;
  } else if (fabs(pi - phi) >= 1e-8) {
    const REAL sinphi = sin(phi);
    const REAL denom = 2.0 * sinphi;
    if (fabs(denom) < 1e-15) {
      fprintf(stderr,
              "ERROR in %s: unstable general branch (|2 sin(phi)| too small). phi=%.17e denom=%.17e\n",
              __func__, phi, denom);
      exit(1);
    }

"""
        + axis_general_codegen
        + r"""
  } else {
    REAL n0, n1, n2;
"""
        + pi_terms_codegen
        + r"""
    if (n0 < 0.0)
      n0 = 0.0;
    if (n1 < 0.0)
      n1 = 0.0;
    if (n2 < 0.0)
      n2 = 0.0;

    n0 = sqrt(n0);
    n1 = sqrt(n1);
    n2 = sqrt(n2);

    int imax = 0;
    if (n1 > n0 && n1 >= n2)
      imax = 1;
    else if (n2 > n0 && n2 > n1)
      imax = 2;

    if (imax == 0) {
      if (n1 > 0.0)
        n1 = copysign(n1, R[0][1] + R[1][0]);
      if (n2 > 0.0)
        n2 = copysign(n2, R[0][2] + R[2][0]);
    } else if (imax == 1) {
      if (n0 > 0.0)
        n0 = copysign(n0, R[0][1] + R[1][0]);
      if (n2 > 0.0)
        n2 = copysign(n2, R[1][2] + R[2][1]);
    } else {
      if (n0 > 0.0)
        n0 = copysign(n0, R[0][2] + R[2][0]);
      if (n1 > 0.0)
        n1 = copysign(n1, R[1][2] + R[2][1]);
    }

    nU_part1[0] = n0;
    nU_part1[1] = n1;
    nU_part1[2] = n2;

    REAL nnorm = sqrt(nU_part1[0] * nU_part1[0] + nU_part1[1] * nU_part1[1] + nU_part1[2] * nU_part1[2]);
    if (nnorm < 1e-15) {
      // Deterministic fallback in fully degenerate pi-branch.
      if (R[0][0] >= R[1][1] && R[0][0] >= R[2][2]) {
        nU_part1[0] = 1.0;
        nU_part1[1] = 0.0;
        nU_part1[2] = 0.0;
      } else if (R[1][1] >= R[0][0] && R[1][1] >= R[2][2]) {
        nU_part1[0] = 0.0;
        nU_part1[1] = 1.0;
        nU_part1[2] = 0.0;
      } else {
        nU_part1[0] = 0.0;
        nU_part1[1] = 0.0;
        nU_part1[2] = 1.0;
      }
      nnorm = 1.0;
    }
    nU_part1[0] /= nnorm;
    nU_part1[1] /= nnorm;
    nU_part1[2] /= nnorm;

    // Deterministic sign convention for pi-branch axis ambiguity.
    if (nU_part1[0] < 0.0 || (fabs(nU_part1[0]) <= 1e-16 && nU_part1[1] < 0.0) ||
        (fabs(nU_part1[0]) <= 1e-16 && fabs(nU_part1[1]) <= 1e-16 && nU_part1[2] < 0.0)) {
      nU_part1[0] = -nU_part1[0];
      nU_part1[1] = -nU_part1[1];
      nU_part1[2] = -nU_part1[2];
    }
    *dphi_part1 = pi;
  }

  REAL nnorm = sqrt(nU_part1[0] * nU_part1[0] + nU_part1[1] * nU_part1[1] + nU_part1[2] * nU_part1[2]);
  if (nnorm < 1e-15) {
    fprintf(stderr, "ERROR in %s: could not recover nonzero rotation axis.\n", __func__);
    exit(1);
  }
  nU_part1[0] /= nnorm;
  nU_part1[1] /= nnorm;
  nU_part1[2] /= nnorm;

  // Debug reconstruction check: R_check(n,phi) should match input R.
  const REAL dphi_local = *dphi_part1;

  REAL Rcheck[3][3];
"""
        + Rcheck_codegen
        + r"""

  REAL max_abs_err = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      const REAL err = fabs(Rcheck[i][j] - R[i][j]);
      if (err > max_abs_err)
        max_abs_err = err;
    }
  }

  if (max_abs_err > 1e-10) {
    fprintf(stderr,
            "ERROR in %s: matrix->axis-angle reconstruction mismatch. max_abs_err=%.17e phi=%.17e "
            "n=(%.17e, %.17e, %.17e)\n",
            __func__, max_abs_err, *dphi_part1, nU_part1[0], nU_part1[1], nU_part1[2]);
    exit(1);
  }

  // Second output slot is deterministic identity.
  nU_part2[0] = 1.0;
  nU_part2[1] = 0.0;
  nU_part2[2] = 0.0;
  *dphi_part2 = 0.0;
"""
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
