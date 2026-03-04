"""
Register SO(3)-matrix-first rotation helpers for BHaH.

Global convention used throughout this module:
- R maps rotating-frame components -> fixed-frame components.
- Columns are rotating-frame basis vectors expressed in fixed coordinates:
  R[:,0] = xhat, R[:,1] = yhat, R[:,2] = zhat.
- v_fixed = R v_rot, and v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, and T_rot = R^T T_fixed R.
- Relative transform from src rotating basis to dst rotating basis:
  DeltaR_dst_from_src = R_dst^T R_src.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import re
from pathlib import Path

import nrpy.c_function as cfc
import nrpy.params as par

SO3_CONVENTION_REQUIRED_SNIPPETS = (
    "R maps rotating-frame components -> fixed-frame components.",
    "R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat",
    "v_fixed = R v_rot, v_rot = R^T v_fixed.",
    "T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.",
    "DeltaR_dst_from_src = R_dst^T R_src.",
    "R[i][j]",
    "row i, column j.",
)


def assert_so3_convention_in_text(text: str, context: str) -> None:
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


def audit_cumulative_hat_usage_in_nrpy() -> dict[str, list[str]]:
    r"""
    Audit cumulative hat read/write usage in active ``nrpy/`` Python sources.

    :return: Mapping with two keys:
        - ``reads``: files containing commondata reads of cumulative hats.
        - ``writes``: files containing direct assignments to cumulative hats.

    Doctests:
    >>> usage = audit_cumulative_hat_usage_in_nrpy()
    >>> len(usage["writes"]) == 0
    True
    >>> any(path.endswith("unrotate_xCart_to_fixed_frame.py") for path in usage["reads"])
    True
    >>> any(path.endswith("unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame.py") for path in usage["reads"])
    True
    """
    nrpy_root = Path(__file__).resolve().parents[3]
    read_pattern = re.compile(r"commondata\.cumulative_regrid_[xyz]hatU\[[0-2]\]")
    write_pattern = re.compile(
        r"(?:commondata\.)?cumulative_regrid_[xyz]hatU\[[0-2]\]\s*="
    )

    reads: list[str] = []
    writes: list[str] = []
    for py_file in nrpy_root.rglob("*.py"):
        text = py_file.read_text(encoding="utf-8")
        if read_pattern.search(text):
            reads.append(str(py_file.relative_to(nrpy_root)))
        if write_pattern.search(text):
            writes.append(str(py_file.relative_to(nrpy_root)))
    return {"reads": sorted(reads), "writes": sorted(writes)}


def verify_so3_convention_lockstep_in_registered_helpers() -> None:
    r"""
    Verify locked convention text appears in all registered SO(3) helper docs.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunctions_so3_matrix_ops()
    >>> names = [
    ...     "build_R_from_cumulative_hats",
    ...     "so3_validate_and_optionally_fix_hats",
    ...     "so3_apply_R_to_vector",
    ...     "so3_apply_RT_to_vector",
    ...     "so3_apply_R_to_tensorDD",
    ...     "so3_apply_RT_to_tensorDD",
    ...     "so3_relative_rotation_dst_from_src",
    ...     "so3_matrix_to_axis_angle",
    ... ]
    >>> for name in names:
    ...     assert_so3_convention_in_text(cfc.CFunction_dict[name].full_function, name)
    """


def register_rotation_commondata_CodeParameters() -> None:
    """
    Register rotation-related commondata CodeParameters.

    Integration contract for future cumulative-hat update sites:
    - Immediately after each authoritative update, call
      ``so3_validate_and_optionally_fix_hats(..., do_fix=1)``.
    - Active read-side consumers use ``do_fix=0`` and hard-fail on invalid hats.
    """
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


def register_CFunction_build_R_from_cumulative_hats() -> None:
    r"""
    Register C function ``build_R_from_cumulative_hats``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_build_R_from_cumulative_hats()
    >>> generated_str = cfc.CFunction_dict["build_R_from_cumulative_hats"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "build_R_from_cumulative_hats_openmp", file_ext="c")
    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Build SO(3) rotation matrix from cumulative rotating-frame hats.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- Columns are hats in fixed basis: R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat.
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout is row-major indexing: R[i][j] is row i, column j.

@param[in] xhatU Rotating-frame xhat in fixed coordinates.
@param[in] yhatU Rotating-frame yhat in fixed coordinates.
@param[in] zhatU Rotating-frame zhat in fixed coordinates.
@param[out] R Rotation matrix mapping rotating -> fixed components.
"""
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name="build_R_from_cumulative_hats",
        params="const REAL xhatU[3], const REAL yhatU[3], const REAL zhatU[3], REAL R[3][3]",
        body=r"""
  // R[i][j] means row i, column j. Fill R by columns from hats.
  R[0][0] = xhatU[0];
  R[1][0] = xhatU[1];
  R[2][0] = xhatU[2];

  R[0][1] = yhatU[0];
  R[1][1] = yhatU[1];
  R[2][1] = yhatU[2];

  R[0][2] = zhatU[0];
  R[1][2] = zhatU[1];
  R[2][2] = zhatU[2];
""",
        include_CodeParameters_h=False,
    )


def register_CFunction_so3_validate_and_optionally_fix_hats() -> None:
    r"""
    Register C function ``so3_validate_and_optionally_fix_hats``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_validate_and_optionally_fix_hats()
    >>> generated_str = cfc.CFunction_dict["so3_validate_and_optionally_fix_hats"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "so3_validate_and_optionally_fix_hats_openmp", file_ext="c")
    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Validate cumulative hats as a right-handed orthonormal basis; optionally repair.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout statement for helper calls: R[i][j] is row i, column j.

Validation triggers if any of:
- |x·y|, |x·z|, |y·z| > 1e-12
- |||x||-1|, |||y||-1|, |||z||-1| > 1e-12
- det(R) < 0.999999999999, where columns of R are (xhat,yhat,zhat)

Deterministic repair policy when @p do_fix != 0:
1) x <- normalize(x)
2) y <- y - (x·y)x; y <- normalize(y)
3) z <- x × y; z <- normalize(z)
4) y <- z × x; y <- normalize(y)
5) if det(R) < 0 then z <- -z and y <- z × x; y <- normalize(y)

If post-fix invariants still fail, this routine aborts with diagnostics.

@param[in,out] xhatU Candidate xhat vector (fixed basis).
@param[in,out] yhatU Candidate yhat vector (fixed basis).
@param[in,out] zhatU Candidate zhat vector (fixed basis).
@param[in] do_fix If nonzero, apply deterministic repair; otherwise report invalidity.
@return 0 if valid (possibly after repair), 1 if invalid and not repaired.
"""
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=includes,
        desc=desc,
        cfunc_type="int",
        name="so3_validate_and_optionally_fix_hats",
        params="REAL xhatU[3], REAL yhatU[3], REAL zhatU[3], const int do_fix",
        body=r"""
  const REAL ortho_tol = 1e-12;
  const REAL norm_tol = 1e-12;
  const REAL det_tol = 0.999999999999;
  const REAL tiny = 1e-300;

  REAL x_dot_y = xhatU[0] * yhatU[0] + xhatU[1] * yhatU[1] + xhatU[2] * yhatU[2];
  REAL x_dot_z = xhatU[0] * zhatU[0] + xhatU[1] * zhatU[1] + xhatU[2] * zhatU[2];
  REAL y_dot_z = yhatU[0] * zhatU[0] + yhatU[1] * zhatU[1] + yhatU[2] * zhatU[2];

  REAL x_norm = sqrt(xhatU[0] * xhatU[0] + xhatU[1] * xhatU[1] + xhatU[2] * xhatU[2]);
  REAL y_norm = sqrt(yhatU[0] * yhatU[0] + yhatU[1] * yhatU[1] + yhatU[2] * yhatU[2]);
  REAL z_norm = sqrt(zhatU[0] * zhatU[0] + zhatU[1] * zhatU[1] + zhatU[2] * zhatU[2]);

  const REAL cross_yz0 = yhatU[1] * zhatU[2] - yhatU[2] * zhatU[1];
  const REAL cross_yz1 = yhatU[2] * zhatU[0] - yhatU[0] * zhatU[2];
  const REAL cross_yz2 = yhatU[0] * zhatU[1] - yhatU[1] * zhatU[0];
  REAL detR = xhatU[0] * cross_yz0 + xhatU[1] * cross_yz1 + xhatU[2] * cross_yz2;

  const int needs_fix =
      (fabs(x_dot_y) > ortho_tol || fabs(x_dot_z) > ortho_tol || fabs(y_dot_z) > ortho_tol ||
       fabs(x_norm - 1.0) > norm_tol || fabs(y_norm - 1.0) > norm_tol || fabs(z_norm - 1.0) > norm_tol || detR < det_tol);

  if (!needs_fix)
    return 0;

  if (!do_fix) {
    fprintf(stderr,
            "ERROR in %s: invalid hats detected and do_fix=0. "
            "dots=(%.17e, %.17e, %.17e) norms=(%.17e, %.17e, %.17e) det=%.17e\n",
            __func__, x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR);
    return 1;
  }

  // Deterministic Gram-Schmidt repair and right-handedness enforcement.
  if (x_norm <= tiny) {
    fprintf(stderr, "ERROR in %s: xhat norm too small during fix (%.17e).\n", __func__, x_norm);
    exit(1);
  }
  for (int i = 0; i < 3; i++)
    xhatU[i] /= x_norm;

  x_dot_y = xhatU[0] * yhatU[0] + xhatU[1] * yhatU[1] + xhatU[2] * yhatU[2];
  for (int i = 0; i < 3; i++)
    yhatU[i] -= x_dot_y * xhatU[i];

  y_norm = sqrt(yhatU[0] * yhatU[0] + yhatU[1] * yhatU[1] + yhatU[2] * yhatU[2]);
  if (y_norm <= tiny) {
    fprintf(stderr, "ERROR in %s: yhat norm too small during fix (%.17e).\n", __func__, y_norm);
    exit(1);
  }
  for (int i = 0; i < 3; i++)
    yhatU[i] /= y_norm;

  zhatU[0] = xhatU[1] * yhatU[2] - xhatU[2] * yhatU[1];
  zhatU[1] = xhatU[2] * yhatU[0] - xhatU[0] * yhatU[2];
  zhatU[2] = xhatU[0] * yhatU[1] - xhatU[1] * yhatU[0];

  z_norm = sqrt(zhatU[0] * zhatU[0] + zhatU[1] * zhatU[1] + zhatU[2] * zhatU[2]);
  if (z_norm <= tiny) {
    fprintf(stderr, "ERROR in %s: zhat norm too small during fix (%.17e).\n", __func__, z_norm);
    exit(1);
  }
  for (int i = 0; i < 3; i++)
    zhatU[i] /= z_norm;

  yhatU[0] = zhatU[1] * xhatU[2] - zhatU[2] * xhatU[1];
  yhatU[1] = zhatU[2] * xhatU[0] - zhatU[0] * xhatU[2];
  yhatU[2] = zhatU[0] * xhatU[1] - zhatU[1] * xhatU[0];

  y_norm = sqrt(yhatU[0] * yhatU[0] + yhatU[1] * yhatU[1] + yhatU[2] * yhatU[2]);
  if (y_norm <= tiny) {
    fprintf(stderr, "ERROR in %s: yhat norm too small after z cross x (%.17e).\n", __func__, y_norm);
    exit(1);
  }
  for (int i = 0; i < 3; i++)
    yhatU[i] /= y_norm;

  const REAL post_cross_yz0 = yhatU[1] * zhatU[2] - yhatU[2] * zhatU[1];
  const REAL post_cross_yz1 = yhatU[2] * zhatU[0] - yhatU[0] * zhatU[2];
  const REAL post_cross_yz2 = yhatU[0] * zhatU[1] - yhatU[1] * zhatU[0];
  detR = xhatU[0] * post_cross_yz0 + xhatU[1] * post_cross_yz1 + xhatU[2] * post_cross_yz2;

  if (detR < 0.0) {
    for (int i = 0; i < 3; i++)
      zhatU[i] = -zhatU[i];

    yhatU[0] = zhatU[1] * xhatU[2] - zhatU[2] * xhatU[1];
    yhatU[1] = zhatU[2] * xhatU[0] - zhatU[0] * xhatU[2];
    yhatU[2] = zhatU[0] * xhatU[1] - zhatU[1] * xhatU[0];

    y_norm = sqrt(yhatU[0] * yhatU[0] + yhatU[1] * yhatU[1] + yhatU[2] * yhatU[2]);
    if (y_norm <= tiny) {
      fprintf(stderr, "ERROR in %s: yhat norm too small in right-handedness fix (%.17e).\n", __func__, y_norm);
      exit(1);
    }
    for (int i = 0; i < 3; i++)
      yhatU[i] /= y_norm;
  }

  x_dot_y = xhatU[0] * yhatU[0] + xhatU[1] * yhatU[1] + xhatU[2] * yhatU[2];
  x_dot_z = xhatU[0] * zhatU[0] + xhatU[1] * zhatU[1] + xhatU[2] * zhatU[2];
  y_dot_z = yhatU[0] * zhatU[0] + yhatU[1] * zhatU[1] + yhatU[2] * zhatU[2];
  x_norm = sqrt(xhatU[0] * xhatU[0] + xhatU[1] * xhatU[1] + xhatU[2] * xhatU[2]);
  y_norm = sqrt(yhatU[0] * yhatU[0] + yhatU[1] * yhatU[1] + yhatU[2] * yhatU[2]);
  z_norm = sqrt(zhatU[0] * zhatU[0] + zhatU[1] * zhatU[1] + zhatU[2] * zhatU[2]);

  const REAL final_cross_yz0 = yhatU[1] * zhatU[2] - yhatU[2] * zhatU[1];
  const REAL final_cross_yz1 = yhatU[2] * zhatU[0] - yhatU[0] * zhatU[2];
  const REAL final_cross_yz2 = yhatU[0] * zhatU[1] - yhatU[1] * zhatU[0];
  detR = xhatU[0] * final_cross_yz0 + xhatU[1] * final_cross_yz1 + xhatU[2] * final_cross_yz2;

  const int still_invalid =
      (fabs(x_dot_y) > ortho_tol || fabs(x_dot_z) > ortho_tol || fabs(y_dot_z) > ortho_tol ||
       fabs(x_norm - 1.0) > norm_tol || fabs(y_norm - 1.0) > norm_tol || fabs(z_norm - 1.0) > norm_tol || detR < det_tol);

  if (still_invalid) {
    fprintf(stderr,
            "ERROR in %s: post-fix hats remain invalid. "
            "dots=(%.17e, %.17e, %.17e) norms=(%.17e, %.17e, %.17e) det=%.17e\n",
            __func__, x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR);
    exit(1);
  }

  return 0;
""",
        include_CodeParameters_h=False,
    )


def register_CFunction_so3_apply_R_to_vector() -> None:
    r"""
    Register C function ``so3_apply_R_to_vector``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_apply_R_to_vector()
    >>> generated_str = cfc.CFunction_dict["so3_apply_R_to_vector"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "so3_apply_R_to_vector_openmp", file_ext="c")
    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Apply R to a vector using v_fixed = R v_rot.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout uses R[i][j] = row i, column j.

@param[in] R Rotation matrix.
@param[in] v_rot Vector components in rotating basis.
@param[out] v_fixed Vector components in fixed basis.
"""
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name="so3_apply_R_to_vector",
        params="const REAL R[3][3], const REAL v_rot[3], REAL v_fixed[3]",
        body=r"""
  // Alias-safe copy: permits v_fixed == v_rot.
  const REAL v_in[3] = {v_rot[0], v_rot[1], v_rot[2]};
  // v_fixed[i] = sum_j R[i][j] v_in[j] with R[i][j] = row i, column j.
  for (int i = 0; i < 3; i++) {
    REAL accum = 0.0;
    for (int j = 0; j < 3; j++) {
      accum += R[i][j] * v_in[j];
    }
    v_fixed[i] = accum;
  }
""",
        include_CodeParameters_h=False,
    )


def register_CFunction_so3_apply_RT_to_vector() -> None:
    r"""
    Register C function ``so3_apply_RT_to_vector``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_apply_RT_to_vector()
    >>> generated_str = cfc.CFunction_dict["so3_apply_RT_to_vector"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "so3_apply_RT_to_vector_openmp", file_ext="c")
    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Apply R^T to a vector using v_rot = R^T v_fixed.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout uses R[i][j] = row i, column j.

@param[in] R Rotation matrix.
@param[in] v_fixed Vector components in fixed basis.
@param[out] v_rot Vector components in rotating basis.
"""
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name="so3_apply_RT_to_vector",
        params="const REAL R[3][3], const REAL v_fixed[3], REAL v_rot[3]",
        body=r"""
  // Alias-safe copy: permits v_rot == v_fixed.
  const REAL v_in[3] = {v_fixed[0], v_fixed[1], v_fixed[2]};
  // v_rot[i] = sum_j (R^T)[i][j] v_in[j] = sum_j R[j][i] v_in[j].
  for (int i = 0; i < 3; i++) {
    REAL accum = 0.0;
    for (int j = 0; j < 3; j++) {
      accum += R[j][i] * v_in[j];
    }
    v_rot[i] = accum;
  }
""",
        include_CodeParameters_h=False,
    )


def register_CFunction_so3_apply_R_to_tensorDD() -> None:
    r"""
    Register C function ``so3_apply_R_to_tensorDD``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_apply_R_to_tensorDD()
    >>> generated_str = cfc.CFunction_dict["so3_apply_R_to_tensorDD"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "so3_apply_R_to_tensorDD_openmp", file_ext="c")
    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Apply R to a rank-2 tensor: T_fixed = R T_rot R^T.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout uses R[i][j] = row i, column j.

@param[in] R Rotation matrix.
@param[in] T_rot Tensor components in rotating basis.
@param[out] T_fixed Tensor components in fixed basis.
"""
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name="so3_apply_R_to_tensorDD",
        params="const REAL R[3][3], const REAL T_rot[3][3], REAL T_fixed[3][3]",
        body=r"""
  REAL tmp[3][3];

  // tmp = R * T_rot: tmp[i][j] = sum_k R[i][k] * T_rot[k][j].
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REAL accum = 0.0;
      for (int k = 0; k < 3; k++) {
        accum += R[i][k] * T_rot[k][j];
      }
      tmp[i][j] = accum;
    }
  }

  // T_fixed = tmp * R^T: T_fixed[i][j] = sum_k tmp[i][k] * R[j][k].
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REAL accum = 0.0;
      for (int k = 0; k < 3; k++) {
        accum += tmp[i][k] * R[j][k];
      }
      T_fixed[i][j] = accum;
    }
  }
""",
        include_CodeParameters_h=False,
    )


def register_CFunction_so3_apply_RT_to_tensorDD() -> None:
    r"""
    Register C function ``so3_apply_RT_to_tensorDD``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_apply_RT_to_tensorDD()
    >>> generated_str = cfc.CFunction_dict["so3_apply_RT_to_tensorDD"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "so3_apply_RT_to_tensorDD_openmp", file_ext="c")
    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Apply R^T to a rank-2 tensor: T_rot = R^T T_fixed R.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout uses R[i][j] = row i, column j.

@param[in] R Rotation matrix.
@param[in] T_fixed Tensor components in fixed basis.
@param[out] T_rot Tensor components in rotating basis.
"""
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name="so3_apply_RT_to_tensorDD",
        params="const REAL R[3][3], const REAL T_fixed[3][3], REAL T_rot[3][3]",
        body=r"""
  REAL tmp[3][3];

  // tmp = R^T * T_fixed: tmp[i][j] = sum_k R[k][i] * T_fixed[k][j].
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REAL accum = 0.0;
      for (int k = 0; k < 3; k++) {
        accum += R[k][i] * T_fixed[k][j];
      }
      tmp[i][j] = accum;
    }
  }

  // T_rot = tmp * R: T_rot[i][j] = sum_k tmp[i][k] * R[k][j].
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REAL accum = 0.0;
      for (int k = 0; k < 3; k++) {
        accum += tmp[i][k] * R[k][j];
      }
      T_rot[i][j] = accum;
    }
  }
""",
        include_CodeParameters_h=False,
    )


def register_CFunction_so3_relative_rotation_dst_from_src() -> None:
    r"""
    Register C function ``so3_relative_rotation_dst_from_src``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_relative_rotation_dst_from_src()
    >>> generated_str = cfc.CFunction_dict["so3_relative_rotation_dst_from_src"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "so3_relative_rotation_dst_from_src_openmp", file_ext="c")
    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Compute relative basis transform DeltaR_dst_from_src = R_dst^T R_src.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout statement for helper calls: R[i][j] is row i, column j.

Given two absolute rotating->fixed rotations, @p R_src and @p R_dst, this
returns the matrix that maps components in the src rotating basis to
components in the dst rotating basis:
v_dst = DeltaR_dst_from_src * v_src.

@param[in] R_dst Absolute destination rotation (rotating->fixed).
@param[in] R_src Absolute source rotation (rotating->fixed).
@param[out] DeltaR_dst_from_src Relative transform from src rotating basis to dst rotating basis.
"""
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name="so3_relative_rotation_dst_from_src",
        params=(
            "const REAL R_dst[3][3], const REAL R_src[3][3], "
            "REAL DeltaR_dst_from_src[3][3]"
        ),
        body=r"""
  // DeltaR_dst_from_src[i][j] = sum_k (R_dst^T)[i][k] R_src[k][j] = sum_k R_dst[k][i] R_src[k][j].
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REAL accum = 0.0;
      for (int k = 0; k < 3; k++) {
        accum += R_dst[k][i] * R_src[k][j];
      }
      DeltaR_dst_from_src[i][j] = accum;
    }
  }
""",
        include_CodeParameters_h=False,
    )


def register_CFunction_so3_matrix_to_axis_angle() -> None:
    r"""
    Register C function ``so3_matrix_to_axis_angle``.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_matrix_to_axis_angle()
    >>> generated_str = cfc.CFunction_dict["so3_matrix_to_axis_angle"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "so3_matrix_to_axis_angle_openmp", file_ext="c")
    """
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Convert SO(3) matrix to axis-angle with robust pi-branch handling.

Convention:
- R maps rotating-frame components -> fixed-frame components.
- R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
- v_fixed = R v_rot, v_rot = R^T v_fixed.
- T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
- DeltaR_dst_from_src = R_dst^T R_src.
- C layout statement for helper calls: R[i][j] is row i, column j.

Algorithm:
1) c = clamp((trace(R)-1)/2, -1, 1), phi = acos(c)
2) Small-angle branch: phi < 1e-12 -> n=(1,0,0), phi=0
3) General branch: |pi-phi| >= 1e-8 -> n from skew(R)/(2 sin(phi))
4) Pi branch: |pi-phi| < 1e-8 -> diagonal-based recovery plus deterministic
   sign resolution from symmetric off-diagonal terms
5) Debug check: reconstruct R(n,phi) and abort if mismatch exceeds tolerance

@param[in] R Input rotation matrix.
@param[out] nU Unit axis vector.
@param[out] dphi Rotation angle in radians.
"""
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name="so3_matrix_to_axis_angle",
        params="const REAL R[3][3], REAL nU[3], REAL *restrict dphi",
        body=r"""
  const REAL pi = acos(-1.0);
  const REAL traceR = R[0][0] + R[1][1] + R[2][2];
  REAL c = 0.5 * (traceR - 1.0);
  if (c > 1.0)
    c = 1.0;
  if (c < -1.0)
    c = -1.0;

  REAL phi = acos(c);
  *dphi = phi;

  if (phi < 1e-12) {
    nU[0] = 1.0;
    nU[1] = 0.0;
    nU[2] = 0.0;
    *dphi = 0.0;
  } else if (fabs(pi - phi) >= 1e-8) {
    const REAL sinphi = sin(phi);
    const REAL denom = 2.0 * sinphi;
    if (fabs(denom) < 1e-15) {
      fprintf(stderr,
              "ERROR in %s: unstable general branch (|2 sin(phi)| too small). phi=%.17e denom=%.17e\n",
              __func__, phi, denom);
      exit(1);
    }

    nU[0] = (R[2][1] - R[1][2]) / denom;
    nU[1] = (R[0][2] - R[2][0]) / denom;
    nU[2] = (R[1][0] - R[0][1]) / denom;
  } else {
    REAL n0 = 0.5 * (R[0][0] + 1.0);
    REAL n1 = 0.5 * (R[1][1] + 1.0);
    REAL n2 = 0.5 * (R[2][2] + 1.0);
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

    nU[0] = n0;
    nU[1] = n1;
    nU[2] = n2;

    REAL nnorm = sqrt(nU[0] * nU[0] + nU[1] * nU[1] + nU[2] * nU[2]);
    if (nnorm < 1e-15) {
      // Deterministic fallback in fully degenerate pi-branch.
      if (R[0][0] >= R[1][1] && R[0][0] >= R[2][2]) {
        nU[0] = 1.0;
        nU[1] = 0.0;
        nU[2] = 0.0;
      } else if (R[1][1] >= R[0][0] && R[1][1] >= R[2][2]) {
        nU[0] = 0.0;
        nU[1] = 1.0;
        nU[2] = 0.0;
      } else {
        nU[0] = 0.0;
        nU[1] = 0.0;
        nU[2] = 1.0;
      }
      nnorm = 1.0;
    }
    nU[0] /= nnorm;
    nU[1] /= nnorm;
    nU[2] /= nnorm;

    // Deterministic sign convention for pi-branch axis ambiguity.
    if (nU[0] < 0.0 || (fabs(nU[0]) <= 1e-16 && nU[1] < 0.0) ||
        (fabs(nU[0]) <= 1e-16 && fabs(nU[1]) <= 1e-16 && nU[2] < 0.0)) {
      nU[0] = -nU[0];
      nU[1] = -nU[1];
      nU[2] = -nU[2];
    }
    *dphi = pi;
  }

  REAL nnorm = sqrt(nU[0] * nU[0] + nU[1] * nU[1] + nU[2] * nU[2]);
  if (nnorm < 1e-15) {
    fprintf(stderr, "ERROR in %s: could not recover nonzero rotation axis.\n", __func__);
    exit(1);
  }
  nU[0] /= nnorm;
  nU[1] /= nnorm;
  nU[2] /= nnorm;

  // Debug reconstruction check: R_check(n,phi) should match input R.
  const REAL cphi = cos(*dphi);
  const REAL sphi = sin(*dphi);
  const REAL one_minus_c = 1.0 - cphi;
  const REAL nx = nU[0];
  const REAL ny = nU[1];
  const REAL nz = nU[2];

  REAL Rcheck[3][3];
  Rcheck[0][0] = cphi + one_minus_c * nx * nx;
  Rcheck[0][1] = one_minus_c * nx * ny - sphi * nz;
  Rcheck[0][2] = one_minus_c * nx * nz + sphi * ny;
  Rcheck[1][0] = one_minus_c * ny * nx + sphi * nz;
  Rcheck[1][1] = cphi + one_minus_c * ny * ny;
  Rcheck[1][2] = one_minus_c * ny * nz - sphi * nx;
  Rcheck[2][0] = one_minus_c * nz * nx - sphi * ny;
  Rcheck[2][1] = one_minus_c * nz * ny + sphi * nx;
  Rcheck[2][2] = cphi + one_minus_c * nz * nz;

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
            __func__, max_abs_err, *dphi, nU[0], nU[1], nU[2]);
    exit(1);
  }
""",
        include_CodeParameters_h=False,
    )


def register_CFunctions_so3_matrix_ops() -> None:
    """Register all SO(3) matrix helper C functions."""
    register_rotation_commondata_CodeParameters()
    if "build_R_from_cumulative_hats" not in cfc.CFunction_dict:
        register_CFunction_build_R_from_cumulative_hats()
    if "so3_validate_and_optionally_fix_hats" not in cfc.CFunction_dict:
        register_CFunction_so3_validate_and_optionally_fix_hats()
    if "so3_apply_R_to_vector" not in cfc.CFunction_dict:
        register_CFunction_so3_apply_R_to_vector()
    if "so3_apply_RT_to_vector" not in cfc.CFunction_dict:
        register_CFunction_so3_apply_RT_to_vector()
    if "so3_apply_R_to_tensorDD" not in cfc.CFunction_dict:
        register_CFunction_so3_apply_R_to_tensorDD()
    if "so3_apply_RT_to_tensorDD" not in cfc.CFunction_dict:
        register_CFunction_so3_apply_RT_to_tensorDD()
    if "so3_relative_rotation_dst_from_src" not in cfc.CFunction_dict:
        register_CFunction_so3_relative_rotation_dst_from_src()
    if "so3_matrix_to_axis_angle" not in cfc.CFunction_dict:
        register_CFunction_so3_matrix_to_axis_angle()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
