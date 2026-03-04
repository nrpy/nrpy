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

import math

import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH.rotation.so3_matrix_ops import (
    register_CFunction_so3_apply_R_to_tensorDD,
    register_CFunction_so3_apply_R_to_vector,
)


def _rodrigues_matrix_from_axis_angle(
    nU: list[float], dphi: float
) -> list[list[float]]:
    """
    Build a 3x3 Rodrigues rotation matrix from axis-angle input.

    :param nU: Rotation axis components.
    :param dphi: Rotation angle in radians.
    :return: Rodrigues rotation matrix.
    """
    nx, ny, nz = nU
    nnorm = math.sqrt(nx * nx + ny * ny + nz * nz)
    if nnorm < 1e-300:
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    nx /= nnorm
    ny /= nnorm
    nz /= nnorm
    c = math.cos(dphi)
    s = math.sin(dphi)
    one_minus_c = 1.0 - c
    return [
        [
            c + one_minus_c * nx * nx,
            one_minus_c * nx * ny - s * nz,
            one_minus_c * nx * nz + s * ny,
        ],
        [
            one_minus_c * ny * nx + s * nz,
            c + one_minus_c * ny * ny,
            one_minus_c * ny * nz - s * nx,
        ],
        [
            one_minus_c * nz * nx - s * ny,
            one_minus_c * nz * ny + s * nx,
            c + one_minus_c * nz * nz,
        ],
    ]


def _apply_R_to_vector(R: list[list[float]], vU: list[float]) -> list[float]:
    """
    Apply v_dst = R v_src.

    :param R: Rotation matrix.
    :param vU: Input vector components.
    :return: Output vector components after applying ``R``.
    """
    return [sum(R[i][j] * vU[j] for j in range(3)) for i in range(3)]


def _apply_R_to_tensorDD(
    R: list[list[float]], tDD: list[list[float]]
) -> list[list[float]]:
    """
    Apply T_dst = R T_src R^T.

    :param R: Rotation matrix.
    :param tDD: Input rank-2 tensor.
    :return: Output tensor after applying ``R`` and ``R^T``.
    """
    tmp = [
        [sum(R[i][k] * tDD[k][j] for k in range(3)) for j in range(3)] for i in range(3)
    ]
    return [
        [sum(tmp[i][k] * R[j][k] for k in range(3)) for j in range(3)] for i in range(3)
    ]


def verify_quaternion_interface_parity() -> None:
    r"""
    Verify axis-angle matrix convention against quaternion tensor rotation interface.

    Doctests:
    >>> from nrpy.equations.quaternion_rotations.tensor_rotation import rotate
    >>> axis = [0.3, -0.4, 0.5]
    >>> dphi = 1.2
    >>> vec = [1.0, -2.0, 0.5]
    >>> ten = [[2.0, 0.1, -0.3], [0.1, 1.5, 0.7], [-0.3, 0.7, 0.8]]
    >>> R = _rodrigues_matrix_from_axis_angle(axis, dphi)
    >>> vec_by_R = _apply_R_to_vector(R, vec)
    >>> vec_by_q = [float(v) for v in rotate(vec, axis, dphi)]
    >>> max(abs(a - b) for a, b in zip(vec_by_R, vec_by_q)) < 1e-13
    True
    >>> ten_by_R = _apply_R_to_tensorDD(R, ten)
    >>> ten_by_q = [[float(v) for v in row] for row in rotate(ten, axis, dphi)]
    >>> max(abs(ten_by_R[i][j] - ten_by_q[i][j]) for i in range(3) for j in range(3)) < 1e-13
    True
    >>> # Near-pi case
    >>> dphi = math.pi - 1e-10
    >>> R = _rodrigues_matrix_from_axis_angle(axis, dphi)
    >>> vec_by_R = _apply_R_to_vector(R, vec)
    >>> vec_by_q = [float(v) for v in rotate(vec, axis, dphi)]
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
    >>> def rand_axis() -> list[float]:
    ...     while True:
    ...         axis = [rng.uniform(-1.0, 1.0) for _ in range(3)]
    ...         if math.sqrt(sum(a * a for a in axis)) > 1e-12:
    ...             return axis
    >>> angles = [0.0, 1e-14, -1e-14, math.pi - 1e-12, -math.pi + 1e-12, 0.7, -1.1]
    >>> max_vec_err = 0.0
    >>> max_ten_err = 0.0
    >>> for idx in range(64):
    ...     axis = rand_axis()
    ...     dphi = angles[idx % len(angles)] if idx < len(angles) else rng.uniform(-math.pi, math.pi)
    ...     vec = [rng.uniform(-2.0, 2.0) for _ in range(3)]
    ...     ten = [[rng.uniform(-2.0, 2.0) for _ in range(3)] for _ in range(3)]
    ...     # Symmetrize to match BSSN tensor usage.
    ...     ten = [[0.5 * (ten[i][j] + ten[j][i]) for j in range(3)] for i in range(3)]
    ...     R = _rodrigues_matrix_from_axis_angle(axis, dphi)
    ...     vec_by_R = _apply_R_to_vector(R, vec)
    ...     vec_by_q = [float(v) for v in rotate(vec, axis, dphi)]
    ...     ten_by_R = _apply_R_to_tensorDD(R, ten)
    ...     ten_by_q = [[float(v) for v in row] for row in rotate(ten, axis, dphi)]
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
    >>> from nrpy.infrastructures.BHaH.rotation.so3_matrix_ops import assert_so3_convention_in_text
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_rotate_BSSN_Cartesian_basis_by_R()
    >>> generated_str = cfc.CFunction_dict["rotate_BSSN_Cartesian_basis_by_R"].full_function
    >>> assert_so3_convention_in_text(generated_str, "rotate_BSSN_Cartesian_basis_by_R")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "rotate_BSSN_Cartesian_basis_by_R_openmp", file_ext="c")
    """
    if "so3_apply_R_to_vector" not in cfc.CFunction_dict:
        register_CFunction_so3_apply_R_to_vector()
    if "so3_apply_R_to_tensorDD" not in cfc.CFunction_dict:
        register_CFunction_so3_apply_R_to_tensorDD()

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
- Vectors: v_dst = DeltaR_dst_from_src * v_src.
- Rank-2 tensors: T_dst = DeltaR_dst_from_src * T_src *
  DeltaR_dst_from_src^T.

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
  // Enforce symmetry by mirroring upper-triangular entries to lower-triangular.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < i; j++) {
      hDD[i][j] = hDD[j][i];
      aDD[i][j] = aDD[j][i];
    } // END LOOP over strictly-lower-triangular entries.
  } // END LOOP over tensor rows.

  REAL vetU_out[3], betU_out[3], lambdaU_out[3];
  REAL hDD_out[3][3], aDD_out[3][3];
  so3_apply_R_to_vector(DeltaR_dst_from_src, vetU, vetU_out);
  so3_apply_R_to_vector(DeltaR_dst_from_src, betU, betU_out);
  so3_apply_R_to_vector(DeltaR_dst_from_src, lambdaU, lambdaU_out);
  so3_apply_R_to_tensorDD(DeltaR_dst_from_src, hDD, hDD_out);
  so3_apply_R_to_tensorDD(DeltaR_dst_from_src, aDD, aDD_out);

  for (int i = 0; i < 3; i++) {
    vetU[i] = vetU_out[i];
    betU[i] = betU_out[i];
    lambdaU[i] = lambdaU_out[i];
    for (int j = 0; j < 3; j++) {
      hDD[i][j] = hDD_out[i][j];
      aDD[i][j] = aDD_out[i][j];
    }
  }

  // Re-enforce symmetry by averaging mirrored tensor entries.
  for (int i = 0; i < 3; i++) {
    for (int j = i + 1; j < 3; j++) {
      const REAL hsym = 0.5 * (hDD[i][j] + hDD[j][i]);
      const REAL asym = 0.5 * (aDD[i][j] + aDD[j][i]);
      hDD[i][j] = hsym;
      hDD[j][i] = hsym;
      aDD[i][j] = asym;
      aDD[j][i] = asym;
    } // END LOOP over strictly-upper-triangular entries.
  } // END LOOP over tensor rows.
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
    >>> from nrpy.infrastructures.BHaH.rotation.so3_matrix_ops import assert_so3_convention_in_text
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_rotate_BSSN_Cartesian_basis()
    >>> generated_str = cfc.CFunction_dict["rotate_BSSN_Cartesian_basis"].full_function
    >>> assert_so3_convention_in_text(generated_str, "rotate_BSSN_Cartesian_basis")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "rotate_BSSN_Cartesian_basis_openmp", file_ext="c")
    """
    if "rotate_BSSN_Cartesian_basis_by_R" not in cfc.CFunction_dict:
        register_CFunction_rotate_BSSN_Cartesian_basis_by_R()

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
    body = r"""
  const REAL n0 = nU[0];
  const REAL n1 = nU[1];
  const REAL n2 = nU[2];
  const REAL nnorm = sqrt(n0 * n0 + n1 * n1 + n2 * n2);

  REAL R[3][3];
  if (nnorm < 1e-300) {
    // Deterministic fallback for degenerate axis input.
    R[0][0] = 1.0; R[0][1] = 0.0; R[0][2] = 0.0;
    R[1][0] = 0.0; R[1][1] = 1.0; R[1][2] = 0.0;
    R[2][0] = 0.0; R[2][1] = 0.0; R[2][2] = 1.0;
  } else {
    const REAL nx = n0 / nnorm;
    const REAL ny = n1 / nnorm;
    const REAL nz = n2 / nnorm;
    const REAL c = cos(dphi);
    const REAL s = sin(dphi);
    const REAL one_minus_c = 1.0 - c;

    // Rodrigues formula: R = c I + (1-c) n n^T + s [n]_x.
    R[0][0] = c + one_minus_c * nx * nx;
    R[0][1] = one_minus_c * nx * ny - s * nz;
    R[0][2] = one_minus_c * nx * nz + s * ny;
    R[1][0] = one_minus_c * ny * nx + s * nz;
    R[1][1] = c + one_minus_c * ny * ny;
    R[1][2] = one_minus_c * ny * nz - s * nx;
    R[2][0] = one_minus_c * nz * nx - s * ny;
    R[2][1] = one_minus_c * nz * ny + s * nx;
    R[2][2] = c + one_minus_c * nz * nz;
  }

  rotate_BSSN_Cartesian_basis_by_R(vetU, betU, lambdaU, hDD, aDD, R);
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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
