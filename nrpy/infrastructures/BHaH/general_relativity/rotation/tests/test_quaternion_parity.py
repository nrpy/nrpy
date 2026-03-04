"""Parity gates for BSSN Cartesian-basis rotation helpers."""

from __future__ import annotations

import math
import random

import sympy as sp

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.equations.quaternion_rotations.tensor_rotation import rotate
from nrpy.infrastructures.BHaH.general_relativity.rotation.rotate_BSSN_Cartesian_basis import (
    register_CFunction_rotate_BSSN_Cartesian_basis,
)
from nrpy.infrastructures.BHaH.general_relativity.rotation.rotate_BSSN_Cartesian_basis_by_R import (
    register_CFunction_rotate_BSSN_Cartesian_basis_by_R,
)


def _as_float(x: sp.Expr | float | int) -> float:
    return float(sp.sympify(x))


def _to_real_float(x: object) -> float:
    return float(sp.re(sp.sympify(x)))


def _rodrigues_matrix_eval(
    nU_input: list[sp.Expr | float | int], dphi_input: sp.Expr | float | int
) -> list[list[float]]:
    n0 = _as_float(nU_input[0])
    n1 = _as_float(nU_input[1])
    n2 = _as_float(nU_input[2])
    axis_sq = n0 * n0 + n1 * n1 + n2 * n2
    if axis_sq < 1.0e-600:
        return [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]

    nnorm = math.sqrt(axis_sq)
    nx = n0 / nnorm
    ny = n1 / nnorm
    nz = n2 / nnorm
    dphi = _as_float(dphi_input)
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


def _apply_rotation_to_vector_eval(
    R_input: list[list[sp.Expr | float | int]],
    v_input: list[sp.Expr | float | int],
) -> list[float]:
    out: list[float] = [0.0, 0.0, 0.0]
    for i in range(3):
        out_i = 0.0
        for j in range(3):
            out_i += _as_float(R_input[i][j]) * _as_float(v_input[j])
        out[i] = out_i
    return out


def _apply_rotation_to_tensor_eval(
    R_input: list[list[sp.Expr | float | int]],
    t_input: list[list[sp.Expr | float | int]],
) -> list[list[float]]:
    out: list[list[float]] = [[0.0, 0.0, 0.0] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            out_ij = 0.0
            for k in range(3):
                for l in range(3):
                    out_ij += (
                        _as_float(R_input[i][k])
                        * _as_float(t_input[k][l])
                        * _as_float(R_input[j][l])
                    )
            out[i][j] = out_ij
    return out


def _verify_quaternion_interface_parity() -> None:
    axis = [sp.Rational(3, 10), -sp.Rational(2, 5), sp.Rational(1, 2)]
    dphi = sp.Rational(6, 5)
    vec = [1, -2, sp.Rational(1, 2)]
    ten = [
        [2, sp.Rational(1, 10), -sp.Rational(3, 10)],
        [sp.Rational(1, 10), sp.Rational(3, 2), sp.Rational(7, 10)],
        [-sp.Rational(3, 10), sp.Rational(7, 10), sp.Rational(4, 5)],
    ]

    R = _rodrigues_matrix_eval(axis, dphi)
    vec_by_R = _apply_rotation_to_vector_eval(R, vec)
    vec_by_q = [_to_real_float(v) for v in rotate(vec, axis, dphi)]
    if max(abs(a - b) for a, b in zip(vec_by_R, vec_by_q)) >= 1e-13:
        raise AssertionError("Vector parity check failed for nominal case")

    ten_by_R = _apply_rotation_to_tensor_eval(R, ten)
    ten_by_q = [[_to_real_float(v) for v in row] for row in rotate(ten, axis, dphi)]
    if (
        max(abs(ten_by_R[i][j] - ten_by_q[i][j]) for i in range(3) for j in range(3))
        >= 1e-13
    ):
        raise AssertionError("Tensor parity check failed for nominal case")

    dphi = sp.pi - sp.sympify("1e-10")
    R = _rodrigues_matrix_eval(axis, dphi)
    vec_by_R = _apply_rotation_to_vector_eval(R, vec)
    vec_by_q = [_to_real_float(v) for v in rotate(vec, axis, dphi)]
    if max(abs(a - b) for a, b in zip(vec_by_R, vec_by_q)) >= 1e-11:
        raise AssertionError("Vector parity check failed for near-pi case")


def _verify_quaternion_interface_parity_randomized() -> None:
    rng = random.Random(20260303)

    def rand_axis() -> list[float]:
        while True:
            axis = [rng.uniform(-1.0, 1.0) for _ in range(3)]
            if axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2] > 1e-24:
                return axis

    angles = [
        0.0,
        1e-14,
        -1e-14,
        sp.pi - sp.sympify("1e-12"),
        -sp.pi + sp.sympify("1e-12"),
        0.7,
        -1.1,
    ]
    max_vec_err = 0.0
    max_ten_err = 0.0

    for idx in range(64):
        axis = rand_axis()
        if idx < len(angles):
            dphi = angles[idx % len(angles)]
        else:
            dphi = rng.uniform(-3.141592653589793, 3.141592653589793)
        vec = [rng.uniform(-2.0, 2.0) for _ in range(3)]
        ten = [[rng.uniform(-2.0, 2.0) for _ in range(3)] for _ in range(3)]
        ten = [[0.5 * (ten[i][j] + ten[j][i]) for j in range(3)] for i in range(3)]

        R = _rodrigues_matrix_eval(axis, dphi)
        vec_by_R = _apply_rotation_to_vector_eval(R, vec)
        vec_by_q = [_to_real_float(v) for v in rotate(vec, axis, dphi)]

        ten_by_R = _apply_rotation_to_tensor_eval(R, ten)
        ten_by_q = [[_to_real_float(v) for v in row] for row in rotate(ten, axis, dphi)]

        max_vec_err = max(
            max_vec_err, max(abs(a - b) for a, b in zip(vec_by_R, vec_by_q))
        )
        max_ten_err = max(
            max_ten_err,
            max(
                abs(ten_by_R[i][j] - ten_by_q[i][j]) for i in range(3) for j in range(3)
            ),
        )

    if max_vec_err >= 2e-11:
        raise AssertionError(
            f"Randomized vector parity exceeded tolerance: {max_vec_err}"
        )
    if max_ten_err >= 5e-11:
        raise AssertionError(
            f"Randomized tensor parity exceeded tolerance: {max_ten_err}"
        )


def run_quaternion_parity_and_convention_gates() -> bool:
    r"""
    Run deterministic quaternion/matrix parity checks and registration sanity.

    Doctests:
    >>> run_quaternion_parity_and_convention_gates()
    True
    """
    _verify_quaternion_interface_parity()
    _verify_quaternion_interface_parity_randomized()

    par.set_parval_from_str("parallelization", "openmp")
    cfc.CFunction_dict.clear()
    register_CFunction_rotate_BSSN_Cartesian_basis_by_R()
    register_CFunction_rotate_BSSN_Cartesian_basis()

    return True


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
