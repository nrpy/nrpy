"""
Symbolic Tensor (Quaternion) Rotation.

This script will perform symbolic tensor rotation on SymPy expressions using Quaternions.

Author: Ken Sible
Email:  ksible *at* outlook *dot* com
"""

from typing import List, Union, cast

import sympy as sp


def rotate(
    tensor: Union[List[sp.Basic], List[List[sp.Basic]]],
    axis: List[sp.Basic],
    angle: sp.Basic,
) -> Union[List[sp.Basic], List[List[sp.Basic]]]:
    """
    Rotate a symbolic vector or tensor about an arbitrary axis, using quaternions.

    :param tensor: 3-vector or (3x3)-matrix to be rotated.
    :param axis: Rotation axis (normal 3-vector).
    :param angle: Rotation angle (in radians).
    :return: Rotated tensor (of original type).
    :raises ValueError: If an unsupported tensor type or invalid size is provided.

    >>> from sympy import pi
    >>> one = sp.sympify(1)
    >>> zero = sp.sympify(0)
    >>> xhat = [one, zero, zero]
    >>> yhat = [zero, one, zero]
    >>> zhat = [zero, zero, one]
    >>> rotate(xhat, zhat, pi / 2)
    [0, 1, 0]
    >>> rotate(yhat, zhat, -pi / 2)
    [1, 0, 0]
    >>> import nrpy.indexedexp as ixp
    >>> gammaDD = ixp.zerorank2(dimension=3)
    >>> gammaDD[0][0] = gammaDD[1][1] = gammaDD[2][2] = one
    >>> for i in range(3):
    ...     for j in range(3):
    ...           gammaDD[i][j] = sp.Rational(1, 2)
    ...           if i==j:
    ...               gammaDD[i][j] = one
    >>> print(gammaDD)
    [[1, 1/2, 1/2], [1/2, 1, 1/2], [1/2, 1/2, 1]]
    >>> rotate(gammaDD, zhat, pi / 2)
    [[1, -1/2, -1/2], [-1/2, 1, 1/2], [-1/2, 1/2, 1]]
    >>> # Tests from nrpytutorial Jupyter notebook Tutorial-Symbolic_Tensor_Rotation.ipynb
    >>> v, angle, M = [1, 0, 1], pi/2, [[1, 2, 1], [0, 1, 0], [2, 1, 2]]
    >>> # vector rotation about x-axis
    >>> expected = sp.rot_axis1(-angle) * sp.Matrix(v)
    >>> received = sp.Matrix(rotate(v, [1, 0, 0], angle))
    >>> assert expected == received; v_ = received
    >>> # matrix rotation about x-axis
    >>> expected = sp.rot_axis1(-angle) * sp.Matrix(M) * sp.transpose(sp.rot_axis1(-angle))
    >>> received = sp.Matrix(rotate(M, [1, 0, 0], angle))
    >>> assert expected == received; M_ = received
    >>> # vector rotation about y-axis
    >>> expected = sp.rot_axis2(-angle) * sp.Matrix(v)
    >>> received = sp.Matrix(rotate(v, [0, 1, 0], angle))
    >>> assert expected == received; v_ = received
    >>> # matrix rotation about y-axis
    >>> expected = sp.rot_axis2(-angle) * sp.Matrix(M) * sp.transpose(sp.rot_axis2(-angle))
    >>> received = sp.Matrix(rotate(M, [0, 1, 0], angle))
    >>> assert expected == received; M_ = received
    >>> # vector rotation about z-axis
    >>> expected = sp.rot_axis3(-angle) * sp.Matrix(v)
    >>> received = sp.Matrix(rotate(v, [0, 0, 1], angle))
    >>> assert expected == received; v_ = received
    >>> # matrix rotation about z-axis
    >>> expected = sp.rot_axis3(-angle) * sp.Matrix(M) * sp.transpose(sp.rot_axis3(-angle))
    >>> received = sp.Matrix(rotate(M, [0, 0, 1], angle))
    >>> assert expected == received; M_ = received
    >>> # Confirm that rotating a symmetric rank-2 tensor hUU does not affect its symmetry.
    >>> angle = sp.symbols('theta', real=True)
    >>> vU    = ixp.declarerank1("vU")
    >>> hUU   = ixp.declarerank2("hUU", symmetry="sym01")
    >>> rotatedhUU, N = rotate(hUU, vU, angle), len(hUU)
    >>> for i in range(N):
    ...     for j in range(N):
    ...         if j >= i: continue
    ...         assert sp.simplify(rotatedhUU[i][j] - rotatedhUU[j][i]) == 0
    ...         print(f'Assertion Passed: rotatedhUU[{i}][{j}] == rotatedhUU[{j}][{i}]')
    Assertion Passed: rotatedhUU[1][0] == rotatedhUU[0][1]
    Assertion Passed: rotatedhUU[2][0] == rotatedhUU[0][2]
    Assertion Passed: rotatedhUU[2][1] == rotatedhUU[1][2]
    """

    def mul(*args: Union[List[sp.Quaternion], sp.Quaternion]) -> List[sp.Quaternion]:
        """
        Perform Quaternion-Matrix Multiplication.

        :param args: Either a Quaternion and a list, or a list and a Quaternion.
        :return: The result of Quaternion-matrix multiplication.
        """
        if isinstance(args[0], list):
            q, M = args[1], args[0]
            for i, col in enumerate(M):
                M[i] = col * q
        else:
            q, M = args[0], cast(List[sp.Quaternion], args[1])
            for i, col in enumerate(M):
                M[i] = q * col
        return M

    # Rotation Quaternion (Axis, Angle)
    q = sp.Quaternion.from_axis_angle(axis, angle)
    if isinstance(tensor[0], list):
        Matrix_from_tensor = sp.Matrix(tensor)
        if Matrix_from_tensor.shape != (3, 3):
            raise ValueError("Invalid sp.Matrix Size: Expected a (3x3) matrix.")
        # Rotation Formula: M' = (q.(q.M.q*)^T.q*)^T
        M = [
            sp.Quaternion(0, *Matrix_from_tensor[:, i])
            for i in range(Matrix_from_tensor.shape[1])
        ]
        M = mul(q, mul(M, q.conjugate()))
        for i in range(Matrix_from_tensor.shape[1]):
            Matrix_from_tensor[:, i] = [M[i].b, M[i].c, M[i].d]
        M = [
            sp.Quaternion(0, *Matrix_from_tensor[i, :])
            for i in range(Matrix_from_tensor.shape[0])
        ]
        M = mul(q, mul(M, q.conjugate()))
        for i in range(Matrix_from_tensor.shape[0]):
            Matrix_from_tensor[i, :] = [[M[i].b, M[i].c, M[i].d]]
        return cast(List[List[sp.Basic]], Matrix_from_tensor.tolist())

    if isinstance(tensor, list):
        if len(tensor) != 3:
            raise ValueError("Invalid Vector Length: Expected a 3-vector.")
        # Rotation Formula: v' = q.v.q*
        v = q * sp.Quaternion(0, *tensor) * q.conjugate()
        return [v.b, v.c, v.d]

    raise ValueError("Unsupported Tensor Type: Expected a 3-vector or a (3x3) matrix.")


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
