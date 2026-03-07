"""
Register SO(3)-native rotation helpers for BHaH.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import List, Tuple

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp
from nrpy.equations.rotation.SO3_rotations import SO3Expressions


def _rank1_named(base_name: str) -> List[sp.Expr]:
    return [sp.Symbol(f"{base_name}[{i}]") for i in range(3)]


def _rank2_named(base_name: str) -> List[List[sp.Expr]]:
    out = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            out[i][j] = sp.Symbol(f"{base_name}[{i}][{j}]")
    return out


def register_CFunction_so3_build_R_from_hats() -> None:
    """Register `so3_build_R_from_hats`."""
    xhat_sym = _rank1_named("xhatU")
    yhat_sym = _rank1_named("yhatU")
    zhat_sym = _rank1_named("zhatU")
    R_expr = SO3Expressions.build_rotation_matrix_from_hats(
        xhat_sym, yhat_sym, zhat_sym
    )
    R_exprs = [R_expr[i][j] for i in range(3) for j in range(3)]
    R_lhses = [f"R[{i}][{j}]" for i in range(3) for j in range(3)]
    body = ccg.c_codegen(R_exprs, R_lhses, include_braces=False, verbose=False)
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=["BHaH_defines.h"],
        desc=r"""
@brief Build `R` from cumulative hats.

Convention:
- `R` maps rotating-frame components -> fixed-frame components.
- `R[:,0]=xhat`, `R[:,1]=yhat`, `R[:,2]=zhat`.
""",
        cfunc_type="void",
        name="so3_build_R_from_hats",
        params=(
            "const REAL xhatU[3], const REAL yhatU[3], "
            "const REAL zhatU[3], REAL R[3][3]"
        ),
        body=body,
        include_CodeParameters_h=False,
    )


def register_CFunction_so3_axis_angle_to_R() -> None:
    """Register `so3_axis_angle_to_R`."""
    n_unit_sym: List[sp.Expr] = [sp.Symbol("nx"), sp.Symbol("ny"), sp.Symbol("nz")]
    R_expr = SO3Expressions.rodrigues_matrix_from_unit_axis(
        n_unit_sym, sp.Symbol("dphi")
    )
    R_exprs = [R_expr[i][j] for i in range(3) for j in range(3)]
    R_lhses = [f"R[{i}][{j}]" for i in range(3) for j in range(3)]
    R_codegen = ccg.c_codegen(R_exprs, R_lhses, include_braces=False, verbose=False)
    body = r"""
  const REAL n0 = nU[0];
  const REAL n1 = nU[1];
  const REAL n2 = nU[2];
  const REAL nnorm = sqrt(n0 * n0 + n1 * n1 + n2 * n2);

  if (fabs(dphi) < 1e-15 || nnorm < 1e-300) {
    R[0][0] = 1.0;
    R[0][1] = 0.0;
    R[0][2] = 0.0;
    R[1][0] = 0.0;
    R[1][1] = 1.0;
    R[1][2] = 0.0;
    R[2][0] = 0.0;
    R[2][1] = 0.0;
    R[2][2] = 1.0;
    return;
  }

  const REAL nx = n0 / nnorm;
  const REAL ny = n1 / nnorm;
  const REAL nz = n2 / nnorm;
""" + R_codegen
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=["BHaH_defines.h"],
        desc=r"""
@brief Build `R` from axis-angle input.

If the angle is effectively zero or the axis norm is negligible, this routine
returns the identity matrix.
""",
        cfunc_type="void",
        name="so3_axis_angle_to_R",
        params="const REAL nU[3], const REAL dphi, REAL R[3][3]",
        body=body,
        include_CodeParameters_h=False,
    )


def _register_vector_map_helper(name: str, apply_transpose: bool) -> None:
    R_sym = _rank2_named("R")
    v_in_sym = _rank1_named("v_in")
    v_expr = (
        SO3Expressions.apply_RT_to_vector(R_sym, v_in_sym)
        if apply_transpose
        else SO3Expressions.apply_R_to_vector(R_sym, v_in_sym)
    )
    body = r"""
  const REAL v_in[3] = {vU[0], vU[1], vU[2]};
""" + ccg.c_codegen(
        [v_expr[i] for i in range(3)],
        [f"vU[{i}]" for i in range(3)],
        include_braces=False,
        verbose=False,
    )
    desc = (
        "@brief Apply `R^T` to a vector in place."
        if apply_transpose
        else "@brief Apply `R` to a vector in place."
    )
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=["BHaH_defines.h"],
        desc=desc,
        cfunc_type="void",
        name=name,
        params="const REAL R[3][3], REAL vU[3]",
        body=body,
        include_CodeParameters_h=False,
    )


def register_CFunction_so3_apply_R_to_vector() -> None:
    """Register `so3_apply_R_to_vector`."""
    _register_vector_map_helper("so3_apply_R_to_vector", apply_transpose=False)


def register_CFunction_so3_apply_RT_to_vector() -> None:
    """Register `so3_apply_RT_to_vector`."""
    _register_vector_map_helper("so3_apply_RT_to_vector", apply_transpose=True)


def register_CFunction_so3_relative_R_dst_from_src() -> None:
    """Register `so3_relative_R_dst_from_src`."""
    R_dst_sym = _rank2_named("R_dst")
    R_src_sym = _rank2_named("R_src")
    delta_expr = SO3Expressions.relative_rotation_dst_from_src(R_dst_sym, R_src_sym)
    body = ccg.c_codegen(
        [delta_expr[i][j] for i in range(3) for j in range(3)],
        [f"DeltaR[{i}][{j}]" for i in range(3) for j in range(3)],
        include_braces=False,
        verbose=False,
    )
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=["BHaH_defines.h"],
        desc=r"""
@brief Compute `DeltaR_dst_from_src = R_dst^T R_src`.
""",
        cfunc_type="void",
        name="so3_relative_R_dst_from_src",
        params=("const REAL R_dst[3][3], const REAL R_src[3][3], REAL DeltaR[3][3]"),
        body=body,
        include_CodeParameters_h=False,
    )


def register_CFunction_so3_left_multiply_hats_with_R() -> None:
    """Register `so3_left_multiply_hats_with_R`."""
    R_sym = _rank2_named("DeltaR")
    xhat_in_sym = _rank1_named("xhat_in")
    yhat_in_sym = _rank1_named("yhat_in")
    zhat_in_sym = _rank1_named("zhat_in")
    xhat_expr = SO3Expressions.apply_R_to_vector(R_sym, xhat_in_sym)
    yhat_expr = SO3Expressions.apply_R_to_vector(R_sym, yhat_in_sym)
    zhat_expr = SO3Expressions.apply_R_to_vector(R_sym, zhat_in_sym)
    exprs = [*xhat_expr, *yhat_expr, *zhat_expr]
    lhses = [f"xhatU[{i}]" for i in range(3)]
    lhses += [f"yhatU[{i}]" for i in range(3)]
    lhses += [f"zhatU[{i}]" for i in range(3)]
    body = r"""
  const REAL xhat_in[3] = {xhatU[0], xhatU[1], xhatU[2]};
  const REAL yhat_in[3] = {yhatU[0], yhatU[1], yhatU[2]};
  const REAL zhat_in[3] = {zhatU[0], zhatU[1], zhatU[2]};
""" + ccg.c_codegen(exprs, lhses, include_braces=False, verbose=False)
    cfc.register_CFunction(
        subdirectory="rotation",
        includes=["BHaH_defines.h"],
        desc=r"""
@brief Update cumulative hats by left-multiplying with `DeltaR`.

This applies `R_new = DeltaR * R_old`, where the hat vectors are the columns of
`R_old`.
""",
        cfunc_type="void",
        name="so3_left_multiply_hats_with_R",
        params=("REAL xhatU[3], REAL yhatU[3], REAL zhatU[3], const REAL DeltaR[3][3]"),
        body=body,
        include_CodeParameters_h=False,
    )


def register_CFunctions() -> Tuple[str, ...]:
    r"""
    Register all public SO(3) helpers in this module.

    :return: Tuple of registered public helper names.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> names = register_CFunctions()
    >>> tuple(names)
    ('so3_build_R_from_hats', 'so3_axis_angle_to_R', 'so3_apply_R_to_vector', 'so3_apply_RT_to_vector', 'so3_relative_R_dst_from_src', 'so3_left_multiply_hats_with_R')
    >>> for func_name in names:
    ...     generated_str = cfc.CFunction_dict[func_name].full_function
    ...     validate_strings(generated_str, f"{func_name}__openmp", file_ext="c")
    """
    register_CFunction_so3_build_R_from_hats()
    register_CFunction_so3_axis_angle_to_R()
    register_CFunction_so3_apply_R_to_vector()
    register_CFunction_so3_apply_RT_to_vector()
    register_CFunction_so3_relative_R_dst_from_src()
    register_CFunction_so3_left_multiply_hats_with_R()
    return (
        "so3_build_R_from_hats",
        "so3_axis_angle_to_R",
        "so3_apply_R_to_vector",
        "so3_apply_RT_to_vector",
        "so3_relative_R_dst_from_src",
        "so3_left_multiply_hats_with_R",
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
