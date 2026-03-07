"""
Generate C code for SO(3) matrix helpers used by the BHaH rotation infrastructure.

This module registers small C helpers that:
1. build rotation matrices from cumulative hat vectors or axis-angle input,
2. apply rotation matrices or their transposes to vectors in place,
3. compute relative rotations between two cumulative rotation matrices, and
4. update cumulative hat vectors after a left multiplication by a new rotation.

SO(3) convention used here:
- ``R`` maps rotating-frame components to fixed-frame components.
- ``R[:,0]=xhat``, ``R[:,1]=yhat``, ``R[:,2]=zhat`` (all in fixed basis).
- ``v_fixed = R v_rot`` and ``v_rot = R^T v_fixed``.
- ``DeltaR_dst_from_src = R_dst^T R_src``.

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
    """
    Construct a symbolic rank-1 object with C-style indexed names.

    :param base_name: Base C identifier, e.g. ``"vU"``.
    :return: Three symbolic entries named ``base_name[i]``.
    """
    return [sp.Symbol(f"{base_name}[{i}]") for i in range(3)]


def _rank2_named(base_name: str) -> List[List[sp.Expr]]:
    """
    Construct a symbolic rank-2 object with C-style indexed names.

    :param base_name: Base C identifier, e.g. ``"R"``.
    :return: ``3 x 3`` symbolic matrix named ``base_name[i][j]``.
    """
    out = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            out[i][j] = sp.Symbol(f"{base_name}[{i}][{j}]")
    return out


def register_CFunction_so3_build_R_from_hats() -> None:
    """
    Register the generated C function ``so3_build_R_from_hats``.

    High-level behavior:
    - Interprets ``xhatU``, ``yhatU``, and ``zhatU`` as the cumulative rotating
      basis vectors expressed in the fixed basis.
    - Packs those hat vectors as the three columns of ``R``.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_build_R_from_hats()
    >>> "so3_build_R_from_hats" in cfc.CFunction_dict
    True
    """
    # Step 1: Set C function metadata and interface.
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Build a cumulative rotation matrix from cumulative hat vectors.

Convention:
- ``R`` maps rotating-frame components to fixed-frame components.
- ``R[:,0]=xhat``, ``R[:,1]=yhat``, ``R[:,2]=zhat``.
- C layout statement: ``R[i][j]`` is row ``i``, column ``j``.

@param[in] xhatU Cumulative rotating-frame x-basis vector in fixed components.
@param[in] yhatU Cumulative rotating-frame y-basis vector in fixed components.
@param[in] zhatU Cumulative rotating-frame z-basis vector in fixed components.
@param[out] R Rotation matrix assembled from the supplied hat vectors.
"""
    cfunc_type = "void"
    name = "so3_build_R_from_hats"
    params = (
        "const REAL xhatU[3], const REAL yhatU[3], const REAL zhatU[3], REAL R[3][3]"
    )

    # Step 2: Build symbolic matrix assignments from equation helpers.
    xhat_sym = _rank1_named("xhatU")
    yhat_sym = _rank1_named("yhatU")
    zhat_sym = _rank1_named("zhatU")
    R_expr = SO3Expressions.build_rotation_matrix_from_hats(
        xhat_sym, yhat_sym, zhat_sym
    )
    exprs = [R_expr[i][j] for i in range(3) for j in range(3)]
    lhses = [f"R[{i}][{j}]" for i in range(3) for j in range(3)]

    # Step 3: Construct the generated C body.
    body = r"""
  // Assemble the matrix column-by-column from cumulative hat vectors.
""" + ccg.c_codegen(
        exprs,
        lhses,
        include_braces=False,
        verbose=False,
    )

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


def register_CFunction_so3_axis_angle_to_R() -> None:
    """
    Register the generated C function ``so3_axis_angle_to_R``.

    High-level behavior:
    - Interprets ``nU`` and ``dphi`` as an axis-angle rotation.
    - Normalizes ``nU`` when possible and builds ``R`` via Rodrigues' formula.
    - Falls back to the identity matrix when the axis norm or angle is
      effectively zero.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_axis_angle_to_R()
    >>> "so3_axis_angle_to_R" in cfc.CFunction_dict
    True
    """
    # Step 1: Set C function metadata and interface.
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Build a rotation matrix from axis-angle input.

Convention:
- ``R`` maps rotating-frame components to fixed-frame components.
- ``v_fixed = R v_rot`` and ``v_rot = R^T v_fixed``.
- C layout statement: ``R[i][j]`` is row ``i``, column ``j``.

If ``|dphi|`` is effectively zero or ``|nU|`` is negligible, this routine
returns the identity matrix as a deterministic no-op rotation.

@param[in] nU Rotation axis, which need not be normalized.
@param[in] dphi Rotation angle in radians.
@param[out] R Rotation matrix computed from the axis-angle pair.
"""
    cfunc_type = "void"
    name = "so3_axis_angle_to_R"
    params = "const REAL nU[3], const REAL dphi, REAL R[3][3]"

    # Step 2: Build symbolic Rodrigues assignments from equation helpers.
    nU_unit_sym: List[sp.Expr] = [sp.Symbol("nx"), sp.Symbol("ny"), sp.Symbol("nz")]
    R_expr = SO3Expressions.rodrigues_matrix_from_unit_axis(
        nU_unit_sym, sp.Symbol("dphi")
    )
    exprs = [R_expr[i][j] for i in range(3) for j in range(3)]
    lhses = [f"R[{i}][{j}]" for i in range(3) for j in range(3)]
    R_codegen = ccg.c_codegen(
        exprs,
        lhses,
        include_braces=False,
        verbose=False,
    )

    # Step 3: Construct the generated C body.
    body = (
        r"""
  const REAL n0 = nU[0];
  const REAL n1 = nU[1];
  const REAL n2 = nU[2];
  const REAL nnorm = sqrt(n0 * n0 + n1 * n1 + n2 * n2);
  const REAL safe_nnorm = (nnorm < 1e-300) ? 1.0 : nnorm;

  // Normalize the input axis before evaluating Rodrigues' formula.
  const REAL nx = n0 / safe_nnorm;
  const REAL ny = n1 / safe_nnorm;
  const REAL nz = n2 / safe_nnorm;
"""
        + R_codegen
        + r"""
  if (fabs(dphi) < 1e-15 || nnorm < 1e-300) {
    // Deterministic identity fallback for degenerate axis-angle input.
    R[0][0] = 1.0;
    R[0][1] = 0.0;
    R[0][2] = 0.0;
    R[1][0] = 0.0;
    R[1][1] = 1.0;
    R[1][2] = 0.0;
    R[2][0] = 0.0;
    R[2][1] = 0.0;
    R[2][2] = 1.0;
  } // END IF: degenerate axis-angle input
"""
    )

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


def register_CFunction_so3_apply_R_to_vector() -> None:
    """
    Register the generated C function ``so3_apply_R_to_vector``.

    High-level behavior:
    - Applies ``R`` to ``vU`` in place.
    - Uses an alias-safe temporary copy of the input vector before overwriting
      ``vU``.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_apply_R_to_vector()
    >>> "so3_apply_R_to_vector" in cfc.CFunction_dict
    True
    """
    # Step 1: Set C function metadata and interface.
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Apply a rotation matrix to a vector in place.

Convention:
- ``R`` maps rotating-frame components to fixed-frame components.
- Vector update in Einstein notation: ``v^i_out = R^i{}_j v^j_in``.
- C layout statement: ``R[i][j]`` is row ``i``, column ``j``.

This routine is alias-safe: it snapshots the input vector before overwriting
``vU``.

@param[in] R Rotation matrix.
@param[in,out] vU Vector updated in place with ``R * vU``.
"""
    cfunc_type = "void"
    name = "so3_apply_R_to_vector"
    params = "const REAL R[3][3], REAL vU[3]"

    # Step 2: Build symbolic vector-rotation assignments.
    R_sym = _rank2_named("R")
    v_in_sym = _rank1_named("v_in")
    v_expr = SO3Expressions.apply_R_to_vector(R_sym, v_in_sym)
    exprs = [v_expr[i] for i in range(3)]
    lhses = [f"vU[{i}]" for i in range(3)]

    # Step 3: Construct the generated C body.
    body = r"""
  // Snapshot the input vector so the in-place update is alias-safe.
  const REAL v_in[3] = {vU[0], vU[1], vU[2]};
""" + ccg.c_codegen(
        exprs,
        lhses,
        include_braces=False,
        verbose=False,
    )

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


def register_CFunction_so3_apply_RT_to_vector() -> None:
    """
    Register the generated C function ``so3_apply_RT_to_vector``.

    High-level behavior:
    - Applies ``R^T`` to ``vU`` in place.
    - Uses an alias-safe temporary copy of the input vector before overwriting
      ``vU``.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_apply_RT_to_vector()
    >>> "so3_apply_RT_to_vector" in cfc.CFunction_dict
    True
    """
    # Step 1: Set C function metadata and interface.
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Apply the transpose of a rotation matrix to a vector in place.

Convention:
- ``R`` maps rotating-frame components to fixed-frame components.
- ``R^T`` maps fixed-frame components to rotating-frame components.
- Vector update in Einstein notation: ``v^i_out = (R^T)^i{}_j v^j_in``.
- C layout statement: ``R[i][j]`` is row ``i``, column ``j``.

This routine is alias-safe: it snapshots the input vector before overwriting
``vU``.

@param[in] R Rotation matrix whose transpose is applied.
@param[in,out] vU Vector updated in place with ``R^T * vU``.
"""
    cfunc_type = "void"
    name = "so3_apply_RT_to_vector"
    params = "const REAL R[3][3], REAL vU[3]"

    # Step 2: Build symbolic transpose-vector assignments.
    R_sym = _rank2_named("R")
    v_in_sym = _rank1_named("v_in")
    v_expr = SO3Expressions.apply_RT_to_vector(R_sym, v_in_sym)
    exprs = [v_expr[i] for i in range(3)]
    lhses = [f"vU[{i}]" for i in range(3)]

    # Step 3: Construct the generated C body.
    body = r"""
  // Snapshot the input vector so the in-place update is alias-safe.
  const REAL v_in[3] = {vU[0], vU[1], vU[2]};
""" + ccg.c_codegen(
        exprs,
        lhses,
        include_braces=False,
        verbose=False,
    )

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


def register_CFunction_so3_relative_R_dst_from_src() -> None:
    """
    Register the generated C function ``so3_relative_R_dst_from_src``.

    High-level behavior:
    - Interprets ``R_dst`` and ``R_src`` as cumulative rotating-to-fixed
      matrices for the destination and source frames.
    - Computes the relative rotation that maps source rotating-basis components
      to destination rotating-basis components.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_relative_R_dst_from_src()
    >>> "so3_relative_R_dst_from_src" in cfc.CFunction_dict
    True
    """
    # Step 1: Set C function metadata and interface.
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Compute the relative rotation from a source rotating basis to a destination rotating basis.

Convention:
- ``R_dst`` and ``R_src`` map rotating-frame components to fixed-frame components.
- ``DeltaR_dst_from_src = R_dst^T R_src``.
- ``DeltaR_dst_from_src`` maps source rotating-basis components to destination
  rotating-basis components.
- C layout statement: ``R[i][j]`` and ``DeltaR[i][j]`` are row-major by index.

@param[in] R_dst Destination cumulative rotation matrix.
@param[in] R_src Source cumulative rotation matrix.
@param[out] DeltaR Relative rotation matrix ``R_dst^T R_src``.
"""
    cfunc_type = "void"
    name = "so3_relative_R_dst_from_src"
    params = "const REAL R_dst[3][3], const REAL R_src[3][3], REAL DeltaR[3][3]"

    # Step 2: Build symbolic relative-rotation assignments.
    R_dst_sym = _rank2_named("R_dst")
    R_src_sym = _rank2_named("R_src")
    delta_expr = SO3Expressions.relative_rotation_dst_from_src(R_dst_sym, R_src_sym)
    exprs = [delta_expr[i][j] for i in range(3) for j in range(3)]
    lhses = [f"DeltaR[{i}][{j}]" for i in range(3) for j in range(3)]

    # Step 3: Construct the generated C body.
    body = r"""
  // Form the destination-from-source relative rotation matrix.
""" + ccg.c_codegen(
        exprs,
        lhses,
        include_braces=False,
        verbose=False,
    )

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


def register_CFunction_so3_left_multiply_hats_with_R() -> None:
    """
    Register the generated C function ``so3_left_multiply_hats_with_R``.

    High-level behavior:
    - Treats ``xhatU``, ``yhatU``, and ``zhatU`` as the columns of the current
      cumulative rotation matrix.
    - Applies ``DeltaR`` on the left so the updated hats represent
      ``R_new = DeltaR * R_old``.
    - Uses alias-safe temporary copies of all three input hat vectors.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_so3_left_multiply_hats_with_R()
    >>> "so3_left_multiply_hats_with_R" in cfc.CFunction_dict
    True
    """
    # Step 1: Set C function metadata and interface.
    includes = ["BHaH_defines.h"]
    desc = r"""
@brief Update cumulative hat vectors by left-multiplying with a rotation matrix.

Convention:
- ``R_old`` is the cumulative rotating-to-fixed matrix whose columns are the
  input hat vectors.
- ``R_new = DeltaR * R_old``.
- Each updated hat vector is one column of ``R_new``.
- C layout statement: ``DeltaR[i][j]`` is row ``i``, column ``j``.

This routine is alias-safe: it snapshots all input hat vectors before
overwriting them in place.

@param[in,out] xhatU Cumulative rotating-frame x-basis vector, updated in place.
@param[in,out] yhatU Cumulative rotating-frame y-basis vector, updated in place.
@param[in,out] zhatU Cumulative rotating-frame z-basis vector, updated in place.
@param[in] DeltaR Rotation matrix that left-multiplies the cumulative basis.
"""
    cfunc_type = "void"
    name = "so3_left_multiply_hats_with_R"
    params = "REAL xhatU[3], REAL yhatU[3], REAL zhatU[3], const REAL DeltaR[3][3]"

    # Step 2: Build symbolic hat-update assignments.
    DeltaR_sym = _rank2_named("DeltaR")
    xhat_in_sym = _rank1_named("xhat_in")
    yhat_in_sym = _rank1_named("yhat_in")
    zhat_in_sym = _rank1_named("zhat_in")
    xhat_expr = SO3Expressions.apply_R_to_vector(DeltaR_sym, xhat_in_sym)
    yhat_expr = SO3Expressions.apply_R_to_vector(DeltaR_sym, yhat_in_sym)
    zhat_expr = SO3Expressions.apply_R_to_vector(DeltaR_sym, zhat_in_sym)
    exprs = [*xhat_expr, *yhat_expr, *zhat_expr]
    lhses = [f"xhatU[{i}]" for i in range(3)]
    lhses += [f"yhatU[{i}]" for i in range(3)]
    lhses += [f"zhatU[{i}]" for i in range(3)]

    # Step 3: Construct the generated C body.
    body = r"""
  // Snapshot the input hat vectors so the in-place update is alias-safe.
  const REAL xhat_in[3] = {xhatU[0], xhatU[1], xhatU[2]};
  const REAL yhat_in[3] = {yhatU[0], yhatU[1], yhatU[2]};
  const REAL zhat_in[3] = {zhatU[0], zhatU[1], zhatU[2]};
""" + ccg.c_codegen(
        exprs,
        lhses,
        include_braces=False,
        verbose=False,
    )

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
