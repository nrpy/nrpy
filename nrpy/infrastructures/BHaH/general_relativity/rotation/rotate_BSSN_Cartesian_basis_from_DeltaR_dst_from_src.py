"""
Generate C code for rotating BSSN Cartesian-basis fields using a matrix.

This module registers the core C routine that rotates BSSN vectors and
symmetric rank-2 tensors using a supplied matrix
``DeltaR_dst_from_src``. The update is performed in place, with explicit
alias-safe temporaries and equation-derived SO(3) expressions.

SO(3) convention used here:
- ``R`` maps rotating-frame components to fixed-frame components.
- ``R[:,0]=xhat``, ``R[:,1]=yhat``, ``R[:,2]=zhat`` (all in fixed basis).
- ``v_fixed = R v_rot`` and ``v_rot = R^T v_fixed``.
- ``T_fixed = R T_rot R^T`` and ``T_rot = R^T T_fixed R``.
- ``DeltaR_dst_from_src = R_dst^T R_src``.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp
from nrpy.equations.rotation.SO3_rotations import SO3Expressions


def register_CFunction_rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src() -> None:
    """
    Register the generated C function
    ``rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src``.

    High-level behavior:
    - Treats ``DeltaR_dst_from_src`` as the basis-change matrix from source
      rotating components to destination rotating components.
    - Rotates vectors ``vetU``, ``betU``, and ``lambdaU`` in place.
    - Rotates symmetric tensors ``hDD`` and ``aDD`` in place.

    Symmetric tensor storage contract:
    - Input: upper-triangular components are authoritative.
    - Internal: full symmetric local tensors are reconstructed.
    - Output: only upper-triangular components are written back.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src()
    >>> generated_str = cfc.CFunction_dict["rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src"].full_function
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     validate_strings(generated_str, "rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src_openmp", file_ext="c")
    """
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
  components before applying equation-derived SO(3) tensor-rotation expressions.

@param[in,out] vetU BSSN rescaled shift vector, rotated in place.
@param[in,out] betU BSSN rescaled driver vector, rotated in place.
@param[in,out] lambdaU BSSN conformal connection vector, rotated in place.
@param[in,out] hDD BSSN conformal metric perturbation tensor, rotated in place.
@param[in,out] aDD BSSN trace-free extrinsic curvature tensor, rotated in place.
@param[in] DeltaR_dst_from_src Relative basis rotation matrix.
"""
    cfunc_type = "void"
    name = "rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src"
    params = (
        "REAL vetU[3], REAL betU[3], REAL lambdaU[3], "
        "REAL hDD[3][3], REAL aDD[3][3], const REAL DeltaR_dst_from_src[3][3]"
    )

    # Step 2: Build equation-driven vector/tensor rotation assignments.
    def _rank1_named(base_name: str) -> List[sp.Expr]:
        return [sp.Symbol(f"{base_name}[{i}]") for i in range(3)]

    def _rank2_named(base_name: str) -> List[List[sp.Expr]]:
        out = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                out[i][j] = sp.Symbol(f"{base_name}[{i}][{j}]")
        return out

    DeltaR_sym = _rank2_named("DeltaR_dst_from_src")

    def _vector_exprs_lhses(
        input_name: str, output_name: str
    ) -> tuple[List[sp.Expr], List[str]]:
        v_sym = _rank1_named(input_name)
        v_rot_expr = SO3Expressions.apply_R_to_vector(DeltaR_sym, v_sym)
        exprs = [v_rot_expr[i] for i in range(3)]
        lhses = [f"{output_name}[{i}]" for i in range(3)]
        return exprs, lhses

    def _tensor_exprs_lhses(
        input_name: str, output_name: str
    ) -> tuple[List[sp.Expr], List[str]]:
        t_sym = _rank2_named(input_name)
        t_rot_expr = SO3Expressions.apply_R_to_tensorDD(DeltaR_sym, t_sym)
        exprs = [t_rot_expr[i][j] for i in range(3) for j in range(3)]
        lhses = [f"{output_name}[{i}][{j}]" for i in range(3) for j in range(3)]
        return exprs, lhses

    all_exprs: List[sp.Expr] = []
    all_lhses: List[str] = []
    for exprs, lhses in (
        _vector_exprs_lhses("vetU_in", "vetU_out"),
        _vector_exprs_lhses("betU_in", "betU_out"),
        _vector_exprs_lhses("lambdaU_in", "lambdaU_out"),
        _tensor_exprs_lhses("hDD_sym", "hDD_out"),
        _tensor_exprs_lhses("aDD_sym", "aDD_out"),
    ):
        all_exprs.extend(exprs)
        all_lhses.extend(lhses)
    rotation_codegen = ccg.c_codegen(
        all_exprs, all_lhses, include_braces=False, verbose=False
    )

    # Step 3: Construct the generated C body.
    body = (
        r"""
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

  // Alias-safe vector rotations: v_dst = DeltaR * v_src.
  const REAL vetU_in[3] = {vetU[0], vetU[1], vetU[2]};
  const REAL betU_in[3] = {betU[0], betU[1], betU[2]};
  const REAL lambdaU_in[3] = {lambdaU[0], lambdaU[1], lambdaU[2]};
"""
        + rotation_codegen
        + r"""
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
