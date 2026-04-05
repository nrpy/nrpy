"""
C function registration for computing the local 3D volume element sqrt(detgammahat) * d3xx at a point.

This module constructs and registers the C helper routine "sqrt_detgammahat_d3xx_volume_element".
The generated C function evaluates the positive volume element used in 3D integrals at the
provided local coordinates (xx0, xx1, xx2). The expression is:
    sqrt(detgammahat(xx0,xx1,xx2)) * abs(dxx0 * dxx1 * dxx2)
The grid spacings dxx{0,1,2} are taken from CodeParameters.h, and reference-metric data are
specialized at registration time for the selected coordinate system.

Function
--------
register_CFunction_sqrt_detgammahat_d3xx_volume_element
    Construct and register the "sqrt_detgammahat_d3xx_volume_element" C function.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)


def register_CFunction_sqrt_detgammahat_d3xx_volume_element(
    CoordSystem: str,
) -> None:
    """
    Construct and register a C function that computes the local 3D volume element at a point.

    This function generates and registers the C helper "sqrt_detgammahat_d3xx_volume_element". The
    generated C code evaluates sqrt(detgammahat) * abs(dxx0 * dxx1 * dxx2) at the provided local
    coordinates (xx0, xx1, xx2), using reference-metric data associated with the requested
    coordinate system. The dxx{0,1,2} grid spacings are declared in CodeParameters.h and pulled
    in via the params_struct definitions. For CUDA builds, the function is decorated with
    __host__ __device__; for OpenMP builds it is a plain C function. The result is returned by
    reference through the dV pointer in C.

    :param CoordSystem: Name of the coordinate system used to select the reference metric data
                        embedded in the generated C function.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "sqrt_detgammahat_d3xx_volume_element"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       register_CFunction_sqrt_detgammahat_d3xx_volume_element(CoordSystem)
    ...       generated_str = cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}"
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    # dxx{0,1,2} must be registered or it won't be declared at the top of the function.
    # fmt: off
    for j in range(3):
        _ = par.CodeParameter("REAL", __name__, f"dxx{j}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    # fmt: on

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
 * @file sqrt_detgammahat_d3xx_volume_element.c
 * @brief Compute the local 3D volume element sqrt(detgammahat) * d3xx at a point.
 *
 * This routine evaluates the positive volume element used when integrating scalar fields on a
 * 3D grid. The evaluation uses reference-metric data and grid spacings stored in the params
 * structure. The computed value is returned by reference via the dV pointer.
 *
 * The expression implemented is:
 *   dV = sqrt(detgammahat(xx0, xx1, xx2)) * abs(dxx0 * dxx1 * dxx2)
 * The absolute value ensures a positive volume element regardless of coordinate orientation.
 *
 * @param[in]  params  Pointer to parameter struct (reference-metric data, grid spacings, and sizes).
 * @param[in]  xx0     Local coordinate 0 at which to evaluate the volume element.
 * @param[in]  xx1     Local coordinate 1 at which to evaluate the volume element.
 * @param[in]  xx2     Local coordinate 2 at which to evaluate the volume element.
 * @param[out] dV      Pointer to the location where the volume element will be stored.
 *
 * @return     void. The result is written to *dV.
 *
 * Note: If a user-editable block is provided in the implementation, users may add custom logic,
 * such as scaling or additional diagnostics, prior to writing the result.
 """
    cfunc_type = "void"
    name = "sqrt_detgammahat_d3xx_volume_element"
    params = "const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict dV"
    rfm = refmetric.reference_metric[CoordSystem]
    # These are set in CodeParameters.h
    dxx0, dxx1, dxx2 = sp.symbols("dxx0 dxx1 dxx2", real=True)
    expr_list = [sp.sqrt(sp.simplify(rfm.detgammahat)) * sp.Abs(dxx0 * dxx1 * dxx2)]
    body = ccg.c_codegen(
        expr_list,
        ["*dV"],
        include_braces=False,
    )

    parallelization = par.parval_from_str("parallelization")
    param_symbols, _ = get_params_commondata_symbols_from_expr_list(expr_list)
    params_definitions = generate_definition_header(
        param_symbols,
        enable_intrinsics=False,
        var_access=parallel_utils.get_params_access("openmp"),
    )

    kernel_body = f"{params_definitions}\n{body}"
    cfunc_decorators = "__host__ __device__" if parallelization == "cuda" else ""
    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=kernel_body,
        cfunc_decorators=cfunc_decorators,
    )
