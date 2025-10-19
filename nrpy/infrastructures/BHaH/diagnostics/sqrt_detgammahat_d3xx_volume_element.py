"""
Register numerical_grids_and_timestep() C function, as well as functions called by this one.

These functions set up numerical grids for use within the BHaH infrastructure.

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
    Register a C function to find the minimum grid spacing ds_min.

    ds_min is the minimum spacing between neighboring gridpoints on a numerical grid.

    :param CoordSystem: The coordinate system of the numerical grid.

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
    Setting up reference_metric[SinhSymTP]...
    Setting up reference_metric[HoleySinhSpherical]...
    Setting up reference_metric[Cartesian]...
    Setting up reference_metric[SinhCylindricalv2n2]...
    """
    # dxx{0,1,2} must be registered or it won't be declared at the top of the function.
    # fmt: off
    for j in range(3):
        _ = par.CodeParameter("REAL", __name__, f"dxx{j}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    # fmt: on

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
 * @brief Compute the local 3D volume element sqrt{detgammahat} d^3xx at a point.
 *
 * Given a set of grid/metric parameters and local coordinates, this function evaluates
 * the positive volume element used for integrating scalar fields over a 3D grid.
 * The implementation uses precomputed parameters in @p params and elementary functions
 * of the local coordinates. The result is stored by reference in @p dV.
 *
 * @param params Pointer to the global/solver parameters required to evaluate the volume element.
 * @param xx0 Local coordinate 0 at which to evaluate the volume element.
 * @param xx1 Local coordinate 1 at which to evaluate the volume element.
 * @param xx2 Local coordinate 2 at which to evaluate the volume element.
 * @param dV  Output pointer set to the absolute value of the local volume element.
 * @return void This function returns no value; the result is written to @p dV.
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
