"""
C function registration for converting between grid coordinate (xx0,xx1,xx2)
(uniform grid spacing) to Cartesian coordinate (x,y,z).

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Set, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
import nrpy.reference_metric as refmetric


def _prepare_sympy_exprs_for_codegen(
    sympy_exprs: List[sp.Expr], local_vars: Set[str]
) -> List[sp.Expr]:
    """
    Substitute symbolic parameters with 'params->' prefix for C code generation.

    This function iterates through a list of SymPy expressions, identifies all unique
    free symbols that are not specified as local variables, and prepends 'params->'
    to their names. This is a common requirement for generating C code that
    accesses parameters from a 'params' struct.

    This implementation is efficient, creating a single substitution dictionary
    for all symbols across all expressions before applying it.

    :param sympy_exprs: A list of SymPy expressions to process.
    :param local_vars: A set of strings, where each string is the name of a
                       symbol that should be treated as a local variable (and thus
                       not be prefixed).
    :return: A new list of SymPy expressions with parameters appropriately substituted.

    Doctests:
    >>> import sympy as sp
    >>> # Define input symbols and expressions
    >>> a, b, c, xx0, xx1 = sp.symbols("a b c xx0 xx1")
    >>> input_exprs = [a + b*xx0, c - xx1]
    >>> local_vars = {"xx0", "xx1"}
    >>>
    >>> # Manually construct the expected SymPy objects for a robust comparison.
    >>> # The symbols "params->a", etc., are valid SymPy symbols.
    >>> params_a, params_b, params_c = sp.symbols("params->a params->b params->c")
    >>> expected_exprs_set = {xx0*params_b + params_a, -xx1 + params_c}
    >>>
    >>> # Generate the expressions and convert to a set for order-independent comparison.
    >>> generated_exprs_set = set(_prepare_sympy_exprs_for_codegen(input_exprs, local_vars))
    >>>
    >>> # The sets of expressions must be mathematically identical.
    >>> generated_exprs_set == expected_exprs_set
    True
    """
    # Create a single substitution dictionary for all unique parameters.
    all_symbols = set().union(*(expr.free_symbols for expr in sympy_exprs))
    substitutions = {
        symbol: sp.symbols(f"params->{symbol.name}")
        for symbol in all_symbols
        if symbol.name not in local_vars
    }

    # Apply the substitution to each expression.
    return [expr.subs(substitutions) for expr in sympy_exprs]


def register_CFunction__Cart_to_xx_and_nearest_i0i1i2(
    CoordSystem: str,
    relative_to: str = "local_grid_center",
    gridding_approach: str = "independent grid(s)",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register a C function that maps Cartesian coordinates to xx and finds the nearest grid indices.

    This function generates a C function which, given Cartesian coordinates (x, y, z),
    computes the corresponding (xx0, xx1, xx2) coordinates and determines the "closest"
    grid indices (i0, i1, i2) for the specified coordinate system. The C function is
    then registered for later use.

    :param CoordSystem: The coordinate system for the local grid patch.
    :param relative_to: Whether the computation is relative to the "local_grid_center"
                        (default) or "global_grid_center".
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :raises ValueError: When the value of `gridding_approach` is not "independent grid(s)"
                        or "multipatch".
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings, clang_format
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "Cart_to_xx_and_nearest_i0i1i2"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       _ = register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
    ...       generated_str = clang_format(cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function)
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}"
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[SinhSymTP]...
    Setting up reference_metric[HoleySinhSpherical]...
    Setting up reference_metric[Cartesian]...
    Setting up reference_metric[SinhCylindricalv2n2]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 1: Basic setup and parameter validation.
    if gridding_approach not in {"independent grid(s)", "multipatch"}:
        raise ValueError(
            "Invalid value for 'gridding_approach'. Must be 'independent grid(s)' or 'multipatch'."
        )

    parallelization = par.parval_from_str("parallelization")
    rfm = refmetric.reference_metric[CoordSystem]
    local_C_vars = {"xx0", "xx1", "xx2", "Cartx", "Carty", "Cartz"}

    namesuffix = f"_{relative_to}" if relative_to == "global_grid_center" else ""
    name = f"Cart_to_xx_and_nearest_i0i1i2{namesuffix}"
    params = "const params_struct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]"
    cfunc_decorators = "__host__ __device__" if parallelization == "cuda" else ""

    desc = "Given Cartesian point (x,y,z), this function "
    if gridding_approach == "multipatch":
        desc += "does stuff needed for multipatch, and then "
    desc += """unshifts the grid back to the origin to output the corresponding
            (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid"""

    # Step 2: Generate the core C-code for the coordinate transformation.
    core_body_list: List[str] = []
    if rfm.requires_NewtonRaphson_for_Cart_to_xx:
        # Part 2a: Handle mixed analytical and Newton-Raphson inversions.
        analytic_exprs: List[sp.Expr] = []
        analytic_names: List[str] = []
        for i in range(3):
            if rfm.NewtonRaphson_f_of_xx[i] == sp.sympify(0):
                analytic_exprs.append(rfm.Cart_to_xx[i])
                analytic_names.append(f"xx[{i}]")

        if analytic_exprs:
            core_body_list.append(
                "  // First compute analytical coordinate inversions:\n"
            )
            processed_analytic_exprs = _prepare_sympy_exprs_for_codegen(
                analytic_exprs, local_C_vars
            )
            core_body_list.append(
                ccg.c_codegen(
                    processed_analytic_exprs, analytic_names, include_braces=False
                )
            )

        core_body_list.append(
            """
  // Next perform Newton-Raphson iterations as needed:
  const REAL XX_TOLERANCE = 1e-12;  // that's 1 part in 1e12 dxxi.
  const REAL F_OF_XX_TOLERANCE = 1e-12;  // tolerance of function for which we're finding the root.
  const int ITER_MAX = 100;
  int iter, tolerance_has_been_met=0;
"""
        )
        for i in range(3):
            if rfm.NewtonRaphson_f_of_xx[i] != sp.sympify(0):
                nr_input_exprs = [
                    rfm.NewtonRaphson_f_of_xx[i],
                    sp.diff(rfm.NewtonRaphson_f_of_xx[i], rfm.xx[i]),
                ]
                nr_processed_exprs = _prepare_sympy_exprs_for_codegen(
                    nr_input_exprs, local_C_vars
                )
                nr_codegen_output = ccg.c_codegen(
                    nr_processed_exprs,
                    [f"f_of_xx{i}", f"fprime_of_xx{i}"],
                    include_braces=True,
                    verbose=False,
                )
                core_body_list.append(
                    f"""
  iter=0;
  REAL xx{i}  = 0.5 * (params->xxmin{i} + params->xxmax{i});
  while(iter < ITER_MAX && !tolerance_has_been_met) {{
    REAL f_of_xx{i}, fprime_of_xx{i};

{nr_codegen_output}
    const REAL xx{i}_np1 = xx{i} - f_of_xx{i} / fprime_of_xx{i};

    if( fabs(xx{i} - xx{i}_np1) <= XX_TOLERANCE * params->dxx{i} && fabs(f_of_xx{i}) <= F_OF_XX_TOLERANCE ) {{
      tolerance_has_been_met = 1;
    }}
    xx{i} = xx{i}_np1;
    iter++;
  }} // END Newton-Raphson iterations to compute xx{i}
  if(iter >= ITER_MAX) {{
    fprintf(stderr, "ERROR: Newton-Raphson failed for {CoordSystem}: xx{i}, x,y,z = %.15e %.15e %.15e\\n", Cartx,Carty,Cartz);
    exit(1);
  }}
  xx[{i}] = xx{i};
"""
                )
    else:
        # Part 2b: Handle purely analytical inversions.
        processed_exprs = _prepare_sympy_exprs_for_codegen(
            list(rfm.Cart_to_xx), local_C_vars
        )
        core_body_list.append(
            ccg.c_codegen(
                processed_exprs,
                ["xx[0]", "xx[1]", "xx[2]"],
                include_braces=False,
            )
        )

    # Part 2c: Add logic to find the nearest grid point index.
    core_body_list.append(
        """
      // Find the nearest grid indices (i0, i1, i2) for the given Cartesian coordinates (x, y, z).
      // Assuming a cell-centered grid, which follows the pattern:
      //   xx0[i0] = params->xxmin0 + ((REAL)(i0 - NGHOSTS) + 0.5) * params->dxx0
      // The index i0 can be derived as:
      //   i0 = (xx0[i0] - params->xxmin0) / params->dxx0 - 0.5 + NGHOSTS
      // Now, including typecasts:
      //   i0 = (int)((xx[0] - params->xxmin0) / params->dxx0 - 0.5 + (REAL)NGHOSTS)
      // The (int) typecast always rounds down, so we add 0.5 inside the outer parenthesis:
      //   i0 = (int)((xx[0] - params->xxmin0) / params->dxx0 - 0.5 + (REAL)NGHOSTS + 0.5)
      // The 0.5 values cancel out:
      //   i0 =           (int)( ( xx[0] - params->xxmin0 ) / params->dxx0 + (REAL)NGHOSTS )
      Cart_to_i0i1i2[0] = (int)( ( xx[0] - params->xxmin0 ) / params->dxx0 + (REAL)NGHOSTS );
      Cart_to_i0i1i2[1] = (int)( ( xx[1] - params->xxmin1 ) / params->dxx1 + (REAL)NGHOSTS );
      Cart_to_i0i1i2[2] = (int)( ( xx[2] - params->xxmin2 ) / params->dxx2 + (REAL)NGHOSTS );
"""
    )
    core_body = "".join(core_body_list)

    # Step 3: Assemble the full C function body.
    body_parts = [
        """
  // Set (Cartx, Carty, Cartz) relative to the global (as opposed to local) grid.
  //   This local grid may be offset from the origin by adjusting
  //   (Cart_originx, Cart_originy, Cart_originz) to nonzero values.
  REAL Cartx = xCart[0];
  REAL Carty = xCart[1];
  REAL Cartz = xCart[2];
"""
    ]
    if relative_to == "local_grid_center":
        body_parts.append(
            """
  // Set the origin, (Cartx, Carty, Cartz) = (0, 0, 0), to the center of the local grid patch.
  Cartx -= params->Cart_originx;
  Carty -= params->Cart_originy;
  Cartz -= params->Cart_originz;
  {
"""
        )
        body_parts.append(core_body)
        body_parts.append("  }\n")
    else:
        body_parts.append(core_body)

    # Step 4: Register the C function.
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc=desc,
        cfunc_type="void",
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body="".join(body_parts),
        cfunc_decorators=cfunc_decorators,
    )
    return pcg.NRPyEnv()


def register_CFunction_xx_to_Cart(
    CoordSystem: str,
    gridding_approach: str = "independent grid(s)",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Convert uniform-grid coordinate (xx[0], xx[1], xx[2]) to the corresponding Cartesian coordinate.

    :param CoordSystem: The coordinate system name as a string.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :raises ValueError: If an invalid gridding_approach is provided.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings, clang_format
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "xx_to_Cart"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       _ = register_CFunction_xx_to_Cart(CoordSystem)
    ...       generated_str = clang_format(cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function)
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}"
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 1: Basic setup and parameter validation.
    if gridding_approach not in {"independent grid(s)", "multipatch"}:
        raise ValueError(
            "Invalid value for 'gridding_approach'. Must be 'independent grid(s)' or 'multipatch'."
        )

    parallelization = par.parval_from_str("parallelization")
    rfm = refmetric.reference_metric[CoordSystem]
    local_C_vars = {"xx0", "xx1", "xx2"}

    # Step 2: Prepare SymPy expressions for C code generation.
    raw_xx_to_Cart_exprs = [rfm.xx_to_Cart[i] + gri.Cart_origin[i] for i in range(3)]
    processed_exprs = _prepare_sympy_exprs_for_codegen(
        raw_xx_to_Cart_exprs, local_C_vars
    )

    # Step 3: Construct the full C function body.
    codegen_results = ccg.c_codegen(
        processed_exprs,
        ["xCart[0]", "xCart[1]", "xCart[2]"],
    )
    body = (
        """
const REAL xx0 = xx[0];
const REAL xx1 = xx[1];
const REAL xx2 = xx[2];
"""
        + codegen_results
    )

    # Step 4: Register the C function.
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc="""Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} in terms of
local grid coordinates {xx[0][i0], xx[1][i1], xx[2][i2]} = {xx0, xx1, xx2},
taking into account the possibility that the origin of this grid is off-center.""",
        cfunc_type="void",
        CoordSystem_for_wrapper_func=CoordSystem,
        name="xx_to_Cart",
        params="const params_struct *restrict params, const REAL xx[3], REAL xCart[3]",
        include_CodeParameters_h=False,
        body=body,
        cfunc_decorators="__host__ __device__" if parallelization == "cuda" else "",
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
