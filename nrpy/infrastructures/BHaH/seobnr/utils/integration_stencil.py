"""
Register CFunction for 8-th order accurate integration stencils.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def get_integration_stencil_uniform_spacing(
    polynomial_order: int,
) -> tuple[list[list[sp.core.numbers.Rational]], list[list[int]]]:
    """
    Evaluate 8-th order accurate numerical integrals.
    This function evaluates the integrals
         /x_i
    I =  | f(x) dx
         /x_{i-1}
    for i = 1,2,...,polynomial_order. The stencil is constructed by
    finding the (polynomial_order + 1)-point Lagrange interpolating polynomial
    that passes through the equally-spaced points (x_i,f_i) and
    integrating it between x_{i-1} and x_i.

    :param polynomial_order: order of the interpolating polynomial.
    :return: tuple of lists of integration weights and offset indices.
    """
    x_0, x = sp.symbols("x_0 x", real=True)
    X = [x_0 + i for i in range(polynomial_order + 1)]
    Y = [sp.Symbol(f"y_{i}", real=True) for i in range(polynomial_order + 1)]
    interpolating_polynomial = 0
    for i in range(polynomial_order + 1):
        l_i = 1
        for j in range(polynomial_order + 1):
            if i != j:
                l_i *= (x - X[j]) / (X[i] - X[j])
        interpolating_polynomial += l_i * Y[i]
    weights: list[list[sp.core.numbers.Rational]] = []
    indices: list[list[int]] = []
    for i in range(polynomial_order):
        integral = sp.simplify(
            sp.integrate(interpolating_polynomial, (x, X[i], X[i + 1]))
        )
        weights_i = [integral.coeff(Y[j]) for j in range(polynomial_order + 1)]
        weights.append(weights_i)
        indices.append([j - i for j in range(polynomial_order + 1)])
    return (weights, indices)


def register_CFunction_integration_stencil() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for 8-th order accurate integration stencils.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
    Evaluate 8-th order accurate integration stencils.
    This function evaluates the integral
           /x_i
    I_i =  | f(x) dx
           /x_{i-1}
    when x and f(x) are given as arrays. The stencil is constructed by
    finding the 8-point Lagrange interpolating polynomial that passes through 
    the point (x_i,f_i) with arbitrary offset (default is centered at x_i)
    and performing integration between x_{i-1} and x_i.
    
    @param offset - offset of the point (determines the forwards/backwards/centralized interpolation).
    @param coeffs - Array of integration weights for the stencil.
    @param indices - Array of indices for the stencil.
    """
    cfunc_type = "void"
    name = "integration_stencil"
    params = "const int offset, REAL *restrict coeffs, int *restrict indices"

    interpolating_polynomial_order = 7
    ws, idxs = get_integration_stencil_uniform_spacing(interpolating_polynomial_order)
    offsets = [3, 2, 1, 0, -1, -2, -3]
    body = """
switch(offset){
"""
    for i, offset in enumerate(offsets):
        body += f"""
  case {offset}:
"""
        for j, coeff in enumerate(ws[i]):
            body += sp.ccode(coeff, assign_to=f"coeffs[{j}]")
        for j, index in enumerate(idxs[i]):
            body += f"indices[{j}] = {index};\n"
        body += """
    break;
"""
    body += """
    default:
        fprintf(stderr,"Error: Invalid offset %d\\nOffset must be in [-3,3]\\n", offset);
        exit(1);
}
"""
    cfc.register_CFunction(
        subdirectory="utils",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
