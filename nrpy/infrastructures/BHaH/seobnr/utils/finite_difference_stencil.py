"""
Register CFunction for 8-th order accurate finite difference stencils.

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
from nrpy.infrastructures import BHaH


def register_CFunction_finite_difference_stencil() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for 8-th order accurate finite difference stencils.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
    Evaluate 8-th order accurate finite difference stencils.

    @param offset - array index (determines the forwards/backwards offset).
    @param coeffs - Array of coefficients for the stencil.
    @param indices - Array of indices for the stencil.
    """
    cfunc_type = "void"
    name = "finite_difference_stencil"
    params = "const int offset, REAL *restrict coeffs, int *restrict indices"

    FDORDER = 8
    offsets = [4, 3, 2, 1, 0, -1, -2, -3, -4]
    body = """
switch(offset){
"""
    for offset in offsets:
        coeffs, indices = (
            BHaH.CurviBoundaryConditions.apply_bcs_outerradiation_and_inner.get_arb_offset_FD_coeffs_indices(
                FDORDER, offset, 1
            )
        )
        body += f"""
  case {offset}:
"""
        for j, coeff in enumerate(coeffs):
            body += sp.ccode(coeff, assign_to=f"coeffs[{j}]")
        for j, index in enumerate(indices):
            body += f"indices[{j}] = {index};\n"
        body += """
    break;
"""
    body += """
    default:
        fprintf(stderr,"Error: Invalid offset %d\\nOffset must be in [-4,4]\\n", offset);
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
