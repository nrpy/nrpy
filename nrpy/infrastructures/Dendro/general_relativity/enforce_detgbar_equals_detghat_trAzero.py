"""
Generate the Dendro algebraic projection for the BSSN variables.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_algebraic_constraints import (
    BSSN_algebraic_constraints,
)
from nrpy.infrastructures.Dendro import state_h


def register_CFunction_enforce_detgbar_equals_detghat_trAzero(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
    CoordSystem: str = "Cartesian",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the owned-node determinant and trace-free projection.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :param enable_fCCZ4: Select the fCCZ4 state layout when true.
    :param CoordSystem: Reference-metric coordinate system.
    :return: Updated NRPy registries, or ``None`` during task collection.
    :raises ValueError: If the Dendro generation settings are invalid.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("Algebraic projection requires Infrastructure='Dendro'.")
    if par.parval_from_str("parallelization") != "none":
        raise ValueError("Dendro point kernels require parallelization='none'.")
    if CoordSystem != "Cartesian":
        raise ValueError("Dendro owned-node projection requires Cartesian coordinates.")

    state_h.validate_registered_state(enable_fCCZ4)
    hprimeDD, aprimeDD = BSSN_algebraic_constraints(CoordSystem, False)
    component_names = tuple(
        f"{tensor_name}{i}{j}"
        for tensor_name in ("hDD", "aDD")
        for i in range(3)
        for j in range(i, 3)
    )
    expressions = [
        tensor[i][j]
        for tensor in (hprimeDD, aprimeDD)
        for i in range(3)
        for j in range(i, 3)
    ]
    scalar_type = gri.DENDRO_SCALAR_TYPE
    kernel = c_codegen(
        expressions,
        [f"out_{name}[pp]" for name in component_names],
        include_braces=False,
        enable_fd_codegen=True,
        enable_fd_functions=False,
        enable_simd=False,
        fp_type=str(par.parval_from_str("fp_type")),
        fp_type_alias=scalar_type,
        cse_sorting="none",
        verbose=False,
    )
    evolved_names = tuple(state_h.evolved_gridfunctions(enable_fCCZ4))
    bindings = []
    for name in component_names:
        index = evolved_names.index(name)
        bindings.append(
            f"const {scalar_type}* {gri.DendroGridFunction.input_pointer(name)} = "
            f"in_gfs[{index}];"
        )
        bindings.append(f"{scalar_type}* out_{name} = in_gfs[{index}];")
    body = "\n".join(
        (
            *bindings,
            "for (unsigned pp = node_begin; pp < node_end; ++pp) {",
            kernel,
            "}  // END LOOP: for pp over node range",
        )
    )
    cfc.register_CFunction(
        subdirectory="generated/src/enforce_detgbar_equals_detghat_trAzero",
        includes=[f"{solver_stem}_defines.h"],
        desc="Enforce det(gammabar)=det(gammahat) and tr(Abar)=0 per owned node.",
        cfunc_type="void",
        name="enforce_detgbar_equals_detghat_trAzero",
        params=f"{scalar_type}* const* in_gfs, unsigned node_begin, unsigned node_end",
        body=body,
    )
    return pcg.NRPyEnv()
