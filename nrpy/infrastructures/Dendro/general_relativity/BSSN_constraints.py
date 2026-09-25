"""
Generate BSSN constraint diagnostics for Dendro.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.finite_difference import stencil_reach_per_axis
from nrpy.infrastructures.Dendro import state_h
from nrpy.infrastructures.Dendro.simple_loop import simple_loop


def register_CFunction_BSSN_constraints(
    solver_stem: str,
    *,
    fd_order: int = 6,
    CoordSystem: str = "Cartesian",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register one order-specific per-block BSSN constraint kernel.

    :param solver_stem: Lowercase solver name used by the generated header.
    :param fd_order: Centered finite-difference order.
    :param CoordSystem: Reference-metric coordinate system.
    :return: Updated NRPy environment, or ``None`` during task collection.
    :raises ValueError: If configuration, layout, or padding is invalid.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("Infrastructure must be 'Dendro' to build BSSN constraints.")
    if par.parval_from_str("parallelization") != "none":
        raise ValueError("Dendro point kernels require parallelization='none'.")
    if fd_order not in (4, 6, 8):
        raise ValueError(f"Unsupported fd_order={fd_order!r}; allowed: (4, 6, 8).")

    old_fd_order = par.parval_from_str("fd_order")
    par.set_parval_from_str("fd_order", fd_order)
    try:
        for name in state_h.BSSN_DIAGNOSTIC_GRIDFUNCTIONS:
            if name not in gri.glb_gridfcs_dict:
                gri.register_gridfunctions(name, group="DIAG", is_basename=False)
        old_register_magnitudes = par.parval_from_str(
            "register_M_and_LAMBDA_CONSTRAINT_gridfunctions"
        )
        old_register_momentum = par.parval_from_str("register_MU_gridfunctions")
        par.set_parval_from_str("register_M_and_LAMBDA_CONSTRAINT_gridfunctions", False)
        par.set_parval_from_str("register_MU_gridfunctions", False)
        try:
            constraints = BSSN_constraints[CoordSystem]
        finally:
            par.set_parval_from_str(
                "register_M_and_LAMBDA_CONSTRAINT_gridfunctions",
                old_register_magnitudes,
            )
            par.set_parval_from_str("register_MU_gridfunctions", old_register_momentum)

        bssn_quantities = BSSN_quantities[CoordSystem]
        momentum_covariant = [
            sum(bssn_quantities.gammabarDD[i][j] * constraints.MU[j] for j in range(3))
            / bssn_quantities.exp_m4phi
            for i in range(3)
        ]
        expressions_by_gridfunction: Dict[str, sp.Expr] = {
            "H": constraints.H,
            "MU0": momentum_covariant[0],
            "MU1": momentum_covariant[1],
            "MU2": momentum_covariant[2],
            "M_CONSTRAINT": sp.sqrt(constraints.Msquared),
            "LAMBDA_CONSTRAINT": constraints.LambdaConstraintMagnitude,
        }
        if tuple(expressions_by_gridfunction) != state_h.BSSN_DIAGNOSTIC_GRIDFUNCTIONS:
            raise ValueError("BSSN diagnostic expressions are not in canonical order.")
        expressions = list(expressions_by_gridfunction.values())
        padding = max(stencil_reach_per_axis(expressions, "unset", fd_order))
        if padding != fd_order // 2:
            raise ValueError(
                f"BSSN constraints FD{fd_order} require padding {padding}, "
                f"expected {fd_order // 2}."
            )

        scalar_type = gri.DENDRO_SCALAR_TYPE
        kernel = c_codegen(
            expressions,
            [f"diagnostic_{name}[pp]" for name in expressions_by_gridfunction],
            include_braces=False,
            enable_fd_codegen=True,
            enable_fd_functions=False,
            enable_simd=False,
            fp_type=str(par.parval_from_str("fp_type")),
            fp_type_alias=scalar_type,
            mem_alloc_style="210",
            rational_const_alias="static const",
            cse_sorting="none",
            verbose=False,
            upwind_control_vec=sp.Symbol("unset"),
        )
        input_bindings = []
        for index, name in enumerate(state_h.BSSN_EVOLVED_GRIDFUNCTIONS):
            input_bindings.append(
                f"const {scalar_type}* {gri.DendroGridFunction.input_pointer(name)} = "
                f"in_gfs[{index}] + offset;"
            )
        diagnostic_bindings = [
            f"{scalar_type}* diagnostic_{name} = diagnostic_gfs[{index}] + offset;"
            for index, name in enumerate(state_h.BSSN_DIAGNOSTIC_GRIDFUNCTIONS)
        ]
        geometry = f"""const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(block.getOffset());
const unsigned nx_block = block.getAllocationSzX();
const unsigned ny_block = block.getAllocationSzY();
const unsigned nz_block = block.getAllocationSzZ();
const unsigned padding_block = block.get1DPadWidth();
if (padding_block < {fd_order // 2}) {{
    throw std::invalid_argument("BSSN_constraints block padding is too small for FD{fd_order}");
}}  // END IF: block padding too small
const {scalar_type} dx_block[3] = {{
    block.computeDx(domain_min, domain_max),
    block.computeDy(domain_min, domain_max),
    block.computeDz(domain_min, domain_max)}};
const {scalar_type} pmin_block[3] = {{
    GRIDX_TO_X(block.getBlockNode().minX()) - padding_block * dx_block[0],
    GRIDY_TO_Y(block.getBlockNode().minY()) - padding_block * dx_block[1],
    GRIDZ_TO_Z(block.getBlockNode().minZ()) - padding_block * dx_block[2]}};"""
        body = "\n".join(
            (
                geometry,
                *input_bindings,
                *diagnostic_bindings,
                simple_loop(
                    kernel,
                    nx="nx_block",
                    ny="ny_block",
                    nz="nz_block",
                    padding="padding_block",
                    pmin_padded="pmin_block",
                    dx="dx_block",
                ),
            )
        )
        cfc.register_CFunction(
            subdirectory="generated/src/BSSN_constraints",
            includes=[f"{solver_stem}_defines.h"],
            desc=f"Per-block BSSN constraints at FD order {fd_order}.",
            cfunc_type="void",
            name=f"BSSN_constraints_order_{fd_order}",
            params=(
                f"const ot::Block& block, const {scalar_type}* const* in_gfs, "
                f"{scalar_type}* const* diagnostic_gfs, const Point& domain_min, "
                "const Point& domain_max"
            ),
            body=body,
        )
    finally:
        par.set_parval_from_str("fd_order", old_fd_order)
    return pcg.NRPyEnv()
