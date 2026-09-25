"""
Generate the separate conformal-Ricci kernel used by Dendro GR RHSs.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.finite_difference import stencil_reach_per_axis
from nrpy.infrastructures.Dendro import state_h
from nrpy.infrastructures.Dendro.simple_loop import simple_loop


def register_CFunction_Ricci_eval(
    solver_stem: str,
    *,
    fd_order: int = 6,
    CoordSystem: str = "Cartesian",
    enable_intrinsics: bool = True,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register one order-specific per-block conformal-Ricci kernel.

    :param solver_stem: Generated application namespace and file-name stem.
    :param fd_order: Centered finite-difference order.
    :param CoordSystem: NRPy reference-metric coordinate system.
    :param enable_intrinsics: Generate SIMD-intrinsic kernels; every vector stays
        inside its own row.
    :return: The parallel-codegen environment outside the registration phase.
    :raises ValueError: If the infrastructure, profile, or Ricci layout is invalid.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("Infrastructure must be 'Dendro' to build Ricci_eval.")
    if fd_order not in (4, 6, 8):
        raise ValueError(f"Unsupported fd_order={fd_order!r}; allowed: (4, 6, 8).")
    if par.parval_from_str("parallelization") != "none":
        raise ValueError("Dendro point kernels require parallelization='none'.")
    if enable_intrinsics and CoordSystem != "Cartesian":
        raise ValueError("Dendro SIMD Ricci kernels require Cartesian coordinates.")

    old_fd_order = par.parval_from_str("fd_order")
    par.set_parval_from_str("fd_order", fd_order)
    try:
        quantities = BSSN_quantities[CoordSystem]
        _ = BSSN_quantities[CoordSystem + "_RbarDD_gridfunctions"]
        if tuple(quantities.Ricci_varnames) != state_h.RICCI_GRIDFUNCTIONS:
            raise ValueError(
                "Ricci_eval must produce exactly the six canonical RbarDD components."
            )
        expressions = list(quantities.Ricci_exprs)
        padding = max(stencil_reach_per_axis(expressions, "unset", fd_order))
        if padding != fd_order // 2:
            raise ValueError(
                f"Ricci_eval FD{fd_order} requires padding {padding}, "
                f"expected {fd_order // 2}."
            )
        scalar_type = gri.DENDRO_SCALAR_TYPE
        kernel = c_codegen(
            expressions,
            [f"ricci_{name}[pp]" for name in state_h.RICCI_GRIDFUNCTIONS],
            enable_simd=enable_intrinsics,
            enable_fd_codegen=True,
            enable_fd_functions=False,
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
            dendro_name = cast(
                gri.DendroGridFunction, gri.glb_gridfcs_dict[name]
            ).dendro_name
            input_bindings.append(
                f"const {scalar_type}* in_{dendro_name} = in_gfs[{index}] + offset;"
            )
        ricci_bindings = [
            f"{scalar_type}* ricci_{name} = ricci_gfs[{index}] + offset;"
            for index, name in enumerate(state_h.RICCI_GRIDFUNCTIONS)
        ]
        geometry = f"""const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(block.getOffset());
const unsigned nx_block = block.getAllocationSzX();
const unsigned ny_block = block.getAllocationSzY();
const unsigned nz_block = block.getAllocationSzZ();
const unsigned padding_block = block.get1DPadWidth();
if (padding_block < {fd_order // 2}) {{
    throw std::invalid_argument("Ricci_eval block padding is too small for FD{fd_order}");
}}  // END IF: block padding too small
const {scalar_type} dx_block[3] = {{
    block.computeDx(domain_min, domain_max),
    block.computeDy(domain_min, domain_max),
    block.computeDz(domain_min, domain_max)}};
[[maybe_unused]] const {scalar_type} pmin_block[3] = {{
    GRIDX_TO_X(block.getBlockNode().minX()) - padding_block * dx_block[0],
    GRIDY_TO_Y(block.getBlockNode().minY()) - padding_block * dx_block[1],
    GRIDZ_TO_Z(block.getBlockNode().minZ()) - padding_block * dx_block[2]}};"""
        body = "\n".join(
            (
                geometry,
                *input_bindings,
                *ricci_bindings,
                simple_loop(
                    kernel,
                    nx="nx_block",
                    ny="ny_block",
                    nz="nz_block",
                    padding="padding_block",
                    pmin_padded="pmin_block",
                    dx="dx_block",
                    enable_intrinsics=enable_intrinsics,
                ),
            )
        )
        cfc.register_CFunction(
            subdirectory="generated/src/Ricci_eval",
            # simd_intrinsics.h must precede the definitions header, which can
            # reach a copy with the same include guard; "./" keeps this order
            # after clang-format sorts the includes.
            includes=(["./simd_intrinsics.h"] if enable_intrinsics else [])
            + [f"{solver_stem}_defines.h"],
            desc=f"Per-block conformal Ricci tensor at FD order {fd_order}.",
            cfunc_type="void",
            name=f"Ricci_eval_order_{fd_order}",
            params=(
                f"const ot::Block& block, const {scalar_type}* const* in_gfs, "
                f"{scalar_type}* const* ricci_gfs, const Point& domain_min, "
                "const Point& domain_max"
            ),
            body=body,
        )
    finally:
        par.set_parval_from_str("fd_order", old_fd_order)
    return pcg.NRPyEnv()
