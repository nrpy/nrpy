"""
Generate the derivative-dependent initial conformal connection for Dendro.

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
import nrpy.reference_metric as refmetric
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.infrastructures.Dendro import state_h
from nrpy.infrastructures.Dendro.simple_loop import simple_loop


def register_CFunction_initial_data_lambdaU(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
    fd_order: int = 6,
    CoordSystem: str = "Cartesian",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register one finite-difference order of the initial connection kernel.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :param enable_fCCZ4: Select the fCCZ4 state layout when true.
    :param fd_order: Centered finite-difference order.
    :param CoordSystem: Reference-metric coordinate system.
    :return: The NRPy registries, or ``None`` during parallel collection.
    :raises ValueError: If the order or Dendro generation settings are invalid.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("initial_data_lambdaU requires Infrastructure='Dendro'.")
    if fd_order not in (4, 6, 8):
        raise ValueError("initial_data_lambdaU supports FD orders 4, 6, and 8.")
    if par.parval_from_str("parallelization") != "none":
        raise ValueError("Dendro point kernels require parallelization='none'.")

    old_fd_order = par.parval_from_str("fd_order")
    par.set_parval_from_str("fd_order", fd_order)
    try:
        state_h.validate_registered_state(enable_fCCZ4)
        quantities = BSSN_quantities[CoordSystem]
        reference_metric = refmetric.reference_metric[CoordSystem]
        expressions = [
            quantities.DGammaU[i] / reference_metric.ReU[i] for i in range(3)
        ]
        scalar_type = gri.DENDRO_SCALAR_TYPE
        kernel = c_codegen(
            expressions,
            [f"out_lambdaU{i}[pp]" for i in range(3)],
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
        for index, name in enumerate(evolved_names):
            bindings.append(
                f"const {scalar_type}* {gri.DendroGridFunction.input_pointer(name)} = "
                f"in_gfs[{index}] + offset;"
            )
        for i in range(3):
            index = evolved_names.index(f"lambdaU{i}")
            bindings.append(
                f"{scalar_type}* out_lambdaU{i} = out_gfs[{index}] + offset;"
            )
        body = "\n".join(
            (
                "const std::ptrdiff_t offset = "
                "static_cast<std::ptrdiff_t>(block.getOffset());",
                "const unsigned nx_block = block.getAllocationSzX();",
                "const unsigned ny_block = block.getAllocationSzY();",
                "const unsigned nz_block = block.getAllocationSzZ();",
                "const unsigned padding_block = block.get1DPadWidth();",
                f"if (padding_block < {fd_order // 2}) {{",
                '    throw std::invalid_argument("initial_data_lambdaU padding '
                f'is too small for FD{fd_order}");',
                "}  // END IF: block padding too small",
                f"const {scalar_type} dx_block[3] = {{",
                "    block.computeDx(domain_min, domain_max),",
                "    block.computeDy(domain_min, domain_max),",
                "    block.computeDz(domain_min, domain_max)};",
                f"const {scalar_type} pmin_block[3] = {{",
                "    GRIDX_TO_X(block.getBlockNode().minX()) - "
                "padding_block * dx_block[0],",
                "    GRIDY_TO_Y(block.getBlockNode().minY()) - "
                "padding_block * dx_block[1],",
                "    GRIDZ_TO_Z(block.getBlockNode().minZ()) - "
                "padding_block * dx_block[2]};",
                *bindings,
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
            subdirectory="generated/src/initial_data_lambdaU",
            includes=[f"{solver_stem}_defines.h"],
            desc=f"Initialize covariant conformal connection with FD{fd_order}.",
            cfunc_type="void",
            name=f"initial_data_lambdaU_order_{fd_order}",
            params=(
                f"const ot::Block& block, const {scalar_type}* const* in_gfs, "
                f"{scalar_type}* const* out_gfs, const Point& domain_min, "
                "const Point& domain_max"
            ),
            body=body,
        )
    finally:
        par.set_parval_from_str("fd_order", old_fd_order)
    return pcg.NRPyEnv()
