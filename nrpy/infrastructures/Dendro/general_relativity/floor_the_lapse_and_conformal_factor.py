"""
Generate the Dendro lapse and conformal-factor floor.

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
from nrpy.infrastructures.Dendro import state_h
from nrpy.infrastructures.Dendro.simple_loop import simple_loop


def register_CFunction_floor_the_lapse_and_conformal_factor(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the per-block lapse and conformal-factor floor.

    Dendro's ``CHI_FLOOR`` is applied directly to the lapse.  Since these
    applications evolve ``W = sqrt(chi)``, the conformal-factor floor is
    ``sqrt(CHI_FLOOR)``.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :param enable_fCCZ4: Select the fCCZ4 state layout when true.
    :return: Updated NRPy registries, or ``None`` during task collection.
    :raises ValueError: If the Dendro generation settings are invalid.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("Field floors require Infrastructure='Dendro'.")
    if par.parval_from_str("parallelization") != "none":
        raise ValueError("Dendro point kernels require parallelization='none'.")

    state_h.validate_registered_state(enable_fCCZ4)
    par.register_CodeParameter(
        "REAL",
        __name__,
        "chi_floor",
        1.0e-4,
        assumption="RealPositive",
        commondata=True,
        add_to_parfile=True,
        description="Floor for chi and alpha; W is floored at sqrt(chi_floor).",
    )
    scalar_type = gri.DENDRO_SCALAR_TYPE
    evolved_names = tuple(state_h.evolved_gridfunctions(enable_fCCZ4))
    alpha_index = evolved_names.index("alpha")
    conformal_factor_index = evolved_names.index("cf")
    loop_body = "\n".join(
        (
            "alpha[pp] = std::max(alpha[pp], chi_floor);",
            "cf_W_or_chi[pp] = std::max(cf_W_or_chi[pp], std::sqrt(chi_floor));",
        )
    )
    body = "\n".join(
        (
            "if (!(chi_floor > 0.0) || !std::isfinite(chi_floor))",
            '    throw std::invalid_argument("CHI_FLOOR must be finite and positive");',
            "const std::ptrdiff_t offset = "
            "static_cast<std::ptrdiff_t>(block.getOffset());",
            "const unsigned nx_block = block.getAllocationSzX();",
            "const unsigned ny_block = block.getAllocationSzY();",
            "const unsigned nz_block = block.getAllocationSzZ();",
            f"const {scalar_type} dx_block[3] = {{",
            "    block.computeDx(domain_min, domain_max),",
            "    block.computeDy(domain_min, domain_max),",
            "    block.computeDz(domain_min, domain_max)};",
            f"const {scalar_type} pmin_block[3] = {{",
            "    GRIDX_TO_X(block.getBlockNode().minX()),",
            "    GRIDY_TO_Y(block.getBlockNode().minY()),",
            "    GRIDZ_TO_Z(block.getBlockNode().minZ())};",
            f"{scalar_type}* alpha = in_gfs[{alpha_index}] + offset;",
            f"{scalar_type}* cf_W_or_chi = in_gfs[{conformal_factor_index}] + offset;",
            simple_loop(
                loop_body,
                nx="nx_block",
                ny="ny_block",
                nz="nz_block",
                padding="0",
                pmin_padded="pmin_block",
                dx="dx_block",
            ),
        )
    )
    cfc.register_CFunction(
        subdirectory="generated/src/floor_the_lapse_and_conformal_factor",
        includes=[f"{solver_stem}_defines.h", "<algorithm>", "<cmath>", "<stdexcept>"],
        desc="Floor alpha and W consistently with Dendro's CHI_FLOOR.",
        cfunc_type="void",
        name="floor_the_lapse_and_conformal_factor",
        params=(
            f"const ot::Block& block, {scalar_type}* const* in_gfs, "
            f"const {scalar_type} chi_floor, const Point& domain_min, "
            "const Point& domain_max"
        ),
        body=body,
        ET_current_thorn_CodeParams_used=["chi_floor"],
    )
    return pcg.NRPyEnv()
