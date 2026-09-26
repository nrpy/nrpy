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


def register_CFunction_floor_the_lapse_and_conformal_factor(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the owned-node lapse and conformal-factor floor.

    Dendro's ``CHI_FLOOR`` is applied directly to the lapse.  Since these
    applications may evolve ``W = sqrt(chi)`` or ``chi``, and the
    conformal-factor floor follows the evolved variable.

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
    conformal_factor = par.parval_from_str("EvolvedConformalFactor_cf")
    if conformal_factor not in ("W", "chi"):
        raise ValueError("Dendro BSSN and fCCZ4 require W or chi.")

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
            "cf_W_or_chi[pp] = std::max(cf_W_or_chi[pp], "
            + ("std::sqrt(chi_floor)" if conformal_factor == "W" else "chi_floor")
            + ");",
        )
    )
    body = "\n".join(
        (
            "if (!(chi_floor > 0.0) || !std::isfinite(chi_floor))",
            '    throw std::invalid_argument("CHI_FLOOR must be finite and positive");',
            f"{scalar_type}* alpha = in_gfs[{alpha_index}];",
            f"{scalar_type}* cf_W_or_chi = in_gfs[{conformal_factor_index}];",
            "for (unsigned pp = node_begin; pp < node_end; ++pp) {",
            loop_body,
            "}  // END LOOP: for pp over node range",
        )
    )
    desc = "Floor alpha and the conformal factor with Dendro's CHI_FLOOR."
    cfunc_type = "void"
    name = "floor_the_lapse_and_conformal_factor"
    params = (
        f"{scalar_type}* const* in_gfs, unsigned node_begin, "
        f"unsigned node_end, const {scalar_type} chi_floor"
    )
    cfc.register_CFunction(
        subdirectory="generated/src/floor_the_lapse_and_conformal_factor",
        includes=[f"{solver_stem}_defines.h", "<algorithm>", "<cmath>", "<stdexcept>"],
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        ET_current_thorn_CodeParams_used=["chi_floor"],
    )
    return pcg.NRPyEnv()
