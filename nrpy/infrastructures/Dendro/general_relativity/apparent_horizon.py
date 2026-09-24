"""
Generate the module-owned BHaHAHA apparent-horizon dispatch for Dendro.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


def register_CFunction_apparent_horizon(
    solver_stem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register dispatch from the evolved BSSN state to BHaHAHA.

    The indices selected when constructing ``AEH_BHaHAHA`` name ``cf``,
    ``trK``, ``aDD00..22``, and ``hDD00..22``. BHaHAHA interpolates these
    evolved fields first, then the callback converts each sample to ADM data.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :return: The NRPy registries, or ``None`` during parallel collection.
    :raises ValueError: If the selected infrastructure is not Dendro.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("apparent_horizon requires Infrastructure='Dendro'.")

    cfc.register_CFunction(
        subdirectory="generated/src/apparent_horizon",
        includes=[f"{solver_stem}_defines.h", "aeh_bhahaha.h"],
        desc="Run BHaHAHA on evolved BSSN fields with post-interpolation ADM conversion.",
        cfunc_type="void",
        name="apparent_horizon",
        params=(
            "dendro_aeh::AEH_BHaHAHA& finder, const ot::Mesh* mesh, "
            "const double** evolved_gfs, unsigned iteration, double time, "
            "const std::vector<Point>& tracked_locations"
        ),
        body="""if (mesh == nullptr)
    throw std::invalid_argument("apparent_horizon received a null mesh");
if (evolved_gfs == nullptr)
    throw std::invalid_argument("apparent_horizon received null evolved fields");
finder.find_horizons(mesh, evolved_gfs, iteration, time, tracked_locations);""",
    )
    return pcg.NRPyEnv()
