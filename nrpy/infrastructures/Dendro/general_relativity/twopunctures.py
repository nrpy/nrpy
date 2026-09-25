"""
Generate Dendro block interpolation from NRPy's TwoPunctures solution.

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
from nrpy.infrastructures.BHaH.general_relativity import (
    NRPyPN_quasicircular_momenta,
)
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import (
    TwoPunctures_lib,
)


def register_CFunction_twopunctures(
    solver_stem: str,
    *,
    tp_orientation: str = "native_cartesian_xy_plane",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register TwoPunctures and interpolate its solved data onto one Dendro block interior.

    The generated solver loads coefficients computed by its single-rank
    ``--tpid`` mode. This function then fills the eighteen ADM scratch fields
    consumed by ``ADM_to_BSSN``. The required TwoPunctures lapse is
    ``alpha=psi**(-2)=W``.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :param tp_orientation: Fixed-frame orientation used by TwoPunctures.
    :return: The NRPy registries, or ``None`` during parallel collection.
    :raises ValueError: If the infrastructure or orientation is unsupported.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("twopunctures requires Infrastructure='Dendro'.")
    if tp_orientation not in {"legacy_swap_xz", "native_cartesian_xy_plane"}:
        raise ValueError(f"Unsupported TwoPunctures orientation: {tp_orientation}.")

    par.set_parval_from_str("Infrastructure", "BHaH")
    TwoPunctures_lib.register_C_functions_explicit(tp_orientation)
    NRPyPN_quasicircular_momenta.register_CFunction_NRPyPN_quasicircular_momenta()
    par.set_parval_from_str("Infrastructure", "Dendro")
    scalar_type = gri.DENDRO_SCALAR_TYPE
    assignments = (
        "sample.gammaSphorCartDD00",
        "sample.gammaSphorCartDD01",
        "sample.gammaSphorCartDD02",
        "sample.gammaSphorCartDD11",
        "sample.gammaSphorCartDD12",
        "sample.gammaSphorCartDD22",
        "sample.KSphorCartDD00",
        "sample.KSphorCartDD01",
        "sample.KSphorCartDD02",
        "sample.KSphorCartDD11",
        "sample.KSphorCartDD12",
        "sample.KSphorCartDD22",
        "sample.betaSphorCartU0",
        "sample.betaSphorCartU1",
        "sample.betaSphorCartU2",
        "sample.BSphorCartU0",
        "sample.BSphorCartU1",
        "sample.BSphorCartU2",
    )
    copy_to_adm = "\n".join(
        f"                adm_gfs[{index}][pp] = {value};"
        for index, value in enumerate(assignments)
    )
    body = f"""if (commondata == nullptr || params == nullptr || punctures == nullptr ||
    adm_gfs == nullptr)
    throw std::invalid_argument("twopunctures received a null argument");
if (std::strcmp(punctures->initial_lapse, "W") != 0)
    throw std::invalid_argument("twopunctures requires alpha=W=(psi+u)^(-2)");
const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(block.getOffset());
const unsigned nx = block.getAllocationSzX();
const unsigned ny = block.getAllocationSzY();
const unsigned nz = block.getAllocationSzZ();
const unsigned padding = block.get1DPadWidth();
const {scalar_type} dx[3] = {{block.computeDx(domain_min, domain_max),
                             block.computeDy(domain_min, domain_max),
                             block.computeDz(domain_min, domain_max)}};
const {scalar_type} pmin[3] = {{
    GRIDX_TO_X(block.getBlockNode().minX()) - padding * dx[0],
    GRIDY_TO_Y(block.getBlockNode().minY()) - padding * dx[1],
    GRIDZ_TO_Z(block.getBlockNode().minZ()) - padding * dx[2]}};
// Interior points only: zip reads only these; the solver refills the padding
// by unzipping the converted evolved state and extrapolating into
// physical-boundary padding.
for (unsigned k = padding; k < nz - padding; ++k) {{
    const {scalar_type} z = pmin[2] + k * dx[2];
    for (unsigned j = padding; j < ny - padding; ++j) {{
        const {scalar_type} y = pmin[1] + j * dx[1];
        for (unsigned i = padding; i < nx - padding; ++i) {{
            const REAL x[3] = {{pmin[0] + i * dx[0], y, z}};
            initial_data_struct sample{{}};
            TP_Interp(commondata, params, x, punctures, &sample);
            const std::size_t pp = offset + i + static_cast<std::size_t>(nx) *
                                               (j + static_cast<std::size_t>(ny) * k);
{copy_to_adm}
            if (!(sample.alpha > 0.0) || !std::isfinite(sample.alpha))
                throw std::runtime_error("TwoPunctures returned invalid W lapse");
        }}  // END LOOP: for i over interior x
    }}  // END LOOP: for j over interior y
}}  // END LOOP: for k over interior z"""
    cfc.register_CFunction(
        subdirectory="generated/src/twopunctures",
        includes=[
            f"{solver_stem}_defines.h",
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
        ],
        desc="Interpolate solved TwoPunctures ADM data onto one Dendro block interior.",
        cfunc_type="void",
        name="twopunctures",
        params=(
            "const ot::Block& block, const commondata_struct* commondata, "
            "const params_struct* params, const ID_persist_struct* punctures, "
            f"{scalar_type}* const* adm_gfs, "
            "const Point& domain_min, const Point& domain_max"
        ),
        body=body,
    )
    return pcg.NRPyEnv()
