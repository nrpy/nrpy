"""
Generate excised, proper-volume constraint reductions for Dendro GR.

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


def register_CFunction_diagnostics(
    solver_stem: str, *, enable_fCCZ4: bool = False
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register one-block constraint accumulation outside spherical excisions.

    The caller zeroes the output accumulators before visiting blocks and performs
    MPI reductions afterward. Squared norms and volume use the physical measure
    ``sqrt(gamma) d^3x = W^-6 d^3x``; maxima are pointwise absolute values.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :param enable_fCCZ4: Select fCCZ4 diagnostics when true.
    :return: The NRPy registries, or ``None`` during parallel collection.
    :raises ValueError: If the infrastructure, W formulation, or state is invalid.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("diagnostics requires Infrastructure='Dendro'.")
    if par.parval_from_str("EvolvedConformalFactor_cf") != "W":
        raise ValueError("Proper-volume diagnostics require W evolution.")

    state_h.validate_registered_state(enable_fCCZ4)
    evolved_names = tuple(state_h.evolved_gridfunctions(enable_fCCZ4))
    conformal_factor_index = evolved_names.index("cf")
    diagnostics_names = (
        state_h.FCCZ4_DIAGNOSTIC_GRIDFUNCTIONS
        if enable_fCCZ4
        else state_h.BSSN_DIAGNOSTIC_GRIDFUNCTIONS
    )
    scalar_type = gri.DENDRO_SCALAR_TYPE
    body = f"""if (diagnostic_gfs == nullptr || state_gfs == nullptr ||
    local_squared_norms == nullptr || local_max_norms == nullptr ||
    local_volume == nullptr)
    throw std::invalid_argument("diagnostics received a null array");
if (num_excision_regions > 0 &&
    (excision_centers == nullptr || excision_radii == nullptr))
    throw std::invalid_argument("diagnostics received null excision data");
const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(block.getOffset());
const unsigned nx = block.getAllocationSzX();
const unsigned ny = block.getAllocationSzY();
const unsigned nz = block.getAllocationSzZ();
const unsigned padding = block.get1DPadWidth();
if (nx <= 2 * padding || ny <= 2 * padding || nz <= 2 * padding)
    throw std::invalid_argument("diagnostics received invalid block padding");
const {scalar_type} dx[3] = {{block.computeDx(domain_min, domain_max),
                             block.computeDy(domain_min, domain_max),
                             block.computeDz(domain_min, domain_max)}};
const {scalar_type} pmin[3] = {{
    GRIDX_TO_X(block.getBlockNode().minX()) - padding * dx[0],
    GRIDY_TO_Y(block.getBlockNode().minY()) - padding * dx[1],
    GRIDZ_TO_Z(block.getBlockNode().minZ()) - padding * dx[2]}};
for (unsigned k = padding; k < nz - padding; ++k) {{
    const {scalar_type} z = pmin[2] + k * dx[2];
    for (unsigned j = padding; j < ny - padding; ++j) {{
        const {scalar_type} y = pmin[1] + j * dx[1];
        for (unsigned i = padding; i < nx - padding; ++i) {{
            const {scalar_type} x = pmin[0] + i * dx[0];
            const {scalar_type} quadrature_weight =
                (i == padding || i + 1 == nx - padding ? 0.5 : 1.0) *
                (j == padding || j + 1 == ny - padding ? 0.5 : 1.0) *
                (k == padding || k + 1 == nz - padding ? 0.5 : 1.0);
            bool excised = false;
            for (unsigned region = 0; region < num_excision_regions; ++region) {{
                const {scalar_type} delta_x = x - excision_centers[region].x();
                const {scalar_type} delta_y = y - excision_centers[region].y();
                const {scalar_type} delta_z = z - excision_centers[region].z();
                const {scalar_type} radius = excision_radii[region];
                excised = excised || delta_x * delta_x + delta_y * delta_y +
                                      delta_z * delta_z < radius * radius;
            }}
            if (excised) continue;
            const std::size_t pp = offset + i + static_cast<std::size_t>(nx) *
                                               (j + static_cast<std::size_t>(ny) * k);
            const {scalar_type} W = state_gfs[{conformal_factor_index}][pp];
            if (!(W > 0.0) || !std::isfinite(W)) {{
                *local_volume = std::numeric_limits<{scalar_type}>::infinity();
                for (unsigned field = 0; field < {len(diagnostics_names)}; ++field) {{
                    local_squared_norms[field] =
                        std::numeric_limits<{scalar_type}>::infinity();
                    local_max_norms[field] =
                        std::numeric_limits<{scalar_type}>::infinity();
                }}
                continue;
            }}
            const {scalar_type} proper_volume = quadrature_weight * dx[0] * dx[1] *
                                                dx[2] / std::pow(W, 6);
            *local_volume += proper_volume;
            for (unsigned field = 0; field < {len(diagnostics_names)}; ++field) {{
                const {scalar_type} value = diagnostic_gfs[field][pp];
                if (!std::isfinite(value)) {{
                    local_squared_norms[field] =
                        std::numeric_limits<{scalar_type}>::infinity();
                    local_max_norms[field] =
                        std::numeric_limits<{scalar_type}>::infinity();
                    continue;
                }}
                local_squared_norms[field] += value * value * proper_volume;
                local_max_norms[field] = std::max(local_max_norms[field], std::abs(value));
            }}
        }}
    }}
}}"""
    cfc.register_CFunction(
        subdirectory="generated/src/diagnostics",
        includes=[f"{solver_stem}_defines.h", "<limits>"],
        desc="Accumulate excised physical-volume constraint norms on one block.",
        cfunc_type="void",
        name="diagnostics",
        params=(
            f"const ot::Block& block, const {scalar_type}* const* diagnostic_gfs, "
            f"const {scalar_type}* const* state_gfs, "
            "const Point& domain_min, const Point& domain_max, "
            f"const Point* excision_centers, const {scalar_type}* excision_radii, "
            f"unsigned num_excision_regions, {scalar_type}* local_squared_norms, "
            f"{scalar_type}* local_max_norms, {scalar_type}* local_volume"
        ),
        body=body,
    )
    return pcg.NRPyEnv()
