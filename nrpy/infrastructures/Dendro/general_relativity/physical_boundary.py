"""
Generate radiative outer-boundary right-hand sides for Dendro GR fields.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, Union, cast

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.infrastructures.Dendro import state_h


def register_CFunction_physical_boundary(
    solver_stem: str, *, enable_fCCZ4: bool = False
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the outgoing-radiation condition on physical block faces.

    The kernel reuses the current padded state. It adds no synchronization and
    evaluates the selected centered derivative at physical boundary nodes.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :param enable_fCCZ4: Select fCCZ4 state layout when true.
    :return: The NRPy registries, or ``None`` during parallel collection.
    :raises ValueError: If the infrastructure or registered state is invalid.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("physical_boundary requires Infrastructure='Dendro'.")

    state_h.validate_registered_state(enable_fCCZ4)
    evolved_names = tuple(state_h.evolved_gridfunctions(enable_fCCZ4))
    falloff: Dict[str, str] = {}
    asymptotic: Dict[str, str] = {}
    for name in evolved_names:
        falloff[name] = "2.0" if name.startswith(("lambdaU", "aDD")) else "1.0"
        asymptotic[name] = "1.0" if name in ("alpha", "cf") else "0.0"

    scalar_type = gri.DENDRO_SCALAR_TYPE
    body = f"""const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(block.getOffset());
const unsigned nx = block.getAllocationSzX();
const unsigned ny = block.getAllocationSzY();
const unsigned nz = block.getAllocationSzZ();
const unsigned padding = block.get1DPadWidth();
const unsigned bflag = block.getBlkNodeFlag();
if (padding == 0 || nx <= 2 * padding || ny <= 2 * padding || nz <= 2 * padding) {{
    throw std::invalid_argument("physical_boundary received invalid block padding");
}}
const {scalar_type} dx[3] = {{block.computeDx(domain_min, domain_max),
                             block.computeDy(domain_min, domain_max),
                             block.computeDz(domain_min, domain_max)}};
const {scalar_type} pmin[3] = {{
    GRIDX_TO_X(block.getBlockNode().minX()) - padding * dx[0],
    GRIDY_TO_Y(block.getBlockNode().minY()) - padding * dx[1],
    GRIDZ_TO_Z(block.getBlockNode().minZ()) - padding * dx[2]}};
const unsigned derivative_radius = padding;
if (derivative_radius < 2 || derivative_radius > 4)
    throw std::invalid_argument("physical_boundary received an unsupported stencil");
const {scalar_type} derivative_coefficients[5][4] = {{
    {{0.0, 0.0, 0.0, 0.0}},
    {{0.0, 0.0, 0.0, 0.0}},
    {{2.0 / 3.0, -1.0 / 12.0, 0.0, 0.0}},
    {{3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0, 0.0}},
    {{4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0}}}};
for (unsigned k = padding; k < nz - padding; ++k) {{
    const {scalar_type} z = pmin[2] + k * dx[2];
    for (unsigned j = padding; j < ny - padding; ++j) {{
        const {scalar_type} y = pmin[1] + j * dx[1];
        for (unsigned i = padding; i < nx - padding; ++i) {{
            const bool on_boundary =
                ((bflag & (1u << OCT_DIR_LEFT)) && i == padding) ||
                ((bflag & (1u << OCT_DIR_RIGHT)) && i + 1 == nx - padding) ||
                ((bflag & (1u << OCT_DIR_DOWN)) && j == padding) ||
                ((bflag & (1u << OCT_DIR_UP)) && j + 1 == ny - padding) ||
                ((bflag & (1u << OCT_DIR_BACK)) && k == padding) ||
                ((bflag & (1u << OCT_DIR_FRONT)) && k + 1 == nz - padding);
            if (!on_boundary) continue;
            const {scalar_type} x = pmin[0] + i * dx[0];
            const {scalar_type} radius = std::sqrt(x * x + y * y + z * z);
            if (!(radius > 0.0))
                throw std::runtime_error("physical boundary intersects the origin");
            const std::size_t pp = i + static_cast<std::size_t>(nx) *
                                           (j + static_cast<std::size_t>(ny) * k);
            for (unsigned field = 0; field < {len(evolved_names)}; ++field) {{
                {scalar_type} gradient[3] = {{0.0, 0.0, 0.0}};
                for (unsigned distance = 1; distance <= derivative_radius;
                     ++distance) {{
                    const {scalar_type} coefficient =
                        derivative_coefficients[derivative_radius][distance - 1];
                    gradient[0] += coefficient *
                        (in_gfs[field][offset + pp + distance] -
                         in_gfs[field][offset + pp - distance]) / dx[0];
                    gradient[1] += coefficient *
                        (in_gfs[field][offset + pp + distance * nx] -
                         in_gfs[field][offset + pp - distance * nx]) / dx[1];
                    gradient[2] += coefficient *
                        (in_gfs[field][offset + pp + distance * nx * ny] -
                         in_gfs[field][offset + pp - distance * nx * ny]) / dx[2];
                }}
                const {scalar_type} radial_derivative =
                    x * gradient[0] + y * gradient[1] + z * gradient[2];
                rhs_gfs[field][offset + pp] =
                    -(radial_derivative + falloff[field] *
                      (in_gfs[field][offset + pp] - asymptotic[field])) / radius;
            }}
        }}
    }}
}}"""
    falloff_values = ", ".join(falloff[name] for name in evolved_names)
    asymptotic_values = ", ".join(asymptotic[name] for name in evolved_names)
    body = f"""const {scalar_type} falloff[{len(evolved_names)}] = {{{falloff_values}}};
const {scalar_type} asymptotic[{len(evolved_names)}] = {{{asymptotic_values}}};
{body}"""
    cfc.register_CFunction(
        subdirectory="generated/src/physical_boundary",
        includes=[f"{solver_stem}_defines.h"],
        desc="Apply state-aware outgoing-radiation RHS data on physical faces.",
        cfunc_type="void",
        name="physical_boundary",
        params=(
            f"const ot::Block& block, const {scalar_type}* const* in_gfs, "
            f"{scalar_type}* const* rhs_gfs, "
            "const Point& domain_min, const Point& domain_max"
        ),
        body=body,
    )
    return pcg.NRPyEnv()
