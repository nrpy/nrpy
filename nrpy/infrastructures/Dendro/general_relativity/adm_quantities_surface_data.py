"""
Generate finite-radius ADM surface-data kernels for Dendro.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.infrastructures.Dendro.simple_loop import simple_loop


def register_CFunction_adm_quantities_surface_data(
    solver_stem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register block-local ADM surface-data kernels.

    The three order-specific block kernels consume the first twelve fields
    emitted by ``BSSN_to_ADM`` (six physical-metric and six extrinsic-curvature
    components). They emit three ADM-mass flux components followed by the six
    symmetric components of ``K_ij - gamma_ij K``. The caller zips and
    interpolates these nine scalar fields with Dendrolib, then passes the
    interpolation ``validIndices`` to the separate ``adm_quantities`` kernel.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :return: The NRPy registries, or ``None`` during parallel collection.
    :raises ValueError: If the selected infrastructure is not Dendro.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "adm_quantities_surface_data requires Infrastructure='Dendro'."
        )

    scalar_type = gri.DENDRO_SCALAR_TYPE
    symmetric_components = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    centered_first_derivatives: Dict[int, Tuple[Tuple[int, int], ...]] = {
        4: ((-2, 1), (-1, -8), (1, 8), (2, -1)),
        6: ((-3, -1), (-2, 9), (-1, -45), (1, 45), (2, -9), (3, 1)),
        8: (
            (-4, 3),
            (-3, -32),
            (-2, 168),
            (-1, -672),
            (1, 672),
            (2, -168),
            (3, 32),
            (4, -3),
        ),
    }
    derivative_denominators = {4: 12, 6: 60, 8: 840}

    for fd_order in (4, 6, 8):
        input_bindings = [
            f"const {scalar_type}* gamma{first}{second} = adm_gfs[{index}] + offset;"
            for index, (first, second) in enumerate(symmetric_components)
        ]
        input_bindings.extend(
            f"const {scalar_type}* curvature{first}{second} = "
            f"adm_gfs[{6 + index}] + offset;"
            for index, (first, second) in enumerate(symmetric_components)
        )
        output_bindings = [
            f"{scalar_type}* mass_flux{direction} = "
            f"surface_gfs[{direction}] + offset;"
            for direction in range(3)
        ]
        output_bindings.extend(
            f"{scalar_type}* momentum{first}{second} = "
            f"surface_gfs[{3 + index}] + offset;"
            for index, (first, second) in enumerate(symmetric_components)
        )

        derivative_lines: List[str] = []
        for first, second in symmetric_components:
            for direction, stride in enumerate(("1", "nx", "nxy")):
                terms = []
                for offset, coefficient in centered_first_derivatives[fd_order]:
                    sign = "+" if coefficient > 0 else "-"
                    magnitude = abs(coefficient)
                    coefficient_text = "" if magnitude == 1 else f"{magnitude} * "
                    terms.append(
                        f" {sign} {coefficient_text}gamma{first}{second}"
                        f"[pp + ({offset}) * {stride}]"
                    )
                derivative_lines.append(
                    f"const {scalar_type} d{direction}_gamma{first}{second} = "
                    f"({''.join(terms).lstrip(' +')}) * invdxx{direction} / "
                    f"static_cast<{scalar_type}>({derivative_denominators[fd_order]});"
                )
        momentum_lines = [
            f"momentum{first}{second}[pp] = curvature{first}{second}[pp] - "
            f"gamma{first}{second}[pp] * trace_curvature;"
            for first, second in symmetric_components
        ]

        point_body = "\n".join(
            (
                *derivative_lines,
                f"const {scalar_type} determinant =",
                "    gamma00[pp] * (gamma11[pp] * gamma22[pp] - "
                "gamma12[pp] * gamma12[pp]) -",
                "    gamma01[pp] * (gamma01[pp] * gamma22[pp] - "
                "gamma02[pp] * gamma12[pp]) +",
                "    gamma02[pp] * (gamma01[pp] * gamma12[pp] - "
                "gamma02[pp] * gamma11[pp]);",
                "if (!(determinant > 0.0) || !std::isfinite(determinant)) {",
                f"    const {scalar_type} invalid = "
                f"std::numeric_limits<{scalar_type}>::quiet_NaN();",
                "    mass_flux0[pp] = invalid;",
                "    mass_flux1[pp] = invalid;",
                "    mass_flux2[pp] = invalid;",
                "    momentum00[pp] = invalid;",
                "    momentum01[pp] = invalid;",
                "    momentum02[pp] = invalid;",
                "    momentum11[pp] = invalid;",
                "    momentum12[pp] = invalid;",
                "    momentum22[pp] = invalid;",
                "    continue;",
                "}  // END IF: invalid metric determinant",
                f"const {scalar_type} inverse_gamma00 = "
                "(gamma11[pp] * gamma22[pp] - gamma12[pp] * gamma12[pp]) / "
                "determinant;",
                f"const {scalar_type} inverse_gamma01 = "
                "(gamma02[pp] * gamma12[pp] - gamma01[pp] * gamma22[pp]) / "
                "determinant;",
                f"const {scalar_type} inverse_gamma02 = "
                "(gamma01[pp] * gamma12[pp] - gamma02[pp] * gamma11[pp]) / "
                "determinant;",
                f"const {scalar_type} inverse_gamma11 = "
                "(gamma00[pp] * gamma22[pp] - gamma02[pp] * gamma02[pp]) / "
                "determinant;",
                f"const {scalar_type} inverse_gamma12 = "
                "(gamma01[pp] * gamma02[pp] - gamma00[pp] * gamma12[pp]) / "
                "determinant;",
                f"const {scalar_type} inverse_gamma22 = "
                "(gamma00[pp] * gamma11[pp] - gamma01[pp] * gamma01[pp]) / "
                "determinant;",
                f"const {scalar_type} trace_curvature =",
                "    inverse_gamma00 * curvature00[pp] +",
                "    2.0 * inverse_gamma01 * curvature01[pp] +",
                "    2.0 * inverse_gamma02 * curvature02[pp] +",
                "    inverse_gamma11 * curvature11[pp] +",
                "    2.0 * inverse_gamma12 * curvature12[pp] +",
                "    inverse_gamma22 * curvature22[pp];",
                "mass_flux0[pp] = d1_gamma01 + d2_gamma02 - d0_gamma11 - "
                "d0_gamma22;",
                "mass_flux1[pp] = d0_gamma01 + d2_gamma12 - d1_gamma00 - "
                "d1_gamma22;",
                "mass_flux2[pp] = d0_gamma02 + d1_gamma12 - d2_gamma00 - "
                "d2_gamma11;",
                *momentum_lines,
            )
        )
        geometry = f"""if (adm_gfs == nullptr || surface_gfs == nullptr) {{
    throw std::invalid_argument(
        "adm_quantities_surface_data_order_{fd_order} received a null field table");
}}  // END IF: null field table
for (unsigned field = 0; field < 12; ++field) {{
    if (adm_gfs[field] == nullptr)
        throw std::invalid_argument(
            "adm_quantities_surface_data_order_{fd_order} received a null ADM field");
}}  // END LOOP: for field over ADM inputs
for (unsigned field = 0; field < 9; ++field) {{
    if (surface_gfs[field] == nullptr)
        throw std::invalid_argument(
            "adm_quantities_surface_data_order_{fd_order} received a null output field");
}}  // END LOOP: for field over surface outputs
const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(block.getOffset());
const unsigned nx_block = block.getAllocationSzX();
const unsigned ny_block = block.getAllocationSzY();
const unsigned nz_block = block.getAllocationSzZ();
const unsigned padding_block = block.get1DPadWidth();
if (padding_block < {fd_order // 2}) {{
    throw std::invalid_argument(
        "ADM surface-data block padding is too small for FD{fd_order}");
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
                *output_bindings,
                simple_loop(
                    point_body,
                    nx="nx_block",
                    ny="ny_block",
                    nz="nz_block",
                    padding="padding_block",
                    pmin_padded="pmin_block",
                    dx="dx_block",
                ),
            )
        )
        desc = (
            "Compute block-local ADM mass flux and momentum surface tensor "
            f"at FD order {fd_order}."
        )
        cfunc_type = "void"
        name = f"adm_quantities_surface_data_order_{fd_order}"
        params = (
            f"const ot::Block& block, const {scalar_type}* const* adm_gfs, "
            f"{scalar_type}* const* surface_gfs, const Point& domain_min, "
            "const Point& domain_max"
        )
        cfc.register_CFunction(
            subdirectory="generated/src/adm_quantities_surface_data",
            includes=[f"{solver_stem}_defines.h", "<limits>"],
            desc=desc,
            cfunc_type=cfunc_type,
            name=name,
            params=params,
            body=body,
        )

    return pcg.NRPyEnv()
