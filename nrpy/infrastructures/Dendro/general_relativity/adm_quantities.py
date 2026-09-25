"""
Generate finite-radius ADM surface quadrature for Dendro.

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


def register_CFunction_adm_quantities(
    solver_stem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register distributed finite-radius ADM quadrature.

    The caller provides the interpolation ownership indices and owns the final
    MPI reduction.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :return: The NRPy registries, or ``None`` during parallel collection.
    :raises ValueError: If the selected infrastructure is not Dendro.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("adm_quantities requires Infrastructure='Dendro'.")

    scalar_type = gri.DENDRO_SCALAR_TYPE
    quadrature_body = f"""if (local_quantities == nullptr)
    throw std::invalid_argument("adm_quantities received a null accumulator");
if (num_valid_points == 0) return;
if (valid_point_indices == nullptr || surface_coordinates == nullptr ||
    surface_normals == nullptr || surface_weights == nullptr ||
    surface_data == nullptr || center == nullptr)
    throw std::invalid_argument("adm_quantities received a null surface array");
for (unsigned field = 0; field < 9; ++field) {{
    if (surface_data[field] == nullptr)
        throw std::invalid_argument("adm_quantities received a null data field");
}}  // END LOOP: for field over surface data
constexpr {scalar_type} inverse_8pi =
    0.039788735772973833942220940843128590508;
for (unsigned valid = 0; valid < num_valid_points; ++valid) {{
    const unsigned point = valid_point_indices[valid];
    if (point >= num_points)
        throw std::out_of_range("ADM quadrature point index is out of range");
    const {scalar_type}* normal = surface_normals + 3 * point;
    const {scalar_type}* coordinate = surface_coordinates + 3 * point;
    const {scalar_type} weighted_inverse_8pi =
        surface_weights[point] * inverse_8pi;
    const {scalar_type} mass_integrand =
        surface_data[0][point] * normal[0] +
        surface_data[1][point] * normal[1] +
        surface_data[2][point] * normal[2];
    local_quantities[0] +=
        static_cast<{scalar_type}>(0.5) * weighted_inverse_8pi * mass_integrand;
    const {scalar_type} momentum_surface[3] = {{
        surface_data[3][point] * normal[0] +
            surface_data[4][point] * normal[1] +
            surface_data[5][point] * normal[2],
        surface_data[4][point] * normal[0] +
            surface_data[6][point] * normal[1] +
            surface_data[7][point] * normal[2],
        surface_data[5][point] * normal[0] +
            surface_data[7][point] * normal[1] +
            surface_data[8][point] * normal[2]}};
    for (unsigned direction = 0; direction < 3; ++direction)
        local_quantities[1 + direction] +=
            weighted_inverse_8pi * momentum_surface[direction];
    const {scalar_type} relative_coordinate[3] = {{
        coordinate[0] - center[0],
        coordinate[1] - center[1],
        coordinate[2] - center[2]}};
    local_quantities[4] += weighted_inverse_8pi *
        (relative_coordinate[1] * momentum_surface[2] -
         relative_coordinate[2] * momentum_surface[1]);
    local_quantities[5] += weighted_inverse_8pi *
        (relative_coordinate[2] * momentum_surface[0] -
         relative_coordinate[0] * momentum_surface[2]);
    local_quantities[6] += weighted_inverse_8pi *
        (relative_coordinate[0] * momentum_surface[1] -
         relative_coordinate[1] * momentum_surface[0]);
}}  // END LOOP: for valid over owned points"""
    cfc.register_CFunction(
        subdirectory="generated/src/adm_quantities",
        includes=[f"{solver_stem}_defines.h"],
        desc="Accumulate rank-local ADM mass and momenta from owned surface samples.",
        cfunc_type="void",
        name="adm_quantities",
        params=(
            "unsigned num_points, unsigned num_valid_points, "
            "const unsigned* valid_point_indices, "
            f"const {scalar_type}* surface_coordinates, "
            f"const {scalar_type}* surface_normals, "
            f"const {scalar_type}* surface_weights, "
            f"const {scalar_type}* const* surface_data, "
            f"const {scalar_type} center[3], "
            f"{scalar_type} local_quantities[7]"
        ),
        body=quadrature_body,
    )
    return pcg.NRPyEnv()
