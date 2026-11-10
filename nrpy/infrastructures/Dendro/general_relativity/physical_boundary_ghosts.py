"""
Generate physical-boundary ghost values for Dendro GR block fields.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


def register_CFunction_physical_boundary_ghosts(
    solver_stem: str,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register extrapolation into physical exterior padding.

    Six interior points give fourth-order boundary second derivatives at
    FD6 and FD8. FD4 uses five points because a one-element block has only
    five interior nodes. The x, y, and z passes successively fill faces,
    edges, and corners without touching interior or inter-block ghost values.
    Call this on the current unzipped state before any centered derivative.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :return: The NRPy registries, or ``None`` during parallel collection.
    :raises ValueError: If the infrastructure is not Dendro.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("physical_boundary_ghosts requires Infrastructure='Dendro'.")

    desc = """Extrapolate current block fields into physical exterior padding.

@param[in] block Dendro block with physical-face flags.
@param[in,out] fields Unzipped fields whose physical ghosts are filled.
@param field_count Number of fields in the pointer table.
"""
    cfunc_type = "void"
    name = "physical_boundary_ghosts"
    params = "const ot::Block& block, DendroScalar* const* fields, unsigned field_count"
    body = r"""const unsigned bflag = block.getBlkNodeFlag();
if (bflag == 0) return;
const unsigned padding = block.get1DPadWidth();
const unsigned dimensions[3] = {block.getAllocationSzX(),
                                block.getAllocationSzY(),
                                block.getAllocationSzZ()};
const unsigned interior_nodes = padding == 2 ? 5 : 6;
if (padding < 2 || padding > 4 ||
    dimensions[0] < 2 * padding + interior_nodes ||
    dimensions[1] < 2 * padding + interior_nodes ||
    dimensions[2] < 2 * padding + interior_nodes)
    throw std::invalid_argument("physical_boundary_ghosts has too few interior nodes");
const std::size_t strides[3] = {
    1, dimensions[0], static_cast<std::size_t>(dimensions[0]) * dimensions[1]};
const unsigned face_flags[3][2] = {
    {OCT_DIR_LEFT, OCT_DIR_RIGHT},
    {OCT_DIR_DOWN, OCT_DIR_UP},
    {OCT_DIR_BACK, OCT_DIR_FRONT}};
const DendroScalar coefficients_fd4[3][6] = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {5.0, -10.0, 10.0, -5.0, 1.0, 0.0},
    {15.0, -40.0, 45.0, -24.0, 5.0, 0.0}};
const DendroScalar coefficients_fd6_fd8[5][6] = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {6.0, -15.0, 20.0, -15.0, 6.0, -1.0},
    {21.0, -70.0, 105.0, -84.0, 35.0, -6.0},
    {56.0, -210.0, 336.0, -280.0, 120.0, -21.0},
    {126.0, -504.0, 840.0, -720.0, 315.0, -56.0}};
const DendroScalar (*coefficients)[6] =
    padding == 2 ? coefficients_fd4 : coefficients_fd6_fd8;
const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(block.getOffset());
for (unsigned field = 0; field < field_count; ++field) {
    DendroScalar* const values = fields[field] + offset;
    for (unsigned axis = 0; axis < 3; ++axis) {
        unsigned lower[3] = {0, 0, 0};
        unsigned upper[3] = {dimensions[0], dimensions[1], dimensions[2]};
        lower[axis] = 0;
        upper[axis] = 1;
        for (unsigned tangent = axis + 1; tangent < 3; ++tangent) {
            lower[tangent] = (bflag & (1u << face_flags[tangent][0]))
                                 ? padding : 0;
            upper[tangent] = (bflag & (1u << face_flags[tangent][1]))
                                 ? dimensions[tangent] - padding
                                 : dimensions[tangent];
        } // END LOOP: for tangent over later axes
        for (unsigned side = 0; side < 2; ++side) {
            if (!(bflag & (1u << face_flags[axis][side]))) continue;
            const unsigned boundary = side == 0 ? padding :
                dimensions[axis] - padding - 1;
            for (unsigned k = lower[2]; k < upper[2]; ++k) {
                for (unsigned j = lower[1]; j < upper[1]; ++j) {
                    for (unsigned i = lower[0]; i < upper[0]; ++i) {
                        const std::size_t base = i + static_cast<std::size_t>(dimensions[0]) *
                            (j + static_cast<std::size_t>(dimensions[1]) * k);
                        for (unsigned distance = 1; distance <= padding; ++distance) {
                            DendroScalar extrapolated = 0.0;
                            for (unsigned inner = 0; inner < interior_nodes; ++inner) {
                                const unsigned source = side == 0 ?
                                    boundary + inner : boundary - inner;
                                extrapolated += coefficients[distance][inner] *
                                    values[base + source * strides[axis]];
                            } // END LOOP: for inner over interpolation nodes
                            const unsigned target = side == 0 ?
                                boundary - distance : boundary + distance;
                            values[base + target * strides[axis]] = extrapolated;
                        } // END LOOP: for distance over ghost points
                    } // END LOOP: for i over tangential points
                } // END LOOP: for j over tangential points
            } // END LOOP: for k over tangential points
        } // END LOOP: for side over physical faces
    } // END LOOP: for axis over Cartesian directions
} // END LOOP: for field over unzipped fields"""
    cfc.register_CFunction(
        subdirectory="generated/src/physical_boundary_ghosts",
        includes=[f"{solver_stem}_defines.h"],
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )
    return pcg.NRPyEnv()
