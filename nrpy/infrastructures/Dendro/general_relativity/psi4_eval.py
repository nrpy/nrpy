"""
Generate Dendro Psi4 evaluation kernels.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.psi4 import Psi4
from nrpy.finite_difference import stencil_reach_per_axis
from nrpy.infrastructures.Dendro import state_h
from nrpy.infrastructures.Dendro.simple_loop import simple_loop


def register_CFunction_psi4_eval(
    solver_stem: str,
    *,
    CoordSystem: str = "Cartesian",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register FD4/6/8 Psi4 evaluation kernels.

    Each point kernel evaluates Psi4 from BSSN fields with the
    Baker-Campanelli-Lousto tetrad.

    :param solver_stem: Lowercase solver name used by the generated header.
    :param CoordSystem: Reference-metric coordinate system; must be Cartesian.
    :return: Updated NRPy environment, or ``None`` during task collection.
    :raises ValueError: If formulation or coordinates are invalid.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("Infrastructure must be 'Dendro' for wave extraction.")
    if par.parval_from_str("parallelization") != "none":
        raise ValueError("Dendro point kernels require parallelization='none'.")
    conformal_factor = par.parval_from_str("EvolvedConformalFactor_cf")
    if conformal_factor not in ("W", "chi"):
        raise ValueError("Dendro wave extraction requires W or chi.")
    if CoordSystem != "Cartesian":
        raise ValueError("Dendro wave extraction requires Cartesian coordinates.")
    missing_fields = [
        name
        for name in state_h.BSSN_EVOLVED_GRIDFUNCTIONS
        if name not in gri.glb_gridfcs_dict
    ]
    if missing_fields:
        raise ValueError(
            "Psi4 evaluation requires canonical BSSN fields; missing "
            + ", ".join(missing_fields)
        )

    scalar_type = gri.DENDRO_SCALAR_TYPE
    psi4 = Psi4(CoordSystem=CoordSystem, enable_rfm_precompute=False)
    contraction = c_codegen(
        [psi4.psi4_re, psi4.psi4_im],
        ["psi4_real[pp]", "psi4_imag[pp]"],
        include_braces=False,
        enable_fd_codegen=False,
        enable_simd=False,
        fp_type=str(par.parval_from_str("fp_type")),
        fp_type_alias=scalar_type,
        cse_sorting="none",
        verbose=False,
    )
    input_bindings: List[str] = []
    for index, name in enumerate(state_h.BSSN_EVOLVED_GRIDFUNCTIONS):
        gridfunction = cast(gri.DendroGridFunction, gri.glb_gridfcs_dict[name])
        input_bindings.append(
            f"const {scalar_type}* in_{gridfunction.dendro_name} = "
            f"in_gfs[{index}] + offset;"
        )
    point_state_lines: List[str] = []
    for name in (
        "cf",
        "trK",
        "hDD00",
        "hDD01",
        "hDD02",
        "hDD11",
        "hDD12",
        "hDD22",
        "aDD00",
        "aDD01",
        "aDD02",
        "aDD11",
        "aDD12",
        "aDD22",
    ):
        gridfunction = cast(gri.DendroGridFunction, gri.glb_gridfcs_dict[name])
        point_state_lines.append(
            f"const {scalar_type} {name} = in_{gridfunction.dendro_name}[pp];"
        )
    point_state = "\n".join(point_state_lines)
    derivative_unpack = "\n".join(
        f"const {scalar_type} {name} = {array_name};"
        for name, array_name in zip(
            psi4.metric_derivs_varname_list,
            psi4.metric_derivs_varname_arr_list,
        )
    )
    inverse_chi_denominator = "cf * cf" if conformal_factor == "W" else "cf"
    tetrad = f"""const {scalar_type} inverse_W2 = 1.0 / ({inverse_chi_denominator});
const {scalar_type} gammaDD00 = (1.0 + hDD00) * inverse_W2;
const {scalar_type} gammaDD01 = hDD01 * inverse_W2;
const {scalar_type} gammaDD02 = hDD02 * inverse_W2;
const {scalar_type} gammaDD11 = (1.0 + hDD11) * inverse_W2;
const {scalar_type} gammaDD12 = hDD12 * inverse_W2;
const {scalar_type} gammaDD22 = (1.0 + hDD22) * inverse_W2;
const {scalar_type} detgamma =
    gammaDD00 * (gammaDD11 * gammaDD22 - gammaDD12 * gammaDD12)
    - gammaDD01 * (gammaDD01 * gammaDD22 - gammaDD02 * gammaDD12)
    + gammaDD02 * (gammaDD01 * gammaDD12 - gammaDD02 * gammaDD11);
const {scalar_type} inverse_detgamma = 1.0 / detgamma;
const {scalar_type} gammaUU00 =
    (gammaDD11 * gammaDD22 - gammaDD12 * gammaDD12) * inverse_detgamma;
const {scalar_type} gammaUU01 =
    (gammaDD02 * gammaDD12 - gammaDD01 * gammaDD22) * inverse_detgamma;
const {scalar_type} gammaUU02 =
    (gammaDD01 * gammaDD12 - gammaDD02 * gammaDD11) * inverse_detgamma;
const {scalar_type} gammaUU11 =
    (gammaDD00 * gammaDD22 - gammaDD02 * gammaDD02) * inverse_detgamma;
const {scalar_type} gammaUU12 =
    (gammaDD01 * gammaDD02 - gammaDD00 * gammaDD12) * inverse_detgamma;
const {scalar_type} gammaUU22 =
    (gammaDD00 * gammaDD11 - gammaDD01 * gammaDD01) * inverse_detgamma;
{scalar_type} v1U[3] = {{-xx1, xx0, 0.0}};
{scalar_type} v2U[3] = {{xx0, xx1, xx2}};
const {scalar_type} radius2 = xx0 * xx0 + xx1 * xx1 + xx2 * xx2;
const {scalar_type} axis2 = xx0 * xx0 + xx1 * xx1;
const {scalar_type} seed_threshold =
    64.0 * std::numeric_limits<{scalar_type}>::epsilon();
if (axis2 <= seed_threshold * std::max(radius2, static_cast<{scalar_type}>(1))) {{
    v1U[0] = 1.0;
    v1U[1] = 0.0;
    if (radius2 <= seed_threshold) {{
        v2U[0] = 0.0;
        v2U[1] = 1.0;
        v2U[2] = 0.0;
    }}
}}
const {scalar_type} crossD[3] = {{
    v1U[1] * v2U[2] - v1U[2] * v2U[1],
    v1U[2] * v2U[0] - v1U[0] * v2U[2],
    v1U[0] * v2U[1] - v1U[1] * v2U[0]}};
const {scalar_type} sqrt_detgamma = std::sqrt(detgamma);
const {scalar_type} v3U[3] = {{
    sqrt_detgamma *
        (gammaUU00 * crossD[0] + gammaUU01 * crossD[1] + gammaUU02 * crossD[2]),
    sqrt_detgamma *
        (gammaUU01 * crossD[0] + gammaUU11 * crossD[1] + gammaUU12 * crossD[2]),
    sqrt_detgamma *
        (gammaUU02 * crossD[0] + gammaUU12 * crossD[1] + gammaUU22 * crossD[2])}};
const auto metric_inner = [&](const {scalar_type}* left,
                              const {scalar_type}* right) {{
    return left[0] * (gammaDD00 * right[0] + gammaDD01 * right[1]
                      + gammaDD02 * right[2])
         + left[1] * (gammaDD01 * right[0] + gammaDD11 * right[1]
                      + gammaDD12 * right[2])
         + left[2] * (gammaDD02 * right[0] + gammaDD12 * right[1]
                      + gammaDD22 * right[2]);
}};
const {scalar_type} inverse_norm_v1 = 1.0 / std::sqrt(metric_inner(v1U, v1U));
{scalar_type} e1U[3] = {{v1U[0] * inverse_norm_v1,
                         v1U[1] * inverse_norm_v1,
                         v1U[2] * inverse_norm_v1}};
const {scalar_type} e1_dot_v2 = metric_inner(e1U, v2U);
{scalar_type} e2U[3] = {{v2U[0] - e1_dot_v2 * e1U[0],
                         v2U[1] - e1_dot_v2 * e1U[1],
                         v2U[2] - e1_dot_v2 * e1U[2]}};
const {scalar_type} inverse_norm_e2 = 1.0 / std::sqrt(metric_inner(e2U, e2U));
for (unsigned component = 0; component < 3; ++component) {{
    e2U[component] *= inverse_norm_e2;
}}
const {scalar_type} e1_dot_v3 = metric_inner(e1U, v3U);
const {scalar_type} e2_dot_v3 = metric_inner(e2U, v3U);
{scalar_type} e3U[3] = {{
    v3U[0] - e1_dot_v3 * e1U[0] - e2_dot_v3 * e2U[0],
    v3U[1] - e1_dot_v3 * e1U[1] - e2_dot_v3 * e2U[1],
    v3U[2] - e1_dot_v3 * e1U[2] - e2_dot_v3 * e2U[2]}};
const {scalar_type} inverse_norm_e3 = 1.0 / std::sqrt(metric_inner(e3U, e3U));
for (unsigned component = 0; component < 3; ++component) {{
    e3U[component] *= inverse_norm_e3;
}}
const {scalar_type} n4U0 = M_SQRT1_2;
const {scalar_type} n4U1 = -M_SQRT1_2 * e2U[0];
const {scalar_type} n4U2 = -M_SQRT1_2 * e2U[1];
const {scalar_type} n4U3 = -M_SQRT1_2 * e2U[2];
const {scalar_type} mre4U0 = 0.0;
const {scalar_type} mre4U1 = M_SQRT1_2 * e3U[0];
const {scalar_type} mre4U2 = M_SQRT1_2 * e3U[1];
const {scalar_type} mre4U3 = M_SQRT1_2 * e3U[2];
const {scalar_type} mim4U0 = 0.0;
const {scalar_type} mim4U1 = M_SQRT1_2 * e1U[0];
const {scalar_type} mim4U2 = M_SQRT1_2 * e1U[1];
const {scalar_type} mim4U3 = M_SQRT1_2 * e1U[2];"""

    old_fd_order = par.parval_from_str("fd_order")
    try:
        for fd_order in (4, 6, 8):
            par.set_parval_from_str("fd_order", fd_order)
            padding = max(
                stencil_reach_per_axis(
                    psi4.metric_derivs_expr_list,
                    "unset",
                    fd_order,
                )
            )
            if padding != fd_order // 2:
                raise ValueError(
                    f"Psi4 FD{fd_order} requires padding {padding}, "
                    f"expected {fd_order // 2}."
                )
            metric_derivatives = c_codegen(
                psi4.metric_derivs_expr_list,
                psi4.metric_derivs_varname_arr_list,
                include_braces=False,
                enable_fd_codegen=True,
                enable_fd_functions=False,
                enable_simd=False,
                fp_type=str(par.parval_from_str("fp_type")),
                fp_type_alias=scalar_type,
                mem_alloc_style="210",
                rational_const_alias="static const",
                cse_sorting="none",
                verbose=False,
                upwind_control_vec=sp.Symbol("unset"),
            )
            geometry = f"""const std::ptrdiff_t offset =
    static_cast<std::ptrdiff_t>(block.getOffset());
const unsigned nx_block = block.getAllocationSzX();
const unsigned ny_block = block.getAllocationSzY();
const unsigned nz_block = block.getAllocationSzZ();
const unsigned padding_block = block.get1DPadWidth();
if (padding_block < {padding}) {{
    throw std::invalid_argument(
        "psi4_eval block padding is too small for FD{fd_order}");
}}
const {scalar_type} dx_block[3] = {{
    block.computeDx(domain_min, domain_max),
    block.computeDy(domain_min, domain_max),
    block.computeDz(domain_min, domain_max)}};
const {scalar_type} pmin_block[3] = {{
    GRIDX_TO_X(block.getBlockNode().minX()) - padding_block * dx_block[0],
    GRIDY_TO_Y(block.getBlockNode().minY()) - padding_block * dx_block[1],
    GRIDZ_TO_Z(block.getBlockNode().minZ()) - padding_block * dx_block[2]}};"""
            point_evaluation = f"""{scalar_type} arr_gammaDDdDD[81] = {{}};
{scalar_type} arr_GammaUDD[27] = {{}};
{scalar_type} arr_KDDdD[27] = {{}};
{{
{metric_derivatives}
}}
{point_state}
{derivative_unpack}
{tetrad}
{contraction}"""
            body = "\n".join(
                (
                    geometry,
                    *input_bindings,
                    f"{scalar_type}* psi4_real = psi4_gfs[0] + offset;",
                    f"{scalar_type}* psi4_imag = psi4_gfs[1] + offset;",
                    simple_loop(
                        point_evaluation,
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
                subdirectory="generated/src/psi4_eval",
                includes=[
                    f"{solver_stem}_defines.h",
                    "<algorithm>",
                    "<cmath>",
                    "<limits>",
                    "<stdexcept>",
                ],
                desc=f"Evaluate Psi4 on one padded block at FD order {fd_order}.",
                cfunc_type="void",
                name=f"psi4_eval_order_{fd_order}",
                params=(
                    f"const ot::Block& block, const {scalar_type}* const* in_gfs, "
                    f"{scalar_type}* const* psi4_gfs, const Point& domain_min, "
                    "const Point& domain_max"
                ),
                body=body,
            )
    finally:
        par.set_parval_from_str("fd_order", old_fd_order)
    return pcg.NRPyEnv()
