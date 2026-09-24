"""
Generate the pointwise BSSN-to-ADM conversion for Dendro.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.infrastructures.Dendro import state_h
from nrpy.infrastructures.Dendro.simple_loop import simple_loop


def register_CFunction_BSSN_to_ADM(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
    CoordSystem: str = "Cartesian",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register conversion of the evolved formulation to Cartesian ADM fields.

    The output order is ``gammaDD``, ``KDD``, ``betaU``, then ``BU``. This is
    the exact inverse storage contract consumed by ``ADM_to_BSSN``. ``BU`` is
    the Gamma-driver auxiliary field, not an ADM variable; retaining it makes
    the scratch representation lossless for restart and service calls.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :param enable_fCCZ4: Validate the fCCZ4 state instead of the BSSN state.
    :param CoordSystem: Reference-metric coordinate system.
    :return: The NRPy registries, or ``None`` during parallel collection.
    :raises ValueError: If generation is not configured for Cartesian data.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("BSSN_to_ADM requires Infrastructure='Dendro'.")
    conformal_factor = par.parval_from_str("EvolvedConformalFactor_cf")
    if conformal_factor not in ("W", "chi"):
        raise ValueError("Dendro BSSN and fCCZ4 require W or chi.")
    if par.parval_from_str("parallelization") != "none":
        raise ValueError("Dendro point kernels require parallelization='none'.")
    if CoordSystem != "Cartesian":
        raise ValueError("Dendro ADM service fields use the Cartesian basis.")

    state_h.validate_registered_state(enable_fCCZ4)
    converted = BSSN_to_ADM(CoordSystem)
    expressions: Dict[str, sp.Expr] = {}
    for tensor_name, tensor in (("gammaDD", converted.gammaDD), ("KDD", converted.KDD)):
        for i in range(3):
            for j in range(i, 3):
                expressions[f"{tensor_name}{i}{j}"] = tensor[i][j]
    for i in range(3):
        expressions[f"betaU{i}"] = converted.betaU[i]
        expressions[f"BU{i}"] = sp.Symbol(f"betU{i}", real=True)

    output_names = tuple(expressions)
    evolved_names = tuple(state_h.evolved_gridfunctions(enable_fCCZ4))
    scalar_type = gri.DENDRO_SCALAR_TYPE
    input_bindings: List[str] = []
    local_inputs: List[str] = []
    for index, name in enumerate(evolved_names):
        gridfunction = gri.glb_gridfcs_dict[name]
        if not isinstance(gridfunction, gri.DendroGridFunction):
            raise ValueError(f"{name} is not registered as a Dendro gridfunction.")
        input_bindings.append(
            f"const {scalar_type}* in_{gridfunction.dendro_name} = "
            f"in_gfs[{index}] + offset;"
        )
        if name != "Theta_fCCZ4" and not name.startswith("lambdaU"):
            local_inputs.append(
                f"const {scalar_type} {name} = in_{gridfunction.dendro_name}[pp];"
            )
    output_bindings = [
        f"{scalar_type}* adm_{name} = adm_gfs[{index}] + offset;"
        for index, name in enumerate(output_names)
    ]
    kernel = c_codegen(
        [expressions[name] for name in output_names],
        [f"adm_{name}[pp]" for name in output_names],
        include_braces=False,
        enable_simd=False,
        fp_type=str(par.parval_from_str("fp_type")),
        fp_type_alias=scalar_type,
        cse_sorting="none",
        verbose=False,
    )
    body = "\n".join(
        (
            "const std::ptrdiff_t offset = "
            "static_cast<std::ptrdiff_t>(block.getOffset());",
            "const unsigned nx_block = block.getAllocationSzX();",
            "const unsigned ny_block = block.getAllocationSzY();",
            "const unsigned nz_block = block.getAllocationSzZ();",
            "const DendroScalar dx_block[3] = {block.computeDx(domain_min, domain_max), "
            "block.computeDy(domain_min, domain_max), block.computeDz(domain_min, domain_max)};",
            "const unsigned padding_block = block.get1DPadWidth();",
            "const DendroScalar pmin_block[3] = {"
            "GRIDX_TO_X(block.getBlockNode().minX()) - padding_block * dx_block[0], "
            "GRIDY_TO_Y(block.getBlockNode().minY()) - padding_block * dx_block[1], "
            "GRIDZ_TO_Z(block.getBlockNode().minZ()) - padding_block * dx_block[2]};",
            *input_bindings,
            *output_bindings,
            simple_loop(
                "\n".join((*local_inputs, kernel)),
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
        subdirectory="generated/src/BSSN_to_ADM",
        includes=[f"{solver_stem}_defines.h"],
        desc="Pointwise conversion from BSSN fields to Cartesian ADM fields.",
        cfunc_type="void",
        name="BSSN_to_ADM",
        params=(
            f"const ot::Block& block, const {scalar_type}* const* in_gfs, "
            f"{scalar_type}* const* adm_gfs, const Point& domain_min, "
            "const Point& domain_max"
        ),
        body=body,
    )
    return pcg.NRPyEnv()
