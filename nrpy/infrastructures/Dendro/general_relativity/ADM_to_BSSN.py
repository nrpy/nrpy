"""
Generate the pointwise ADM-to-BSSN conversion for Dendro.

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
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.infrastructures.Dendro import state_h
from nrpy.infrastructures.Dendro.simple_loop import simple_loop


def register_CFunction_ADM_to_BSSN(
    solver_stem: str,
    *,
    fd_order: int,
    enable_fCCZ4: bool = False,
    CoordSystem: str = "Cartesian",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the pointwise conversion from puncture ADM data to evolved data.

    The input component order is ``gammaDD``, ``KDD``, ``betaU``, then ``BU``.
    The lapse is initialized to the converted conformal factor, so W evolution
    starts with ``alpha=W``. The derivative-dependent connection is written
    later by :mod:`initial_data_lambdaU`.

    :param solver_stem: Lowercase formulation name used by generated headers.
    :param fd_order: Centered finite-difference order naming this kernel variant.
    :param enable_fCCZ4: Include the fCCZ4 constraint scalar when true.
    :param CoordSystem: Reference-metric coordinate system.
    :return: The NRPy registries, or ``None`` during parallel collection.
    :raises ValueError: If the selected infrastructure, conformal factor,
        point-loop mode, or registered state is incompatible with this kernel.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError("ADM_to_BSSN requires Infrastructure='Dendro'.")
    if par.parval_from_str("EvolvedConformalFactor_cf") != "W":
        raise ValueError("The generated Dendro applications evolve W.")
    if par.parval_from_str("parallelization") != "none":
        raise ValueError("Dendro point kernels require parallelization='none'.")
    if fd_order not in (4, 6, 8):
        raise ValueError("Dendro ADM conversion supports FD orders 4, 6, and 8.")

    state_h.validate_registered_state(enable_fCCZ4)
    gammaDD = ixp.declarerank2("gammaDD", symmetry="sym01")
    KDD = ixp.declarerank2("KDD", symmetry="sym01")
    betaU = ixp.declarerank1("betaU")
    BU = ixp.declarerank1("BU")
    converted = ADM_to_BSSN(
        gammaDD,
        KDD,
        betaU,
        BU,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=False,
    )
    quantities = BSSN_quantities[CoordSystem]
    expressions: Dict[str, sp.Expr] = {
        "alpha": converted.cf,
        "cf": converted.cf,
        "trK": converted.trK,
    }
    for i in range(3):
        expressions[f"vetU{i}"] = converted.vetU[i]
        expressions[f"betU{i}"] = converted.betU[i]
    for tensor_name, tensor in (("hDD", converted.hDD), ("aDD", converted.aDD)):
        for i in range(3):
            for j in range(i, 3):
                expressions[f"{tensor_name}{i}{j}"] = tensor[i][j]
    if enable_fCCZ4:
        expressions["Theta_fCCZ4"] = sp.sympify(0)

    connection_names = {str(quantities.lambdaU[i]) for i in range(3)}
    evolved_names = tuple(state_h.evolved_gridfunctions(enable_fCCZ4))
    output_names = tuple(name for name in evolved_names if name not in connection_names)
    if set(output_names) != set(expressions):
        raise ValueError(
            "ADM_to_BSSN outputs differ from the canonical non-connection state: "
            f"outputs={output_names}, expressions={tuple(expressions)}."
        )

    scalar_type = gri.DENDRO_SCALAR_TYPE
    output_lvalues: List[str] = []
    for name in output_names:
        dendro_name = cast(
            gri.DendroGridFunction, gri.glb_gridfcs_dict[name]
        ).dendro_name
        output_lvalues.append(f"out_{dendro_name}[pp]")
    kernel = c_codegen(
        [expressions[name] for name in output_names],
        output_lvalues,
        include_braces=False,
        enable_simd=False,
        fp_type=str(par.parval_from_str("fp_type")),
        fp_type_alias=scalar_type,
        cse_sorting="none",
        verbose=False,
    )

    input_names = (
        "gammaDD00",
        "gammaDD01",
        "gammaDD02",
        "gammaDD11",
        "gammaDD12",
        "gammaDD22",
        "KDD00",
        "KDD01",
        "KDD02",
        "KDD11",
        "KDD12",
        "KDD22",
        "betaU0",
        "betaU1",
        "betaU2",
        "BU0",
        "BU1",
        "BU2",
    )
    input_bindings = [
        f"const {scalar_type}* adm_{name} = adm_gfs[{index}] + offset;"
        for index, name in enumerate(input_names)
    ]
    output_bindings = []
    for name in output_names:
        index = evolved_names.index(name)
        dendro_name = cast(
            gri.DendroGridFunction, gri.glb_gridfcs_dict[name]
        ).dendro_name
        output_bindings.append(
            f"{scalar_type}* out_{dendro_name} = out_gfs[{index}] + offset;"
        )
    local_inputs = [
        f"const {scalar_type} {name} = adm_{name}[pp];" for name in input_names
    ]
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
        subdirectory="generated/src/ADM_to_BSSN",
        includes=[f"{solver_stem}_defines.h"],
        desc="Pointwise conversion from TwoPunctures ADM fields to W-BSSN fields.",
        cfunc_type="void",
        name=f"ADM_to_BSSN_order_{fd_order}",
        params=(
            f"const ot::Block& block, const {scalar_type}* const* adm_gfs, "
            f"{scalar_type}* const* out_gfs, const Point& domain_min, "
            "const Point& domain_max"
        ),
        body=body,
    )
    return pcg.NRPyEnv()
