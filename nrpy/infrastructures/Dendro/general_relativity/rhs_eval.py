"""
Generate direct finite-difference BSSN or fCCZ4 right-hand sides for Dendro.

The formulation equations are assembled here and lowered directly into one
``ot::Block`` kernel per finite-difference order. The solver calls the
separate Ricci kernel immediately before this kernel without an intervening
communication step.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from collections import OrderedDict
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, Mapping, Tuple, Union, cast

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.equations.general_relativity.BSSN_gauge_RHSs import BSSN_gauge_RHSs
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs
from nrpy.equations.general_relativity.fCCZ4_RHSs import fCCZ4_RHSs
from nrpy.equations.general_relativity.kreiss_oliger_terms import (
    add_KreissOliger_dissipation_terms,
)
from nrpy.finite_difference import stencil_reach_per_axis
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list
from nrpy.infrastructures.Dendro import CodeParameters, state_h
from nrpy.infrastructures.Dendro.simple_loop import simple_loop

# (KO base order, required padding). NRPy's dKOD construction adds two to the
# base order, giving effective KO differences 4, 6, and 8.
DENDRO_FD_PROFILES: Mapping[int, Tuple[int, int]] = {
    4: (2, 2),
    6: (4, 3),
    8: (6, 4),
}


def register_CFunction_rhs_eval(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
    fd_order: int = 6,
    enable_KreissOliger_dissipation: bool = True,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
    enable_SSL: bool = False,
    enable_CAHD: bool = False,
    capture_validation_expressions: bool = False,
    enable_intrinsics: bool = True,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register one order-specific BSSN or fCCZ4 right-hand-side kernel.

    :param solver_stem: Lowercase formulation stem used by the generated header.
    :param enable_fCCZ4: Generate fCCZ4 instead of BSSN.
    :param fd_order: Centered finite-difference order (4, 6, or 8).
    :param enable_KreissOliger_dissipation: Add matching-order KO dissipation.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse gauge condition.
    :param ShiftEvolutionOption: Shift gauge condition.
    :param enable_SSL: Add slow-start lapse.
    :param enable_CAHD: Add formulation-specific Hamiltonian damping.
    :param capture_validation_expressions: Retain expanded expressions for validation.
    :param enable_intrinsics: Generate SIMD-intrinsic kernels; every vector stays
        inside its own row.
    :return: Updated NRPy registries, or ``None`` while collecting parallel work.
    :raises ValueError: If the requested configuration or state layout is invalid.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro RHS, got "
            f"{par.parval_from_str('Infrastructure')!r}."
        )
    if fd_order not in DENDRO_FD_PROFILES:
        raise ValueError(f"Unsupported fd_order={fd_order!r}; allowed: (4, 6, 8).")
    conformal_factor = par.parval_from_str("EvolvedConformalFactor_cf")
    if conformal_factor not in ("W", "chi"):
        raise ValueError("Dendro BSSN and fCCZ4 require W or chi.")
    if par.parval_from_str("parallelization") != "none":
        raise ValueError("Dendro point kernels require parallelization='none'.")
    if enable_intrinsics and CoordSystem != "Cartesian":
        raise ValueError("Dendro SIMD RHS kernels require Cartesian coordinates.")

    old_fd_order = par.parval_from_str("fd_order")
    par.set_parval_from_str("fd_order", fd_order)
    try:
        ko_fd_order, expected_padding = DENDRO_FD_PROFILES[fd_order]
        staged_coord_system = CoordSystem + "_RbarDD_gridfunctions"
        quantities = BSSN_quantities[staged_coord_system]

        if enable_fCCZ4:
            fccz4_rhs = fCCZ4_RHSs.get_rhs(
                staged_coord_system,
                enable_YBS_Gamma_constraint_adjustment=False,
                enable_YBS_momentum_constraint_adjustment=False,
            )
            rhs_by_symbol_name: Dict[str, sp.Expr] = OrderedDict(
                sorted(fccz4_rhs.fCCZ4_RHSs_varname_to_expr_dict.items())
            )
            alpha_rhs, vet_rhsU, bet_rhsU = BSSN_gauge_RHSs(
                CoordSystem=CoordSystem,
                enable_rfm_precompute=False,
                enable_T4munu=False,
                LapseEvolutionOption=LapseEvolutionOption,
                ShiftEvolutionOption=ShiftEvolutionOption,
                enable_YBS_Gamma_constraint_adjustment=False,
                evolved_connection_rhsU=fccz4_rhs.Lambdatilde_rhsU,
            )
            if LapseEvolutionOption == "OnePlusLog":
                alpha_rhs += 4 * quantities.alpha * fccz4_rhs.Theta
        else:
            bssn_rhs = BSSN_RHSs.get_rhs(
                staged_coord_system,
                enable_YBS_Gamma_constraint_adjustment=False,
                enable_YBS_momentum_constraint_adjustment=False,
            )
            rhs_by_symbol_name = OrderedDict(
                sorted(bssn_rhs.BSSN_RHSs_varname_to_expr_dict.items())
            )
            alpha_rhs, vet_rhsU, bet_rhsU = BSSN_gauge_RHSs(
                CoordSystem=CoordSystem,
                enable_rfm_precompute=False,
                enable_T4munu=False,
                LapseEvolutionOption=LapseEvolutionOption,
                ShiftEvolutionOption=ShiftEvolutionOption,
            )

        rhs_by_symbol_name["alpha_rhs"] = alpha_rhs
        for i in range(3):
            rhs_by_symbol_name[f"vet_rhsU{i}"] = vet_rhsU[i]
            rhs_by_symbol_name[f"bet_rhsU{i}"] = bet_rhsU[i]
        rhs_by_symbol_name = OrderedDict(sorted(rhs_by_symbol_name.items()))

        if enable_KreissOliger_dissipation:
            add_KreissOliger_dissipation_terms(
                rhs_by_symbol_name,
                CoordSystem=CoordSystem,
                enable_rfm_precompute=False,
                registering_module=__name__,
                ShiftEvolutionOption=ShiftEvolutionOption,
                KreissOliger_strength_gauge=0.3,
                KreissOliger_strength_nongauge=0.3,
                enable_CAKO=False,
                W=(
                    sp.sqrt(quantities.cf)
                    if conformal_factor == "chi"
                    else quantities.cf
                ),
                include_Theta_fCCZ4=enable_fCCZ4,
            )

        ssl_exponent = sp.Integer(0)
        if enable_SSL:
            SSL_h, SSL_sigma = par.register_CodeParameters(
                "REAL",
                __name__,
                ["SSL_h", "SSL_sigma"],
                [0.6, 20.0],
                add_to_parfile=True,
            )
            W = sp.sqrt(quantities.cf) if conformal_factor == "chi" else quantities.cf
            ssl_exponent = -(sp.Symbol("stage_time") ** 2) / (2 * SSL_sigma**2)
            rhs_by_symbol_name["alpha_rhs"] -= (
                W * SSL_h * sp.Symbol("SSL_exp_factor") * (quantities.alpha - W)
            )

        if enable_CAHD:
            gate = "register_M_and_LAMBDA_CONSTRAINT_gridfunctions"
            previous_gate = par.parval_from_str(gate)
            par.set_parval_from_str(gate, False)
            try:
                hamiltonian_constraint = BSSN_constraints[staged_coord_system].H
            finally:
                par.set_parval_from_str(gate, previous_gate)
            if enable_fCCZ4:
                c_cahd = par.register_CodeParameter(
                    "REAL", __name__, "C_CAHD", 0.15, add_to_parfile=True
                )
                rhs_by_symbol_name["cf_rhs"] += (
                    (4 if conformal_factor == "chi" else 2)
                    * c_cahd
                    * quantities.cf
                    * hamiltonian_constraint
                    * sp.Symbol("time_step")
                )
            else:
                c_cahd = par.register_CodeParameter(
                    "REAL", __name__, "C_CAHD", 0.06, add_to_parfile=True
                )
                rhs_by_symbol_name["cf_rhs"] += (
                    (1 if conformal_factor == "chi" else sp.Rational(1, 2))
                    * c_cahd
                    * quantities.cf
                    * hamiltonian_constraint
                    * sp.Symbol("grid_spacing") ** 2
                    / (
                        sp.Symbol("time_step")
                        * (1 + 10 * sp.Symbol("grid_spacing") ** 2)
                    )
                )

        rhs_by_gridfunction_name: Dict[str, sp.Expr] = OrderedDict()
        rhs_by_gridfunction_name["alpha"] = rhs_by_symbol_name["alpha_rhs"]
        rhs_by_gridfunction_name["cf"] = rhs_by_symbol_name["cf_rhs"]
        rhs_by_gridfunction_name["trK"] = rhs_by_symbol_name["trK_rhs"]
        for i in range(3):
            rhs_by_gridfunction_name[f"lambdaU{i}"] = rhs_by_symbol_name[
                f"lambda_rhsU{i}"
            ]
        for i in range(3):
            rhs_by_gridfunction_name[f"vetU{i}"] = rhs_by_symbol_name[f"vet_rhsU{i}"]
        for i in range(3):
            rhs_by_gridfunction_name[f"betU{i}"] = rhs_by_symbol_name[f"bet_rhsU{i}"]
        for tensor_name, rhs_name in (("hDD", "h_rhsDD"), ("aDD", "a_rhsDD")):
            for i in range(3):
                for j in range(i, 3):
                    rhs_by_gridfunction_name[f"{tensor_name}{i}{j}"] = (
                        rhs_by_symbol_name[f"{rhs_name}{i}{j}"]
                    )
        if enable_fCCZ4:
            rhs_by_gridfunction_name["Theta_fCCZ4"] = rhs_by_symbol_name[
                "Theta_fCCZ4_rhs"
            ]

        # Dendro uses centered derivatives. Build the finite substitution set
        # from known state families; never infer names by walking expressions.
        centered_advection_substitutions: Dict[sp.Basic, sp.Basic] = {}
        for scalar_name in ("alpha", "cf", "trK", "Theta_fCCZ4"):
            for direction in range(3):
                for directional_operator in ("dupD", "ddnD"):
                    centered_advection_substitutions[
                        sp.Symbol(f"{scalar_name}_{directional_operator}{direction}")
                    ] = sp.Symbol(f"{scalar_name}_dD{direction}")
        for vector_name in ("lambdaU", "vetU", "betU"):
            for component in range(3):
                for direction in range(3):
                    for directional_operator in ("dupD", "ddnD"):
                        centered_advection_substitutions[
                            sp.Symbol(
                                f"{vector_name}_{directional_operator}{component}{direction}"
                            )
                        ] = sp.Symbol(f"{vector_name}_dD{component}{direction}")
        for tensor_name in ("hDD", "aDD"):
            for i in range(3):
                for j in range(i, 3):
                    for direction in range(3):
                        for directional_operator in ("dupD", "ddnD"):
                            centered_advection_substitutions[
                                sp.Symbol(
                                    f"{tensor_name}_{directional_operator}{i}{j}{direction}"
                                )
                            ] = sp.Symbol(f"{tensor_name}_dD{i}{j}{direction}")
        rhs_by_gridfunction_name = OrderedDict(
            (name, expression.xreplace(centered_advection_substitutions))
            for name, expression in rhs_by_gridfunction_name.items()
        )

        evol_order = tuple(state_h.evolved_gridfunctions(enable_fCCZ4))
        state_h.validate_registered_state(enable_fCCZ4)
        if tuple(rhs_by_gridfunction_name) != evol_order:
            raise ValueError(
                "RHS keys must equal the canonical EVOL order: "
                f"got={tuple(rhs_by_gridfunction_name)!r} expected={evol_order!r}."
            )

        lvalues = []
        for name in evol_order:
            dendro_name = cast(
                gri.DendroGridFunction, gri.glb_gridfcs_dict[name]
            ).dendro_name
            lvalues.append(f"rhs_{dendro_name}[pp]")
        if len(lvalues) != len(set(lvalues)):
            raise ValueError("RHS outputs must map bijectively onto EVOL storage.")

        kernel_expressions = list(rhs_by_gridfunction_name.values())
        fp_type = str(par.parval_from_str("fp_type"))
        scalar_type = gri.DENDRO_SCALAR_TYPE
        prefix = "NOSIMD" if enable_intrinsics else ""
        kernel = c_codegen(
            kernel_expressions,
            lvalues,
            enable_simd=enable_intrinsics,
            enable_fd_codegen=True,
            enable_fd_functions=False,
            fp_type=fp_type,
            fp_type_alias=scalar_type,
            mem_alloc_style="210",
            rational_const_alias="static const",
            cse_sorting="none",
            verbose=False,
            upwind_control_vec=sp.Symbol("unset"),
            ko_fd_order=ko_fd_order,
        )
        padding = max(
            stencil_reach_per_axis(
                kernel_expressions,
                "unset",
                fd_order,
                ko_fd_order=ko_fd_order,
            )
        )
        if padding != expected_padding:
            raise ValueError(
                f"Dendro fd_order={fd_order} requires padding {expected_padding}, "
                f"but emitted operators reach {padding}."
            )

        parameter_symbols, commondata_symbols = (
            get_params_commondata_symbols_from_expr_list(
                kernel_expressions + [ssl_exponent]
            )
        )
        used_codeparameters = tuple(
            sorted(
                name
                for name in parameter_symbols + commondata_symbols
                if name in par.glb_code_params_dict
            )
        )
        geometry = f"""const std::ptrdiff_t offset = static_cast<std::ptrdiff_t>(block.getOffset());
const unsigned nx_block = block.getAllocationSzX();
const unsigned ny_block = block.getAllocationSzY();
const unsigned nz_block = block.getAllocationSzZ();
const unsigned padding_block = block.get1DPadWidth();
if (padding_block < {expected_padding}) {{
    throw std::invalid_argument("rhs_eval block padding is too small for FD{fd_order}");
}}
const {scalar_type} dx_block[3] = {{
    block.computeDx(domain_min, domain_max),
    block.computeDy(domain_min, domain_max),
    block.computeDz(domain_min, domain_max)}};
const {scalar_type} {prefix}grid_spacing = dx_block[0];
[[maybe_unused]] const {scalar_type} pmin_block[3] = {{
    GRIDX_TO_X(block.getBlockNode().minX()) - padding_block * dx_block[0],
    GRIDY_TO_Y(block.getBlockNode().minY()) - padding_block * dx_block[1],
    GRIDZ_TO_Z(block.getBlockNode().minZ()) - padding_block * dx_block[2]}};"""
        input_bindings = []
        rhs_bindings = []
        for index, name in enumerate(evol_order):
            dendro_name = cast(
                gri.DendroGridFunction, gri.glb_gridfcs_dict[name]
            ).dendro_name
            input_bindings.append(
                f"const {scalar_type}* in_{dendro_name} = in_gfs[{index}] + offset;"
            )
            rhs_bindings.append(
                f"{scalar_type}* rhs_{dendro_name} = rhs_gfs[{index}] + offset;"
            )
        ricci_bindings = [
            f"const {scalar_type}* in_{name} = ricci_gfs[{index}] + offset;"
            for index, name in enumerate(state_h.RICCI_GRIDFUNCTIONS)
        ]
        preloop = []
        if enable_SSL:
            preloop.append(
                f"const {scalar_type} {prefix}SSL_exp_factor = "
                f"std::exp(-({prefix}stage_time * {prefix}stage_time) / "
                f"(2 * {prefix}SSL_sigma * {prefix}SSL_sigma));"
            )
        if enable_intrinsics:
            broadcast_names = (
                *used_codeparameters,
                "stage_time",
                "time_step",
                "grid_spacing",
                *(("SSL_exp_factor",) if enable_SSL else ()),
            )
            preloop += [
                f"[[maybe_unused]] const REAL_SIMD_ARRAY {name} = ConstSIMD(NOSIMD{name});"
                for name in broadcast_names
            ]
        block_body = "\n".join(
            (
                geometry,
                *input_bindings,
                *ricci_bindings,
                *rhs_bindings,
                *preloop,
                simple_loop(
                    kernel,
                    nx="nx_block",
                    ny="ny_block",
                    nz="nz_block",
                    padding="padding_block",
                    pmin_padded="pmin_block",
                    dx="dx_block",
                    enable_intrinsics=enable_intrinsics,
                ),
            )
        )
        cparam_args = ", ".join(
            "const "
            f"{CodeParameters.c_type(par.glb_code_params_dict[name].cparam_type)} "
            f"{prefix}{name}"
            for name in used_codeparameters
        )
        block_params = (
            f"const ot::Block& block, const {scalar_type}* const* in_gfs, "
            f"const {scalar_type}* const* ricci_gfs, {scalar_type}* const* rhs_gfs, "
            f"const Point& domain_min, const Point& domain_max, "
            f"const {scalar_type} {prefix}stage_time, "
            f"const {scalar_type} {prefix}time_step"
            + (f", {cparam_args}" if cparam_args else "")
        )
        formulation = "fCCZ4" if enable_fCCZ4 else "BSSN"
        cfc.register_CFunction(
            subdirectory="generated/src/rhs_eval",
            # simd_intrinsics.h must precede the definitions header, which can
            # reach a copy with the same include guard; "./" keeps this order
            # after clang-format sorts the includes.
            includes=(["./simd_intrinsics.h"] if enable_intrinsics else [])
            + [f"{solver_stem}_defines.h"],
            desc=f"Per-block direct-FD {formulation} RHS ({len(evol_order)} fields).",
            cfunc_type="void",
            name=f"rhs_eval_order_{fd_order}",
            params=block_params,
            body=block_body,
            ET_current_thorn_CodeParams_used=list(used_codeparameters),
        )

        if capture_validation_expressions:
            validation = cast(
                Dict[int, Dict[str, Dict[str, object]]],
                par.glb_extras_dict.setdefault("Dendro", {}).setdefault(
                    "validation_candidates", {}
                ),
            )
            candidates = validation.setdefault(fd_order, {}).setdefault("rhs", {})
            owner = (
                f"rhs_eval_order_{fd_order}"
                f"|fccz4={int(enable_fCCZ4)}|coord={CoordSystem}"
                f"|lapse={LapseEvolutionOption}|shift={ShiftEvolutionOption}"
                f"|ko={int(enable_KreissOliger_dissipation)}"
                f"|ssl={int(enable_SSL)}|cahd={int(enable_CAHD)}"
            )
            if owner in candidates:
                raise ValueError(
                    "Validation RHS expressions already exist for exact worker "
                    f"configuration {owner}."
                )
            expanded_quantities = BSSN_quantities[CoordSystem]
            ricci_substitutions = {
                sp.Symbol(name): expression
                for name, expression in zip(
                    expanded_quantities.Ricci_varnames,
                    expanded_quantities.Ricci_exprs,
                )
            }
            candidates[owner] = {
                name: expression.xreplace(
                    {sp.Symbol("SSL_exp_factor"): sp.exp(ssl_exponent)}
                ).xreplace(ricci_substitutions)
                for name, expression in rhs_by_gridfunction_name.items()
            }
    finally:
        par.set_parval_from_str("fd_order", old_fd_order)
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
