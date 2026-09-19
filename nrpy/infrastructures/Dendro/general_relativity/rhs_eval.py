"""
Direct finite-difference BSSN or fCCZ4 right-hand sides for Dendro.

One module, one entry point, and the formulation chosen by the
``enable_fCCZ4`` argument, as ``BHaH/general_relativity/rhs_eval.py`` does for
the same two formulations.  The two expression-assembly routes stay separate
below: fCCZ4 comes from the shared ``fCCZ4_system`` factory, BSSN from
``BSSN_RHSs`` plus the gauge terms and optional Kreiss-Oliger dissipation, and
BHaH's own divergence guard records that unifying those branches is not wanted.

Everything downstream of the expressions is shared: the point kernel is lowered
by ``c_codegen`` with direct finite differences, wrapped in the Dendro interior
point loop, and registered as a per-block, an all-block and a flat-block
CFunction.  This module authors no field name, no physics default and no
finite-difference coefficient.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from collections import OrderedDict
from dataclasses import dataclass
from typing import Dict, Iterable, Mapping, Tuple

import sympy as sp

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.BSSN_gauge_RHSs import BSSN_gauge_RHSs
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs
from nrpy.equations.general_relativity.fCCZ4_system import (
    build_fccz4_expression_bundle,
)
from nrpy.equations.general_relativity.kreiss_oliger_terms import (
    add_KreissOliger_dissipation_terms,
)
from nrpy.finite_difference import (
    extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars,
)
from nrpy.infrastructures.Dendro import CFunction_roles as roles
from nrpy.infrastructures.Dendro import block_kernel_helpers as bkh
from nrpy.infrastructures.Dendro import gridfunction_name_decorations as gf_names
from nrpy.infrastructures.Dendro.general_relativity import generation_parameters
from nrpy.infrastructures.Dendro.simple_loop import (
    block_loop,
    require_serial_parallelization,
)

# The emitted CFunction names are the established NRPy operation name,
# prefixed with the solver stem the caller supplies, exactly as initial_data.py
# and enforce_detgbar_equals_detghat_trAzero.py spell theirs.  The formulation
# is carried by the stem, so one suffix serves both formulations.  The LTS
# flat-block adapter uses the same numerical body with a flat layout.
RHS_EVAL_BLOCK_SUFFIX = "rhs_eval_block"
RHS_EVAL_ALL_BLOCKS_SUFFIX = "rhs_eval"
RHS_EVAL_FLAT_BLOCK_SUFFIX = "rhs_eval_flat_block"

# (KO base order, required padding).  NRPy's dKOD construction adds two to
# the base order, so these profiles produce effective KO differences 4/6/8.
DENDRO_FD_PROFILES: Mapping[int, Tuple[int, int]] = {
    4: (2, 2),
    6: (4, 3),
    8: (6, 4),
}


@dataclass(frozen=True)
class RHSBuild:
    """
    Immutable result of building the direct-FD RHS for one profile.

    :param evol_order: The EVOL gridfunction names, in registry order.
    :param fd_order: Centered finite-difference order.
    :param ko_fd_order: Base order supplied to NRPy's ``dKOD`` construction.
    :param ko_enabled: Whether the emitted RHS contains KO dissipation.
    :param lvalues: The output lvalues (``rhs_<name>[pp]``).
    :param padding: Ghost points the emitted operators reach, widest axis.
    :param block_body: The per-block CFunction body (point loop + bindings).
    :param block_params: The per-block CFunction parameter list.
    :param all_blocks_body: The all-block CFunction body (NRPy block loop).
    :param all_blocks_params: The all-block CFunction parameter list.
    :param flat_block_body: The flat-block adapter CFunction body.
    :param flat_block_params: The flat-block adapter CFunction parameter list.
    :param used_codeparameters: The CodeParameters the signatures forward, in
        signature order, computed from the expression free symbols.
    :param rhs_by_symbol_name: The assembled symbolic right-hand sides, kept so
        the module's ``__main__`` can pin them against trusted values without
        reassembling them.
    """

    evol_order: Tuple[str, ...]
    fd_order: int
    ko_fd_order: int
    ko_enabled: bool
    lvalues: Tuple[str, ...]
    padding: int
    block_body: str
    block_params: str
    all_blocks_body: str
    all_blocks_params: str
    flat_block_body: str
    flat_block_params: str
    used_codeparameters: Tuple[str, ...]
    rhs_by_symbol_name: Mapping[str, sp.Expr]


# The BSSN evolved state: hDD (6), aDD (6), cf, trK, lambdaU (3), alpha,
# vetU (3), betU (3).
BSSN_EVOL_COUNT = 24


def _directional_operators(expressions: Iterable[sp.Expr]) -> Tuple[str, ...]:
    """
    Return directional finite-difference operators in an expression sequence.

    :param expressions: Symbolic expressions to inspect.
    :return: Sorted unique directional finite-difference operator names.
    """
    directional_families = ("dupD", "ddnD", "dfullupD", "dfulldnD")
    candidates = []
    for expression in expressions:
        for symbol in expression.free_symbols:
            suffix = str(symbol).rsplit("_", maxsplit=1)[-1]
            if any(
                suffix.startswith(family) and suffix[len(family) :].isdigit()
                for family in directional_families
            ):
                candidates.append(symbol)
    if not candidates:
        return ()
    _, operators = extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars(
        candidates
    )
    return tuple(sorted(set(operators)))


def BSSN_rhs_expressions(
    *,
    CoordSystem: str,
    LapseEvolutionOption: str,
    ShiftEvolutionOption: str,
    enable_KreissOliger_dissipation: bool,
) -> Dict[str, sp.Expr]:
    """
    Assemble the BSSN RHS expression set.

    The assembly order matches ETLegacy's and BHaH's: the non-gauge RHSs come
    from the cached ``BSSN_RHSs`` object, the gauge RHSs are added to a copy of
    its dictionary, Kreiss-Oliger terms are applied through the shared helper,
    and the Dendro lowering later centers directional advection derivatives.

    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :return: RHS expressions keyed by their symbolic output names.
    """
    rhs = BSSN_RHSs[CoordSystem]
    alpha_rhs, vet_rhsU, bet_rhsU = BSSN_gauge_RHSs(
        CoordSystem=CoordSystem,
        enable_rfm_precompute=False,
        enable_T4munu=False,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    # Copy before adding: BSSN_RHSs caches its dictionary, and mutating it
    # would leak this profile's gauge and dissipation choices into every later
    # generated RHS in the same process.
    rhs_by_symbol_name: Dict[str, sp.Expr] = OrderedDict(
        sorted(rhs.BSSN_RHSs_varname_to_expr_dict.items())
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
            # W is consumed only under enable_CAKO, which this profile does
            # not offer, so no conformal factor needs deriving here.
            W=sp.sympify(1),
            # BSSN has no Z4 scalar; the fCCZ4 profile is the only one that
            # dissipates Theta_fCCZ4.
            include_Theta_fCCZ4=False,
        )

    return rhs_by_symbol_name


def build_rhs_eval(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
    fd_order: int = 6,
    enable_KreissOliger_dissipation: bool = True,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
) -> RHSBuild:
    """
    Build the Dendro right-hand-side kernels for one formulation.

    One entry point for both formulations, taking the formulation as an
    argument, exactly as ``BHaH/general_relativity/rhs_eval.py`` does.  Expression assembly remains formulation-specific; lowering and registration
    share one implementation.

    :param solver_stem: Lowercase formulation stem, prefixed onto every
        emitted CFunction name.
    :param enable_fCCZ4: Build fCCZ4 instead of BSSN.
    :param fd_order: The centered finite-difference order (4, 6, or 8).
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :return: The formulation's build record.
    :raises ValueError: If the registered fields do not match the formulation.

    Doctests:
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fp_type", "double")
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    >>> par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    >>> import nrpy.grid as _gri
    >>> _gri.glb_gridfcs_dict.clear()
    >>> par.glb_extras_dict.pop("Dendro", None) and None
    >>> try:
    ...     build_rhs_eval("fccz4", enable_fCCZ4=True, fd_order=2, enable_KreissOliger_dissipation=False)
    ... except ValueError as error:
    ...     print(str(error).splitlines()[0])
    Unsupported fd_order=2; allowed: (4, 6, 8).
    >>> import contextlib, io
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _build = build_rhs_eval(
    ...         "fccz4", enable_fCCZ4=True, fd_order=4, enable_KreissOliger_dissipation=False
    ...     )
    >>> len(_build.evol_order), "Theta_fCCZ4" in _build.evol_order
    (25, True)
    >>> sorted(_build.lvalues)[:2]
    ['rhs_Theta_fCCZ4[pp]', 'rhs_aDD00[pp]']
    >>> _build.padding
    2
    >>> _build.ko_fd_order
    2
    >>> roles.required_padding()
    2


    >>> import contextlib, io
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fp_type", "double")
    >>> par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    >>> import nrpy.grid as _gri
    >>> from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities as _bq
    >>> from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs as _brhs
    >>> _gri.glb_gridfcs_dict.clear()
    >>> par.glb_extras_dict.pop("Dendro", None) and None
    >>> # The equations layer memoizes per CoordSystem and registers the evolved
    >>> # state in its constructor, so the memo must go with the registry or the
    >>> # rebuild finds nothing registered.
    >>> _bq.clear() or _brhs.clear()
    >>> try:
    ...     build_rhs_eval("bssn", fd_order=2, enable_KreissOliger_dissipation=False)
    ... except ValueError as error:
    ...     print(str(error).splitlines()[0])
    Unsupported fd_order=2; allowed: (4, 6, 8).
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _build = build_rhs_eval(
    ...         "bssn", fd_order=4, enable_KreissOliger_dissipation=False
    ...     )
    >>> len(_build.evol_order), "Theta_fCCZ4" in _build.evol_order
    (24, False)
    >>> sorted(_build.lvalues)[:2]
    ['rhs_aDD00[pp]', 'rhs_aDD01[pp]']
    >>> _build.padding
    2
    >>> _build.ko_fd_order
    2
    >>> roles.required_padding()
    2

    The mapping failure check is exercised through this public builder before
    lowering.  The original ``BSSN_rhs_expressions`` binding is restored even if
    an assertion fails.

    >>> _owner_globals = build_rhs_eval.__globals__
    >>> _original_bssn_expressions = _owner_globals["BSSN_rhs_expressions"]
    >>> _good_rhs = OrderedDict(_build.rhs_by_symbol_name)
    >>> _fake_rhs = _good_rhs
    >>> def _temporary_bssn_expressions(**_kwargs):
    ...     return _fake_rhs
    >>> try:
    ...     _owner_globals["BSSN_rhs_expressions"] = _temporary_bssn_expressions
    ...     _missing = OrderedDict(_good_rhs)
    ...     _ = _missing.pop(next(iter(_missing)))
    ...     _fake_rhs = _missing
    ...     try:
    ...         build_rhs_eval("bssn", fd_order=4, enable_KreissOliger_dissipation=False)
    ...     except ValueError as error:
    ...         assert "missing=" in str(error)
    ...     else:
    ...         raise AssertionError("missing RHS target was accepted")
    ...     _extra = OrderedDict(_good_rhs)
    ...     _extra["not_registered_rhs"] = sp.Integer(0)
    ...     _fake_rhs = _extra
    ...     try:
    ...         build_rhs_eval("bssn", fd_order=4, enable_KreissOliger_dissipation=False)
    ...     except ValueError as error:
    ...         assert "extra=['not_registered']" in str(error)
    ...     else:
    ...         raise AssertionError("extra RHS target was accepted")
    ...     _repeated = OrderedDict(_good_rhs)
    ...     _repeated["aDD00_rhs"] = _repeated["a_rhsDD00"]
    ...     _fake_rhs = _repeated
    ...     try:
    ...         build_rhs_eval("bssn", fd_order=4, enable_KreissOliger_dissipation=False)
    ...     except ValueError as error:
    ...         assert "do not map bijectively" in str(error)
    ...     else:
    ...         raise AssertionError("repeated RHS target was accepted")
    ... finally:
    ...     _owner_globals["BSSN_rhs_expressions"] = _original_bssn_expressions

    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro RHS, got "
            f"{par.parval_from_str('Infrastructure')!r}."
        )
    if fd_order not in DENDRO_FD_PROFILES:
        raise ValueError(f"Unsupported fd_order={fd_order!r}; allowed: (4, 6, 8).")
    ko_fd_order, expected_padding = DENDRO_FD_PROFILES[fd_order]
    # Current GR initial-data kernels require a conformal-factor representation
    # whose Minkowski value matches the registered asymptotic field value.
    generation_parameters.validate_generation_parameters()
    # Dendro owns the outer block traversal. The qualified point kernel is
    # serial so it does not introduce nested inner-loop parallelism. Reject a
    # different request instead of silently generating an unqualified kernel.
    require_serial_parallelization()
    formulation = "fCCZ4" if enable_fCCZ4 else "BSSN"
    expected_evol_count = 25 if enable_fCCZ4 else BSSN_EVOL_COUNT
    if enable_fCCZ4:
        bundle = build_fccz4_expression_bundle(
            CoordSystem=CoordSystem,
            LapseEvolutionOption=LapseEvolutionOption,
            ShiftEvolutionOption=ShiftEvolutionOption,
            enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        )
        rhs_by_symbol_name = bundle.rhs_by_symbol_name
    else:
        rhs_by_symbol_name = BSSN_rhs_expressions(
            CoordSystem=CoordSystem,
            LapseEvolutionOption=LapseEvolutionOption,
            ShiftEvolutionOption=ShiftEvolutionOption,
            enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        )
    centered: Dict[str, sp.Expr] = OrderedDict()
    for rhs_name, expression in rhs_by_symbol_name.items():
        replacements: Dict[sp.Basic, sp.Basic] = {}
        for symbol in expression.free_symbols:
            name = str(symbol)
            separator = name.rfind("_")
            if separator < 0:
                continue
            suffix_position = separator + 1
            suffix = name[suffix_position:]
            prefix = next(
                (
                    candidate
                    for candidate in ("dupD", "ddnD")
                    if suffix.startswith(candidate)
                    and suffix[len(candidate) :].isdigit()
                ),
                None,
            )
            if prefix is None:
                continue
            _, derivative_operators = (
                extract_base_gfs_and_deriv_ops_lists__from_list_of_deriv_vars([symbol])
            )
            if len(derivative_operators) != 1 or not derivative_operators[0].startswith(
                ("dupD", "ddnD")
            ):
                continue
            # Preserve every tensor-component digit following the derivative
            # suffix. For example, aDD_dupD000 becomes aDD_dD000, not
            # aDD00_dD0. The parser above validates the derivative symbol;
            # this replacement changes only its final operator suffix.
            replacement_name = name[:suffix_position] + "dD" + suffix[len(prefix) :]
            replacements[symbol] = sp.Symbol(replacement_name, **symbol.assumptions0)
        centered[rhs_name] = expression.xreplace(replacements)
    remaining = _directional_operators(centered.values())
    if remaining:
        raise ValueError(
            "Dendro RHS contains directional derivative operators after centered "
            f"normalization: {remaining}."
        )
    rhs_by_symbol_name = centered
    kernel_expressions = list(rhs_by_symbol_name.values())
    remaining_directional = _directional_operators(kernel_expressions)
    if remaining_directional:
        raise ValueError(
            "Dendro RHS contains directional derivative operators before code "
            f"generation: {remaining_directional}."
        )
    rhs_symbols = tuple(rhs_by_symbol_name)
    lvalues = tuple(
        f"rhs_{gf_names.rhs_symbol_to_gridfunction_name(symbol)}[pp]"
        for symbol in rhs_symbols
    )
    evol_order = roles.registered_evol_order()
    if len(evol_order) != expected_evol_count:
        raise ValueError(
            f"{formulation} EVOL registry must hold exactly {expected_evol_count} fields, found "
            f"{len(evol_order)}."
        )
    # The RHS symbols must map bijectively onto the registered EVOL fields
    # (under their exact registered names).
    mapped = {
        gf_names.rhs_symbol_to_gridfunction_name(symbol) for symbol in rhs_symbols
    }
    if len(lvalues) != len(set(lvalues)) or mapped != set(evol_order):
        raise ValueError(
            f"{formulation} RHS symbols do not map bijectively onto the registered "
            f"EVOL fields: missing={sorted(set(evol_order) - mapped)} "
            f"extra={sorted(mapped - set(evol_order))}"
        )
    par.set_parval_from_str("fd_order", fd_order)
    # Both scalar spellings come from the registries,
    # and every codegen option is passed explicitly.
    fp_type = str(par.parval_from_str("fp_type"))
    scalar_type = gri.DENDRO_SCALAR_TYPE
    kernel = c_codegen(
        list(rhs_by_symbol_name.values()),
        list(lvalues),
        enable_fd_codegen=True,
        enable_fd_functions=False,
        enable_simd=False,
        fp_type=fp_type,
        fp_type_alias=scalar_type,
        mem_alloc_style="210",
        rational_const_alias="static const",
        verbose=False,
        upwind_control_vec=sp.Symbol("unset"),
        ko_fd_order=ko_fd_order,
    )
    # Padding is the widest reach of the derivative operators the emitted
    # kernel actually contains, taken per axis from the same coefficient
    # source used to generate the kernel. The Dendro profiles make the KO and
    # centered regular derivatives fit the same 2-, 3-, or 4-point padding.
    # The consumed CodeParameters are the expression free symbols that are
    # registered CodeParameters.  Reading the symbols the equations actually
    # contain is exact; scanning the emitted C text for names is not.  Sorted
    # for a deterministic CFunction parameter order (the caller forwards the
    # values in this order).
    used_codeparameters = bkh.used_codeparameters(rhs_by_symbol_name.values())
    point_loop_body = bkh.point_loop(kernel)
    block_body = (
        bkh.block_pointer_bindings(evol_order, scalar_type) + "\n" + point_loop_body
    )
    all_blocks_body = block_loop(
        f"{solver_stem}_{RHS_EVAL_BLOCK_SUFFIX}(mesh.geom[blk], in_gfs, rhs_gfs"
        + (
            f", {bkh.cparam_arguments(used_codeparameters)}"
            if used_codeparameters
            else ""
        )
        + ");",
        num_blocks="mesh.num_blocks",
    )
    # Single-body rule: the flat adapter derives flat-layout
    # pointers, packs the per-component pointer arrays, and calls the
    # registered block kernel.  There is exactly one numerical body, so the
    # two paths cannot diverge.
    flat_call_args = (
        f", {bkh.cparam_arguments(used_codeparameters)}" if used_codeparameters else ""
    )
    flat_block_body = (
        bkh.flat_block_pointer_bindings(evol_order, scalar_type)
        + f"\nconst {scalar_type}* const in_gfs_call[] = {{"
        + ", ".join(gf_names.input_pointer(name) for name in evol_order)
        + f"}};\n{scalar_type}* rhs_gfs_call[] = {{"
        + ", ".join(gf_names.rhs_pointer(name) for name in evol_order)
        + "};\n"
        + f"{solver_stem}_{RHS_EVAL_BLOCK_SUFFIX}"
        + f"(geom, in_gfs_call, rhs_gfs_call{flat_call_args});\n"
    )
    cparam_args = bkh.cparam_declarations(used_codeparameters)
    block_params = (
        f"const block_geometry_struct& geom, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* rhs_gfs" + (f", {cparam_args}" if cparam_args else "")
    )
    # The all-block wrapper retains its standalone mesh interface. Real host
    # contexts normalize each ot::Block into the shared block_geometry_struct ABI and
    # call the registered block kernel directly.
    all_blocks_params = (
        f"const standalone_host_mesh_struct& mesh, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* rhs_gfs" + (f", {cparam_args}" if cparam_args else "")
    )
    flat_block_params = (
        f"const block_geometry_struct& geom, const {scalar_type}* const in_gfs_flat, "
        f"{scalar_type}* const rhs_gfs_flat"
        + (f", {cparam_args}" if cparam_args else "")
    )
    operators = bkh.emitted_derivative_operators(kernel_expressions)
    padding = bkh.padding_from_derivative_operators(
        kernel_expressions,
        fd_order,
        ko_fd_order=ko_fd_order,
    )
    remaining_directional = _directional_operators(kernel_expressions)
    if remaining_directional:
        raise ValueError(
            "Dendro RHS contains directional derivative operators after centered "
            f"lowering: {remaining_directional}."
        )
    if padding != expected_padding:
        raise ValueError(
            f"Dendro fd_order={fd_order} requires padding {expected_padding}, "
            f"but emitted operators reach {padding}."
        )
    # Add dKOD operators to the emitted kernel exactly once if
    # and only if Kreiss-Oliger dissipation was requested.
    if (
        any(operator.startswith("dKOD") for operator in operators)
        != enable_KreissOliger_dissipation
    ):
        raise ValueError(
            "dKOD-operator presence does not match "
            f"enable_KreissOliger_dissipation={enable_KreissOliger_dissipation!r}."
        )
    roles.set_required_padding(padding)
    return RHSBuild(
        evol_order=evol_order,
        fd_order=fd_order,
        ko_fd_order=ko_fd_order,
        ko_enabled=enable_KreissOliger_dissipation,
        lvalues=lvalues,
        padding=padding,
        block_body=block_body,
        block_params=block_params,
        all_blocks_body=all_blocks_body,
        all_blocks_params=all_blocks_params,
        flat_block_body=flat_block_body,
        flat_block_params=flat_block_params,
        used_codeparameters=tuple(used_codeparameters),
        rhs_by_symbol_name=dict(rhs_by_symbol_name),
    )


def register_CFunctions_rhs_eval(
    solver_stem: str,
    *,
    enable_fCCZ4: bool = False,
    fd_order: int = 6,
    enable_KreissOliger_dissipation: bool = True,
    CoordSystem: str = "Cartesian",
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant__Hatted",
) -> RHSBuild:
    """
    Register the right-hand-side CFunctions for one formulation.

    :param solver_stem: Lowercase formulation stem, prefixed onto every
        emitted CFunction name and onto the emitted include.
    :param enable_fCCZ4: Register fCCZ4 instead of BSSN.
    :param fd_order: The centered finite-difference order (4, 6, or 8).
    :param enable_KreissOliger_dissipation: Enable Kreiss-Oliger dissipation.
    :param CoordSystem: Reference-metric coordinate system.
    :param LapseEvolutionOption: Lapse evolution option.
    :param ShiftEvolutionOption: Shift evolution option.
    :return: Registered RHS build record with canonical expressions and order.
    """
    build = build_rhs_eval(
        solver_stem,
        enable_fCCZ4=enable_fCCZ4,
        fd_order=fd_order,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        CoordSystem=CoordSystem,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    formulation = "fCCZ4" if enable_fCCZ4 else "BSSN"
    subdirectory = "generated/src/rhs_eval"
    includes = [f"{solver_stem}_defines.h"]
    cfunc_type = "void"
    block_name = f"{solver_stem}_{RHS_EVAL_BLOCK_SUFFIX}"
    block_desc = (
        f"Per-block direct-FD {formulation} RHS ({len(build.evol_order)} fields)."
    )
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=block_desc,
        cfunc_type=cfunc_type,
        name=block_name,
        params=build.block_params,
        body=build.block_body,
    )
    roles.set_CFunction_role(block_name, "rhs_eval_block")
    roles.set_CFunction_codeparameters(block_name, build.used_codeparameters)
    all_blocks_name = f"{solver_stem}_{RHS_EVAL_ALL_BLOCKS_SUFFIX}"
    all_blocks_desc = f"All-block direct-FD {formulation} RHS (NRPy block loop)."
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=all_blocks_desc,
        cfunc_type=cfunc_type,
        name=all_blocks_name,
        params=build.all_blocks_params,
        body=build.all_blocks_body,
    )
    roles.set_CFunction_role(all_blocks_name, "rhs_eval")
    roles.set_CFunction_codeparameters(all_blocks_name, build.used_codeparameters)
    flat_block_name = f"{solver_stem}_{RHS_EVAL_FLAT_BLOCK_SUFFIX}"
    flat_block_desc = "LTS flat-block adapter (same numerical body, flat layout)."
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=flat_block_desc,
        cfunc_type=cfunc_type,
        name=flat_block_name,
        params=build.flat_block_params,
        body=build.flat_block_body,
    )
    roles.set_CFunction_role(flat_block_name, "rhs_eval_flat_block")
    roles.set_CFunction_codeparameters(flat_block_name, build.used_codeparameters)
    return build


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
