# nrpy/infrastructures/Dendro/CFunction_roles.py
"""
Dendro CFunction roles, plus the padding, upwind-control and registry-order records the emitters read.

A Dendro kernel is an ordinary NRPy CFunction plus scheduling metadata keyed
by the registered function name.  The sidecar (stored in
``par.glb_extras_dict["Dendro"]["CFunction_roles"]``) contains no function
body, signature, parameter default, field declaration, or source path:
the CFunction registry is the only body/signature store.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Any, Dict, Sequence, Tuple, cast

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par


def registered_evol_order() -> Tuple[str, ...]:
    """
    Return the registered EVOL gridfunction names in NRPy list order.

    The NRPy gridfunction registry is the sole ordering and naming authority;
    every builder reads the order from here so no second ordering can appear.
    An empty EVOL set raises: a builder that reached this point with nothing
    registered would emit a well-formed kernel with an empty body rather than
    fail, which is the one failure a generated solver cannot report.

    :return: The ordered EVOL names.
    :raises ValueError: If no EVOL gridfunction is registered.

    Doctests:
    >>> gri.glb_gridfcs_dict.clear()
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> try:
    ...     registered_evol_order()
    ... except ValueError as error:
    ...     print(error)
    No EVOL gridfunction is registered; register the evolved state before building a Dendro kernel.
    >>> _ = gri.register_gridfunctions(["bXX", "aYY"], group="EVOL")
    >>> registered_evol_order()
    ('aYY', 'bXX')
    """
    evol, _auxevol, _diag, _aux = gri.GridFunction.gridfunction_lists()
    if not evol:
        raise ValueError(
            "No EVOL gridfunction is registered; register the evolved state "
            "before building a Dendro kernel."
        )
    return tuple(evol)


def registered_auxevol_order() -> Tuple[str, ...]:
    """
    Return the registered AUXEVOL gridfunction names in NRPy list order.

    :return: The ordered AUXEVOL names.
    """
    _evol, auxevol, _diag, _aux = gri.GridFunction.gridfunction_lists()
    return tuple(auxevol)


def registered_diag_order() -> Tuple[str, ...]:
    """
    Return the registered DIAG gridfunction names in NRPy list order.

    :return: The ordered DIAG names.
    """
    _evol, _auxevol, diag, _aux = gri.GridFunction.gridfunction_lists()
    return tuple(diag)


def registered_aux_order() -> Tuple[str, ...]:
    """
    Return the registered AUX gridfunction names in NRPy list order.

    :return: The ordered AUX names.
    """
    _evol, _auxevol, _diag, aux = gri.GridFunction.gridfunction_lists()
    return tuple(aux)


def _extras() -> Dict[str, Any]:
    """
    Return (creating if needed) the Dendro section of the NRPy extras dict.

    :return: The mutable Dendro extras dictionary.
    """
    return par.glb_extras_dict.setdefault("Dendro", {})


def set_required_padding(padding: int) -> None:
    """
    Record the ghost points the registered kernels need on every axis.

    The right-hand-side builder derives this from the derivative operators it
    actually emitted.  Recording it beside the role metadata keeps the value
    in the NRPy registries, so the emitters read it at the point of use
    instead of having it threaded through every caller.  It is one number, not
    one per axis, for the reason
    :func:`nrpy.infrastructures.Dendro.block_kernel_helpers.padding_from_derivative_operators`
    records.

    :param padding: Ghost points required on every axis.
    """
    _extras()["required_padding"] = int(padding)


def required_padding() -> int:
    """
    Return the recorded ghost points required on every axis.

    :return: The recorded padding.
    :raises ValueError: If no kernel has recorded a padding requirement.
    """
    padding = _extras().get("required_padding")
    if padding is None:
        raise ValueError(
            "No Dendro kernel has recorded a padding requirement; register the "
            "right-hand-side CFunctions before emitting the project."
        )
    return cast(int, padding)


def set_upwind_control_fields(names: Tuple[str, ...]) -> None:
    """
    Record the EVOL fields that drive the emitted upwind stencil selection.

    The right-hand-side builder derives these from the shared expression
    factory's upwind control vector.  The state header renders their registry
    positions, so recording them here keeps the emitted table and the emitted
    kernel derived from one value.

    :param names: Exact registered EVOL names, in registry order.
    """
    _extras()["upwind_control_fields"] = tuple(names)


def upwind_control_fields() -> Tuple[str, ...]:
    """
    Return the recorded upwind control field names.

    :return: The recorded EVOL names, in registry order.
    :raises ValueError: If no kernel has recorded an upwind control set.
    """
    names = _extras().get("upwind_control_fields")
    if names is None:
        raise ValueError(
            "No Dendro kernel has recorded an upwind control set; register the "
            "right-hand-side CFunctions before emitting the state header."
        )
    return cast(Tuple[str, ...], names)


def set_CFunction_role(name: str, role: str) -> None:
    """
    Record one registered CFunction's Dendro scheduling role.

    Builders call :func:`nrpy.c_function.register_CFunction` themselves, as every
    established infrastructure does, and then record the role here.  The role is
    non-authoritative scheduling metadata, so the host-adapter emitters can ask
    for "the all-block RHS entry point" rather than taking a dozen name
    arguments.

    :param name: The registered CFunction name.
    :param role: Scheduling role, e.g. ``"rhs_eval_block"``.
    :raises ValueError: If no CFunction of that name is registered, which would
        leave the sidecar naming something the registry does not carry.

    Doctests:
    >>> cfc.CFunction_dict.clear()
    >>> par.glb_extras_dict.pop("Dendro", None) and None
    >>> cfc.register_CFunction(
    ...     desc="Per-block RHS.",
    ...     name="bssn_rhs_eval_block",
    ...     params="int n",
    ...     body="(void)n;",
    ... )
    >>> set_CFunction_role("bssn_rhs_eval_block", "rhs_eval_block")
    >>> CFunction_name_for_role("rhs_eval_block")
    'bssn_rhs_eval_block'
    >>> try:
    ...     set_CFunction_role("not_registered", "rhs_eval_block")
    ... except ValueError as error:
    ...     print(error)
    Cannot record a Dendro role for 'not_registered': no CFunction of that name is registered.
    >>> try:
    ...     CFunction_name_for_role("constraints_eval")
    ... except ValueError as error:
    ...     print(error)
    Expected exactly one registered CFunction with Dendro role 'constraints_eval', found [].
    """
    if name not in cfc.CFunction_dict:
        raise ValueError(
            f"Cannot record a Dendro role for {name!r}: no CFunction of that "
            "name is registered."
        )
    _CFunction_roles()[name] = role


def set_CFunction_codeparameters(name: str, names: Sequence[str]) -> None:
    """
    Record the CodeParameters one registered CFunction's signature forwards.

    The builder computed this set exactly, from the expression free symbols, so
    it is recorded here rather than recovered by splitting the emitted signature
    on commas: reading the free symbols the equations contain is exact, scanning
    emitted C text for names is not.

    :param name: The registered CFunction name.
    :param names: CodeParameter names, in signature order.
    :raises ValueError: If no CFunction of that name is registered.

    Doctests:
    >>> cfc.CFunction_dict.clear()
    >>> par.glb_extras_dict.pop("Dendro", None) and None
    >>> cfc.register_CFunction(
    ...     desc="k", name="bssn_k", params="int n", body="(void)n;"
    ... )
    >>> set_CFunction_codeparameters("bssn_k", ("eta", "kappa1"))
    >>> CFunction_codeparameters("bssn_k")
    ('eta', 'kappa1')
    >>> CFunction_codeparameters("bssn_unrecorded")
    ()
    """
    if name not in cfc.CFunction_dict:
        raise ValueError(
            f"Cannot record CodeParameters for {name!r}: no CFunction of that "
            "name is registered."
        )
    _CFunction_codeparameters()[name] = tuple(names)


def _CFunction_codeparameters() -> Dict[str, Tuple[str, ...]]:
    """
    Return the mutable CodeParameter sidecar, creating it on first use.

    :return: Mapping from registered CFunction name to CodeParameter names.
    """
    extras = _extras()
    table = extras.setdefault("CFunction_codeparameters", {})
    return cast(Dict[str, Tuple[str, ...]], table)


def CFunction_codeparameters(name: str) -> Tuple[str, ...]:
    """
    Return the CodeParameters one registered CFunction's signature forwards.

    :param name: The registered CFunction name.
    :return: CodeParameter names in signature order, empty when none were
        recorded.
    """
    return _CFunction_codeparameters().get(name, ())


def _CFunction_roles() -> Dict[str, str]:
    """
    Return the Dendro role sidecar keyed by registered CFunction name.

    :return: Mapping of CFunction name to its scheduling role.
    """
    return cast(Dict[str, str], _extras().setdefault("CFunction_roles", {}))


def CFunction_name_for_role(role: str) -> str:
    """
    Return the registered CFunction name that carries one Dendro role.

    The host-adapter emitters need the name of, say, the all-block RHS entry
    point.  The role sidecar already records it, so they read it from there
    rather than take a dozen name arguments or rebuild the naming convention.

    :param role: Dendro scheduling role, e.g. ``"rhs_eval_block"``.
    :return: The registered CFunction name carrying that role.
    :raises ValueError: If no registered CFunction, or more than one, carries
        the role.
    """
    matches = sorted(
        name for name, recorded in _CFunction_roles().items() if recorded == role
    )
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one registered CFunction with Dendro role {role!r}, "
            f"found {matches}."
        )
    return matches[0]


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
