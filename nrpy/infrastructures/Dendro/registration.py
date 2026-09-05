# nrpy/infrastructures/Dendro/registration.py
"""
Dendro CFunction registration with non-authoritative role metadata.

A Dendro kernel is an ordinary NRPy CFunction plus scheduling metadata keyed
by the registered function name.  The sidecar (stored in
``par.glb_extras_dict["Dendro"]["CFunction_roles"]``) contains no function
body, signature, parameter default, field declaration, or source path:
the CFunction registry is the only body/signature store.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Any, Dict, Tuple, cast

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
    >>> import nrpy.infrastructures.Dendro.generation_parameters  # noqa: F401
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


def _Dendro_extras() -> Dict[str, Any]:
    """
    Return (creating if needed) the Dendro section of the NRPy extras dict.

    :return: The mutable Dendro extras dictionary.
    """
    return par.glb_extras_dict.setdefault("Dendro", {})


def set_required_padding(padding: Tuple[int, int, int]) -> None:
    """
    Record the ghost points the registered kernels need on each axis.

    The right-hand-side builder derives this from the derivative operators it
    actually emitted.  Recording it beside the role metadata keeps the value
    in the NRPy registries, so the emitters read it at the point of use
    instead of having it threaded through every caller.

    :param padding: Ghost points required on the x, y and z axes.
    """
    _Dendro_extras()["required_padding"] = tuple(int(p) for p in padding)


def required_padding() -> Tuple[int, int, int]:
    """
    Return the recorded ghost points required on the x, y and z axes.

    :return: The recorded (px, py, pz).
    :raises ValueError: If no kernel has recorded a padding requirement.
    """
    padding = _Dendro_extras().get("required_padding")
    if padding is None:
        raise ValueError(
            "No Dendro kernel has recorded a padding requirement; register the "
            "right-hand-side CFunctions before emitting the project."
        )
    return cast(Tuple[int, int, int], padding)


def set_upwind_control_fields(names: Tuple[str, ...]) -> None:
    """
    Record the EVOL fields that drive the emitted upwind stencil selection.

    The right-hand-side builder derives these from the shared expression
    factory's upwind control vector.  The state header renders their registry
    positions, so recording them here keeps the emitted table and the emitted
    kernel derived from one value.

    :param names: Exact registered EVOL names, in registry order.
    """
    _Dendro_extras()["upwind_control_fields"] = tuple(names)


def upwind_control_fields() -> Tuple[str, ...]:
    """
    Return the recorded upwind control field names.

    :return: The recorded EVOL names, in registry order.
    :raises ValueError: If no kernel has recorded an upwind control set.
    """
    names = _Dendro_extras().get("upwind_control_fields")
    if names is None:
        raise ValueError(
            "No Dendro kernel has recorded an upwind control set; register the "
            "right-hand-side CFunctions before emitting the state header."
        )
    return cast(Tuple[str, ...], names)


def register_Dendro_CFunction(*, role: str, **cfunction_kwargs: Any) -> None:
    """
    Register a CFunction in the NRPy registry and record its Dendro role.

    The CFunction body and signature are registered exactly once in
    ``cfc.CFunction_dict``; the sidecar records only the scheduling role, so
    the host-adapter emitters can ask for "the all-block RHS entry point"
    instead of taking a dozen name arguments.

    :param role: Non-authoritative scheduling role (e.g., ``"rhs_block"``).
    :param cfunction_kwargs: Arguments forwarded to
        :func:`nrpy.c_function.register_CFunction`.

    Doctests:
    >>> cfc.CFunction_dict.clear()
    >>> par.glb_extras_dict.pop("Dendro", None) and None
    >>> register_Dendro_CFunction(
    ...     role="rhs_block", desc="Per-block RHS.", name="bssn_rhs_block",
    ...     params="int n", body="(void)n;")
    >>> CFunction_name_for_role("rhs_block")
    'bssn_rhs_block'
    >>> "bssn_rhs_block" in cfc.CFunction_dict
    True
    >>> try:
    ...     CFunction_name_for_role("diagnostics")
    ... except ValueError as error:
    ...     print(error)
    Expected exactly one registered CFunction with Dendro role 'diagnostics', found [].
    """
    # ``cfc.register_CFunction`` already rejects a duplicate name, so the
    # sidecar cannot acquire two entries for one CFunction.
    cfc.register_CFunction(**cfunction_kwargs)
    _CFunction_roles()[str(cfunction_kwargs["name"])] = role


def _CFunction_roles() -> Dict[str, str]:
    """
    Return the Dendro role sidecar keyed by registered CFunction name.

    :return: Mapping of CFunction name to its scheduling role.
    """
    return cast(Dict[str, str], _Dendro_extras().setdefault("CFunction_roles", {}))


def CFunction_name_for_role(role: str) -> str:
    """
    Return the registered CFunction name that carries one Dendro role.

    The host-adapter emitters need the name of, say, the all-block RHS entry
    point.  The role sidecar already records it, so they read it from there
    rather than take a dozen name arguments or rebuild the naming convention.

    :param role: Dendro scheduling role, e.g. ``"rhs_block"``.
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
