# nrpy/infrastructures/Dendro/registration.py
"""
Dendro CFunction registration with non-authoritative role metadata.

A Dendro kernel is an ordinary NRPy CFunction plus scheduling metadata keyed
by the registered function name.  The sidecar (stored in
``par.glb_extras_dict["Dendro"]["CFunction_roles"]``) contains no function
body, signature, parameter default, field declaration, or source path:
the CFunction registry is the only body/signature store.
"""

from typing import Any, Dict, Tuple, cast

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par

# Metadata state machine: registration is open until freeze() seals it.
# Stored in a mutable module dict so toggling it never needs a global.
_registration_state: Dict[str, bool] = {"open": True}


def registered_evol_order() -> Tuple[str, ...]:
    """
    Return the registered EVOL gridfunction names in NRPy list order.

    The NRPy gridfunction registry is the sole ordering and naming authority
    (whitepaper section 5.1); every builder reads the order from here so no
    second ordering can appear.

    :return: The ordered EVOL names.
    """
    evol, _auxevol, _diag, _aux = gri.GridFunction.gridfunction_lists()
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


def assert_dendro_registration_open() -> None:
    """
    Raise if Dendro registration has been sealed by a freeze.

    :raises RuntimeError: If a Dendro role is registered after freeze.
    """
    if not _registration_state["open"]:
        raise RuntimeError(
            "Dendro registration is sealed; freeze has already run. "
            "Start a new generation (fresh process) to register CFunctions."
        )


def _dendro_extras() -> Dict[str, Any]:
    """
    Return (creating if needed) the Dendro section of the NRPy extras dict.

    :return: The mutable Dendro extras dictionary.
    """
    return par.glb_extras_dict.setdefault("Dendro", {})


def _seal_registration(sealed: bool) -> None:
    """
    Set the registration-open flag (used by freeze and transactions).

    :param sealed: True to seal (block) role registration, False to open it.
    """
    _registration_state["open"] = not sealed


def set_registration_open(open_flag: bool) -> None:
    """
    Explicitly open or seal Dendro registration.

    :param open_flag: True to allow role registration, False to seal.
    """
    _seal_registration(not open_flag)


def register_dendro_CFunction(
    *,
    role: str,
    entry_point: bool = False,
    calls: Tuple[str, ...] = (),
    lifecycle_hook: str = "",
    **cfunction_kwargs: Any,
) -> None:
    """
    Register a CFunction in the NRPy registry and record its Dendro role.

    The CFunction body and signature are registered exactly once in
    ``cfc.CFunction_dict``; the role sidecar only records scheduling
    metadata (role, entry point, callees, lifecycle hook).

    :param role: Non-authoritative scheduling role (e.g., ``"rhs_block"``).
    :param entry_point: Whether Dendro invokes this CFunction directly.
    :param calls: Names of registered CFunctions (or allowed host APIs)
        that this CFunction invokes.
    :param lifecycle_hook: Dendro lifecycle point at which this CFunction
        runs (e.g., ``"rhs"``, ``"post_timestep"``), or empty.
    :param cfunction_kwargs: Arguments forwarded to
        :func:`nrpy.c_function.register_CFunction`.
    :raises ValueError: If a Dendro role was already registered for this
        CFunction name.
    """
    assert_dendro_registration_open()
    cfc.register_CFunction(**cfunction_kwargs)
    name = str(cfunction_kwargs["name"])
    roles = _dendro_extras().setdefault("CFunction_roles", {})
    if name in roles:
        raise ValueError(f"Duplicate Dendro role for {name}")
    roles[name] = {
        "role": role,
        "entry_point": bool(entry_point),
        "calls": tuple(sorted(set(calls), key=lambda c: (c.lower(), c))),
        "lifecycle_hook": lifecycle_hook,
    }


def get_CFunction_roles() -> Dict[str, Dict[str, Any]]:
    """
    Return the Dendro role sidecar keyed by registered CFunction name.

    :return: Mapping of CFunction name to its scheduling metadata.
    """
    return cast(
        Dict[str, Dict[str, Any]], _dendro_extras().setdefault("CFunction_roles", {})
    )


def _clear_role_state() -> None:
    """Remove the Dendro role sidecar (used by generation transactions)."""
    _dendro_extras().pop("CFunction_roles", None)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
