# nrpy/infrastructures/Dendro/access_capture.py
"""
Scoped capture of Dendro gridfunction memory accesses.

During NRPy code generation, :class:`gri.DendroGridFunction` records every
one-point access (gridfunction name plus signed offsets) into the active
capture context.  The capture context is keyed by the CFunction being built so
that the padding required by that function can be derived from the exact
accesses it emits, rather than from a second stencil model.

The recorded data lives in ``par.glb_extras_dict["Dendro"]`` and contains only
names and signed integer offsets: no function bodies, no field declarations,
and no parallel state.
"""

from contextlib import contextmanager
from contextvars import ContextVar
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple

import nrpy.params as par

# The CFunction name of the currently active capture context, if any.
_ACTIVE_CAPTURE: ContextVar[Optional[str]] = ContextVar(
    "dendro_active_capture", default=None
)

# Single authority for access records: par.glb_extras_dict["Dendro"]
# ["accesses_by_cfunction"].  There is no module-global mirror; every helper
# below reads/writes the extras store directly (whitepaper sections 8.3/9.4).


def _dendro_extras() -> Dict[str, Any]:
    """
    Return (creating if needed) the Dendro section of the NRPy extras dict.

    :return: The mutable Dendro extras dictionary.
    """
    return par.glb_extras_dict.setdefault("Dendro", {})


def _access_store() -> Dict[str, Set[Tuple[str, int, int, int]]]:
    """
    Return (creating if needed) the single-authority access record store.

    :return: Mapping of CFunction name to its recorded access set.
    """
    extras = _dendro_extras()
    store = extras.get("accesses_by_cfunction")
    if not isinstance(store, dict):
        store = {}
        extras["accesses_by_cfunction"] = store
    return store


def record_dendro_access(
    gridfunction_name: str, i0_offset: int, i1_offset: int, i2_offset: int
) -> None:
    """
    Record one gridfunction memory access for the active capture context.

    If no capture context is active the access is simply not recorded; Dendro
    gridfunction reads remain fully functional outside a capture.

    :param gridfunction_name: Exact registered NRPy gridfunction name.
    :param i0_offset: Signed offset in the fastest (x) direction.
    :param i1_offset: Signed offset in the middle (y) direction.
    :param i2_offset: Signed offset in the slowest (z) direction.
    :raises ValueError: If two capture contexts are active at once.
    """
    active = _ACTIVE_CAPTURE.get()
    if active is None:
        return
    store = _access_store()
    if active not in store:
        raise ValueError(
            f"Active Dendro capture {active!r} has no recorded access set."
        )
    store[active].add((gridfunction_name, i0_offset, i1_offset, i2_offset))


@contextmanager
def capture_gridfunction_accesses(cfunction_name: str) -> Iterator[None]:
    """
    Capture all Dendro gridfunction accesses while building one CFunction.

    On success the recorded accesses remain available through
    :func:`get_captured_accesses`.  On failure the partial record is removed.

    :param cfunction_name: Name of the CFunction whose code is being built.
    :yields: Nothing (the context is a null context manager).
    :raises ValueError: If accesses were already captured for this CFunction
        or a capture context is already active.
    :raises Exception: If the body raises, the partial record is removed and
        the exception is re-raised.
    """
    if _ACTIVE_CAPTURE.get() is not None:
        raise ValueError(
            "Nested Dendro access captures are not supported; "
            f"a capture for {_ACTIVE_CAPTURE.get()!r} is already active."
        )
    store = _access_store()
    if cfunction_name in store:
        raise ValueError(f"Accesses already captured for {cfunction_name}")
    store[cfunction_name] = set()
    token = _ACTIVE_CAPTURE.set(cfunction_name)
    try:
        yield
    except Exception:
        store.pop(cfunction_name, None)
        raise
    finally:
        _ACTIVE_CAPTURE.reset(token)


def record_empty_capture(cfunction_name: str) -> None:
    """
    Record an explicit empty capture for a neighbor-free numerical CFunction.

    Writers such as the Minkowski fill touch only the current point, so their
    correct access set is empty and their padding is (0, 0, 0).  Registering
    the empty set explicitly satisfies the freeze leaf-capture gate (section
    8.6) without inventing reads.

    :param cfunction_name: CFunction name to record as neighbor-free.
    :raises ValueError: If a capture already exists for this CFunction.
    """
    store = _access_store()
    if cfunction_name in store:
        raise ValueError(f"Accesses already captured for {cfunction_name}")
    store[cfunction_name] = set()


def get_captured_accesses(cfunction_name: str) -> Tuple[Tuple[str, int, int, int], ...]:
    """
    Return the deterministic (sorted) accesses captured for one CFunction.

    :param cfunction_name: Name of the CFunction whose accesses were captured.
    :return: Tuple of (gridfunction name, i0, i1, i2) access tuples, sorted.
    :raises ValueError: If no capture exists for this CFunction.
    """
    store = _access_store()
    if cfunction_name not in store:
        raise ValueError(f"No accesses captured for {cfunction_name}")
    return tuple(sorted(store[cfunction_name]))


def accessed_gridfunction_names(cfunction_name: str) -> Tuple[str, ...]:
    """
    Return the distinct gridfunction names one CFunction actually reads.

    Binding every registered pointer in a kernel that reads a subset trips
    ``-Werror=unused-variable``, so the builders bind exactly the fields the
    capture recorded.  The names are sorted for determinism; callers that need
    registry order re-order them against the registry.

    :param cfunction_name: Name of the CFunction whose accesses were captured.
    :return: Sorted tuple of distinct gridfunction names.
    """
    return tuple(
        sorted({name for name, _i0, _i1, _i2 in get_captured_accesses(cfunction_name)})
    )


def required_padding(cfunction_name: str) -> Tuple[int, int, int]:
    """
    Reduce the captured accesses of one CFunction to per-axis padding.

    The padding along each axis is the maximum absolute signed offset reached
    by any access of that CFunction.

    :param cfunction_name: Name of the CFunction whose accesses were captured.
    :return: Tuple of (padding_x, padding_y, padding_z).

    Doctests:
    >>> import nrpy.grid as gri
    >>> import nrpy.params as par
    >>> import nrpy.infrastructures.Dendro.generation_parameters  # noqa: F401
    >>> from nrpy.infrastructures.Dendro import access_capture as cap

    The capture is exercised through the imported package module, not through
    this file's ``__main__`` copy: ``nrpy.grid`` records into the package
    module, and the two hold separate context variables.

    >>> cap._reset_state()
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> gri.glb_gridfcs_dict.clear()
    >>> _ = gri.register_gridfunctions("probe")

    Whitepaper section 8.3: a read records its exact signed offsets, and
    section 8.6 reduces them to the padding that read requires.

    >>> with cap.capture_gridfunction_accesses("probe_kernel"):
    ...     _ = gri.glb_gridfcs_dict["probe"].read_gf_from_memory_Ccode_onept(2, 0, 0)
    ...     _ = gri.glb_gridfcs_dict["probe"].read_gf_from_memory_Ccode_onept(0, -3, 1)
    >>> cap.get_captured_accesses("probe_kernel")
    (('probe', 0, -3, 1), ('probe', 2, 0, 0))
    >>> cap.required_padding("probe_kernel")
    (2, 3, 1)
    >>> cap.accessed_gridfunction_names("probe_kernel")
    ('probe',)

    A neighbour-free kernel records an explicit empty capture, which reduces
    to zero padding rather than to a missing record.

    >>> cap.record_empty_capture("pointwise_kernel")
    >>> cap.required_padding("pointwise_kernel")
    (0, 0, 0)
    >>> cap._reset_state()
    >>> gri.glb_gridfcs_dict.clear()
    """
    accesses = get_captured_accesses(cfunction_name)
    if not accesses:
        return (0, 0, 0)
    padding = [0, 0, 0]
    for _name, i0, i1, i2 in accesses:
        padding[0] = max(padding[0], abs(i0))
        padding[1] = max(padding[1], abs(i1))
        padding[2] = max(padding[2], abs(i2))
    return (padding[0], padding[1], padding[2])


def _dendro_extras_state() -> Tuple[Dict[str, Any], Any]:
    """
    Return the pieces of Dendro state owned by this module, for snapshots.

    :return: (extras dict, active capture token value).
    """
    return (_dendro_extras(), _ACTIVE_CAPTURE.get())


def _reset_state() -> None:
    """
    Remove all Dendro access-capture state.

    Intended for generation transactions: called before capturing and, on
    rollback, to restore the pre-transaction (empty) state.
    """
    par.glb_extras_dict.get("Dendro", {}).pop("accesses_by_cfunction", None)


def captured_cfunction_names() -> List[str]:
    """
    Return the names of all CFunctions with recorded accesses, sorted.

    :return: Sorted list of CFunction names with captures.
    """
    return sorted(_access_store())


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
