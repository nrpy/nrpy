# nrpy/infrastructures/Dendro/transaction.py
"""
Generation transaction for NRPy Dendro generation.

NRPy uses process-global registries.  The production command runs one
generation profile in a fresh Python process; this transaction adds rollback
and test isolation on top: it snapshots the relevant global registries,
verifies a safe clean starting state, and restores everything on failure.
"""

import copy
import hashlib
import json
from contextlib import contextmanager
from typing import Any, Dict, Iterator, List, Tuple

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.Dendro import registration as dendro_reg

# Ordered generation phases; transitions are forward-only.
_GENERATION_PHASES = (
    "CREATED",
    "CONFIGURED",
    "EXPRESSIONS_BUILT",
    "GRIDFUNCTIONS_AND_CODEPARAMETERS_REGISTERED",
    "CFUNCTIONS_REGISTERED",
    "ROLE_AND_ACCESS_METADATA_COMPLETE",
    "FROZEN",
    "EMITTED",
    "VERIFIED",
)


class GenerationTransaction:
    """
    Snapshot and rollback of the NRPy registries for one generation run.

    :param captured_state: State captured at construction (used internally).
    """

    def __init__(self, captured_state: Dict[str, Any]) -> None:
        self._params_values: Dict[str, Any] = dict(captured_state["params_values"])
        self._params_keys: List[str] = list(captured_state["params_keys"])
        self._code_params_snapshot: Dict[str, Any] = dict(
            captured_state["code_params_snapshot"]
        )
        self._gridfcs_keys: List[str] = list(captured_state["gridfcs_keys"])
        self._cfunctions_keys: List[str] = list(captured_state["cfunctions_keys"])
        self._dendro_extras_snapshot: Dict[str, Any] = copy.deepcopy(
            captured_state["dendro_extras_snapshot"]
        )
        self._registration_was_open: bool = bool(
            captured_state["registration_was_open"]
        )
        self._phase_index: int = 0

    @classmethod
    def capture(cls) -> "GenerationTransaction":
        """
        Capture the current state of the relevant NRPy registries.

        :return: A GenerationTransaction holding the captured state.
        """
        return cls(
            {
                "params_values": {
                    name: param.value for name, param in par.glb_params_dict.items()
                },
                "params_keys": list(par.glb_params_dict),
                "code_params_snapshot": {
                    name: copy.deepcopy(cp.defaultvalue)
                    for name, cp in par.glb_code_params_dict.items()
                },
                "gridfcs_keys": list(gri.glb_gridfcs_dict),
                "cfunctions_keys": list(cfc.CFunction_dict),
                "dendro_extras_snapshot": copy.deepcopy(
                    par.glb_extras_dict.get("Dendro", {})
                ),
                "registration_was_open": bool(dendro_reg._registration_state["open"]),
            }
        )

    @property
    def phase(self) -> str:
        """
        Return the current generation phase name.

        :return: The current phase.
        """
        return _GENERATION_PHASES[self._phase_index]

    def advance_to(self, phase: str) -> None:
        """
        Advance the generation phase machine to a later phase.

        :param phase: Target phase name.
        :raises ValueError: If the target phase is unknown or not later than
            the current phase.
        """
        if phase not in _GENERATION_PHASES:
            raise ValueError(f"Unknown generation phase: {phase}")
        target = _GENERATION_PHASES.index(phase)
        if target <= self._phase_index:
            raise ValueError(
                f"Cannot move from {self.phase} to {phase}; phases are forward-only."
            )
        self._phase_index = target

    def require_clean(self, allowlist: Tuple[str, ...] = ()) -> None:
        """
        Require a safe clean gridfunction and CFunction registry.

        :param allowlist: CFunction names permitted to pre-exist (e.g., from
            imported modules that register at import time).
        :raises ValueError: If the registry is not clean.
        """
        if gri.glb_gridfcs_dict:
            raise ValueError(
                "Expected an empty gridfunction registry for a Dendro "
                f"generation, found: {sorted(gri.glb_gridfcs_dict)}"
            )
        unexpected = [name for name in cfc.CFunction_dict if name not in allowlist]
        if unexpected:
            raise ValueError(
                "Expected an empty CFunction registry for a Dendro "
                f"generation, found: {sorted(unexpected)}"
            )

    def rollback(self) -> None:
        """
        Restore the registries to the captured pre-transaction state.

        Parameter *values* are restored (import-time registrations are
        reused), and any gridfunction, CodeParameter, or CFunction entry
        created during the transaction is removed.  Dendro-owned extras
        (role sidecar, access captures) and the registration seal are reset.
        """
        for name, value in self._params_values.items():
            if name in par.glb_params_dict:
                par.glb_params_dict[name].value = value
        for name in list(par.glb_params_dict):
            if name not in self._params_keys:
                del par.glb_params_dict[name]
        for name in list(par.glb_code_params_dict):
            if name not in self._code_params_snapshot:
                del par.glb_code_params_dict[name]
        for name, defaultvalue in self._code_params_snapshot.items():
            if name in par.glb_code_params_dict:
                par.glb_code_params_dict[name].defaultvalue = copy.deepcopy(
                    defaultvalue
                )
        for name in list(gri.glb_gridfcs_dict):
            if name not in self._gridfcs_keys:
                del gri.glb_gridfcs_dict[name]
        for name in list(cfc.CFunction_dict):
            if name not in self._cfunctions_keys:
                del cfc.CFunction_dict[name]
        if self._dendro_extras_snapshot:
            par.glb_extras_dict["Dendro"] = copy.deepcopy(self._dendro_extras_snapshot)
        else:
            par.glb_extras_dict.pop("Dendro", None)
        dendro_reg.set_registration_open(self._registration_was_open)
        self._phase_index = 0


@contextmanager
def dendro_generation_transaction(
    require_clean: bool = True,
    allowlist: Tuple[str, ...] = (),
) -> Iterator[GenerationTransaction]:
    """
    Context manager for one Dendro generation profile.

    On any exception, the transaction rolls back the registries and the
    exception is re-raised.

    :param require_clean: Require empty gridfunction/CFunction registries
        (outside the allowlist) at transaction start.
    :param allowlist: CFunction names permitted to pre-exist.
    :yields: The active GenerationTransaction.
    :raises Exception: On any exception the registries are rolled back and the
        exception is re-raised.

    Doctests:
    >>> import nrpy.grid as gri
    >>> import nrpy.params as par
    >>> import nrpy.infrastructures.Dendro.generation_parameters  # noqa: F401
    >>> from nrpy.infrastructures.Dendro import transaction as tx_mod
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> gri.glb_gridfcs_dict.clear()

    Whitepaper section 4.3: a failed generation leaves the registries exactly
    as it found them, so a later run in the same process is not poisoned.

    >>> before = tx_mod.registry_digest()
    >>> try:
    ...     with tx_mod.dendro_generation_transaction() as tx:
    ...         _ = gri.register_gridfunctions("scratch_probe")
    ...         raise RuntimeError("generation failed")
    ... except RuntimeError:
    ...     pass
    >>> "scratch_probe" in gri.glb_gridfcs_dict
    False
    >>> tx_mod.registry_digest() == before
    True

    Section 4.4: the phase machine is forward-only.

    >>> with tx_mod.dendro_generation_transaction() as tx:
    ...     tx.advance_to("FROZEN")
    ...     print(tx.phase)
    ...     try:
    ...         tx.advance_to("FROZEN")
    ...     except ValueError:
    ...         print("Backward or repeated phase rejected. Good.")
    FROZEN
    Backward or repeated phase rejected. Good.
    """
    tx = GenerationTransaction.capture()
    if require_clean:
        tx.require_clean(allowlist)
    try:
        yield tx
    except Exception:
        tx.rollback()
        raise


def registry_digest() -> str:
    """
    Compute a deterministic digest of the mutable NRPy registries.

    Used before and after export to prove the frozen snapshot still matches
    the mutable registries.

    :return: Hex SHA-256 digest of the canonical registry contents.
    """
    canonical = {
        "params": {
            name: par.glb_params_dict[name].value
            for name in sorted(par.glb_params_dict)
        },
        "code_params": {
            name: repr(par.glb_code_params_dict[name].defaultvalue)
            for name in sorted(par.glb_code_params_dict)
        },
        "gridfcs": sorted(gri.glb_gridfcs_dict),
        "cfunctions": {
            name: {
                "prototype": cfc.CFunction_dict[name].function_prototype,
                "body": cfc.CFunction_dict[name].full_function,
                "subdirectory": cfc.CFunction_dict[name].subdirectory,
            }
            for name in sorted(cfc.CFunction_dict)
        },
    }
    payload = json.dumps(canonical, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
