# nrpy/infrastructures/Dendro/generation_parameters.py
"""
Dendro target generation parameters.

These are ordinary NRPy parameters registered through ``par.register_param``.
They hold the *non-scientific* generation choices for lowering a registered NRPy
environment into a Dendro-GR module.  No parallel configuration object is kept:
the values live only in the NRPy parameter registry and are read back from it.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, Tuple

import nrpy.params as par

# Register the Dendro target generation choices.  register_param() is
# idempotent, so importing this module (or calling register again) is safe.
par.register_param(
    str,
    __name__,
    "Dendro_scalar_type",
    "DendroScalar",
    description="C/C++ scalar type alias used for all generated numerical code.",
)
par.register_param(
    bool,
    __name__,
    "Dendro_enable_KreissOliger_dissipation",
    True,
    description="Enable NRPy-owned Kreiss-Oliger dissipation in the canonical backend.",
)


_ALLOWED_VALUES: Dict[str, Tuple[object, ...]] = {
    # The registered `f_infinity` of `cf` is 1 regardless of the selected
    # representation, so only representations whose Minkowski value is 1 may
    # be generated: `phi` has Minkowski value 0 and would make the generated
    # initial data wrong.
    "EvolvedConformalFactor_cf": ("chi", "W"),
}


def validate_generation_parameters() -> None:
    """
    Validate the Dendro generation parameters against allowed values.

    Every registered Dendro generation parameter is validated here.  Only the
    qualified profiles are accepted; anything else fails generation
    instead of silently producing an unqualified configuration.

    :raises ValueError: If any parameter holds a disallowed value.

    Doctests:
    >>> par.set_parval_from_str("Dendro_scalar_type", "DendroScalar")
    >>> validate_generation_parameters()
    >>> par.set_parval_from_str("Dendro_scalar_type", "not an identifier")
    >>> try:
    ...     validate_generation_parameters()
    ... except ValueError as error:
    ...     print(error)
    Invalid Dendro_scalar_type: 'not an identifier'.
    >>> par.set_parval_from_str("Dendro_scalar_type", "DendroScalar")
    """
    scalar_type = par.parval_from_str("Dendro_scalar_type")
    if not isinstance(scalar_type, str) or not scalar_type.isidentifier():
        raise ValueError(f"Invalid Dendro_scalar_type: {scalar_type!r}.")
    for name, allowed in _ALLOWED_VALUES.items():
        # Shared NRPy parameters (e.g. EvolvedConformalFactor_cf) are
        # registered by the equation modules, so a profile that builds no
        # equations simply has nothing to validate for them.
        if name not in par.glb_params_dict:
            continue
        value = par.parval_from_str(name)
        if value not in allowed:
            raise ValueError(f"Unsupported {name}={value!r}; allowed: {list(allowed)}.")
    enable_ko = par.parval_from_str("Dendro_enable_KreissOliger_dissipation")
    if not isinstance(enable_ko, bool):
        raise ValueError(
            f"Invalid Dendro_enable_KreissOliger_dissipation: {enable_ko!r}; expected bool."
        )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
