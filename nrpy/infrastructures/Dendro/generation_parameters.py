# nrpy/infrastructures/Dendro/generation_parameters.py
"""
Validation of the Dendro generation choices.

Dendro registers no NRPy parameters of its own.  The scalar alias is the core
constant :data:`nrpy.grid.DENDRO_SCALAR_TYPE`, hardcoded as the three sibling
gridfunction classes hardcode theirs, and Kreiss-Oliger dissipation is a
per-call builder argument as it is in BHaH and ETLegacy.  What remains is the
qualified-value check on the core parameters a Dendro profile depends on, so an
unqualified configuration fails generation instead of being emitted.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, Tuple

import nrpy.params as par

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

    Only the qualified profiles are accepted; anything else fails generation
    instead of silently producing an unqualified configuration.  The scalar
    alias is not validated here because it is no longer a parameter: it is the
    constant :data:`nrpy.grid.DENDRO_SCALAR_TYPE`, hardcoded exactly as the
    three sibling gridfunction classes hardcode theirs.

    :raises ValueError: If any parameter holds a disallowed value.

    Doctests:
    >>> import nrpy.equations.general_relativity.BSSN_quantities  # noqa: F401
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    >>> validate_generation_parameters()
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "phi")
    >>> try:
    ...     validate_generation_parameters()
    ... except ValueError as error:
    ...     print(error)
    Unsupported EvolvedConformalFactor_cf 'phi'; qualified values: chi, W.
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    """
    for name, allowed in _ALLOWED_VALUES.items():
        value = par.parval_from_str(name)
        if value not in allowed:
            raise ValueError(
                f"Unsupported {name} {value!r}; qualified values: "
                + ", ".join(str(item) for item in allowed)
                + "."
            )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
