# nrpy/infrastructures/Dendro/general_relativity/generation_parameters.py
"""
Validate generation choices required by Dendro GR applications.

Dendro registers no NRPy parameters of its own. The scalar alias is the core
constant :data:`nrpy.grid.DENDRO_SCALAR_TYPE`, and Kreiss-Oliger dissipation is
a per-call builder argument. This module owns the qualified conformal-factor
representations required by the current GR initial-data paths.

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
    Validate the GR generation parameters against qualified values.

    :raises ValueError: If any parameter holds a disallowed value.

    Doctests:
    >>> import nrpy.equations.general_relativity.BSSN_quantities  # noqa: F401
    >>> for conformal_factor in ("chi", "W"):
    ...     par.set_parval_from_str("EvolvedConformalFactor_cf", conformal_factor)
    ...     validate_generation_parameters()
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
