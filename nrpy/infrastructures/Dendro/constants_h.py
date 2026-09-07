# nrpy/infrastructures/Dendro/constants_h.py
"""
Emit the ``<stem>_constants.h`` header of generated compile-time constants.

The constants are the finite-difference order the kernels were generated at, the
ghost points their derivative operators reach, and whether those operators
include Kreiss-Oliger dissipation.  All three are read from the NRPy registries
or from the caller that built the kernels, rather than restated here, so the
generated solver cannot disagree with the kernels it was generated with.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List

import nrpy.params as par
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.header_guards import header_guard

BANNER = generated_file_banner()


def output_constants_h(
    solver_stem: str,
    solver_namespace: str,
    required_padding: int,
    enable_KreissOliger_dissipation: bool,
) -> str:
    """
    Emit the generated constants header.

    It carries the finite-difference order, the required block padding and the
    Kreiss-Oliger switch.  The padding is supplied by the kernel builder, which
    takes it from the widest reach of the derivative
    operators the emitted kernel actually contains.  It is not ``fd_order // 2``:
    the upwinded and Kreiss-Oliger families reach one point further than the
    centred ones, and a host that sized its ghost zones from the radius would
    have the kernel read past the end of a block.

    :param solver_stem: Lowercase formulation stem for the emitted header name.
    :param solver_namespace: Solver namespace, following Dendro's lowercase
        formulation habit (``namespace bssn``).
    :param required_padding: Ghost points required on every axis.
    :param enable_KreissOliger_dissipation: Whether the emitted kernel carries
        Kreiss-Oliger dissipation, passed by the caller that built the kernels.
        The generated padding self-test needs it: the Kreiss-Oliger family
        reaches one point past the centred radius, exactly as the upwinded one
        does, and no other emitted constant records its presence.
    :return: The complete C++ header text.
    :raises ValueError: If the padding is below the centred-plus-Kreiss-Oliger
        floor.  This is a floor, not the exact reach: an upwinded kernel reaches
        one point further still, and whether the kernel upwinds is recorded in
        the state header rather than here, so the generated ``test_padding``
        case is what checks the exact reach.

    Doctests:
    >>> import nrpy.finite_difference  # noqa: F401
    >>> par.set_parval_from_str("fd_order", 4)
    >>> header = output_constants_h("bssn", "bssn", 3, False)
    >>> "#ifndef BSSN_CONSTANTS_H" in header
    True
    >>> header.rstrip().endswith("#endif  // BSSN_CONSTANTS_H")
    True
    >>> from nrpy.helpers.generic import clang_format
    >>> "}  // END NAMESPACE: bssn::generated" in clang_format(header)
    True
    >>> "inline constexpr bool KO_ENABLED = false;" in header
    True
    >>> try:
    ...     output_constants_h("bssn", "bssn", 2, True)
    ... except ValueError as error:
    ...     print(error)
    required_padding 2 is below the stencil floor 3 at fd_order 4.
    """
    fd_order = int(par.parval_from_str("fd_order"))
    padding = int(required_padding)
    # The Kreiss-Oliger family reaches one point past the centred radius, so the
    # floor is not fd_order // 2 when it is present.  The upwinded family
    # reaches that far too, but its presence is recorded in the state header, so
    # the generated test_padding case is what compares against the exact reach.
    floor = fd_order // 2 + (1 if enable_KreissOliger_dissipation else 0)
    if padding < floor:
        raise ValueError(
            f"required_padding {padding} is below the stencil floor {floor} "
            f"at fd_order {fd_order}."
        )
    opening, closing = header_guard(f"{solver_stem}_constants.h")
    lines: List[str] = [BANNER.rstrip("\n"), opening, ""]
    lines.append(f"namespace {solver_namespace}::generated {{")
    lines.append("")
    lines.append(f"inline constexpr unsigned FD_ORDER = {fd_order};")
    lines.append(f"inline constexpr unsigned REQUIRED_PADDING = {padding};")
    ko = "true" if enable_KreissOliger_dissipation else "false"
    lines.append(f"inline constexpr bool KO_ENABLED = {ko};")
    lines.append("")
    lines.append("// clang-format off")
    lines.append(f"}}  // END NAMESPACE: {solver_namespace}::generated")
    lines.append("// clang-format on")
    lines.append("")
    lines.append(closing)
    lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
