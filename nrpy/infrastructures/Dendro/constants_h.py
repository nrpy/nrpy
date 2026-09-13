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
    takes it from the widest reach of the derivative operators the emitted
    kernel actually contains.  It is not ``fd_order // 2``: the upwinded and
    Kreiss-Oliger families reach one point further than the centered ones, and a
    host that sized its ghost zones from the radius would have the kernel read
    past the end of a block.  That is exactly why this emitter states the
    recorded reach and derives nothing: a floor rebuilt here from the order and
    the dissipation switch would be a second stencil model, and it would be
    wrong in the very configuration both applications ship -- at ``fd_order``
    4 with Kreiss-Oliger off it gives 2, while the emitted kernel's upwinded
    operators reach 3.  Only the coefficients know which families a kernel
    actually contains.  The ``dfullupD``/``dfulldnD`` families reach
    ``fd_order``.  They have stencils and no C-code path: ``c_codegen`` classifies
    them as ordinary symbols rather than rejecting them, so a kernel containing
    one would emit an undeclared identifier and fail to compile rather than
    lower incorrectly.

    :param solver_stem: Lowercase formulation stem for the emitted header name.
    :param solver_namespace: Solver namespace, following Dendro's lowercase
        formulation habit (``namespace bssn``).
    :param required_padding: Ghost points required on every axis.  Pass the
        exact reach the kernel builder recorded through
        ``CFunction_roles.set_required_padding``, which every production
        caller reads back with ``CFunction_roles.required_padding()``.
        Nothing here re-derives or bounds it: the coefficient-derived reach in
        ``nrpy.finite_difference.stencil_reach_per_axis`` is its only
        authority, so a caller that supplies some other number is emitting a
        padding this generator cannot vouch for.
    :param enable_KreissOliger_dissipation: Whether the emitted kernel carries
        Kreiss-Oliger dissipation, passed by the caller that built the kernels.
        It is recorded as ``KO_ENABLED`` so a reader of the generated header
        can tell which dissipation the kernel was lowered with.  Real-host
        profile validation compares ``ko_enabled`` with this constant, but no
        emitted numerical code reads it or derives a stencil reach from it.
        The emitted parameter file records the same switch as ``ko_enabled``.
    :return: The complete C++ header text.

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
    >>> "inline constexpr unsigned REQUIRED_PADDING = 3;" in header
    True
    """
    fd_order = int(par.parval_from_str("fd_order"))
    padding = int(required_padding)
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
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
