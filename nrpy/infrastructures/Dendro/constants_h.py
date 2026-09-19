# nrpy/infrastructures/Dendro/constants_h.py
"""
Emit the ``<stem>_constants.h`` header of generated compile-time constants.

The constants state the centered and Kreiss-Oliger finite-difference orders,
the ghost points their derivative operators reach, and whether the kernel
includes Kreiss-Oliger dissipation.  The caller passes the values recorded by
the kernel builder so the generated solver cannot disagree with its kernels.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List

from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.header_guards import header_guard

BANNER = generated_file_banner()


def output_constants_h(
    solver_stem: str,
    solver_namespace: str,
    fd_order: int,
    ko_fd_order: int,
    ko_effective_difference_order: int,
    required_padding: int,
    enable_KreissOliger_dissipation: bool,
) -> str:
    """
    Emit the generated constants header.

    It carries the numerical profile supplied by the kernel builder.  The
    builder computes the padding from the operators actually present, then
    checks it against the selected Dendro profile.  This emitter repeats that
    consistency check before writing the values.

    :param solver_stem: Lowercase formulation stem for the emitted header name.
    :param solver_namespace: NRPy-qualified solver namespace, e.g.
        ``nrpy::bssn``.
    :param fd_order: Centered finite-difference order.
    :param ko_fd_order: Base order supplied to NRPy's ``dKOD`` construction.
    :param ko_effective_difference_order: Actual even KO difference order.
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
    >>> header = output_constants_h("bssn", "bssn", 4, 2, 4, 2, False)
    >>> "#ifndef BSSN_CONSTANTS_H" in header
    True
    >>> header.rstrip().endswith("#endif  // BSSN_CONSTANTS_H")
    True
    >>> from nrpy.helpers.generic import clang_format
    >>> "}  // END NAMESPACE: bssn::generated" in clang_format(header)
    True
    >>> "inline constexpr bool KO_ENABLED = false;" in header
    True
    >>> "inline constexpr unsigned REQUIRED_PADDING = 2;" in header
    True
    """
    fd_order = int(fd_order)
    padding = int(required_padding)
    ko_order = int(ko_fd_order)
    ko_effective_order = int(ko_effective_difference_order)
    if fd_order not in (4, 6, 8):
        raise ValueError(f"Unsupported Dendro FD_ORDER={fd_order}; allowed: (4, 6, 8).")
    if (
        ko_order != fd_order - 2
        or ko_effective_order != fd_order
        or padding != fd_order // 2
    ):
        raise ValueError(
            "Dendro profile must satisfy KO_FD_ORDER=FD_ORDER-2, "
            "KO_EFFECTIVE_DIFFERENCE_ORDER=FD_ORDER, and "
            "REQUIRED_PADDING=FD_ORDER/2."
        )
    opening, closing = header_guard(f"{solver_stem}_constants.h")
    lines: List[str] = [BANNER.rstrip("\n"), opening, ""]
    lines.append(f"namespace {solver_namespace}::generated {{")
    lines.append("")
    lines.append(f"inline constexpr unsigned FD_ORDER = {fd_order};")
    lines.append(f"inline constexpr unsigned KO_FD_ORDER = {ko_order};")
    lines.append(
        f"inline constexpr unsigned KO_EFFECTIVE_DIFFERENCE_ORDER = {ko_effective_order};"
    )
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
