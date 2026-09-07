# nrpy/infrastructures/Dendro/header_guards.py
"""
Build the traditional header guard every generated Dendro header carries.

``coding_style.md`` section 4 requires the ``#ifndef``/``#define`` form with an
UPPER_SNAKE_CASE macro, and BHaH's own emitted headers use it, so the generated
Dendro headers do too.  The macro is derived from the header's file name, which
already carries the solver stem, so the headers this builds guards for cannot
collide when two generated solvers sit in one Dendro-GR checkout.  The copied
``dendro_standalone_host.h`` is not one of them: it carries no stem, it is
byte-identical in every project, and each target's include path holds only its
own copy.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Tuple


def header_guard(header_name: str) -> Tuple[str, str]:
    """
    Return the opening and closing header-guard text for one header.

    :param header_name: File name of the generated header, e.g. ``bssn_types.h``.
    :return: The two lines that open the guard, and the line that closes it.

    Doctests:
    >>> opening, closing = header_guard("bssn_types.h")
    >>> print(opening)
    #ifndef BSSN_TYPES_H
    #define BSSN_TYPES_H
    >>> print(closing)
    #endif  // BSSN_TYPES_H
    >>> header_guard("bssnCtx.h")[0].splitlines()[0]
    '#ifndef BSSNCTX_H'
    """
    macro = "".join(
        character if character.isalnum() else "_" for character in header_name
    ).upper()
    return f"#ifndef {macro}\n#define {macro}", f"#endif  // {macro}"


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
