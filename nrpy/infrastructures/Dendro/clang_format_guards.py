# nrpy/infrastructures/Dendro/clang_format_guards.py
"""
Find the emitted ``END NAMESPACE`` markers clang-format would rewrite away.

``coding_style.md`` section 10 requires ``}  // END NAMESPACE: <what>`` on every
namespace closer, and clang-format's ``FixNamespaceComments`` rewrites that
comment to ``} // namespace <name>``.  Both NRPy's own formatter options and
Dendro-GR's ``.clang-format`` turn that option on, so a closer keeps its marker
only inside a ``// clang-format off`` / ``// clang-format on`` guard.  The
emitters that build a header line by line can assert on formatted output
directly; the ones carrying whole C++ files as template strings cannot be run
without a fully registered solver, so they assert on their templates through
this checker instead.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List

MARKER = "// END NAMESPACE:"
GUARD_OFF = "// clang-format off"
GUARD_ON = "// clang-format on"


def unguarded_end_namespace_markers(text: str) -> List[str]:
    r"""
    Return the ``END NAMESPACE`` marker lines that no clang-format guard covers.

    This proves the guards, not the markers: a namespace closer carrying no
    marker at all is not reported, because the marker itself is what
    ``coding_style.md`` section 10 requires and what the emitters' own doctests
    assert.

    :param text: The emitted C++ text, or a template for it.
    :return: The offending marker lines, in the order they appear.

    Doctests:
    >>> guarded = "// clang-format off\n}  // END NAMESPACE: bssn\n// clang-format on\n"
    >>> unguarded_end_namespace_markers(guarded)
    []
    >>> unguarded_end_namespace_markers("}  // END NAMESPACE: bssn\n")
    ['}  // END NAMESPACE: bssn']
    >>> unguarded_end_namespace_markers("  }  // END NAMESPACE: inner\n")
    ['  }  // END NAMESPACE: inner']
    >>> unguarded_end_namespace_markers("} // END NAMESPACE: bssn\n")
    ['} // END NAMESPACE: bssn']
    >>> two = ("// clang-format off\n}  // END NAMESPACE: a\n"
    ...        "}  // END NAMESPACE: b\n// clang-format on\n")
    >>> unguarded_end_namespace_markers(two)
    []
    """
    offenders: List[str] = []
    guarded = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped == GUARD_OFF:
            guarded = True
            continue
        if stripped == GUARD_ON:
            guarded = False
            continue
        if stripped.startswith("}") and MARKER in stripped and not guarded:
            offenders.append(line)
    return offenders


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
