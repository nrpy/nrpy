# nrpy/infrastructures/Dendro/gridfunction_name_decorations.py
"""
Reversible syntactic name transformations for Dendro generated code.

Only syntactic decorations that map one-to-one onto the exact NRPy name are
allowed here.  No semantic aliases (e.g., ``cf`` -> ``chi``) are permitted:
every machine-readable identity preserves the registered name byte-for-byte.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Optional, Tuple

import nrpy.grid as gri

# NRPy RHS naming convention: <base>_rhs<suffix>, where the optional suffix
# carries the tensor-variance and component indices (e.g., "DD12", "U0").


_VARIANCE_LETTERS = frozenset({"U", "D"})
_INDEX_DIGITS = "0123456789"


def tensor_family_of(name: str) -> Optional[Tuple[str, int]]:
    """
    Split an exact NRPy name into its tensor family and rank, if it has one.

    A component name ends with a variance run (``U``/``D``) followed by one
    index digit per variance letter, so ``Z4constraintU0`` is component 0 of the
    rank-1 family ``Z4constraintU`` while ``H_Z4`` is a scalar whose name merely
    ends in a digit.  Plain string methods decide this, so no regular
    expression is needed.

    :param name: Exact registered or factory-supplied name.
    :return: ``(family, rank)`` for a tensor component, or None for a scalar.

    Doctests:
    >>> tensor_family_of("Z4constraintU0")
    ('Z4constraintU', 1)
    >>> tensor_family_of("hDD01")
    ('hDD', 2)
    >>> print(tensor_family_of("H_Z4"))
    None
    >>> print(tensor_family_of("alpha"))
    None
    """
    index_run = len(name) - len(name.rstrip(_INDEX_DIGITS))
    if index_run == 0:
        return None
    family = name[: len(name) - index_run]
    if len(family) <= index_run:
        return None
    if not set(family[-index_run:]) <= _VARIANCE_LETTERS:
        return None
    return family, index_run


def _is_variance_suffix(suffix: str) -> bool:
    """
    Return whether a suffix is empty or a variance run plus its index digits.

    :param suffix: The text following ``_rhs`` in an NRPy RHS symbol name.
    :return: True for ``""``, ``"U0"``, ``"DD12"``; False for ``"x"``, ``"U"``,
        ``"0"``, ``"DD1"``.

    Doctests:
    >>> [_is_variance_suffix(s) for s in ("", "U0", "DD12")]
    [True, True, True]
    >>> [_is_variance_suffix(s) for s in ("x", "U", "0", "DD1")]
    [False, False, False, False]
    """
    if not suffix:
        return True
    variance = suffix.rstrip(_INDEX_DIGITS)
    digits = suffix[len(variance) :]
    return (
        bool(variance)
        and len(digits) == len(variance)
        and set(variance) <= _VARIANCE_LETTERS
    )


def validate_cpp_identifier(name: str) -> str:
    """
    Validate that a string is a legal C/C++ identifier.

    :param name: String to validate.
    :return: The validated string, unchanged.
    :raises ValueError: If the string is not a legal C/C++ identifier.
    """
    # ``str.isidentifier()`` is Python's identifier grammar, which is the C
    # grammar once non-ASCII names are excluded.  No pattern is needed.
    if not (name.isidentifier() and name.isascii()):
        raise ValueError(f"Not a valid C/C++ identifier: {name!r}")
    return name


def enum_member(name: str) -> str:
    """
    Return the enum member spelling for a registered gridfunction name.

    :param name: Exact registered NRPy gridfunction name.
    :return: The validated name, unchanged.
    """
    return validate_cpp_identifier(name)


def input_pointer(name: str) -> str:
    """
    Return the Dendro input-role pointer name for a gridfunction.

    The ``in_`` spelling itself lives on
    :class:`nrpy.grid.DendroGridFunction`, which is what the shared one-point
    memory read uses, so the two cannot drift.

    :param name: Exact registered NRPy gridfunction name.
    :return: ``in_<name>``.
    """
    return gri.DendroGridFunction.input_pointer(validate_cpp_identifier(name))


def rhs_pointer(name: str) -> str:
    """
    Return the Dendro RHS-role pointer name for a gridfunction.

    :param name: Exact registered NRPy gridfunction name.
    :return: ``rhs_<name>``.
    """
    return f"rhs_{validate_cpp_identifier(name)}"


def out_pointer(name: str) -> str:
    """
    Return the Dendro output-role pointer name for a gridfunction.

    Writers that produce *state* (initial data, constraint enforcement) use this role;
    ``rhs_`` is reserved for right-hand-side output, so a reader of a
    generated signature can tell the two apart.

    The other role prefixes (``in_``, ``rhs_``, ``diag_``) all denote
    something other than "state written by this kernel".  Reusing
    ``rhs_`` for the initial-data writer made its generated signature say it
    fills the right-hand-side vector, which is how a caller silently zeroes
    the state.  ``out_`` is therefore a fifth reversible decoration rather
    than an overload of ``rhs_``; it renames no scientific object, so the
    exact-name rule is preserved.

    :param name: Exact registered NRPy gridfunction name.
    :return: ``out_<name>``.
    """
    return f"out_{validate_cpp_identifier(name)}"


def diag_pointer(name: str) -> str:
    """
    Return the Dendro diagnostic-role pointer name for a gridfunction.

    :param name: Exact registered NRPy gridfunction name.
    :return: ``diag_<name>``.
    """
    return f"diag_{validate_cpp_identifier(name)}"


def rhs_symbol_to_gridfunction_name(rhs_name: str) -> str:
    """
    Map an NRPy RHS symbol name to its evolved gridfunction name.

    The NRPy RHS naming convention is ``<base>_rhs<suffix>`` where the
    optional suffix carries the tensor-variance and component indices
    (e.g., ``DD12``, ``U0``), so ``h_rhsDD00`` maps to ``hDD00`` and
    ``alpha_rhs`` maps to ``alpha``.  The transformation is algorithmic and
    reversible; callers validate the resulting set against the active EVOL
    registry (bijection).

    :param rhs_name: NRPy RHS symbol name (e.g., ``"lambda_rhsU2"``).
    :return: The corresponding evolved gridfunction name (e.g., ``"lambdaU2"``).
    :raises ValueError: If the name does not follow the NRPy RHS convention.

    Doctests:
    >>> rhs_symbol_to_gridfunction_name("h_rhsDD00")
    'hDD00'
    >>> rhs_symbol_to_gridfunction_name("a_rhsDD12")
    'aDD12'
    >>> rhs_symbol_to_gridfunction_name("lambda_rhsU2")
    'lambdaU2'
    >>> rhs_symbol_to_gridfunction_name("alpha_rhs")
    'alpha'
    >>> rhs_symbol_to_gridfunction_name("Theta_fCCZ4_rhs")
    'Theta_fCCZ4'
    >>> try:
    ...     rhs_symbol_to_gridfunction_name("notarhs")
    ... except ValueError:
    ...     print("Unrecognized RHS symbol rejected. Good.")
    Unrecognized RHS symbol rejected. Good.
    """
    # Walk the "_rhs" occurrences from the right and take the first split whose
    # base is non-empty and whose suffix is either empty or a variance run
    # followed by one digit per variance letter -- the same shape
    # ``tensor_family_of`` above tests, and what the pattern this replaced
    # matched greedily.
    index = rhs_name.rfind("_rhs")
    while index > 0:
        base = rhs_name[:index]
        suffix = rhs_name[index + len("_rhs") :]
        if _is_variance_suffix(suffix):
            return f"{base}{suffix}"
        index = rhs_name.rfind("_rhs", 0, index)
    raise ValueError(f"Unrecognized NRPy RHS symbol: {rhs_name}")


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
