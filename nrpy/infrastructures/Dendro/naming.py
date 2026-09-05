# nrpy/infrastructures/Dendro/naming.py
"""
Reversible syntactic name transformations for Dendro generated code.

Only syntactic decorations that map one-to-one onto the exact NRPy name are
allowed here.  No semantic aliases (e.g., ``cf`` -> ``chi``) are permitted:
every machine-readable identity preserves the registered name byte-for-byte.
"""

import re

_C_IDENTIFIER_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

# NRPy RHS naming convention: <base>_rhs<suffix>, where the optional suffix
# carries the tensor-variance and component indices (e.g., "DD12", "U0").
_RHS_RE = re.compile(r"^(?P<base>.+)_rhs(?P<suffix>(?:[UD]+[0-9]+)?)$")


def validate_cpp_identifier(name: str) -> str:
    """
    Validate that a string is a legal C/C++ identifier.

    :param name: String to validate.
    :return: The validated string, unchanged.
    :raises ValueError: If the string is not a legal C/C++ identifier.
    """
    if not _C_IDENTIFIER_RE.match(name):
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

    :param name: Exact registered NRPy gridfunction name.
    :return: ``in_<name>``.
    """
    return f"in_{validate_cpp_identifier(name)}"


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

    Writers that produce *state* (initial data, projection) use this role;
    ``rhs_`` is reserved for right-hand-side output, so a reader of a
    generated signature can tell the two apart.

    DEVIATION FROM SECTION 5.2: the whitepaper enumerates four role prefixes
    (``in_``, ``rhs_``, ``diag_``, ``aux_``) and none of them denotes "state
    written by this kernel".  Reusing ``rhs_`` for the initial-data writer
    made its generated signature say it fills the right-hand-side vector,
    which is how a caller silently zeroes the state.  ``out_`` is added as a
    fifth reversible decoration rather than overloading ``rhs_``; it renames
    no scientific object, and the section 5.2 exact-name rule is preserved.

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


def aux_pointer(name: str) -> str:
    """
    Return the Dendro auxiliary-role pointer name for a gridfunction.

    :param name: Exact registered NRPy gridfunction name.
    :return: ``aux_<name>``.
    """
    return f"aux_{validate_cpp_identifier(name)}"


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
    match = _RHS_RE.fullmatch(rhs_name)
    if match is None:
        raise ValueError(f"Unrecognized NRPy RHS symbol: {rhs_name}")
    return f"{match.group('base')}{match.group('suffix')}"


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
