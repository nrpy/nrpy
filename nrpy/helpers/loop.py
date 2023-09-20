""" NRPy+ Loop Generation

The following script generate a single or nested loop of arbitrary
dimension in C, and has support for cache blocking (loop tiling).

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
"""

from typing import List, Union, Tuple


def loop1D(
    idx_var: str = "i",
    lower_bound: str = "0",
    upper_bound: str = "N",
    increment: str = "1",
    pragma: str = "#pragma omp parallel for",
) -> Tuple[str, str]:
    """
    Generate a one-dimensional loop in C.

    Args:
        idx_var: str, optional: The index variable for the loop.
        lower_bound: str, optional: The lower bound on the index variable.
        upper_bound: str, optional: The upper bound on the index variable.
        increment: str, optional: The increment for the index variable.
        pragma: str, optional: The OpenMP pragma (https://en.wikipedia.org/wiki/OpenMP).

    Raises:
        ValueError: If any argument has a type other than string, an error is thrown.

    Returns:
        tuple[str, str]: Returns a string header and a string footer.

    >>> header, footer = loop1D(pragma='')
    >>> print(header)
    for (int i = 0; i < N; i++) {
    <BLANKLINE>

    >>> print(footer)
    } // END LOOP: for (int i = 0; i < N; i++)
    <BLANKLINE>

    >>> header, footer = loop1D(increment='2', pragma='')
    >>> print(header)
    for (int i = 0; i < N; i += 2) {
    <BLANKLINE>
    """
    # If some argument has a different type other than string, then throw an error
    if any(
        not isinstance(i, str)
        for i in (idx_var, lower_bound, upper_bound, increment, pragma)
    ):
        raise ValueError("all parameters must have type string.")
    # Generate header and footer for a one-dimensional loop with optional parallelization using OpenMP
    pragma = pragma + "\n" if pragma else ""
    increment = " += " + increment if increment != "1" else "++"
    header = "for (int {i0} = {i1}; {i0} < {i2}; {i0}{i3})".format(
        i0=idx_var, i1=lower_bound, i2=upper_bound, i3=increment
    )
    footer = "} // END LOOP: " + header.strip() + "\n"
    return pragma + header + " {\n", footer


def loop(
    idx_var: Union[str, List[str]],
    lower_bound: Union[str, List[str]],
    upper_bound: Union[str, List[str]],
    increment: Union[str, List[str]],
    pragma: Union[str, List[str]],
    loop_body: str = "",
    tile_size: Union[str, List[str]] = "",
) -> Union[Tuple[str, str], str]:
    """
    Generate a nested loop of arbitrary dimension in C.

    Args:
        idx_var: Union[str, List[str]]: The index variable for the loop.
        lower_bound: Union[str, List[str]]: The lower bound on the index variable.
        upper_bound: Union[str, List[str]]: The upper bound on the index variable.
        increment: Union[str, List[str]]: The increment for the index variable.
        pragma: Union[str, List[str]]: The OpenMP pragma (https://en.wikipedia.org/wiki/OpenMP).
        loop_body: str, optional: The body of the loop.
        tile_size: Union[str, List[str]], optional: The tile size for cache blocking.

    Raises:
        ValueError: If all list parameters do not have the same length, an error is thrown.

    Returns:
        str: Returns a string of the loop or (header, footer).

    >>> from nrpy.helpers.generic import clang_format
    >>> header, footer = loop('i', '0', 'N', '1', '')
    >>> print(clang_format(header))
    for (int i = 0; i < N; i++) {
    <BLANKLINE>

    >>> print(clang_format(footer))
    } // END LOOP: for (int i = 0; i < N; i++)
    <BLANKLINE>

    >>> print(clang_format(loop('i', '0', 'N', '1', '', loop_body='// <INTERIOR>')))
    for (int i = 0; i < N; i++) {
      // <INTERIOR>
    } // END LOOP: for (int i = 0; i < N; i++)
    <BLANKLINE>

    >>> print(clang_format(loop('i', '0', 'N', '1', '', loop_body='// <INTERIOR>', tile_size='16')))
    for (int iB = 0; iB < N; iB += 16) {
      for (int i = iB; i < MIN(N, iB + 16); i++) {
        // <INTERIOR>
      } // END LOOP: for (int i = iB; i < MIN(N, iB + 16); i++)
    } // END LOOP: for (int iB = 0; iB < N; iB += 16)
    <BLANKLINE>
    """
    # Convert all parameters to lists for consistency
    idx_var = [idx_var] if isinstance(idx_var, str) else idx_var
    lower_bound = [lower_bound] if isinstance(lower_bound, str) else lower_bound
    upper_bound = [upper_bound] if isinstance(upper_bound, str) else upper_bound
    increment = [increment] if isinstance(increment, str) else increment
    pragma = [pragma] if isinstance(pragma, str) else pragma
    tile_size = (
        []
        if tile_size == ""
        else [tile_size]
        if isinstance(tile_size, str)
        else tile_size
    )

    if len(set(map(len, [idx_var, lower_bound, upper_bound, increment, pragma]))) != 1:
        raise ValueError("All list parameters must have the same length.")

    headers: List[str] = []
    footers: List[str] = []
    for i, var in enumerate(idx_var):
        if tile_size:
            ext_header, ext_footer = loop1D(
                f"{var}B",
                lower_bound[i],
                upper_bound[i],
                tile_size[i],
                "",
            )
            header, footer = loop1D(
                var,
                f"{var}B",
                f"MIN({upper_bound[i]}, {var}B + {tile_size[i]})",
                increment[i],
                pragma[i],
            )
            headers.insert(i, ext_header)
            footers.insert(i, ext_footer)
        else:
            header, footer = loop1D(
                var,
                lower_bound[i],
                upper_bound[i],
                increment[i],
                pragma[i],
            )
        headers.append(header)
        footers.append(footer)

    if loop_body:
        loop_body_lines = [f"{line}\n" for line in loop_body.split("\n")]
        loop_body = "".join(loop_body_lines)

    header = "".join(headers)
    footer = "".join(footers[::-1])

    return header + loop_body + footer if loop_body else (header, footer)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
