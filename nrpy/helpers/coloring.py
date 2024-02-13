"""
The coloring package allows you to create output using one of six basic colors.
It automatically shuts off coloring if output is not to a terminal.
Exception: It will still color output in a Jupyter notebook.

Author: Steven Brandt
"""

import sys
from typing import Any, Dict
from typing_extensions import Literal

# Define the type for color names
ColorNames = Literal["red", "green", "yellow", "blue", "magenta", "cyan"]

# Dictionary mapping color names to terminal color codes
colors: Dict[ColorNames, str] = {
    "red": "\033[31m",
    "green": "\033[32m",
    "yellow": "\033[33m",
    "blue": "\033[34m",
    "magenta": "\033[35m",
    "cyan": "\033[36m",
}

# Reset color code
reset: str = "\033[0m"


def coloring_is_disabled(
    arg: Any, c: ColorNames  # pylint: disable=unused-argument
) -> str:
    """
    Provide a stringified version of the argument with coloring disabled.

    :param arg: The object to convert to a plain string.
    :param c: A valid color name, unused here but kept for interface consistency.
    :return: A stringified and uncolored version of `arg`.

    >>> coloring_is_disabled('fish', 'blue')
    'fish'
    >>> coloring_is_disabled('fish', 'green')
    'fish'
    """
    return str(arg)


def coloring_is_enabled(arg: Any, c: ColorNames) -> str:
    r"""
    Return the stringified version of the argument with the specified color.
    Coloring will be disabled if the output is not being sent to a notebook cell
    or a plain old console.

    :param arg: The object to convert to a string.
    :param c: The name of the color to use.
    :return: A stringified and colored version of `arg`.
    :raises AssertionError: If `c` is not a string or not a valid color name.

    >>> import re
    >>> re.sub(r'\033\[', 'ESC', coloring_is_enabled('fish', 'blue'))
    'ESC34mfishESC0m'
    >>> re.sub(r'\033\[', 'ESC', coloring_is_enabled('fish', 'green'))
    'ESC32mfishESC0m'
    """
    assert isinstance(c, str), "Color name must be a string"
    assert c in colors, f"Invalid color name: {c}"
    return colors[c] + str(arg) + reset


# Determine if output is to a terminal or a Jupyter notebook
is_tty: bool = sys.stdout.isatty() if hasattr(sys.stdout, "isatty") else False
is_jupyter: bool = (
    type(sys.stdout).__name__ == "OutStream"
    and type(sys.stdout).__module__ == "ipykernel.iostream"
)

# Choose the appropriate coloring function based on the output destination
is_colored = (
    coloring_is_disabled if not is_tty and not is_jupyter else coloring_is_enabled
)

if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
