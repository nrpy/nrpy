"""
The colored package allows you to create output using one of six basic colors. It automatically shuts off coloring if output is not to a terminal. Exception: It will still color output in a Jupyter notebook.

Author: Steven Brandt
"""

import sys
from typing import Any, Dict
from typing_extensions import Literal

color_names = Literal["red", "green", "yellow", "blue", "magenta", "cyan"]


def coloring_is_disabled(
    arg: Any, c: color_names
) -> str:  # pylint: disable=unused-argument
    """
    Provide stringified version of the argument with coloring disabled.

    :param arg: the object to convert to a plain string
    :param c: a valid color name

    :return: a stringified version of arg
    >>> coloring_is_disabled('fish','blue')
    'fish'
    >>> coloring_is_disabled('fish','green')
    'fish'
    """
    return str(arg)


colors: Dict[color_names, str] = {
    "red": "\033[31m",
    "green": "\033[32m",
    "yellow": "\033[33m",
    "blue": "\033[34m",
    "magenta": "\033[35m",
    "cyan": "\033[36m",
}
reset = "\033[0m"


def coloring_is_enabled(arg: Any, c: color_names) -> str:
    r"""
    Return the stringified version of the argument with the specified color.
    Coloring will be disabled if the output is not being sent to a notebook cell
    or a plain old console.

    :param arg: The object to convert to a string
    :param c: The name of the color to use
    :return: A stringified and colored version of `arg`
    :raises AssertionError: If `c` is not a string or not a valid color name
    >>> import re
    >>> re.sub(r'.\[','ESC',coloring_is_enabled('fish','blue'))
    'ESC34mfishESC0m'
    >>> re.sub(r'.\[','ESC',coloring_is_enabled('fish','green'))
    'ESC32mfishESC0m'
    """
    assert isinstance(c, str), "Color name must be a string"
    assert c in colors, f"Invalid color name: {c}"
    return colors[c] + str(arg) + reset


if hasattr(sys.stdout, "isatty"):
    is_tty = sys.stdout.isatty()
else:
    is_tty = False

is_jupyter = (
    type(sys.stdout).__name__ == "OutStream"
    and type(sys.stdout).__module__ == "ipykernel.iostream"
)
if (not is_tty) and (not is_jupyter):
    colored = coloring_is_disabled
else:
    colored = coloring_is_enabled

if __name__ == "__main__":
    print(colored("red", "red"), colored("green", "green"))
