"""
Module to print colored text in terminal or Jupyter notebook.
"""

import sys


def not_colored(text, _):
    """
    Returns the text without any color codes.
    """
    return repr(text)


# Define the color codes
colors = {
    "red": "\033[31m",
    "green": "\033[32m",
    "yellow": "\033[33m",
    "blue": "\033[34m",
    "magenta": "\033[35m",
    "cyan": "\033[36m",
}
reset = "\033[0m"


def colored(text, color):
    """
    Returns the text string colored according to the color argument.
    """
    assert isinstance(color, str), "Color argument must be a string"
    assert color in colors, "Invalid color argument"

    return colors[color] + str(text) + reset


# Check if the standard output is a terminal
is_tty = hasattr(sys.stdout, "isatty") and sys.stdout.isatty()

# Check if running in Jupyter notebook
is_jupyter = (
    type(sys.stdout).__name__ == "OutStream"
    and type(sys.stdout).__module__ == "ipykernel.iostream"
)

# If not running in a terminal or Jupyter notebook, disable coloring
if not is_tty and not is_jupyter:
    colored = not_colored

if __name__ == "__main__":
    print(colored("This text should be green in color", "green"))
