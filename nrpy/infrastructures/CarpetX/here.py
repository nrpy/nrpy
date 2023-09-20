import os
import re
import sys
import inspect
from colored import colored

# Get the current working directory
_here = os.path.realpath(os.getcwd())


def here(*args: str) -> None:
    """
    Function to print the stack frame at the call site along with a string.
    """
    herell(False, *args)


def herecc(*args: str) -> None:
    """
    Function to print the stack frame at the call site along with a string.
    Additionally, it also includes code context.
    """
    herell(True, *args)


def herell(usecc: bool, *args: str) -> None:
    """
    Internal function to print the stack frame at the call site. If 'usecc' is True,
    the function will also include code context.
    """
    # Get the current stack frame
    stack = inspect.stack()
    frame = stack[2][0]
    herestr = "HERE:"

    # Include code context if usecc is True
    if usecc:
        herestr = re.sub(
            r"^herecc\((.*)\)$", r"HERE: \1:", frame.code_context[0].strip()
        )

    # Get the filename and line number
    fname = os.path.realpath(frame.f_code.co_filename)
    line = frame.f_lineno

    # Remove the current directory from the filename if present
    if fname.startswith(_here):
        fname = fname[len(_here) + 1 :]

    # Prepare the arguments for print
    nargs = [colored(herestr, "cyan"), f"{fname}:{colored(line, 'yellow')}"] + list(
        args
    )

    # Print the stack trace with the colored output
    print(*nargs)

    # Flush the stdout and stderr to ensure that the output is displayed
    sys.stdout.flush()
    sys.stderr.flush()


if __name__ == "__main__":
    # Test the 'here' and 'herecc' functions
    here(_here)
    herecc(_here)
