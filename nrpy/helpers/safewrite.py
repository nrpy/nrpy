"""
The Safewrite class helps users see and understand changes to generated output files.

Author: Steven R. Brandt
"""
from typing import Any, Optional, Union
import io
import os
from difflib import context_diff
import sys
from subprocess import Popen, PIPE
from pathlib import Path
from nrpy.helpers.colored import colored

try:
    from clang_format import _get_executable as get_executable  # type: ignore[import-not-found]
except ModuleNotFoundError:

    def get_executable(_: str) -> str:
        """Use /bin/cat if clang_format is not available."""
        return "/bin/cat"


clang_formatter = get_executable("clang-format")

verbose = False
nochange = False


class SafeWrite:
    """
    SafeWrite ensures that files are only written if their contents would change.
    In addition, it prints a message to the screen notifying the user if the file has changed.
    If verbose is set, it will print a diff of the change.
    If nochange is set, then it will avoid changing file contents at all. This setting
    is useful when debugging.
    """

    # black insists on indenting this way, pylint does not allow it
    # pylint: disable=C0330
    def __init__(
        self,
        fname: Union[Path, str],
        encoding: Optional[str] = None,
        do_format: bool = False,
    ) -> None:
        self.fname = os.path.abspath(str(fname))
        self.fd: io.StringIO
        self.do_format = do_format
        self.encoding = encoding

    def __enter__(self) -> io.TextIOWrapper:
        self.fd = io.StringIO()
        return self.fd

    def __exit__(self, ty: Any, val: Any, tb: Any) -> None:
        print("Checking", self.fname, end="...")
        newcontent = self.fd.getvalue()
        if self.do_format:
            pipe = Popen(
                [clang_formatter], stdout=PIPE, stdin=PIPE, universal_newlines=True
            )
            out, _err = pipe.communicate(newcontent)
            newcontent = out
        if os.path.exists(self.fname):
            if self.encoding is not None:
                with open(self.fname, encoding=self.encoding) as fd:
                    oldcontent = fd.read()
            else:
                with open(self.fname) as fd:
                    oldcontent = fd.read()
            do_write = newcontent.strip() != oldcontent.strip()
            if do_write and verbose:
                print("Diff for:", self.fname)
                oldlines = [line + "\n" for line in oldcontent.strip().split("\n")]
                newlines = [line + "\n" for line in newcontent.strip().split("\n")]
                sys.stdout.writelines(
                    context_diff(oldlines, newlines, fromfile="before", tofile="after")
                )
        else:
            do_write = True
        if do_write:
            assert not nochange
            if self.encoding is not None:
                with open(self.fname, "w", encoding=self.encoding) as fd:
                    fd.write(newcontent)
            else:
                with open(self.fname, "w") as fd:
                    fd.write(newcontent)
            print(colored("[written]", "red"))
        else:
            print(colored("[no changes]", "green"))
