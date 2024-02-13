"""
The ConditionalFileUpdater class helps users see and understand changes to generated output files.

Author: Steven R. Brandt
"""

from typing import Any, Optional, Union
import io
import os
from difflib import context_diff
import sys
from pathlib import Path
from nrpy.helpers.generic import clang_format
import nrpy.helpers.coloring as c

verbose = False
nochange = False


class ConditionalFileUpdater:
    """
    Content Comparison: It compares the existing file content with the new content to determine if an update is necessary.
    Optional Formatting: It can format the new content before the comparison and potential update, using tools like clang-format.
    Change Notification: It notifies the user of any changes or the lack thereof, with an option to display a detailed diff if in verbose mode.
    Safe Writing: It avoids writing to the file if there are no changes, reducing unnecessary disk writes and preserving file timestamps.
    >>> import os
    >>> c.colored = c.coloring_is_disabled
    >>> testf = "/tmp/_test_.txt"
    >>> os.unlink(testf)
    >>> with ConditionalFileUpdater(testf) as fd:
    ...   print("Testing", file=fd)
    Checking /tmp/_test_.txt...[written]
    >>> with ConditionalFileUpdater(testf) as fd:
    ...   print("Testing", file=fd)
    Checking /tmp/_test_.txt...[no changes]
    """

    # black insists on indenting this way, pylint does not allow it
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
        """
        Create a StringIO proxy for the file.
        We will store the output
        here so that we can later compare it with the contents of the file.

        :return: A StringIO object.
        """
        self.fd = io.StringIO()
        return self.fd

    def __exit__(self, ty: Any, val: Any, tb: Any) -> None:
        """
        Output is finished. Open the file and check whether its contents match.
        Only if they do not should the file be updated.

        :param ty: unused
        :param val: unused
        :param tb: unused

        :return: None
        """
        print("Checking", self.fname, end="...")
        newcontent = self.fd.getvalue()
        encoding = self.encoding
        if encoding is None:
            encoding = "utf-8"
        if self.do_format:
            newcontent = clang_format(newcontent)
        if os.path.exists(self.fname):
            with open(self.fname, encoding=encoding) as fd:
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
            with open(self.fname, "w", encoding=encoding) as fd:
                fd.write(newcontent)
            print(c.colored("[written]", "red"))
        else:
            print(c.colored("[no changes]", "green"))
