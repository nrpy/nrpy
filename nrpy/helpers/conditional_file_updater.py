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
    A class to conditionally update files based on content comparison, with optional formatting.

    This class compares the existing file content with the new content to determine if an update is necessary.
    It can optionally format the new content before the comparison and potential update using tools like clang-format.
    It notifies the user of any changes or the lack thereof, with an option to display a detailed diff if in verbose mode.
    It avoids writing to the file if there are no changes, reducing unnecessary disk writes and preserving file timestamps.

    :param fname: The name or path of the file to be potentially updated.
    :param encoding: The encoding to use when reading and writing the file. Defaults to 'utf-8'.
    :param do_format: Flag indicating whether the new content should be formatted before comparison. Defaults to False.

    >>> import os
    >>> c.is_colored = c.coloring_is_disabled
    >>> test_file = "/tmp/_test_.txt"
    >>> with open(test_file, 'w') as f: # Ensure the file exists for the unlink example
    ...     _ = f.write('Temporary file content') # Ignore the return value
    >>> os.unlink(test_file)  # Now safe to unlink, as the file is guaranteed to exist
    >>> with ConditionalFileUpdater(test_file) as fd:
    ...     print("Testing", file=fd)
    Checking /tmp/_test_.txt...[written]
    >>> with ConditionalFileUpdater(test_file) as fd:
    ...     print("Testing", file=fd)
    Checking /tmp/_test_.txt...[no changes]
    """

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

    def __exit__(
        self, exception_type: Any, exception_value: Any, traceback: Any
    ) -> None:
        """
        Exit the runtime context and compare the file's new content with its existing content, updating if necessary.

        :param exception_type: Exception type if an exception was raised within the context.
        :param exception_value: Exception value if an exception was raised within the context.
        :param traceback: Traceback object if an exception was raised within the context.
        """
        print(f"Checking {self.fname}...", end="")
        new_content = self.fd.getvalue()
        if self.do_format:
            new_content = clang_format(new_content)

        # Determine if the file exists and read its content if it does.
        old_content = ""
        if os.path.exists(self.fname):
            with open(self.fname, encoding=self.encoding) as file_descriptor:
                old_content = file_descriptor.read()

        # Decide whether to write based on a comparison of stripped content.
        should_write = new_content.strip() != old_content.strip()

        if should_write:
            if verbose:
                print("Diff for:", self.fname)
                old_lines = [line + "\n" for line in old_content.strip().split("\n")]
                new_lines = [line + "\n" for line in new_content.strip().split("\n")]
                sys.stdout.writelines(
                    context_diff(
                        old_lines, new_lines, fromfile="before", tofile="after"
                    )
                )

            if not nochange:
                with open(self.fname, "w", encoding=self.encoding) as file_descriptor:
                    file_descriptor.write(new_content)
                print(c.is_colored("[written]", "red"))
        else:
            print(c.is_colored("[no changes]", "green"))


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
