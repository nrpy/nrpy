# nrpy/helpers/generic.py
"""
Generic helper functions used throughout NRPy.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Dave Kirby (super-fast uniq function)
"""

import inspect
import pkgutil
import platform
import subprocess
import sys
from difflib import ndiff
from pathlib import Path
from types import FrameType, ModuleType
from typing import Any, List, Optional, cast

import nrpy.params as par
from nrpy.helpers.cached_functions import is_cached, read_cached, write_cached


def superfast_uniq(seq: List[Any]) -> List[Any]:
    """
    Remove duplicate elements from a list while preserving their original order.

    :param seq: List of elements.
    :return: List with unique elements in the order they first appear in the original list.

    Example from https://www.peterbe.com/plog/uniqifiers-benchmark by Dave Kirby.
    """
    # Order preserving
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]  # type: ignore


def clang_format(c_code_str: str) -> str:
    r"""
    Format a given C code string using clang-format.

    Leverages the clang-format tool to automatically format C code according to specified style options.
    Ensures consistent code formatting across the project.

    :param c_code_str: The C code string to be formatted.
    :return: Formatted C code string.
    :raises RuntimeError: If clang-format encounters any error.
    :raises OSError: If clang-format is not installed / not found on PATH.

    Doctest:
    >>> print(clang_format(r'''int main() { printf("Hello, World!"); for(int i=0;i<10;i++) for(int j=i;j<10;j++) printf("%d %d\\n",i,j); return 0; }'''))
    int main() {
      printf("Hello, World!");
      for (int i = 0; i < 10; i++)
        for (int j = i; j < 10; j++)
          printf("%d %d\\n", i, j);
      return 0;
    }
    """
    clang_format_options = par.parval_from_str("clang_format_options")  # always a str

    # Returned cached output if inputs are exactly the same
    unique_id = __name__ + c_code_str + clang_format_options
    if is_cached(unique_id):
        return cast(str, read_cached(unique_id))

    # Build the command. Since options is a STRING, pass it as ONE argv element
    cmd = ["clang-format"]
    if clang_format_options.strip():
        cmd.append(clang_format_options)

    try:
        # Python 3.6-compatible: use Popen with PIPEs (no capture_output)
        with subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ) as process:
            stdout, stderr = process.communicate(input=c_code_str.encode("utf-8"))

            if process.returncode == 0:
                formatted = stdout.decode("utf-8", errors="replace")
                write_cached(unique_id, formatted)
                return formatted

            err = stderr.decode("utf-8", errors="replace").strip()
            raise RuntimeError(
                f"clang-format exited with code {process.returncode}.\n"
                f"Command: {' '.join(cmd)}\n\nstderr:\n{err or '<no stderr>'}"
            )

    except FileNotFoundError as exc:
        # Detailed, platform-aware installation message
        arch = platform.machine().lower()
        system = platform.system()
        details = [
            "clang-format was not found on your system (executable 'clang-format' is not in PATH).",
            "",
            "How to install:",
        ]
        if arch in ("x86_64", "amd64"):
            details.append(
                "- On x86_64, you can install via pip: `pip install clang-format` "
                "(this provides a `clang-format` executable in your Python environment)."
            )
        details.append(
            "- Otherwise, install it from your OS package manager (it is usually a standard package):\n"
            "  * macOS (Homebrew): `brew install clang-format`\n"
            "  * Ubuntu/Debian:    `sudo apt-get install clang-format`\n"
            "  * Fedora:           `sudo dnf install clang-tools-extra`\n"
            "  * Arch:             `sudo pacman -S clang`\n"
            "  * Windows (winget): `winget install LLVM.LLVM` (ensure LLVM/bin is on PATH)"
        )
        details.extend(
            [
                "",
                f"Detected platform: {system} / {arch}",
                "After installation, ensure `clang-format` is on your PATH and re-run the command.",
            ]
        )
        raise OSError("\n".join(details)) from exc


def validate_strings(
    to_check: str,
    string_desc: str,
    file_ext: str = "c",
) -> None:
    """
    Validate a string against a trusted value stored in a .c file; manage trusted file creation if file not found.

    Compare the provided string "to_check" to a trusted value stored in
    [caller module's directory]/tests/[caller module]_{string_desc}.c. Create the file with the provided string if
    it is missing. Report mismatches with detailed differences and provide instructions for updating the trusted file.

    :param to_check: Specify the string to validate, representing the expected value or output.
    :param string_desc: Provide a short, non-empty, whitespace-free label to use in the trusted file's name.
    :param file_ext: Specify the file extension for the trusted file, defaulting to "c".
    :raises ValueError: Raise if:
        - `string_desc` is empty or contains whitespace.
        - The caller's frame or filename cannot be identified.
        - The provided string does not match the trusted value.
    :raises RuntimeError: Raise if:
        - The caller's frame lacks code information, preventing file determination.
        - System or environment errors occur during file creation or access.
    """
    if not string_desc or " " in string_desc:
        raise ValueError(
            "String description cannot be blank or have whitespace inside."
        )

    # Get the caller's frame
    caller_frame: Optional[FrameType] = inspect.currentframe()
    if caller_frame is None or caller_frame.f_back is None:
        raise RuntimeError("Unable to retrieve caller frame.")
    caller_frame = caller_frame.f_back

    # Get the caller's filename
    try:
        caller_filename: str = inspect.getfile(caller_frame)
    except TypeError as exc:
        raise RuntimeError("Unable to determine the caller's filename.") from exc

    # Handle the case when called from a doctest
    if caller_filename.startswith("<doctest"):
        module: Optional[ModuleType] = inspect.getmodule(caller_frame)
        if module and hasattr(module, "__file__") and module.__file__:
            caller_filename = module.__file__
        else:
            # Fallback to sys.argv[0] if module filename is not available
            caller_filename = sys.argv[0]

    caller_directory = Path(caller_filename).parent

    outdir = caller_directory / "tests"

    # Determine the output directory and create it if it doesn't exist
    outdir.mkdir(parents=True, exist_ok=True)

    # Get the function name or script name
    if caller_frame.f_code is None:
        raise RuntimeError("Caller frame does not have code information.")
    function_name: str = caller_frame.f_code.co_name

    if function_name.startswith("<"):
        # Use the script's filename if function name is not meaningful
        function_name = Path(caller_filename).stem

    outfile_path = outdir / f"{function_name}_{string_desc}.{file_ext}"

    if outfile_path.exists():
        trusted_string = outfile_path.read_text(encoding="utf-8")
        if trusted_string != to_check:
            raise ValueError(f"""\
*** FAILURE ***
****{string_desc}**** mismatch in {outfile_path.name}!
Computed result differs from trusted value in file {outfile_path}
String output does not match the trusted version.
Differences:

{diff_strings(trusted_string, to_check)}

If you trust the new version, then delete {outfile_path} and rerun to generate a new version.
BE SURE TO INDICATE THE REASON FOR UPDATING THE TRUSTED FILE IN THE COMMIT MESSAGE.
***************
""")
    else:
        # Write the trusted string to the file
        outfile_path.write_text(to_check, encoding="utf-8")
        print(f"Trusted file '{outfile_path}' created with the provided string.")


def diff_strings(str1: str, str2: str) -> str:
    r"""
    Generate a side-by-side diff between two strings, highlighting only added or removed lines.

    Compares two multi-line strings and indicates lines that are removed or added, ignoring changes within lines
    for a cleaner overview of differences.

    :param str1: First string for comparison.
    :param str2: Second string for comparison.
    :return: Side-by-side diff of the two strings.

    Doctests:
    >>> str1 = "line1\nline2\nline3"
    >>> str2 = "line1\nline2\nline4"
    >>> print(diff_strings(str1, str2))
    - line3
    + line4
    >>> print(diff_strings(str1, str1))
    <BLANKLINE>
    """
    lines1 = str1.splitlines()
    lines2 = str2.splitlines()
    diff = ndiff(lines1, lines2)

    # Filter out lines that start with '?' or ' '
    clean_diff = [
        line for line in diff if not (line.startswith("?") or line.startswith(" "))
    ]

    return "\n".join(clean_diff)


def copy_files(
    package: str, filenames_list: List[str], project_dir: str, subdirectory: str
) -> None:
    """
    Copy specified files into a designated subdirectory within the project directory.

    Retrieves the specified files from a given Python package and copies them into a specified subdirectory
    inside the project directory. Ensures that the target subdirectory exists, creating it if necessary.

    :param package: The package name where the files are located.
    :param filenames_list: A list of filenames to be copied.
    :param project_dir: The path of the project directory where the files will be copied.
    :param subdirectory: The name of the subdirectory within the project directory.
    :raises FileNotFoundError: If a specified file is not found within the package.
    """
    # Ensure the subdirectory exists or create it if necessary
    target_subdir = Path(project_dir) / subdirectory
    target_subdir.mkdir(parents=True, exist_ok=True)

    for filename in filenames_list:
        # Retrieve the file data from the package
        data = pkgutil.get_data(package, filename)

        if data is not None:
            target_file = target_subdir / filename
            # Write the data to the target file in binary mode
            with open(target_file, "wb") as f:
                f.write(data)
        else:
            raise FileNotFoundError(
                f"Cannot find file '{filename}' in package '{package}'"
            )


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
