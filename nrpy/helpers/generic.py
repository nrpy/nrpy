"""
Generic helper functions used throughout NRPy.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Dave Kirby (super-fast uniq function)
"""

import base64
import hashlib
import lzma
import pkgutil
import subprocess
from difflib import ndiff
from pathlib import Path
from typing import Any, List, cast

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


# Used within c_function to create multi-line comments.
def prefix_with_star(input_string: str) -> str:
    r"""
    Prefix every line in the input string with "   * ".

    Useful for formatting multi-line comments or documentation blocks by adding a consistent prefix to each line.

    :param input_string: The input multi-line string to be prefixed.
    :return: The modified string with "* " prefixed on every line.

    Example:
    >>> s = "   Hello\n   World   "
    >>> prefix_with_star(s)
    '   * Hello\n   * World'
    """
    # Splitting the string by line breaks, trimming each line, and prefixing
    lines = input_string.split("\n")
    prefixed_lines = ["   * " + line.strip() for line in lines]

    # Joining the prefixed lines back into a single string
    result = "\n".join(prefixed_lines)
    return result


def clang_format(
    c_code_str: str,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> str:
    r"""
    Format a given C code string using clang-format.

    Leverages the clang-format tool to automatically format C code according to specified style options.
    Ensures consistent code formatting across the project.

    :param c_code_str: The C code string to be formatted.
    :param clang_format_options: Formatting options for clang-format.
    :return: Formatted C code string.
    :raises RuntimeError: If clang-format encounters any error.

    Doctest:
    >>> print(clang_format(r'''int main() { printf("Hello, World!"); for(int i=0;i<10;i++) for(int j=i;j<10;j++) printf("%d %d\n",i,j); return 0; }'''))
    int main() {
      printf("Hello, World!");
      for (int i = 0; i < 10; i++)
        for (int j = i; j < 10; j++)
          printf("%d %d\n", i, j);
      return 0;
    }
    """
    unique_id = __name__ + c_code_str + clang_format_options
    if is_cached(unique_id):
        return cast(str, read_cached(unique_id))
    # For Python 3.6 compatibility, use subprocess.PIPE instead of capture_output
    with subprocess.Popen(
        ["clang-format", clang_format_options],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ) as process:
        # Send C code string to clang-format and fetch the result
        stdout, stderr = process.communicate(input=c_code_str.encode())

        # If the process exited without errors, return the formatted code
        if process.returncode == 0:
            write_cached(unique_id, stdout.decode())
            return stdout.decode()

        raise RuntimeError(f"Error using clang-format: {stderr.decode()}")


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


def hash_to_signed_32bit(s: str) -> int:
    """
    Compute the SHA-256 hash of a string and convert it to a signed 32-bit integer.

    Generates a consistent 32-bit integer hash from an input string by:
    1. Computing the SHA-256 hash.
    2. Reducing the hash to a 32-bit unsigned integer.
    3. Converting it to a signed integer by interpreting the highest bit as the sign bit.

    :param s: Input string.
    :return: Signed 32-bit integer representation of the hash.

    >>> hash_to_signed_32bit("Hello")
    641210729
    """
    sha256 = hashlib.sha256()
    sha256.update(s.encode())
    hash_value = int(sha256.hexdigest(), 16)

    # Reduce the hash to a 32-bit unsigned integer
    reduced_hash = hash_value % (2**32)

    # Interpret the high bit as the sign bit to convert to a signed integer
    if reduced_hash >= 2**31:
        reduced_hash -= 2**32

    return reduced_hash


def compress_string_to_base64(input_string: str) -> str:
    """
    Compress the input string and return it as a Base64 encoded string.

    Uses LZMA compression to reduce the size of the input string and encodes the compressed data using Base64
    for safe transport or storage.

    :param input_string: String to be compressed.
    :return: Base64 encoded compressed string.

    Doctests:
    >>> original_string = "This is a test string that will be compressed."
    >>> compressed_string = compress_string_to_base64(original_string)
    >>> print(compressed_string)  # Output will be a Base64 encoded compressed string with line breaks
    /Td6WFoAAATm1rRGAgAhARwAAAAQz1jMAQAtVGhpcyBpcyBhIHRlc3Qgc3RyaW5nIHRoYXQgd2lsbCBiZSBjb21wcmVzc2VkLgAAAMZJbwnhOK2qAAFGLmdQc1oftvN9AQAAAAAEWVo=
    """
    # Compress the input string using LZMA with maximum compression
    compressed_data = lzma.compress(input_string.encode(), preset=9)

    # Encode the compressed data as Base64
    base64_encoded = base64.b64encode(compressed_data)

    return base64_encoded.decode()


def decompress_base64_to_string(input_base64: str) -> str:
    """
    Decompress a Base64 encoded string back to its original form.

    Decodes the Base64 string to obtain the compressed binary data and decompresses it using LZMA to retrieve the original string.

    :param input_base64: Base64 encoded compressed string.
    :return: Original decompressed string.

    Doctests:
    >>> original_string = "This is a test string that will be compressed."
    >>> compressed_string = compress_string_to_base64(original_string)
    >>> decompressed_string = decompress_base64_to_string(compressed_string)
    >>> decompressed_string == original_string
    True
    """
    # Decode the Base64 encoded string
    base64_decoded = base64.b64decode(input_base64)

    # Decompress the data using LZMA
    decompressed_data = lzma.decompress(base64_decoded)

    return decompressed_data.decode()


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
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
