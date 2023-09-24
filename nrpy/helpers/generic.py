"""
Generic helper functions used throughout NRPy+.
"""
from typing import List, Any, cast
import subprocess
import hashlib
import lzma
import base64
from difflib import ndiff
from nrpy.helpers.cached_functions import is_cached, read_cached, write_cached


#
def superfast_uniq(seq: List[Any]) -> List[Any]:
    """
    super fast 'uniq' function, that preserves order.
    f8() function from https://www.peterbe.com/plog/uniqifiers-benchmark
    Author: Dave Kirby
    """
    # Order preserving
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]  # type: ignore


# Used within c_function to create multi-line comments.
def prefix_with_star(input_string: str) -> str:
    """
    Prefixes every line in the given multi-line string with "* ".
    First, removes leading and trailing whitespaces from each line.

    Args:
    - input_string (str): The input multi-line string to be prefixed.

    Returns:
    - str: The modified string with "* " prefixed on every line.

    Example:
    >>> s = "   Hello   \\n   World   "
    >>> prefix_with_star(s)
    '   * Hello\\n   * World'
    """

    # Splitting the string by line breaks, trimming each line, and prefixing
    lines = input_string.split("\n")
    prefixed_lines = ["   * " + line.strip() for line in lines]

    # Joining the prefixed lines back into a single string
    result = "\n".join(prefixed_lines)
    return result


def clang_format(
    c_code_str: str,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 0}",
) -> str:
    """
    Formats a given C code string using clang-format.

    Parameters:
    - c_code_str (str): The C code string that needs to be formatted.

    Returns:
    - str: Formatted C code string.

    Raises:
    - RuntimeError: If clang-format encounters any error.

    Examples:
    >>> print(clang_format(r'''int main() { printf("Hello, World!"); for(int i=0;i<10;i++) for(int j=i;j<10;j++) printf("%d %d\\n",i,j); return 0; }'''))
    int main() {
      printf("Hello, World!");
      for (int i = 0; i < 10; i++)
        for (int j = i; j < 10; j++)
          printf("%d %d\\n", i, j);
      return 0;
    }
    """
    unique_id = __name__ + c_code_str + clang_format_options
    if is_cached(unique_id):
        return cast(str, read_cached(unique_id))
    with subprocess.Popen(
        ["clang-format", clang_format_options],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ) as process:
        # Send your C code string to clang-format and fetch the result
        stdout, stderr = process.communicate(input=c_code_str.encode())

        # If the process exited without errors, return the formatted code
        if process.returncode == 0:
            write_cached(unique_id, stdout.decode())
            return stdout.decode()

        raise RuntimeError(f"Error using clang-format: {stderr.decode()}")


def diff_strings(str1: str, str2: str) -> str:
    """
    Generates a side-by-side diff between two strings using difflib, excluding intra-line details.

    :param str1: The first string to compare.
    :param str2: The second string to compare.
    :return: The side-by-side diff between the two strings.

    Doctests:
    >>> str1 = "line1\\nline2\\nline3"
    >>> str2 = "line1\\nline2\\nline4"
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
    Compute the SHA-256 hash of a given string and convert it to a signed 32-bit integer.

    :param s: The input string.
    :return: The signed 32-bit integer representation.

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
    Compresses the given input string using LZMA with maximum compression and returns the result as a Base64 encoded string.

    :param input_string: The string to be compressed.
    :return: The Base64 encoded compressed string.

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
    Decompresses the given Base64 encoded input string that was compressed using LZMA with maximum compression and returns the original string.

    :param input_base64: The Base64 encoded compressed string.
    :return: The original decompressed string.

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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
