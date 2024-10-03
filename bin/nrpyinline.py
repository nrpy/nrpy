"""
Extract and execute Python code from a text file marked by 'NRPYSTART' and 'NRPYEND'.

This version ignores leading //, *, /*, and trailing */ along with any whitespace before or after these comment markers.

Example usage is provided in the `run_script_from_file` function's docstring.

Author: Zachariah Etienne
Email: zachetie *at** gmail **dot * com
"""

import argparse
import doctest
import sys
import textwrap
import traceback
from typing import Any, Dict, List


def strip_comment_markers(line: str) -> str:
    """
    Remove leading //, *, /*, and trailing */ from a line, preserving code indentation.

    :param line: The line of text to process.
    :return: The line with comment markers removed.
    """
    # Remove leading whitespace to find comment markers
    stripped_line = line.lstrip()
    leading_spaces = line[: len(line) - len(stripped_line)]

    # Remove leading comment markers
    if stripped_line.startswith("//"):
        stripped_line = stripped_line[2:].lstrip()
    elif stripped_line.startswith("/*"):
        stripped_line = stripped_line[2:].lstrip()
    elif stripped_line.startswith("*"):
        stripped_line = stripped_line[1:].lstrip()

    # Remove trailing comment markers
    if stripped_line.rstrip().endswith("*/"):
        stripped_line = stripped_line.rstrip()[:-2].rstrip()

    # Return the line with original leading spaces and processed content
    return leading_spaces + stripped_line


def run_script_from_file(file_path: str) -> None:
    """
    Extract and execute Python code between 'NRPYSTART' and 'NRPYEND' markers.

    Read a text file line by line, extract Python code between 'NRPYSTART' and 'NRPYEND' markers,
    and execute each script found between these markers.

    The script ignores leading //, *, /*, and trailing */ along with any whitespace before or after these comment markers.

    :param file_path: The path to the text file containing the Python scripts.
    :raises ValueError: If 'NRPYEND' is encountered without a preceding 'NRPYSTART' or
                        if 'NRPYSTART' is encountered while already capturing a script block.
    :raises Exception: Any exception raised during the execution of the extracted script.

    Example usage:
    >>> import os
    >>> from appdirs import user_cache_dir
    >>> cache_dir = user_cache_dir('run_script_from_file')
    >>> os.makedirs(cache_dir, exist_ok=True)
    >>> file_path = os.path.join(cache_dir, 'dummy_path.txt')

    # Example 1: File with code blocks inside comments
    >>> file_content = '''Some random text...
    ...    // NRPYSTART
    ...    //     print("Hello from script 1")
    ...    // NRPYEND
    ... Some more text...
    ...    /* NRPYSTART */
    ...    /*     for i in range(2): */
    ...    /*         print(f"Loop iteration {i}") */
    ...    /* NRPYEND */
    ... Even more text...
    ... Some more text...
    ...    /* NRPYSTART
    ...    *     for i in range(2): */
    ...    *         print(f"Loop2 iteration {i}")
    ...    NRPYEND */
    ... Even more text...
    ... /*
    ...   NRPYSTART
    ... print("hello3")
    ... NRPYEND
    ... */
    ... this text should also be ignored /**/ //'''
    >>> with open(file_path, 'w', encoding='utf-8') as f:
    ...     _ = f.write(file_content)
    >>> from io import StringIO
    >>> from contextlib import redirect_stdout
    >>> f_output = StringIO()
    >>> with redirect_stdout(f_output):
    ...     run_script_from_file(file_path)
    >>> print(f_output.getvalue(), end='')
    Hello from script 1
    Loop iteration 0
    Loop iteration 1
    Loop2 iteration 0
    Loop2 iteration 1
    hello3
    >>> os.remove(file_path)

    # Example 2: 'NRPYEND' appears without 'NRPYSTART'
    >>> file_content = '''Some text...
    ... NRPYEND'''
    >>> with open(file_path, 'w', encoding='utf-8') as f:
    ...     _ = f.write(file_content)
    >>> try:
    ...     run_script_from_file(file_path)
    ... except ValueError as e:
    ...     print(f"ValueError: {e}")
    ValueError: Found 'NRPYEND' without encountering 'NRPYSTART'.
    >>> os.remove(file_path)

    # Example 3: 'NRPYSTART' appears twice without 'NRPYEND'
    >>> file_content = '''NRPYSTART
    ... NRPYSTART'''
    >>> with open(file_path, 'w', encoding='utf-8') as f:
    ...     _ = f.write(file_content)
    >>> try:
    ...     run_script_from_file(file_path)
    ... except ValueError as e:
    ...     print(f"ValueError: {e}")
    ValueError: Found 'NRPYSTART' while already capturing a script block.
    >>> os.remove(file_path)

    # Example 4: Missing 'NRPYEND' at the end of the file
    >>> file_content = '''NRPYSTART
    ... print("This script never ends")'''
    >>> with open(file_path, 'w', encoding='utf-8') as f:
    ...     _ = f.write(file_content)
    >>> try:
    ...     run_script_from_file(file_path)
    ... except ValueError as e:
    ...     print(f"ValueError: {e}")
    ValueError: File ended while still capturing a script block. Missing 'NRPYEND'.
    >>> os.remove(file_path)
    """
    script_lines: List[str] = []
    capturing: bool = False
    namespace: Dict[str, Any] = {
        "__builtins__": __builtins__
    }  # Include built-ins for exec()

    with open(file_path, "r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, 1):
            # Process the line to strip comment markers
            processed_line = strip_comment_markers(line)

            # Check for the start marker
            if processed_line.strip() == "NRPYSTART":
                if capturing:
                    raise ValueError(
                        "Found 'NRPYSTART' while already capturing a script block."
                    )
                capturing = True
                script_lines = []
                continue  # Move to the next line

            # Check for the end marker
            if processed_line.strip() == "NRPYEND":
                if not capturing:
                    raise ValueError(
                        "Found 'NRPYEND' without encountering 'NRPYSTART'."
                    )
                # Execute the captured script
                script_code = "".join(script_lines)
                script_code = textwrap.dedent(script_code)  # Remove common indentation
                try:
                    # Warning: this will exec the Python script.
                    exec(script_code, namespace)  # pylint: disable=exec-used
                except Exception:  # Removed 'as e' since it's unused
                    # Print the entire traceback along with the error message
                    print(
                        f"Error executing script ending at line {line_number}:\n{traceback.format_exc()}",
                        file=sys.stderr,
                    )
                    raise
                capturing = False
                continue  # Move to the next line

            # Capture lines if within a script block
            if capturing:
                # Strip comment markers from the line before adding to script_lines
                processed_line = strip_comment_markers(line)
                script_lines.append(processed_line)

    # After file processing, check if a script was left unclosed
    if capturing:
        raise ValueError(
            "File ended while still capturing a script block. Missing 'NRPYEND'."
        )


def main() -> None:
    """Parse command-line arguments and invoke the script runner."""
    parser = argparse.ArgumentParser(
        description="Run Python scripts from a file with NRPYSTART/NRPYEND markers."
    )
    parser.add_argument(
        "file_path", help="Path to the text file containing Python scripts."
    )
    args = parser.parse_args()

    try:
        run_script_from_file(args.file_path)
    except ValueError as ve:
        print(f"Error: {ve}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError as fnfe:
        print(f"Error: {fnfe}", file=sys.stderr)
        sys.exit(1)
    except SyntaxError as se:
        print(f"Syntax error in the script: {se}", file=sys.stderr)
        sys.exit(1)
    # Let other exceptions propagate without catching broad exceptions


if __name__ == "__main__":
    results = doctest.testmod(
        optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS
    )

    if results.failed > 0:
        print(
            f"Doctest failed: {results.failed} of {results.attempted} test(s)",
            file=sys.stderr,
        )
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
