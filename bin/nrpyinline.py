"""
Extract and execute Python code from a text file marked by 'NRPYSTART' and 'NRPYEND'.

Example usage is provided in the `run_script_from_file` function's docstring.

Author: Zachariah Etienne
Email: zachetie *at** gmail **dot * com
"""

import argparse
import doctest
import sys
from typing import Any, Dict, List


def run_script_from_file(file_path: str) -> None:
    """
    Extract and execute Python code between 'NRPYSTART' and 'NRPYEND' markers.

    Read a text file line by line, extract Python code between 'NRPYSTART' and 'NRPYEND' markers,
    and execute each script found between these markers.

    :param file_path: The path to the text file containing the Python scripts.
    :raises ValueError: If 'NRPYEND' is encountered without a preceding 'NRPYSTART' or
                        if 'NRPYSTART' is encountered while already capturing a script block.
    :raises Exception: Any exception raised during the execution of the extracted script.

    Example usage:
    >>> from unittest.mock import mock_open, patch

    # Example 1: File with two NRPYSTART/NRPYEND code blocks
    >>> file_content = '''Some random text...
    ...    NRPYSTART
    ... print("Hello from script 1")
    ... NRPYEND
    ... Some more text...
    ... NRPYSTART
    ... for i in range(2):
    ...     print(f"Loop iteration {i}")
    ...    NRPYEND
    ... Even more text...'''
    >>> m = mock_open(read_data=file_content)
    >>> with patch('builtins.open', m):
    ...     run_script_from_file('dummy_path.txt')
    Hello from script 1
    Loop iteration 0
    Loop iteration 1

    # Example 2: 'NRPYEND' appears without 'NRPYSTART'
    >>> file_content = '''Some text...
    ... NRPYEND'''
    >>> m = mock_open(read_data=file_content)
    >>> with patch('builtins.open', m):
    ...     run_script_from_file('dummy_path.txt')
    Traceback (most recent call last):
    ...
    ValueError: Found 'NRPYEND' without encountering 'NRPYSTART'.

    # Example 3: 'NRPYSTART' appears twice without 'NRPYEND'
    >>> file_content = '''NRPYSTART
    ... NRPYSTART'''
    >>> m = mock_open(read_data=file_content)
    >>> with patch('builtins.open', m):
    ...     run_script_from_file('dummy_path.txt')
    Traceback (most recent call last):
    ...
    ValueError: Found 'NRPYSTART' while already capturing a script block.

    # Example 4: Missing 'NRPYEND' at the end of the file
    >>> file_content = '''NRPYSTART
    ... print("This script never ends")'''
    >>> m = mock_open(read_data=file_content)
    >>> with patch('builtins.open', m):
    ...     run_script_from_file('dummy_path.txt')
    Traceback (most recent call last):
    ...
    ValueError: File ended while still capturing a script block. Missing 'NRPYEND'.
    """
    script_lines: List[str] = []
    capturing: bool = False
    namespace: Dict[str, Any] = {}  # Separate namespace for executed scripts

    with open(file_path, "r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, 1):
            # Strip leading and trailing whitespace
            stripped_line = line.strip()

            # Check for the start marker
            if stripped_line == "NRPYSTART":
                if capturing:
                    raise ValueError(
                        "Found 'NRPYSTART' while already capturing a script block."
                    )
                capturing = True
                script_lines = []
                continue  # Move to the next line

            # Check for the end marker
            if stripped_line == "NRPYEND":
                if not capturing:
                    raise ValueError(
                        "Found 'NRPYEND' without encountering 'NRPYSTART'."
                    )
                # Execute the captured script
                script_code = "".join(script_lines)
                try:
                    # Warning: this will exec the Python script.
                    exec(script_code, namespace)  # pylint: disable=exec-used
                except Exception as e:  # pylint: disable=broad-except
                    # Catching Exception is acceptable here because we want to handle any error from exec
                    print(
                        f"Error executing script ending at line {line_number}: {e}",
                        file=sys.stderr,
                    )
                    raise
                capturing = False
                continue  # Move to the next line

            # Capture lines if within a script block
            if capturing:
                script_lines.append(line)

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
    results = doctest.testmod()

    if results.failed > 0:
        print(
            f"Doctest failed: {results.failed} of {results.attempted} test(s)",
            file=sys.stderr,
        )
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
