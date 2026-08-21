# nrpy/helpers/cached_functions.py
"""
Provide caching functionality to accelerate the execution of repeated, time-consuming tasks.

Offer utility functions that allow for:
- Generating and checking cache files based on unique IDs.
- Storing and retrieving cached data.
- Efficiently simplifying SymPy expressions by leveraging cached results.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import hashlib
import pickle
import tempfile
from pathlib import Path
from typing import Any, cast

import sympy as sp
from appdirs import user_cache_dir  # type: ignore


def get_hash(unique_id: str) -> str:
    """
    Generate a SHA-256 hash string for a given unique ID.

    :param unique_id: A unique identifier to be hashed.
    :return: The SHA-256 hash string.

    Doctests:
    >>> get_hash("test_id") == hashlib.sha256("test_id".encode("utf-8")).hexdigest()
    True
    """
    return hashlib.sha256(unique_id.encode("utf-8")).hexdigest()


def cache_file(unique_id: str) -> Path:
    """
    Generate a cache file path using a unique ID.

    :param unique_id: A unique identifier for generating the file path.
    :return: The cache file path.

    Doctests:
    >>> cache_file("test_id").name == f"{get_hash('test_id')}.nrpycache"
    True
    """
    if not Path(user_cache_dir("nrpy")).exists():
        Path(user_cache_dir("nrpy")).mkdir(parents=True, exist_ok=True)
    return Path(user_cache_dir("nrpy")) / f"{get_hash(unique_id)}.nrpycache"


def is_cached(unique_id: str) -> bool:
    """
    Check if the file with the given unique ID exists in cache.

    :param unique_id: A unique identifier to check for.
    :return: True if the file exists, False otherwise.

    Doctests:
    >>> # Assuming "nonexistent_id" is not cached
    >>> is_cached("nonexistent_id")
    False
    """
    return cache_file(unique_id).exists()


def read_cached(unique_id: str) -> Any:
    """
    Read the cached data associated with a unique ID.

    :param unique_id: A unique identifier to read data for.
    :return: The data read from the cache file.

    Doctests:
    >>> write_cached('test_read', {'data': 123})
    >>> read_cached('test_read') == {'data': 123}
    True
    """
    # print(f"Reading " + str(unique_id[:80]).replace("\n", ""))
    with open(cache_file(unique_id), "rb") as file:
        # print(f"Reading cached file {file.name}.")
        return pickle.load(file)


def _write_pickle_via_temporary_file(cache_path: Path, data: Any) -> None:
    """
    Write pickled data via a unique temporary file in the cache directory.

    :param cache_path: Destination path for the pickled data.
    :param data: Data to pickle.
    """
    temporary_file = tempfile.NamedTemporaryFile(
        mode="wb", dir=str(cache_path.parent), delete=False
    )
    temporary_path = Path(temporary_file.name)
    try:
        with temporary_file as file:
            pickle.dump(data, file)
        temporary_path.replace(cache_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def write_cached(unique_id: str, data: Any) -> None:
    """
    Write data to the cache file associated with a unique ID.

    :param unique_id: A unique identifier to write data for.
    :param data: The data to be written to the cache.

    Doctests:
    >>> write_cached('test_write', {'value': 456})
    >>> read_cached('test_write') == {'value': 456}
    True
    >>> import os
    >>> from threading import Event, Thread
    >>> from unittest.mock import patch
    >>> concurrent_write_test_id = f"test_concurrent_write_{os.getpid()}"
    >>> write_cached(concurrent_write_test_id, "old")
    >>> dump_started = Event()
    >>> allow_dump = Event()
    >>> original_dump = pickle.dump
    >>> def paused_dump(data, file):
    ...     dump_started.set()
    ...     if not allow_dump.wait(timeout=5):
    ...         raise RuntimeError("Timed out waiting to finish cache dump.")
    ...     original_dump(data, file)
    >>> with patch.object(pickle, "dump", paused_dump):
    ...     writer = Thread(target=write_cached, args=(concurrent_write_test_id, "new"))
    ...     writer.start()
    ...     started = dump_started.wait(timeout=5)
    ...     try:
    ...         value_during_write = read_cached(concurrent_write_test_id)
    ...     finally:
    ...         allow_dump.set()
    ...         writer.join(timeout=5)
    >>> started
    True
    >>> writer.is_alive()
    False
    >>> value_during_write
    'old'
    >>> read_cached(concurrent_write_test_id)
    'new'
    """
    # print(f"Writing " + str(unique_id[:80]).replace("\n", ""))
    _write_pickle_via_temporary_file(cache_file(unique_id), data)


def cached_simplify(expr: sp.Basic) -> sp.Expr:
    r"""
    Simplify a given sympy expression using a cache to speed up repeated simplifications.

    :param expr: SymPy expression to be simplified.
    :return: Simplified SymPy expression.

    Doctests:
    >>> x = sp.symbols('x')
    >>> expr = sp.sympify("x**2 + 2*x + 1")
    >>> simplified_expr = cached_simplify(expr)
    >>> simplified_expr.equals((x + 1)**2)
    True
    >>> import sys
    >>> from threading import Event, Thread
    >>> from unittest.mock import patch
    >>> concurrent_expr = sp.sympify("x + x")
    >>> expected_expr = sp.sympify("2*x")
    >>> module = sys.modules[cached_simplify.__module__]
    >>> dump_started = Event()
    >>> allow_dump = Event()
    >>> original_dump = pickle.dump
    >>> def paused_first_dump(data, file):
    ...     if not dump_started.is_set():
    ...         dump_started.set()
    ...         if not allow_dump.wait(timeout=5):
    ...             raise RuntimeError("Timed out waiting to finish cache dump.")
    ...     original_dump(data, file)
    >>> worker_results = []
    >>> def cached_write():
    ...     worker_results.append(cached_simplify(concurrent_expr))
    >>> with tempfile.TemporaryDirectory(prefix="nrpy-cache-") as cache_root:
    ...     with patch.object(module, "user_cache_dir", return_value=cache_root):
    ...         with patch.object(pickle, "dump", paused_first_dump):
    ...             writer = Thread(target=cached_write)
    ...             writer.start()
    ...             started = dump_started.wait(timeout=5)
    ...             try:
    ...                 concurrent_result = cached_simplify(concurrent_expr)
    ...             finally:
    ...                 allow_dump.set()
    ...                 writer.join(timeout=5)
    ...             if writer.is_alive():
    ...                 raise RuntimeError("Timed out waiting for cache writer.")
    >>> started and concurrent_result == expected_expr
    True
    >>> worker_results == [expected_expr]
    True
    """
    if expr == sp.sympify(0):
        return cast(sp.Expr, sp.sympify(0))
    cache_dir = Path(user_cache_dir("nrpy"))
    try:
        pickle_expr = pickle.dumps(expr)
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_fle = cache_dir / (hashlib.sha256(pickle_expr).hexdigest() + ".nrpycache")
        if cache_fle.exists():
            with open(cache_fle, "rb") as file:
                return cast(sp.Expr, pickle.load(file))
        else:
            simplified = sp.simplify(expr)
            _write_pickle_via_temporary_file(cache_fle, simplified)
            return cast(sp.Expr, simplified)
    except pickle.PicklingError:
        return cast(sp.Expr, sp.simplify(expr))


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
