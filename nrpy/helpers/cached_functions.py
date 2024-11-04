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
from pathlib import Path
from typing import Any, cast

import sympy as sp
from appdirs import user_cache_dir  # type: ignore

import nrpy.params as par


def get_hash(unique_id: str) -> str:
    """
    Generate a SHA-256 hash string for a given unique ID.

    :param unique_id: A unique identifier to be hashed.
    :return: The SHA-256 hash string.

    DocTests:
    >>> get_hash("test_id") == hashlib.sha256("test_id".encode("utf-8")).hexdigest()
    True
    """
    return hashlib.sha256(unique_id.encode("utf-8")).hexdigest()


def cache_file(unique_id: str) -> Path:
    """
    Generate a cache file path using a unique ID.

    :param unique_id: A unique identifier for generating the file path.
    :return: The cache file path.

    DocTests:
    >>> cache_file("test_id").name == f"{get_hash('test_id')}.nrpycache"
    True
    """
    if not Path(user_cache_dir("nrpy")).exists():
        Path(user_cache_dir("nrpy")).mkdir(parents=True, exist_ok=True)
    return Path(user_cache_dir("nrpy")) / f"{get_hash(unique_id)}.nrpycache"


def NRPy_params_checksum() -> str:
    """
    Generate a checksum of NRPy+ parameters stored in par.glb_params_dict.

    :return: The checksum string.

    DocTests:
    >>> isinstance(NRPy_params_checksum(), str)
    True
    """
    return get_hash(str(pickle.dumps(dict(sorted(par.glb_params_dict.items())))))


def is_cached(unique_id: str) -> bool:
    """
    Check if the file with the given unique ID exists in cache.

    :param unique_id: A unique identifier to check for.
    :return: True if the file exists, False otherwise.

    DocTests:
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

    DocTests:
    >>> write_cached('test_read', {'data': 123})
    >>> read_cached('test_read') == {'data': 123}
    True
    """
    # print(f"Reading " + str(unique_id[:80]).replace("\n", ""))
    with open(cache_file(unique_id), "rb") as file:
        # print(f"Reading cached file {file.name}.")
        return pickle.load(file)


def write_cached(unique_id: str, data: Any) -> None:
    """
    Write data to the cache file associated with a unique ID.

    :param unique_id: A unique identifier to write data for.
    :param data: The data to be written to the cache.

    DocTests:
    >>> write_cached('test_write', {'value': 456})
    >>> read_cached('test_write') == {'value': 456}
    True
    """
    # print(f"Writing " + str(unique_id[:80]).replace("\n", ""))
    with open(cache_file(unique_id), "wb") as file:
        # print(f"Writing cached file {file.name}.")
        pickle.dump(data, file)


def cached_simplify(expr: sp.Basic) -> sp.Expr:
    r"""
    Simplify a given sympy expression using a cache to speed up repeated simplifications.

    :param expr: SymPy expression to be simplified.
    :return: Simplified SymPy expression.

    DocTests:
    >>> x = sp.symbols('x')
    >>> expr = sp.sympify("x**2 + 2*x + 1")
    >>> simplified_expr = cached_simplify(expr)
    >>> simplified_expr.equals((x + 1)**2)
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
            with open(cache_fle, "wb") as file:
                pickle.dump(simplified, file)
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
