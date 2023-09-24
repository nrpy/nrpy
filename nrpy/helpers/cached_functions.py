"""
Cache-enable versions of functions that are
* normally slow, and
* called with the same inputs many times during a development cycle.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
import hashlib
import pickle
from pathlib import Path
from typing import cast, Any
from appdirs import user_cache_dir  # type: ignore
import sympy as sp


def cache_file(unique_id: str) -> Path:
    """
    Generate a cache file path using a unique ID.

    :param unique_id: A unique identifier to be used for generating the file path.
    :return: The cache file path.
    """

    def get_hash(unique_id: str) -> str:
        """
        Generate a hash string for a given unique ID.

        :param unique_id: A unique identifier to be hashed.
        :return: The hash string.
        """
        return hashlib.sha256(unique_id.encode("utf-8")).hexdigest()

    return Path(user_cache_dir("nrpy")) / f"{get_hash(unique_id)}.nrpycache"


def is_cached(unique_id: str) -> bool:
    """
    Check if the file with the given unique ID exists in cache.

    :param unique_id: A unique identifier to check for.
    :return: True if the file exists, False otherwise.
    """
    return cache_file(unique_id).exists()


def read_cached(unique_id: str) -> Any:
    """
    Read the cached data associated with a unique ID.

    :param unique_id: A unique identifier to read data for.
    :return: The data read from the cache file.
    """
    with open(cache_file(unique_id), "rb") as file:
        # print(f"Reading cached file {file.name}.")
        return pickle.load(file)


def write_cached(unique_id: str, data: Any) -> None:
    """
    Write data to the cache file associated with a unique ID.

    :param unique_id: A unique identifier to write data for.
    :param data: The data to be written to the cache.
    """
    with open(cache_file(unique_id), "wb") as file:
        # print(f"Writing cached file {file.name}.")
        pickle.dump(data, file)


def cached_simplify(expr: sp.Basic) -> sp.Expr:
    """
    Simplify a given sympy expression using a cache to speed up repeated simplifications.

    :param expr: Sympy expression to be simplified.
    :return: Simplified sympy expression.
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
