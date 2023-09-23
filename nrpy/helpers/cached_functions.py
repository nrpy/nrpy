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
from typing import cast
from appdirs import user_cache_dir  # type: ignore
import sympy as sp


def cached_simplify(expr: sp.Basic) -> sp.Basic:
    if expr == sp.sympify(0):
        return cast(sp.Basic, sp.sympify(0))
    cache_dir = Path(user_cache_dir("nrpy"))
    try:
        pickle_expr = pickle.dumps(expr)
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file = cache_dir / (
            hashlib.sha256(pickle_expr).hexdigest() + ".nrpycache"
        )
        if cache_file.exists():
            with open(cache_file, "rb") as file:
                return cast(sp.Basic, pickle.load(file))
        else:
            simplified = sp.simplify(expr)
            with open(cache_file, "wb") as file:
                pickle.dump(simplified, file)
            return simplified
    except pickle.PicklingError:
        return sp.simplify(expr)
