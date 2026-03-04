"""Structural contracts for rotation C-function generator modules."""

from __future__ import annotations

import ast
from pathlib import Path

MODULE_CONTRACTS = {
    "rotate_BSSN_Cartesian_basis_from_axis_angle.py": "register_CFunction_rotate_BSSN_Cartesian_basis_from_axis_angle",
    "rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src.py": "register_CFunction_rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src",
}


def _top_level_function_names(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return [node.name for node in tree.body if isinstance(node, ast.FunctionDef)]


def _assert_module_contract(path: Path, expected_register_name: str) -> None:
    fn_names = _top_level_function_names(path)
    register_names = [
        name for name in fn_names if name.startswith("register_CFunction_")
    ]
    if register_names != [expected_register_name]:
        raise AssertionError(
            f"{path.name} must define exactly one register_CFunction_* named "
            f"{expected_register_name}; got {register_names}"
        )

    parity_names = [name for name in fn_names if "parity" in name]
    if parity_names:
        raise AssertionError(
            f"{path.name} must not define parity test helpers; found {parity_names}"
        )


def run_rotation_module_contracts() -> bool:
    r"""
    Enforce one-file-per-register-function and no test helpers in production files.

    Doctests:
    >>> run_rotation_module_contracts()
    True
    """
    rotation_dir = Path(__file__).resolve().parent.parent
    for module_name, register_name in MODULE_CONTRACTS.items():
        _assert_module_contract(rotation_dir / module_name, register_name)
    return True


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
