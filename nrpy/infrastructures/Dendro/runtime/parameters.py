# nrpy/infrastructures/Dendro/runtime/parameters.py
"""
Generated runtime parameter CFunctions (whitepaper section 6.4).

After all scientific CFunctions have established the complete use closure,
this module registers the parameter default/validation/print CFunctions.
Their bodies are generated from the frozen CodeParameter records: every name,
type, default, and description is read from the live NRPy CodeParameter
registry (no hand-maintained parameter table).

The CFunctions are registered *without* a Dendro role: they are called by the
host context (the PR 7 mock vehicle and, after the I0-1 gates, the real
Dendro-GR context), not scheduled as numerical kernels.
"""

from typing import Any, List, Tuple

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.Dendro.CodeParameters_output import c_type
from nrpy.infrastructures.Dendro.freeze import body_uses_symbol

PARAMS_SUBDIRECTORY = "generated/src/parameters"


def _used_cparams_in_order() -> Tuple[str, ...]:
    """
    Return the used CodeParameter names in deterministic order.

    The use closure is computed token-aware over the registered CFunction
    bodies (section 6.2), exactly as freeze computes it, so the parameter
    CFunctions and the frozen snapshot describe the same set.

    :return: Sorted tuple of used CodeParameter names.
    """
    bodies = {
        name: cfc.CFunction_dict[name].full_function
        for name in sorted(cfc.CFunction_dict)
    }
    used = set()
    for body in bodies.values():
        for cp_name in par.glb_code_params_dict:
            if body_uses_symbol(body, cp_name):
                used.add(cp_name)
    return tuple(sorted(used))


def _cparam_cxx_type(cp_name: str, scalar_type: str) -> str:
    """
    Map one registered CodeParameter type to its generated C++ type.

    :param cp_name: Registered CodeParameter name.
    :param scalar_type: Registered Dendro scalar alias.
    :return: The generated C++ type spelling (section 6.5).
    """
    return c_type(par.glb_code_params_dict[cp_name].cparam_type, scalar_type)


def _default_literal(cp_name: str, scalar_type: str) -> str:
    """
    Render the C++ literal for a registered CodeParameter default.

    :param cp_name: CodeParameter name.
    :param scalar_type: Registered Dendro scalar alias.
    :return: A C++ literal expression (REAL-family defaults are rendered with
        a double literal plus cast to the scalar alias, so the generated text
        does not depend on the alias spelling for formatting).
    """
    cp = par.glb_code_params_dict[cp_name]
    value: Any = cp.defaultvalue
    cparam_type = cp.cparam_type
    if cparam_type == "bool":
        return "true" if value else "false"
    if cparam_type == "int":
        return str(int(value))
    if cparam_type == "#define":
        return str(value)
    # REAL / float / double / DendroScalar and array types: float literal.
    if cparam_type in ("REAL", "DendroScalar"):
        return f"static_cast<{scalar_type}>({float(value)!r})"
    return repr(float(value))


def register_parameter_CFunctions_last() -> None:
    """
    Register the generated parameter CFunctions.

    Every parameter CFunction takes the generated parameter table by
    reference (``fccz4::generated::GeneratedParams& params``); the table
    object itself is owned by the host context (mock vehicle) or the real
    Dendro-GR context, so the generated bodies contain no host storage.

    Must be called after every scientific CFunction is registered (so the
    token-aware use closure is complete) and before freeze.

    """
    scalar_type = str(par.parval_from_str("Dendro_scalar_type"))
    used = _used_cparams_in_order()

    # --- set_defaults -----------------------------------------------------
    set_lines: List[str] = ["params = fccz4::generated::GeneratedParams{};"]
    for cp_name in used:
        cxx_type = _cparam_cxx_type(cp_name, scalar_type)
        if cxx_type.startswith("std::array"):
            # Element defaults are unavailable from a registered scalar
            # default; leave the array value-initialized (zero-initialized).
            continue
        set_lines.append(
            f"params.{cp_name} = {_default_literal(cp_name, scalar_type)};"
        )
    cfc.register_CFunction(
        subdirectory=PARAMS_SUBDIRECTORY,
        desc="Generated fccz4 parameter defaults (from the frozen CodeParameter records).",
        name="fccz4_params_set_defaults",
        params="fccz4::generated::GeneratedParams& params",
        body="\n".join(set_lines),
    )

    # --- validate -----------------------------------------------------------
    validate_lines: List[str] = ["bool ok = true;"]
    for cp_name in used:
        cxx_type = _cparam_cxx_type(cp_name, scalar_type)
        if cxx_type == "bool" or cxx_type == "int" or cxx_type.startswith("std::array"):
            continue
        # Finite-check for every floating point parameter (section 13.7:
        # all required runtime parameters must be parsed and validated).
        validate_lines.append(
            f"if (!std::isfinite(params.{cp_name})) {{ ok = false; }}"
        )
    validate_lines.append("return ok;")
    cfc.register_CFunction(
        subdirectory=PARAMS_SUBDIRECTORY,
        desc="Generated fccz4 parameter validation (finite checks for used floating point parameters).",
        cfunc_type="bool",
        includes=["<cmath>"],
        name="fccz4_params_validate",
        params="const fccz4::generated::GeneratedParams& params",
        body="\n".join(validate_lines),
    )

    # --- print_effective -----------------------------------------------------
    print_lines: List[str] = ['std::printf("fccz4 effective parameters:\\n");']
    for cp_name in used:
        cxx_type = _cparam_cxx_type(cp_name, scalar_type)
        if cxx_type.startswith("std::array"):
            continue
        if cxx_type in ("bool", "int"):
            print_lines.append(
                f'std::printf("  {cp_name} = %d\\n", (int) params.{cp_name});'
            )
        else:
            print_lines.append(
                f'std::printf("  {cp_name} = %g\\n", (double) params.{cp_name});'
            )
    cfc.register_CFunction(
        subdirectory=PARAMS_SUBDIRECTORY,
        desc="Generated fccz4 effective-parameter printout.",
        includes=["<cstdio>"],
        name="fccz4_params_print_effective",
        params="const fccz4::generated::GeneratedParams& params",
        body="\n".join(print_lines),
    )

    # --- parse_toml -----------------------------------------------------------
    # The production TOML parser binds the generated descriptors through the
    # host TOML library after the I0-1 gates.  Until then there is no binding,
    # so this returns nonzero: a caller that supplied a parameter file must be
    # told its values are not being applied, never left believing they were
    # (section 6.4; section 13.7 requires every runtime parameter to be parsed
    # and validated before the first step).
    parse_lines: List[str] = [
        "// No TOML binding exists in this generated profile: the real host",
        "// parser binds the same generated descriptor table after the pinned",
        "// Dendrolib gates (section 6.4).  Returning nonzero keeps a supplied",
        "// parameter file from being silently ignored.",
        "(void) params;",
        "(void) toml_path;",
        "return 1;",
    ]
    cfc.register_CFunction(
        subdirectory=PARAMS_SUBDIRECTORY,
        desc=(
            "Generated fccz4 TOML binding (unimplemented in this profile: "
            "refuses a supplied parameter file)."
        ),
        cfunc_type="int",
        name="fccz4_params_parse_toml",
        params="fccz4::generated::GeneratedParams& params, const char* toml_path",
        body="\n".join(parse_lines),
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
