# nrpy/infrastructures/Dendro/runtime/parameters.py
"""
Register the generated runtime-parameter CFunctions for a Dendro solver.

Every name, type, default and description is read from the NRPy CodeParameter
registry at the point of use, as BHaH's ``BHaH_defines_h`` reads
``par.glb_code_params_dict``.  There is no hand-maintained parameter table and
no use closure: the whole registry is emitted, exactly as BHaH emits it.

The CFunctions are registered *without* a Dendro role: the host context calls
them, and they are not scheduled as numerical kernels.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Any, List

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.Dendro.CodeParameters import array_size, c_type

PARAMETERS_SUBDIRECTORY = "generated/src/parameters"


def emitted_parameter_names() -> List[str]:
    """
    Return the registered CodeParameter names the generated struct carries.

    The order and the ``#define`` exclusion match
    :func:`nrpy.infrastructures.Dendro.CodeParameters.output_Dendro_parameters_h`,
    so a generated body and the generated struct cannot disagree about which
    members exist.

    :return: Sorted CodeParameter names, ``#define`` parameters excluded.
    """
    return [
        cp_name
        for cp_name, code_param in sorted(par.glb_code_params_dict.items())
        if code_param.cparam_type != "#define"
    ]


def _default_statement(cp_name: str, scalar_type: str) -> str:
    """
    Render the statement that sets one registered CodeParameter's default.

    A ``char`` array is filled with ``snprintf``, which guarantees null
    termination, exactly as BHaH's ``CodeParameters`` emitter does.  Numeric
    array parameters are left value-initialized: a registered scalar default
    supplies no element values, and inventing them would be worse than zero.

    :param cp_name: CodeParameter name.
    :param scalar_type: Registered Dendro scalar alias.
    :return: One C++ statement, or the empty string when nothing is set.
    """
    code_param = par.glb_code_params_dict[cp_name]
    value: Any = code_param.defaultvalue
    cparam_type = code_param.cparam_type
    size = array_size(cparam_type)
    if size is not None:
        if c_type(cparam_type, scalar_type) != "char":
            return ""
        return f'std::snprintf(params.{cp_name}, {size}, "%s", "{value}");'
    if cparam_type == "bool":
        return f"params.{cp_name} = {'true' if value else 'false'};"
    if cparam_type == "int":
        return f"params.{cp_name} = {int(value)};"
    if cparam_type in ("REAL", "DendroScalar"):
        return f"params.{cp_name} = static_cast<{scalar_type}>({float(value)!r});"
    return f"params.{cp_name} = {float(value)!r};"


def register_CFunctions_parameters(solver_namespace: str, solver_stem: str) -> None:
    """
    Register the generated parameter CFunctions.

    Every parameter CFunction takes the generated parameter table by reference
    (``<solver_namespace>::generated::GeneratedParams& params``).  The table
    object is owned by the host context, so the generated bodies carry no host
    storage.

    Call this after every scientific CFunction is registered.

    :param solver_namespace: Solver namespace, following Dendro's lowercase
        formulation habit (``namespace bssn``).
    :param solver_stem: Lowercase formulation stem for the emitted CFunction
        names, following Dendro's habit of naming solver symbols for the
        formulation.
    """
    scalar_type = str(par.parval_from_str("Dendro_scalar_type"))
    params_type = f"{solver_namespace}::generated::GeneratedParams"
    names = emitted_parameter_names()

    set_lines: List[str] = [f"params = {params_type}{{}};"]
    for cp_name in names:
        statement = _default_statement(cp_name, scalar_type)
        if statement:
            set_lines.append(statement)
    cfc.register_CFunction(
        subdirectory=PARAMETERS_SUBDIRECTORY,
        desc="Generated parameter defaults, from the registered CodeParameters.",
        includes=["<cstdio>"],
        name=f"{solver_stem}_params_set_defaults",
        params=f"{params_type}& params",
        body="\n".join(set_lines),
    )

    validate_lines: List[str] = ["bool ok = true;"]
    for cp_name in names:
        cparam_type = par.glb_code_params_dict[cp_name].cparam_type
        if array_size(cparam_type) is not None:
            continue
        if c_type(cparam_type, scalar_type) in ("bool", "int", "char"):
            continue
        validate_lines.append(
            f"if (!std::isfinite(params.{cp_name})) {{ ok = false; }}"
        )
    validate_lines.append("return ok;")
    cfc.register_CFunction(
        subdirectory=PARAMETERS_SUBDIRECTORY,
        desc="Generated parameter validation: finite checks for floating point parameters.",
        cfunc_type="bool",
        includes=["<cmath>"],
        name=f"{solver_stem}_params_validate",
        params=f"const {params_type}& params",
        body="\n".join(validate_lines),
    )

    print_lines: List[str] = [f'std::printf("{solver_stem} effective parameters:\\n");']
    for cp_name in names:
        cparam_type = par.glb_code_params_dict[cp_name].cparam_type
        base = c_type(cparam_type, scalar_type)
        if array_size(cparam_type) is not None:
            if base == "char":
                print_lines.append(
                    f'std::printf("  {cp_name} = %s\\n", params.{cp_name});'
                )
            continue
        if base in ("bool", "int"):
            print_lines.append(
                f'std::printf("  {cp_name} = %d\\n", (int) params.{cp_name});'
            )
        else:
            print_lines.append(
                f'std::printf("  {cp_name} = %g\\n", (double) params.{cp_name});'
            )
    cfc.register_CFunction(
        subdirectory=PARAMETERS_SUBDIRECTORY,
        desc="Generated effective-parameter printout.",
        includes=["<cstdio>"],
        name=f"{solver_stem}_params_print_effective",
        params=f"const {params_type}& params",
        body="\n".join(print_lines),
    )

    # No parameter-file binding exists in this profile: the host parser binds
    # the same generated struct once Dendrolib is pinned.  Returning nonzero
    # keeps a supplied parameter file from being silently ignored.
    parse_lines: List[str] = [
        "// This profile has no parameter-file binding: the host parser binds",
        "// the same generated struct once Dendrolib is pinned.  Refusing the",
        "// file keeps a caller from believing its values were applied.",
        "(void) params;",
        "(void) parfile_path;",
        "return 1;",
    ]
    cfc.register_CFunction(
        subdirectory=PARAMETERS_SUBDIRECTORY,
        desc=(
            "Generated parameter-file binding, unimplemented in this profile: "
            "refuses a supplied parameter file."
        ),
        cfunc_type="int",
        name=f"{solver_stem}_params_parse_file",
        params=f"{params_type}& params, const char* parfile_path",
        body="\n".join(parse_lines),
    )
