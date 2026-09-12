# nrpy/infrastructures/Dendro/CodeParameters.py
"""
Emit the parameter header and parameter CFunctions.

Every name, type and default is read from the NRPy CodeParameter registry, as
BHaH's ``CodeParameters.py`` does.  Every registered parameter whose
``cparam_type`` is not ``#define`` is emitted -- there is no use closure, so the
generated struct is the whole registry and a caller that registers a parameter
gets it in the table; no physics parameter table is authored here and no
equation module is imported.

The parameter CFunctions (set defaults, parse a file, validate, print) are
registered from here as well, as BHaH registers ``params_struct_set_to_default``
from its own ``CodeParameters``.  They carry no Dendro role: the host context
calls them, and they are not scheduled as numerical kernels.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.header_guards import header_guard

BANNER = generated_file_banner()


def c_type(cparam_type: str) -> str:
    """
    Map a registered CodeParameter type to its generated C++ base type.

    ``REAL`` maps to the registered Dendro scalar alias; ``float``, ``double``,
    ``int``, ``bool`` and ``char`` map to the same spelled built-in types.  For
    an array type the base is mapped and the length is carried separately, as
    BHaH's ``BHaH_defines_h`` does, because a C array declares its length after
    the member name.

    :param cparam_type: Registered ``cparam_type`` string.
    :return: The generated C++ base-type spelling.
    :raises ValueError: If the type has no mapping.

    Doctests:
    >>> c_type("REAL")
    'DendroScalar'
    >>> c_type("char[100]")
    'char'
    """
    base, _size, _is_array = par.parse_cparam_type(cparam_type)
    if base in ("REAL", gri.DENDRO_SCALAR_TYPE):
        return gri.DENDRO_SCALAR_TYPE
    if base in ("float", "double", "int", "bool", "char"):
        return base
    if base == "#define":
        raise ValueError(
            "`#define` CodeParameters are compile-time constants and are never "
            "`params_struct` members."
        )
    raise ValueError(f"CodeParameter type {cparam_type!r} has no generated mapping.")


def member_declaration(cp_name: str, cparam_type: str) -> str:
    """
    Render one ``params_struct`` member declaration.

    :param cp_name: Registered CodeParameter name.
    :param cparam_type: Registered ``cparam_type`` string.
    :return: The C++ member declaration, without indentation.

    Doctests:
    >>> member_declaration("eta", "REAL")
    'DendroScalar eta;'
    >>> member_declaration("steps", "int")
    'int steps;'
    >>> member_declaration("enable_filter", "bool")
    'bool enable_filter;'
    >>> member_declaration("CoordSystemName", "char[100]")
    'char CoordSystemName[100];'
    """
    _base, size, _is_array = par.parse_cparam_type(cparam_type)
    base = c_type(cparam_type)
    return f"{base} {cp_name}[{size}];" if size else f"{base} {cp_name};"


def output_parameters_h(solver_stem: str, solver_namespace: str) -> str:
    """
    Emit the generated parameter struct header.

    One member per registered CodeParameter in sorted name order, as BHaH's
    ``BHaH_defines_h`` iterates ``par.glb_code_params_dict``.  No default is
    baked into the struct: defaults live in the registered parameter-defaults
    CFunction and in the sample parameter file.

    :param solver_stem: Lowercase formulation stem for emitted header names.
    :param solver_namespace: Solver namespace, following Dendro's lowercase
        formulation habit (``namespace bssn``).
    :return: The complete C++ header text.

    Doctests:
    >>> header = output_parameters_h("bssn", "bssn")
    >>> "#ifndef BSSN_PARAMETERS_H" in header
    True
    >>> header.rstrip().endswith("#endif  // BSSN_PARAMETERS_H")
    True
    >>> "namespace bssn::generated {" in header
    True
    >>> from nrpy.helpers.generic import clang_format
    >>> "}  // END NAMESPACE: bssn::generated" in clang_format(header)
    True
    """
    opening, closing = header_guard(f"{solver_stem}_parameters.h")
    lines: List[str] = [BANNER.rstrip("\n"), opening, ""]
    lines += ["", f'#include "{solver_stem}_types.h"', ""]
    lines.append(f"namespace {solver_namespace}::generated {{")
    lines += ["", "struct params_struct {"]
    emitted = 0
    for cp_name, code_param in sorted(par.glb_code_params_dict.items()):
        if code_param.cparam_type == "#define":
            continue
        lines.append("    " + member_declaration(cp_name, code_param.cparam_type))
        emitted += 1
    lines += ["};  // END STRUCT: params_struct", ""]
    lines.append(f"inline constexpr unsigned NUM_CODE_PARAMETERS = {emitted};")
    lines += [
        "",
        "// clang-format off",
        f"}}  // END NAMESPACE: {solver_namespace}::generated",
        "// clang-format on",
        "",
        closing,
        "",
    ]
    return "\n".join(lines)


PARAMETERS_SUBDIRECTORY = "generated/src/parameters"


def emitted_parameter_names() -> List[str]:
    """
    Return the registered CodeParameter names the generated struct carries.

    The order and the ``#define`` exclusion match
    :func:`nrpy.infrastructures.Dendro.CodeParameters.output_parameters_h`,
    so a generated body and the generated struct cannot disagree about which
    members exist.

    :return: Sorted CodeParameter names, ``#define`` parameters excluded.
    """
    return [
        cp_name
        for cp_name, code_param in sorted(par.glb_code_params_dict.items())
        if code_param.cparam_type != "#define"
    ]


def register_CFunctions_parameters(solver_stem: str, solver_namespace: str) -> None:
    """
    Register the generated parameter CFunctions.

    Every parameter CFunction takes the generated parameter table by reference
    (``<solver_namespace>::generated::params_struct& params``).  The table
    object is owned by the host context, so the generated bodies carry no host
    storage.

    Call this after every scientific CFunction is registered.

    :param solver_stem: Lowercase formulation stem for the emitted CFunction
        names, following Dendro's habit of naming solver symbols for the
        formulation.
    :param solver_namespace: Solver namespace, following Dendro's lowercase
        formulation habit (``namespace bssn``).
    """
    params_type = f"{solver_namespace}::generated::params_struct"
    names = emitted_parameter_names()
    subdirectory = PARAMETERS_SUBDIRECTORY
    includes = [f"{solver_stem}_defines.h"]

    # A char array is filled with snprintf, which guarantees null termination,
    # exactly as BHaH's CodeParameters emitter does.  A numeric array is set
    # element by element: core broadcasts a scalar default across the length at
    # registration (nrpy/params.py, CodeParameter.__init__), so the registered
    # default is already a list and leaving the member value-initialized would
    # emit zeros where core holds the default.
    set_lines: List[str] = [f"params = {params_type}{{}};"]
    for cp_name in names:
        code_param = par.glb_code_params_dict[cp_name]
        value = code_param.defaultvalue
        cparam_type = code_param.cparam_type
        # Core permits "unset" for CodeParameters that a later setup stage
        # must supply.  Leave those members value-initialized here instead of
        # trying to coerce the sentinel into a numeric literal.
        if value == "unset" or (
            isinstance(value, (list, tuple))
            and all(element == "unset" for element in value)
        ):
            continue
        base_type, size, _is_array = par.parse_cparam_type(cparam_type)
        if size is not None:
            if c_type(cparam_type) == "char":
                set_lines.append(
                    f'std::snprintf(params.{cp_name}, {size}, "%s", "{value}");'
                )
            else:
                elements = value
                for index, element in enumerate(elements[: int(size)]):
                    if base_type == "int":
                        set_lines.append(f"params.{cp_name}[{index}] = {int(element)};")
                    else:
                        set_lines.append(
                            f"params.{cp_name}[{index}] = "
                            f"static_cast<{gri.DENDRO_SCALAR_TYPE}>({float(element)!r});"
                        )
        elif cparam_type == "bool":
            set_lines.append(f"params.{cp_name} = {'true' if value else 'false'};")
        elif cparam_type == "int":
            set_lines.append(f"params.{cp_name} = {int(value)};")
        elif cparam_type in ("REAL", gri.DENDRO_SCALAR_TYPE):
            set_lines.append(
                f"params.{cp_name} = static_cast<{gri.DENDRO_SCALAR_TYPE}>({float(value)!r});"
            )
        else:
            set_lines.append(f"params.{cp_name} = {float(value)!r};")
    set_to_default_desc = (
        "Generated parameter defaults, from the registered CodeParameters."
    )
    set_to_default_cfunc_type = "void"
    set_to_default_name = f"{solver_stem}_params_struct_set_to_default"
    set_to_default_params = f"{params_type}& params"
    set_to_default_body = "\n".join(set_lines)
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=set_to_default_desc,
        cfunc_type=set_to_default_cfunc_type,
        name=set_to_default_name,
        params=set_to_default_params,
        body=set_to_default_body,
    )

    validate_lines: List[str] = ["bool ok = true;"]
    for cp_name in names:
        cparam_type = par.glb_code_params_dict[cp_name].cparam_type
        if par.parse_cparam_type(cparam_type)[2]:
            continue
        if c_type(cparam_type) in ("bool", "int", "char"):
            continue
        validate_lines.append(
            f"if (!std::isfinite(params.{cp_name})) {{ ok = false; }}"
        )
    validate_lines.append("return ok;")
    validate_desc = (
        "Generated parameter validation: finite checks for floating point parameters."
    )
    validate_cfunc_type = "bool"
    validate_name = f"{solver_stem}_params_validate"
    validate_params = f"const {params_type}& params"
    validate_body = "\n".join(validate_lines)
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=validate_desc,
        cfunc_type=validate_cfunc_type,
        name=validate_name,
        params=validate_params,
        body=validate_body,
    )

    print_lines: List[str] = [f'std::printf("{solver_stem} effective parameters:\\n");']
    for cp_name in names:
        cparam_type = par.glb_code_params_dict[cp_name].cparam_type
        base = c_type(cparam_type)
        if par.parse_cparam_type(cparam_type)[2]:
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
    print_effective_desc = "Generated effective-parameter printout."
    print_effective_cfunc_type = "void"
    print_effective_name = f"{solver_stem}_params_print_effective"
    print_effective_params = f"const {params_type}& params"
    print_effective_body = "\n".join(print_lines)
    cfc.register_CFunction(
        subdirectory=subdirectory,
        includes=includes,
        desc=print_effective_desc,
        cfunc_type=print_effective_cfunc_type,
        name=print_effective_name,
        params=print_effective_params,
        body=print_effective_body,
    )


def output_toml_bindings() -> str:
    """
    Emit checked TOML assignments for registered runtime parameters only.

    :return: C++ assignments for the generated parameter table.
    :raises ValueError: If a character parameter has no supported TOML binding.
    """
    lines: List[str] = []
    for name, parameter in sorted(par.glb_code_params_dict.items()):
        if not parameter.add_to_parfile:
            continue
        base, size, is_array = par.parse_cparam_type(parameter.cparam_type)
        mapped = c_type(parameter.cparam_type)
        if base == "char":
            raise ValueError("Dendro TOML character parameters are not supported.")
        lines.append(f'if (item.first == "{name}") {{')
        if is_array:
            lines.append(
                f"const auto values = toml::get<std::vector<{mapped}>>(item.second);"
            )
            lines.append(
                f'if (values.size() != {size}) throw std::runtime_error("wrong parameter array length: {name}");'
            )
            lines.append(f"for (unsigned i = 0; i < {size}; ++i) {{")
            if base not in ("int", "bool"):
                lines.append(
                    'if (!std::isfinite(values[i])) throw std::runtime_error("nonfinite parameter");'
                )
            lines.append(f"params.{name}[i] = values[i];")
            lines.append("} // END LOOP: assign parameter array elements")
        else:
            lines.append(f"params.{name} = toml::get<{mapped}>(item.second);")
        lines += ["continue;", "} // END IF: bind registered runtime parameter"]
    return "\n".join(lines)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
