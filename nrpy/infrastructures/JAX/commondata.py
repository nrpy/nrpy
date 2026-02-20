"""
Manages the registration of common data parameters and generation of a Commondata dataclass for JAX projects.

Author: Siddharth Mahesh
        sidmahesh **at** gmail **dot* com
"""

import textwrap
from typing import Any, Dict, List

# Dictionary to store registered Commondata parameters
commondata_params_dict: Dict[str, Dict[str, Any]] = {}


def register_commondata_param(
    name: str,
    dtype: str,
    default: Any,
    description: str = "",
) -> None:
    """
    Register a parameter for the Commondata dataclass.

    :param name: The name of the parameter.
    :param dtype: The data type of the parameter.
    :param default: The default value of the parameter.
    :param description: A description of the parameter.

    :raises ValueError: If the parameter name is already registered.

    Doctests:
    >>> register_commondata_param("param1", "float", 0.0, "An example parameter.")
    >>> "param1" in commondata_params_dict
    True
    >>> commondata_params_dict["param1"]['dtype']
    'float'
    >>> register_commondata_param("param1", "float", 0.0) # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
        ...
    ValueError: Parameter 'param1' is already registered in commondata_params_dict.
    """
    if name in commondata_params_dict:
        raise ValueError(
            f"Parameter '{name}' is already registered in commondata_params_dict."
        )
    commondata_params_dict[name] = {
        "dtype": dtype,
        "default": default,
        "description": description,
    }


def register_commondata_params(
    names: List[str],
    dtypes: List[str],
    defaults: List[Any],
    descriptions: List[str],
) -> None:
    """
    Register multiple parameters for the Commondata dataclass.

    :param names: A list of parameter names to register.
    :param dtypes: A list of data types for the parameters.
    :param defaults: A list of default values for the parameters.
    :param descriptions: A list of descriptions for the parameters.
    """
    for name, dtype, default, description in zip(names, dtypes, defaults, descriptions):
        register_commondata_param(name, dtype, default, description)


def generate_commondata_dataclass() -> str:
    """
    Generate the string content for the Commondata dataclass module.

    :return: A string containing the full content of the Commondata module.
    """
    imports = "from dataclasses import dataclass\n"

    fields = []
    for name, props in commondata_params_dict.items():
        field_str = f"    {name}: {props['dtype']} = {props['default']}"
        if props["description"]:
            field_str += f"  # {props['description']}"
        fields.append(field_str)

    class_body = textwrap.dedent("""
        @dataclass
        class Commondata:
            \"\"\"Dataclass to store common simulation parameters.\"\"\"
        """)

    if fields:
        class_body += "\n".join(fields)
    else:
        class_body += "    pass\n"

    return f"{imports}\n{class_body}"


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
