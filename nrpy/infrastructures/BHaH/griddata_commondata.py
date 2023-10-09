"""
Manage registration and storage of data stored within griddata_struct and commondata_struct.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import Dict, List

import nrpy.params as par


class GridCommonData:
    """
    Represent grid data with a module and declaration information.

    :param module: Module associated with the grid data
    :param c_declaration: C code declaration; e.g., "struct params_struct"
    :param description: Description of the module (default is empty string)

    :raises ValueError: If the declaration contains a semicolon
    """

    def __init__(self, module: str, c_declaration: str, description: str = "") -> None:
        self.module = module
        self.c_declaration = c_declaration
        self.description = description
        if ";" in c_declaration:
            raise ValueError("GridData.c_declaration cannot have a semicolon inside.")


def register_griddata_commondata(
    module: str, c_declaration: str, description: str = "", is_commondata: bool = False
) -> None:
    """
    Register grid data into the global dictionary `glb_griddata_struct_dict`.

    :param module: Module associated with the grid data
    :param c_declaration: C code declaration; e.g., "struct params_struct"
    :param description: Description of the module (default is empty string)
    :param is_commondata: Whether to register as commondata (default is False)

    :raises ValueError: If the same declaration is registered into the same module twice

    Doctest:
    >>> register_griddata_commondata("my_module", "struct my_module", "my_module's description")
    >>> par.glb_extras_dict["griddata_struct"]["my_module"][0].c_declaration
    'struct my_module'
    >>> print(par.glb_extras_dict["griddata_struct"]["my_module"][0].description)
    my_module's description
    >>> register_griddata_commondata("my_module", "struct my_module", "my description")  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    ValueError: Cannot register the same declaration (struct my_module) into glb_griddata_struct_dict[my_module] twice.
    """

    def register_griddata_or_commondata(
        dictionary: Dict[str, List[GridCommonData]]
    ) -> None:
        if module in dictionary:
            if any(gd.c_declaration == c_declaration for gd in dictionary[module]):
                raise ValueError(
                    f"Cannot register the same declaration ({c_declaration}) into glb_griddata_struct_dict[{module}] twice."
                )
            dictionary[module].append(
                GridCommonData(module, c_declaration, description)
            )
        else:
            dictionary[module] = [GridCommonData(module, c_declaration, description)]

    if is_commondata:
        # commondata_struct stores data common to all grids
        if "commondata_struct" not in par.glb_extras_dict:
            par.glb_extras_dict["commondata_struct"] = {}
        register_griddata_or_commondata(par.glb_extras_dict["commondata_struct"])
    else:
        # griddata_struct stores data specific to each grid
        if "griddata_struct" not in par.glb_extras_dict:
            par.glb_extras_dict["griddata_struct"] = {}
        register_griddata_or_commondata(par.glb_extras_dict["griddata_struct"])


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
