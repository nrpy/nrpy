"""
Module for constructing param.ccl for Cactus thorns.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import List
from pathlib import Path
import nrpy.params as par
import nrpy.c_function as cfc


def construct_param_ccl(
    project_dir: str,
    thorn_name: str,
    shares_extends_str: str = "",
) -> List[str]:
    """
    Generate the content for the param.ccl file and write it to disk.

    :param project_dir: The directory where the project resides.
    :param thorn_name: Name of the thorn.
    :param shares_extends_str: String that specifies any shared or extended Cactus thorns.
    :return: A list of parameter names that are registered to the param.ccl file.
    """
    paramccl_str = f"""# This param.ccl file was automatically generated by NRPy+.
#   You are advised against modifying it directly; instead
#   modify the Python code that generates it.
{shares_extends_str}

restricted:
"""
    CParams_registered_to_params_ccl: List[str] = []

    for CFunction in cfc.CFunction_dict.values():
        if (
            CFunction.ET_thorn_name == thorn_name
            and CFunction.ET_current_thorn_CodeParams_used
        ):
            for CPname in CFunction.ET_current_thorn_CodeParams_used:
                # only declare parameters once
                if CPname not in CParams_registered_to_params_ccl:
                    CParam = par.glb_code_params_dict[CPname]
                    paramccl_str += f'{CParam.c_type_alias} {CParam.name} "(see NRPy+ for parameter definition)"\n'
                    paramccl_str += "{\n"
                    paramccl_str += ' *:* :: "All values accepted. NRPy+ does not restrict the allowed ranges of parameters yet."\n'
                    paramccl_str += f"}} {CParam.defaultvalue}\n\n"
                    CParams_registered_to_params_ccl += [CPname]
    output_Path = Path(project_dir) / thorn_name
    output_Path.mkdir(parents=True, exist_ok=True)
    with open(output_Path / "param.ccl", "w", encoding="utf-8") as file:
        file.write(paramccl_str)
    return CParams_registered_to_params_ccl
