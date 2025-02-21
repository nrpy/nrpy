"""
Manage registration and storage of data stored within griddata_struct and commondata_struct.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, List

import nrpy.c_function as cfc
import nrpy.helpers.gpu.utilities as gpu_utils
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
        dictionary: Dict[str, List[GridCommonData]],
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


def register_CFunction_griddata_free(
    enable_rfm_precompute: bool,
    enable_CurviBCs: bool,
    enable_bhahaha: bool = False,
    parallelization: str = "openmp",
) -> None:
    """
    Register the C function griddata_free() to free all memory within the griddata struct.

    :param enable_rfm_precompute: A flag to enable/disable rfm_precompute_free within the C function body.
    :param enable_CurviBCs: Whether to free CurviBCs within the C function body.
    :param enable_bhahaha: Whether to enable freeing of BHaHAHA memory.
    :param parallelization: Parallelization method to use. Default is "openmp".

    :raises ValueError: If BHaHAHA is not supported in the parallelization mode.
    """
    desc = """Free all memory within the griddata struct,
except perhaps non_y_n_gfs (e.g., after a regrid, in which non_y_n_gfs are freed first)."""
    cfunc_type = "void"
    name = "griddata_free"
    params = (
        "const commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_host, const bool free_non_y_n_gfs_and_core_griddata_pointers"
        if parallelization == "cuda"
        else "const commondata_struct *restrict commondata, griddata_struct *restrict griddata, const bool free_non_y_n_gfs_and_core_griddata_pointers"
    )
    body = ""
    free_func = gpu_utils.get_memory_free_function(parallelization)
    if enable_bhahaha and parallelization == "openmp":
        body += r"""  // Free BHaHAHA memory.
  for (int which_horizon = 0; which_horizon < commondata->bah_max_num_horizons; which_horizon++) {
    free(commondata->bhahaha_params_and_data[which_horizon].prev_horizon_m1);
    free(commondata->bhahaha_params_and_data[which_horizon].prev_horizon_m2);
    free(commondata->bhahaha_params_and_data[which_horizon].prev_horizon_m3);
  }
""".replace(
            "free(", f"{free_func}("
        )
    elif enable_bhahaha and parallelization != "openmp":
        raise ValueError(
            "BHaHAHA is not yet supported in parallelization mode: " + parallelization
        )
    body += r"""  // Free memory allocated inside griddata[].
  for(int grid=0;grid<commondata->NUMGRIDS;grid++) {
"""
    if enable_rfm_precompute:
        body += "  rfm_precompute_free(commondata, &griddata[grid].params, griddata[grid].rfmstruct);\n"
        body += "  free(griddata[grid].rfmstruct);\n".replace("free(", f"{free_func}(")
        body += f'  {gpu_utils.get_check_errors_str(parallelization, free_func, opt_msg="Free: rfmstruct failed")}'
    if enable_CurviBCs:
        body += rf"""
  free(griddata[grid].bcstruct.inner_bc_array);
  {gpu_utils.get_check_errors_str(parallelization, free_func, opt_msg="Free: bcstruct.inner_bc_array failed")}
  for(int ng=0;ng<NGHOSTS*3;ng++) {{
      free(griddata[grid].bcstruct.pure_outer_bc_array[ng]);
      {gpu_utils.get_check_errors_str(parallelization, free_func, opt_msg="Free: bcstruct.pure_outer_bc_array failed")}
}}
""".replace(
            "free(", f"{free_func}("
        )
    body += rf"""

  MoL_free_memory_y_n_gfs(&griddata[grid].gridfuncs);
  if(free_non_y_n_gfs_and_core_griddata_pointers) {{
    MoL_free_memory_non_y_n_gfs(&griddata[grid].gridfuncs);
  }}
  for(int i=0;i<3;i++) {{
    free(griddata[grid].xx[i]);
    {gpu_utils.get_check_errors_str(parallelization, free_func, opt_msg="Free: grid.XX failed")}
""".replace(
        "free(", f"{free_func}("
    )
    body += "free(griddata_host[grid].xx[i]);\n" if parallelization == "cuda" else ""
    body += """}
} // END for(int grid=0;grid<commondata->NUMGRIDS;grid++)
"""
    body += rf"""if(free_non_y_n_gfs_and_core_griddata_pointers) {{
        free(griddata);
        {"free(griddata_host);" if parallelization == "cuda" else ""}
    }}"""
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
