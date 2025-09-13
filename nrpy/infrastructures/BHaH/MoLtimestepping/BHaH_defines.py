"""
A one-stop function to register MoL's contribution to the BHaH_defines.h file.

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot** com
"""

from typing import Dict, List, Tuple, Union

import nrpy.params as par
from nrpy.infrastructures import BHaH


def register_BHaH_defines_h(
    Butcher_dict: Dict[str, Tuple[List[List[Union[int, str]]], int]],
    MoL_method: str,
) -> None:
    """
    Register MoL's contribution to the BHaH_defines.h file.

    :param Butcher_dict: Dictionary of Butcher tables for the MoL method.
    :param MoL_method: Method for the Method of Lines.
    """
    parallelization = par.parval_from_str("parallelization")
    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        _diag_pt,
        _diag_pt2,
    ) = BHaH.MoLtimestepping.gridfunction_names.generate_gridfunction_names(
        Butcher_dict, MoL_method=MoL_method
    )

    BHaH_MoL_body: str = (
        "typedef struct __MoL_gridfunctions_struct__ {\n"
        + f"REAL *restrict {y_n_gridfunctions};\n"
        + "".join(f"REAL *restrict {gfs};\n" for gfs in non_y_n_gridfunctions_list)
        + r"""REAL *restrict diagnostic_output_gfs;
REAL *restrict diagnostic_output_gfs2;
} MoL_gridfunctions_struct;
"""
    )
    if parallelization != "openmp":
        BHaH_MoL_body = BHaH_MoL_body.replace("REAL *restrict ", "REAL * ")

    BHaH.BHaH_defines_h.register_BHaH_defines(__name__, BHaH_MoL_body)
