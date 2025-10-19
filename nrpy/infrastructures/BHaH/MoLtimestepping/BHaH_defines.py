"""
A one-stop function to register MoL's contribution to the BHaH_defines.h file.

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot** com
"""

from typing import Dict, List, Tuple, Union

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
    intermediate_stage_gfs_gridfunctions_list = BHaH.MoLtimestepping.rk_butcher_table_dictionary.intermediate_stage_gf_names_list(
        Butcher_dict, MoL_method=MoL_method
    )
    BHaH_MoL_body: str = (
        "typedef struct __MoL_gridfunctions_struct__ {\n"
        + "REAL *y_n_gfs;\n"
        + "".join(
            f"REAL *{gfs};\n" for gfs in intermediate_stage_gfs_gridfunctions_list
        )
        + """REAL *auxevol_gfs;
} MoL_gridfunctions_struct;
"""
    )

    BHaH.BHaH_defines_h.register_BHaH_defines(__name__, BHaH_MoL_body)
