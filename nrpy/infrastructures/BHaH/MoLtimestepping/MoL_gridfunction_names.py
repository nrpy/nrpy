# nrpy/infrastructures/BHaH/MoLtimestepping/MoL_gridfunction_names.py
"""
Functions for naming the gridfunctions (y_n_gfs, k_i gfs, etc.) based on the chosen Method of Lines method.
Also checks whether a Butcher table is diagonal.

Authors: Zachariah B. Etienne (lead maintainer)
         zachetie **at** gmail **dot* com
         Samuel Tootle (GPU support in NRPy2)
         Brandon Clark (original, NRPy1 version)
"""

from typing import Dict, List, Tuple, Union

import sympy as sp


def is_diagonal_Butcher(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    key: str,
) -> bool:
    """
    Check if a given Butcher table (from Butcher_dict) is diagonal.

    :param Butcher_dict: A dictionary containing Butcher tables. Each key maps to a table (list of lists).
    :param key: The key to retrieve the Butcher table from the Butcher_dict.

    :return: True if the table is diagonal, False otherwise.

    Note:
    A Butcher table is diagonal if all its non-diagonal elements are zero.
    """
    Butcher = Butcher_dict[key][0]
    L = len(Butcher) - 1

    for i in range(L):
        for j in range(1, i):
            if Butcher[i][j] != sp.sympify(0):
                return False
    return True


def generate_gridfunction_names(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    MoL_method: str = "RK4",
) -> Tuple[str, List[str], str, Union[str, None]]:
    """
    Generate gridfunction names for the specified Method of Lines (MoL) method.
    Used for setting up MoL_malloc, MoL_step_forward_in_time, MoL_free, and BHaH_defines.h

    :param Butcher_dict: Dictionary of Butcher tables for the MoL method.
    :param MoL_method: The MoL method to generate gridfunction names for.
    :return: A tuple containing y_n_gridfunctions, non_y_n_gridfunctions_list,
             diagnostic_gridfunctions_point_to, and diagnostic_gridfunctions2_point_to.

    Doctests:
    >>> from nrpy.infrastructures.BHaH.MoLtimestepping.RK_Butcher_Table_Dictionary import generate_Butcher_tables
    >>> Butcher_dict = generate_Butcher_tables()
    >>> generate_gridfunction_names(Butcher_dict, "RK2 Heun")
    ('y_n_gfs', ['y_nplus1_running_total_gfs', 'k_odd_gfs', 'k_even_gfs', 'auxevol_gfs'], 'y_nplus1_running_total_gfs', 'k_odd_gfs')
    >>> generate_gridfunction_names(Butcher_dict, "RK3")
    ('y_n_gfs', ['next_y_input_gfs', 'k1_gfs', 'k2_gfs', 'k3_gfs', 'auxevol_gfs'], 'k1_gfs', 'k2_gfs')
    >>> generate_gridfunction_names(Butcher_dict, "RK3 Ralston")
    ('y_n_gfs', ['k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs', 'k2_or_y_nplus_a32_k2_gfs', 'auxevol_gfs'], 'k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs', 'k2_or_y_nplus_a32_k2_gfs')
    >>> generate_gridfunction_names(Butcher_dict, "RK4")
    ('y_n_gfs', ['y_nplus1_running_total_gfs', 'k_odd_gfs', 'k_even_gfs', 'auxevol_gfs'], 'y_nplus1_running_total_gfs', 'k_odd_gfs')
    """
    y_n_gridfunctions = "y_n_gfs"
    non_y_n_gridfunctions_list = []

    if is_diagonal_Butcher(Butcher_dict, MoL_method) and "RK3" in MoL_method:
        non_y_n_gridfunctions_list.extend(
            [
                "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs",
                "k2_or_y_nplus_a32_k2_gfs",
            ]
        )
        diagnostic_gridfunctions_point_to = (
            "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs"
        )
        diagnostic_gridfunctions2_point_to = "k2_or_y_nplus_a32_k2_gfs"
    else:
        if not is_diagonal_Butcher(Butcher_dict, MoL_method):
            num_k = len(Butcher_dict[MoL_method][0]) - 1
            non_y_n_gridfunctions_list.append("next_y_input_gfs")
            for i in range(num_k):
                non_y_n_gridfunctions_list.append(f"k{i + 1}_gfs")
            diagnostic_gridfunctions_point_to = "k1_gfs"
            if "k2_gfs" in non_y_n_gridfunctions_list:
                diagnostic_gridfunctions2_point_to = "k2_gfs"
            else:
                print(
                    "MoL WARNING: No gridfunction group available for diagnostic_output_gfs2"
                )
                diagnostic_gridfunctions2_point_to = None
        else:
            non_y_n_gridfunctions_list.append("y_nplus1_running_total_gfs")
            if MoL_method != "Euler":
                non_y_n_gridfunctions_list.extend(["k_odd_gfs", "k_even_gfs"])
            diagnostic_gridfunctions_point_to = "y_nplus1_running_total_gfs"
            diagnostic_gridfunctions2_point_to = "k_odd_gfs"

    non_y_n_gridfunctions_list.append("auxevol_gfs")

    return (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        diagnostic_gridfunctions_point_to,
        diagnostic_gridfunctions2_point_to,
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
