"""
Library of C functions for solving the BSSN equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Tuple

import nrpy.c_codegen as ccg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy.equations.general_relativity import psi4_tetrads
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)


def generate_CFunction_psi4_tetrad(
    CoordSystem: str,
    tetrad: str = "quasiKinnersley",
    use_metric_to_construct_unit_normal: bool = False,
) -> Tuple[str, str]:
    """
    Register C function for psi4 tetrad computations.

    :param CoordSystem: The coordinate system to be used.
    :param tetrad: The type of tetrad. Defaults to "quasiKinnersley".
    :param use_metric_to_construct_unit_normal: Whether to use the metric to construct the unit normal. Defaults to False.

    :return: Tuple containing the prefunction and prefunction launch code.
    """
    parallelization = par.parval_from_str("parallelization")
    desc = f"Compute {tetrad} tetrad for psi4, with use_metric_to_construct_unit_normal={use_metric_to_construct_unit_normal}"
    name = "psi4_tetrad"
    cfunc_type = "void"
    cfunc_decorators = "__device__" if parallelization == "cuda" else ""

    # Initialize psi4 tetrad
    psi4tet_class = psi4_tetrads.Psi4Tetrads(
        CoordSystem,
        enable_rfm_precompute=False,
        tetrad=tetrad,
        use_metric_to_construct_unit_normal=use_metric_to_construct_unit_normal,
    )
    list_of_vrnms = []
    list_of_exprs = []

    for i in range(4):
        list_of_vrnms.append("*mre4U" + str(i))
        list_of_exprs.append(psi4tet_class.mre4U[i])
    for i in range(4):
        list_of_vrnms.append("*mim4U" + str(i))
        list_of_exprs.append(psi4tet_class.mim4U[i])
    for i in range(4):
        list_of_vrnms.append("*n4U" + str(i))
        list_of_exprs.append(psi4tet_class.n4U[i])

    list_of_metricvarnames = ["cf"] + [
        f"hDD{i}{j}" for i in range(3) for j in range(i, 3)
    ]

    arg_dict_cuda = {
        **{var: "const REAL" for var in list_of_metricvarnames},
        **{var.replace("*", ""): "REAL *" for var in list_of_vrnms},
        **{f"xx{i}": "const REAL" for i in range(3)},
    }

    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }

    # Sort the lhss list alphabetically, and rhss to match:
    lhss, rhss = [
        list(x)
        for x in zip(
            *sorted(zip(list_of_vrnms, list_of_exprs), key=lambda pair: pair[0])
        )
    ]

    # Find symbols stored in params
    param_symbols, _ = get_params_commondata_symbols_from_expr_list(
        rhss, exclude=[f"xx{j}" for j in range(3)]
    )
    loop_params = ""

    params_definitions = generate_definition_header(
        param_symbols,
        var_access=parallel_utils.get_params_access(parallelization),
    )
    kernel_body = f"{loop_params}\n{params_definitions}\n"
    kernel_body += "  // Compute tetrads:\n"

    kernel_body += ccg.c_codegen(
        rhss,
        lhss,
        verbose=False,
        enable_cse=True,
        include_braces=False,
    )

    prefunc, prefunc_launch = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=desc,
        cfunc_type=cfunc_type,
        launchblock_with_braces=False,
        cfunc_decorators=cfunc_decorators,
    )

    for var in list_of_metricvarnames:
        prefunc_launch = prefunc_launch.replace(
            var, f"in_gfs[IDX4pt({var.upper()}GF, idx3)]"
        )
    for pvar in list_of_vrnms:
        var = pvar.replace("*", "")
        prefunc_launch = prefunc_launch.replace(var, rf"&{var}")
    return prefunc, prefunc_launch
