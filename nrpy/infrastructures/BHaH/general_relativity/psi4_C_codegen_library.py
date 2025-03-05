"""
Library of C functions for solving the BSSN equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp
from nrpy.equations.general_relativity import psi4, psi4_tetrads


def register_CFunction_psi4(
    CoordSystem: str,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Add psi4 to Cfunction dictionary.

    :param CoordSystem: Coordinate system to be used.
    :param OMP_collapse: OpenMP collapse clause integer value.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Set up the C function for psi4
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = "Compute psi4 at all interior gridpoints"
    name = "psi4"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict in_gfs, REAL *restrict diagnostic_output_gfs"

    gri.register_gridfunctions(["psi4_re", "psi4_im"], group="AUX")

    psi4_class = psi4.Psi4(CoordSystem=CoordSystem, enable_rfm_precompute=False)
    body = r"""if(! (Cart_originx == 0 &&Cart_originy == 0 &&Cart_originz == 0) ) {
  fprintf(stderr, "Error: psi4 output assumes that the grid is centered on the origin.\n");
  fprintf(stderr, "       Good news: check out the C code for proposed modifications.\n");
  exit(1);
}
"""

    loop_prefix = rf"""
REAL xx0, xx1, xx2;
{{
  xx0 = xx[0][i0];
  xx1 = xx[1][i1];
  xx2 = xx[2][i2];

/* PROPOSED MODIFICATIONS FOR COMPUTING PSI4 ON GRIDS NOT CENTERED ON THE ORIGIN
    REAL xCart_rel_to_globalgrid_center[3];
    xx_to_Cart(commondata, params, xx, i0, i1, i2,  xCart_rel_to_globalgrid_center);
    int ignore_Cart_to_i0i1i2[3];  REAL xx_rel_to_globalgridorigin[3];
    Cart_to_xx_and_nearest_i0i1i2_global_grid_center(commondata, params, xCart_rel_to_globalgrid_center,xx_rel_to_globalgridorigin,ignore_Cart_to_i0i1i2);
    xx0=xx_rel_to_globalgridorigin[0];
    xx1=xx_rel_to_globalgridorigin[1];
    xx2=xx_rel_to_globalgridorigin[2];
*/
}}

// Metric derivative quantities.
REAL arr_gammaDDdDD[3*3*3*3], arr_GammaUDD[3*3*3], arr_KDDdD[3*3*3];
// Tetrad quantities.
REAL mre4U0,mre4U1,mre4U2,mre4U3,mim4U0,mim4U1,mim4U2,mim4U3,n4U0,n4U1,n4U2,n4U3;
const int idx3 = IDX3(i0, i1, i2);

psi4_tetrad(commondata, params,
    in_gfs[IDX4pt(CFGF, idx3)],
    in_gfs[IDX4pt(HDD00GF, idx3)],
    in_gfs[IDX4pt(HDD01GF, idx3)],
    in_gfs[IDX4pt(HDD02GF, idx3)],
    in_gfs[IDX4pt(HDD11GF, idx3)],
    in_gfs[IDX4pt(HDD12GF, idx3)],
    in_gfs[IDX4pt(HDD22GF, idx3)],
    &mre4U0,&mre4U1,&mre4U2,&mre4U3,&mim4U0,&mim4U1,&mim4U2,&mim4U3,&n4U0,&n4U1,&n4U2,&n4U3,
    xx, i0,i1,i2);

psi4_metric_deriv_quantities(commondata, params, in_gfs, xx0, xx1, xx2, i0, i1, i2, arr_gammaDDdDD, arr_GammaUDD, arr_KDDdD);
// Next, unpack gammaDDdDD, GammaUDD, KDDdD from their arrays:
MAYBE_UNUSED REAL {psi4_class.metric_deriv_var_list_str};
{psi4_class.metric_deriv_unpack_arrays}
"""
    body += lp.simple_loop(
        loop_body=loop_prefix
        + ccg.c_codegen(
            [psi4_class.psi4_re, psi4_class.psi4_im],
            [
                gri.BHaHGridFunction.access_gf(
                    "psi4_re",
                    gf_array_name="diagnostic_output_gfs",
                ),
                gri.BHaHGridFunction.access_gf(
                    "psi4_im",
                    gf_array_name="diagnostic_output_gfs",
                ),
            ],
            enable_fd_codegen=True,
        ),
        loop_region="interior",
        enable_intrinsics=False,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=False,
        read_xxs=False,
        OMP_collapse=OMP_collapse,
    )

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_psi4_metric_deriv_quantities(
    CoordSystem: str,
    enable_fd_functions: bool,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for psi4 metric derivative quantity computations.

    :param CoordSystem: The coordinate system to be used.
    :param enable_fd_functions: Flag to enable or disable the finite difference functions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Initialize psi4 tetrad
    psi4_class = psi4.Psi4(
        CoordSystem,
        enable_rfm_precompute=False,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Compute metric derivative quantities gamma_{ij,kl}, Gamma^i_{jk}, and K_{ij,k} needed for psi4."
    name = "psi4_metric_deriv_quantities"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL *restrict in_gfs, const REAL xx0, const REAL xx1, const REAL xx2, const int i0, const int i1, const int i2, REAL arr_gammaDDdDD[81], REAL arr_GammaUDD[27], REAL arr_KDDdD[27]"""

    body = ccg.c_codegen(
        psi4_class.metric_derivs_expr_list,
        psi4_class.metric_derivs_varname_arr_list,
        verbose=False,
        enable_cse=True,
        include_braces=False,
        enable_fd_codegen=True,
        enable_fd_functions=enable_fd_functions,
    )

    cfc.register_CFunction(
        includes=includes,
        prefunc=fin.construct_FD_functions_prefunc() if enable_fd_functions else "",
        desc=desc,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_psi4_tetrad(
    CoordSystem: str,
    tetrad: str = "quasiKinnersley",
    use_metric_to_construct_unit_normal: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for psi4 tetrad computations.

    :param CoordSystem: The coordinate system to be used.
    :param tetrad: The type of tetrad. Defaults to "quasiKinnersley".
    :param use_metric_to_construct_unit_normal: Whether to use the metric to construct the unit normal. Defaults to False.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"Compute {tetrad} tetrad for psi4, with use_metric_to_construct_unit_normal={use_metric_to_construct_unit_normal}"
    name = "psi4_tetrad"

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

    params = "const commondata_struct *restrict commondata, const params_struct *restrict params,"
    list_of_metricvarnames = ["cf"] + [
        f"hDD{i}{j}" for i in range(3) for j in range(i, 3)
    ]
    params += ", ".join([f"const REAL {var}" for var in list_of_metricvarnames])
    params += ", " + ", ".join([f"REAL {var}" for var in list_of_vrnms])
    params += ", REAL *restrict xx[3], const int i0, const int i1, const int i2"

    body = ""

    for i in range(3):
        body += f"  const REAL xx{i} = xx[{i}][i{i}];\n"
    body += "  // Compute tetrads:\n"

    # Sort the lhss list alphabetically, and rhss to match:
    lhss, rhss = [
        list(x)
        for x in zip(
            *sorted(zip(list_of_vrnms, list_of_exprs), key=lambda pair: pair[0])
        )
    ]
    body += ccg.c_codegen(
        rhss,
        lhss,
        verbose=False,
        enable_cse=True,
        include_braces=False,
    )

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
