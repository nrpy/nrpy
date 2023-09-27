"""
Library of C functions for solving the BSSN equations in
 ***curvilinear*** coordinates, using a reference-metric
 formalism

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import List, Union, cast
from pathlib import Path
from inspect import currentframe as cfr
from types import FrameType as FT
import sympy as sp

import nrpy.grid as gri
import nrpy.params as par
import nrpy.c_codegen as ccg
import nrpy.finite_difference as fin
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric

import nrpy.helpers.parallel_codegen as pcg

from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs
from nrpy.equations.general_relativity.BSSN_gauge_RHSs import BSSN_gauge_RHSs
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.equations.general_relativity.InitialData_Cartesian import (
    InitialData_Cartesian,
)
import nrpy.equations.general_relativity.psi4 as psifour
import nrpy.equations.general_relativity.psi4_tetrads as psifourtet
import nrpy.infrastructures.BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter as admid
import nrpy.infrastructures.BHaH.simple_loop as lp


def register_CFunction_initial_data(
    CoordSystem: str, IDtype: str, IDCoordSystem: str
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C functions for converting ADM initial data to BSSN variables and applying boundary conditions.
    The function performs the following operations:
    1. Registers the exact ADM initial data function.
    2. Registers a function for converting ADM initial data to BSSN variables in the specified coordinate system.
    3. Generates C code for setting initial data and applying boundary conditions.

    :param CoordSystem: The coordinate system for the calculation
    :param IDtype: The type of initial data
    :param IDCoordSystem: The native coordinate system of the initial data
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    ID = InitialData_Cartesian(IDtype=IDtype)

    admid.register_CFunction_exact_ADM_ID_function(
        IDCoordSystem, IDtype, ID.alpha, ID.betaU, ID.BU, ID.gammaDD, ID.KDD
    )
    admid.register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(
        CoordSystem, IDCoordSystem=IDCoordSystem
    )

    desc = "Set initial data."
    c_type = "void"
    name = "initial_data"
    params = "const commondata_struct *restrict commondata, griddata_struct *restrict griddata"

    # Unpack griddata
    body = """
ID_persist_struct ID_persist;
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
"""
    body += f"initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN(commondata, griddata, &ID_persist, {IDtype});"
    body += """
  apply_bcs_outerextrap_and_inner(commondata, params, &griddata->bcstruct, griddata->gridfuncs.y_n_gfs);
}
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_rhs_eval(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    enable_fd_functions: bool,
    LapseEvolutionOption: str,
    ShiftEvolutionOption: str,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side evaluation function for the BSSN equations.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD (Single Instruction, Multiple Data).
    :param LapseEvolutionOption: Lapse evolution equation choice.
    :param ShiftEvolutionOption: Lapse evolution equation choice.
    :return: None
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    if enable_simd:
        includes += [str(Path("simd") / "simd_intrinsics.h")]
    desc = r"""Set RHSs for the BSSN evolution equations."""
    c_type = "void"
    name = "rhs_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs"
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    # Populate BSSN rhs variables
    rhs = BSSN_RHSs[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    BSSN_RHSs_varnames = rhs.BSSN_RHSs_varnames
    BSSN_RHSs_exprs = rhs.BSSN_RHSs_exprs
    alpha_rhs, vet_rhsU, bet_rhsU = BSSN_gauge_RHSs(
        CoordSystem,
        enable_rfm_precompute,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    BSSN_RHSs_varnames += ["alpha_rhs"]
    BSSN_RHSs_exprs += [alpha_rhs]
    for i in range(3):
        BSSN_RHSs_varnames += [f"vet_rhsU{i}", f"bet_rhsU{i}"]
        BSSN_RHSs_exprs += [vet_rhsU[i], bet_rhsU[i]]
    sorted_list = sorted(zip(BSSN_RHSs_varnames, BSSN_RHSs_exprs))
    BSSN_RHSs_varnames, BSSN_RHSs_exprs = [list(t) for t in zip(*sorted_list)]

    BSSN_RHSs_access_gf: List[str] = []
    for var in BSSN_RHSs_varnames:
        BSSN_RHSs_access_gf += [
            gri.BHaHGridFunction.access_gf(
                var.replace("_rhs", ""),
                0,
                0,
                0,
                gf_array_name="rhs_gfs",
            )
        ]
    # Set up upwind control vector (betaU)
    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    betaU = ixp.zerorank1()
    vetU = ixp.declarerank1("vetU")
    for i in range(3):
        # self.lambda_rhsU[i] = self.Lambdabar_rhsU[i] / rfm.ReU[i]
        betaU[i] = vetU[i] * rfm.ReU[i]
    body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            BSSN_RHSs_exprs,
            BSSN_RHSs_access_gf,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            upwind_control_vec=betaU,
            enable_fd_functions=enable_fd_functions,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        prefunc=fin.construct_FD_functions_prefunc() if enable_fd_functions else "",
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        enable_simd=enable_simd,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_Ricci_eval(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    enable_fd_functions: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the Ricci evaluation function.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD instructions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :return: None
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    orig_enable_RbarDD_gridfunctions = par.parval_from_str(
        "enable_RbarDD_gridfunctions"
    )
    if orig_enable_RbarDD_gridfunctions:
        # try/except in case BSSN_quantities hasn't been set yet.
        try:
            del BSSN_quantities[
                CoordSystem + "_rfm_precompute"
                if enable_rfm_precompute
                else CoordSystem
            ]
        except KeyError:
            pass
        par.set_parval_from_str("enable_RbarDD_gridfunctions", False)
    Bq = BSSN_quantities[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    includes = ["BHaH_defines.h"]
    if enable_simd:
        includes += [str(Path("simd") / "simd_intrinsics.h")]
    desc = r"""Set Ricci tensor."""
    c_type = "void"
    name = "Ricci_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict in_gfs, REAL *restrict auxevol_gfs"
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    # Populate Ricci tensor
    Ricci_access_gfs: List[str] = []
    for var in Bq.Ricci_varnames:
        Ricci_access_gfs += [
            gri.BHaHGridFunction.access_gf(var, 0, 0, 0, gf_array_name="auxevol_gfs")
        ]
    body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            Bq.Ricci_exprs,
            Ricci_access_gfs,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            enable_fd_functions=enable_fd_functions,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    if orig_enable_RbarDD_gridfunctions:
        par.set_parval_from_str(
            "enable_RbarDD_gridfunctions", orig_enable_RbarDD_gridfunctions
        )
        del BSSN_quantities[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        _ = BSSN_quantities[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        prefunc=fin.construct_FD_functions_prefunc() if enable_fd_functions else "",
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        enable_simd=enable_simd,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_constraints(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    enable_fd_functions: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the BSSN constraints evaluation function.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD instructions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :return: None
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    Bcon = BSSN_constraints[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    includes = ["BHaH_defines.h"]
    if enable_simd:
        includes += [str(Path("simd") / "simd_intrinsics.h")]
    desc = r"""Evaluate BSSN constraints."""
    c_type = "void"
    name = "constraints_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict in_gfs, const REAL *restrict auxevol_gfs, REAL *restrict diagnostic_output_gfs"
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    Constraints_access_gfs: List[str] = []
    for var in ["H", "MSQUARED"]:
        Constraints_access_gfs += [
            gri.BHaHGridFunction.access_gf(
                var, 0, 0, 0, gf_array_name="diagnostic_output_gfs"
            )
        ]
    body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            [Bcon.H, Bcon.Msquared],
            Constraints_access_gfs,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            enable_fd_functions=enable_fd_functions,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    cfc.register_CFunction(
        includes=includes,
        prefunc=fin.construct_FD_functions_prefunc() if enable_fd_functions else "",
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
        enable_simd=enable_simd,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_enforce_detgammabar_equals_detgammahat(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_fd_functions: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the function that enforces the det(gammabar) = det(gammahat) constraint.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :return: None
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    Bq = BSSN_quantities[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    includes = ["BHaH_defines.h"]
    desc = r"""Enforce det(gammabar) = det(gammahat) constraint. Required for strong hyperbolicity."""
    c_type = "void"
    name = "enforce_detgammabar_equals_detgammahat"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs"
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )

    # First define the Kronecker delta:
    KroneckerDeltaDD = ixp.zerorank2()
    for i in range(3):
        KroneckerDeltaDD[i][i] = sp.sympify(1)

    # The detgammabar in BSSN_RHSs is set to detgammahat when BSSN_RHSs::detgbarOverdetghat_equals_one=True (default),
    #    so we manually compute it here:
    dummygammabarUU, detgammabar = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)

    # Next apply the constraint enforcement equation above.
    nrpyAbs = sp.Function("nrpyAbs")
    hprimeDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            hprimeDD[i][j] = (nrpyAbs(rfm.detgammahat) / detgammabar) ** (
                sp.Rational(1, 3)
            ) * (KroneckerDeltaDD[i][j] + Bq.hDD[i][j]) - KroneckerDeltaDD[i][j]

    hDD_access_gfs: List[str] = []
    hprimeDD_expr_list: List[sp.Expr] = []
    for i in range(3):
        for j in range(i, 3):
            hDD_access_gfs += [
                gri.BHaHGridFunction.access_gf(
                    f"hDD{i}{j}", 0, 0, 0, gf_array_name="in_gfs"
                )
            ]
            hprimeDD_expr_list += [hprimeDD[i][j]]

    # To evaluate the cube root, SIMD support requires e.g., SLEEF.
    #   Also need to be careful to not access memory out of bounds!
    #   After all this is a loop over ALL POINTS.
    #   Exercise to the reader: prove that for any reasonable grid,
    #   SIMD loops over grid interiors never write data out of bounds
    #   and are threadsafe for any reasonable number of threads.
    body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            hprimeDD_expr_list,
            hDD_access_gfs,
            enable_fd_codegen=True,
            enable_simd=False,
            enable_fd_functions=enable_fd_functions,
        ),
        loop_region="all points",
        enable_simd=False,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        prefunc=fin.construct_FD_functions_prefunc() if enable_fd_functions else "",
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        enable_simd=False,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


#
def register_CFunction_psi4_part(
    CoordSystem: str,
    which_part: int,
    enable_fd_functions: bool,
    OMP_collapse: int,
    output_empty_function: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Add psi4 to Cfunction dictionary. psi4 is a really huge expression, so we output
    it in 3 parts: psi4_part0, psi4_part1, and psi4_part2

    :param CoordSystem: Coordinate system to be used.
    :param enable_rfm_precompute: Flag to enable or disable the precomputation of reference metric.
    :param which_part: Specifies which part of psi4 to compute.
    :param OMP_collapse: OpenMP collapse clause integer value.
    :param output_empty_function: If True, psi4 will be set to zero.

    :return: Returns the modified NRPy environment, if applicable.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Set up the C function for psi4
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = f"Compute psi4 at all interior gridpoints, part {which_part}"
    name = f"psi4_part{which_part}"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict in_gfs, REAL *restrict diagnostic_output_gfs"

    gri.register_gridfunctions(
        [f"psi4_part{which_part}re", f"psi4_part{which_part}im"], group="AUX"
    )

    if output_empty_function:
        body = " "
        enable_fd_functions = False
    else:
        psi4 = psifour.Psi4(CoordSystem=CoordSystem, enable_rfm_precompute=False)
        body = r"""if(! (Cart_originx == 0 &&Cart_originy == 0 &&Cart_originz == 0) ) {
  fprintf(stderr, "Error: psi4 output assumes that the grid is centered on the origin.\n");
  fprintf(stderr, "       Good news: check out the C code for proposed modifications.\n");
  exit(1);
}
"""
        loop_prefix = """
REAL mre4U0,mre4U1,mre4U2,mre4U3,mim4U0,mim4U1,mim4U2,mim4U3,n4U0,n4U1,n4U2,n4U3;
REAL xx0,xx1,xx2;
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
{
"""
        rfm = refmetric.reference_metric[CoordSystem]
        found_xx = [False] * 3
        for i in range(3):
            for expr in [psi4.psi4_re_pt[which_part], psi4.psi4_im_pt[which_part]]:
                if rfm.xx[i] in expr.free_symbols:
                    found_xx[i] = True
                    break
        for i in range(3):
            if found_xx[i]:
                loop_prefix += f"    xx{i} = xx[{i}][i{i}];\n"
        loop_prefix += r"""
    /* PROPOSED MODIFICATIONS FOR COMPUTING PSI4 ON GRIDS NOT CENTERED ON THE ORIGIN
        REAL xCart_rel_to_globalgrid_center[3];
        xx_to_Cart(commondata, params, xx, i0, i1, i2,  xCart_rel_to_globalgrid_center);
        int ignore_Cart_to_i0i1i2[3];  REAL xx_rel_to_globalgridorigin[3];
        Cart_to_xx_and_nearest_i0i1i2_global_grid_center(commondata, params, xCart_rel_to_globalgrid_center,xx_rel_to_globalgridorigin,ignore_Cart_to_i0i1i2);
        xx0=xx_rel_to_globalgridorigin[0];
        xx1=xx_rel_to_globalgridorigin[1];
        xx2=xx_rel_to_globalgridorigin[2];
    */
}
"""
        body += lp.simple_loop(
            loop_body=loop_prefix
            + ccg.c_codegen(
                [psi4.psi4_re_pt[which_part], psi4.psi4_im_pt[which_part]],
                [
                    gri.BHaHGridFunction.access_gf(
                        f"psi4_part{which_part}re",
                        gf_array_name="diagnostic_output_gfs",
                    ),
                    gri.BHaHGridFunction.access_gf(
                        f"psi4_part{which_part}im",
                        gf_array_name="diagnostic_output_gfs",
                    ),
                ],
                enable_fd_codegen=True,
                enable_fd_functions=enable_fd_functions,
            ),
            loop_region="interior",
            enable_simd=False,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=False,
            read_xxs=False,
            OMP_collapse=OMP_collapse,
        )

    cfc.register_CFunction(
        includes=includes,
        prefunc=fin.construct_FD_functions_prefunc() if enable_fd_functions else "",
        desc=desc,
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
    output_empty_function: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for psi4 tetrad computations.

    :param CoordSystem: The coordinate system to be used.
    :param tetrad: The type of tetrad. Defaults to "quasiKinnersley".
    :param use_metric_to_construct_unit_normal: Whether to use the metric to construct the unit normal. Defaults to False.
    :param output_empty_function: If True, output an empty function body. Defaults to False.
    :return: A NRPy Environment object or None.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"Compute {tetrad} tetrad for psi4, with use_metric_to_construct_unit_normal={use_metric_to_construct_unit_normal}"
    name = "psi4_tetrad"

    # Initialize psi4 tetrad
    psi4tet = psifourtet.Psi4Tetrads(
        CoordSystem,
        enable_rfm_precompute=False,
        tetrad=tetrad,
        use_metric_to_construct_unit_normal=use_metric_to_construct_unit_normal,
    )
    list_of_vrnms = []
    list_of_exprs = []

    for i in range(4):
        list_of_vrnms.append("*mre4U" + str(i))
        list_of_exprs.append(psi4tet.mre4U[i])
    for i in range(4):
        list_of_vrnms.append("*mim4U" + str(i))
        list_of_exprs.append(psi4tet.mim4U[i])
    for i in range(4):
        list_of_vrnms.append("*n4U" + str(i))
        list_of_exprs.append(psi4tet.n4U[i])

    # for i in range(4):
    #     list_of_varnames.extend([f"*mre4U{i}", f"*mim4U{i}", f"*n4U{i}"])
    #     list_of_symbvars.extend([psi4tet.mre4U[i], psi4tet.mim4U[i], psi4tet.n4U[i]])
    #     print(i, sp.count_ops(psi4tet.mre4U[i]))
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params,"
    list_of_metricvarnames = ["cf"] + [
        f"hDD{i}{j}" for i in range(3) for j in range(i, 3)
    ]
    params += ", ".join([f"const REAL {var}" for var in list_of_metricvarnames])
    params += ", " + ", ".join([f"REAL {var}" for var in list_of_vrnms])
    params += ", REAL *restrict xx[3], const int i0, const int i1, const int i2"

    body = ""

    if output_empty_function:
        body += " "
    else:
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
            rhss, lhss, verbose=False, enable_cse=True, include_braces=False
        )

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
