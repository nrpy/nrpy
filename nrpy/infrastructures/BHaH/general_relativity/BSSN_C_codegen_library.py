"""
Library of C functions for solving the BSSN equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from collections import OrderedDict as ODict
from typing import List, Union, cast, Tuple, Dict
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
from nrpy.equations.general_relativity.InitialData_Spherical import (
    InitialData_Spherical,
)
import nrpy.equations.general_relativity.psi4 as psifour
import nrpy.equations.general_relativity.psi4_tetrads as psifourtet
import nrpy.infrastructures.BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter as admid
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.infrastructures.BHaH.diagnostics.output_0d_1d_2d_slices as out012d


def register_CFunction_initial_data(
    CoordSystem: str,
    IDtype: str,
    IDCoordSystem: str,
    ID_persist_struct_str: str,
    enable_checkpointing: bool = False,
    populate_ID_persist_struct_str: str = "",
    free_ID_persist_struct_str: str = "",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C functions for converting ADM initial data to BSSN variables and applying boundary conditions.

    The function performs the following operations:
    1. Registers the exact ADM initial data function.
    2. Registers a function for converting ADM initial data to BSSN variables in the specified coordinate system.
    3. Generates C code for setting initial data and applying boundary conditions.

    :param CoordSystem: The coordinate system for the calculation.
    :param IDtype: The type of initial data.
    :param IDCoordSystem: The native coordinate system of the initial data.
    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial data.
    :param ID_persist_struct_str: A string representing the persistent structure for the initial data.
    :param populate_ID_persist_struct_str: Optional string to populate the persistent structure for initial data.
    :param free_ID_persist_struct_str: Optional string to free the persistent structure for initial data.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    try:
        ID: Union[InitialData_Cartesian, InitialData_Spherical]
        if IDCoordSystem == "Cartesian":
            ID = InitialData_Cartesian(IDtype=IDtype)
        else:
            ID = InitialData_Spherical(IDtype=IDtype)

        admid.register_CFunction_exact_ADM_ID_function(
            IDCoordSystem, IDtype, ID.alpha, ID.betaU, ID.BU, ID.gammaDD, ID.KDD
        )
    except (ValueError, RuntimeError):
        print(
            f"Warning: {IDtype} does not correspond to an implemented exact initial data type."
        )
        print("Assuming initial data functionality is implemented elsewhere.")

    admid.register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(
        CoordSystem,
        IDCoordSystem=IDCoordSystem,
        ID_persist_struct_str=ID_persist_struct_str,
    )

    desc = "Set initial data."
    c_type = "void"
    name = "initial_data"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    body = ""
    if enable_checkpointing:
        body += """// Attempt to read checkpoint file. If it doesn't exist, then continue. Otherwise return.
if( read_checkpoint(commondata, griddata) ) return;
"""
    body += "ID_persist_struct ID_persist;\n"
    if populate_ID_persist_struct_str:
        body += populate_ID_persist_struct_str
    body += """
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
"""
    body += f"""initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN(commondata, params,
griddata[grid].xx, &griddata[grid].bcstruct, &griddata[grid].gridfuncs, &ID_persist, {IDtype});"""
    body += """
  apply_bcs_outerextrap_and_inner(commondata, params, &griddata[grid].bcstruct, griddata[grid].gridfuncs.y_n_gfs);
}
"""
    if free_ID_persist_struct_str:
        body += free_ID_persist_struct_str

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


def register_CFunction_diagnostics(
    list_of_CoordSystems: List[str],
    default_diagnostics_out_every: float,
    enable_psi4_diagnostics: bool = False,
    enable_progress_indicator: bool = True,
    grid_center_filename_tuple: Tuple[str, str] = (
        "out0d-conv_factor%.2f.txt",
        "convergence_factor",
    ),
    axis_filename_tuple: Tuple[str, str] = (
        "out1d-AXIS-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
    plane_filename_tuple: Tuple[str, str] = (
        "out2d-PLANE-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
    out_quantities_dict: Union[str, Dict[Tuple[str, str], str]] = "default",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for simulation diagnostics.

    :param list_of_CoordSystems: Lists of unique CoordSystems used.
    :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
    :param enable_psi4_diagnostics: Whether or not to enable psi4 diagnostics.
    :param enable_progress_indicator: Whether or not to enable the progress indicator.
    :param grid_center_filename_tuple: Tuple containing filename and variables for grid center output.
    :param axis_filename_tuple: Tuple containing filename and variables for axis output.
    :param plane_filename_tuple: Tuple containing filename and variables for plane output.
    :param out_quantities_dict: Dictionary or string specifying output quantities.

    :return: None if in registration phase, else the updated NRPy environment.
    :raises TypeError: If `out_quantities_dict` is not a dictionary and not set to "default".
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    _ = par.CodeParameter(
        "REAL",
        __name__,
        "diagnostics_output_every",
        default_diagnostics_out_every,
        commondata=True,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    # fmt: off
    if out_quantities_dict == "default":
        out_quantities_dict = {
            ("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))",
            ("REAL", "log10sqrtM2L"): "log10(sqrt(diagnostic_output_gfs[IDX4pt(MSQUAREDGF, idx3)]) + 1e-16)",
            ("REAL", "cfL"): "y_n_gfs[IDX4pt(CFGF, idx3)]",
            ("REAL", "alphaL"): "y_n_gfs[IDX4pt(ALPHAGF, idx3)]",
            ("REAL", "trKL"): "y_n_gfs[IDX4pt(TRKGF, idx3)]",
        }
    if not isinstance(out_quantities_dict, dict):
        raise TypeError(f"out_quantities_dict was initialized to {out_quantities_dict}, which is not a dictionary!")
    # fmt: on

    for CoordSystem in list_of_CoordSystems:
        out012d.register_CFunction_diagnostics_grid_center(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            filename_tuple=grid_center_filename_tuple,
        )
        for axis in ["y", "z"]:
            out012d.register_CFunction_diagnostics_1d_axis(
                CoordSystem=CoordSystem,
                out_quantities_dict=out_quantities_dict,
                filename_tuple=axis_filename_tuple,
                axis=axis,
            )
        for plane in ["xy", "yz"]:
            out012d.register_CFunction_diagnostics_2d_plane(
                CoordSystem=CoordSystem,
                out_quantities_dict=out_quantities_dict,
                filename_tuple=plane_filename_tuple,
                plane=plane,
            )

    desc = r"""Diagnostics."""
    c_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    body = r"""
const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->diagnostics_output_every;
// Explanation of the if() below:
// Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime.
// Step 2: Multiplying by outevery yields the exact time we should output again, t_out.
// Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible!
if(fabs(round(currtime / outevery) * outevery - currtime) < 0.5*currdt) {
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    // Unpack griddata struct:
    const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;
    REAL *restrict xx[3];
    {
      for (int ww = 0; ww < 3; ww++)
        xx[ww] = griddata[grid].xx[ww];
    }
    const params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

    // Constraint output
    {
      Ricci_eval(commondata, params, &griddata[grid].rfmstruct, y_n_gfs, auxevol_gfs);
      constraints_eval(commondata, params, &griddata[grid].rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_output_gfs);
    }

    // 0D output
    diagnostics_grid_center(commondata, params, &griddata[grid].gridfuncs);

    // 1D output
    diagnostics_1d_y_axis(commondata, params, xx, &griddata[grid].gridfuncs);
    diagnostics_1d_z_axis(commondata, params, xx, &griddata[grid].gridfuncs);

    // 2D output
    diagnostics_2d_xy_plane(commondata, params, xx, &griddata[grid].gridfuncs);
    diagnostics_2d_yz_plane(commondata, params, xx, &griddata[grid].gridfuncs);
"""
    if enable_psi4_diagnostics:
        body += r"""      // Do psi4 output, but only if the grid is spherical-like.
      if (strstr(CoordSystemName, "Spherical") != NULL) {

        // Adjusted to match Tutorial-Start_to_Finish-BSSNCurvilinear-Two_BHs_Collide-Psi4.ipynb
        const int psi4_spinweightm2_sph_harmonics_max_l = 2;
#define num_of_R_exts 24
        const REAL list_of_R_exts[num_of_R_exts] =
        { 10.0, 20.0, 21.0, 22.0, 23.0,
          24.0, 25.0, 26.0, 27.0, 28.0,
          29.0, 30.0, 31.0, 32.0, 33.0,
          35.0, 40.0, 50.0, 60.0, 70.0,
          80.0, 90.0, 100.0, 150.0 };

        // Set psi4.
        psi4_part0(commondata, params, xx, y_n_gfs, diagnostic_output_gfs);
        psi4_part1(commondata, params, xx, y_n_gfs, diagnostic_output_gfs);
        psi4_part2(commondata, params, xx, y_n_gfs, diagnostic_output_gfs);
        // Decompose psi4 into spin-weight -2  spherical harmonics & output to files.
        psi4_spinweightm2_decomposition_on_sphlike_grids(commondata, params, diagnostic_output_gfs, list_of_R_exts, num_of_R_exts, psi4_spinweightm2_sph_harmonics_max_l, xx);
      }
"""
    body += r"""
  }
}
"""
    if enable_progress_indicator:
        body += "progress_indicator(commondata, griddata);"
    body += r"""
if(commondata->time + commondata->dt > commondata->t_final) printf("\n");
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
    enable_KreissOliger_dissipation: bool,
    LapseEvolutionOption: str,
    ShiftEvolutionOption: str,
    KreissOliger_strength_mult_by_W: bool = False,
    # when mult by W, strength_gauge=0.99 & strength_nongauge=0.3 is best.
    KreissOliger_strength_gauge: float = 0.3,
    KreissOliger_strength_nongauge: float = 0.3,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side evaluation function for the BSSN equations.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD (Single Instruction, Multiple Data).
    :param enable_fd_functions: Whether or not to enable finite difference functions.
    :param enable_KreissOliger_dissipation: Whether or not to enable Kreiss-Oliger dissipation.
    :param LapseEvolutionOption: Lapse evolution equation choice.
    :param ShiftEvolutionOption: Lapse evolution equation choice.
    :param KreissOliger_strength_mult_by_W: Whether to multiply Kreiss-Oliger strength by W.
    :param KreissOliger_strength_gauge: Gauge strength for Kreiss-Oliger dissipation.
    :param KreissOliger_strength_nongauge: Non-gauge strength for Kreiss-Oliger dissipation.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
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
    alpha_rhs, vet_rhsU, bet_rhsU = BSSN_gauge_RHSs(
        CoordSystem,
        enable_rfm_precompute,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    rhs.BSSN_RHSs_varname_to_expr_dict["alpha_rhs"] = alpha_rhs
    for i in range(3):
        rhs.BSSN_RHSs_varname_to_expr_dict[f"vet_rhsU{i}"] = vet_rhsU[i]
        rhs.BSSN_RHSs_varname_to_expr_dict[f"bet_rhsU{i}"] = bet_rhsU[i]

    rhs.BSSN_RHSs_varname_to_expr_dict = ODict(
        sorted(rhs.BSSN_RHSs_varname_to_expr_dict.items())
    )

    # Add Kreiss-Oliger dissipation to the BSSN RHSs:
    if enable_KreissOliger_dissipation:
        diss_strength_gauge, diss_strength_nongauge = par.register_CodeParameters(
            "REAL",
            __name__,
            ["KreissOliger_strength_gauge", "KreissOliger_strength_nongauge"],
            [KreissOliger_strength_gauge, KreissOliger_strength_nongauge],
            commondata=True,
        )

        if KreissOliger_strength_mult_by_W:
            Bq = BSSN_quantities[
                CoordSystem + "_rfm_precompute"
                if enable_rfm_precompute
                else CoordSystem
            ]
            EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
            if EvolvedConformalFactor_cf == "W":
                diss_strength_gauge *= Bq.cf
                diss_strength_nongauge *= Bq.cf
            elif EvolvedConformalFactor_cf == "chi":
                diss_strength_gauge *= sp.sqrt(Bq.cf)
                diss_strength_nongauge *= sp.sqrt(Bq.cf)
            elif EvolvedConformalFactor_cf == "phi":
                diss_strength_gauge *= sp.exp(-2 * Bq.cf)
                diss_strength_nongauge *= sp.exp(-2 * Bq.cf)

        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        alpha_dKOD = ixp.declarerank1("alpha_dKOD")
        cf_dKOD = ixp.declarerank1("cf_dKOD")
        trK_dKOD = ixp.declarerank1("trK_dKOD")
        betU_dKOD = ixp.declarerank2("betU_dKOD", symmetry="nosym")
        vetU_dKOD = ixp.declarerank2("vetU_dKOD", symmetry="nosym")
        lambdaU_dKOD = ixp.declarerank2("lambdaU_dKOD", symmetry="nosym")
        aDD_dKOD = ixp.declarerank3("aDD_dKOD", symmetry="sym01")
        hDD_dKOD = ixp.declarerank3("hDD_dKOD", symmetry="sym01")
        for k in range(3):
            rhs.BSSN_RHSs_varname_to_expr_dict["alpha_rhs"] += (
                diss_strength_gauge * alpha_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.BSSN_RHSs_varname_to_expr_dict["cf_rhs"] += (
                diss_strength_nongauge * cf_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.BSSN_RHSs_varname_to_expr_dict["trK_rhs"] += (
                diss_strength_nongauge * trK_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            for i in range(3):
                if "2ndOrder" in ShiftEvolutionOption:
                    rhs.BSSN_RHSs_varname_to_expr_dict[f"bet_rhsU{i}"] += (
                        diss_strength_gauge * betU_dKOD[i][k] * rfm.ReU[k]
                    )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                rhs.BSSN_RHSs_varname_to_expr_dict[f"vet_rhsU{i}"] += (
                    diss_strength_gauge * vetU_dKOD[i][k] * rfm.ReU[k]
                )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                rhs.BSSN_RHSs_varname_to_expr_dict[f"lambda_rhsU{i}"] += (
                    diss_strength_nongauge * lambdaU_dKOD[i][k] * rfm.ReU[k]
                )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                for j in range(i, 3):
                    rhs.BSSN_RHSs_varname_to_expr_dict[f"a_rhsDD{i}{j}"] += (
                        diss_strength_nongauge * aDD_dKOD[i][j][k] * rfm.ReU[k]
                    )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                    rhs.BSSN_RHSs_varname_to_expr_dict[f"h_rhsDD{i}{j}"] += (
                        diss_strength_nongauge * hDD_dKOD[i][j][k] * rfm.ReU[k]
                    )  # ReU[k] = 1/scalefactor_orthog_funcform[k]

    BSSN_RHSs_access_gf: List[str] = []
    for var in rhs.BSSN_RHSs_varname_to_expr_dict.keys():
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
            list(rhs.BSSN_RHSs_varname_to_expr_dict.values()),
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
        CoordSystem_for_wrapper_func=CoordSystem,
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
    :param enable_fd_functions: Whether or not to enable finite difference functions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
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
        CoordSystem_for_wrapper_func=CoordSystem,
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
    :param enable_fd_functions: Whether or not to enable finite difference functions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
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
        CoordSystem_for_wrapper_func=CoordSystem,
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

    :return: None if in registration phase, else the updated NRPy environment.
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
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=body,
        enable_simd=False,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_psi4_part(
    CoordSystem: str,
    which_part: int,
    enable_fd_functions: bool,
    OMP_collapse: int,
    output_empty_function: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Add psi4 to Cfunction dictionary.

    psi4 is a really huge expression, so we output it in 3 parts:
      psi4_part0, psi4_part1, and psi4_part2.

    :param CoordSystem: Coordinate system to be used.
    :param which_part: Specifies which part of psi4 to compute.
    :param enable_fd_functions: Flag to enable or disable the finite difference functions.
    :param OMP_collapse: OpenMP collapse clause integer value.
    :param output_empty_function: If True, psi4 will be set to zero.

    :return: None if in registration phase, else the updated NRPy environment.
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
    output_empty_function: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for psi4 tetrad computations.

    :param CoordSystem: The coordinate system to be used.
    :param tetrad: The type of tetrad. Defaults to "quasiKinnersley".
    :param use_metric_to_construct_unit_normal: Whether to use the metric to construct the unit normal. Defaults to False.
    :param output_empty_function: If True, output an empty function body. Defaults to False.

    :return: None if in registration phase, else the updated NRPy environment.
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
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_psi4_spinweightm2_decomposition_on_sphlike_grids() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register C function for decomposing psi4 into spin-weighted spherical harmonics.

    :param None: No parameters for this function.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    prefunc = r"""
static void lowlevel_decompose_psi4_into_swm2_modes(const int Nxx_plus_2NGHOSTS1,const int Nxx_plus_2NGHOSTS2,
                                                    const REAL dxx1, const REAL dxx2,
                                                    const int swm2sh_maximum_l_mode_to_compute, const REAL curr_time, const REAL R_ext,
                                                    const REAL *restrict th_array, const REAL *restrict sinth_array, const REAL *restrict ph_array,
                                                    const REAL *restrict psi4r_at_R_ext, const REAL *restrict psi4i_at_R_ext) {
  char filename[100];  FILE *outpsi4_l_m;
  // Output header at t=0:
  if(curr_time==0) {
    for(int l=2;l<=swm2sh_maximum_l_mode_to_compute;l++) {
      sprintf(filename,"Rpsi4_l%d-r%06.1f.txt",l,(double)R_ext);
      outpsi4_l_m = fopen(filename, "w");
      fprintf(outpsi4_l_m, "# column 1: t-R_ext = [retarded time]\n");
      int col=2;
      for(int m=-l;m<=l;m++) {
        fprintf(outpsi4_l_m, "# column %d: Re(psi4_{l=%d,m=%d}) * R_ext\n", col,l,m);  col++;
        fprintf(outpsi4_l_m, "# column %d: Im(psi4_{l=%d,m=%d}) * R_ext\n", col,l,m);  col++;
      }
      fclose(outpsi4_l_m);
    }
  }

  // Output one file per l mode; each column represents a unique complex component of l,m
  for(int l=2;l<=swm2sh_maximum_l_mode_to_compute;l++) {
    sprintf(filename,"Rpsi4_l%d-r%06.1f.txt",l,(double)R_ext);
    outpsi4_l_m = fopen(filename, "a");
    char oneline[10000];
    sprintf(oneline, "%e", (double)(curr_time - R_ext));
    for(int m=-l;m<=l;m++) {
      // Parallelize the integration loop:
      REAL psi4r_l_m = 0.0;
      REAL psi4i_l_m = 0.0;
#pragma omp parallel for reduction(+:psi4r_l_m,psi4i_l_m)
      for(int i1=0;i1<Nxx_plus_2NGHOSTS1-2*NGHOSTS;i1++) {
        const REAL th    = th_array[i1];
        const REAL sinth = sinth_array[i1];
        for(int i2=0;i2<Nxx_plus_2NGHOSTS2-2*NGHOSTS;i2++) {
          const REAL ph = ph_array[i2];
          // Construct integrand for psi4 spin-weight s=-2 spherical harmonic
          REAL ReY_sm2_l_m,ImY_sm2_l_m;
          spin_weight_minus2_sph_harmonics(l,m, th,ph,  &ReY_sm2_l_m,&ImY_sm2_l_m);

          const int idx2d = i1*(Nxx_plus_2NGHOSTS2-2*NGHOSTS)+i2;
          const REAL a = psi4r_at_R_ext[idx2d];
          const REAL b = psi4i_at_R_ext[idx2d];
          const REAL c = ReY_sm2_l_m;
          const REAL d = ImY_sm2_l_m;
          psi4r_l_m += (a*c + b*d) * dxx2  * sinth*dxx1;
          psi4i_l_m += (b*c - a*d) * dxx2  * sinth*dxx1;
        }
      }
      sprintf(oneline + strlen(oneline), " %.15e %.15e", (double)(R_ext*psi4r_l_m), (double)(R_ext*psi4i_l_m));
    }
    fprintf(outpsi4_l_m, "%s\n", oneline);
    fclose(outpsi4_l_m);
  }
}
"""

    desc = "Decompose psi4 across all l,m modes from l=2 up to and including L_MAX (global variable)"
    name = "psi4_spinweightm2_decomposition_on_sphlike_grids"
    params = r"""const commondata_struct *restrict commondata,const params_struct *restrict params,
    REAL *restrict diagnostic_output_gfs,
    const REAL *restrict list_of_R_exts, const int num_of_R_exts,
    const int psi4_spinweightm2_sph_harmonics_max_l, REAL *restrict xx[3]"""

    body = r"""  // Step 1: Allocate memory for 2D arrays used to store psi4, theta, sin(theta), and phi.
  const int sizeof_2Darray = sizeof(REAL)*(Nxx_plus_2NGHOSTS1-2*NGHOSTS)*(Nxx_plus_2NGHOSTS2-2*NGHOSTS);
  REAL *restrict psi4r_at_R_ext = (REAL *restrict)malloc(sizeof_2Darray);
  REAL *restrict psi4i_at_R_ext = (REAL *restrict)malloc(sizeof_2Darray);
  //         ... also store theta, sin(theta), and phi to corresponding 1D arrays.
  REAL *restrict sinth_array = (REAL *restrict)malloc(sizeof(REAL)*(Nxx_plus_2NGHOSTS1-2*NGHOSTS));
  REAL *restrict th_array    = (REAL *restrict)malloc(sizeof(REAL)*(Nxx_plus_2NGHOSTS1-2*NGHOSTS));
  REAL *restrict ph_array    = (REAL *restrict)malloc(sizeof(REAL)*(Nxx_plus_2NGHOSTS2-2*NGHOSTS));

  const int NinterpGHOSTS = MIN(2, NGHOSTS-1);
  const int N0 = 2*NinterpGHOSTS; // Interp stencil is 2*NinterpGHOSTS+1 in size;
  //                                 reaches NinterpGHOSTS to the left & right of
  //                                 central point.
  const REAL pow_dxx0__N0 = pow(params->dxx0, N0);

  // Step 2: Loop over all extraction indices:
  for(int which_R_ext=0;which_R_ext<num_of_R_exts;which_R_ext++) {
    // Step 2.a: Set the extraction radius R_ext based on the radial index R_ext_idx
    const REAL R_ext = list_of_R_exts[which_R_ext];
    const REAL xCart_R_ext[3] = { R_ext, 0.0, 0.0 }; // just put a point on the x-axis.

    int Cart_to_i0i1i2[3]; REAL closest_xx[3];
    Cart_to_xx_and_nearest_i0i1i2(commondata,params, xCart_R_ext, closest_xx, Cart_to_i0i1i2);

    const int closest_i0=Cart_to_i0i1i2[0];

    // We want a src grid point inside the source grid (duh) with
    //  mask=+0, as all mask=+0 points will have at least
    //  NGHOSTS>=NinterpGHOSTS of filled neighbor pts.
    if(IS_IN_GRID_INTERIOR(Cart_to_i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NGHOSTS)) {

      // Step 2.a.i: Set radial interpolation coefficients for r=R_ext.
      REAL l0i__times__w0i_inv[2*NinterpGHOSTS+1];
      {
        for(int i=0;i<=N0;i++) {
          REAL prod_numer_i = 1.0;
          int  prod_denom_i = 1;
          for(int l=0;  l<i;  l++) { prod_denom_i *= i-l; prod_numer_i *= closest_xx[0] - xx[0][closest_i0 - N0/2 + l]; }
          for(int l=i+1;l<=N0;l++) { prod_denom_i *= i-l; prod_numer_i *= closest_xx[0] - xx[0][closest_i0 - N0/2 + l]; }
          l0i__times__w0i_inv[i] = prod_numer_i / ( pow_dxx0__N0 * (REAL)prod_denom_i );
        }
      }

      // Step 2.b: Compute psi_4 at this extraction radius and store to a local 2D array.
#pragma omp parallel for
      for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS1-NGHOSTS;i1++) {
        th_array[i1-NGHOSTS]    =     xx[1][i1];
        sinth_array[i1-NGHOSTS] = sin(xx[1][i1]);
        for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS2-NGHOSTS;i2++) {
          ph_array[i2-NGHOSTS] = xx[2][i2];

          REAL sum_psi4r=0;
          REAL sum_psi4i=0;
          // Perform radial interpolation to get psi4 at desired extraction radius R_ext.
          for(int i=0;i<=N0;i++) {
            // psi4r and psi4i in fixed frame have been stored to diagnostic_output_gfs.
            //  Here we interpolate to specific radius.
            sum_psi4r += ( diagnostic_output_gfs[IDX4(PSI4_PART0REGF, i + closest_i0 - N0/2, i1, i2)] +
                           diagnostic_output_gfs[IDX4(PSI4_PART1REGF, i + closest_i0 - N0/2, i1, i2)] +
                           diagnostic_output_gfs[IDX4(PSI4_PART2REGF, i + closest_i0 - N0/2, i1, i2)] ) * l0i__times__w0i_inv[i];
            sum_psi4i += ( diagnostic_output_gfs[IDX4(PSI4_PART0IMGF, i + closest_i0 - N0/2, i1, i2)] +
                           diagnostic_output_gfs[IDX4(PSI4_PART1IMGF, i + closest_i0 - N0/2, i1, i2)] +
                           diagnostic_output_gfs[IDX4(PSI4_PART2IMGF, i + closest_i0 - N0/2, i1, i2)] ) * l0i__times__w0i_inv[i];
          }
          // Store result to "2D" array (actually 1D array with 2D storage):
          const int idx2d = (i1-NGHOSTS)*(Nxx_plus_2NGHOSTS2-2*NGHOSTS)+(i2-NGHOSTS);
          psi4r_at_R_ext[idx2d] = sum_psi4r;
          psi4i_at_R_ext[idx2d] = sum_psi4i;
        }
      }
      // Step 3: Perform integrations across all l,m modes from l=2 up to and including L_MAX (global variable):
      lowlevel_decompose_psi4_into_swm2_modes(Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2,
                                              dxx1,dxx2, swm2sh_maximum_l_mode_to_compute,
                                              time, R_ext, th_array, sinth_array, ph_array,
                                              psi4r_at_R_ext,psi4i_at_R_ext);
    }
  }
  // Step 4: Free all allocated memory:
  free(psi4r_at_R_ext); free(psi4i_at_R_ext);
  free(sinth_array); free(th_array); free(ph_array);
"""

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        prefunc=prefunc,
        desc=desc,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
