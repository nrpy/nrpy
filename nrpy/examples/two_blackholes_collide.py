"""
Sets up a complete C code for solving the general relativistic
  field equations in curvilinear coordinates on a cell-centered
  grid, using a reference metric.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
from typing import List
from pathlib import Path
import shutil
import os
import sympy as sp

import nrpy.params as par
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.c_codegen as ccg
import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric

from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs
from nrpy.equations.general_relativity.BSSN_gauge_RHSs import BSSN_gauge_RHSs
from nrpy.equations.general_relativity.InitialData_Cartesian import (
    InitialData_Cartesian,
)
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints

from nrpy.infrastructures.BHaH.MoLtimestepping import MoL
from nrpy.infrastructures.BHaH import rfm_precompute
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.main_c as main
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions as cbc
import nrpy.infrastructures.BHaH.numerical_grids_and_timestep as numgrids
import nrpy.infrastructures.BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter as admid
import nrpy.infrastructures.BHaH.xx_tofrom_Cart as xxCart

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "two_blackholes_collide"
CoordSystem = "Spherical"
IDtype = "BrillLindquist"
IDCoordSystem = "Cartesian"
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_Covariant"
GammaDriving_eta = 1.0
grid_physical_size = 7.5
t_final = 1.0 * grid_physical_size
Nxx_dict = {
    "Spherical": [72, 12, 2],
    "SinhSpherical": [72, 12, 2],
    "Cartesian": [64, 64, 64],
}
default_BH1_mass = default_BH2_mass = 0.5
default_BH1_z_posn = +0.5
default_BH2_z_posn = -0.5
enable_rfm_precompute = True
MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 2
enable_simd = True
separate_Ricci_and_BSSN_RHS = True

OMP_collapse = 1
if "Spherical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "2")
    OMP_collapse = 2  # about 2x faster
project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("LeaveRicciSymbolic", separate_Ricci_and_BSSN_RHS)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)

# In the following function we overwrite t_final (registered as a CodeParameter in MoL),
#    so let's remove it from the parfile
par.adjust_CodeParam_default("t_final", t_final)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
def register_CFunction_initial_data() -> None:
    """
    Register the initial data function for the wave equation with specific parameters.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    ID = InitialData_Cartesian(IDtype=IDtype)
    admid.register_CFunction_exact_ADM_ID_function(
        IDCoordSystem, IDtype, ID.alpha, ID.betaU, ID.BU, ID.gammaDD, ID.KDD
    )
    admid.register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(
        CoordSystem, IDCoordSystem=IDCoordSystem
    )

    desc = r"""Set initial data."""
    c_type = "void"
    name = "initial_data"
    params = "const commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    # Unpack griddata
    body = r"""
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


def register_CFunction_diagnostics(in_CoordSystem: str, plane: str = "yz") -> None:
    """
    Register C function for simulation diagnostics.

    :param in_CoordSystem: Specifies the coordinate system for the diagnostics.
    :param plane: The default plane for diagnostics; defaults to "yz".
    :return: None
    """
    _ = par.CodeParameter(
        "REAL", __name__, "diagnostics_output_every", 0.25, commondata=True
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

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
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    const params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

    // Constraint output
    {
      Ricci_eval(commondata, params, &griddata->rfmstruct, y_n_gfs, auxevol_gfs);
      constraints_eval(commondata, params, &griddata->rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_output_gfs);
    }
    // 0D output
    {
      char filename[256];
      sprintf(filename, "out0d-conv_factor%.2f.txt", convergence_factor);
      FILE *outfile;
      if (nn == 0)
        outfile = fopen(filename, "w");
      else
        outfile = fopen(filename, "a");
      if (outfile == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
        exit(1);
      }

      const int center_of_grid_idx = IDX3(NGHOSTS, Nxx_plus_2NGHOSTS1 / 2, Nxx_plus_2NGHOSTS2 / 2);

      const REAL H_at_center = diagnostic_output_gfs[IDX4pt(HGF, center_of_grid_idx)];
      const REAL M2_at_center = diagnostic_output_gfs[IDX4pt(MSQUAREDGF, center_of_grid_idx)];
      const REAL alpha_at_center = y_n_gfs[IDX4pt(ALPHAGF, center_of_grid_idx)];
      const REAL trK_at_center = y_n_gfs[IDX4pt(TRKGF, center_of_grid_idx)];

      fprintf(outfile, "%e %e %e %e %e\n", time, log10(fabs(H_at_center + 1e-16)), log10(fabs(M2_at_center + 1e-16)), alpha_at_center, trK_at_center);

      fclose(outfile);
    }

    // 1D output
    {
      char filename[256];
      sprintf(filename, "out1d-conv_factor%.2f-t%.2f.txt", convergence_factor, time);
      FILE *outfile;
      outfile = fopen(filename, "w");
      if (outfile == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
        exit(1);
      }
      for (int i1 = Nxx_plus_2NGHOSTS1 - NGHOSTS - 1; i1 >= NGHOSTS; i1 -= Nxx1-1) {
        int i0_start, i0_end, i0_step;

        if (i1 == (Nxx_plus_2NGHOSTS1 - NGHOSTS - 1)) {
          i0_start = Nxx_plus_2NGHOSTS0 - NGHOSTS - 1;
          i0_end = NGHOSTS - 1;
          i0_step = -1;
        } else if (i1 == NGHOSTS) {
          i0_start = NGHOSTS;
          i0_end = Nxx_plus_2NGHOSTS0 - NGHOSTS;
          i0_step = 1;
        } else {
          continue; // Skip this iteration if i1 is not one of the special values
        }

        for (int i0 = i0_start; i0 != i0_end; i0 += i0_step) {
          const int i2 = Nxx_plus_2NGHOSTS2 / 2;
          const int idx3 = IDX3(i0, i1, i2);
          REAL xCart[3];
          xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

          const REAL alphaL = y_n_gfs[IDX4pt(ALPHAGF, idx3)];
          const REAL trKL = y_n_gfs[IDX4pt(TRKGF, idx3)];
          const REAL HL = diagnostic_output_gfs[IDX4pt(HGF, idx3)];
          const REAL M2L = diagnostic_output_gfs[IDX4pt(MSQUAREDGF, idx3)];

          fprintf(outfile, "%e %e %e %e %e\n", xCart[2], log10(fabs(HL + 1e-16)), log10(fabs(M2L + 1e-16)), alphaL, trKL);
        }
      }
      fclose(outfile);
    }
"""
    body += rf"""
    // 2D output:
    {{
      char filename[256];
      sprintf(filename, "out2d-{plane}_plane-conv_factor%.2f-t%.2f.txt", convergence_factor, time);
      FILE *outfile;
      outfile = fopen(filename, "w");
      if (outfile == NULL) {{
        fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
        exit(1);
      }}
"""
    body += lp.simple_loop_2D(
        loop_body=r"""
const int idx3 = IDX3(i0, i1, i2);
REAL xCart[3];
xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

const REAL alphaL = y_n_gfs[IDX4pt(ALPHAGF, idx3)];
const REAL trKL = y_n_gfs[IDX4pt(TRKGF, idx3)];
const REAL HL = diagnostic_output_gfs[IDX4pt(HGF, idx3)];
const REAL M2L = diagnostic_output_gfs[IDX4pt(MSQUAREDGF, idx3)];

fprintf(outfile, "%e %e %e %e %e %e %e\n", xCart[0], xCart[1], xCart[2], log10(fabs(HL + 1e-16)), log10(fabs(M2L + 1e-16)), alphaL, trKL);
""",
        CoordSystem=in_CoordSystem,
        plane=plane,
    )
    body += r"""
      fclose(outfile);
    }
  }
}
progress_indicator(commondata, griddata);
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


def register_CFunction_rhs_eval(
    in_CoordSystem: str,
    in_enable_rfm_precompute: bool,
    in_enable_simd: bool,
    in_LapseEvolutionOption: str,
    in_ShiftEvolutionOption: str,
) -> None:
    """
    Register the right-hand side evaluation function for the BSSN equations.

    :param in_CoordSystem: The coordinate system to be used.
    :param in_enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param in_enable_simd: Whether or not to enable SIMD (Single Instruction, Multiple Data).
    :param in_LapseEvolutionOption: The lapse function evolution option.
    :param in_ShiftEvolutionOption: The shift vector evolution option.
    :return: None
    """
    includes = ["BHaH_defines.h"]
    if in_enable_simd:
        includes += [os.path.join("simd", "simd_intrinsics.h")]
    desc = r"""Set RHSs for the BSSN evolution equations."""
    c_type = "void"
    name = "rhs_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs"
    if in_enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    # Populate BSSN rhs variables
    rhs = BSSN_RHSs[
        CoordSystem + "_rfm_precompute" if in_enable_rfm_precompute else CoordSystem
    ]
    BSSN_RHSs_varnames = rhs.BSSN_RHSs_varnames
    BSSN_RHSs_exprs = rhs.BSSN_RHSs_exprs
    alpha_rhs, vet_rhsU, bet_rhsU = BSSN_gauge_RHSs(
        in_CoordSystem,
        in_enable_rfm_precompute,
        LapseEvolutionOption=in_LapseEvolutionOption,
        ShiftEvolutionOption=in_ShiftEvolutionOption,
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
        in_CoordSystem + "_rfm_precompute"
        if in_enable_rfm_precompute
        else in_CoordSystem
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
            enable_simd=in_enable_simd,
            upwind_control_vec=betaU,
        ),
        loop_region="interior",
        enable_simd=in_enable_simd,
        CoordSystem=in_CoordSystem,
        enable_rfm_precompute=in_enable_rfm_precompute,
        read_xxs=not in_enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        enable_simd=in_enable_simd,
    )


def register_CFunction_Ricci_eval(
    in_CoordSystem: str,
    in_enable_rfm_precompute: bool,
    in_enable_simd: bool,
    in_OMP_collapse: int,
) -> None:
    """
    Register the Ricci evaluation function.

    :param in_CoordSystem: The coordinate system to be used.
    :param in_enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param in_enable_simd: Whether or not to enable SIMD instructions.
    :param in_OMP_collapse: Degree of OpenMP loop collapsing.
    :return: None
    """
    orig_LeaveRicciSymbolic = par.parval_from_str("LeaveRicciSymbolic")
    if orig_LeaveRicciSymbolic:
        del BSSN_quantities[
            in_CoordSystem + "_rfm_precompute"
            if in_enable_rfm_precompute
            else in_CoordSystem
        ]

        par.set_parval_from_str("LeaveRicciSymbolic", False)
    Bq = BSSN_quantities[
        in_CoordSystem + "_rfm_precompute"
        if in_enable_rfm_precompute
        else in_CoordSystem
    ]

    includes = ["BHaH_defines.h"]
    if in_enable_simd:
        includes += [os.path.join("simd", "simd_intrinsics.h")]
    desc = r"""Set Ricci tensor."""
    c_type = "void"
    name = "Ricci_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict in_gfs, REAL *restrict auxevol_gfs"
    if in_enable_rfm_precompute:
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
            enable_simd=in_enable_simd,
        ),
        loop_region="interior",
        enable_simd=in_enable_simd,
        CoordSystem=in_CoordSystem,
        enable_rfm_precompute=in_enable_rfm_precompute,
        read_xxs=not in_enable_rfm_precompute,
        OMP_collapse=in_OMP_collapse,
    )

    if orig_LeaveRicciSymbolic:
        par.set_parval_from_str("LeaveRicciSymbolic", orig_LeaveRicciSymbolic)
        del BSSN_quantities[
            in_CoordSystem + "_rfm_precompute"
            if in_enable_rfm_precompute
            else in_CoordSystem
        ]
        _ = BSSN_quantities[
            in_CoordSystem + "_rfm_precompute"
            if in_enable_rfm_precompute
            else in_CoordSystem
        ]

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        enable_simd=in_enable_simd,
    )


def register_CFunction_constraints(
    in_CoordSystem: str,
    in_enable_rfm_precompute: bool,
    in_enable_simd: bool,
    in_OMP_collapse: int,
) -> None:
    """
    Register the BSSN constraints evaluation function.

    :param in_CoordSystem: The coordinate system to be used.
    :param in_enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param in_enable_simd: Whether or not to enable SIMD instructions.
    :param in_OMP_collapse: Degree of OpenMP loop collapsing.
    :return: None
    """
    Bcon = BSSN_constraints[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    includes = ["BHaH_defines.h"]
    if enable_simd:
        includes += [os.path.join("simd", "simd_intrinsics.h")]
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
            enable_simd=in_enable_simd,
        ),
        loop_region="interior",
        enable_simd=in_enable_simd,
        CoordSystem=in_CoordSystem,
        enable_rfm_precompute=in_enable_rfm_precompute,
        read_xxs=not in_enable_rfm_precompute,
        OMP_collapse=in_OMP_collapse,
    )

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        enable_simd=enable_simd,
    )


def register_CFunction_enforce_detgammabar_equals_detgammahat(
    in_CoordSystem: str,
    in_enable_rfm_precompute: bool,
    in_OMP_collapse: int,
) -> None:
    """
    Register the function that enforces the det(gammabar) = det(gammahat) constraint.

    :param in_CoordSystem: The coordinate system to be used.
    :param in_enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param in_OMP_collapse: Degree of OpenMP loop collapsing.
    :return: None
    """
    Bq = BSSN_quantities[
        in_CoordSystem + "_rfm_precompute"
        if in_enable_rfm_precompute
        else in_CoordSystem
    ]
    rfm = refmetric.reference_metric[
        in_CoordSystem + "_rfm_precompute"
        if in_enable_rfm_precompute
        else in_CoordSystem
    ]

    includes = ["BHaH_defines.h"]
    desc = r"""Enforce det(gammabar) = det(gammahat) constraint. Required for strong hyperbolicity."""
    c_type = "void"
    name = "enforce_detgammabar_equals_detgammahat"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs"
    if in_enable_rfm_precompute:
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
    body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            hprimeDD_expr_list,
            hDD_access_gfs,
            enable_fd_codegen=True,
            enable_simd=False,
        ),
        loop_region="all points",
        enable_simd=False,
        CoordSystem=in_CoordSystem,
        enable_rfm_precompute=in_enable_rfm_precompute,
        read_xxs=not in_enable_rfm_precompute,
        OMP_collapse=in_OMP_collapse,
    )

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        enable_simd=False,
    )


register_CFunction_initial_data()
numgrids.register_CFunction_numerical_grids_and_timestep_setup(
    CoordSystem, grid_physical_size, Nxx_dict
)
register_CFunction_diagnostics(in_CoordSystem=CoordSystem)
if enable_rfm_precompute:
    rfm_precompute.register_CFunctions_rfm_precompute(CoordSystem)
print("Constructing rhs_eval C function...")
register_CFunction_rhs_eval(
    in_CoordSystem=CoordSystem,
    in_enable_rfm_precompute=enable_rfm_precompute,
    in_enable_simd=enable_simd,
    in_LapseEvolutionOption=LapseEvolutionOption,
    in_ShiftEvolutionOption=ShiftEvolutionOption,
)
print("Finished constructing rhs_eval C function.")
print("Started constructing Ricci C function.")
register_CFunction_Ricci_eval(
    in_CoordSystem=CoordSystem,
    in_enable_rfm_precompute=enable_rfm_precompute,
    in_enable_simd=enable_simd,
    in_OMP_collapse=OMP_collapse,
)
print("Finished constructing Ricci C function.")
register_CFunction_enforce_detgammabar_equals_detgammahat(
    in_CoordSystem=CoordSystem,
    in_enable_rfm_precompute=enable_rfm_precompute,
    in_OMP_collapse=OMP_collapse,
)
print("Started constructing constraints C function.")
register_CFunction_constraints(
    in_CoordSystem=CoordSystem,
    in_enable_rfm_precompute=enable_rfm_precompute,
    in_enable_simd=enable_simd,
    in_OMP_collapse=OMP_collapse,
)
print("Finished constructing constraints C function.")

cbc.CurviBoundaryConditions_register_C_functions(
    CoordSystem, radiation_BC_fd_order=radiation_BC_fd_order
)
rhs_string = """
Ricci_eval(commondata, params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
if (strncmp(commondata->outer_bc_type, "radiation", 50) == 0)
  apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata->xx,
                                     gridfunctions_wavespeed,gridfunctions_f_infinity,
                                     RK_INPUT_GFS, RK_OUTPUT_GFS);"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")
MoL.MoL_register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_string="""if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);
  enforce_detgammabar_equals_detgammahat(commondata, params, rfmstruct, RK_OUTPUT_GFS);""",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
)
xxCart.register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
xxCart.register_CFunction_xx_to_Cart(CoordSystem)
progress.register_CFunction_progress_indicator()

#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
if CoordSystem == "SinhSpherical":
    par.adjust_CodeParam_default("SINHW", 0.4)
par.adjust_CodeParam_default("eta", GammaDriving_eta)

CPs.write_CodeParameters_h_files(project_dir=project_dir)
CPs.register_CFunctions_params_commondata_struct_set_to_default()
cmdpar.generate_default_parfile(project_dir=project_dir, project_name=project_name)
cmdpar.register_CFunction_cmdline_input_and_parfile_parser(
    project_name=project_name, cmdline_inputs=["convergence_factor"]
)
Bdefines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    enable_simd=enable_simd,
    fin_NGHOSTS_add_one_for_upwinding=True,
    MoL_method=MoL_method,
    CoordSystem=CoordSystem,
)
main.register_CFunction_main_c(
    MoL_method=MoL_method,
    # initial_data_desc=WaveType,
    # initial_data_function_call="initial_data(&griddata);",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
    boundary_conditions_desc="Quadratic extrapolation, manually defined",
)

if enable_simd:
    project_path = Path("project") / project_name
    simd_path = project_path / "simd"
    simd_path.mkdir(parents=True, exist_ok=True)

    try:
        import importlib.resources  # pylint: disable=E1101  # for Python 3.7+

        source_path = importlib.resources.files("nrpy.helpers") / "simd_intrinsics.h"
        shutil.copy(str(source_path), str(simd_path))
    except ImportError:  # Fallback to resource_filename for older Python versions
        from pkg_resources import resource_filename  # type: ignore

        source_path = resource_filename("nrpy.helpers", "simd_intrinsics.h")
        shutil.copy(source_path, str(simd_path))  # type: ignore

Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir, project_name=project_name, exec_name=project_name
)
print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")

# print(cfc.CFunction_dict["initial_data"].full_function)
# print(cfc.CFunction_dict["rhs_eval"].full_function)
# print(cfc.CFunction_dict["apply_bcs"].full_function)
# print(cfc.CFunction_dict["parameter_file_read_and_parse"].full_function)
# print(cfc.CFunction_dict["main"].full_function)
