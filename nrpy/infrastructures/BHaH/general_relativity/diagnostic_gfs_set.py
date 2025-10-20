"""
Register and emit the C `diagnostic_gfs_set()` routine for interpolation and integration diagnostics.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


def register_CFunction_diagnostic_gfs_set(enable_psi4: bool) -> Union[None, pcg.NRPyEnv_type]:
    """
    Generate and register the C `diagnostic_gfs_set()` routine that fills per-grid diagnostic arrays.

    :returns: None during the registration phase; otherwise the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "diagnostics/diagnostic_gfs.h",
    ]
    desc = """
 * @brief Populate diagnostic_gfs[][] for interpolation & integration diagnostics.
"""
    cfunc_type = "void"
    name = "diagnostic_gfs_set"
    params = "const commondata_struct *restrict commondata, const griddata_struct *restrict griddata, REAL *restrict diagnostic_gfs[MAXNUMGRIDS]"

    diagnostic_gfs_names_dict = {
        "DIAG_HAMILTONIAN": "H_constraint",
        "DIAG_MSQUAREDGF": "M^2",
        "DIAG_LAPSE": "Lapse",
        "DIAG_W": "Conformal_factor_W",
    }
    if enable_psi4:
        diagnostic_gfs_names_dict.update({
            "DIAG_PSI4_RE": "Psi4_Re",
            "DIAG_PSI4_IM": "Psi4_Im",
        })

    body = """
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    const params_struct *restrict params = &griddata[grid].params;
    SET_NXX_PLUS_2NGHOSTS_VARS(grid);
    const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;

    // Poison diagnostic_gfs (for debugging purposes only; WARNING: this might make valgrind ineffective)
    // #pragma omp parallel for
    //     for (int ii = 0; ii < TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2; ii++) {
    //       diagnostic_gfs[grid][ii] = NAN;
    //     } // END LOOP over all points & gridfunctions, poisoning diagnostic_gfs
    """
    if enable_psi4:
        body += """
"""
    body += """

    // Apply inner bcs to DIAG_RESIDUAL, as it depends on finite differences.
    {
      const bc_struct bcstruct = griddata[grid].bcstruct;
      const innerpt_bc_struct *restrict inner_bc_array = bcstruct.inner_bc_array;
      const int num_inner_boundary_points = bcstruct.bc_info.num_inner_boundary_points;
#pragma omp parallel for
      for (int pt = 0; pt < num_inner_boundary_points; pt++) {
        const int dstpt = inner_bc_array[pt].dstpt;
        const int srcpt = inner_bc_array[pt].srcpt;
        const int evol_gf_with_same_parity = UUGF; // <- IMPORTANT
        diagnostic_gfs[grid][IDX4pt(DIAG_RESIDUAL, dstpt)] =
            inner_bc_array[pt].parity[evol_gf_parity[evol_gf_with_same_parity]] * diagnostic_gfs[grid][IDX4pt(DIAG_RESIDUAL, srcpt)];
      } // END LOOP over inner boundary points
    } // END applying inner bcs to DIAG_RESIDUAL

    LOOP_OMP("omp parallel for", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      const int idx3 = IDX3(i0, i1, i2);
      diagnostic_gfs[grid][IDX4pt(DIAG_UUGF, idx3)] = y_n_gfs[IDX4pt(UUGF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_VVGF, idx3)] = y_n_gfs[IDX4pt(VVGF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_GRIDINDEX, idx3)] = (REAL)grid;
    } // END LOOP over all gridpoints to set diagnostic_gfs
  } // END LOOP over grids
"""
    par.glb_extras_dict["diagnostic_gfs_names_dict"] = diagnostic_gfs_names_dict

    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
