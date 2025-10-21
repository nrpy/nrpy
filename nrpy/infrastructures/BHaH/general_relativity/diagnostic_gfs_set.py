"""
C function registration for populating per-grid diagnostic arrays used by interpolation and integration routines.

This module constructs and registers the C routine "diagnostic_gfs_set".
The generated C function iterates over all grids to fill per-grid diagnostic arrays for
downstream interpolation and integration routines: it computes a residual-type diagnostic;
optionally applies inner boundary conditions using parity-consistent signs; copies selected
evolved gridfunctions from y_n_gfs into designated diagnostic channels; and records a
per-point grid identifier.

Function
--------
register_CFunction_diagnostic_gfs_set
    Construct and register the "diagnostic_gfs_set" C function.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


def register_CFunction_diagnostic_gfs_set(
    enable_interp_diagnostics: bool, enable_psi4: bool
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register a C function that populates per-grid diagnostic arrays used by interpolation and integration routines.

    This function generates and registers the C helper "diagnostic_gfs_set", which loops over all
    grids and fills per-grid diagnostic arrays as follows: it computes a residual-type diagnostic at
    all points; if enable_interp_diagnostics is True, it applies an inner-boundary pass using parity
    to ensure consistent values near symmetry or excision boundaries; it copies selected current-time
    level evolved gridfunctions from y_n_gfs into designated diagnostic channels; it sets a per-point
    grid-identifier channel; and, if enable_psi4 is True, it also registers additional waveform-related
    diagnostic channels produced by separate code paths. Each per-grid output buffer is assumed to hold
    TOTAL_NUM_DIAG_GFS times the number of points in that grid.

    :param enable_interp_diagnostics: If True, apply an inner-boundary parity pass to the residual-type
                                      diagnostic; if False, skip that step.
    :param enable_psi4:              If True, include additional waveform-related diagnostic channels
                                      in the registration; if False, omit them.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    TBD
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
 * @file diagnostic_gfs_set.c
 * @brief Populate per-grid diagnostic arrays used by interpolation and integration routines.
 *
 * The function "diagnostic_gfs_set" loops over all grids and fills per-grid diagnostic arrays:
 *   1) Compute a residual-type diagnostic at all points using a helper that evaluates the
 *      finite-difference residual.
 *   2) If enabled at code-generation time, apply inner boundary conditions to that residual by
 *      copying from a source point to a destination point with a sign determined by the relevant
 *      parity, ensuring parity-consistent values near symmetry or excision boundaries.
 *   3) Copy selected evolved gridfunctions from the current time level (y_n_gfs) into designated
 *      diagnostic channels for downstream consumers.
 *   4) Set a per-point grid identifier channel to the grid index (converted to REAL).
 *
 * The routine assumes each per-grid output buffer is contiguous and large enough to store all
 * diagnostic channels:
 *     TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2
 * Loops over grid points may be parallelized with OpenMP if available.
 *
 * If a user-editable block is present in the implementation, users may add custom logic such as
 * extra diagnostics or filtering before finalizing values.
 *
 * @param[in]  commondata
 *     Pointer to global simulation metadata (e.g., counters and configuration) accessed by helpers
 *     and used to determine the number of grids to process.
 * @param[in]  griddata
 *     Pointer to an array of per-grid data. For each grid, this provides parameters, coordinates,
 *     boundary condition metadata, and gridfunctions (including y_n_gfs and any auxiliary data)
 *     referenced by this routine and its helpers.
 * @param[out] diagnostic_gfs
 *     Array of per-grid output buffers. For each grid, diagnostic_gfs[grid] must point to a buffer
 *     of size TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2.
 *
 * @return void.
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
        diagnostic_gfs_names_dict.update(
            {
                "DIAG_PSI4_RE": "Psi4_Re",
                "DIAG_PSI4_IM": "Psi4_Im",
            }
        )

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
    if enable_interp_diagnostics:
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
"""
    body += """
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
