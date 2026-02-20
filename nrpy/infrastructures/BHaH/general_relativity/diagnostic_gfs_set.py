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
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


def register_CFunction_diagnostic_gfs_set(
    enable_interp_diagnostics: bool, enable_psi4: bool, enable_T4munu: bool = False
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
    :param enable_T4munu:            If True, copy latest T4munu to DIAG_T4UU gridfunctions.
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

    gri.register_gridfunctions(
        names="DIAG_HAMILTONIAN", desc="H_constraint", group="DIAG"
    )
    gri.register_gridfunctions(names="DIAG_MSQUARED", desc="M^2", group="DIAG")
    gri.register_gridfunctions(names="DIAG_LAPSE", desc="Lapse", group="DIAG")
    gri.register_gridfunctions(names="DIAG_W", desc="Conformal_factor_W", group="DIAG")
    if enable_psi4:
        gri.register_gridfunctions(names="DIAG_PSI4_RE", desc="Psi4_Re", group="DIAG")
        gri.register_gridfunctions(names="DIAG_PSI4_IM", desc="Psi4_Im", group="DIAG")
    gri.register_gridfunctions_for_single_rank2(
        "DIAG_RBARDD",
        desc="Ricci_tensor_component_RbarDD",
        symmetry="sym01",
        dimension=3,
        group="DIAG",
    )
    if enable_T4munu:
        gri.register_gridfunctions_for_single_rank2(
            "DIAG_T4UU",
            desc="Stress_energy_tensor_component_T4UU",
            symmetry="sym01",
            dimension=4,
            group="DIAG",
        )
        
    diag_gf_parity_types = gri.BHaHGridFunction.set_parity_types(
        sorted([v.name for v in gri.glb_gridfcs_dict.values() if v.group == "DIAG"])
    )
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

    body = f"MAYBE_UNUSED const int8_t diag_gf_parities[{len(diag_gf_parity_types)}] = {{ {', '.join(map(str, diag_gf_parity_types))} }};\n"
    body += """
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    const params_struct *restrict params = &griddata[grid].params;
    const rfm_struct *restrict rfmstruct = griddata[grid].rfmstruct;
    const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;

    // Poison diagnostic_gfs (for debugging purposes only; WARNING: this might make valgrind ineffective)
    // #pragma omp parallel for
    //     for (int ii = 0; ii < TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2; ii++) {
    //       diagnostic_gfs[grid][ii] = NAN;
    //     } // END LOOP over all points & gridfunctions, poisoning diagnostic_gfs
"""
    parallelization = par.parval_from_str("parallelization")
    ricci_call = (
        "Ricci_eval_host(params, rfmstruct, y_n_gfs, diagnostic_gfs[grid]);"
        if parallelization == "cuda"
        else "Ricci_eval(params, rfmstruct, y_n_gfs, auxevol_gfs);"
    )
    body += f"""
    // Set Ricci and constraints gridfunctions
    {ricci_call}
    constraints_eval(commondata, params, rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_gfs[grid]);
    """
    if enable_psi4:
        body += """
    // NOTE: Inner boundary conditions must be set before any interpolations are performed, whether for psi4 decomp. or interp diags.
    // Set psi4 gridfunctions
    psi4(commondata, params, (REAL * restrict*)griddata[grid].xx, y_n_gfs, diagnostic_gfs[grid]);
    const int inner_bc_apply_gfs[] = {DIAG_PSI4_REGF, DIAG_PSI4_IMGF};
    const int num_inner_bc_apply_gfs = (int)(sizeof(inner_bc_apply_gfs) / sizeof(inner_bc_apply_gfs[0]));
    apply_bcs_inner_only_specific_gfs(commondata, params, &griddata[grid].bcstruct, diagnostic_gfs[grid], num_inner_bc_apply_gfs, diag_gf_parities,
                                      inner_bc_apply_gfs);
"""
    if enable_interp_diagnostics:
        body += """
    {
      // NOTE: Inner boundary conditions must be set before any interpolations are performed, whether for psi4 decomp. or interp diags.
      // Apply inner bcs to constraints needed to do interpolation correctly
      const int inner_bc_apply_gfs[] = {DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF};
      const int num_inner_bc_apply_gfs = (int)(sizeof(inner_bc_apply_gfs) / sizeof(inner_bc_apply_gfs[0]));
      apply_bcs_inner_only_specific_gfs(commondata, params, &griddata[grid].bcstruct, diagnostic_gfs[grid], num_inner_bc_apply_gfs, diag_gf_parities,
                                        inner_bc_apply_gfs);
    } // END set inner BCs on desired GFs
"""
    body += "  } // END LOOP over grids\n"

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
