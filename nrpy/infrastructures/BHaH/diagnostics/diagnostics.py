"""
C function registration for the top-level diagnostics() driver.

This module constructs and registers the C diagnostics() driver and the commondata
parameter diagnostics_output_every. The generated C code determines whether the
current step is an output step from time, dt, and diagnostics_output_every; allocates
and initializes diagnostic_gfs; conditionally calls enabled diagnostics routines
(nearest, interpolation, volume integration, and optional extensions); releases
temporary storage; and advances the progress indicator. Hooks to free and later
restore MoL scratch arrays are present but currently disabled in the emitted code.

Preconditions / constraints (documentation only; not enforced here):
- diagnostics_output_every > 0 (the scheduling rule divides by this value).
- NUMGRIDS <= MAXNUMGRIDS (the generated driver uses a fixed-size pointer array
  sized by MAXNUMGRIDS and loops over NUMGRIDS).

Functions
---------
register_all_diagnostics
    Public entry point. Stages required helper headers (nearest + volume integration
    only, in this module), registers the "diagnostics()" C driver, and registers
    family-specific helpers for the families that require per-coordinate registration
    here.
_register_CFunction_diagnostics
    Internal helper that generates and registers the "diagnostics()" C driver and
    the CodeParameter "diagnostics_output_every". This helper participates in the
    parallel-codegen registration phase.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Set, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures import BHaH


def _register_CFunction_diagnostics(  # pylint: disable=unused-argument
    default_diagnostics_out_every: float,
    enable_nearest_diagnostics: bool,
    enable_interp_diagnostics: bool,
    enable_volume_integration_diagnostics: bool,
    enable_free_auxevol: bool = True,
    enable_psi4_diagnostics: bool = False,
    enable_bhahaha: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register a C function that drives all scheduled diagnostics.

    Two-phase behavior:
    - During the parallel-codegen registration phase, this function records the call
      (including argument values) and returns None.
    - Outside the registration phase, this function generates and registers the C
      driver "diagnostics" and returns the updated NRPy environment.

    The generated driver is invoked once per timestep. It decides whether the current
    time is within half a timestep of the nearest multiple of diagnostics_output_every,
    and on output steps it allocates and populates diagnostic_gfs via
    diagnostic_gfs_set(...), invokes the enabled diagnostics routines
    (diagnostics_nearest(...), diagnostics_interp(...), diagnostics_volume_integration(...),
    and optional extensions), then frees temporary storage.

    Scheduling note:
    The rule uses a tolerance window of 0.5 * dt around the nominal output times.
    With floating-point time accumulation and/or variable dt, an output may occur up to
    that tolerance early/late relative to exact multiples of diagnostics_output_every.

    CUDA note (generation-time vs compile-time):
    - If the code is generated with parallelization="cuda", the driver signature includes
      griddata_device. Separately, device-to-host transfers in the body are compiled only
      under __CUDACC__.
    - On output steps in the CUDA path, the driver copies all NUM_EVOL_GFS evolution
      gridfunctions from device to host prior to diagnostics that require host-side I/O,
      with additional optional transfers guarded by feature macros (e.g., T4UU00GF).

    Preconditions / constraints (documentation only; not enforced here):
    - diagnostics_output_every > 0 (the scheduling rule divides by this value).
    - NUMGRIDS <= MAXNUMGRIDS (the driver uses a fixed-size pointer array sized by
      MAXNUMGRIDS and loops over NUMGRIDS).

    :param default_diagnostics_out_every: Default value for the commondata parameter
        "diagnostics_output_every", which controls the diagnostics output cadence.
        Must satisfy diagnostics_output_every > 0.
    :param enable_nearest_diagnostics: If True, include a call to diagnostics_nearest(...)
        on output steps.
    :param enable_interp_diagnostics: If True, include a call to diagnostics_interp(...)
        on output steps.
    :param enable_volume_integration_diagnostics: If True, include a call to
        diagnostics_volume_integration(...) on output steps.
    :param enable_free_auxevol: Currently ignored and has no effect on emitted C code.
        The MoL scratch free/restore hooks remain commented out in the generated driver.
    :param enable_psi4_diagnostics: If True, include a call to
        psi4_spinweightm2_decomposition(...) on output steps.
    :param enable_bhahaha: If True, include a call to bhahaha_find_horizons(...) on output steps.
    :return: None if in registration phase (after recording the requested registration),
        else the updated NRPy environment.

    Doctests:
    None.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # --- C Function Registration ---
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "diagnostics/diagnostic_gfs.h",
    ]
    desc = """
 * @file diagnostics.c
 * @brief Top-level driver that schedules and runs all enabled diagnostics.
 *
 * The function "diagnostics" is generated by NRPy and invoked once per timestep.
 * It checks whether the current time falls on a diagnostics output step. On output
 * steps it allocates per-grid temporary arrays (diagnostic_gfs), initializes them
 * via diagnostic_gfs_set(...), and runs the enabled diagnostics routines
 * (for example diagnostics_nearest(...), diagnostics_interp(...),
 * diagnostics_volume_integration(...), and optional extensions). Temporary storage is
 * freed before returning. Independently of output steps, a progress indicator is
 * advanced every call, and a trailing newline is printed when the run is about to finish.
 *
 * Scheduling rule (time-based with tolerance window):
 *   fabs(round(time / diagnostics_output_every) * diagnostics_output_every - time) < 0.5 * dt
 *
 * CUDA note (generation-time vs compile-time):
 * - If this file is generated with parallelization="cuda", the signature includes
 *   griddata_device. Separately, device-to-host transfers are compiled only under
 *   __CUDACC__.
 * - On output steps in the CUDA path, the driver copies all NUM_EVOL_GFS evolution
 *   gridfunctions from device to host prior to diagnostics that require host-side I/O,
 *   with additional optional transfers guarded by feature macros (e.g., T4UU00GF).
 *
 * @pre
 * - commondata and griddata are non-null and initialized.
 * - commondata->NUMGRIDS >= 1 and commondata->NUMGRIDS <= MAXNUMGRIDS.
 * - diagnostics_output_every > 0.
 * - grid dimensions are valid, and backing arrays exist.
 * - Generated diagnostics interfaces and indices are consistent with the build.
 *
 * @post
 * - On output steps, all enabled diagnostics execute and may write output.
 * - On non-output steps, no diagnostic I/O occurs and the solution state is unchanged.
 * - The progress indicator advances every call; a newline is printed if time + dt > t_final.
 *
 * @param[in,out] commondata  Global simulation metadata and run-time parameters
 *                            (e.g., time, dt, diagnostics_output_every, t_final, NUMGRIDS).
 * @param[in,out] griddata_device  Device-side per-grid data used for device-to-host
 *                                 synchronization when generated with parallelization="cuda";
 *                                 omitted otherwise.
 * @param[in,out] griddata    Host-side per-grid data (parameters, fields, and workspace).
 *
 * @warning
 * - Diagnostics that encounter allocation or I/O failures may abort the program.
 * - The set of diagnostics compiled in is fixed at code generation time; manual changes
 *   must remain consistent with generated headers and prototypes.
 *
 * @return void
"""
    parallelization = par.parval_from_str("parallelization")
    _ = par.CodeParameter(
        "REAL",
        __name__,
        "diagnostics_output_every",
        default_diagnostics_out_every,
        commondata=True,
    )
    cfunc_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    if parallelization == "cuda":
        params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata_device, griddata_struct *restrict griddata"
    newline = "\n"  # Keep newline sequences as named constants to avoid escape/brace pitfalls inside f-strings.
    rnewline = "\\n"  # Keep newline sequences as named constants to avoid escape/brace pitfalls inside f-strings.
    body = f"""
  const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->diagnostics_output_every;
  // Explanation of the if() below:                                                                                                                                                                                                           
  // Step 1: round(currtime / outevery) gives the nearest integer n to the ratio currtime/outevery.                                                                                                                                           
  // Step 2: Multiplying by outevery yields the nearest output time t_out = n * outevery.                                                                                                                                                     
  // Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible.                                                                                                                                        
  // Note: This rule divides by outevery; outevery must be > 0.
  if (fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt) {{
    // Optional hook: free MoL scratch storage before diagnostics and restore afterward.
    // Currently disabled in emitted code (calls are commented out below).
    // for (int grid = 0; grid < commondata->NUMGRIDS; grid++)
    //   MoL_free_intermediate_stage_gfs(&griddata[grid].gridfuncs);

    {"// Find apparent horizon(s)." if enable_bhahaha else ""}
    {"bhahaha_find_horizons(commondata, griddata);" + newline if enable_bhahaha else ""}

    // Allocate temporary storage for diagnostic_gfs.
    REAL *diagnostic_gfs[MAXNUMGRIDS];
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
      SET_NXX_PLUS_2NGHOSTS_VARS(grid);
      const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
      BHAH_MALLOC(diagnostic_gfs[grid], TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS_tot * sizeof(REAL));

#ifdef __CUDACC__
      // This does not leverage async memory transfers using multiple streams at the moment
      // given the current intent is one cuda stream per grid. This could be leveraged
      // in the future by increasing NUM_STREAMS such that a diagnostic stream is included per grid
      const params_struct *restrict params = &griddata_device[grid].params;
      size_t streamid = params->grid_idx % NUM_STREAMS;
      cpyHosttoDevice_params__constant(&griddata_device[grid].params, streamid);
      // Copy solution to host (all evolution gridfunctions).
      for(int gf=0;gf<NUM_EVOL_GFS;gf++)
        cpyDevicetoHost__gf(commondata, params, griddata[grid].gridfuncs.y_n_gfs, griddata_device[grid].gridfuncs.y_n_gfs, gf,gf, streamid);
#ifdef T4UU00GF
      // Transfer the 10 T4UU gridfunctions from device to host (optional; requires T4UU00GF).
      const int idx0_host = IDX4pt(DIAG_T4UU00GF, 0), idx0_device = IDX4pt(T4UU00GF, 0);
      cudaMemcpyAsync(&diagnostic_gfs[grid][idx0_host], &griddata_device[grid].gridfuncs.auxevol_gfs[idx0_device],
                      10 * Nxx_plus_2NGHOSTS_tot * sizeof(REAL), cudaMemcpyDeviceToHost, streams[streamid]);
#endif // T4UU00GF
      // Sync data before attempting to write to file
      cudaStreamSynchronize(streams[streamid]);
#endif // __CUDACC__
    }} // END LOOP over grids

    // Set diagnostic_gfs; see generated diagnostics/diagnostic_gfs.h for the interface.
    diagnostic_gfs_set(commondata, griddata, diagnostic_gfs);

    {"// Nearest-point diagnostics, at center, along y,z axes (1D) and xy and yz planes (2D)." if enable_nearest_diagnostics else ""}
    {"diagnostics_nearest(commondata, griddata, (const REAL **)diagnostic_gfs);" + newline if enable_nearest_diagnostics else ""}
    {"// Interpolation diagnostics, at center, along x,y,z axes (1D) and xy and yz planes (2D)." if enable_interp_diagnostics else ""}
    {"diagnostics_interp(commondata, griddata, (const REAL **)diagnostic_gfs);" + newline if enable_interp_diagnostics else ""}
    {"// Volume-integration diagnostics." if enable_volume_integration_diagnostics else ""}
    {"diagnostics_volume_integration(commondata, griddata, (const REAL **)diagnostic_gfs);" + newline if enable_volume_integration_diagnostics else ""}
    {"// Decompose psi4 into modes." if enable_psi4_diagnostics else ""}
    {"psi4_spinweightm2_decomposition(commondata, griddata, (const REAL **)diagnostic_gfs);" + newline if enable_psi4_diagnostics else ""}

    // Free temporary storage allocated to diagnostic_gfs.
    for(int grid=0; grid<commondata->NUMGRIDS; grid++)
      free(diagnostic_gfs[grid]);

    // Optional hook: restore MoL scratch storage after diagnostics (currently disabled).
    // for (int grid = 0; grid < commondata->NUMGRIDS; grid++)
    //   MoL_malloc_intermediate_stage_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
  }} // END if output step

  progress_indicator(commondata, griddata);
  if (commondata->time + commondata->dt > commondata->t_final) printf("{rnewline}");
"""
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


def register_all_diagnostics(
    project_dir: str,
    set_of_CoordSystems: Set[str],
    default_diagnostics_out_every: float,
    enable_nearest_diagnostics: bool,
    enable_interp_diagnostics: bool,
    enable_volume_integration_diagnostics: bool,
    enable_free_auxevol: bool = True,
    enable_psi4_diagnostics: bool = False,
    enable_bhahaha: bool = False,
) -> None:
    """
    Register and stage diagnostics-related C code and helper headers.

    What this function does in this module:
    - Copies helper headers for the nearest and volume-integration diagnostics families
      into the project's diagnostics/ subdirectory (when enabled).
    - Registers the top-level diagnostics() C driver with the requested default output cadence.
    - Registers per-coordinate helper C functions for:
        - nearest diagnostics (several sampling helpers), and
        - volume integration (volume-element helper).

    What this function does NOT do here:
    - It does not copy/register any interpolation-specific helper headers or per-coordinate
      helper functions in this module.
    - It does not stage/register psi4 or bhahaha helper code in this module; enabling those
      flags only affects whether the top-level driver emits calls on output steps.

    :param project_dir: Target project directory for emitted diagnostics assets.
    :param set_of_CoordSystems: Set of coordinate-system names to generate helpers for.
    :param default_diagnostics_out_every: Default value for the commondata parameter
        "diagnostics_output_every" controlling output cadence. Must satisfy
        diagnostics_output_every > 0.
    :param enable_nearest_diagnostics: If True, include sampling-based "nearest" diagnostics.
    :param enable_interp_diagnostics: If True, include interpolation-based diagnostics (driver call only).
    :param enable_volume_integration_diagnostics: If True, include volume-integration diagnostics.
    :param enable_free_auxevol: Currently ignored and has no effect on emitted C code.
        The MoL scratch free/restore hooks remain commented out in the generated driver.
    :param enable_psi4_diagnostics: If True, decompose psi4 into spin-weight -2 spherical harmonics
        (driver call only).
    :param enable_bhahaha: If True, include a call to bhahaha_find_horizons(...) on output steps.

    Doctests:
    None.
    """
    filenames_list_to_copy = []
    if enable_nearest_diagnostics:
        filenames_list_to_copy += ["diagnostics_nearest_common.h"]
    if enable_volume_integration_diagnostics:
        filenames_list_to_copy += ["diagnostics_volume_integration_helpers.h"]
    if filenames_list_to_copy:
        copy_files(
            package="nrpy.infrastructures.BHaH.diagnostics",
            filenames_list=filenames_list_to_copy,
            project_dir=project_dir,
            subdirectory="diagnostics",
        )

    _register_CFunction_diagnostics(
        default_diagnostics_out_every=default_diagnostics_out_every,
        enable_nearest_diagnostics=enable_nearest_diagnostics,
        enable_interp_diagnostics=enable_interp_diagnostics,
        enable_volume_integration_diagnostics=enable_volume_integration_diagnostics,
        enable_free_auxevol=enable_free_auxevol,
        enable_psi4_diagnostics=enable_psi4_diagnostics,
        enable_bhahaha=enable_bhahaha,
    )
    if enable_nearest_diagnostics:
        for CoordSystem in set_of_CoordSystems:
            BHaH.diagnostics.diagnostics_nearest_grid_center.register_CFunction_diagnostics_nearest_grid_center(
                CoordSystem=CoordSystem
            )
            BHaH.diagnostics.diagnostics_nearest_1d_y_and_z_axes.register_CFunction_diagnostics_nearest_1d_y_and_z_axes(
                CoordSystem=CoordSystem
            )
            BHaH.diagnostics.diagnostics_nearest_2d_xy_and_yz_planes.register_CFunction_diagnostics_nearest_2d_xy_and_yz_planes(
                CoordSystem=CoordSystem
            )
    if enable_volume_integration_diagnostics:
        for CoordSystem in set_of_CoordSystems:
            BHaH.diagnostics.sqrt_detgammahat_d3xx_volume_element.register_CFunction_sqrt_detgammahat_d3xx_volume_element(
                CoordSystem=CoordSystem
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
