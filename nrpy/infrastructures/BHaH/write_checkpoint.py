# nrpy/infrastructures/BHaH/write_checkpoint.py
"""
Register the write_checkpoint CFunction.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc
import nrpy.params as par


def register_CFunction_write_checkpoint(
    default_checkpoint_every: float = 2.0,
    enable_multipatch: bool = False,
    enable_bhahaha: bool = False,
) -> None:
    """
    Register write_checkpoint CFunction for writing checkpoints.

    :param default_checkpoint_every: The default checkpoint interval in physical time units.
    :param enable_multipatch: Whether to enable multipatch support.
    :param enable_bhahaha: Whether to enable BHaHAHA.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name="write_checkpoint"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    cfc.CFunction_dict.clear()
    ...    register_CFunction_write_checkpoint()
    ...    generated_str = cfc.CFunction_dict[name].full_function
    ...    validation_desc = f"{name}__{parallelization}".replace(" ", "_")
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    par.register_CodeParameter(
        "REAL", __name__, "checkpoint_every", default_checkpoint_every, commondata=True
    )
    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Write a checkpoint file.

@param[in] commondata Global state to serialize.
@param[in,out] griddata Per-grid metadata and evolved gridfunctions."""
    if parallelization == "cuda":
        desc += "\n@param[in,out] griddata_device CUDA device hierarchy."
    cfunc_type = "void"
    name = "write_checkpoint"
    params = "const commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    params += (
        ", griddata_struct *restrict griddata_device"
        if parallelization in ["cuda"]
        else ""
    )
    prefunc = r"""
#ifdef __CUDACC__
#define BHAH_CPY_DEVICE_TO_HOST_PARAMS() memcpy(&griddata[grid].params, &griddata_device[grid].params, sizeof(params_struct));
#define BHAH_CPY_DEVICE_TO_HOST_ALL_GFS() \
for (int gf = 0; gf < NUM_EVOL_GFS; ++gf) { \
  size_t stream_id = griddata_device[grid].params.grid_idx % NUM_STREAMS; \
  cpyDevicetoHost__gf(commondata, &griddata[grid].params, griddata[grid].gridfuncs.y_n_gfs, griddata_device[grid].gridfuncs.y_n_gfs, gf, gf, stream_id); \
}
#endif // __CUDACC__

/**
 * Write one complete checkpoint payload or terminate on failure.
 *
 * @param[in] ptr Payload address.
 * @param size Size of one payload element.
 * @param nmemb Number of payload elements.
 * @param[in,out] fp Open checkpoint stream.
 * @param[in] what Payload description for diagnostics.
 * @param[in] file Caller source filename.
 * @param line Caller source line.
 * @param[in] func Caller function name.
 */
static inline void BHAH_safe_write_impl(const void *ptr, size_t size, size_t nmemb, FILE *fp, const char *what, const char *file, int line,
                                        const char *func) {
  clearerr(fp);
  errno = 0;

  const size_t wrote = fwrite(ptr, size, nmemb, fp);
  if (wrote != nmemb) {
    const int err = errno;
    fprintf(stderr, "%s:%d: %s: fwrite failed writing %s (wanted %zu, wrote %zu)\n", file, line, func, what, nmemb, wrote);

    if (err != 0)
      fprintf(stderr, "%s:%d: %s: errno=%d (%s)\n", file, line, func, err, strerror(err));
    else if (ferror(fp))
      fprintf(stderr, "%s:%d: %s: stream error set but errno not set\n", file, line, func);

    // Fatal: checkpoint output must not silently continue on partial write.
    exit(1);
  } // END IF: checkpoint write is incomplete
} // END FUNCTION: BHAH_safe_write_impl

#define FWRITE(ptr, size, nmemb, fp, what)                                                                                                  \
  BHAH_safe_write_impl((ptr), (size_t)(size), (size_t)(nmemb), (fp), (what), __FILE__, __LINE__, __func__)

"""

    body = r"""
  char filename[256];
  snprintf(filename, 256, "checkpoint-conv_factor%.2f.dat", commondata->convergence_factor);

  const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->checkpoint_every;
  if (outevery <= (REAL)0.0) return; // outevery <= 0 means do not checkpoint.
  // Explanation of the if() below:
  // Step 1: round(currtime / outevery) gives the nearest integer n to the ratio currtime/outevery.
  // Step 2: Multiplying by outevery yields the nearest output time t_out = n * outevery.
  // Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible.
  if (fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt) {
    FILE *cp_file = fopen(filename, "wb");
    if (cp_file == NULL) {
      perror("write_checkpoint: Failed to open checkpoint file. Check permissions and disk space availability.");
      exit(1);
    } // END IF: checkpoint file cannot be opened
"""
    if enable_bhahaha:
        body += r"""
    commondata_struct checkpoint_commondata = *commondata;
    for (int i = 0; i < checkpoint_commondata.bah_max_num_horizons; i++) {
      checkpoint_commondata.bhahaha_params_and_data[i].input_metric_data = NULL;
      checkpoint_commondata.bhahaha_params_and_data[i].prev_horizon_m1 = NULL;
      checkpoint_commondata.bhahaha_params_and_data[i].prev_horizon_m2 = NULL;
      checkpoint_commondata.bhahaha_params_and_data[i].prev_horizon_m3 = NULL;
    } // END LOOP: for i over apparent horizons in sanitized commondata copy
"""
    body += r"""
"""
    if enable_bhahaha:
        body += r"""
    FWRITE(&checkpoint_commondata, sizeof(commondata_struct), 1, cp_file, "commondata");
"""
    else:
        body += r"""
    FWRITE(commondata, sizeof(commondata_struct), 1, cp_file, "commondata");
"""
    body += r"""
    fprintf(stderr, "WRITING CHECKPOINT: cd struct size = %zu time=%e\n", sizeof(commondata_struct), commondata->time);
"""
    if enable_bhahaha:
        body += r"""
    for (int i = 0; i < commondata->bah_max_num_horizons; i++) {
      const bhahaha_params_and_data_struct *restrict horizon_params = &commondata->bhahaha_params_and_data[i];
      const int has_m1 = horizon_params->prev_horizon_m1 != NULL;
      const int has_m2 = horizon_params->prev_horizon_m2 != NULL;
      const int has_m3 = horizon_params->prev_horizon_m3 != NULL;
      const int any_prev_horizon_shapes = has_m1 || has_m2 || has_m3;
      const int all_prev_horizon_shapes = has_m1 && has_m2 && has_m3;
      if (any_prev_horizon_shapes && !all_prev_horizon_shapes) {
        fprintf(stderr,
                "write_checkpoint: FATAL: inconsistent BHaHAHA horizon-shape allocation state for horizon %d.\n",
                i);
        exit(EXIT_FAILURE);
      } // END IF: horizon-shape allocation state is inconsistent

      const uint8_t has_prev_horizon_shapes = all_prev_horizon_shapes ? (uint8_t)1 : (uint8_t)0;
      FWRITE(&has_prev_horizon_shapes, sizeof(uint8_t), 1, cp_file, "has_prev_horizon_shapes");
      if (has_prev_horizon_shapes != (uint8_t)0) {
        const int n = commondata->bah_num_resolutions_multigrid - 1;
        const int ntheta_max = commondata->bah_Ntheta_array_multigrid[n];
        const int nphi_max = commondata->bah_Nphi_array_multigrid[n];
        const size_t npts = (size_t)ntheta_max * (size_t)nphi_max;
        if (n < 0 || ntheta_max <= 0 || nphi_max <= 0) {
          fprintf(stderr, "write_checkpoint: FATAL: invalid BHaHAHA horizon-shape dimensions for horizon %d.\n", i);
          exit(EXIT_FAILURE);
        } // END IF: horizon-shape dimensions are invalid
        FWRITE(horizon_params->prev_horizon_m1, sizeof(REAL), npts, cp_file, "bhahaha_prev_horizon_m1");
        FWRITE(horizon_params->prev_horizon_m2, sizeof(REAL), npts, cp_file, "bhahaha_prev_horizon_m2");
        FWRITE(horizon_params->prev_horizon_m3, sizeof(REAL), npts, cp_file, "bhahaha_prev_horizon_m3");
      } // END IF: previous horizon shapes exist
    } // END LOOP: for i over all apparent horizons
"""
    body += r"""
  for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
      const int ntot_grid =
          (griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2);

#ifdef __CUDACC__
    BHAH_CPY_DEVICE_TO_HOST_PARAMS();
    BHAH_CPY_DEVICE_TO_HOST_ALL_GFS();
#endif // __CUDACC__
      FWRITE(&griddata[grid].params, sizeof(params_struct), 1, cp_file, "params_struct");

      // First we free up memory so we can malloc more.
      MoL_free_intermediate_stage_gfs(&griddata[grid].gridfuncs);

      int count = 0;
#pragma omp parallel for reduction(+ : count)
      for (int i = 0; i < ntot_grid; i++) {
"""
    if enable_multipatch:
        body += "const int maskval = griddata[grid].mask[i];\n"
        body += "const int owned_interior_maskval = griddata[grid].params.grid_idx;\n"
        mask_condition = "maskval == owned_interior_maskval || maskval == BUFFER_ZONE || maskval == OUTER_BOUNDARY"
    else:
        body += "const int maskval = 1;\n"
        mask_condition = "maskval >= +0"
    body += (
        r"""
        if ("""
        + mask_condition
        + r""")
          count++;
      } // END LOOP: for i over all grid points
      FWRITE(&count, sizeof(int), 1, cp_file, "gridpoint_count");

      int *restrict out_data_indices;
      BHAH_MALLOC(out_data_indices, sizeof(int) * count);
      REAL *restrict compact_out_data;
      BHAH_MALLOC(compact_out_data, sizeof(REAL) * NUM_EVOL_GFS * count);
      int which_el = 0;

      for (int i = 0; i < ntot_grid; i++) {
"""
    )
    if enable_multipatch:
        body += "const int maskval = griddata[grid].mask[i];\n"
        body += "const int owned_interior_maskval = griddata[grid].params.grid_idx;\n"
        mask_condition = "maskval == owned_interior_maskval || maskval == BUFFER_ZONE || maskval == OUTER_BOUNDARY"
    else:
        body += "const int maskval = 1;\n"
        mask_condition = "maskval >= +0"
    body += (
        r"""
        if ("""
        + mask_condition
        + r""") {
          out_data_indices[which_el] = i;
          for (int gf = 0; gf < NUM_EVOL_GFS; gf++)
            compact_out_data[which_el * NUM_EVOL_GFS + gf] = griddata[grid].gridfuncs.y_n_gfs[ntot_grid * gf + i];
          which_el++;
        } // END IF: grid point is selected for checkpoint output
      } // END LOOP: for i over all grid points

      FWRITE(out_data_indices, sizeof(int), count, cp_file, "out_data_indices");
      FWRITE(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file, "compact_out_data");
      free(out_data_indices);
      free(compact_out_data);

      // Re-allocate memory for intermediate-stage scratch storage for MoL.
      MoL_malloc_intermediate_stage_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
    } // END LOOP: for grid over active grids
    fclose(cp_file);
    fprintf(stderr, "FINISHED WRITING CHECKPOINT\n");
  } // END IF: current time is ready for checkpoint output
"""
    )
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        prefunc=prefunc,
        body=body,
        include_CodeParameters_h=False,
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
