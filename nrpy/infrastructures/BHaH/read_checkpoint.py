# nrpy/infrastructures/BHaH/read_checkpoint.py
"""
Register the read_checkpoint CFunction.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc
import nrpy.params as par


def register_CFunction_read_checkpoint(enable_bhahaha: bool = False) -> None:
    """
    Register read_checkpoint CFunction for reading checkpoints.

    :param enable_bhahaha: Whether to enable BHaHAHA.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name="read_checkpoint"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    cfc.CFunction_dict.clear()
    ...    register_CFunction_read_checkpoint()
    ...    generated_str = cfc.CFunction_dict[name].full_function
    ...    validation_desc = f"{name}__{parallelization}".replace(" ", "_")
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_read_checkpoint(enable_bhahaha=True)
    >>> generated_str = cfc.CFunction_dict[name].full_function
    >>> generated_str.index("bah_max_num_horizons < 0") < generated_str.index("sanitize_checkpoint_commondata_pointers(commondata);")
    True
    >>> generated_str.index("n < 0 || n >= ntheta_capacity") < generated_str.index("bah_Ntheta_array_multigrid[n]")
    True
    >>> generated_str.index("count < 0 || count > ntot_grid") < generated_str.index("BHAH_MALLOC(out_data_indices")
    True
    >>> generated_str.index("if (count > 0)") < generated_str.index("FREAD(out_data_indices")
    True
    >>> generated_str.index("out_data_indices[i] < 0") < generated_str.index("compact_out_data[(size_t)i")
    True
    """
    parallelization = par.parval_from_str("parallelization")
    includes = [
        "<limits.h>",
        "<unistd.h>",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = """Read a checkpoint file.

If griddata == NULL, read only commondata metadata and return 1 when the
checkpoint exists. This supports rebuilding grids from restart state before
loading the full checkpoint payload.

@param[in,out] commondata Global state loaded from the checkpoint.
@param[in,out] griddata Rebuilt grid hierarchy, or NULL for metadata-only loading.
@return 1 when a checkpoint is read, or 0 when no checkpoint exists."""
    cfunc_type = "int"
    name = "read_checkpoint"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    params += (
        ", griddata_struct *restrict griddata_device"
        if parallelization in ["cuda"]
        else ""
    )
    prefunc = r"""
// Safe fread wrapper: verifies that the expected number of items were read.
// On mismatch, prints a descriptive error to stderr and aborts the program.
// This is intentionally fatal: partial/failed checkpoint reads are treated
// as unrecoverable corruption, not a soft "no restart".
#define FREAD(ptr, size, nmemb, stream, filename, context)                                                   \
  do {                                                                                                       \
    size_t _expected = (size_t)(nmemb);                                                                      \
    size_t _got = fread((ptr), (size), (nmemb), (stream));                                                   \
    if (_got != _expected) {                                                                                 \
      fprintf(stderr,                                                                                        \
              "read_checkpoint: FATAL: error while reading %s (%s): expected %zu items, got %zu.\n",       \
              (filename), (context), _expected, _got);                                                       \
      exit(EXIT_FAILURE);                                                                                    \
    }                                                                                                        \
  } while (0) // END DO-WHILE: verify checkpoint read count

"""
    if enable_bhahaha:
        prefunc += r"""
/**
 * Clear non-serializable BHaHAHA pointers after loading commondata.
 *
 * @param[in,out] commondata Loaded global state to sanitize.
 */
static void sanitize_checkpoint_commondata_pointers(commondata_struct *restrict commondata) {
  for (int i = 0; i < commondata->bah_max_num_horizons; i++) {
    commondata->bhahaha_params_and_data[i].input_metric_data = NULL;
    commondata->bhahaha_params_and_data[i].prev_horizon_m1 = NULL;
    commondata->bhahaha_params_and_data[i].prev_horizon_m2 = NULL;
    commondata->bhahaha_params_and_data[i].prev_horizon_m3 = NULL;
  } // END LOOP: for i over horizons
} // END FUNCTION: sanitize_checkpoint_commondata_pointers

"""

    body = r"""
  char filename[256];
  snprintf(filename, 256, "checkpoint-conv_factor%.2f.dat", commondata->convergence_factor);

  // If the checkpoint doesn't exist then return 0; if it does exist and can't be read, then error out.
  FILE *cp_file = fopen(filename, "rb");
  if (cp_file == NULL) {
    if (errno == ENOENT) return 0;  // checkpoint doesn't exist
    fprintf(stderr,
            "read_checkpoint: FATAL: could not open %s for reading: %s\n",
            filename, strerror(errno));
    exit(EXIT_FAILURE);
  } // END IF: checkpoint file cannot be opened

  FREAD(commondata, sizeof(commondata_struct), 1, cp_file, filename, "commondata_struct");

  if (commondata->NUMGRIDS <= 0 || commondata->NUMGRIDS > MAXNUMGRIDS) {
    fprintf(stderr,
            "read_checkpoint: FATAL: invalid NUMGRIDS=%d (expected 1..%d) in %s.\n",
            commondata->NUMGRIDS, MAXNUMGRIDS, filename);
    exit(EXIT_FAILURE);
  } // END IF: checkpoint grid count is invalid
"""
    if enable_bhahaha:
        body += r"""
  const int horizon_capacity =
      (int)(sizeof(commondata->bhahaha_params_and_data) / sizeof(commondata->bhahaha_params_and_data[0]));
  const int ntheta_capacity =
      (int)(sizeof(commondata->bah_Ntheta_array_multigrid) / sizeof(commondata->bah_Ntheta_array_multigrid[0]));
  const int nphi_capacity =
      (int)(sizeof(commondata->bah_Nphi_array_multigrid) / sizeof(commondata->bah_Nphi_array_multigrid[0]));
  if (commondata->bah_max_num_horizons < 0 ||
      commondata->bah_max_num_horizons > horizon_capacity ||
      commondata->bah_num_resolutions_multigrid <= 0 ||
      commondata->bah_num_resolutions_multigrid > ntheta_capacity ||
      commondata->bah_num_resolutions_multigrid > nphi_capacity) {
    fprintf(stderr,
            "read_checkpoint: FATAL: invalid BHaHAHA capacities in %s: horizons=%d/%d resolutions=%d/%d/%d.\n",
            filename, commondata->bah_max_num_horizons, horizon_capacity,
            commondata->bah_num_resolutions_multigrid, ntheta_capacity, nphi_capacity);
    exit(EXIT_FAILURE);
  } // END IF: checkpoint BHaHAHA capacities are invalid

  sanitize_checkpoint_commondata_pointers(commondata);
"""
    body += r"""
  fprintf(stderr, "cd struct size = %zu time=%e\n", sizeof(commondata_struct), commondata->time);
  if (griddata == NULL) {
    fclose(cp_file);
    return 1;
  }
"""
    if enable_bhahaha:
        body += r"""
  for (int i = 0; i < commondata->bah_max_num_horizons; i++) {
    bhahaha_params_and_data_struct *restrict horizon_params = &commondata->bhahaha_params_and_data[i];

    uint8_t has_prev_horizon_shapes = 0;
    FREAD(&has_prev_horizon_shapes, sizeof(uint8_t), 1, cp_file, filename, "has_prev_horizon_shapes");
    if (has_prev_horizon_shapes != (uint8_t)0) {
      const int n = commondata->bah_num_resolutions_multigrid - 1;
      if (n < 0 || n >= ntheta_capacity || n >= nphi_capacity) {
        fprintf(stderr,
                "read_checkpoint: FATAL: invalid horizon-shape resolution index %d in %s.\n",
                n, filename);
        exit(EXIT_FAILURE);
      } // END IF: horizon-shape resolution index is invalid

      const int ntheta_max = commondata->bah_Ntheta_array_multigrid[n];
      const int nphi_max = commondata->bah_Nphi_array_multigrid[n];
      if (ntheta_max <= 0 || nphi_max <= 0 ||
          (size_t)ntheta_max > SIZE_MAX / (size_t)nphi_max) {
        fprintf(stderr,
                "read_checkpoint: FATAL: invalid horizon-shape dimensions in commondata for horizon %d: (%d,%d).\n",
                i, ntheta_max, nphi_max);
        exit(EXIT_FAILURE);
      } // END IF: horizon-shape dimensions are invalid
      const size_t npts = (size_t)ntheta_max * (size_t)nphi_max;
      if (npts > SIZE_MAX / sizeof(REAL)) {
        fprintf(stderr,
                "read_checkpoint: FATAL: horizon-shape payload is too large for horizon %d in %s.\n",
                i, filename);
        exit(EXIT_FAILURE);
      } // END IF: horizon-shape allocation size overflows

      BHAH_MALLOC(horizon_params->prev_horizon_m1, sizeof(REAL) * npts);
      BHAH_MALLOC(horizon_params->prev_horizon_m2, sizeof(REAL) * npts);
      BHAH_MALLOC(horizon_params->prev_horizon_m3, sizeof(REAL) * npts);
      FREAD(horizon_params->prev_horizon_m1, sizeof(REAL), npts, cp_file, filename, "prev_horizon_m1");
      FREAD(horizon_params->prev_horizon_m2, sizeof(REAL), npts, cp_file, filename, "prev_horizon_m2");
      FREAD(horizon_params->prev_horizon_m3, sizeof(REAL), npts, cp_file, filename, "prev_horizon_m3");
    } // END IF: previous horizon shapes exist
  } // END LOOP: for i over horizons
"""

    body += r"""
  for(int grid=0; grid < commondata->NUMGRIDS; grid++) {
    params_struct checkpoint_params;
    FREAD(&checkpoint_params, sizeof(params_struct), 1, cp_file, filename, "griddata[grid].params");

    int count;
    FREAD(&count, sizeof(int), 1, cp_file, filename, "count");

    const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
    if (Nxx_plus_2NGHOSTS0 <= 0 || Nxx_plus_2NGHOSTS1 <= 0 || Nxx_plus_2NGHOSTS2 <= 0 ||
        (size_t)Nxx_plus_2NGHOSTS0 > SIZE_MAX / (size_t)Nxx_plus_2NGHOSTS1) {
      fprintf(stderr,
              "read_checkpoint: FATAL: invalid rebuilt dimensions on grid %d: (%d,%d,%d).\n",
              grid, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2);
      exit(EXIT_FAILURE);
    } // END IF: first rebuilt-grid dimension product is invalid
    const size_t Nxx_plus_2NGHOSTS01 = (size_t)Nxx_plus_2NGHOSTS0 * (size_t)Nxx_plus_2NGHOSTS1;
    if (Nxx_plus_2NGHOSTS01 > SIZE_MAX / (size_t)Nxx_plus_2NGHOSTS2) {
      fprintf(stderr, "read_checkpoint: FATAL: rebuilt grid %d size overflows size_t.\n", grid);
      exit(EXIT_FAILURE);
    } // END IF: full rebuilt-grid dimension product overflows
    const size_t ntot_grid_size = Nxx_plus_2NGHOSTS01 * (size_t)Nxx_plus_2NGHOSTS2;
    if (ntot_grid_size > (size_t)INT_MAX) {
      fprintf(stderr, "read_checkpoint: FATAL: rebuilt grid %d has too many points: %zu.\n", grid, ntot_grid_size);
      exit(EXIT_FAILURE);
    } // END IF: rebuilt-grid point count cannot fit its int representation
    const int ntot_grid = (int)ntot_grid_size;

    if (checkpoint_params.Nxx_plus_2NGHOSTS0 != Nxx_plus_2NGHOSTS0 ||
        checkpoint_params.Nxx_plus_2NGHOSTS1 != Nxx_plus_2NGHOSTS1 ||
        checkpoint_params.Nxx_plus_2NGHOSTS2 != Nxx_plus_2NGHOSTS2 ||
        checkpoint_params.CoordSystem_hash != griddata[grid].params.CoordSystem_hash) {
      fprintf(stderr,
              "read_checkpoint: FATAL: checkpoint/grid rebuild mismatch on grid %d: checkpoint dims=(%d,%d,%d) hash=%d, rebuilt dims=(%d,%d,%d) hash=%d\n",
              grid,
              checkpoint_params.Nxx_plus_2NGHOSTS0, checkpoint_params.Nxx_plus_2NGHOSTS1, checkpoint_params.Nxx_plus_2NGHOSTS2,
              checkpoint_params.CoordSystem_hash,
              Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
              griddata[grid].params.CoordSystem_hash);
      exit(EXIT_FAILURE);
    } // END IF: checkpoint parameters disagree with rebuilt grid metadata

    if (count < 0 || count > ntot_grid) {
      fprintf(stderr,
              "read_checkpoint: FATAL: invalid count=%d for grid=%d (ntot_grid=%d) in %s\n",
              count, grid, ntot_grid, filename);
      exit(EXIT_FAILURE);
    } // END IF: serialized grid-point count is invalid
    const size_t num_evol_gfs = (size_t)NUM_EVOL_GFS;
    if (num_evol_gfs == 0 ||
        (size_t)count > SIZE_MAX / num_evol_gfs ||
        ntot_grid_size > SIZE_MAX / num_evol_gfs) {
      fprintf(stderr, "read_checkpoint: FATAL: gridfunction count overflows on grid %d in %s.\n", grid, filename);
      exit(EXIT_FAILURE);
    } // END IF: gridfunction element count overflows
    const size_t compact_count = (size_t)count * num_evol_gfs;
    const size_t numpts_in_all_evol_gfs = ntot_grid_size * num_evol_gfs;
    if ((size_t)count > SIZE_MAX / sizeof(int) ||
        compact_count > SIZE_MAX / sizeof(REAL) || numpts_in_all_evol_gfs > SIZE_MAX / sizeof(REAL)) {
      fprintf(stderr, "read_checkpoint: FATAL: gridfunction allocation size overflows on grid %d in %s.\n", grid, filename);
      exit(EXIT_FAILURE);
    } // END IF: gridfunction byte count overflows

    int *restrict out_data_indices = NULL;
    REAL *restrict compact_out_data = NULL;
    if (count > 0) {
      BHAH_MALLOC(out_data_indices, sizeof(int) * (size_t)count);
      FREAD(out_data_indices, sizeof(int), (size_t)count, cp_file, filename, "out_data_indices");
      for (int i = 0; i < count; i++) {
        if (out_data_indices[i] < 0 || out_data_indices[i] >= ntot_grid) {
          fprintf(stderr,
                  "read_checkpoint: FATAL: invalid gridpoint index=%d at element=%d for grid=%d (ntot_grid=%d) in %s.\n",
                  out_data_indices[i], i, grid, ntot_grid, filename);
          exit(EXIT_FAILURE);
        } // END IF: serialized grid-point index is invalid
      } // END LOOP: for i over serialized grid-point indices
      BHAH_MALLOC(compact_out_data, sizeof(REAL) * compact_count);
      FREAD(compact_out_data, sizeof(REAL), compact_count, cp_file, filename, "compact_out_data");
    } // END IF: checkpoint stores points for this grid

    fprintf(stderr, "Reading checkpoint: grid = %d | pts = %d / %d | %d\n", grid, count, ntot_grid, Nxx_plus_2NGHOSTS2);
#ifdef __CUDACC__
    memcpy(&griddata_device[grid].params, &griddata[grid].params, sizeof(params_struct));
    BHAH_FREE_DEVICE(griddata_device[grid].gridfuncs.y_n_gfs);
    BHAH_MALLOC_DEVICE(griddata_device[grid].gridfuncs.y_n_gfs, sizeof(REAL) * numpts_in_all_evol_gfs);
    BHAH_FREE_PINNED(griddata[grid].gridfuncs.y_n_gfs);
    BHAH_MALLOC_PINNED(griddata[grid].gridfuncs.y_n_gfs, sizeof(REAL) * numpts_in_all_evol_gfs);
#else
    BHAH_FREE(griddata[grid].gridfuncs.y_n_gfs);
    BHAH_MALLOC(griddata[grid].gridfuncs.y_n_gfs, sizeof(REAL) * numpts_in_all_evol_gfs);
#endif // __CUDACC__

#pragma omp parallel for
    for (size_t i = 0; i < numpts_in_all_evol_gfs; i++)
      griddata[grid].gridfuncs.y_n_gfs[i] = 0.0;

#pragma omp parallel for
    for (int i = 0; i < count; i++)
      for (int gf = 0; gf < NUM_EVOL_GFS; gf++)
        griddata[grid].gridfuncs.y_n_gfs[(size_t)gf * ntot_grid_size + (size_t)out_data_indices[i]] =
            compact_out_data[(size_t)i * num_evol_gfs + (size_t)gf];
    free(out_data_indices);
    free(compact_out_data);

    IFCUDARUN(for (int gf = 0; gf < NUM_EVOL_GFS; gf++)
                  cpyHosttoDevice__gf(commondata, &griddata[grid].params, griddata[grid].gridfuncs.y_n_gfs, griddata_device[grid].gridfuncs.y_n_gfs,
                                      gf, gf, griddata[grid].params.grid_idx % NUM_STREAMS););

  } // END LOOP: for grid over active grids
  fclose(cp_file);
  fprintf(stderr, "FINISHED WITH READING\n");

  // Next set t_0 and n_0
  commondata->t_0 = commondata->time;
  commondata->nn_0 = commondata->nn;

  IFCUDARUN(BHAH_DEVICE_SYNC());

  return 1;
"""
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
