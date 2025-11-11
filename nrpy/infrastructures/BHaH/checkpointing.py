"""
Register CFunctions read_checkpoint and write_checkpoint.

Provides checkpointing capabilities to BHaH simulations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc
import nrpy.params as par


def register_CFunction_read_checkpoint(enable_bhahaha: bool = False) -> None:
    """
    Register read_checkpoint CFunction for reading checkpoints.

    :param enable_bhahaha: Whether to enable BHaHAHA.

    Doctest:
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
    """
    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "unistd.h"]
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
              "read_checkpoint: FATAL: error while reading %s (%s): expected %zu items, got %zu.\n",         \
              (filename), (context), _expected, _got);                                                       \
      exit(EXIT_FAILURE);                                                                                    \
    }                                                                                                        \
  } while (0)

#define BHAH_CHKPT_CPY_HOST_TO_DEVICE_ALL_GFS() \
for (int gf = 0; gf < NUM_EVOL_GFS; gf++) { \
  cpyHosttoDevice__gf(commondata, &griddata[grid].params, griddata[grid].gridfuncs.y_n_gfs, griddata_device[grid].gridfuncs.y_n_gfs, gf, gf, griddata[grid].params.grid_idx % NUM_STREAMS); \
} // END LOOP over gridfunctions
"""
    desc = "Read a checkpoint file"
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
    body = r"""
  char filename[256];
  snprintf(filename, 256, "checkpoint-conv_factor%.2f.dat", commondata->convergence_factor);

  // If the checkpoint doesn't exist then return 0.
  if (access(filename, F_OK) != 0)
    return 0;

  FILE *cp_file = fopen(filename, "r");
  FREAD(commondata, sizeof(commondata_struct), 1, cp_file, filename, "commondata_struct");
  fprintf(stderr, "cd struct size = %ld time=%e\n", sizeof(commondata_struct), commondata->time);
"""
    if enable_bhahaha:
        body += r"""
  for (int i = 0; i < commondata->bah_max_num_horizons; i++) {
    FREAD(&commondata->bhahaha_params_and_data[i], sizeof(bhahaha_params_and_data_struct), 1, cp_file, filename, "bhahaha_params_and_data[i]");
    commondata->bhahaha_params_and_data[i].prev_horizon_m1 = malloc(sizeof(REAL) * 64 * 32);
    commondata->bhahaha_params_and_data[i].prev_horizon_m2 = malloc(sizeof(REAL) * 64 * 32);
    commondata->bhahaha_params_and_data[i].prev_horizon_m3 = malloc(sizeof(REAL) * 64 * 32);
    FREAD(commondata->bhahaha_params_and_data[i].prev_horizon_m1, sizeof(REAL), 64 * 32, cp_file, filename, "bhahaha_params_and_data[i].prev_horizon_m1");
    FREAD(commondata->bhahaha_params_and_data[i].prev_horizon_m2, sizeof(REAL), 64 * 32, cp_file, filename, "bhahaha_params_and_data[i].prev_horizon_m2");
    FREAD(commondata->bhahaha_params_and_data[i].prev_horizon_m3, sizeof(REAL), 64 * 32, cp_file, filename, "bhahaha_params_and_data[i].prev_horizon_m3");
  } // END LOOP over horizons
"""

    body += r"""
  for(int grid=0; grid < commondata->NUMGRIDS; grid++) {
    FREAD(&griddata[grid].params, sizeof(params_struct), 1, cp_file, filename, "griddata[grid].params");

    int count;
    FREAD(&count, sizeof(int), 1, cp_file, filename, "count");

    int *out_data_indices = (int *)malloc(sizeof(int) * count);
    REAL *compact_out_data = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);

    const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
    const int ntot_grid =
        griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2;
    fprintf(stderr, "Reading checkpoint: grid = %d | pts = %d / %d | %d\n", grid, count, ntot_grid, Nxx_plus_2NGHOSTS2);
    FREAD(out_data_indices, sizeof(int), count, cp_file, filename, "out_data_indices");
    FREAD(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file, filename, "compact_out_data");
#ifdef __CUDACC__
    memcpy(&griddata_device[grid].params, &griddata[grid].params, sizeof(params_struct));
    BHAH_FREE_DEVICE(griddata_device[grid].gridfuncs.y_n_gfs);
    BHAH_MALLOC_DEVICE(griddata_device[grid].gridfuncs.y_n_gfs, sizeof(REAL) * ntot_grid * NUM_EVOL_GFS);
    BHAH_FREE_PINNED(griddata[grid].gridfuncs.y_n_gfs);
    BHAH_MALLOC_PINNED(griddata[grid].gridfuncs.y_n_gfs, sizeof(REAL) * ntot_grid * NUM_EVOL_GFS);
#else
    BHAH_FREE(griddata[grid].gridfuncs.y_n_gfs);
    BHAH_MALLOC(griddata[grid].gridfuncs.y_n_gfs, sizeof(REAL) * ntot_grid * NUM_EVOL_GFS);
#endif // __CUDACC__

#pragma omp parallel for
    for (int i = 0; i < count; i++) {
      for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {
        griddata[grid].gridfuncs.y_n_gfs[IDX4pt(gf, out_data_indices[i])] = compact_out_data[i * NUM_EVOL_GFS + gf];
      } // END LOOP over gfs
    } // END LOOP over points
    free(out_data_indices);
    free(compact_out_data);

    IFCUDARUN(BHAH_CHKPT_CPY_HOST_TO_DEVICE_ALL_GFS());

  } // END LOOP over grids
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
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def register_CFunction_write_checkpoint(
    default_checkpoint_every: float = 2.0,
    enable_bhahaha: bool = False,
) -> None:
    """
    Register write_checkpoint CFunction for writing checkpoints.

    :param default_checkpoint_every: The default checkpoint interval in physical time units.
    :param enable_bhahaha: Whether to enable BHaHAHA.

    Doctest:
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
    desc = "Write a checkpoint file"
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
"""
    body = r"""
  char filename[256];
  snprintf(filename, 256, "checkpoint-conv_factor%.2f.dat", commondata->convergence_factor);

  const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->checkpoint_every;
  // Explanation of the if() below:
  // Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime/outevery.
  // Step 2: Multiplying by outevery yields the exact time we should output again, t_out.
  // Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible!
  if (fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt) {
    FILE *cp_file = fopen(filename, "w+");
    fwrite(commondata, sizeof(commondata_struct), 1, cp_file);
    fprintf(stderr, "WRITING CHECKPOINT: cd struct size = %ld time=%e\n", sizeof(commondata_struct), commondata->time);
"""
    if enable_bhahaha:
        body += r"""
    for (int i = 0; i < commondata->bah_max_num_horizons; i++) {
      fwrite(&commondata->bhahaha_params_and_data[i], sizeof(bhahaha_params_and_data_struct), 1, cp_file);
      fwrite(commondata->bhahaha_params_and_data[i].prev_horizon_m1, sizeof(REAL), 64 * 32, cp_file);
      fwrite(commondata->bhahaha_params_and_data[i].prev_horizon_m2, sizeof(REAL), 64 * 32, cp_file);
      fwrite(commondata->bhahaha_params_and_data[i].prev_horizon_m3, sizeof(REAL), 64 * 32, cp_file);
    } // END LOOP over all apparent horizons
"""
    body += r"""
  for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
      const int ntot_grid =
          (griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2);

#ifdef __CUDACC__
    BHAH_CPY_DEVICE_TO_HOST_PARAMS();
    BHAH_CPY_DEVICE_TO_HOST_ALL_GFS();
#endif // __CUDACC__
      fwrite(&griddata[grid].params, sizeof(params_struct), 1, cp_file);

      // First we free up memory so we can malloc more.
      MoL_free_intermediate_stage_gfs(&griddata[grid].gridfuncs);

      int count = 0;
      const int maskval = 1; // to be replaced with griddata[grid].mask[i].
#pragma omp parallel for reduction(+ : count)
      for (int i = 0; i < ntot_grid; i++) {
        if (maskval >= +0)
          count++;
      } // END LOOP over all gridpoints
      fwrite(&count, sizeof(int), 1, cp_file);

      int *out_data_indices = (int *)malloc(sizeof(int) * count);
      REAL *compact_out_data = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);
      int which_el = 0;

      for (int i = 0; i < ntot_grid; i++) {
        if (maskval >= +0) {
          out_data_indices[which_el] = i;
          for (int gf = 0; gf < NUM_EVOL_GFS; gf++)
            compact_out_data[which_el * NUM_EVOL_GFS + gf] = griddata[grid].gridfuncs.y_n_gfs[ntot_grid * gf + i];
          which_el++;
        } // END IF maskval >= +0
      } // END LOOP over all gridpoints

      fwrite(out_data_indices, sizeof(int), count, cp_file);
      fwrite(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);
      free(out_data_indices);
      free(compact_out_data);

      // Re-allocate memory for intermediate-stage scratch storage for MoL.
      MoL_malloc_intermediate_stage_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
    } // END LOOP over grids
    fclose(cp_file);
    fprintf(stderr, "FINISHED WRITING CHECKPOINT\n");
  } // END IF ready to write checkpoint
"""
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def register_CFunctions(
    default_checkpoint_every: float = 2.0,
    enable_bhahaha: bool = False,
) -> None:
    """
    Register CFunctions for checkpointing.

    :param default_checkpoint_every: The default checkpoint interval in physical time units.
    :param enable_bhahaha: Whether to enable BHaHAHA.
    """
    register_CFunction_read_checkpoint(enable_bhahaha=enable_bhahaha)
    register_CFunction_write_checkpoint(
        default_checkpoint_every=default_checkpoint_every,
        enable_bhahaha=enable_bhahaha,
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
