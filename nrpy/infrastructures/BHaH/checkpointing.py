"""
Register CFunctions read_checkpoint and write_checkpoint.

Provides checkpointing capabilities to BHaH simulations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Tuple

import nrpy.c_function as cfc
import nrpy.params as par


def register_CFunction_read_checkpoint(
    filename_tuple: Tuple[str, str] = (
        r"checkpoint-conv_factor%.2f.dat",
        "commondata->convergence_factor",
    ),
    enable_bhahaha: bool = False,
) -> None:
    """
    Register read_checkpoint CFunction for reading checkpoints.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
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
#define FREAD(ptr, size, nmemb, stream) \
  MAYBE_UNUSED const int numitems=fread((ptr), (size), (nmemb), (stream));
"""
    prefunc += (
        r"""
    #define BHAH_CHKPT_HOST_MOL_GF_FREE(gf_ptr) \
      CUDA__free_host_gfs(gf_ptr); \
      CUDA__free_host_diagnostic_gfs(gf_ptr);
    #define BHAH_CHKPT_HOST_MOL_GF_MALLOC(cd, params_ptr, gf_ptr) \
      CUDA__malloc_host_gfs(cd, params_ptr, gf_ptr); \
      CUDA__malloc_host_diagnostic_gfs(cd, params_ptr, gf_ptr);
    #define BHAH_CHKPT_DEVICE_MOL_GF_FREE(gf_ptr) \
      MoL_free_memory_y_n_gfs(gf_ptr);
    #define BHAH_CHKPT_DEVICE_MOL_GF_MALLOC(cd, params_ptr, gf_ptr) \
      MoL_malloc_y_n_gfs(cd, params_ptr, gf_ptr);
    #define BHAH_CHKPT_CPY_HOST_TO_DEVICE_PARAMS() \
      memcpy(&griddata_device[grid].params, &griddata[grid].params, sizeof(params_struct));
    #define BHAH_CHKPT_CPY_HOST_TO_DEVICE_ALL_GFS() \
    for (int gf = 0; gf < NUM_EVOL_GFS; ++gf) { \
      cpyHosttoDevice__gf(commondata, &griddata[grid].params, griddata[grid].gridfuncs.y_n_gfs, griddata_device[grid].gridfuncs.y_n_gfs, gf, gf, griddata[grid].params.grid_idx % NUM_STREAMS); \
    }
    """
        if parallelization == "cuda"
        else r"""
    #define BHAH_CHKPT_HOST_MOL_GF_FREE(gf_ptr)
      MoL_free_memory_y_n_gfs(gf_ptr);
    #define BHAH_CHKPT_HOST_MOL_GF_MALLOC(cd, params_ptr, gf_ptr)
      MoL_malloc_y_n_gfs(cd, params_ptr, gf_ptr);
    """
    )
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
    body = rf"""
  char filename[256];
  snprintf(filename, 256, "{filename_tuple[0]}", {filename_tuple[1]});
"""
    body += r"""  // If the checkpoint doesn't exist then return 0.
  if (access(filename, F_OK) != 0)
    return 0;

  FILE *cp_file = fopen(filename, "r");
  FREAD(commondata, sizeof(commondata_struct), 1, cp_file);
  fprintf(stderr, "cd struct size = %ld time=%e\n", sizeof(commondata_struct), commondata->time);
"""
    if enable_bhahaha:
        body += r"""
  for (int i = 0; i < commondata->bah_max_num_horizons; i++) {
    FREAD(&commondata->bhahaha_params_and_data[i], sizeof(bhahaha_params_and_data_struct), 1, cp_file);
    commondata->bhahaha_params_and_data[i].prev_horizon_m1 = malloc(sizeof(REAL) * 64 * 32);
    commondata->bhahaha_params_and_data[i].prev_horizon_m2 = malloc(sizeof(REAL) * 64 * 32);
    commondata->bhahaha_params_and_data[i].prev_horizon_m3 = malloc(sizeof(REAL) * 64 * 32);
    FREAD(commondata->bhahaha_params_and_data[i].prev_horizon_m1, sizeof(REAL), 64 * 32, cp_file);
    FREAD(commondata->bhahaha_params_and_data[i].prev_horizon_m2, sizeof(REAL), 64 * 32, cp_file);
    FREAD(commondata->bhahaha_params_and_data[i].prev_horizon_m3, sizeof(REAL), 64 * 32, cp_file);
  }
"""

    body += r"""
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    FREAD(&griddata[grid].params, sizeof(params_struct), 1, cp_file);

    int count;
    FREAD(&count, sizeof(int), 1, cp_file);

    int * out_data_indices = (int *)malloc(sizeof(int) * count);
    REAL * compact_out_data = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);

    const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
    const int ntot = griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2;
    fprintf(stderr, "Reading checkpoint: grid = %d | pts = %d / %d | %d\n", grid, count, ntot, Nxx_plus_2NGHOSTS2);
    FREAD(out_data_indices, sizeof(int), count, cp_file);
    FREAD(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);
"""
    if parallelization in ["cuda"]:
        body += r"""
    BHAH_CHKPT_CPY_HOST_TO_DEVICE_PARAMS();
    BHAH_CHKPT_DEVICE_MOL_GF_FREE(&griddata_device[grid].gridfuncs);
    BHAH_CHKPT_DEVICE_MOL_GF_MALLOC(commondata, &griddata_device[grid].params, &griddata_device[grid].gridfuncs);
"""
    body += r"""
    BHAH_CHKPT_HOST_MOL_GF_FREE(&griddata[grid].gridfuncs);
    BHAH_CHKPT_HOST_MOL_GF_MALLOC(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
#pragma omp parallel for
    for (int i = 0; i < count; i++) {
      for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {
        griddata[grid].gridfuncs.y_n_gfs[IDX4pt(gf, out_data_indices[i])] = compact_out_data[i * NUM_EVOL_GFS + gf];
      }
    }
    free(out_data_indices);
    free(compact_out_data);
    BHAH_CHKPT_CPY_HOST_TO_DEVICE_ALL_GFS();
  }
  fclose(cp_file);
  fprintf(stderr, "FINISHED WITH READING\n");

  // Next set t_0 and n_0
  commondata->t_0 = commondata->time;
  commondata->nn_0 = commondata->nn;
""".replace(
        "BHAH_CHKPT_CPY_HOST_TO_DEVICE_ALL_GFS();",
        (
            ""
            if parallelization not in ["cuda"]
            else "BHAH_CHKPT_CPY_HOST_TO_DEVICE_ALL_GFS();"
        ),
    )
    if parallelization in ["cuda"]:
        body += "BHAH_DEVICE_SYNC();\n"
    body += "return 1;\n"
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
    filename_tuple: Tuple[str, str] = (
        "checkpoint-conv_factor%.2f.dat",
        "commondata->convergence_factor",
    ),
    enable_bhahaha: bool = False,
) -> None:
    """
    Register write_checkpoint CFunction for writing checkpoints.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
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
    prefunc = (
        r"""
    #define BHAH_CHKPT_HOST_MOL_GF_FREE(gf_ptr) \
      CUDA__free_host_gfs(gf_ptr);
    #define BHAH_CHKPT_HOST_MOL_GF_MALLOC(cd, params_ptr, gf_ptr) \
      CUDA__malloc_host_gfs(cd, params_ptr, gf_ptr);
    #define BHAH_CPY_DEVICE_TO_HOST_PARAMS() \
      memcpy(&griddata[grid].params, &griddata_device[grid].params, sizeof(params_struct));
    #define BHAH_CPY_DEVICE_TO_HOST_ALL_GFS() \
    for (int gf = 0; gf < NUM_EVOL_GFS; ++gf) { \
      size_t stream_id = griddata_device[grid].params.grid_idx % NUM_STREAMS; \
      cpyDevicetoHost__gf(commondata, &griddata[grid].params, griddata[grid].gridfuncs.y_n_gfs, griddata_device[grid].gridfuncs.y_n_gfs, gf, gf, stream_id); \
    }
    """
        if parallelization == "cuda"
        else r"""
    #define BHAH_CHKPT_HOST_MOL_GF_FREE(gf_ptr) \
      MoL_free_memory_y_n_gfs(gf_ptr);
    #define BHAH_CHKPT_HOST_MOL_GF_MALLOC(cd, params_ptr, gf_ptr) \
      MoL_malloc_y_n_gfs(cd, params_ptr, gf_ptr);
    """
    )
    body = rf"""
  char filename[256];
  snprintf(filename, 256, "{filename_tuple[0]}", {filename_tuple[1]});
"""
    body += r"""const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->checkpoint_every;
// Explanation of the if() below:
// Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime.
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
  }
"""
    body += r"""
  for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
"""
    if parallelization in ["cuda"]:
        body += r"""
    BHAH_CPY_DEVICE_TO_HOST_PARAMS();
    BHAH_CPY_DEVICE_TO_HOST_ALL_GFS();
"""
    body += r"""
    fwrite(&griddata[grid].params, sizeof(params_struct), 1, cp_file);
    const int ntot = ( griddata[grid].params.Nxx_plus_2NGHOSTS0*
                       griddata[grid].params.Nxx_plus_2NGHOSTS1*
                       griddata[grid].params.Nxx_plus_2NGHOSTS2 );
    // First we free up memory so we can malloc more: copy y_n_gfs to diagnostic_output_gfs & then free y_n_gfs.
#pragma omp parallel for
    for(int i=0;i<ntot*NUM_EVOL_GFS;i++) griddata[grid].gridfuncs.diagnostic_output_gfs[i] = griddata[grid].gridfuncs.y_n_gfs[i];
    BHAH_CHKPT_HOST_MOL_GF_FREE(&griddata[grid].gridfuncs);

    int count = 0;
    const int maskval = 1; // to be replaced with griddata[grid].mask[i].
#pragma omp parallel for reduction(+:count)
    for(int i=0;i<ntot;i++) {
      if(maskval >= +0) count++;
    }
    fwrite(&count, sizeof(int), 1, cp_file);

    int  * out_data_indices = (int  *)malloc(sizeof(int)                 * count);
    REAL * compact_out_data = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);
    int which_el = 0;
    BHAH_DEVICE_SYNC();
    for(int i=0;i<ntot;i++) {
      if(maskval >= +0) {
        out_data_indices[which_el] = i;
        for(int gf=0; gf<NUM_EVOL_GFS; gf++)
          compact_out_data[which_el*NUM_EVOL_GFS + gf] = griddata[grid].gridfuncs.diagnostic_output_gfs[ntot*gf + i];
        which_el++;
      }
    }
    //printf("HEY which_el = %d | count = %d\n", which_el, count);
    fwrite(out_data_indices, sizeof(int) , count               , cp_file);
    fwrite(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);
    free(out_data_indices); free(compact_out_data);
    BHAH_CHKPT_HOST_MOL_GF_MALLOC(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
#pragma omp parallel for
    for(int i=0;i<ntot*NUM_EVOL_GFS;i++) griddata[grid].gridfuncs.y_n_gfs[i] = griddata[grid].gridfuncs.diagnostic_output_gfs[i];
  }
  fclose(cp_file);
  fprintf(stderr, "FINISHED WRITING CHECKPOINT\n");
}
""".replace(
        "BHAH_DEVICE_SYNC();",
        ("" if parallelization not in ["cuda"] else "BHAH_DEVICE_SYNC();"),
    )
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
    filename_tuple: Tuple[str, str] = (
        "checkpoint-conv_factor%.2f.dat",
        "commondata->convergence_factor",
    ),
    default_checkpoint_every: float = 2.0,
    enable_bhahaha: bool = False,
) -> None:
    """
    Register CFunctions for checkpointing.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
    :param default_checkpoint_every: The default checkpoint interval in physical time units.
    :param enable_bhahaha: Whether to enable BHaHAHA.
    """
    register_CFunction_read_checkpoint(
        filename_tuple=filename_tuple, enable_bhahaha=enable_bhahaha
    )
    register_CFunction_write_checkpoint(
        filename_tuple=filename_tuple,
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
