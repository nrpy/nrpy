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
) -> None:
    """
    Register read_checkpoint CFunction for reading checkpoints.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
    :return: None
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "unistd.h"]
    prefunc = r"""
#define FREAD(ptr, size, nmemb, stream) { const int numitems=fread((ptr), (size), (nmemb), (stream)); }
"""
    desc = "Read a checkpoint file"
    c_type = "int"
    name = "read_checkpoint"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
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
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    FREAD(&griddata[grid].params, sizeof(params_struct), 1, cp_file);

    int count;
    FREAD(&count, sizeof(int), 1, cp_file);

    int *restrict out_data_indices = (int *restrict)malloc(sizeof(int) * count);
    REAL *restrict compact_out_data = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);

    const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
    const int ntot = griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2;
    fprintf(stderr, "Reading checkpoint: grid = %d | pts = %d / %d | %d\n", grid, count, ntot, Nxx_plus_2NGHOSTS2);
    FREAD(out_data_indices, sizeof(int), count, cp_file);
    FREAD(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);

    MoL_malloc_y_n_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
    int which_el = 0;
#pragma omp parallel for
    for (int i = 0; i < count; i++) {
      for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {
        griddata[grid].gridfuncs.y_n_gfs[IDX4pt(gf, out_data_indices[i])] = compact_out_data[i * NUM_EVOL_GFS + gf];
      }
    }
    free(out_data_indices);
    free(compact_out_data);
  }
  fclose(cp_file);
  fprintf(stderr, "FINISHED WITH READING\n");

  // Next set t_0 and n_0
  commondata->t_0 = commondata->time;
  commondata->nn_0 = commondata->nn;

  return 1;
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type,
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
) -> None:
    """
    Register write_checkpoint CFunction for writing checkpoints.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
    :param default_checkpoint_every: The default checkpoint interval in physical time units.
    :return: None
    """
    par.register_CodeParameter(
        "REAL", __name__, "checkpoint_every", default_checkpoint_every, commondata=True
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Write a checkpoint file"
    c_type = "void"
    name = "write_checkpoint"
    params = "const commondata_struct *restrict commondata, griddata_struct *restrict griddata"

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
  for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
    fwrite(&griddata[grid].params, sizeof(params_struct), 1, cp_file);
    const int ntot = ( griddata[grid].params.Nxx_plus_2NGHOSTS0*
                       griddata[grid].params.Nxx_plus_2NGHOSTS1*
                       griddata[grid].params.Nxx_plus_2NGHOSTS2 );
    // First we free up memory so we can malloc more: copy y_n_gfs to diagnostic_output_gfs & then free y_n_gfs.
#pragma omp parallel for
    for(int i=0;i<ntot*NUM_EVOL_GFS;i++) griddata[grid].gridfuncs.diagnostic_output_gfs[i] = griddata[grid].gridfuncs.y_n_gfs[i];
    MoL_free_memory_y_n_gfs(&griddata[grid].gridfuncs);

    int count = 0;
    const int maskval = 1; // to be replaced with griddata[grid].mask[i].
#pragma omp parallel for reduction(+:count)
    for(int i=0;i<ntot;i++) {
      if(maskval >= +0) count++;
    }
    fwrite(&count, sizeof(int), 1, cp_file);

    int  *restrict out_data_indices = (int  *restrict)malloc(sizeof(int)                 * count);
    REAL *restrict compact_out_data = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);
    int which_el = 0;
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
    MoL_malloc_y_n_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
#pragma omp parallel for
    for(int i=0;i<ntot*NUM_EVOL_GFS;i++) griddata[grid].gridfuncs.y_n_gfs[i] = griddata[grid].gridfuncs.diagnostic_output_gfs[i];
  }
  fclose(cp_file);
  fprintf(stderr, "FINISHED WRITING CHECKPOINT\n");
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


def register_CFunctions(
    filename_tuple: Tuple[str, str] = (
        "checkpoint-conv_factor%.2f.dat",
        "commondata->convergence_factor",
    ),
    default_checkpoint_every: float = 2.0,
) -> None:
    """
    Register CFunctions for checkpointing.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
    :param default_checkpoint_every: The default checkpoint interval in physical time units.
    :return: None
    """
    register_CFunction_read_checkpoint(filename_tuple=filename_tuple)
    register_CFunction_write_checkpoint(
        filename_tuple=filename_tuple, default_checkpoint_every=default_checkpoint_every
    )
