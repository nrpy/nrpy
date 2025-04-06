#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "unistd.h"

#define FREAD(ptr, size, nmemb, stream) MAYBE_UNUSED const int numitems = fread((ptr), (size), (nmemb), (stream));

#define BHAH_CHKPT_HOST_MOL_GF_FREE(gf_ptr)
MoL_free_memory_y_n_gfs(gf_ptr);
#define BHAH_CHKPT_HOST_MOL_GF_MALLOC(cd, params_ptr, gf_ptr)
MoL_malloc_y_n_gfs(cd, params_ptr, gf_ptr);

/**
 * Read a checkpoint file
 */
int read_checkpoint(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {

  char filename[256];
  snprintf(filename, 256, "checkpoint-conv_factor%.2f.dat", commondata->convergence_factor);
  // If the checkpoint doesn't exist then return 0.
  if (access(filename, F_OK) != 0)
    return 0;

  FILE *cp_file = fopen(filename, "r");
  FREAD(commondata, sizeof(commondata_struct), 1, cp_file);
  fprintf(stderr, "cd struct size = %ld time=%e\n", sizeof(commondata_struct), commondata->time);

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    FREAD(&griddata[grid].params, sizeof(params_struct), 1, cp_file);

    int count;
    FREAD(&count, sizeof(int), 1, cp_file);

    int *out_data_indices = (int *)malloc(sizeof(int) * count);
    REAL *compact_out_data = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);

    const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
    const int ntot = griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2;
    fprintf(stderr, "Reading checkpoint: grid = %d | pts = %d / %d | %d\n", grid, count, ntot, Nxx_plus_2NGHOSTS2);
    FREAD(out_data_indices, sizeof(int), count, cp_file);
    FREAD(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);

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
  }
  fclose(cp_file);
  fprintf(stderr, "FINISHED WITH READING\n");

  // Next set t_0 and n_0
  commondata->t_0 = commondata->time;
  commondata->nn_0 = commondata->nn;
  return 1;
} // END FUNCTION read_checkpoint
