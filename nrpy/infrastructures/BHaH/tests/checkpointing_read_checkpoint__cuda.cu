#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "unistd.h"

// clang formatting disabled due to issues with brace placement:
//   placing the opening brace on the same line and other times on a new line, which causes CI failures.
// clang-format off
#define FREAD(ptr, size, nmemb, stream) \
  { MAYBE_UNUSED const int numitems=fread((ptr), (size), (nmemb), (stream)); }
// clang-format on

#define BHAH_CHKPT_CPY_HOST_TO_DEVICE_ALL_GFS()                                                                                                      \
  for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {                                                                                                        \
    cpyHosttoDevice__gf(commondata, &griddata[grid].params, griddata[grid].gridfuncs.y_n_gfs, griddata_device[grid].gridfuncs.y_n_gfs, gf, gf,       \
                        griddata[grid].params.grid_idx % NUM_STREAMS);                                                                               \
  } // END LOOP over gridfunctions

/**
 * Read a checkpoint file
 */
int read_checkpoint(commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_device) {
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
    const int ntot_grid =
        griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2;
    fprintf(stderr, "Reading checkpoint: grid = %d | pts = %d / %d | %d\n", grid, count, ntot_grid, Nxx_plus_2NGHOSTS2);
    FREAD(out_data_indices, sizeof(int), count, cp_file);
    FREAD(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);
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
    BHAH_CHKPT_CPY_HOST_TO_DEVICE_ALL_GFS();
  } // END LOOP over grids
  fclose(cp_file);
  fprintf(stderr, "FINISHED WITH READING\n");

  // Next set t_0 and n_0
  commondata->t_0 = commondata->time;
  commondata->nn_0 = commondata->nn;
  BHAH_DEVICE_SYNC();
  return 1;
} // END FUNCTION read_checkpoint
