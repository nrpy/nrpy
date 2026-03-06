#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "unistd.h"

// Safe fread wrapper: verifies that the expected number of items were read.
// On mismatch, prints a descriptive error to stderr and aborts the program.
// This is intentionally fatal: partial/failed checkpoint reads are treated
// as unrecoverable corruption, not a soft "no restart".
#define FREAD(ptr, size, nmemb, stream, filename, context)                                                                                           \
  do {                                                                                                                                               \
    size_t _expected = (size_t)(nmemb);                                                                                                              \
    size_t _got = fread((ptr), (size), (nmemb), (stream));                                                                                           \
    if (_got != _expected) {                                                                                                                         \
      fprintf(stderr, "read_checkpoint: FATAL: error while reading %s (%s): expected %zu items, got %zu.\n", (filename), (context), _expected,       \
              _got);                                                                                                                                 \
      exit(EXIT_FAILURE);                                                                                                                            \
    }                                                                                                                                                \
  } while (0)

/**
 * Read a checkpoint file
 */
int read_checkpoint(commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_device) {
  char filename[256];
  snprintf(filename, 256, "checkpoint-conv_factor%.2f.dat", commondata->convergence_factor);

  // If the checkpoint doesn't exist then return 0; if it does exist and can't be read, then error out.
  FILE *cp_file = fopen(filename, "r");
  if (cp_file == NULL) {
    if (errno == ENOENT)
      return 0; // checkpoint doesn't exist
    fprintf(stderr, "read_checkpoint: FATAL: could not open %s for reading: %s\n", filename, strerror(errno));
    exit(EXIT_FAILURE);
  } // END IF fopen failed

  FREAD(commondata, sizeof(commondata_struct), 1, cp_file, filename, "commondata_struct");
  fprintf(stderr, "cd struct size = %zu time=%e\n", sizeof(commondata_struct), commondata->time);

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    FREAD(&griddata[grid].params, sizeof(params_struct), 1, cp_file, filename, "griddata[grid].params");

    int count;
    FREAD(&count, sizeof(int), 1, cp_file, filename, "count");

    int *restrict out_data_indices;
    BHAH_MALLOC(out_data_indices, sizeof(int) * count);
    REAL *restrict compact_out_data;
    BHAH_MALLOC(compact_out_data, sizeof(REAL) * NUM_EVOL_GFS * count);

    const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
    const int ntot_grid = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

    if (count < 0 || count > ntot_grid) {
      fprintf(stderr, "read_checkpoint: FATAL: invalid count=%d for grid=%d (ntot_grid=%d) in %s\n", count, grid, ntot_grid, filename);
      exit(EXIT_FAILURE);
    } // END IF count is invalid, error out.
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

    IFCUDARUN(for (int gf = 0; gf < NUM_EVOL_GFS; gf++)
                  cpyHosttoDevice__gf(commondata, &griddata[grid].params, griddata[grid].gridfuncs.y_n_gfs, griddata_device[grid].gridfuncs.y_n_gfs,
                                      gf, gf, griddata[grid].params.grid_idx % NUM_STREAMS););

  } // END LOOP over grids
  fclose(cp_file);
  fprintf(stderr, "FINISHED WITH READING\n");

  // Next set t_0 and n_0
  commondata->t_0 = commondata->time;
  commondata->nn_0 = commondata->nn;

  IFCUDARUN(BHAH_DEVICE_SYNC());

  return 1;
} // END FUNCTION read_checkpoint
