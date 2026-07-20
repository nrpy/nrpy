#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include <limits.h>
#include <unistd.h>

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
  } while (0) // END DO-WHILE: verify checkpoint read count

// Allocate reader-owned host memory or terminate before first use.
#define READ_CHECKPOINT_MALLOC_OR_EXIT(ptr, num_bytes, filename, context)                                                                            \
  do {                                                                                                                                               \
    const size_t _read_checkpoint_num_bytes = (size_t)(num_bytes);                                                                                   \
    BHAH_MALLOC(ptr, _read_checkpoint_num_bytes);                                                                                                    \
    if ((ptr) == NULL) {                                                                                                                             \
      fprintf(stderr, "read_checkpoint: FATAL: allocation failed for %zu bytes while reading %s (%s).\n", _read_checkpoint_num_bytes, (filename),    \
              (context));                                                                                                                            \
      exit(EXIT_FAILURE);                                                                                                                            \
    }                                                                                                                                                \
  } while (0) // END DO-WHILE: allocate reader-owned host memory

/**
 * Read a checkpoint file.
 *
 * If griddata == NULL, read only commondata metadata and return 1 when the
 * checkpoint exists. This supports rebuilding grids from restart state before
 * loading the full checkpoint payload.
 *
 * Serialized grid-point indices must be strictly increasing, as emitted by the
 * current writer. Noncanonical indices are treated as fatal checkpoint
 * corruption; reader-owned host allocation failures also terminate immediately.
 *
 * @param[in,out] commondata Global state loaded from the checkpoint.
 * @param[in,out] griddata Rebuilt grid hierarchy, or NULL for metadata-only loading.
 * @return 1 when a checkpoint is read, or 0 when no checkpoint exists.
 */
int read_checkpoint(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  char filename[256];
  snprintf(filename, 256, "checkpoint-conv_factor%.2f.dat", commondata->convergence_factor);

  // If the checkpoint doesn't exist then return 0; if it does exist and can't be read, then error out.
  FILE *cp_file = fopen(filename, "rb");
  if (cp_file == NULL) {
    if (errno == ENOENT)
      return 0; // checkpoint doesn't exist
    fprintf(stderr, "read_checkpoint: FATAL: could not open %s for reading: %s\n", filename, strerror(errno));
    exit(EXIT_FAILURE);
  } // END IF: checkpoint file cannot be opened

  FREAD(commondata, sizeof(commondata_struct), 1, cp_file, filename, "commondata_struct");

  if (commondata->NUMGRIDS <= 0 || commondata->NUMGRIDS > MAXNUMGRIDS) {
    fprintf(stderr, "read_checkpoint: FATAL: invalid NUMGRIDS=%d (expected 1..%d) in %s.\n", commondata->NUMGRIDS, MAXNUMGRIDS, filename);
    exit(EXIT_FAILURE);
  } // END IF: checkpoint grid count is invalid

  fprintf(stderr, "cd struct size = %zu time=%e\n", sizeof(commondata_struct), commondata->time);
  if (griddata == NULL) {
    fclose(cp_file);
    return 1;
  }

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    params_struct checkpoint_params;
    FREAD(&checkpoint_params, sizeof(params_struct), 1, cp_file, filename, "griddata[grid].params");

    int count;
    FREAD(&count, sizeof(int), 1, cp_file, filename, "count");

    const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
    if (Nxx_plus_2NGHOSTS0 <= 0 || Nxx_plus_2NGHOSTS1 <= 0 || Nxx_plus_2NGHOSTS2 <= 0 ||
        (size_t)Nxx_plus_2NGHOSTS0 > SIZE_MAX / (size_t)Nxx_plus_2NGHOSTS1) {
      fprintf(stderr, "read_checkpoint: FATAL: invalid rebuilt dimensions on grid %d: (%d,%d,%d).\n", grid, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1,
              Nxx_plus_2NGHOSTS2);
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

    if (checkpoint_params.Nxx_plus_2NGHOSTS0 != Nxx_plus_2NGHOSTS0 || checkpoint_params.Nxx_plus_2NGHOSTS1 != Nxx_plus_2NGHOSTS1 ||
        checkpoint_params.Nxx_plus_2NGHOSTS2 != Nxx_plus_2NGHOSTS2 || checkpoint_params.CoordSystem_hash != griddata[grid].params.CoordSystem_hash) {
      fprintf(stderr,
              "read_checkpoint: FATAL: checkpoint/grid rebuild mismatch on grid %d: checkpoint dims=(%d,%d,%d) hash=%d, rebuilt dims=(%d,%d,%d) "
              "hash=%d\n",
              grid, checkpoint_params.Nxx_plus_2NGHOSTS0, checkpoint_params.Nxx_plus_2NGHOSTS1, checkpoint_params.Nxx_plus_2NGHOSTS2,
              checkpoint_params.CoordSystem_hash, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, griddata[grid].params.CoordSystem_hash);
      exit(EXIT_FAILURE);
    } // END IF: checkpoint parameters disagree with rebuilt grid metadata

    if (count < 0 || count > ntot_grid) {
      fprintf(stderr, "read_checkpoint: FATAL: invalid count=%d for grid=%d (ntot_grid=%d) in %s\n", count, grid, ntot_grid, filename);
      exit(EXIT_FAILURE);
    } // END IF: serialized grid-point count is invalid
    const size_t num_evol_gfs = (size_t)NUM_EVOL_GFS;
    if (num_evol_gfs == 0 || (size_t)count > SIZE_MAX / num_evol_gfs || ntot_grid_size > SIZE_MAX / num_evol_gfs) {
      fprintf(stderr, "read_checkpoint: FATAL: gridfunction count overflows on grid %d in %s.\n", grid, filename);
      exit(EXIT_FAILURE);
    } // END IF: gridfunction element count overflows
    const size_t compact_count = (size_t)count * num_evol_gfs;
    const size_t numpts_in_all_evol_gfs = ntot_grid_size * num_evol_gfs;
    if ((size_t)count > SIZE_MAX / sizeof(int) || compact_count > SIZE_MAX / sizeof(REAL) || numpts_in_all_evol_gfs > SIZE_MAX / sizeof(REAL)) {
      fprintf(stderr, "read_checkpoint: FATAL: gridfunction allocation size overflows on grid %d in %s.\n", grid, filename);
      exit(EXIT_FAILURE);
    } // END IF: gridfunction byte count overflows

    int *restrict out_data_indices = NULL;
    REAL *restrict compact_out_data = NULL;
    if (count > 0) {
      READ_CHECKPOINT_MALLOC_OR_EXIT(out_data_indices, sizeof(int) * (size_t)count, filename, "out_data_indices");
      FREAD(out_data_indices, sizeof(int), (size_t)count, cp_file, filename, "out_data_indices");
      for (int i = 0; i < count; i++) {
        const int idx = out_data_indices[i];
        if (idx < 0 || idx >= ntot_grid) {
          fprintf(stderr, "read_checkpoint: FATAL: invalid gridpoint index=%d at element=%d for grid=%d (ntot_grid=%d) in %s.\n", idx, i, grid,
                  ntot_grid, filename);
          exit(EXIT_FAILURE);
        } // END IF: serialized grid-point index is invalid
        if (i > 0 && idx <= out_data_indices[i - 1]) {
          fprintf(stderr, "read_checkpoint: FATAL: non-increasing gridpoint index=%d after previous=%d at element=%d for grid=%d in %s.\n", idx,
                  out_data_indices[i - 1], i, grid, filename);
          exit(EXIT_FAILURE);
        } // END IF: serialized grid-point indices are not strictly increasing
      } // END LOOP: for i over serialized grid-point indices
      READ_CHECKPOINT_MALLOC_OR_EXIT(compact_out_data, sizeof(REAL) * compact_count, filename, "compact_out_data");
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
    READ_CHECKPOINT_MALLOC_OR_EXIT(griddata[grid].gridfuncs.y_n_gfs, sizeof(REAL) * numpts_in_all_evol_gfs, filename,
                                   "griddata[grid].gridfuncs.y_n_gfs");
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
} // END FUNCTION: read_checkpoint
