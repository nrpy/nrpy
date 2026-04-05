#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

#ifdef __CUDACC__
#define BHAH_CPY_DEVICE_TO_HOST_PARAMS() memcpy(&griddata[grid].params, &griddata_device[grid].params, sizeof(params_struct));
#define BHAH_CPY_DEVICE_TO_HOST_ALL_GFS()                                                                                                            \
  for (int gf = 0; gf < NUM_EVOL_GFS; ++gf) {                                                                                                        \
    size_t stream_id = griddata_device[grid].params.grid_idx % NUM_STREAMS;                                                                          \
    cpyDevicetoHost__gf(commondata, &griddata[grid].params, griddata[grid].gridfuncs.y_n_gfs, griddata_device[grid].gridfuncs.y_n_gfs, gf, gf,       \
                        stream_id);                                                                                                                  \
  }
#endif // __CUDACC__

static inline void BHAH_safe_write_impl(const void *ptr, size_t size, size_t nmemb, FILE *fp, const char *what, const char *file, int line,
                                        const char *func) {
  clearerr(fp);
  errno = 0;

  const size_t wrote = fwrite(ptr, size, nmemb, fp);
  if (wrote != nmemb) {
    const int err = errno;
    fprintf(stderr, "%s:%d: %s: fwrite failed writing %s (wanted %zu, wrote %zu)\n", file, line, func, what, nmemb, wrote);

    if (err != 0) {
      fprintf(stderr, "%s:%d: %s: errno=%d (%s)\n", file, line, func, err, strerror(err));
    } else if (ferror(fp)) {
      fprintf(stderr, "%s:%d: %s: stream error set but errno not set\n", file, line, func);
    } // END IF/ELSE no error

    // Fatal: checkpoint output must not silently continue on partial write.
    exit(1);
  } // END IF wrote != nmemb
} // END FUNCTION BHAH_safe_write_impl()

#define FWRITE(ptr, size, nmemb, fp, what) BHAH_safe_write_impl((ptr), (size_t)(size), (size_t)(nmemb), (fp), (what), __FILE__, __LINE__, __func__)

/**
 * Write a checkpoint file
 */
void write_checkpoint(const commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  char filename[256];
  snprintf(filename, 256, "checkpoint-conv_factor%.2f.dat", commondata->convergence_factor);

  const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->checkpoint_every;
  if (outevery <= (REAL)0.0)
    return; // outevery <= 0 means do not checkpoint.
  // Explanation of the if() below:
  // Step 1: round(currtime / outevery) gives the nearest integer n to the ratio currtime/outevery.
  // Step 2: Multiplying by outevery yields the nearest output time t_out = n * outevery.
  // Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible.
  if (fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt) {
    FILE *cp_file = fopen(filename, "wb");
    if (cp_file == NULL) {
      perror("write_checkpoint: Failed to open checkpoint file. Check permissions and disk space availability.");
      exit(1);
    } // END IF cp_file == NULL
    FWRITE(commondata, sizeof(commondata_struct), 1, cp_file, "commondata");
    fprintf(stderr, "WRITING CHECKPOINT: cd struct size = %zu time=%e\n", sizeof(commondata_struct), commondata->time);

    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
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
        const int maskval = 1;

        if (maskval >= +0)
          count++;
      } // END LOOP over all gridpoints
      FWRITE(&count, sizeof(int), 1, cp_file, "gridpoint_count");

      int *restrict out_data_indices;
      BHAH_MALLOC(out_data_indices, sizeof(int) * count);
      REAL *restrict compact_out_data;
      BHAH_MALLOC(compact_out_data, sizeof(REAL) * NUM_EVOL_GFS * count);
      int which_el = 0;

      for (int i = 0; i < ntot_grid; i++) {
        const int maskval = 1;

        if (maskval >= +0) {
          out_data_indices[which_el] = i;
          for (int gf = 0; gf < NUM_EVOL_GFS; gf++)
            compact_out_data[which_el * NUM_EVOL_GFS + gf] = griddata[grid].gridfuncs.y_n_gfs[ntot_grid * gf + i];
          which_el++;
        } // END IF maskval >= +0
      } // END LOOP over all gridpoints

      FWRITE(out_data_indices, sizeof(int), count, cp_file, "out_data_indices");
      FWRITE(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file, "compact_out_data");
      free(out_data_indices);
      free(compact_out_data);

      // Re-allocate memory for intermediate-stage scratch storage for MoL.
      MoL_malloc_intermediate_stage_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
    } // END LOOP over grids
    fclose(cp_file);
    fprintf(stderr, "FINISHED WRITING CHECKPOINT\n");
  } // END IF ready to write checkpoint
} // END FUNCTION write_checkpoint
