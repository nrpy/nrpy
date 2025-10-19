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

/**
 * Write a checkpoint file
 */
void write_checkpoint(const commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_device) {
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

    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
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
} // END FUNCTION write_checkpoint
