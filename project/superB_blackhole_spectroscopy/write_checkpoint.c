#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Write a checkpoint file
 */
void write_checkpoint(const commondata_struct *restrict commondata, griddata_struct *restrict griddata) {

  char filename[256];
  snprintf(filename, 256, "checkpoint-conv_factor%.2f.dat", commondata->convergence_factor);
  const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->checkpoint_every;
  // Explanation of the if() below:
  // Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime.
  // Step 2: Multiplying by outevery yields the exact time we should output again, t_out.
  // Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible!
  if (fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt) {
    FILE *cp_file = fopen(filename, "w+");
    fwrite(commondata, sizeof(commondata_struct), 1, cp_file);
    fprintf(stderr, "WRITING CHECKPOINT: cd struct size = %ld time=%e\n", sizeof(commondata_struct), commondata->time);
    for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
      fwrite(&griddata[grid].params, sizeof(params_struct), 1, cp_file);
      const int ntot =
          (griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2);
      // First we free up memory so we can malloc more: copy y_n_gfs to diagnostic_output_gfs & then free y_n_gfs.
#pragma omp parallel for
      for (int i = 0; i < ntot * NUM_EVOL_GFS; i++)
        griddata[grid].gridfuncs.diagnostic_output_gfs[i] = griddata[grid].gridfuncs.y_n_gfs[i];
      MoL_free_memory_y_n_gfs(&griddata[grid].gridfuncs);

      int count = 0;
      const int maskval = 1; // to be replaced with griddata[grid].mask[i].
#pragma omp parallel for reduction(+ : count)
      for (int i = 0; i < ntot; i++) {
        if (maskval >= +0)
          count++;
      }
      fwrite(&count, sizeof(int), 1, cp_file);

      int *restrict out_data_indices = (int *restrict)malloc(sizeof(int) * count);
      REAL *restrict compact_out_data = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);
      int which_el = 0;
      for (int i = 0; i < ntot; i++) {
        if (maskval >= +0) {
          out_data_indices[which_el] = i;
          for (int gf = 0; gf < NUM_EVOL_GFS; gf++)
            compact_out_data[which_el * NUM_EVOL_GFS + gf] = griddata[grid].gridfuncs.diagnostic_output_gfs[ntot * gf + i];
          which_el++;
        }
      }
      // printf("HEY which_el = %d | count = %d\n", which_el, count);
      fwrite(out_data_indices, sizeof(int), count, cp_file);
      fwrite(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);
      free(out_data_indices);
      free(compact_out_data);
      MoL_malloc_y_n_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
#pragma omp parallel for
      for (int i = 0; i < ntot * NUM_EVOL_GFS; i++)
        griddata[grid].gridfuncs.y_n_gfs[i] = griddata[grid].gridfuncs.diagnostic_output_gfs[i];
    }
    fclose(cp_file);
    fprintf(stderr, "FINISHED WRITING CHECKPOINT\n");
  }
}
