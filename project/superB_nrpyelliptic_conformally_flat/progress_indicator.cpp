#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "time.h"

/**
 * Output progress indicator, including elapsed time and other metrics to give the user an
 *     idea of the status of the computation.
 */
void progress_indicator(commondata_struct *restrict commondata, const griddata_struct *restrict griddata) {
  // Proceed only if progress output is enabled (output_progress_every > 0)
  // and the current iteration (nn) is a multiple of the output frequency (output_progress_every)
  if (!(commondata->output_progress_every >= 0 && (commondata->nn % commondata->output_progress_every == 0)))
    return;
  // Output simulation progress:

  fprintf(stderr, "nn / nn_max = %d / %d ; log10(residual) / log10(residual_target) =  %.4f / %.4f \r", commondata->nn, commondata->nn_max,
          commondata->log10_current_residual, commondata->log10_residual_tolerance);
  fflush(stderr); // Flush the stderr buffer
} // END FUNCTION progress_indicator
