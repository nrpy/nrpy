#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "time.h"
/*
 * Output progress indicator, including elapsed time and other metrics to give the user an
 * idea of the status of the the computation.
 */
void progress_indicator(commondata_struct *restrict commondata, const griddata_struct *restrict griddata) {

  if (commondata->nn == commondata->nn_0) {
    CURRTIME_FUNC(&commondata->start_wallclock_time);
  }
  TIMEVAR currtime;
  CURRTIME_FUNC(&currtime);
  const REAL time_in_ns = TIME_IN_NS(commondata->start_wallclock_time, currtime);
  const REAL seconds_elapsed = time_in_ns / 1e9;
  const REAL phys_time_per_sec = (commondata->time - commondata->t_0) / seconds_elapsed;
  // const REAL RHS_pt_evals_per_sec = num_RHS_pt_evals / (time_in_ns / 1.0e9);
  const REAL seconds_remaining = (int)((commondata->t_final - commondata->time) / phys_time_per_sec);
  int time_remaining__hrs = (int)(seconds_remaining / 3600.0);
  int time_remaining__mins = (int)(seconds_remaining / 60.0) % 60;
  int time_remaining__secs = (int)(seconds_remaining)-time_remaining__mins * 60 - time_remaining__hrs * 3600;
  if (commondata->nn == 0) {
    // Just display zeros for ETA at the zeroth iteration
    time_remaining__hrs = 0;
    time_remaining__mins = 0;
    time_remaining__secs = 0;
  }
  // Step 2: Output simulation progress to stderr
  fprintf(stderr, "%c[2K", 27); // Clear the line
  fprintf(stderr, "It: %d t=%.3f / %.1f = %.2f%% dt=1/%.1f | t/h=%.2f ETA %dh%02dm%02ds\r", commondata->nn, commondata->time, commondata->t_final,
          100.0 * commondata->time / commondata->t_final, 1.0 / (double)commondata->dt, (double)(phys_time_per_sec * 3600.0), time_remaining__hrs,
          time_remaining__mins, time_remaining__secs);
  fflush(stderr); // Flush the stderr buffer
}
