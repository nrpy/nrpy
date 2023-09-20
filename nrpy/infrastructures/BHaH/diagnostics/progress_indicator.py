"""
The BHaH progress indicator, indicating the current
  iteration, time, timestep, benchmark (t/hr), etc.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc
import nrpy.params as par

_ = par.register_CodeParameter(
    "TIMEVAR",
    __name__,
    "start_wallclock_time",
    0.0,
    commondata=True,
    add_to_parfile=False,
    add_to_set_CodeParameters_h=False,
)


def register_CFunction_progress_indicator() -> None:
    """
    Registers a C function that serves as a progress indicator for the simulation.

    The function updates the elapsed time and other metrics to give the user an
    idea of the progress being made in the computation.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "time.h"]
    desc = """Output progress indicator, including elapsed time and other metrics to give the user an
    idea of the status of the the computation."""
    c_type = "void"
    name = "progress_indicator"
    params = """commondata_struct *restrict commondata, const griddata_struct *restrict griddata"""

    body = r"""
  if (commondata->nn == 0) {
    CURRTIME_FUNC(&commondata->start_wallclock_time);
  }
  TIMEVAR currtime;
  CURRTIME_FUNC(&currtime);
  const REAL time_in_ns = TIME_IN_NS(commondata->start_wallclock_time, currtime);
  const REAL seconds_elapsed = time_in_ns / 1e9;
  const REAL phys_time_per_sec = (commondata->time - 0) / seconds_elapsed;
  // const REAL RHS_pt_evals_per_sec = num_RHS_pt_evals / (time_in_ns / 1.0e9);

  // Step 2: Output simulation progress to stderr
  fprintf(stderr,"%c[2K", 27); // Clear the line
  fprintf(stderr, "It: %d t=%.3f / %.3f = %.2f%% dt=1/%.1f | t/h=%.2f\r", commondata->nn,
          commondata->time, commondata->t_final, 100.0*commondata->time / commondata->t_final,
          1.0 / (double)commondata->dt, (double)(phys_time_per_sec * 3600.0));
  fflush(stderr); // Flush the stderr buffer
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
