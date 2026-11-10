"""
Check the 10M reference and inspiral-to-ringdown matching in SEBOBv2.

This script runs a fixed set of spin configurations and compares the
recovered orbital and spin reference with a separate numerical calculation.
It also checks that attachment times, ringing frequencies, and damping
times remain consistent, and that the waveform samples are finite and
ordered in time. An optional baseline executable checks that the
aligned-spin waveform output remains unchanged.

The production calculation uses the eighth-order Runge-Kutta
Prince-Dormand method (RK8PD, an embedded 8/9 pair). The reference calculation
uses the fourth/fifth-order Runge-Kutta-Fehlberg pair (RKF45). These checks
assess numerical consistency within the existing model; accuracy against
numerical-relativity waveforms requires separate validation.

Use a project generated with ``--projected-attachment-diagnostics``.
Each external program has an execution-time limit and an output-file-size
limit. Temporary inputs and outputs are removed even when a check fails.

Author: Suchindram Dasgupta
        sd00113 **at** mix **dot** wvu **dot** edu
"""

import bisect
import cmath
import math
import os
import signal
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

NumericRecord = Dict[str, float]
Waveform = List[Tuple[float, complex]]
FineSamples = List[Tuple[float, complex, float]]
FineDynamics = List[Tuple[float, complex, float, float, float]]
CASE_DATA = (
    ("default", 1.0, (0.0, 0.1, 0.4), (0.0, -0.2, -0.3), 0.01118),
    ("antialigned", 1.0, (0.0, 0.0, -0.8), (0.0, 0.0, -0.8), 0.01118),
    (
        "tilted",
        1.0,
        (0.4749999999999999, 0.0, -0.8227241335952168),
        (0.0, 0.4749999999999999, -0.8227241335952168),
        0.01118,
    ),
    ("unequal", 2.0, (0.1, 0.0, -0.8), (0.0, 0.1, -0.8), 0.01118),
    ("short", 1.0, (0.0, 0.1, 0.4), (0.0, -0.2, -0.3), 0.04),
    ("zero", 1.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.01118),
    ("equalspin", 1.0, (0.0, 0.0, 0.4), (0.0, 0.0, 0.4), 0.01118),
    ("q20", 20.0, (0.0, 0.8, 0.0), (0.0, -0.6, -0.3), 0.01118),
    ("short_allfine", 1.0, (0.0, 0.1, 0.4), (0.0, -0.2, -0.3), 0.05),
)

# Integration-only numerical reference. It links unchanged generated equation
# functions, never the production reference routine, and stops PN time evolution
# at the same orbital-frequency index used by the decoupled model.
REFERENCE_SOURCE = r"""
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include <float.h>

typedef struct {
  gsl_odeiv2_driver *driver;
  size_t dimension;
  REAL initial[NUMVARS_SPIN];
} independent_evolution; // END STRUCT: independent evolution data

/**
 * Evolve a saved state with independent RKF45 tolerances.
 * @param[in,out] data Solver and saved initial state.
 * @param duration Signed integration interval.
 * @param[out] state Recovered state.
 */
static void independent_state(independent_evolution *data, REAL duration, REAL *state) {
  memcpy(state, data->initial, data->dimension * sizeof(REAL));
  if (duration != 0.0) {
    int status = gsl_odeiv2_driver_reset_hstart(data->driver, copysign(0.01, duration));
    REAL time = 0.0;
    if (status == GSL_SUCCESS)
      status = gsl_odeiv2_driver_apply(data->driver, &time, duration, state);
    if (status != GSL_SUCCESS) {
      fprintf(stderr, "Independent evolution failed: %s\n", gsl_strerror(status));
      exit(EXIT_FAILURE);
    } // END IF: independent evolution failed
  } // END IF: nonzero independent evolution interval
  for (size_t i = 0; i < data->dimension; i++) {
    if (!isfinite(state[i])) {
      fprintf(stderr, "Nonfinite independent state component %zu\n", i);
      exit(EXIT_FAILURE);
    } // END IF: independent state is nonfinite
  } // END LOOP: independent state finite checks
} // END FUNCTION: evolve independent reference state

/**
 * Evaluate the shared orbital equations and reject nonfinite reference rates.
 * @param time Orbital evolution time.
 * @param[in] state Orbital state.
 * @param[out] rates Orbital derivatives.
 * @param[in,out] parameters Initialized model parameters.
 * @return GSL_SUCCESS, GSL_EBADFUNC for nonfinite rates, or the equation status.
 */
static int independent_orbit_rhs(REAL time, const REAL *state, REAL *rates, void *parameters) {
  const int status = SEOBNRv5_aligned_spin_right_hand_sides(time, state, rates, parameters);
  if (status != GSL_SUCCESS)
    return status;
  for (size_t i = 0; i < 4; i++)
    if (!isfinite(rates[i]))
      return GSL_EBADFUNC;
  return GSL_SUCCESS;
} // END FUNCTION: check independent orbital derivatives

/**
 * Locate a component value by independent bisection.
 * @param[in,out] data Solver and saved state.
 * @param component State component defining the event.
 * @param target Required physical value.
 * @param low Lower local-time endpoint.
 * @param high Upper local-time endpoint.
 * @return Local time of the recovered event.
 */
static REAL bisect_event(independent_evolution *data, size_t component, REAL target, REAL low, REAL high) {
  REAL state[NUMVARS_SPIN];
  independent_state(data, low, state);
  REAL f_low = state[component] - target;
  independent_state(data, high, state);
  REAL f_high = state[component] - target;
  if (f_low == 0.0)
    return low;
  if (f_high == 0.0)
    return high;
  if (f_low * f_high > 0.0) {
    fprintf(stderr, "Independent event is not bracketed\n");
    exit(EXIT_FAILURE);
  } // END IF: independent event not bracketed
  for (size_t iteration = 0; iteration < 80; iteration++) {
    const REAL midpoint = low + 0.5 * (high - low);
    independent_state(data, midpoint, state);
    const REAL value = state[component] - target;
    if (value == 0.0 || high - low < 2e-10)
      return midpoint;
    if (f_low * value <= 0.0) {
      high = midpoint;
    } // END IF: root in lower half
    else {
      low = midpoint;
      f_low = value;
    } // END ELSE: root in upper half
  } // END LOOP: independent event bisection
  fprintf(stderr, "Independent event did not converge\n");
  exit(EXIT_FAILURE);
} // END FUNCTION: bisect independent event

typedef struct {
  commondata_struct *parameters;
  REAL min_omega_dot;
} independent_spin_parameters; // END STRUCT: independent PN observations

static int independent_spin_rhs(REAL time, const REAL *state, REAL *rates, void *parameters) {
  independent_spin_parameters *data = (independent_spin_parameters *)parameters;
  const int status = SEOBNRv5_quasi_precessing_spin_equations(time, state, rates, data->parameters);
  if (status != GSL_SUCCESS)
    return status;
  for (size_t i = 0; i < NUMVARS_SPIN; i++)
    if (!isfinite(rates[i]))
      return GSL_EBADFUNC;
  data->min_omega_dot = fmin(data->min_omega_dot, rates[OMEGA_PN]);
  return status;
} // END FUNCTION: observe time-domain PN evolution

/**
 * Recover a reference using independent time-domain evolution.
 * @param argc Argument count.
 * @param[in] argv Parameter-file argument.
 * @return EXIT_SUCCESS on recovery, EXIT_FAILURE on input or evolution failure.
 */
int main(int argc, const char *argv[]) {
  if (argc != 2)
    return EXIT_FAILURE;
  commondata_struct parameters;
  commondata_struct_set_to_default(&parameters);
  cmdline_input_and_parfile_parser(&parameters, argc, argv);
  parameters.chi1 = parameters.chi1_z;
  parameters.chi2 = parameters.chi2_z;
  SEOBNRv5_quasi_precessing_spin_coefficients(&parameters);
  SEOBNRv5_quasi_precessing_spin_dynamics(&parameters);
  SEOBNRv5_aligned_spin_initial_conditions_conservative(&parameters);
  SEOBNRv5_aligned_spin_pa_integration(&parameters);

  // Recover the raw time origin from a saved interface sample, independently
  // of the newly registered metadata used by production.
  REAL origin = 0.0;
  int found_origin = 0;
  for (size_t segment = 0; segment < 2; segment++) {
    const REAL *trajectory = segment == 0 ? parameters.dynamics_low : parameters.dynamics_fine;
    const size_t count = segment == 0 ? parameters.nsteps_low : parameters.nsteps_fine;
    for (size_t i = 0; i < count; i++) {
      int identical = 1;
      for (size_t component = R; component <= PPHI; component++)
        identical = identical && trajectory[IDX(i,component)] == parameters.dynamics_raw[IDX(0,component)];
      if (identical) {
        origin = trajectory[IDX(i,TIME)];
        found_origin = 1;
      } // END IF: raw interface state recovered
    } // END LOOP: saved interface candidates
  } // END LOOP: low and fine trajectories
  if (!found_origin) {
    fprintf(stderr, "Independent PA/ODE time origin not observed\n");
    return EXIT_FAILURE;
  } // END IF: time origin not observed

  size_t left = 0;
  const REAL initial_radius = parameters.dynamics_raw[IDX(0,R)];
  if (initial_radius > 10.0) {
    while (left + 1 < parameters.nsteps_raw && parameters.dynamics_raw[IDX(left + 1,R)] > 10.0)
      left++;
    if (left + 1 == parameters.nsteps_raw)
      return EXIT_FAILURE;
    if (parameters.dynamics_raw[IDX(left + 1,R)] == 10.0)
      left++;
  } // END IF: raw inspiral crosses 10M
  gsl_odeiv2_system orbit_system = {independent_orbit_rhs, NULL, 4, &parameters};
  independent_evolution orbit;
  orbit.dimension = 4;
  orbit.driver = gsl_odeiv2_driver_alloc_y_new(&orbit_system, gsl_odeiv2_step_rkf45, 0.01, 1e-14, 1e-13);
  if (orbit.driver == NULL)
    return EXIT_FAILURE;
  gsl_odeiv2_driver_set_nmax(orbit.driver, 100000);
  for (size_t i = 0; i < 4; i++)
    orbit.initial[i] = parameters.dynamics_raw[IDX(left,i + 1)];
  REAL low = 0.0, high = 0.0, duration = 0.0;
  REAL reference[NUMVARS_SPIN];
  REAL original[4];
  memcpy(original, orbit.initial, sizeof(original));
  REAL base_time = parameters.dynamics_raw[IDX(left,TIME)];
  if (orbit.initial[0] < 10.0) {
    REAL elapsed = 0.0;
    for (size_t iteration = 0; iteration < 100000; iteration++) {
      independent_state(&orbit, -1.0, reference);
      if (reference[0] >= 10.0)
        break;
      memcpy(orbit.initial, reference, 4 * sizeof(REAL));
      elapsed -= 1.0;
    } // END LOOP: earlier orbital reference interval
    base_time += elapsed;
    low = -1.0;
  } // END IF: initialized orbit starts inside 10M
  else if (orbit.initial[0] > 10.0) {
    high = parameters.dynamics_raw[IDX(left + 1,TIME)] - parameters.dynamics_raw[IDX(left,TIME)];
  } // END ELSE IF: reference follows raw sample
  if (orbit.initial[0] != 10.0)
    duration = bisect_event(&orbit, 0, 10.0, low, high);
  independent_state(&orbit, duration, reference);
  REAL rates[4];
  if (independent_orbit_rhs(0.0, reference, rates, &parameters) != GSL_SUCCESS)
    return EXIT_FAILURE;
  const REAL omega_reference = rates[1];
  REAL continued[4];
  memcpy(orbit.initial, original, sizeof(original));
  const REAL recovery_interval = initial_radius < 10.0 ? base_time + duration : duration;
  REAL maximum_rdot = -DBL_MAX;
  for (size_t i = 0; i <= 32; i++) {
    REAL sampled[NUMVARS_SPIN], derivative[4];
    independent_state(&orbit, recovery_interval * i / 32.0, sampled);
    if (independent_orbit_rhs(0.0, sampled, derivative, &parameters) != GSL_SUCCESS)
      return EXIT_FAILURE;
    maximum_rdot = fmax(maximum_rdot, derivative[0]);
  } // END LOOP: reference continuation samples
  memcpy(orbit.initial, reference, sizeof(original));
  independent_state(&orbit, -recovery_interval, continued);
  REAL roundtrip = 0.0;
  for (size_t i = 0; i < 4; i++)
    roundtrip = fmax(roundtrip, fabs(continued[i] - original[i]));
  gsl_odeiv2_driver_free(orbit.driver);

  // Integrate the untransformed PN time equations one adaptive step at a time.
  // This stops before merger-side PN evolution, then bisects only the final
  // step that encloses the specified orbital-frequency index.
  REAL spin[NUMVARS_SPIN] = {0,0,1,parameters.chi1_x,parameters.chi1_y,parameters.chi1_z,
                           parameters.chi2_x,parameters.chi2_y,parameters.chi2_z,parameters.initial_omega};
  independent_spin_parameters spin_parameters = {&parameters, DBL_MAX};
  gsl_odeiv2_system spin_system = {independent_spin_rhs, NULL, NUMVARS_SPIN, &spin_parameters};
  gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, NUMVARS_SPIN);
  gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(1e-14, 1e-13);
  gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(NUMVARS_SPIN);
  if (step == NULL || control == NULL || evolve == NULL)
    return EXIT_FAILURE;
  const REAL direction = copysign(1.0, omega_reference - parameters.initial_omega);
  REAL time = 0.0, previous_time = 0.0, h = direction * 0.01;
  REAL previous[NUMVARS_SPIN];
  memcpy(previous, spin, sizeof(spin));
  size_t iterations = 0;
  while (direction * (spin[OMEGA_PN] - omega_reference) < 0.0 && iterations < 100000) {
    memcpy(previous, spin, sizeof(spin));
    previous_time = time;
    const int status = gsl_odeiv2_evolve_apply(evolve, control, step, &spin_system, &time, direction * 2e9, &h, spin);
    if (status != GSL_SUCCESS) {
      fprintf(stderr, "Independent PN evolution failed: %s\n", gsl_strerror(status));
      return EXIT_FAILURE;
    } // END IF: time-domain spin evolution failed
    iterations++;
  } // END WHILE: approach reference spin frequency
  gsl_odeiv2_evolve_free(evolve);
  gsl_odeiv2_control_free(control);
  gsl_odeiv2_step_free(step);
  independent_evolution spin_evolution;
  spin_evolution.dimension = NUMVARS_SPIN;
  memcpy(spin_evolution.initial, previous, sizeof(previous));
  spin_evolution.driver = gsl_odeiv2_driver_alloc_y_new(&spin_system, gsl_odeiv2_step_rkf45, direction * 0.01, 1e-14, 1e-13);
  if (spin_evolution.driver == NULL)
    return EXIT_FAILURE;
  const REAL elapsed = time - previous_time;
  const REAL spin_duration = bisect_event(&spin_evolution, OMEGA_PN, omega_reference, fmin(0.0, elapsed), fmax(0.0, elapsed));
  independent_state(&spin_evolution, spin_duration, spin);
  gsl_odeiv2_driver_free(spin_evolution.driver);
  const REAL ln_norm = sqrt(spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2]);
  for (size_t i = 0; i < 3; i++)
    spin[i] /= ln_norm;
  printf("NRPY_REFERENCE r=%.17g t=%.17g omega=%.17g raw_origin=%.17g phi=%.17g prstar=%.17g pphi=%.17g roundtrip=%.17g maximum_rdot=%.17g reference_rdot=%.17g min_omega_dot=%.17g ln_norm=%.17g raw_start_r=%.17g\n",
      reference[0], origin + base_time + duration, omega_reference, origin,
      reference[1], reference[2], reference[3], roundtrip, maximum_rdot, rates[0], spin_parameters.min_omega_dot, ln_norm, initial_radius);
  printf("NRPY_SPINS lnx=%.17g lny=%.17g lnz=%.17g c1x=%.17g c1y=%.17g c1z=%.17g c2x=%.17g c2y=%.17g c2z=%.17g p1=%.17g p2=%.17g\n",
      spin[0],spin[1],spin[2],spin[3],spin[4],spin[5],spin[6],spin[7],spin[8],
      spin[0]*spin[3]+spin[1]*spin[4]+spin[2]*spin[5],spin[0]*spin[6]+spin[1]*spin[7]+spin[2]*spin[8]);
  spline_data *spin_splines[] = {
      &parameters.chi1_lnhat, &parameters.chi2_lnhat, &parameters.chi1_l, &parameters.chi2_l,
      &parameters.chi1_x_spline, &parameters.chi1_y_spline, &parameters.chi1_z_spline,
      &parameters.chi2_x_spline, &parameters.chi2_y_spline, &parameters.chi2_z_spline,
      &parameters.lnhat_x, &parameters.lnhat_y, &parameters.lnhat_z,
      &parameters.L_x, &parameters.L_y, &parameters.L_z};
  for (size_t i = 0; i < sizeof(spin_splines) / sizeof(spin_splines[0]); i++) {
    if (spin_splines[i]->spline != NULL)
      gsl_spline_free(spin_splines[i]->spline);
    if (spin_splines[i]->acc != NULL)
      gsl_interp_accel_free(spin_splines[i]->acc);
  } // END LOOP: free owned reference spin splines
  free(parameters.dynamics_low);
  free(parameters.dynamics_fine);
  free(parameters.dynamics_raw);
  return EXIT_SUCCESS;
} // END FUNCTION: independent numerical reference
"""


def run_checked(
    arguments: Sequence[str],
    cwd: Path,
    stdout: Path,
    stderr: Path,
    expected_status: int = 0,
) -> None:
    """
    Run a bounded generated-product or compiler subprocess.

    :param arguments: Executable and argument vector.
    :param cwd: Owned subprocess working directory.
    :param stdout: Owned stdout file.
    :param stderr: Owned stderr file.
    :param expected_status: Required subprocess return status.
    :raises RuntimeError: If the subprocess fails or exceeds its execution limit.
    """
    environment = dict(os.environ, OMP_NUM_THREADS="1")
    # prlimit applies the cap before exec without Python's thread-unsafe
    # preexec_fn hook. Compiler children inherit the same file-size limit.
    command = ["prlimit", "--fsize=67108864:67108864", "--", *arguments]
    with stdout.open("wb") as out, stderr.open("wb") as err:
        with subprocess.Popen(
            command,
            cwd=cwd,
            env=environment,
            stdout=out,
            stderr=err,
            shell=False,
            start_new_session=True,
        ) as process:
            try:
                returncode = process.wait(timeout=180)
            except subprocess.TimeoutExpired as exc:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
                raise RuntimeError(
                    f"Subprocess exceeded 180 seconds: {arguments[0]}"
                ) from exc
    if returncode != expected_status:
        with stderr.open("rb") as stream:
            stream.seek(max(0, stderr.stat().st_size - 4096))
            detail = stream.read(4096).decode("utf-8", errors="replace")
        raise RuntimeError(f"{arguments[0]} returned {returncode}: {detail}")


def parse_records(
    path: Path,
) -> Tuple[Dict[str, NumericRecord], FineSamples, FineDynamics]:
    """
    Parse finite named observations and corrected fine samples.

    :param path: Generated-executable diagnostic output.
    :return: Named observations, corrected fine waveform, and original fine dynamics.
    :raises ValueError: If a record is duplicated, malformed, or nonfinite.
    """
    records: Dict[str, NumericRecord] = {}
    corrected: FineSamples = []
    original: FineDynamics = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.startswith("NRPY_"):
                continue
            label, *items = line.split()
            values = {
                key: float(value)
                for key, value in (item.split("=", 1) for item in items)
            }
            if not all(math.isfinite(value) for value in values.values()):
                raise ValueError(f"Nonfinite {label} observation")
            if label == "NRPY_FINE":
                original.append(
                    (
                        values["t"],
                        complex(values["re"], values["im"]),
                        values["r"],
                        values["prstar"],
                        values["omega"],
                    )
                )
                continue
            if label == "NRPY_CORRECTED":
                corrected.append(
                    (values["t"], complex(values["re"], values["im"]), values["phi"])
                )
            elif label in records:
                raise ValueError(f"Duplicate {label} observation")
            else:
                records[label] = values
    return records, corrected, original


def read_waveform(path: Path) -> Waveform:
    """
    Read the complete finite waveform and require strictly increasing times.

    :param path: Generated waveform stdout.
    :return: Time and complex strain samples.
    :raises ValueError: If output is empty, malformed, nonfinite, or unordered.
    """
    samples: Waveform = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            values = [float(item) for item in line.split()]
            if len(values) != 3 or not all(math.isfinite(value) for value in values):
                raise ValueError("Malformed or nonfinite waveform row")
            time, real, imaginary = values
            if samples and time <= samples[-1][0]:
                raise ValueError("Waveform times are not strictly increasing")
            samples.append((time, complex(real, imaginary)))
    if not samples:
        raise ValueError("Generated waveform is empty")
    return samples


def check_handoff(
    before: float,
    after: float,
    omega: float,
    omega22: float,
    checked_omega: float,
    tau: float,
    tau22: float,
    checked_tau: float,
) -> None:
    """
    Require the selected attachment time and QNM aliases to remain identical.

    :param before: Time selected by special-amplitude coefficients.
    :param after: Time retained by NQC.
    :param omega: Generic QNM frequency.
    :param omega22: Mode-specific frequency.
    :param checked_omega: Frequency returned independently by the shared helper.
    :param tau: Generic damping time.
    :param tau22: Mode-specific damping time.
    :param checked_tau: Damping time returned by the shared helper.
    :raises ValueError: If a handoff or QNM equality fails.

    Doctests:
    >>> check_handoff(3.0, 3.1, 0.4, 0.4, 0.4, 12.0, 12.0, 12.0)
    Traceback (most recent call last):
    ...
    ValueError: NQC changed the selected attachment time
    >>> check_handoff(3.0, 3.0, 0.4, 0.5, 0.4, 12.0, 12.0, 12.0)
    Traceback (most recent call last):
    ...
    ValueError: QNM fields disagree with the shared helper
    """
    if before != after:
        raise ValueError("NQC changed the selected attachment time")
    if not omega == omega22 == checked_omega or not tau == tau22 == checked_tau:
        raise ValueError("QNM fields disagree with the shared helper")


def cubic_values(
    times: Sequence[float], values: Sequence[float], time: float
) -> Tuple[float, float, float]:
    """
    Independently evaluate a natural cubic spline and its first two derivatives.

    :param times: Ordered sample times.
    :param values: Sampled amplitude or unwrapped phase.
    :param time: Evaluation time within the sample interval.
    :return: Value, first derivative, and second derivative.
    :raises ValueError: If the interpolation data or requested time is invalid.

    Doctests:
    >>> cubic_values((0.0, 1.0, 2.0), (0.0, 1.0, 2.0), 0.5)
    (0.5, 1.0, 0.0)
    >>> cubic_values((0.0, 0.01, 0.02), (1e307, 2e307, 3e307), 0.015)
    Traceback (most recent call last):
    ...
    ValueError: Nonfinite cubic-spline value or derivative
    """
    count = len(times)
    if count < 3 or count != len(values) or not times[0] <= time <= times[-1]:
        raise ValueError("Invalid cubic-spline data or evaluation time")
    spacing = [times[i + 1] - times[i] for i in range(count - 1)]
    if min(spacing) <= 0:
        raise ValueError("Cubic-spline times must increase")
    diagonal = [1.0] * count
    rhs = [0.0] * count
    upper = [0.0] * count
    for i in range(1, count - 1):
        lower = spacing[i - 1]
        diagonal[i] = 2 * (lower + spacing[i])
        upper[i] = spacing[i]
        rhs[i] = 6 * (
            (values[i + 1] - values[i]) / spacing[i]
            - (values[i] - values[i - 1]) / lower
        )
        factor = lower / diagonal[i - 1]
        diagonal[i] -= factor * upper[i - 1]
        rhs[i] -= factor * rhs[i - 1]
    second = [0.0] * count
    for i in range(count - 2, 0, -1):
        second[i] = (rhs[i] - upper[i] * second[i + 1]) / diagonal[i]
    index = min(bisect.bisect_right(times, time) - 1, count - 2)
    width = spacing[index]
    right = (time - times[index]) / width
    left = 1 - right
    value = left * values[index] + right * values[index + 1]
    value += (
        ((left**3 - left) * second[index] + (right**3 - right) * second[index + 1])
        * width**2
        / 6
    )
    first = (values[index + 1] - values[index]) / width
    first += (
        ((1 - 3 * left**2) * second[index] + (3 * right**2 - 1) * second[index + 1])
        * width
        / 6
    )
    result = value, first, left * second[index] + right * second[index + 1]
    if not all(math.isfinite(component) for component in result):
        raise ValueError("Nonfinite cubic-spline value or derivative")
    return result


def check_matching(
    samples: FineSamples,
    original: FineDynamics,
    attachment: float,
    targets: NumericRecord,
) -> List[float]:
    """
    Compare the corrected attachment derivatives with the BOB NQC targets.

    :param samples: Corrected fine waveform before IMR resampling.
    :param original: Original fine strain and dynamics used in the NQC basis.
    :param attachment: Common attachment time.
    :param targets: Amplitude and waveform-frequency targets.
    :return: Scaled residuals of the five NQC matching equations.
    :raises ValueError: If samples fail to cover attachment or matching fails.
    """
    times = [sample[0] for sample in samples]
    peak = min(range(len(times)), key=lambda i: abs(times[i] - attachment))
    selected = samples[max(peak, 5) - 5 : min(peak + 5, len(samples))]
    cropped_times = [sample[0] for sample in selected]
    amplitudes = [abs(sample[1]) for sample in selected]
    phases = [cmath.phase(sample[1]) for sample in selected]
    for i in range(1, len(phases)):
        phases[i] = phases[i - 1] + math.remainder(
            phases[i] - phases[i - 1], 2 * math.pi
        )
    amp = cubic_values(cropped_times, amplitudes, attachment)
    phase = cubic_values(cropped_times, phases, attachment)
    scale = targets["A"]
    frequency_scale = targets["w"]
    actual = (amp[0], amp[1], amp[2], -phase[1], -phase[2])
    expected = tuple(targets[key] for key in ("A", "Adot", "Addot", "w", "wdot"))
    scales = (
        scale,
        scale * frequency_scale,
        scale * frequency_scale**2,
        frequency_scale,
        frequency_scale**2,
    )
    if not all(math.isfinite(value) and value > 0 for value in scales):
        raise ValueError("Invalid NQC matching normalization")
    residuals = [
        abs(value - target) / normalization
        for value, target, normalization in zip(actual, expected, scales)
    ]
    if not all(math.isfinite(value) for value in residuals):
        raise ValueError("Nonfinite NQC matching residual")
    raw = original[max(peak, 5) - 5 : min(peak + 5, len(original))]
    amplitude_condition: List[float] = []
    phase_condition: List[float] = []
    for _, strain, radius, prstar, omega in raw:
        q1 = prstar**2 / (radius**2 * omega**2)
        q2, q3 = q1 / radius, q1 / radius**1.5
        p1 = -prstar / (radius * omega)
        p2 = -p1 * prstar**2
        amplitude_condition.append(
            abs(strain)
            * (
                1
                + abs(targets["a1"] * q1)
                + abs(targets["a2"] * q2)
                + abs(targets["a3"] * q3)
            )
        )
        phase_condition.append(
            max(2 * math.pi, max(abs(value) for value in phases))
            + abs(targets["b1"] * p1)
            + abs(targets["b2"] * p2)
        )
    # Natural-cubic interpolation is linear in its sample values. Its basis
    # weights quantify amplification of floating-point cancellation in the NQC
    # coefficient sums and their first/second derivatives. The factor 64 covers
    # the small LU systems, basis evaluation, and complex correction arithmetic.
    weights = [
        cubic_values(
            cropped_times, [float(i == j) for i in range(len(raw))], attachment
        )
        for j in range(len(raw))
    ]
    bounds = []
    for index, normalization in enumerate(scales):
        derivative = index if index < 3 else index - 2
        condition = amplitude_condition if index < 3 else phase_condition
        amplification = sum(
            abs(weight[derivative]) * value for weight, value in zip(weights, condition)
        )
        if not math.isfinite(amplification):
            raise ValueError("Nonfinite NQC arithmetic error estimate")
        bounds.append(
            max(1e-11, 64 * sys.float_info.epsilon * amplification / normalization)
        )
    if not all(math.isfinite(value) for value in bounds) or max(bounds) > 1e-6:
        raise ValueError(
            f"NQC system requires excessive arithmetic error allowance: {bounds}"
        )
    if any(residual > bound for residual, bound in zip(residuals, bounds)):
        raise ValueError(
            f"NQC matching residuals exceed their numerical bound: {residuals}"
        )
    return residuals


def check_phase_overlap(
    waveform: Waveform, corrected: FineSamples, attachment: float, split: int
) -> float:
    """
    Check IMR phase alignment against inspiral at the identical overlap time.

    :param waveform: Complete assembled IMR output.
    :param corrected: Corrected fine inspiral before IMR resampling.
    :param attachment: Time used to shift the printed waveform.
    :param split: First ringdown output index observed from IMR assembly.
    :return: Wrapped phase-alignment residual.
    :raises ValueError: If the split is invalid or phase alignment fails.
    """
    if not 0 < split < len(waveform):
        raise ValueError("Invalid IMR split")
    time, strain = waveform[split]
    times = [sample[0] for sample in corrected]
    phases = [sample[2] for sample in corrected]
    envelopes = [sample[1] * cmath.exp(2j * sample[2]) for sample in corrected]
    phase = cubic_values(times, phases, time + attachment)[0]
    real = cubic_values(times, [value.real for value in envelopes], time + attachment)[
        0
    ]
    imaginary = cubic_values(
        times, [value.imag for value in envelopes], time + attachment
    )[0]
    expected = cmath.exp(-2j * phase) * complex(real, imaginary)
    if not math.isfinite(expected.real) or not math.isfinite(expected.imag):
        raise ValueError("Nonfinite IMR overlap strain")
    residual = abs(
        math.remainder(cmath.phase(strain) - cmath.phase(expected), 2 * math.pi)
    )
    # The product resamples demodulated complex strain, not phase itself.
    # This bound covers roundoff of the orbital phase and 15-digit stdout.
    bound = 1024 * sys.float_info.epsilon * max(1.0, abs(2 * phase)) + 5e-14
    if not math.isfinite(residual) or not math.isfinite(bound) or residual > bound:
        raise ValueError(f"IMR phase overlap failed: {residual}")
    return residual


def main(project: Path, baseline: Optional[Path]) -> None:
    """
    Compile the independent reference and check deterministic projected inputs.

    :param project: Generated diagnostic SEBOBv2 project.
    :param baseline: Optional unchanged scalar-control project.
    :raises RuntimeError: If numerical reference, generated output, or a control fails.
    """
    project = project.resolve()
    if baseline is not None:
        baseline = baseline.resolve()
    with tempfile.TemporaryDirectory(prefix="sebobv2-projected-") as directory:
        work = Path(directory)
        source = work / "independent_reference.c"
        source.write_text(REFERENCE_SOURCE, encoding="utf-8")
        independent = work / "independent_reference"
        objects = sorted(
            str(path) for path in project.rglob("*.o") if path.name != "main.o"
        )
        run_checked(
            [
                os.environ.get("CC", "gcc"),
                "-std=gnu99",
                "-O2",
                "-Wall",
                "-I",
                str(project),
                str(source),
                *objects,
                "-lgsl",
                "-lgslcblas",
                "-lm",
                "-fopenmp",
                "-o",
                str(independent),
            ],
            work,
            work / "compile.out",
            work / "compile.err",
        )
        for name, q, spin1, spin2, omega in CASE_DATA:
            fields: Dict[str, Union[float, bool]] = {
                "mass_ratio": q,
                "chi1": spin1[2],
                "chi2": spin2[2],
                "initial_omega": omega,
                "total_mass": 50.0,
                "dt": 2.4627455127717882e-05,
                "use_projected_attachment": True,
            }
            for body, spin in ((1, spin1), (2, spin2)):
                for axis, value in zip("xyz", spin):
                    fields[f"chi{body}_{axis}"] = value
            parameter_file = work / f"{name}.par"
            parameter_file.write_text(
                "".join(f"{key} = {value}\n" for key, value in fields.items()),
                encoding="utf-8",
            )
            output, diagnostic = work / f"{name}.out", work / f"{name}.err"
            run_checked(
                [str(project / "sebobv2"), str(parameter_file)],
                work,
                output,
                diagnostic,
            )
            waveform = read_waveform(output)
            records, corrected, original = parse_records(diagnostic)
            reference_output = work / f"{name}.reference"
            run_checked(
                [str(independent), str(parameter_file)],
                work,
                reference_output,
                work / "reference.err",
            )
            expected, _, _ = parse_records(reference_output)
            observed_ref, independent_ref = (
                records["NRPY_REFERENCE"],
                expected["NRPY_REFERENCE"],
            )
            for key, tolerance in (
                ("r", 2e-9),
                ("t", 2e-7),
                ("omega", 2e-10),
                ("raw_origin", 1e-10),
                ("raw_start_r", 1e-10),
                ("prstar", 2e-10),
                ("pphi", 2e-10),
                ("phi", 2e-9),
            ):
                if abs(observed_ref[key] - independent_ref[key]) > tolerance:
                    raise RuntimeError(f"{name}: independent {key} disagreement")
            if independent_ref["raw_start_r"] == 10.0:
                if (
                    observed_ref["t"] != observed_ref["raw_origin"]
                    or observed_ref["r"] != 10.0
                ):
                    raise RuntimeError(
                        f"{name}: exact 10M interface state was replaced"
                    )
            if (
                independent_ref["reference_rdot"] >= 0
                or independent_ref["roundtrip"] > 2e-9
            ):
                raise RuntimeError(
                    f"{name}: reference continuation failed its inward/recovery criterion: {independent_ref}"
                )
            observed_spins, independent_spins = (
                records["NRPY_SPINS"],
                expected["NRPY_SPINS"],
            )
            spin_error = max(
                abs(observed_spins[key] - value)
                for key, value in independent_spins.items()
            )
            # The forward calculation intentionally retains the shipped cubic
            # sampler (ODE atol=1e-12, rtol=1e-11). Its interpolation error is
            # distinct from earlier direct evolution. Refining those sampling
            # tolerances establishes convergence toward the time-domain result.
            needs_earlier_spin = observed_ref["omega"] < omega
            spin_bound = 5e-10 if needs_earlier_spin else 5e-6
            if spin_error > spin_bound or independent_ref["min_omega_dot"] <= 0:
                raise RuntimeError(
                    f"{name}: independent spin reference disagreement {spin_error}"
                )
            if observed_ref["initial_omega"] != omega:
                raise RuntimeError(
                    f"{name}: requested spin-reference frequency changed"
                )
            for body, initial_spin in ((1, spin1), (2, spin2)):
                recovered_norm = math.sqrt(
                    sum(observed_spins[f"c{body}{axis}"] ** 2 for axis in "xyz")
                )
                initial_norm = math.sqrt(sum(value**2 for value in initial_spin))
                if abs(recovered_norm - initial_norm) > 3 * spin_bound:
                    raise RuntimeError(f"{name}: reference spin magnitude changed")
            if needs_earlier_spin and (
                observed_ref["t"] >= 0 or observed_ref["min_omega_dot"] <= 0
            ):
                raise RuntimeError(
                    "Short-start reference did not preserve the earlier spin state"
                )
            nqc = records["NRPY_NQC"]
            if name in ("zero", "equalspin") and any(
                records["NRPY_ATTACH"][key] != 0 for key in ("c21", "c43", "c55")
            ):
                raise RuntimeError(f"{name}: equal-spin coefficient defaults changed")
            check_handoff(
                records["NRPY_ATTACH"]["before"],
                nqc["after"],
                nqc["omega"],
                nqc["omega22"],
                nqc["checked_omega"],
                nqc["tau"],
                nqc["tau22"],
                nqc["checked_tau"],
            )
            residuals = check_matching(corrected, original, nqc["after"], nqc)
            phase_residual = check_phase_overlap(
                waveform, corrected, nqc["after"], int(records["NRPY_IMR"]["split"])
            )
            if not waveform[0][0] < 0 < waveform[-1][0]:
                raise RuntimeError(f"{name}: IMR output does not span attachment")
            start_bound = 32 * sys.float_info.epsilon * max(1.0, abs(nqc["after"]))
            if abs(waveform[0][0] + nqc["after"]) > start_bound:
                raise RuntimeError(f"{name}: requested waveform start changed")
            print(
                f"{name}: {len(waveform)} finite increasing IMR rows; spin error {spin_error:.3e}; NQC residual {max(residuals):.3e}; IMR phase residual {phase_residual:.3e}",
                flush=True,
            )
            if name == "default":
                rejection = work / "cli-rejection.err"
                run_checked(
                    [
                        str(project / "sebobv2"),
                        str(parameter_file),
                        "1",
                        "0.4",
                        "-0.3",
                        "0.01118",
                        "50",
                        "2.4627455127717882e-05",
                    ],
                    work,
                    work / "cli-rejection.out",
                    rejection,
                    expected_status=1,
                )
                expected_error = "Error: projected attachment requires vector spin inputs from a parameter file; scalar chi1 and chi2 command-line overrides are not supported.\n"
                if rejection.read_text(encoding="utf-8") != expected_error:
                    raise RuntimeError("Projected CLI rejection changed")
            if baseline is not None:
                fields["use_projected_attachment"] = False
                parameter_file.write_text(
                    "".join(f"{key} = {value}\n" for key, value in fields.items()),
                    encoding="utf-8",
                )
                scalar = work / f"{name}.scalar"
                trusted = work / f"{name}.trusted"
                run_checked(
                    [str(project / "sebobv2"), str(parameter_file)],
                    work,
                    scalar,
                    work / "scalar.err",
                )
                run_checked(
                    [str(baseline / "sebobv2"), str(parameter_file)],
                    work,
                    trusted,
                    work / "trusted.err",
                )
                read_waveform(scalar)
                if scalar.read_bytes() != trusted.read_bytes():
                    raise RuntimeError(
                        f"{name}: scalar output changed from the unchanged executable"
                    )
                print(f"{name}: scalar control byte-identical", flush=True)


if __name__ == "__main__":
    import argparse
    import doctest

    results = doctest.testmod()
    if results.failed:
        sys.exit(1)
    argument_parser = argparse.ArgumentParser(
        description="Check projected SEBOBv2 generated-executable numerics."
    )
    argument_parser.add_argument("--project", type=Path)
    argument_parser.add_argument("--baseline", type=Path)
    cli_args = argument_parser.parse_args()
    if cli_args.project is not None:
        main(cli_args.project, cli_args.baseline)
