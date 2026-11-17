"""
Generate C code for the 10M reference used in projected-spin attachment.

The calculation uses the stored inspiral trajectory to find the 10M state.
For waveforms that start inside 10M, it recovers the earlier reference by
continuing the initialized equations backward. The supplied spins remain
defined at the requested initial orbital frequency, and the requested
waveform start is preserved.

Author: Suchindram Dasgupta
        sd00113 **at** mix **dot** wvu **dot** edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


def register_CFunction_SEOBNRv5_projected_attachment_reference() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register the numerical 10M reference calculation.

    The raw EOB trajectory supplies an inward crossing without assuming that
    merger-side radius samples are monotone. For an orbit starting inside 10M,
    continue its actual initialized state backward with the same EOB equations.
    PN spins below their original spline interval are evolved backward from the
    supplied initial frequency, in the original frame. EOB and PN time are not
    identified: the existing model evaluates PN spins at the EOB frequency.

    This function owns reference-state parameters. Register the orbital and PN
    right-hand sides, spin dynamics, and PA integration before generating the
    complete executable. Call it after PA/ODE evolution and before projected
    special-amplitude coefficients. Solver failures terminate the calculation.

    :return: None during discovery, otherwise the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "t_r10M",
            "omega_r10M",
            "r_r10M",
            "phi_r10M",
            "prstar_r10M",
            "pphi_r10M",
            "chi1_r10_x",
            "chi1_r10_y",
            "chi1_r10_z",
            "chi2_r10_x",
            "chi2_r10_y",
            "chi2_r10_z",
            "lnhat_r10_x",
            "lnhat_r10_y",
            "lnhat_r10_z",
            "chi1_projected_r10",
            "chi2_projected_r10",
            "min_reference_omega_dot",
        ],
        commondata=True,
        add_to_parfile=False,
        add_to_set_CodeParameters_h=False,
    )
    par.register_CodeParameters(
        "bool",
        __name__,
        ["projected_reference_ready"],
        [False],
        commondata=True,
        add_to_parfile=False,
        add_to_set_CodeParameters_h=False,
    )

    prefunc = r"""
#include <float.h>

typedef struct {
  gsl_odeiv2_driver *driver;
  REAL initial[4];
} r10_orbit_context; // END STRUCT: orbital continuation data

/**
 * Continue a saved EOB state by a signed time interval.
 * @param[in,out] data Orbital driver and saved state.
 * @param duration Signed interval in total-mass units.
 * @param[out] state Continued EOB state.
 */
static void r10_orbit_state(r10_orbit_context *data, REAL duration, REAL state[4]) {
  memcpy(state, data->initial, 4 * sizeof(REAL));
  if (duration == 0.0)
    return;
  // reset_hstart resets both evolve and step state before setting the step size.
  int status = gsl_odeiv2_driver_reset_hstart(data->driver, copysign(0.01, duration));
  REAL time = 0.0;
  if (status == GSL_SUCCESS)
    status = gsl_odeiv2_driver_apply(data->driver, &time, duration, state);
  if (status != GSL_SUCCESS) {
    fprintf(stderr, "Error: 10M EOB continuation failed at duration %.17g: %s\n", duration, gsl_strerror(status));
    exit(EXIT_FAILURE);
  } // END IF: EOB continuation failed
} // END FUNCTION: continue saved EOB state

/**
 * Evaluate the local 10M crossing equation for GSL.
 * @param duration Signed time interval from the saved state.
 * @param[in,out] parameters Orbital continuation data.
 * @return Continued radius minus 10M.
 */
static double r10_radius_event(double duration, void *parameters) {
  REAL state[4];
  r10_orbit_state((r10_orbit_context *)parameters, duration, state);
  return state[0] - 10.0;
} // END FUNCTION: evaluate 10M crossing equation

typedef struct {
  commondata_struct *parameters;
  REAL min_omega_dot;
} r10_spin_context; // END STRUCT: PN reference evolution data

/**
 * Evaluate existing PN spin equations with orbital frequency as independent variable.
 * @param omega PN orbital frequency.
 * @param[in] state Spin and orbital-plane variables in the original frame.
 * @param[out] derivatives Derivatives with respect to frequency.
 * @param[in,out] parameters Mass parameters and recorded minimum frequency derivative.
 * @return GSL_SUCCESS or the PN evolution failure.
 */
static int r10_spin_frequency(double omega, const double state[], double derivatives[], void *parameters) {
  r10_spin_context *data = (r10_spin_context *)parameters;
  REAL full[NUMVARS_SPIN], rhs[NUMVARS_SPIN];
  memcpy(full, state, 9 * sizeof(REAL));
  full[OMEGA_PN] = omega;
  const int status = SEOBNRv5_quasi_precessing_spin_equations(0.0, full, rhs, data->parameters);
  if (status != GSL_SUCCESS)
    return status;
  data->min_omega_dot = fmin(data->min_omega_dot, rhs[OMEGA_PN]);
  if (!isfinite(rhs[OMEGA_PN]) || rhs[OMEGA_PN] <= 0.0)
    return GSL_EBADFUNC;
  for (size_t i = 0; i < 9; i++)
    derivatives[i] = rhs[i] / rhs[OMEGA_PN];
  return GSL_SUCCESS;
} // END FUNCTION: evaluate frequency-parametrized PN evolution
"""
    body = r"""
commondata->projected_reference_ready = false;
if (commondata->nsteps_raw == 0 || commondata->dynamics_raw == NULL) {
  fprintf(stderr, "Error: 10M reference requires the initialized raw EOB trajectory\n");
  exit(EXIT_FAILURE);
} // END IF: raw EOB trajectory absent

// Step 1: Select the first inward crossing, or the initialized short-start state.
// PA ends at max(10M, ...), so a preceding PA segment cannot cross below 10M.
size_t left = 0;
const REAL first_radius = commondata->dynamics_raw[IDX(0, R)];
if (first_radius > 10.0) {
  while (left + 1 < commondata->nsteps_raw && commondata->dynamics_raw[IDX(left + 1, R)] > 10.0)
    left++;
  if (left + 1 == commondata->nsteps_raw) {
    fprintf(stderr, "Error: raw EOB trajectory has no inward 10M crossing\n");
    exit(EXIT_FAILURE);
  } // END IF: inward 10M crossing absent
  if (commondata->dynamics_raw[IDX(left + 1, R)] == 10.0)
    left++;
} // END IF: reference follows EOB initialization

REAL state[4];
for (size_t i = 0; i < 4; i++)
  state[i] = commondata->dynamics_raw[IDX(left, i + 1)];
REAL base_time = commondata->dynamics_raw[IDX(left, TIME)];
REAL duration = 0.0;
commondata_struct orbital_parameters = *commondata;

// Step 2: Refine the physical crossing with local EOB evolution. An exact node
// is retained exactly, including the PA endpoint at 10M.
if (state[0] != 10.0) {
  gsl_odeiv2_system system = {SEOBNRv5_aligned_spin_right_hand_sides, NULL, 4, &orbital_parameters};
  gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rk8pd, 0.01, 1e-13, 1e-12);
  if (driver == NULL) {
    fprintf(stderr, "Error: 10M orbital driver allocation failed\n");
    exit(EXIT_FAILURE);
  } // END IF: orbital driver allocation failed
  gsl_odeiv2_driver_set_nmax(driver, 100000);
  r10_orbit_context data;
  data.driver = driver;
  memcpy(data.initial, state, sizeof(state));
  REAL lower = 0.0;
  REAL upper = 0.0;
  if (first_radius < 10.0) {
    gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, 4);
    gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(1e-13, 1e-12);
    gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(4);
    if (step == NULL || control == NULL || evolve == NULL) {
      fprintf(stderr, "Error: earlier EOB continuation allocation failed\n");
      exit(EXIT_FAILURE);
    } // END IF: earlier EOB allocation failed
    REAL time = 0.0, previous_time = 0.0, h = -0.01;
    size_t steps = 0;
    while (state[0] < 10.0 && steps < 100000) {
      memcpy(data.initial, state, sizeof(state));
      previous_time = time;
      const int status = gsl_odeiv2_evolve_apply(evolve, control, step, &system, &time, -2e9, &h, state);
      if (status != GSL_SUCCESS) {
        fprintf(stderr, "Error: earlier EOB continuation failed: %s\n", gsl_strerror(status));
        exit(EXIT_FAILURE);
      } // END IF: earlier EOB continuation failed
      steps++;
    } // END WHILE: approach first earlier 10M crossing
    gsl_odeiv2_evolve_free(evolve);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);
    if (state[0] < 10.0) {
      fprintf(stderr, "Error: backward EOB continuation did not reach 10M\n");
      exit(EXIT_FAILURE);
    } // END IF: earlier reference not constructed
    if (state[0] == 10.0) {
      memcpy(data.initial, state, sizeof(state));
      base_time += time;
    } // END IF: exact earlier 10M sample
    else {
      base_time += previous_time;
      lower = time - previous_time;
    } // END ELSE: bracket earlier 10M event
  } // END IF: requested orbit starts inside 10M
  else {
    upper = commondata->dynamics_raw[IDX(left + 1, TIME)] - base_time;
    size_t extensions = 0;
    // Re-integration can differ from a saved endpoint by its integration error.
    // Extend the evolution interval, never the radius or an interpolated state.
    while (r10_radius_event(upper, &data) > 0.0 && extensions < 32) {
      upper *= 2.0;
      extensions++;
    } // END WHILE: refine forward crossing interval
    if (upper <= 0.0 || r10_radius_event(upper, &data) > 0.0) {
      fprintf(stderr, "Error: local EOB evolution did not bracket 10M\n");
      exit(EXIT_FAILURE);
    } // END IF: local crossing not bracketed
  } // END ELSE: reference lies within inspiral
  gsl_function event = {r10_radius_event, &data};
  if (lower != upper)
    duration = root_finding_1d(lower, upper, &event);
  r10_orbit_state(&data, duration, state);
  gsl_odeiv2_driver_free(driver);
} // END IF: 10M event requires refinement

REAL rates[4];
const int orbital_status = SEOBNRv5_aligned_spin_right_hand_sides(0.0, state, rates, &orbital_parameters);
if (orbital_status != GSL_SUCCESS || !isfinite(rates[0]) || !isfinite(rates[1]) || rates[1] <= 0.0 || rates[0] >= 0.0) {
  fprintf(stderr, "Error: 10M reference is not a regular inward EOB state\n");
  exit(EXIT_FAILURE);
} // END IF: reference evolution is not inward
commondata->t_r10M = commondata->t_dynamics_raw_origin + base_time + duration;
commondata->r_r10M = state[0];
commondata->phi_r10M = state[1];
commondata->prstar_r10M = state[2];
commondata->pphi_r10M = state[3];
commondata->omega_r10M = rates[1];

// Step 3: Preserve the existing forward spin splines where they cover the
// reference. Only the missing earlier point is continued from the supplied
// spin state. The requested forward evolution is never restarted or replaced.
REAL spin[9];
if (commondata->omega_r10M >= commondata->omega_spin_min && commondata->omega_r10M <= commondata->omega_spin_max) {
  spin[0] = gsl_spline_eval(commondata->lnhat_x.spline, commondata->omega_r10M, commondata->lnhat_x.acc);
  spin[1] = gsl_spline_eval(commondata->lnhat_y.spline, commondata->omega_r10M, commondata->lnhat_y.acc);
  spin[2] = gsl_spline_eval(commondata->lnhat_z.spline, commondata->omega_r10M, commondata->lnhat_z.acc);
  spin[3] = gsl_spline_eval(commondata->chi1_x_spline.spline, commondata->omega_r10M, commondata->chi1_x_spline.acc);
  spin[4] = gsl_spline_eval(commondata->chi1_y_spline.spline, commondata->omega_r10M, commondata->chi1_y_spline.acc);
  spin[5] = gsl_spline_eval(commondata->chi1_z_spline.spline, commondata->omega_r10M, commondata->chi1_z_spline.acc);
  spin[6] = gsl_spline_eval(commondata->chi2_x_spline.spline, commondata->omega_r10M, commondata->chi2_x_spline.acc);
  spin[7] = gsl_spline_eval(commondata->chi2_y_spline.spline, commondata->omega_r10M, commondata->chi2_y_spline.acc);
  spin[8] = gsl_spline_eval(commondata->chi2_z_spline.spline, commondata->omega_r10M, commondata->chi2_z_spline.acc);
  commondata->chi1_projected_r10 = gsl_spline_eval(commondata->chi1_lnhat.spline, commondata->omega_r10M, commondata->chi1_lnhat.acc);
  commondata->chi2_projected_r10 = gsl_spline_eval(commondata->chi2_lnhat.spline, commondata->omega_r10M, commondata->chi2_lnhat.acc);
  commondata->min_reference_omega_dot = 0.0;
} // END IF: forward spin splines cover reference
else if (commondata->omega_r10M < commondata->omega_spin_min) {
  const REAL initial[9] = {0.0, 0.0, 1.0, commondata->chi1_x, commondata->chi1_y, commondata->chi1_z,
                         commondata->chi2_x, commondata->chi2_y, commondata->chi2_z};
  memcpy(spin, initial, sizeof(spin));
  r10_spin_context data = {commondata, DBL_MAX};
  gsl_odeiv2_system system = {r10_spin_frequency, NULL, 9, &data};
  gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rk8pd, -1e-6, 1e-13, 1e-12);
  if (driver == NULL) {
    fprintf(stderr, "Error: 10M spin driver allocation failed\n");
    exit(EXIT_FAILURE);
  } // END IF: spin driver allocation failed
  gsl_odeiv2_driver_set_nmax(driver, 100000);
  REAL omega = commondata->initial_omega;
  const int status = gsl_odeiv2_driver_apply(driver, &omega, commondata->omega_r10M, spin);
  gsl_odeiv2_driver_free(driver);
  if (status != GSL_SUCCESS) {
    fprintf(stderr, "Error: earlier PN spin reference failed on [%.17g, %.17g], minimum omega_dot %.17g: %s\n",
            commondata->omega_r10M, commondata->initial_omega, data.min_omega_dot, gsl_strerror(status));
    exit(EXIT_FAILURE);
  } // END IF: earlier PN reference failed
  commondata->min_reference_omega_dot = data.min_omega_dot;
  const REAL ln_inverse = 1.0 / sqrt(spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2]);
  for (size_t i = 0; i < 3; i++)
    spin[i] *= ln_inverse;
  commondata->chi1_projected_r10 = spin[0]*spin[3] + spin[1]*spin[4] + spin[2]*spin[5];
  commondata->chi2_projected_r10 = spin[0]*spin[6] + spin[1]*spin[7] + spin[2]*spin[8];
} // END ELSE IF: reference precedes supplied spin frequency
else {
  fprintf(stderr, "Error: 10M reference frequency %.17g exceeds the evolved PN spin interval [%.17g, %.17g]\n",
          commondata->omega_r10M, commondata->omega_spin_min, commondata->omega_spin_max);
  exit(EXIT_FAILURE);
} // END ELSE: reference exceeds forward spin evolution
commondata->lnhat_r10_x = spin[0];
commondata->lnhat_r10_y = spin[1];
commondata->lnhat_r10_z = spin[2];
commondata->chi1_r10_x = spin[3];
commondata->chi1_r10_y = spin[4];
commondata->chi1_r10_z = spin[5];
commondata->chi2_r10_x = spin[6];
commondata->chi2_r10_y = spin[7];
commondata->chi2_r10_z = spin[8];
commondata->projected_reference_ready = true;
"""
    desc = """
Prepare the 10M EOB and PN spin state for projected attachment.

The reference preserves the initialized orbital trajectory and spins at their
supplied frequency, including a consistent earlier continuation for short starts.
Times are measured from the requested waveform start in total-mass units.

@param[in,out] commondata Initialized orbital/spin dynamics and prepared reference.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_projected_attachment_reference"
    params = "commondata_struct *restrict commondata"
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
