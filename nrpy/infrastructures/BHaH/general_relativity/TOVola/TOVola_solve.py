"""
Register TOVola code TOVola_solve.c.

TOVola creates Tolman-Oppenheimer-Volkoff spherically symmetric initial data,
 typically for single neutron stars.

Authors: David Boyer
         Zachariah B. Etienne
         zachetie **at** gmail **dot* com
"""

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.tov.TOV_equations as tov_eqs


def register_CFunction_TOVola_solve() -> None:
    """Register C function TP_solve(), main driver function for pseudospectral solve."""
    includes = ["BHaH_defines.h", "gsl/gsl_errno.h", "gsl/gsl_odeiv2.h"]
    prefunc = r"""
#define ODE_SOLVER_DIM 4
#define TOVOLA_PRESSURE 0
#define TOVOLA_NU 1
#define TOVOLA_MASS 2
#define TOVOLA_R_ISO 3

#define NEGATIVE_R_INTERP_BUFFER 11

/* Structure to hold raw TOV data */
typedef struct {
  // Current state variables
  REAL rho_baryon;
  REAL rho_energy;
  REAL r_lengthscale;

  REAL *restrict r_Schw_arr;
  REAL *restrict rho_energy_arr;
  REAL *restrict rho_baryon_arr;
  REAL *restrict P_arr;
  REAL *restrict M_arr;
  REAL *restrict nu_arr;
  REAL *restrict r_iso_arr;
  int numels_alloced_TOV_arr;

  int numpoints_actually_saved;

  const commondata_struct *restrict commondata;
} TOVola_data_struct;

/* Mock GRHayLib Functions for Simple Polytrope */
void ghl_hybrid_get_K_and_Gamma(const TOVola_data_struct *TOVdata, REAL rho_baryon, REAL *restrict K, REAL *restrict Gamma) {
  // For Simple Polytrope, K and Gamma are constants
  *K = TOVdata->commondata->poly_eos_K;         // Retrieved from TOVdata
  *Gamma = TOVdata->commondata->poly_eos_Gamma; // Retrieved from TOVdata
}

void ghl_hybrid_compute_P_cold_and_eps_cold(const TOVola_data_struct *TOVdata, REAL rho_baryon, REAL *restrict P, REAL *restrict eps) {
  // For Simple Polytrope: P = K * rho_baryon^Gamma
  REAL K = TOVdata->commondata->poly_eos_K;
  REAL Gamma = TOVdata->commondata->poly_eos_Gamma;
  *P = K * pow(rho_baryon, Gamma);
  // Assuming eps = P / ((Gamma - 1) * rho_baryon)
  *eps = *P / ((Gamma - 1.0) * rho_baryon);
}

/* Exception handler to prevent negative pressures */
void TOVola_exception_handler(REAL r, REAL y[]) {
  // Ensure pressure does not become negative due to numerical errors
  if (y[TOVOLA_PRESSURE] < 0) {
    y[TOVOLA_PRESSURE] = 0;
  }
}

/* Termination condition for the integration */
int TOVola_do_we_terminate(REAL r, REAL y[], TOVola_data_struct *TOVdata) {
  /* if (TOVdata->eos_type == TABULATED_EOS) { */
  /*     // Not implemented in this standalone version */
  /*     // Return 0 to continue integration */
  /*     return 0; */
  /* } */
  /* else { */
  if (y[TOVOLA_PRESSURE] <= 0.0) { // For Simple and Piecewise Polytrope
    return 1;
  }
  /* } */

  return 0; // Continue integration
}

/* Evaluate rho_baryon and rho_energy based on the EOS type */
void TOVola_evaluate_rho_and_eps(REAL r, const REAL y[], TOVola_data_struct *TOVdata) {
  // Simple Polytrope
  /* if (TOVdata->eos_type == SIMPLE_POLYTROPE) { */
  REAL aK, aGamma;
  REAL aRho_baryon = TOVdata->rho_baryon;
  REAL eps, aPress;

  // Retrieve K and Gamma from GRHayL
  ghl_hybrid_get_K_and_Gamma(TOVdata, aRho_baryon, &aK, &aGamma);
  TOVdata->rho_baryon = pow(y[TOVOLA_PRESSURE] / aK, 1.0 / aGamma);
  aRho_baryon = TOVdata->rho_baryon;
  ghl_hybrid_compute_P_cold_and_eps_cold(TOVdata, aRho_baryon, &aPress, &eps);
  TOVdata->rho_energy = TOVdata->rho_baryon * (1.0 + eps);
  /* } */
  /*
  // Piecewise Polytrope (Not implemented in this standalone version)
  else if (TOVdata->eos_type == PIECEWISE_POLYTROPE) {
  fprintf(stderr, "PIECEWISE_POLYTROPE EOS not implemented.\n");
  exit(EXIT_FAILURE);
  }
  */
  /*
  // Tabulated EOS (Not implemented in this standalone version)
  else if (TOVdata->eos_type == TABULATED_EOS) {
  // Not implemented
  fprintf(stderr, "TABULATED_EOS not implemented.\n");
  exit(EXIT_FAILURE);
  }
  */
}

/* The main ODE function for GSL */
int TOVola_ODE(REAL r_Schw, const REAL y[], REAL dydr_Schw[], void *params) {
  // Cast params to TOVdata_struct
  TOVola_data_struct *TOVdata = (TOVola_data_struct *)params;

  // Evaluate rho_baryon and rho_energy based on current state
  TOVola_evaluate_rho_and_eps(r_Schw, y, TOVdata);

  // Dereference the struct to use rho_energy
  REAL rho_energy = TOVdata->rho_energy;

  if (isnan(rho_energy)) {
    // Outside the star gives NaNs from the pow function, but we know they
    // should be zeros.
    rho_energy = 0.0;
  }

  // At the center of the star (r_Schw == 0), the TOV equations diverge, so we set reasonable values here.
  if (r_Schw == 0) {
    dydr_Schw[TOVOLA_PRESSURE] = 0.0; // dP/dr
    dydr_Schw[TOVOLA_NU] = 0.0;       // dnu/dr
    dydr_Schw[TOVOLA_MASS] = 0.0;     // dM/dr
    dydr_Schw[TOVOLA_R_ISO] = 1.0;    // dr_iso/dr
    return GSL_SUCCESS;
  }

  {
    // TOV Equations
"""
    tov = tov_eqs.TOV_Equations()
    prefunc += ccg.c_codegen(
        [tov.dP_dr, tov.dnu_dr, tov.dM_dr, tov.dr_iso_dr],
        [
            "dydr_Schw[TOVOLA_PRESSURE]",
            "dydr_Schw[TOVOLA_NU]",
            "dydr_Schw[TOVOLA_MASS]",
            "dydr_Schw[TOVOLA_R_ISO]",
        ],
    )
    prefunc += r"""
  }
  if (y[TOVOLA_R_ISO] > 0 && fabs(dydr_Schw[TOVOLA_R_ISO]) > 0) {
    TOVdata->r_lengthscale = fabs(y[TOVOLA_R_ISO] / dydr_Schw[TOVOLA_R_ISO]);
  }

  return GSL_SUCCESS;
}

/* Placeholder Jacobian function required by GSL */
int TOVola_jacobian_placeholder(REAL t, const REAL y[], REAL *restrict dfdy, REAL dfdt[], void *params) {
  // Jacobian is not necessary for the TOV solution, but GSL requires some
  // function Leave it empty as it does not affect the final results
  return GSL_SUCCESS;
}

/* Initialize the ODE variables */
void TOVola_get_initial_condition(REAL y[], TOVola_data_struct *TOVdata) {
  // Simple Polytrope
  /* if (TOVdata->eos_type == SIMPLE_POLYTROPE) { */
  REAL aK, aGamma;
  REAL rhoC_baryon = TOVdata->commondata->initial_central_density;

  // Retrieve K and Gamma from GRHayL
  ghl_hybrid_get_K_and_Gamma(TOVdata, rhoC_baryon, &aK, &aGamma);
  y[TOVOLA_PRESSURE] = aK * pow(rhoC_baryon, aGamma); // Pressure
  y[TOVOLA_NU] = 0.0;                                 // nu
  y[TOVOLA_MASS] = 0.0;                               // Mass
  y[TOVOLA_R_ISO] = 0.0;                              // r_iso

  // Assign initial conditions
  TOVdata->rho_baryon = rhoC_baryon;
  TOVdata->rho_energy = pow(y[TOVOLA_PRESSURE] / aK, 1.0 / aGamma) + y[TOVOLA_PRESSURE] / (aGamma - 1.0);
  // Pinitial is no longer needed as it's part of y[TOVOLA_PRESSURE]
  /* } */
  /*
  // Piecewise Polytrope (Not implemented in this standalone version)
  else if (TOVdata->eos_type == PIECEWISE_POLYTROPE) {
  // Implement similarly to Simple Polytrope if needed
  fprintf(stderr, "PIECEWISE_POLYTROPE EOS not implemented.\n");
  exit(EXIT_FAILURE);
  }
  */
  /*
  // Tabulated EOS (Not implemented in this standalone version)
  else if (TOVdata->eos_type == TABULATED_EOS) {
  // Not implemented
  fprintf(stderr, "TABULATED_EOS not implemented.\n");
  exit(EXIT_FAILURE);
  }
  */

  printf("Initial Conditions Set: P = %.6e, nu = %.6e, M = %.6e, r_iso = %.6e\n", y[TOVOLA_PRESSURE], y[TOVOLA_NU], y[TOVOLA_MASS], y[TOVOLA_R_ISO]);
}

/* Assign constants after each integration step */
void TOVola_assign_constants(REAL c[], TOVola_data_struct *TOVdata) {
  // Assign the densities
  c[0] = TOVdata->rho_energy; // Total energy density
  c[1] = TOVdata->rho_baryon; // Baryon density

  // Handle NaN cases
  if (isnan(TOVdata->rho_energy)) {
    c[0] = 0.0;
  }
}

/* Function to set up the GSL ODE system and driver */
static int setup_ode_system(const char *ode_method, gsl_odeiv2_system *system, gsl_odeiv2_driver **driver, TOVola_data_struct *TOVdata) {
  const commondata_struct *restrict commondata = TOVdata->commondata;
  system->function = TOVola_ODE;
  system->jacobian = TOVola_jacobian_placeholder;
  system->dimension = 4; // Hardcoded as per requirements
  system->params = TOVdata;

  if (strcmp(ode_method, "ARKF") == 0) {
    *driver = gsl_odeiv2_driver_alloc_y_new(system, gsl_odeiv2_step_rkf45, commondata->initial_ode_step_size, commondata->ode_error_limit,
                                            commondata->ode_error_limit);
  } else if (strcmp(ode_method, "ADP8") == 0) {
    *driver = gsl_odeiv2_driver_alloc_y_new(system, gsl_odeiv2_step_rk8pd, commondata->initial_ode_step_size, commondata->ode_error_limit,
                                            commondata->ode_error_limit);
  } else {
    fprintf(stderr, "Invalid ODE method. Use 'ARKF' or 'ADP8'.\n");
    return -1;
  }

  if (*driver == NULL) {
    fprintf(stderr, "Failed to allocate GSL ODE driver.\n");
    return -1;
  }

  /* Set minimum and maximum step sizes */
  gsl_odeiv2_driver_set_hmin(*driver, commondata->min_step_size);
  gsl_odeiv2_driver_set_hmax(*driver, commondata->max_step_size);

  return 0;
}

/* Initialize TOVola_data_struct structure with initial allocation */
static int initialize_tovola_data(TOVola_data_struct *TOVdata) {
  TOVdata->r_Schw_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->rho_energy_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->rho_baryon_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->P_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->M_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->nu_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numels_alloced_TOV_arr);
  TOVdata->r_iso_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numels_alloced_TOV_arr);

  if (!TOVdata->r_Schw_arr || !TOVdata->rho_energy_arr || !TOVdata->rho_baryon_arr || !TOVdata->P_arr || !TOVdata->M_arr || !TOVdata->nu_arr ||
      !TOVdata->r_iso_arr) {
    fprintf(stderr, "Memory allocation failed for TOVola_data_struct.\n");
    return -1;
  }
  return 0;
}

/* Free TOVola_data_struct structure */
static void free_tovola_data(TOVola_data_struct *TOVdata) {
  free(TOVdata->r_Schw_arr);
  free(TOVdata->rho_energy_arr);
  free(TOVdata->rho_baryon_arr);
  free(TOVdata->P_arr);
  free(TOVdata->M_arr);
  free(TOVdata->nu_arr);
  free(TOVdata->r_iso_arr);
  TOVdata->numels_alloced_TOV_arr = 0;
}

/* Normalize and set data */
void TOVola_Normalize_and_set_data_integrated(TOVola_data_struct *TOVdata, REAL *restrict r_Schw, REAL *restrict rho_energy,
                                              REAL *restrict rho_baryon, REAL *restrict P, REAL *restrict M, REAL *restrict expnu,
                                              REAL *restrict exp4phi, REAL *restrict r_iso) {
  printf("TOVola Normalizing raw TOV data...\n");

  /* Check if there are enough points to normalize */
  if (TOVdata->numpoints_actually_saved < 2) {
    fprintf(stderr, "Not enough data points to normalize.\n");
    exit(EXIT_FAILURE);
  }

  /* Copy raw data to normalized arrays */
  for (int i = 0; i < TOVdata->numpoints_actually_saved; i++) {
    r_Schw[i] = TOVdata->r_Schw_arr[i];
    rho_energy[i] = TOVdata->rho_energy_arr[i];
    rho_baryon[i] = TOVdata->rho_baryon_arr[i];
    P[i] = TOVdata->P_arr[i];
    M[i] = TOVdata->M_arr[i];
    expnu[i] = TOVdata->nu_arr[i];
    r_iso[i] = TOVdata->r_iso_arr[i];
  }

  /* Surface values for normalization */
  const REAL R_Schw_surface = r_Schw[TOVdata->numpoints_actually_saved - 1];
  const REAL M_surface = M[TOVdata->numpoints_actually_saved - 1];
  const REAL r_iso_surface = r_iso[TOVdata->numpoints_actually_saved - 1];
  const REAL nu_surface = expnu[TOVdata->numpoints_actually_saved - 1];

  const REAL normalize = 0.5 * (sqrt(R_Schw_surface * (R_Schw_surface - 2.0 * M_surface)) + R_Schw_surface - M_surface) / r_iso_surface;

  /* Normalize r_iso and calculate expnu and exp4phi */
  for (int i = 0; i < TOVdata->numpoints_actually_saved; i++) {
    r_iso[i] *= normalize;
    expnu[i] = exp(expnu[i] - nu_surface + log(1.0 - 2.0 * M_surface / R_Schw_surface));
    exp4phi[i] = (r_Schw[i] / r_iso[i]) * (r_Schw[i] / r_iso[i]);
  }
  printf("Normalization of raw data complete!\n");
}

/* Extend data to r<0, to ensure we can interpolate to r=0 */
void extend_to_negative_r(REAL *restrict arr, const REAL parity, REAL *restrict tmp, const TOVola_data_struct *restrict TOVdata) {
  for(int i=0;i<NEGATIVE_R_INTERP_BUFFER; i++) tmp[i] = parity * arr[NEGATIVE_R_INTERP_BUFFER - i - 1];
  for(int i=0;i<TOVdata->numpoints_actually_saved; i++) tmp[i+NEGATIVE_R_INTERP_BUFFER] = arr[i];
  memcpy(arr, tmp, sizeof(REAL) * (TOVdata->numpoints_actually_saved+NEGATIVE_R_INTERP_BUFFER));
}
"""
    desc = "Driver routine for TOV solve."
    name = "TOVola_solve"
    params = (
        "const commondata_struct *restrict commondata, ID_persist_struct *ID_persist"
    )
    body = r"""
  printf("Starting TOV Integration using GSL for TOVola...\n");

  REAL current_position = 0;
  const char *ode_method = "ADP8"; // Choose between "ARKF" and "ADP8"

  /* Set up ODE system and driver */
  TOVola_data_struct TOVdata_tmp; // allocates memory for the pointer below.
  TOVola_data_struct *restrict TOVdata = &TOVdata_tmp;
  gsl_odeiv2_system system;
  gsl_odeiv2_driver *driver;
  TOVdata->commondata = commondata;
  TOVdata->numpoints_actually_saved = 0;
  if (setup_ode_system(ode_method, &system, &driver, TOVdata) != 0) {
    fprintf(stderr, "Failed to set up ODE system.\n");
    exit(EXIT_FAILURE);
  }

  /* Initialize ODE variables */
  REAL y[ODE_SOLVER_DIM];
  REAL c[2];
  TOVola_get_initial_condition(y, TOVdata);
  TOVola_assign_constants(c, TOVdata);

  /* Initial memory allocation */
  TOVdata->numels_alloced_TOV_arr = 1024;
  if (initialize_tovola_data(TOVdata) != 0) {
    gsl_odeiv2_driver_free(driver);
    fprintf(stderr, "Failed to initialize TOVola_data_struct.\n");
    exit(EXIT_FAILURE);
  }

  /* Integration loop */
  TOVdata->r_lengthscale = TOVdata->commondata->initial_ode_step_size; // initialize dr to a crazy small value in double precision.
  for (int i = 0; i < TOVdata->commondata->ode_max_steps; i++) {
    REAL dr = 0.01 * TOVdata->r_lengthscale;
    if (TOVdata->rho_baryon < 0.05 * TOVdata->commondata->initial_central_density) {
      // To get a super-accurate mass, reduce the dr sampling near the surface of the star.
      dr = 1e-6 * TOVdata->r_lengthscale;
    }
    /* Exception handling */
    TOVola_exception_handler(current_position, y);

    /* Apply ODE step */
    int status = gsl_odeiv2_driver_apply(driver, &current_position, current_position + dr, y);
    if (status != GSL_SUCCESS) {
      fprintf(stderr, "GSL ODE solver failed with status %d.\n", status);
      gsl_odeiv2_driver_free(driver);
      exit(EXIT_FAILURE);
    }

    /* Post-step exception handling */
    TOVola_exception_handler(current_position, y);

    /* Evaluate densities */
    TOVola_evaluate_rho_and_eps(current_position, y, TOVdata);
    TOVola_assign_constants(c, TOVdata);
    /* Check if reallocation is needed */
    if (TOVdata->numpoints_actually_saved + NEGATIVE_R_INTERP_BUFFER + 1 >= TOVdata->numels_alloced_TOV_arr) {
      // Update arr_size instead of modifying the macro
      const int new_arr_size = 1.5 * TOVdata->numels_alloced_TOV_arr;
      TOVdata->numels_alloced_TOV_arr = new_arr_size;
      TOVdata->r_Schw_arr = realloc(TOVdata->r_Schw_arr, sizeof(REAL) * new_arr_size);
      TOVdata->rho_energy_arr = realloc(TOVdata->rho_energy_arr, sizeof(REAL) * new_arr_size);
      TOVdata->rho_baryon_arr = realloc(TOVdata->rho_baryon_arr, sizeof(REAL) * new_arr_size);
      TOVdata->P_arr = realloc(TOVdata->P_arr, sizeof(REAL) * new_arr_size);
      TOVdata->M_arr = realloc(TOVdata->M_arr, sizeof(REAL) * new_arr_size);
      TOVdata->nu_arr = realloc(TOVdata->nu_arr, sizeof(REAL) * new_arr_size);
      TOVdata->r_iso_arr = realloc(TOVdata->r_iso_arr, sizeof(REAL) * new_arr_size);

      if (!TOVdata->r_Schw_arr || !TOVdata->rho_energy_arr || !TOVdata->rho_baryon_arr || !TOVdata->P_arr || !TOVdata->M_arr || !TOVdata->nu_arr ||
          !TOVdata->r_iso_arr) {
        fprintf(stderr, "Memory reallocation failed during integration.\n");
        gsl_odeiv2_driver_free(driver);
        exit(EXIT_FAILURE);
      }
    }

    /* Store data */
    TOVdata->r_Schw_arr[TOVdata->numpoints_actually_saved] = current_position;
    TOVdata->rho_energy_arr[TOVdata->numpoints_actually_saved] = c[0];
    TOVdata->rho_baryon_arr[TOVdata->numpoints_actually_saved] = c[1];
    TOVdata->P_arr[TOVdata->numpoints_actually_saved] = y[TOVOLA_PRESSURE];
    TOVdata->M_arr[TOVdata->numpoints_actually_saved] = y[TOVOLA_MASS];
    TOVdata->nu_arr[TOVdata->numpoints_actually_saved] = y[TOVOLA_NU];
    TOVdata->r_iso_arr[TOVdata->numpoints_actually_saved] = y[TOVOLA_R_ISO];
    TOVdata->numpoints_actually_saved++;

    // r_SchwArr_np,rhoArr_np,rho_baryonArr_np,PArr_np,mArr_np,exp2phiArr_np,confFactor_exp4phi_np,r_isoArr_np),
    // printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e soln\n", current_position, dr, c[0], c[1], y[TOVOLA_PRESSURE], y[TOVOLA_MASS], y[TOVOLA_NU]);

    /* Termination condition */
    if (TOVola_do_we_terminate(current_position, y, TOVdata)) {
      printf("Finished Integration at position %.6e with Mass %.14e\n", current_position, y[TOVOLA_MASS]);
      break;
    }
  }

  /* Cleanup */
  gsl_odeiv2_driver_free(driver);
  printf("ODE Solver using GSL for TOVola Shutting Down...\n");

  {
    // Data in TOVdata->*_arr are stored at r=TOVdata->commondata->initial_ode_step_size > 0 up to the stellar surface.
    // However, we may need data at r=0, which would require extrapolation.
    // To prevent that, we copy INTERP_BUFFER data points from r>0 to r<0 so that we can always interpolate.
    REAL *restrict tmp = malloc(sizeof(REAL) * (TOVdata->numpoints_actually_saved + NEGATIVE_R_INTERP_BUFFER));

    //printf("Rbefor = %.15e\n", TOVdata->r_Schw_arr[TOVdata->numpoints_actually_saved-1]);

    extend_to_negative_r(TOVdata->r_Schw_arr, -1.0, tmp, TOVdata);
    extend_to_negative_r(TOVdata->rho_energy_arr, +1.0, tmp, TOVdata);
    extend_to_negative_r(TOVdata->rho_baryon_arr, +1.0, tmp, TOVdata);
    extend_to_negative_r(TOVdata->P_arr, +1.0, tmp, TOVdata);
    extend_to_negative_r(TOVdata->M_arr, +1.0, tmp, TOVdata);
    extend_to_negative_r(TOVdata->nu_arr, +1.0, tmp, TOVdata);
    extend_to_negative_r(TOVdata->r_iso_arr, -1.0, tmp, TOVdata);

    free(tmp);
    TOVdata->numpoints_actually_saved += NEGATIVE_R_INTERP_BUFFER;

    //printf("Rafter = %.15e\n", TOVdata->r_Schw_arr[TOVdata->numpoints_actually_saved-1]);
    //for(int i=0;i<TOVdata->numpoints_actually_saved; i++) printf("%e %e %e iiii\n", TOVdata->r_Schw_arr[i], TOVdata->r_iso_arr[i], TOVdata->rho_energy_arr[i]);;
  }

  /* Allocate and populate ID_persist_struct arrays */
  ID_persist->r_Schw_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numpoints_actually_saved);
  ID_persist->rho_energy_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numpoints_actually_saved);
  ID_persist->rho_baryon_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numpoints_actually_saved);
  ID_persist->P_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numpoints_actually_saved);
  ID_persist->M_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numpoints_actually_saved);
  ID_persist->expnu_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numpoints_actually_saved);
  ID_persist->exp4phi_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numpoints_actually_saved);
  ID_persist->r_iso_arr = (REAL *restrict)malloc(sizeof(REAL) * TOVdata->numpoints_actually_saved);

  if (!ID_persist->r_Schw_arr || !ID_persist->rho_energy_arr || !ID_persist->rho_baryon_arr || !ID_persist->P_arr || !ID_persist->M_arr ||
      !ID_persist->expnu_arr || !ID_persist->exp4phi_arr || !ID_persist->r_iso_arr) {
    free_tovola_data(TOVdata);
    fprintf(stderr, "Memory allocation failed for ID_persist_struct arrays.\n");
    exit(EXIT_FAILURE);
  }

  ID_persist->numpoints_arr = TOVdata->numpoints_actually_saved;

  /* Normalize and set data */
  TOVola_Normalize_and_set_data_integrated(TOVdata, ID_persist->r_Schw_arr, ID_persist->rho_energy_arr, ID_persist->rho_baryon_arr, ID_persist->P_arr,
                                           ID_persist->M_arr, ID_persist->expnu_arr, ID_persist->exp4phi_arr, ID_persist->r_iso_arr);

  /* Free raw data as it's no longer needed */
  free_tovola_data(TOVdata);
"""
    cfc.register_CFunction(
        subdirectory="TOVola",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )
