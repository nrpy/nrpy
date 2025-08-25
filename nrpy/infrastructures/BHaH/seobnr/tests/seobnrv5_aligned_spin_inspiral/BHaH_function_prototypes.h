void BOB_aligned_spin_NQC_rhs(commondata_struct *restrict commondata, REAL *restrict amps, REAL *restrict omegas);
void BOB_aligned_spin_waveform(const REAL t, commondata_struct *restrict commondata, REAL *restrict waveform);
void BOB_aligned_spin_waveform_from_times(REAL *restrict times, REAL *restrict amps, REAL *restrict phases, const size_t nsteps_BOB,
                                          commondata_struct *restrict commondata);
void cmdline_input_and_parfile_parser(commondata_struct *restrict commondata, int argc, const char *argv[]);
int commondata_io(commondata_struct *restrict data, const char *restrict filename);
void commondata_struct_set_to_default(commondata_struct *restrict commondata);
double eval_abs_deriv(double t, void *params);
size_t find_local_minimum_index(REAL *restrict arr, size_t size, int order);
void handle_gsl_return_status(int status, int status_desired[], int num_desired, const char *restrict function_name);
int main(int argc, const char *argv[]);
void params_struct_set_to_default(commondata_struct *restrict commondata, griddata_struct *restrict griddata);
void SEOBNRv5_aligned_multidimensional_root_wrapper(gsl_multiroot_function_fdf f, const REAL *restrict x_guess, const size_t n,
                                                    REAL *restrict x_result);
size_t SEOBNRv5_aligned_spin_argrelmin(REAL *restrict arr, size_t nsteps_arr);
void SEOBNRv5_aligned_spin_augments(commondata_struct *restrict commondata);
void SEOBNRv5_aligned_spin_coefficients(commondata_struct *restrict commondata);
int SEOBNRv5_aligned_spin_flux(const REAL *restrict y, const REAL Hreal, const REAL Omega, const REAL Omega_circ, REAL *restrict f,
                               void *restrict params);
double complex SEOBNRv5_aligned_spin_gamma_wrapper(const REAL z_real, const REAL z_imag);
int SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit(const gsl_vector *restrict x, void *restrict params, gsl_vector *restrict f);
int SEOBNRv5_aligned_spin_IMR_waveform(commondata_struct *restrict commondata);
int SEOBNRv5_aligned_spin_initial_conditions_conservative(commondata_struct *restrict commondata);
int SEOBNRv5_aligned_spin_initial_conditions_dissipative(commondata_struct *restrict commondata);
int SEOBNRv5_aligned_spin_interpolate_dynamics(commondata_struct *restrict commondata, REAL *restrict dynamics_fine_prelim,
                                               const size_t nsteps_fine_prelim, const REAL t_peak, const int stop);
void SEOBNRv5_aligned_spin_interpolate_modes(commondata_struct *restrict commondata, const REAL dT);
REAL SEOBNRv5_aligned_spin_iterative_refinement(spline_data *sdata, double initial_left, double initial_right, int levels, double dt_initial,
                                                bool pr);
int SEOBNRv5_aligned_spin_NQC_corrections(commondata_struct *restrict commondata);
int SEOBNRv5_aligned_spin_ode_integration(commondata_struct *restrict commondata);
REAL SEOBNRv5_aligned_spin_radial_momentum_conditions(REAL x, void *restrict params);
int SEOBNRv5_aligned_spin_right_hand_sides(REAL t, const REAL *restrict y, REAL *restrict f, void *restrict params);
void SEOBNRv5_aligned_spin_unwrap(REAL *restrict angles_in, REAL *restrict angles_out, size_t nsteps_arr);
double complex SEOBNRv5_aligned_spin_waveform(REAL *restrict dynamics, commondata_struct *restrict commondata);
int SEOBNRv5_aligned_spin_waveform_from_dynamics(commondata_struct *restrict commondata);