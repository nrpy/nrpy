"""
Generate the C orchestrator for a single-particle debug run.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc


def integrate_single_particle_DEBUG() -> None:
    """
    Generate the C orchestrator for a single-particle debug run.

    This function traces the full path of a single particle, writing its state
    at each internal step of the integrator to a text file for validation.
    """
    name = "integrate_single_particle_DEBUG"
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "gsl/gsl_errno.h",
        "gsl/gsl_odeiv2.h",
    ]
    desc = "DEBUG integrator for a single massive particle. Writes the full trajectory to a text file."
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double start_y[8],
                  double final_y_state[8]"""
    body = r"""
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step * step = gsl_odeiv2_step_alloc(T, 8);
    gsl_odeiv2_control * control = gsl_odeiv2_control_yp_new(1e-14, 1e-14);
    gsl_odeiv2_evolve * evol = gsl_odeiv2_evolve_alloc(8);
    gsl_params gsl_parameters = {commondata, params, metric};
    gsl_odeiv2_system sys = {ode_gsl_wrapper_massive, NULL, 8, &gsl_parameters};
    double y_c[8];
    double t = 0.0, dt = 0.01;
    for (int j = 0; j < 8; j++) { y_c[j] = start_y[j]; }
    FILE *fp = fopen("massive_particle_path.txt", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open massive_particle_path.txt for writing.n");
        exit(1);
    }
    fprintf(fp, "# ProperTime_tau	CoordTime_t	x	y	z	u^t	u^x	u^y	u^z\n");
    printf("Starting debug trace for single massive particle...\n");
    printf("Step | Proper Time (τ) | Coord Time (t) |      x     |      y     |      z     |      u^t   \n");
    printf("-------------------------------------------------------------------------------------------\n");
    for (int i = 0; i < 2000000; i++) {
        int status = gsl_odeiv2_evolve_apply(evol, control, step, &sys, &t, 1e10, &dt, y_c);
        fprintf(fp, "%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	%.6e\n",
                t, y_c[0], y_c[1], y_c[2], y_c[3], y_c[4], y_c[5], y_c[6], y_c[7]);
        if (i % 500 == 0) {
            printf("%4d | %15.4e | %14.4f | %10.4f | %10.4f | %10.4f | %10.4f\n",
                   i, t, y_c[0], y_c[1], y_c[2], y_c[3], y_c[4]);
        }
        const double r_sq = y_c[1]*y_c[1] + y_c[2]*y_c[2] + y_c[3]*y_c[3];
        const double r_horizon = commondata->M_scale * (1.0 + sqrt(1.0 - commondata->a_spin*commondata->a_spin));
        if (status != GSL_SUCCESS) { printf("Termination: GSL ERROR (status = %d)\n", status); break; }
        if (r_sq < r_horizon*r_horizon) { printf("Termination: Fell below event horizon (r=%.2f)\n", sqrt(r_sq)); break; }
        if (r_sq > commondata->r_escape*commondata->r_escape) { printf("Termination: Escaped to r > %.1f\n", commondata->r_escape); break; }
        if (fabs(y_c[4]) > commondata->ut_max) { printf("Termination: Runaway u^t > %.1e\n", commondata->ut_max); break; }
        if (fabs(y_c[0]) > commondata->t_max_integration) { printf("Termination: Exceeded max integration time t > %.1f\n", commondata->t_max_integration); break; }
    }
    for(int j=0; j<8; j++) { final_y_state[j] = y_c[j]; }
    fclose(fp);
    gsl_odeiv2_evolve_free(evol);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);
    """

    cfc.register_CFunction(
        name=name,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
        include_CodeParameters_h=True,
    )