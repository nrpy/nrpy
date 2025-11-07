"""
Generates the main C integration loop for a single massive particle.

This high-performance "production" version uses the GSL driver to ensure
the state is returned at the exact requested snapshot times.

Author: Dalton J. Moone
"""

# Step 0: Import core nrpy modules
import nrpy.c_function as cfc

def integrate_single_particle() -> None:
    """
    Generate and register the C function for integrating a single particle.
    """
    # Step 1: Define C function metadata
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "gsl/gsl_errno.h", "gsl/gsl_odeiv2.h", "<math.h>"]
    desc = r"""@brief Integrates a single massive particle path between two times.
    
    This function uses the GSL driver, which internally uses an adaptive
    step-size algorithm (RKF45) to evolve the particle's state vector y_in_out
    from t_start to t_end, returning the state at the precise t_end.
    
    @param[in]      commondata  Pointer to commondata struct.
    @param[in]      params      Pointer to params_struct.
    @param[in]      metric      Pointer to metric_params struct.
    @param[in]      t_start     The starting coordinate time (t) for the integration.
    @param[in]      t_end       The ending coordinate time (t) for the integration.
    @param[in,out]  y_in_out    The 8-component state vector. Input is the state at t_start, output is the state at t_end.
    
    @return 0 on success, 1 on GSL failure or if particle is terminated.
    """
    cfunc_type = "int"
    name = "integrate_single_particle"
    params = """const commondata_struct *restrict commondata,
    const params_struct *restrict params,
    const metric_params *restrict metric,
    const double t_start, const double t_end,
    double y_in_out[8]"""

    # Step 2: Define C function body
    body = r"""
    // Define the GSL ODE system
    gsl_params gsl_parameters = {commondata, params, metric};
    gsl_odeiv2_system sys = {ode_gsl_wrapper_massive, NULL, 8, &gsl_parameters};
    
    // Set up the GSL driver with adaptive step-size control (RKF45)
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-11, 1e-11);
    
    // The GSL driver takes coordinate time 't' as its independent variable.
    // We use the first component of the state vector, y[0], for this.
    double t = t_start;
    
    // The driver will take internal steps to reach t_end precisely.
    int status = gsl_odeiv2_driver_apply(d, &t, t_end, y_in_out);

    // Free the driver memory *before* any early returns.
    gsl_odeiv2_driver_free(d);

    if (status != GSL_SUCCESS) {
        // Don't print an error here; the orchestrator will check the status.
        return 1; // Return failure code
    }

    // Robustness check after the step to terminate unphysical/escaped particles.
    const double r_sq = y_in_out[1]*y_in_out[1] + y_in_out[2]*y_in_out[2] + y_in_out[3]*y_in_out[3];
    const double r_horizon = commondata->M_scale * (1.0 + sqrt(1.0 - commondata->a_spin*commondata->a_spin));

    if (r_sq < r_horizon*r_horizon || r_sq > commondata->r_escape*commondata->r_escape || fabs(y_in_out[4]) > commondata->ut_max) {
        return 1; // Return failure/termination code
    }

    return 0; // Return success code
    """
    
    # Step 3: Register the C function
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body, 
        include_CodeParameters_h=True
    )