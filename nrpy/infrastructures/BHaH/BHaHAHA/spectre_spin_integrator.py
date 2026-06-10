"""
Register a C function to perform numerical integrations over the apparent horizon surface for spin and mass diagnostics.
Needs: Integrands built in equations/general_relativity/bhahaha/SpECTRESpinEstimate.py

Author: Ralston Graves
"""

from typing import Union

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg

# Import the SpECTRESpinEstimate factory from the other module
from nrpy.equations.general_relativity.bhahaha.SpECTRESpinEstimate import (
    SpECTRESpinEstimate,
)
from nrpy.helpers.generic import clang_format


def register_CFunction_diagnostics_spectre_spin(
    CoordSystem: str = "Spherical",
    enable_rfm_precompute: bool = False,
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a C function that computes the SpECTRE-style spin vector.

    :param CoordSystem: The coordinate system to use, defaults to "Spherical".
    :param enable_rfm_precompute: Whether to enable RFM precompute.
    :param enable_fd_functions: Whether to enable finite-difference functions.
    :return: An NRPyEnv_type object if registration is successful, otherwise None.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(
            f"{__name__}.{register_CFunction_diagnostics_spectre_spin.__name__}",
            locals(),
        )
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""
    Compute the SpECTRE-style spin vector diagnostic and store it in the diagnostics struct.
    """
    cfunc_type = "int"
    cfunc_name = "diagnostics_spectre_spin"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    # Step 1: Get an instance of the symbolic calculator.
    spin_calc = SpECTRESpinEstimate[
        CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    ]

    # Step 2: Register the memory-backed gridfunctions needed by the diagnostic.
    _ = gri.register_gridfunctions_for_single_rank2(
        "SE_qDD",
        symmetry="sym01",
        dimension=2,
        group="AUXEVOL",
        gf_array_name="auxevol_gfs",
    )
    _ = gri.register_gridfunctions_for_single_rank1(
        "SE_XD",
        dimension=2,
        group="AUXEVOL",
        gf_array_name="auxevol_gfs",
    )
    _ = gri.register_gridfunctions_for_single_rank1(
        "zU",
        dimension=3,
        group="AUXEVOL",
        gf_array_name="auxevol_gfs",
    )

    gf_assignments = spin_calc.get_gridfunction_assignments(include_flux_density=False)
    gf_names = [str(sym) for sym in gf_assignments.keys()]

    lhss_precompute = [
        gri.BHaHGridFunction.access_gf(gf_name, gf_array_name="auxevol_gfs")
        for gf_name in gf_names
    ]
    rhss_precompute = list(gf_assignments.values())

    # Step 3: Generate C code for intermediate surface fields.
    precompute_c_code = ccg.c_codegen(
        rhss_precompute,
        lhss_precompute,
        enable_fd_codegen=True,
        enable_fd_functions=enable_fd_functions,
    )

    # Step 4: Retrieve the dictionary of all per-point integrands.
    integrands_dict = spin_calc.get_public_integrands()

    # Step 5: Extract the symbolic expressions we need to integrate.
    # According to the SpECTRESpinEstimate documentation, we need to integrate:
    # - 1 (for the Area A)
    # - x^i (for the centroid XU)
    # - R (for R0)
    # - x^i * R (for XRU)
    # - Omega (for O0)
    # - x^i * Omega (for XOU)
    # - z_alpha * Omega (for ZOU)
    # - |Omega| (for Oabs, used in the near-zero policy)

    # Note: The 'integrand' is the quantity 'f' in ∮ f dA.
    # The differential area element is dA = sqrt(q) * weights * dθ * dφ.
    # We will pass sqrt(q) to c_codegen as 'area_density' and handle the
    # weights and coordinate steps in the C code loop.

    area_density = integrands_dict["area_density"]

    # List of all symbolic quantities to be evaluated inside the loop
    integrand_c_vars = [
        "const REAL area_density",
        "const REAL A_integrand",
        "const REAL XU0_integrand",
        "const REAL XU1_integrand",
        "const REAL XU2_integrand",
        "const REAL R0_integrand",
        "const REAL XRU0_integrand",
        "const REAL XRU1_integrand",
        "const REAL XRU2_integrand",
        "const REAL O0_integrand",
        "const REAL XOU0_integrand",
        "const REAL XOU1_integrand",
        "const REAL XOU2_integrand",
        "const REAL ZOU0_integrand",
        "const REAL ZOU1_integrand",
        "const REAL ZOU2_integrand",
        "const REAL Oabs_integrand",
    ]

    sympy_expressions = [
        area_density,
        integrands_dict["area_integrand"],  # This is just 1.
        *integrands_dict["measurement_frame_xU"],
        integrands_dict["ricci_scalar"],
        *integrands_dict["x_times_R_integrand"],
        integrands_dict["spin_function"],
        *integrands_dict["xOmega_momentU"],
        *integrands_dict["zOmegaU"],
        integrands_dict["abs_omega_integrand"],
    ]

    # Step 6: Construct the body of the C function.
    body = r"""
    const int grid=0;
    const params_struct *restrict params = &griddata[grid].params;
    REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    const REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs; // for hh and its time-derivs
    REAL *restrict xx[3];
    for(int ww=0;ww<3;ww++) xx[ww] = griddata[grid].xx[ww];
#include "set_CodeParameters.h"

    // This diagnostic assumes that all required inputs, including the results
    // of any elliptic solves (e.g., for the spin potential 'z' or test modes 'z_alpha'),
    // are already populated in the auxevol gridfunctions.
    
    // Initialize all RunSums accumulators to zero.
    REAL A_sum = 0.0;
    REAL XU_sum[3] = {0.0, 0.0, 0.0};
    REAL R0_sum = 0.0;
    REAL XRU_sum[3] = {0.0, 0.0, 0.0};
    REAL O0_sum = 0.0;
    REAL XOU_sum[3] = {0.0, 0.0, 0.0};
    REAL ZOU_sum[3] = {0.0, 0.0, 0.0};
    REAL Oabs_sum = 0.0;
    
    // Set integration weights (e.g., for Simpson's rule).
    // This external function provides the 1D weights array.
    const REAL *restrict weights;
    int weight_stencil_size;
    bah_diagnostics_integration_weights(Nxx1, Nxx2, &weights, &weight_stencil_size);

#pragma omp parallel
{
#pragma omp for
    for (int i2 = 0; i2 < Nxx2 + 2 * NGHOSTS; i2++) {
        for (int i1 = 0; i1 < Nxx1 + 2 * NGHOSTS; i1++) {
            const REAL xx1 = xx[1][i1];
            for (int i0 = NGHOSTS; i0 < NGHOSTS + 1; i0++) {
"""
    body += precompute_c_code
    body += r"""
            } // END LOOP: for i0 over radial horizon-grid points
        } // END LOOP: for i1 over theta horizon-grid points
    } // END LOOP: for i2 over phi horizon-grid points
    // Private accumulators for each thread
    REAL A_sum_private = 0.0;
    REAL XU_sum_private[3] = {0.0, 0.0, 0.0};
    REAL R0_sum_private = 0.0;
    REAL XRU_sum_private[3] = {0.0, 0.0, 0.0};
    REAL O0_sum_private = 0.0;
    REAL XOU_sum_private[3] = {0.0, 0.0, 0.0};
    REAL ZOU_sum_private[3] = {0.0, 0.0, 0.0};
    REAL Oabs_sum_private = 0.0;

#pragma omp for
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
        const REAL weight2 = weights[(i2 - NGHOSTS) % weight_stencil_size];
        const REAL xx2 = xx[2][i2];
        for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
            const REAL weight1 = weights[(i1 - NGHOSTS) % weight_stencil_size];
            const REAL xx1 = xx[1][i1];
            // The horizon is a 2D surface, so we loop over a single radial index i0.
            // All quantities are evaluated on this surface.
            for (int i0 = NGHOSTS; i0 < NGHOSTS + 1; i0++) {
"""
    # Step 7: Generate C code for all integrands and the area density.
    # enable_fd_codegen=True tells c_codegen to automatically handle all
    # finite difference derivatives of gridfunctions.
    body += ccg.c_codegen(
        sympy_expressions,
        integrand_c_vars,
        enable_fd_codegen=True,
        enable_fd_functions=enable_fd_functions,
    )

    body += r"""
                // The differential area element, excluding coordinate steps (dθ, dφ)
                const REAL dA_unscaled = area_density * weight1 * weight2;
                
                // Accumulate into thread-private variables
                A_sum_private += A_integrand * dA_unscaled;
                XU_sum_private[0] += XU0_integrand * dA_unscaled;
                XU_sum_private[1] += XU1_integrand * dA_unscaled;
                XU_sum_private[2] += XU2_integrand * dA_unscaled;
                XRU_sum_private[0] += XRU0_integrand * dA_unscaled;
                XRU_sum_private[1] += XRU1_integrand * dA_unscaled;
                XRU_sum_private[2] += XRU2_integrand * dA_unscaled;
                XOU_sum_private[0] += XOU0_integrand * dA_unscaled;
                XOU_sum_private[1] += XOU1_integrand * dA_unscaled;
                XOU_sum_private[2] += XOU2_integrand * dA_unscaled;
                ZOU_sum_private[0] += ZOU0_integrand * dA_unscaled;
                ZOU_sum_private[1] += ZOU1_integrand * dA_unscaled;
                ZOU_sum_private[2] += ZOU2_integrand * dA_unscaled;
                R0_sum_private   += R0_integrand * dA_unscaled;
                O0_sum_private   += O0_integrand * dA_unscaled;
                Oabs_sum_private += Oabs_integrand * dA_unscaled;
            } // END LOOP: for i0 over radial horizon-grid points
        } // END LOOP: for i1 over theta horizon-grid points
    } // END LOOP: for i2 over phi horizon-grid points

    // Use a critical section for the final reduction from private to shared sums
    #pragma omp critical
    {
        A_sum += A_sum_private;
        for(int i=0; i<3; ++i) {
            XU_sum[i]  += XU_sum_private[i];
            XRU_sum[i] += XRU_sum_private[i];
            XOU_sum[i] += XOU_sum_private[i];
            ZOU_sum[i] += ZOU_sum_private[i];
        }
        R0_sum   += R0_sum_private;
        O0_sum   += O0_sum_private;
        Oabs_sum += Oabs_sum_private;
    } // END OMP CRITICAL: update shared spin-diagnostic sums
} // END OMP PARALLEL: precompute and integrate spin diagnostic

// Step 8: Compute the spin vector from the integrated quantities.
bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;

const REAL A = A_sum * dxx1 * dxx2;
const REAL R0 = R0_sum * dxx1 * dxx2;
const REAL XRU[3] = {XRU_sum[0] * dxx1 * dxx2, XRU_sum[1] * dxx1 * dxx2, XRU_sum[2] * dxx1 * dxx2};
const REAL XOU[3] = {XOU_sum[0] * dxx1 * dxx2, XOU_sum[1] * dxx1 * dxx2, XOU_sum[2] * dxx1 * dxx2};
const REAL ZOU[3] = {ZOU_sum[0] * dxx1 * dxx2, ZOU_sum[1] * dxx1 * dxx2, ZOU_sum[2] * dxx1 * dxx2};

const REAL spin_U[3] = {
    (XOU[0] - ZOU[0] * (XRU[0] / R0)) / A,
    (XOU[1] - ZOU[1] * (XRU[1] / R0)) / A,
    (XOU[2] - ZOU[2] * (XRU[2] / R0)) / A
};

bhahaha_diags->spin_chi_x_spectre = spin_U[0];
bhahaha_diags->spin_chi_y_spectre = spin_U[1];
bhahaha_diags->spin_chi_z_spectre = spin_U[2];

return BHAHAHA_SUCCESS;
"""
    # Format and register the C function using the standard helper.
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=cfunc_name,
        params=params,
        include_CodeParameters_h=False,
        body=clang_format(body),
    )
    return pcg.NRPyEnv()
