"""
This module registers a C function to perform numerical integrations
over the apparent horizon surface for spin and mass diagnostics.
"""

import sympy as sp
from typing import Union

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
import nrpy.helpers.parallel_codegen as pcg
from nrpy.helpers.generic import clang_format

# Import the SpECTRESpinEstimate factory from the other module
from nrpy.equations.general_relativity.bhahaha.SpECTRESpinEstimate import (
    SpECTRESpinEstimate,
)


def register_CFunction_SpECTRE_diagnostics_integration(
  CoordSystem: str = "Spherical",
  enable_rfm_precompute: bool = False,
  enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
  """
  Register a C function that integrates SpECTRE-style spin integrands over
  the apparent-horizon 2-surface and stores RunSums in the diagnostics struct.

  This follows the same registration pattern as other BHaH diagnostics.
  """
  if pcg.pcg_registration_phase():
    pcg.register_func_call(f"{__name__}.{register_CFunction_SpECTRE_diagnostics_integration.__name__}", locals())
    return None

  includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
  desc = r"""
  Perform all surface integrations for the SpECTRE-style spin diagnostic.
  This function computes the raw RunSums and stores them. A separate
  function should be called to reduce these sums to a final spin vector.
  """
  cfunc_type = "void"
  name = "BHaHAHA_SpECTRE_diagnostics_integration"
  params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata"

  # Step 1: Get an instance of the symbolic calculator.
  suffix = ("_rfm_precompute" if enable_rfm_precompute else "")
  spin_calc = SpECTRESpinEstimate[CoordSystem + suffix]

  # Step 2: Retrieve the dictionary of all per-point integrands.
  integrands_dict = spin_calc.get_public_integrands()

    # Step 3: Extract the symbolic expressions we need to integrate.
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
      "const REAL XU0_integrand", "const REAL XU1_integrand", "const REAL XU2_integrand",
      "const REAL R0_integrand",
      "const REAL XRU0_integrand", "const REAL XRU1_integrand", "const REAL XRU2_integrand",
      "const REAL O0_integrand",
      "const REAL XOU0_integrand", "const REAL XOU1_integrand", "const REAL XOU2_integrand",
      "const REAL ZOU0_integrand", "const REAL ZOU1_integrand", "const REAL ZOU2_integrand",
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

  # Step 4: Construct the body of the C function.
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
    for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
      const REAL weight1 = weights[(i1 - NGHOSTS) % weight_stencil_size];
      // The horizon is a 2D surface, so we loop over a single radial index i0.
      // All quantities are evaluated on this surface.
      for (int i0 = NGHOSTS; i0 < NGHOSTS + 1; i0++) {
"""
  # Step 5: Generate C code for all integrands and the area density.
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
      for(int i=0; i<3; ++i) {
          XU_sum_private[i]  += (&XU0_integrand)[i] * dA_unscaled;
          XRU_sum_private[i] += (&XRU0_integrand)[i] * dA_unscaled;
          XOU_sum_private[i] += (&XOU0_integrand)[i] * dA_unscaled;
          ZOU_sum_private[i] += (&ZOU0_integrand)[i] * dA_unscaled;
      }
      R0_sum_private   += R0_integrand * dA_unscaled;
      O0_sum_private   += O0_integrand * dA_unscaled;
      Oabs_sum_private += Oabs_integrand * dA_unscaled;
    } // END LOOP i0
  } // END LOOP i1
} // END LOOP i2

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
} // END OMP CRITICAL
} // END OMP PARALLEL

// Step 6: Finalize sums by multiplying by coordinate steps and store results.
bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;
const REAL dxx1 = params->dxx1;
const REAL dxx2 = params->dxx2;

bhahaha_diags->RunSums.A = A_sum * dxx1 * dxx2;
for(int i=0; i<3; ++i) {
    bhahaha_diags->RunSums.XU[i]  = XU_sum[i]  * dxx1 * dxx2;
    bhahaha_diags->RunSums.XRU[i] = XRU_sum[i] * dxx1 * dxx2;
    bhahaha_diags->RunSums.XOU[i] = XOU_sum[i] * dxx1 * dxx2;
    bhahaha_diags->RunSums.ZOU[i] = ZOU_sum[i] * dxx1 * dxx2;
}
bhahaha_diags->RunSums.R0   = R0_sum   * dxx1 * dxx2;
bhahaha_diags->RunSums.O0   = O0_sum   * dxx1 * dxx2;
bhahaha_diags->RunSums.Oabs = Oabs_sum * dxx1 * dxx2;

return; // Success
"""
  # Format and register the C function using the standard helper.
  cfc.register_CFunction(
    subdirectory="",
    includes=includes,
    desc=desc,
    cfunc_type=cfunc_type,
    name=name,
    params=params,
    include_CodeParameters_h=False,
    body=clang_format(body),
  )
  return pcg.NRPyEnv()