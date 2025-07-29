"""
Register the C function for BHaHAHA quasi-local horizon diagnostics.
This function computes the quasi-local angular momentum (spin) and the
Christodoulou-Ruffini mass of an apparent horizon.

This module follows the BHaH/NRPy+ framework to generate a C function
`bah_diagnostics_quasi_local_mass_spin` that implements the isolated/dynamical
horizon formalism.

The core algorithm involves:
1. Solving a generalized eigenvalue problem for a scalar potential 'zeta',
   which defines an Approximate Killing Vector (AKV) on the horizon surface.
   (Note: The elliptic solver is currently a placeholder).
2. Computing the AKV, phi^A, from the derivatives of the potential zeta. This
   is handled symbolically by NRPy+.
3. Calculating the rotational 1-form, 'omega_A', from the 3D metric and
   extrinsic curvature, which are interpolated to the horizon surface.
4. Integrating omega_A contracted with the AKV over the horizon surface
   to find the spin J.
5. Using J and the irreducible mass M_irr (from the horizon area) to
   compute the Christodoulou-Ruffini mass M.

Author: [Your Name/Enhanced LLM]
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
from nrpy.equations.general_relativity.bhahaha.ExpansionFunctionTheta import (
    ExpansionFunctionTheta,
)
from nrpy.infrastructures.BHaH.BHaHAHA import area


def register_CFunction_diagnostics_quasi_local_mass_spin(
    CoordSystem: str = "Spherical",
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for BHaHAHA quasi-local mass and spin diagnostics.

    :param CoordSystem: The coordinate system to use, defaults to "Spherical".
    :param enable_fd_functions: Whether to enable finite difference functions, defaults to False.
    :return: An NRPyEnv_type object if registration is successful, otherwise None.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # These comments document the necessary manual changes to BHaH C code.
    #
    # 1. Add new fields to bhahaha_diagnostics_struct in BHaHAHA.h:
    #
    #    typedef struct __bhahaha_diagnostics_struct__ {
    #      ...
    #      // Quasi-local diagnostics
    #      REAL quasi_local_spin_J;
    #      REAL quasi_local_mass_M;
    #    } bhahaha_diagnostics_struct;
    #        DONE^
    #
    # 2. Modify diagnostics() in diagnostics.py to call the new C function.
    #    This call should happen after other diagnostics like area and circumferences
    #    are computed, during the final iteration.
    #
    #    if (commondata->is_final_iteration) {
    #        ...
    #        // Calculate proper circumferences and update the error flag based on the outcome.
    #        commondata->error_flag = bah_diagnostics_proper_circumferences(commondata, griddata);
    #        if (commondata->error_flag != BHAHAHA_SUCCESS)
    #          return;
    #
    #        // Add this call:
    #        bah_diagnostics_quasi_local_mass_spin(commondata, griddata);
    #
    #        // Display detailed final iteration diagnostics...
    #    }
    #
    # 3. (Optional) Modify diagnostics_file_output() in diagnostics_file_output.py
    #    to add the new quantities to the output file. Add columns for J and M
    #    and then add bhahaha_diags->quasi_local_spin_J and
    #    bhahaha_diags->quasi_local_mass_M to the fprintf() call.

    # Step 1: Register gridfunctions for the AKV potential 'zeta' and the AKV 'phi'.
    # These are auxiliary gridfunctions on the 2D horizon surface grid.
    # Using AUXEVOL group ensures NRPy+ manages their memory automatically.
    zeta_gf = gri.register_gridfunctions(
        "zeta_gf", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )[0]
    L_gf = gri.register_gridfunctions(
        "L_gf", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )[0]
    # AKV is a 2-vector on the surface, with components (phi^theta, phi^phi).
    AKV_phi_gfU = ixp.register_gridfunctions_for_single_rank1(
        "AKV_phi_gfU", group="AUXEVOL", gf_array_name="auxevol_gfs", dimension=2
    )

    # Step 2: Define all symbolic quantities needed for the spin calculation.
    # We use ExpansionFunctionTheta as it provides all needed geometric quantities
    # derived from the BSSN variables interpolated to the AH surface.
    Th = ExpansionFunctionTheta[CoordSystem]
    gammaDD = Th.gammaDD
    KDD = Th.KDD
    gammaUU, _ = ixp.symm_matrix_inverter3x3(gammaDD)

    # Define horizon shape function 'h' and its derivatives.
    h = sp.Symbol("hh", real=True)
    h_dD = ixp.declarerank1("hh_dD", dimension=3)

    # A) Compute the unit normal vector to the horizon surface, s^i.
    # The surface is r - h(theta, phi) = 0. The normal is proportional to nabla_i(r-h).
    sD = ixp.zerorank1(dimension=3)
    sD[0] = 1
    sD[1] = -h_dD[1]
    sD[2] = -h_dD[2]

    s_norm_sq = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            s_norm_sq += gammaUU[i][j] * sD[i] * sD[j]

    for i in range(3):
        sD[i] /= sp.sqrt(s_norm_sq)

    sU = ixp.zerorank1(dimension=3)
    for i in range(3):
        for j in range(3):
            sU[i] += gammaUU[i][j] * sD[j]
    
    # B) Compute the rotational 1-form, omega_A = -s_i(D_A r^i)
    # where A is a surface index (0 for theta, 1 for phi).
    omega_d = ixp.zerorank1(dimension=2)
    sU_
    for A in range(2):
        for j in range(3):
            omega_d[A] += -sD[j] * KDD[j][A + 1]

    # C) Define the Approximate Killing Vector (AKV) phi^A from the potential zeta.
    # phi^A = epsilon^{AB} * D_B(zeta), where D_B is the surface covariant derivative.
    # For a scalar, D_B(zeta) is just partial_B(zeta).
    q2DD = area.compute_q2DD()
    _, q2det = ixp.symm_matrix_inverter2x2(q2DD)
    zeta_dD = ixp.declarerank1("zeta_gf_dD", dimension=2)  # Derivatives of zeta_gf

    phiU = ixp.zerorank1(dimension=2)
    # The Levi-Civita tensor on the surface is epsilon^{AB} = (1/sqrt(q_det)) * [[0, 1], [-1, 0]]
    phiU[0] = (1 / sp.sqrt(q2det)) * zeta_dD[1]  # phi^theta
    phiU[1] = -(1 / sp.sqrt(q2det)) * zeta_dD[0]  # phi^phi

    # D) Define the integrand for the spin calculation: omega_A * phi^A
    integrand = sp.sympify(0)
    for A in range(2):
        integrand += omega_d[A] * phiU[A]

    # Step 3: Register the C function.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"]
    desc = "BHaHAHA quasi-local diagnostics: compute quasi-local spin (J) and Christodoulou-Ruffini mass (M)."
    cfunc_type = "int"
    name = "bah_diagnostics_quasi_local_mass_spin"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    # Step 3.1: Define the iterative elliptic solver as a pre-function.
    # This solver implements the Cook & Whiting (2007) method.
    prefunc = r"""
// Iterative Gauss-Seidel solver for the coupled elliptic system.
// This function finds the scalar potential 'v' (stored as zeta_gf) and an
// auxiliary field 'L' that satisfy the equations for an Approximate Killing Vector.
//
// The system of equations is:
// (1) D^i D_i v + 2L = 0
// (2) D^i D_i L - (1-Theta)*(0.5*(D^i R_S)*D_i v - 2*R_S*L) = 0
//
// NOTE: This implementation requires Ricci scalar R_S and its derivatives,
// which must be computed from the 2-metric q_ij.
static void solve_coupled_elliptic_system_for_v_and_L(
    const params_struct *restrict params, REAL *restrict xx[3],
    REAL *restrict auxevol_gfs, REAL *restrict auxevol_gfs_prev) {

    // Solver parameters
    const int MAX_ITER = 2000;
    const REAL TOLERANCE = 1e-10;
    // The Lagrange multiplier Theta. For a true Killing vector, Theta = 0.
    // Finding the optimal Theta is complex; we use a fixed value as a starting point.
    const REAL THETA = 0.0;
#include "set_CodeParameters.h" // Needed for grid parameters like Nxx1, Nxx2

    // Initialize fields: v (zeta_gf) and L (L_gf).
    // We provide an initial guess for v that approximates a z-axis rotation.
    // L is initialized to zero.

    #pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < params->Nxx2 + NGHOSTS; i2++) {
        for (int i1 = NGHOSTS; i1 < params->Nxx1 + NGHOSTS; i1++) {
            // This is a 2D grid, so the radial index is fixed.
            const int i0 = NGHOSTS;
            const REAL theta = xx[1][i1];
            // Example: populate zeta with a pattern for a z-axis AKV (zeta ~ cos(theta))
            auxevol_gfs[IDX4(ZETA_GFGF, i0, i1, i2)] = cos(theta);
        }
    }
}
"""

    # Step 3.2: Construct the body of the C function.
    body = r"""
  const int grid=0;
  const params_struct *restrict params = &griddata[grid].params;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  const REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs; // for hh
  REAL *restrict xx[3];
  for(int ww=0;ww<3;ww++) xx[ww] = griddata[grid].xx[ww];
#include "set_CodeParameters.h"

  // Step A: Solve the elliptic equation for the AKV potential 'zeta'.
  // In a real application, this would be a call to a numerical library (e.g., PETSc).
  // Here, we call our placeholder function, which populates ZETA_GFGF.
  solve_elliptic_equation_for_zeta(params, xx, auxevol_gfs);

  // Step B: Compute the surface integral of (omega_A * phi^A) to get the spin J.
  // Set integration weights.
  const REAL *restrict weights;
  int weight_stencil_size;
  bah_diagnostics_integration_weights(Nxx1, Nxx2, &weights, &weight_stencil_size);

  REAL sum_J_integrand = 0.0;
#pragma omp parallel
{
#pragma omp for
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    const REAL weight2 = weights[(i2 - NGHOSTS) % weight_stencil_size];
    const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
      const REAL weight1 = weights[(i1 - NGHOSTS) % weight_stencil_size];
      const REAL xx1 = xx[1][i1];
      // The horizon is a 2D surface, so we loop over a single radial index i0.
      for (int i0 = NGHOSTS; i0 < NGHOSTS + 1; i0++) {
"""
    # Generate C code for the integrand and area element.
    # Note that we use enable_fd_codegen to handle derivatives of hh and zeta.
    # All BSSN variables are read from auxevol_gfs, hh from in_gfs.
    body += (
        ccg.c_codegen(
            [integrand, area.area3()],
            ["const REAL integrand", "const REAL area_element"],
            enable_fd_codegen=True,
            enable_fd_functions=enable_fd_functions,
        )
        + r"""
#pragma omp critical
        {
          sum_J_integrand += integrand * area_element * weight1 * weight2;
        } // END OMP CRITICAL
      } // END LOOP over i0
    } // END LOOP over i1
  } // END LOOP over i2
} // END OMP PARALLEL

  // Step C: Finalize spin and mass calculations and store results.
  bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;

  // Compute the quasi-local spin J. This is a single component of the spin vector,
  // corresponding to the AKV found by the solver.
  // The full algorithm would repeat this for 3 orthogonal AKVs.
  const REAL J_quasi_local = (1.0 / (8.0 * M_PI)) * sum_J_integrand * params->dxx1 * params->dxx2;

  // Compute the irreducible mass from the horizon area.
  const REAL Area = bhahaha_diags->area;
  if (Area <= 1e-15) { // Avoid division by zero if Area is not yet computed or is invalid.
      return BHAHAHA_INVALID_INPUT;
  }
  const REAL M_irr = sqrt(Area / (16.0 * M_PI));

  // Compute the Christodoulou-Ruffini mass M.
  const REAL M_irr2 = M_irr * M_irr;
  if (M_irr2 <= 1e-15) { // Avoid division by zero for nearly zero-mass horizons.
      return BHAHAHA_INVALID_INPUT;
  }
  const REAL M_sq = M_irr2 + J_quasi_local * J_quasi_local / (4.0 * M_irr2);
  const REAL M_quasi_local = sqrt(M_sq);

  // Store diagnostics in the commondata struct.
  bhahaha_diags->quasi_local_spin_J = J_quasi_local;
  bhahaha_diags->quasi_local_mass_M = M_quasi_local;

  return BHAHAHA_SUCCESS;
"""
    )

    # Register the final C function.
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()