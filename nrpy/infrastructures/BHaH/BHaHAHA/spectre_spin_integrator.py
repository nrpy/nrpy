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
from nrpy.equations.general_relativity.bhahaha.SpECTRESpinEstimate import (
    SpECTRESpinEstimate,
)
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH.CurviBoundaryConditions.apply_bcs_inner_only import (
    APPLY_PARITY_BRANCHLESS_PREFUNC,
)


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
    :raises ValueError: If a precompute gridfunction has an unsupported rank.

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
    gf_macros = [f"{gf_name.upper()}GF" for gf_name in gf_names]

    parity_conditions_rank2 = {
        (0, 0): 4,
        (0, 1): 5,
        (0, 2): 6,
        (1, 1): 7,
        (1, 2): 8,
        (2, 2): 9,
    }
    parity_entries = []
    for gf_name, gf_macro in zip(gf_names, gf_macros):
        gf = gri.glb_gridfcs_dict[gf_name]
        if gf.rank == 0:
            parity_value = 0
        elif gf.rank == 1:
            parity_value = int(gf.name[-1]) + 1
        elif gf.rank == 2:
            parity_value = parity_conditions_rank2[(int(gf.name[-2]), int(gf.name[-1]))]
        else:
            raise ValueError(
                f"Unsupported spin-diagnostic precompute gridfunction rank: {gf.name}, rank={gf.rank}"
            )
        parity_entries.append(f"  [{gf_macro}] = {parity_value},")
    parity_table_entries = "\n".join(parity_entries)
    selected_precompute_gfs = ", ".join(gf_macros)

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

    prefunc = APPLY_PARITY_BRANCHLESS_PREFUNC + rf"""
static const int8_t spectre_spin_auxevol_gf_parity[NUM_AUXEVOL_GFS] = {{
{parity_table_entries}
}};

/**
 * Apply inner boundary conditions to selected AUXEVOL gridfunctions.
 *
 * The spin diagnostic precomputes SE_qDD and SE_XD on physical angular points
 * before generated finite-difference code differentiates those fields. This
 * helper fills their angular ghost zones for every active radial horizon slab.
 *
 * @param[in] bcstruct Boundary metadata for inner curvilinear points.
 * @param[in,out] auxevol_gfs AUXEVOL gridfunction storage.
 * @param Nxx0 Number of active radial horizon slabs.
 * @param Nxx_plus_2NGHOSTS0 Radial grid size including ghost zones.
 * @param Nxx_plus_2NGHOSTS1 Theta grid size including ghost zones.
 * @param Nxx_plus_2NGHOSTS2 Phi grid size including ghost zones.
 * @param[in] which_gfs Selected AUXEVOL gridfunction indices.
 * @param num_gfs Number of selected gridfunctions.
 * @param[in] auxevol_gf_parity Parity table indexed by AUXEVOL gridfunction.
 */
static void apply_inner_bc_for_selected_auxevol_gfs(const bc_struct *restrict bcstruct, REAL *restrict auxevol_gfs, const int Nxx0,
                                                    const int Nxx_plus_2NGHOSTS0, const int Nxx_plus_2NGHOSTS1,
                                                    const int Nxx_plus_2NGHOSTS2, const int *restrict which_gfs, const int num_gfs,
                                                    const int8_t *restrict auxevol_gf_parity) {{
  const bc_info_struct *bc_info = &bcstruct->bc_info;

#pragma omp parallel for collapse(2)
  for (int gf_idx = 0; gf_idx < num_gfs; gf_idx++) {{
    for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {{
      const int which_gf = which_gfs[gf_idx];
      const innerpt_bc_struct *restrict bc = &bcstruct->inner_bc_array[pt];
      const int dstpt = bc->dstpt;
      const int dst_i0 = dstpt % Nxx_plus_2NGHOSTS0;

      if (dst_i0 < NGHOSTS || dst_i0 >= NGHOSTS + Nxx0)
        continue;

      const int8_t p = bc->parity[auxevol_gf_parity[which_gf]];
      auxevol_gfs[IDX4pt(which_gf, dstpt)] = apply_parity_branchless(auxevol_gfs[IDX4pt(which_gf, bc->srcpt)], p);
    }} // END LOOP: for pt over inner boundary points
  }} // END LOOP: for selected AUXEVOL gridfunctions
}} // END FUNCTION: apply_inner_bc_for_selected_auxevol_gfs
"""

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

    // Precompute SE_qDD and SE_XD only where generated finite-difference
    // stencils are valid.
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
        const REAL xx2 = xx[2][i2];
        for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
            const REAL xx1 = xx[1][i1];
            for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
"""
    body += precompute_c_code
    body += rf"""
            }} // END LOOP: for i0 over active radial horizon-grid slabs
        }} // END LOOP: for i1 over physical theta horizon-grid points
    }} // END LOOP: for i2 over physical phi horizon-grid points

    {{
        const int spectre_spin_precompute_gfs[{len(gf_macros)}] = {{{selected_precompute_gfs}}};
        apply_inner_bc_for_selected_auxevol_gfs(
            &griddata[grid].bcstruct, auxevol_gfs, Nxx0, Nxx_plus_2NGHOSTS0,
            Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
            spectre_spin_precompute_gfs, {len(gf_macros)}, spectre_spin_auxevol_gf_parity);
    }} // END BLOCK: fill SE_qDD and SE_XD ghost zones before differentiating them

#pragma omp parallel
{{
    // Private accumulators for each thread
    REAL A_sum_private = 0.0;
    REAL XU_sum_private[3] = {{0.0, 0.0, 0.0}};
    REAL R0_sum_private = 0.0;
    REAL XRU_sum_private[3] = {{0.0, 0.0, 0.0}};
    REAL O0_sum_private = 0.0;
    REAL XOU_sum_private[3] = {{0.0, 0.0, 0.0}};
    REAL ZOU_sum_private[3] = {{0.0, 0.0, 0.0}};
    REAL Oabs_sum_private = 0.0;

#pragma omp for
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {{
        const REAL weight2 = weights[(i2 - NGHOSTS) % weight_stencil_size];
        const REAL xx2 = xx[2][i2];
        for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {{
            const REAL weight1 = weights[(i1 - NGHOSTS) % weight_stencil_size];
            const REAL xx1 = xx[1][i1];
            for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {{
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
            } // END LOOP: for i0 over active radial horizon-grid slabs
        } // END LOOP: for i1 over physical theta horizon-grid points
    } // END LOOP: for i2 over physical phi horizon-grid points

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
} // END OMP PARALLEL: integrate spin diagnostic

// Step 8: Compute the spin vector from the integrated quantities.
bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;

const REAL surface_weight = dxx1 * dxx2;
const REAL spin_norm_tolerance = 1.0e-14;
const REAL A = A_sum * surface_weight;
const REAL XU[3] = {
    XU_sum[0] * surface_weight,
    XU_sum[1] * surface_weight,
    XU_sum[2] * surface_weight};
const REAL R0 = R0_sum * surface_weight;
const REAL XRU[3] = {
    XRU_sum[0] * surface_weight,
    XRU_sum[1] * surface_weight,
    XRU_sum[2] * surface_weight};
const REAL O0 = O0_sum * surface_weight;
const REAL XOU[3] = {
    XOU_sum[0] * surface_weight,
    XOU_sum[1] * surface_weight,
    XOU_sum[2] * surface_weight};
const REAL ZOU[3] = {
    ZOU_sum[0] * surface_weight,
    ZOU_sum[1] * surface_weight,
    ZOU_sum[2] * surface_weight};
const REAL Oabs = Oabs_sum * surface_weight;

REAL spin_U[3] = {0.0, 0.0, 0.0};

if (fabs(A) > spin_norm_tolerance) {
    REAL x0U[3], xRcorrU[3], IU[3], SalphaU[3];

    for (int i = 0; i < 3; i++) {
        x0U[i] = XU[i] / A;
        xRcorrU[i] = (XRU[i] - x0U[i] * R0) / (8.0 * M_PI);
        IU[i] = XOU[i] - (x0U[i] + xRcorrU[i]) * O0;
        SalphaU[i] = ZOU[i] / (8.0 * M_PI);
    } // END LOOP: for i over spatial spin-vector components

    const REAL normI = sqrt(IU[0] * IU[0] + IU[1] * IU[1] + IU[2] * IU[2]);
    const REAL Salpha_norm = sqrt(SalphaU[0] * SalphaU[0] + SalphaU[1] * SalphaU[1] + SalphaU[2] * SalphaU[2]);
    const REAL S = Salpha_norm;

    if (Oabs > spin_norm_tolerance && normI > spin_norm_tolerance) {
        for (int i = 0; i < 3; i++)
            spin_U[i] = S * IU[i] / normI;
    } else if (Salpha_norm > spin_norm_tolerance) {
        for (int i = 0; i < 3; i++)
            spin_U[i] = S * SalphaU[i] / Salpha_norm;
    } // END ELSE IF: use z_alpha fallback direction when nominal direction is unsafe
} // END IF: area is safe for centroid reduction

bhahaha_diags->spin_chi_x_spectre = spin_U[0];
bhahaha_diags->spin_chi_y_spectre = spin_U[1];
bhahaha_diags->spin_chi_z_spectre = spin_U[2];

return BHAHAHA_SUCCESS;
"""
    # Format and register the C function using the standard helper.
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=cfunc_name,
        params=params,
        include_CodeParameters_h=False,
        body=clang_format(body),
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
