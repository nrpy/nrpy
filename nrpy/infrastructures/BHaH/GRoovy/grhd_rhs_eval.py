"""
C function responsible for computing the right-hand sides of the GRHD equations.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_grhd_rhs_eval(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    evolving_entropy: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function for complete GRHD right hand side evaluation.

    This function orchestrates the evaluation of partial_t U = -partial_i F^i + S.
    It calculates source terms (S), then loops over spatial directions to:
    1. Interpolate metric terms to cell faces.
    2. Reconstruct primitive variables to cell faces (Left/Right states).
    3. Compute u^0 at cell faces.
    4. Solve the Riemann problem (HLL) to get fluxes (F^i).
    5. Compute flux divergences (-partial_i F^i) and add to the RHS.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param evolving_entropy: whether we're using GRHayL to evolve entropy or not

    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up C function headers, signature, and parameters
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    cfunc_type = "void"
    name = "grhd_rhs_eval"
    desc = "Calculate RHSs of GRHD equations"
    params = r"""
const commondata_struct *restrict commondata,
const params_struct *restrict params,
const ghl_parameters *restrict ghl_params,
const ghl_eos_parameters *restrict eos,
REAL *restrict xx[3],
const REAL *restrict evol_gfs,
REAL *restrict auxevol_gfs,
REAL *restrict rhs_gfs"""

    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )

    # Step 3: Source Terms
    # Calculate source terms (S) for all equations and initialize RHS gridfunctions.
    body = r"""

    // Function pointer to allow for loop over fluxes
  void (*calculate_HLL_fluxes)(const commondata_struct *, const params_struct *,
  const ghl_eos_parameters *restrict, REAL *) = NULL;

  calculate_all_source_terms(commondata, params, xx, eos, auxevol_gfs, evol_gfs, rhs_gfs);

  //Loop over flux directions (x, y, z)
  for(int flux_dir=0; flux_dir<3; flux_dir++) {
    //Set function pointer to specific HLL flux function for the current direction
    switch(flux_dir) {
      case 0:
        calculate_HLL_fluxes = &calculate_HLL_fluxes_direction_0;
        break;
      case 1:
        calculate_HLL_fluxes = &calculate_HLL_fluxes_direction_1;
        break;
      case 2:
        calculate_HLL_fluxes = &calculate_HLL_fluxes_direction_2;
        break;
    }

    //Interpolate metric quantities to cell faces
    interpolate_metric_gfs_to_cell_faces(commondata,
                                         params,
                                         flux_dir,
                                         evol_gfs,
                                         auxevol_gfs);


    //Reconstruct primitive variables (rho, P, etc.,) to cell faces.
    //          This will be mc, PPM or WenoZ.
    reconstruction_loop(commondata,
                        params,
                        ghl_params,
                        flux_dir,
                        eos,
                        auxevol_gfs);

"""
    if evolving_entropy:
        body += r"""
    //Reconstruct entropy if necessary
    reconstruction_entropy_loop(commondata,
                               params,
                               ghl_params,
                               flux_dir,
                               eos,
                               auxevol_gfs);
"""
    body += r"""

    //Compute u^0 (time component of 4-velocity) at cell faces
    //          We need this for the Riemann solver. We calculate it for both
    //          Right (R) and Left (L) states at the interface.
    #pragma omp parallel for
    LOOP_REGION(NGHOSTS, (Nxx_plus_2NGHOSTS0-NGHOSTS) + (flux_dir == 0),
                NGHOSTS, (Nxx_plus_2NGHOSTS1-NGHOSTS) + (flux_dir == 1),
                NGHOSTS, (Nxx_plus_2NGHOSTS2-NGHOSTS) + (flux_dir == 2)){

      const int index = IDX3(i0, i1, i2);

      // Compute u^0 for the Right state
      compute_up_index_velocity_time_component_pointwise(commondata, params, ghl_params,
      auxevol_gfs[IDX4pt(ALPHA_FACEGF, index)], auxevol_gfs[IDX4pt(VET_FACEU0GF, index)],
      auxevol_gfs[IDX4pt(VET_FACEU1GF, index)], auxevol_gfs[IDX4pt(VET_FACEU2GF, index)],
      auxevol_gfs[IDX4pt(H_FACEDD00GF, index)], auxevol_gfs[IDX4pt(H_FACEDD01GF, index)],
      auxevol_gfs[IDX4pt(H_FACEDD02GF, index)],
      auxevol_gfs[IDX4pt(H_FACEDD11GF, index)], auxevol_gfs[IDX4pt(H_FACEDD12GF, index)],
      auxevol_gfs[IDX4pt(H_FACEDD22GF, index)], auxevol_gfs[IDX4pt(CF_FACEGF, index)],
      &auxevol_gfs[IDX4pt(RESCALEDVRU0GF, index)], &auxevol_gfs[IDX4pt(RESCALEDVRU1GF, index)],
      &auxevol_gfs[IDX4pt(RESCALEDVRU2GF, index)], &auxevol_gfs[IDX4pt(U4RUTGF, index)]);

      // Compute u^0 for the Left state
      compute_up_index_velocity_time_component_pointwise(commondata, params, ghl_params,
      auxevol_gfs[IDX4pt(ALPHA_FACEGF, index)], auxevol_gfs[IDX4pt(VET_FACEU0GF, index)],
      auxevol_gfs[IDX4pt(VET_FACEU1GF, index)], auxevol_gfs[IDX4pt(VET_FACEU2GF, index)],
      auxevol_gfs[IDX4pt(H_FACEDD00GF, index)], auxevol_gfs[IDX4pt(H_FACEDD01GF, index)],
      auxevol_gfs[IDX4pt(H_FACEDD02GF, index)],
      auxevol_gfs[IDX4pt(H_FACEDD11GF, index)], auxevol_gfs[IDX4pt(H_FACEDD12GF, index)],
      auxevol_gfs[IDX4pt(H_FACEDD22GF, index)], auxevol_gfs[IDX4pt(CF_FACEGF, index)],
      &auxevol_gfs[IDX4pt(RESCALEDVLU0GF, index)], &auxevol_gfs[IDX4pt(RESCALEDVLU1GF, index)],
      &auxevol_gfs[IDX4pt(RESCALEDVLU2GF, index)], &auxevol_gfs[IDX4pt(U4LUTGF, index)]);
    }

    //Calculate HLL Fluxes at interfaces
    calculate_HLL_fluxes(commondata, params, eos, auxevol_gfs);

    //Calculate Flux Divergences and update RHS
    //          RHS -= \partial_i F^i
    calculate_flux_divergences(flux_dir, commondata, params, xx, eos, auxevol_gfs, evol_gfs, rhs_gfs);
  }
"""
    if enable_rfm_precompute:
        body = body.replace("xx,", "rfmstruct,")

    # Step 5: Register the final C function
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=body,
    )
    return pcg.NRPyEnv()
