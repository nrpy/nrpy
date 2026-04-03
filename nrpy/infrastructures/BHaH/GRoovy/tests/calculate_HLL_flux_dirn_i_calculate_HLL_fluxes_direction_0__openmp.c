#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Compute HLL fluxes on the 0-face for the GRHD evolution variables.
 */
void calculate_HLL_fluxes_direction_0(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                      const ghl_eos_parameters *restrict eos, REAL *restrict auxevol_gfs) {
#include "set_CodeParameters.h"
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS + 1; i2++) {
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS + 1; i1++) {
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS + 1; i0++) {

        // Initialize GRHayL primitive states on the right and left sides.
        ghl_primitive_quantities prims_r, prims_l;
        ghl_initialize_primitives(auxevol_gfs[IDX4(RHOB_RGF, i0, i1, i2)], auxevol_gfs[IDX4(P_RGF, i0, i1, i2)], NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                                  NAN, NAN, NAN, &prims_r);

        ghl_initialize_primitives(auxevol_gfs[IDX4(RHOB_LGF, i0, i1, i2)], auxevol_gfs[IDX4(P_LGF, i0, i1, i2)], NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                                  NAN, NAN, NAN, &prims_l);

        double h_r, h_l, cs2_r, cs2_l;
        ghl_compute_h_and_cs2(eos, &prims_r, &h_r, &cs2_r);
        ghl_compute_h_and_cs2(eos, &prims_l, &h_l, &cs2_l);
        /*
         *  Original SymPy expression:

         */
        {
        }

        /*
         * NRPy-Generated GF Access/FD Code, Step 1 of 1:
         * Evaluate SymPy expressions and write to main memory.
         */
        /*
         *  Original SymPy expressions:
         *  "[auxevol_gfs[IDX4(RHO_STAR_HLL_FLUXGF, i0, i1, i2)] = (alpha_face*rescaledvlU0*rhob_l*u4lUt*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 - (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*Abs(cf_face)/cf_face**4 +
         * alpha_face*rescaledvrU0*rhob_r*u4rUt*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 -
         * cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) - 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) +
         * cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 +
         * 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) -
         * vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*Abs(cf_face)/cf_face**4 -
         * (-alpha_face*rhob_l*u4lUt*Abs(cf_face)/cf_face**4 + alpha_face*rhob_r*u4rUt*Abs(cf_face)/cf_face**4)*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1
         * - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 - (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))))/(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2)]"
         *  "[auxevol_gfs[IDX4(TAU_TILDE_HLL_FLUXGF, i0, i1, i2)] = ((-alpha_face*rescaledvlU0*rhob_l*u4lUt*Abs(cf_face)/cf_face**4 + (P_l*vet_faceU0
         * + alpha_face**2*h_l*rescaledvlU0*rhob_l*u4lUt**2)*Abs(cf_face)/cf_face**4)*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 - (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))) +
         * (-alpha_face*rescaledvrU0*rhob_r*u4rUt*Abs(cf_face)/cf_face**4 + (P_r*vet_faceU0 +
         * alpha_face**2*h_r*rescaledvrU0*rhob_r*u4rUt**2)*Abs(cf_face)/cf_face**4)*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))) - (alpha_face*rhob_l*u4lUt*Abs(cf_face)/cf_face**4 -
         * alpha_face*rhob_r*u4rUt*Abs(cf_face)/cf_face**4 - (-P_l + alpha_face**2*h_l*rhob_l*u4lUt**2)*Abs(cf_face)/cf_face**4 + (-P_r +
         * alpha_face**2*h_r*rhob_r*u4rUt**2)*Abs(cf_face)/cf_face**4)*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 - (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))))/(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2)]"
         *  "[auxevol_gfs[IDX4(RESCALEDSTILDE_FLUX_HLLD0GF, i0, i1, i2)] =
         * (alpha_face*(h_faceDD01*(P_l*(alpha_face**2*cf_face**2*(-h_faceDD01*(h_faceDD22 + 1) + h_faceDD02*h_faceDD12) -
         * vet_faceU0*vet_faceU1*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_l*rescaledvlU0*rescaledvlU1*rhob_l*u4lUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * h_faceDD02*(P_l*(alpha_face**2*cf_face**2*(h_faceDD01*h_faceDD12 - h_faceDD02*(h_faceDD11 + 1)) -
         * vet_faceU0*vet_faceU2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_l*rescaledvlU0*rescaledvlU2*rhob_l*u4lUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (h_faceDD00 +
         * 1)*(P_l*(alpha_face**2*cf_face**2*(-h_faceDD12**2 + (h_faceDD11 + 1)*(h_faceDD22 + 1)) - vet_faceU0**2*(-h_faceDD01**2*(h_faceDD22 + 1) +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 +
         * 1)*(h_faceDD22 + 1))) + alpha_face**2*h_l*rescaledvlU0**2*rhob_l*u4lUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 +
         * 1)*(h_faceDD22 + 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (P_l*vet_faceU0 +
         * alpha_face**2*h_l*rescaledvlU0*rhob_l*u4lUt**2)*(h_faceDD01*vet_faceU1 + h_faceDD02*vet_faceU2 + vet_faceU0*(h_faceDD00 +
         * 1))/(alpha_face**2*cf_face**2))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l)
         * + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) - 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) +
         * cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 +
         * 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) -
         * vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 - (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*Abs(cf_face)/cf_face**4 +
         * alpha_face*(h_faceDD01*(P_r*(alpha_face**2*cf_face**2*(-h_faceDD01*(h_faceDD22 + 1) + h_faceDD02*h_faceDD12) -
         * vet_faceU0*vet_faceU1*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_r*rescaledvrU0*rescaledvrU1*rhob_r*u4rUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * h_faceDD02*(P_r*(alpha_face**2*cf_face**2*(h_faceDD01*h_faceDD12 - h_faceDD02*(h_faceDD11 + 1)) -
         * vet_faceU0*vet_faceU2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_r*rescaledvrU0*rescaledvrU2*rhob_r*u4rUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (h_faceDD00 +
         * 1)*(P_r*(alpha_face**2*cf_face**2*(-h_faceDD12**2 + (h_faceDD11 + 1)*(h_faceDD22 + 1)) - vet_faceU0**2*(-h_faceDD01**2*(h_faceDD22 + 1) +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 +
         * 1)*(h_faceDD22 + 1))) + alpha_face**2*h_r*rescaledvrU0**2*rhob_r*u4rUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 +
         * 1)*(h_faceDD22 + 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (P_r*vet_faceU0 +
         * alpha_face**2*h_r*rescaledvrU0*rhob_r*u4rUt**2)*(h_faceDD01*vet_faceU1 + h_faceDD02*vet_faceU2 + vet_faceU0*(h_faceDD00 +
         * 1))/(alpha_face**2*cf_face**2))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l)
         * + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) - 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) +
         * cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 +
         * 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) -
         * vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*Abs(cf_face)/cf_face**4 -
         * (-alpha_face*(h_faceDD01*(P_l*vet_faceU1 + alpha_face**2*h_l*rescaledvlU1*rhob_l*u4lUt**2)/(alpha_face**2*cf_face**2) +
         * h_faceDD02*(P_l*vet_faceU2 + alpha_face**2*h_l*rescaledvlU2*rhob_l*u4lUt**2)/(alpha_face**2*cf_face**2) + (-P_l +
         * alpha_face**2*h_l*rhob_l*u4lUt**2)*(h_faceDD01*vet_faceU1 + h_faceDD02*vet_faceU2 + vet_faceU0*(h_faceDD00 + 1))/(alpha_face**2*cf_face**2)
         * + (h_faceDD00 + 1)*(P_l*vet_faceU0 + alpha_face**2*h_l*rescaledvlU0*rhob_l*u4lUt**2)/(alpha_face**2*cf_face**2))*Abs(cf_face)/cf_face**4 +
         * alpha_face*(h_faceDD01*(P_r*vet_faceU1 + alpha_face**2*h_r*rescaledvrU1*rhob_r*u4rUt**2)/(alpha_face**2*cf_face**2) +
         * h_faceDD02*(P_r*vet_faceU2 + alpha_face**2*h_r*rescaledvrU2*rhob_r*u4rUt**2)/(alpha_face**2*cf_face**2) + (-P_r +
         * alpha_face**2*h_r*rhob_r*u4rUt**2)*(h_faceDD01*vet_faceU1 + h_faceDD02*vet_faceU2 + vet_faceU0*(h_faceDD00 + 1))/(alpha_face**2*cf_face**2)
         * + (h_faceDD00 + 1)*(P_r*vet_faceU0 +
         * alpha_face**2*h_r*rescaledvrU0*rhob_r*u4rUt**2)/(alpha_face**2*cf_face**2))*Abs(cf_face)/cf_face**4)*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1
         * - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 - (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))))/(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2)]"
         *  "[auxevol_gfs[IDX4(RESCALEDSTILDE_FLUX_HLLD1GF, i0, i1, i2)] = (alpha_face*(h_faceDD01*(P_l*(alpha_face**2*cf_face**2*(-h_faceDD12**2 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)) - vet_faceU0**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_l*rescaledvlU0**2*rhob_l*u4lUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * h_faceDD12*(P_l*(alpha_face**2*cf_face**2*(h_faceDD01*h_faceDD12 - h_faceDD02*(h_faceDD11 + 1)) -
         * vet_faceU0*vet_faceU2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_l*rescaledvlU0*rescaledvlU2*rhob_l*u4lUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (h_faceDD11 +
         * 1)*(P_l*(alpha_face**2*cf_face**2*(-h_faceDD01*(h_faceDD22 + 1) + h_faceDD02*h_faceDD12) -
         * vet_faceU0*vet_faceU1*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_l*rescaledvlU0*rescaledvlU1*rhob_l*u4lUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (P_l*vet_faceU0 +
         * alpha_face**2*h_l*rescaledvlU0*rhob_l*u4lUt**2)*(h_faceDD01*vet_faceU0 + h_faceDD12*vet_faceU2 + vet_faceU1*(h_faceDD11 +
         * 1))/(alpha_face**2*cf_face**2))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l)
         * + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) - 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) +
         * cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 +
         * 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) -
         * vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 - (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*Abs(cf_face)/cf_face**4 +
         * alpha_face*(h_faceDD01*(P_r*(alpha_face**2*cf_face**2*(-h_faceDD12**2 + (h_faceDD11 + 1)*(h_faceDD22 + 1)) -
         * vet_faceU0**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_r*rescaledvrU0**2*rhob_r*u4rUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * h_faceDD12*(P_r*(alpha_face**2*cf_face**2*(h_faceDD01*h_faceDD12 - h_faceDD02*(h_faceDD11 + 1)) -
         * vet_faceU0*vet_faceU2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_r*rescaledvrU0*rescaledvrU2*rhob_r*u4rUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (h_faceDD11 +
         * 1)*(P_r*(alpha_face**2*cf_face**2*(-h_faceDD01*(h_faceDD22 + 1) + h_faceDD02*h_faceDD12) -
         * vet_faceU0*vet_faceU1*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_r*rescaledvrU0*rescaledvrU1*rhob_r*u4rUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (P_r*vet_faceU0 +
         * alpha_face**2*h_r*rescaledvrU0*rhob_r*u4rUt**2)*(h_faceDD01*vet_faceU0 + h_faceDD12*vet_faceU2 + vet_faceU1*(h_faceDD11 +
         * 1))/(alpha_face**2*cf_face**2))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l)
         * + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) - 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) +
         * cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 +
         * 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) -
         * vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*Abs(cf_face)/cf_face**4 -
         * (-alpha_face*(h_faceDD01*(P_l*vet_faceU0 + alpha_face**2*h_l*rescaledvlU0*rhob_l*u4lUt**2)/(alpha_face**2*cf_face**2) +
         * h_faceDD12*(P_l*vet_faceU2 + alpha_face**2*h_l*rescaledvlU2*rhob_l*u4lUt**2)/(alpha_face**2*cf_face**2) + (-P_l +
         * alpha_face**2*h_l*rhob_l*u4lUt**2)*(h_faceDD01*vet_faceU0 + h_faceDD12*vet_faceU2 + vet_faceU1*(h_faceDD11 + 1))/(alpha_face**2*cf_face**2)
         * + (h_faceDD11 + 1)*(P_l*vet_faceU1 + alpha_face**2*h_l*rescaledvlU1*rhob_l*u4lUt**2)/(alpha_face**2*cf_face**2))*Abs(cf_face)/cf_face**4 +
         * alpha_face*(h_faceDD01*(P_r*vet_faceU0 + alpha_face**2*h_r*rescaledvrU0*rhob_r*u4rUt**2)/(alpha_face**2*cf_face**2) +
         * h_faceDD12*(P_r*vet_faceU2 + alpha_face**2*h_r*rescaledvrU2*rhob_r*u4rUt**2)/(alpha_face**2*cf_face**2) + (-P_r +
         * alpha_face**2*h_r*rhob_r*u4rUt**2)*(h_faceDD01*vet_faceU0 + h_faceDD12*vet_faceU2 + vet_faceU1*(h_faceDD11 + 1))/(alpha_face**2*cf_face**2)
         * + (h_faceDD11 + 1)*(P_r*vet_faceU1 +
         * alpha_face**2*h_r*rescaledvrU1*rhob_r*u4rUt**2)/(alpha_face**2*cf_face**2))*Abs(cf_face)/cf_face**4)*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1
         * - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 - (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))))/(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2)]"
         *  "[auxevol_gfs[IDX4(RESCALEDSTILDE_FLUX_HLLD2GF, i0, i1, i2)] = (alpha_face*(h_faceDD02*(P_l*(alpha_face**2*cf_face**2*(-h_faceDD12**2 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)) - vet_faceU0**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_l*rescaledvlU0**2*rhob_l*u4lUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * h_faceDD12*(P_l*(alpha_face**2*cf_face**2*(-h_faceDD01*(h_faceDD22 + 1) + h_faceDD02*h_faceDD12) -
         * vet_faceU0*vet_faceU1*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_l*rescaledvlU0*rescaledvlU1*rhob_l*u4lUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (h_faceDD22 +
         * 1)*(P_l*(alpha_face**2*cf_face**2*(h_faceDD01*h_faceDD12 - h_faceDD02*(h_faceDD11 + 1)) - vet_faceU0*vet_faceU2*(-h_faceDD01**2*(h_faceDD22
         * + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11
         * + 1)*(h_faceDD22 + 1))) + alpha_face**2*h_l*rescaledvlU0*rescaledvlU2*rhob_l*u4lUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 +
         * 1)*(h_faceDD22 + 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (P_l*vet_faceU0 +
         * alpha_face**2*h_l*rescaledvlU0*rhob_l*u4lUt**2)*(h_faceDD02*vet_faceU0 + h_faceDD12*vet_faceU1 + vet_faceU2*(h_faceDD22 +
         * 1))/(alpha_face**2*cf_face**2))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l)
         * + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) - 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) +
         * cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 +
         * 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) -
         * vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 - (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*Abs(cf_face)/cf_face**4 +
         * alpha_face*(h_faceDD02*(P_r*(alpha_face**2*cf_face**2*(-h_faceDD12**2 + (h_faceDD11 + 1)*(h_faceDD22 + 1)) -
         * vet_faceU0**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_r*rescaledvrU0**2*rhob_r*u4rUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * h_faceDD12*(P_r*(alpha_face**2*cf_face**2*(-h_faceDD01*(h_faceDD22 + 1) + h_faceDD02*h_faceDD12) -
         * vet_faceU0*vet_faceU1*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) +
         * alpha_face**2*h_r*rescaledvrU0*rescaledvrU1*rhob_r*u4rUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) -
         * h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (h_faceDD22 +
         * 1)*(P_r*(alpha_face**2*cf_face**2*(h_faceDD01*h_faceDD12 - h_faceDD02*(h_faceDD11 + 1)) - vet_faceU0*vet_faceU2*(-h_faceDD01**2*(h_faceDD22
         * + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11
         * + 1)*(h_faceDD22 + 1))) + alpha_face**2*h_r*rescaledvrU0*rescaledvrU2*rhob_r*u4rUt**2*(-h_faceDD01**2*(h_faceDD22 + 1) +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12 - h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 +
         * 1)*(h_faceDD22 + 1)))/(alpha_face**2*cf_face**2*(-h_faceDD01**2*(h_faceDD22 + 1) + 2*h_faceDD01*h_faceDD02*h_faceDD12 -
         * h_faceDD02**2*(h_faceDD11 + 1) - h_faceDD12**2*(h_faceDD00 + 1) + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1))) + (P_r*vet_faceU0 +
         * alpha_face**2*h_r*rescaledvrU0*rhob_r*u4rUt**2)*(h_faceDD02*vet_faceU0 + h_faceDD12*vet_faceU1 + vet_faceU2*(h_faceDD22 +
         * 1))/(alpha_face**2*cf_face**2))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l)
         * + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) - 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) +
         * cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 +
         * 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) -
         * vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*Abs(cf_face)/cf_face**4 -
         * (-alpha_face*(h_faceDD02*(P_l*vet_faceU0 + alpha_face**2*h_l*rescaledvlU0*rhob_l*u4lUt**2)/(alpha_face**2*cf_face**2) +
         * h_faceDD12*(P_l*vet_faceU1 + alpha_face**2*h_l*rescaledvlU1*rhob_l*u4lUt**2)/(alpha_face**2*cf_face**2) + (-P_l +
         * alpha_face**2*h_l*rhob_l*u4lUt**2)*(h_faceDD02*vet_faceU0 + h_faceDD12*vet_faceU1 + vet_faceU2*(h_faceDD22 + 1))/(alpha_face**2*cf_face**2)
         * + (h_faceDD22 + 1)*(P_l*vet_faceU2 + alpha_face**2*h_l*rescaledvlU2*rhob_l*u4lUt**2)/(alpha_face**2*cf_face**2))*Abs(cf_face)/cf_face**4 +
         * alpha_face*(h_faceDD02*(P_r*vet_faceU0 + alpha_face**2*h_r*rescaledvrU0*rhob_r*u4rUt**2)/(alpha_face**2*cf_face**2) +
         * h_faceDD12*(P_r*vet_faceU1 + alpha_face**2*h_r*rescaledvrU1*rhob_r*u4rUt**2)/(alpha_face**2*cf_face**2) + (-P_r +
         * alpha_face**2*h_r*rhob_r*u4rUt**2)*(h_faceDD02*vet_faceU0 + h_faceDD12*vet_faceU1 + vet_faceU2*(h_faceDD22 + 1))/(alpha_face**2*cf_face**2)
         * + (h_faceDD22 + 1)*(P_r*vet_faceU2 +
         * alpha_face**2*h_r*rescaledvrU2*rhob_r*u4rUt**2)/(alpha_face**2*cf_face**2))*Abs(cf_face)/cf_face**4)*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1
         * - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))*(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/8 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/8 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 - (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(16*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(16*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))))/(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 - nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 -
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + nrpyAbs(nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/4 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/4 + nrpyAbs(-nrpyAbs((-2*rescaledvlU0*u4lUt**2*(1 -
         * cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(2*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + sqrt(-(-cs2_l*((-h_faceDD12**2/cf_face**4 +
         * (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2)/2 +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2/2 + nrpyAbs(-(-cs2_l*((-h_faceDD12**2/cf_face**4 + (h_faceDD11
         * + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 -
         * h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 + (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 +
         * 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvlU0**2*u4lUt**2*(1 - cs2_l))*(4*u4lUt**2*(1 - cs2_l) + 4*cs2_l/alpha_face**2) +
         * (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) + 2*cs2_l*vet_faceU0/alpha_face**2)**2)/2)/(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2))/2 +
         * nrpyAbs((-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) - 2*cs2_r*vet_faceU0/alpha_face**2)/(2*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) +
         * sqrt(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6 +
         * 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2)/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) + 2*cs2_r*vet_faceU0/alpha_face**2)**2/2
         * + nrpyAbs(-(-cs2_r*((-h_faceDD12**2/cf_face**4 + (h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**4)/(-h_faceDD01**2*(h_faceDD22 + 1)/cf_face**6
         * + 2*h_faceDD01*h_faceDD02*h_faceDD12/cf_face**6 - h_faceDD02**2*(h_faceDD11 + 1)/cf_face**6 - h_faceDD12**2*(h_faceDD00 + 1)/cf_face**6 +
         * (h_faceDD00 + 1)*(h_faceDD11 + 1)*(h_faceDD22 + 1)/cf_face**6) - vet_faceU0**2/alpha_face**2) + rescaledvrU0**2*u4rUt**2*(1 -
         * cs2_r))*(4*u4rUt**2*(1 - cs2_r) + 4*cs2_r/alpha_face**2) + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)**2)/2)/(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(4*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) + (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(4*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2 + (-2*rescaledvrU0*u4rUt**2*(1 - cs2_r) +
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) - (2*rescaledvrU0*u4rUt**2*(1 - cs2_r) -
         * 2*cs2_r*vet_faceU0/alpha_face**2)/(8*(u4rUt**2*(1 - cs2_r) + cs2_r/alpha_face**2)) + (-2*rescaledvlU0*u4lUt**2*(1 - cs2_l) +
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)) - (2*rescaledvlU0*u4lUt**2*(1 - cs2_l) -
         * 2*cs2_l*vet_faceU0/alpha_face**2)/(8*(u4lUt**2*(1 - cs2_l) + cs2_l/alpha_face**2)))/2)]"
         */
        {
          const REAL FDPart3tmp0 = ((alpha_face) * (alpha_face));
          const REAL FDPart3tmp2 = ((u4lUt) * (u4lUt));
          const REAL FDPart3tmp3 = 1 - cs2_l;
          const REAL FDPart3tmp13 = (1.0 / ((cf_face) * (cf_face) * (cf_face) * (cf_face)));
          const REAL FDPart3tmp14 = ((h_faceDD12) * (h_faceDD12));
          const REAL FDPart3tmp15 = h_faceDD11 + 1;
          const REAL FDPart3tmp16 = h_faceDD22 + 1;
          const REAL FDPart3tmp17 = pow(cf_face, -6);
          const REAL FDPart3tmp20 = h_faceDD00 + 1;
          const REAL FDPart3tmp25 = ((u4rUt) * (u4rUt));
          const REAL FDPart3tmp26 = 1 - cs2_r;
          const REAL FDPart3tmp60 = fabs(cf_face);
          const REAL FDPart3tmp1 = (1.0 / (FDPart3tmp0));
          const REAL FDPart3tmp4 = FDPart3tmp2 * FDPart3tmp3;
          const REAL FDPart3tmp18 = FDPart3tmp16 * ((h_faceDD01) * (h_faceDD01));
          const REAL FDPart3tmp19 = FDPart3tmp15 * ((h_faceDD02) * (h_faceDD02));
          const REAL FDPart3tmp27 = FDPart3tmp25 * FDPart3tmp26;
          const REAL FDPart3tmp61 = FDPart3tmp13 * FDPart3tmp60;
          const REAL FDPart3tmp70 = FDPart3tmp0 * FDPart3tmp2 * h_l * rhob_l;
          const REAL FDPart3tmp73 = FDPart3tmp0 * FDPart3tmp25 * h_r * rhob_r;
          const REAL FDPart3tmp76 = FDPart3tmp0 * FDPart3tmp25 * h_r * rhob_r - P_r;
          const REAL FDPart3tmp77 = FDPart3tmp0 * FDPart3tmp2 * h_l * rhob_l - P_l;
          const REAL FDPart3tmp83 = FDPart3tmp0 * ((cf_face) * (cf_face));
          const REAL FDPart3tmp5 = -2 * FDPart3tmp1 * cs2_l * vet_faceU0 + 2 * FDPart3tmp4 * rescaledvlU0;
          const REAL FDPart3tmp22 =
              -FDPart3tmp1 * ((vet_faceU0) * (vet_faceU0)) +
              (-FDPart3tmp13 * FDPart3tmp14 + FDPart3tmp13 * FDPart3tmp15 * FDPart3tmp16) /
                  (-FDPart3tmp14 * FDPart3tmp17 * FDPart3tmp20 + FDPart3tmp15 * FDPart3tmp16 * FDPart3tmp17 * FDPart3tmp20 -
                   FDPart3tmp17 * FDPart3tmp18 - FDPart3tmp17 * FDPart3tmp19 + 2 * FDPart3tmp17 * h_faceDD01 * h_faceDD02 * h_faceDD12);
          const REAL FDPart3tmp28 = -2 * FDPart3tmp1 * cs2_r * vet_faceU0 + 2 * FDPart3tmp27 * rescaledvrU0;
          const REAL FDPart3tmp62 = FDPart3tmp61 * alpha_face;
          const REAL FDPart3tmp79 = FDPart3tmp1 / ((cf_face) * (cf_face));
          const REAL FDPart3tmp81 = -FDPart3tmp14 * FDPart3tmp20 + FDPart3tmp15 * FDPart3tmp16 * FDPart3tmp20 - FDPart3tmp18 - FDPart3tmp19 +
                                    2 * h_faceDD01 * h_faceDD02 * h_faceDD12;
          const REAL FDPart3tmp104 = FDPart3tmp70 * rescaledvlU1 + P_l * vet_faceU1;
          const REAL FDPart3tmp105 = FDPart3tmp70 * rescaledvlU2 + P_l * vet_faceU2;
          const REAL FDPart3tmp106 = FDPart3tmp73 * rescaledvrU1 + P_r * vet_faceU1;
          const REAL FDPart3tmp107 = FDPart3tmp73 * rescaledvrU2 + P_r * vet_faceU2;
          const REAL FDPart3tmp8 = (1.0 / (FDPart3tmp1 * cs2_l + FDPart3tmp4));
          const REAL FDPart3tmp23 =
              (4 * FDPart3tmp1 * cs2_l + 4 * FDPart3tmp4) * (FDPart3tmp2 * FDPart3tmp3 * ((rescaledvlU0) * (rescaledvlU0)) - FDPart3tmp22 * cs2_l);
          const REAL FDPart3tmp31 = (1.0 / (FDPart3tmp1 * cs2_r + FDPart3tmp27));
          const REAL FDPart3tmp35 = (4 * FDPart3tmp1 * cs2_r + 4 * FDPart3tmp27) *
                                    (-FDPart3tmp22 * cs2_r + FDPart3tmp25 * FDPart3tmp26 * ((rescaledvrU0) * (rescaledvrU0)));
          const REAL FDPart3tmp63 = FDPart3tmp62 * rhob_l * u4lUt;
          const REAL FDPart3tmp66 = FDPart3tmp62 * rhob_r * u4rUt;
          const REAL FDPart3tmp72 = FDPart3tmp70 * rescaledvlU0 + P_l * vet_faceU0;
          const REAL FDPart3tmp75 = FDPart3tmp73 * rescaledvrU0 + P_r * vet_faceU0;
          const REAL FDPart3tmp80 = FDPart3tmp79 * (FDPart3tmp20 * vet_faceU0 + h_faceDD01 * vet_faceU1 + h_faceDD02 * vet_faceU2);
          const REAL FDPart3tmp82 = FDPart3tmp70 * FDPart3tmp81 * rescaledvlU0;
          const REAL FDPart3tmp87 = (1.0 / (FDPart3tmp81));
          const REAL FDPart3tmp88 = FDPart3tmp79 * h_faceDD01;
          const REAL FDPart3tmp92 = FDPart3tmp79 * h_faceDD02;
          const REAL FDPart3tmp94 = -FDPart3tmp81 * ((vet_faceU0) * (vet_faceU0)) + FDPart3tmp83 * (-FDPart3tmp14 + FDPart3tmp15 * FDPart3tmp16);
          const REAL FDPart3tmp96 = FDPart3tmp20 * FDPart3tmp79;
          const REAL FDPart3tmp99 = FDPart3tmp73 * FDPart3tmp81 * rescaledvrU0;
          const REAL FDPart3tmp108 = FDPart3tmp79 * (FDPart3tmp15 * vet_faceU1 + h_faceDD01 * vet_faceU0 + h_faceDD12 * vet_faceU2);
          const REAL FDPart3tmp109 = FDPart3tmp79 * h_faceDD12;
          const REAL FDPart3tmp111 = FDPart3tmp15 * FDPart3tmp79;
          const REAL FDPart3tmp113 = FDPart3tmp79 * (FDPart3tmp16 * vet_faceU2 + h_faceDD02 * vet_faceU0 + h_faceDD12 * vet_faceU1);
          const REAL FDPart3tmp114 = FDPart3tmp16 * FDPart3tmp79;
          const REAL FDPart3tmp68 = FDPart3tmp63 - FDPart3tmp66;
          const REAL FDPart3tmp85 = -FDPart3tmp81 * vet_faceU0 * vet_faceU1 + FDPart3tmp83 * (-FDPart3tmp16 * h_faceDD01 + h_faceDD02 * h_faceDD12);
          const REAL FDPart3tmp89 = FDPart3tmp87 * FDPart3tmp88;
          const REAL FDPart3tmp90 = -FDPart3tmp81 * vet_faceU0 * vet_faceU2 + FDPart3tmp83 * (-FDPart3tmp15 * h_faceDD02 + h_faceDD01 * h_faceDD12);
          const REAL FDPart3tmp93 = FDPart3tmp87 * FDPart3tmp92;
          const REAL FDPart3tmp95 = FDPart3tmp70 * FDPart3tmp81 * ((rescaledvlU0) * (rescaledvlU0)) + FDPart3tmp94 * P_l;
          const REAL FDPart3tmp102 = FDPart3tmp73 * FDPart3tmp81 * ((rescaledvrU0) * (rescaledvrU0)) + FDPart3tmp94 * P_r;
          const REAL FDPart3tmp110 = FDPart3tmp109 * FDPart3tmp87;
          const REAL FDPart3tmp24 = fabs(FDPart3tmp8 * sqrt(-1.0 / 2.0 * FDPart3tmp23 + (1.0 / 2.0) * ((FDPart3tmp5) * (FDPart3tmp5)) +
                                                            (1.0 / 2.0) * fabs(-FDPart3tmp23 + ((FDPart3tmp5) * (FDPart3tmp5)))));
          const REAL FDPart3tmp36 = fabs(FDPart3tmp31 * sqrt((1.0 / 2.0) * ((FDPart3tmp28) * (FDPart3tmp28)) - 1.0 / 2.0 * FDPart3tmp35 +
                                                             (1.0 / 2.0) * fabs(((FDPart3tmp28) * (FDPart3tmp28)) - FDPart3tmp35)));
          const REAL FDPart3tmp39 = -1.0 / 4.0 * FDPart3tmp5 * FDPart3tmp8;
          const REAL FDPart3tmp41 = (1.0 / 4.0) * FDPart3tmp28 * FDPart3tmp31;
          const REAL FDPart3tmp45 = (1.0 / 8.0) * FDPart3tmp5 * FDPart3tmp8;
          const REAL FDPart3tmp47 = (1.0 / 8.0) * FDPart3tmp28 * FDPart3tmp31;
          const REAL FDPart3tmp53 = -1.0 / 16.0 * FDPart3tmp5 * FDPart3tmp8;
          const REAL FDPart3tmp55 = -1.0 / 16.0 * FDPart3tmp28 * FDPart3tmp31;
          const REAL FDPart3tmp56 = (1.0 / 16.0) * FDPart3tmp5 * FDPart3tmp8;
          const REAL FDPart3tmp57 = (1.0 / 16.0) * FDPart3tmp28 * FDPart3tmp31;
          const REAL FDPart3tmp86 = FDPart3tmp82 * rescaledvlU1 + FDPart3tmp85 * P_l;
          const REAL FDPart3tmp91 = FDPart3tmp82 * rescaledvlU2 + FDPart3tmp90 * P_l;
          const REAL FDPart3tmp100 = FDPart3tmp85 * P_r + FDPart3tmp99 * rescaledvrU1;
          const REAL FDPart3tmp101 = FDPart3tmp90 * P_r + FDPart3tmp99 * rescaledvrU2;
          const REAL FDPart3tmp37 = (1.0 / 4.0) * FDPart3tmp24 + (1.0 / 4.0) * FDPart3tmp36;
          const REAL FDPart3tmp42 = (1.0 / 2.0) * FDPart3tmp24 - 1.0 / 2.0 * FDPart3tmp36;
          const REAL FDPart3tmp58 = (1.0 / 8.0) * FDPart3tmp24 + (1.0 / 8.0) * FDPart3tmp36;
          const REAL FDPart3tmp43 =
              fabs(-1.0 / 4.0 * FDPart3tmp28 * FDPart3tmp31 - FDPart3tmp39 - FDPart3tmp41 - FDPart3tmp42 + (1.0 / 4.0) * FDPart3tmp5 * FDPart3tmp8);
          const REAL FDPart3tmp49 =
              fabs(-1.0 / 4.0 * FDPart3tmp28 * FDPart3tmp31 - FDPart3tmp39 - FDPart3tmp41 + FDPart3tmp42 + (1.0 / 4.0) * FDPart3tmp5 * FDPart3tmp8);
          const REAL FDPart3tmp48 =
              (1.0 / 4.0) * FDPart3tmp43 + (1.0 / 2.0) * fabs(-1.0 / 8.0 * FDPart3tmp28 * FDPart3tmp31 + FDPart3tmp37 + (1.0 / 2.0) * FDPart3tmp43 -
                                                              FDPart3tmp45 - FDPart3tmp47 - 1.0 / 8.0 * FDPart3tmp5 * FDPart3tmp8);
          const REAL FDPart3tmp50 =
              (1.0 / 4.0) * FDPart3tmp49 + (1.0 / 2.0) * fabs(-1.0 / 8.0 * FDPart3tmp28 * FDPart3tmp31 - FDPart3tmp37 - FDPart3tmp45 - FDPart3tmp47 -
                                                              1.0 / 2.0 * FDPart3tmp49 - 1.0 / 8.0 * FDPart3tmp5 * FDPart3tmp8);
          const REAL FDPart3tmp51 = (1.0 / (FDPart3tmp37 + FDPart3tmp48 + FDPart3tmp50));
          const REAL FDPart3tmp59 = FDPart3tmp50 - FDPart3tmp53 - FDPart3tmp55 + FDPart3tmp56 + FDPart3tmp57 + FDPart3tmp58;
          const REAL FDPart3tmp65 = FDPart3tmp48 + FDPart3tmp53 + FDPart3tmp55 - FDPart3tmp56 - FDPart3tmp57 + FDPart3tmp58;
          const REAL FDPart3tmp69 = FDPart3tmp59 * FDPart3tmp65;
          const REAL FDPart3tmp98 = FDPart3tmp59 * FDPart3tmp62;
          const REAL FDPart3tmp103 = FDPart3tmp62 * FDPart3tmp65;
          auxevol_gfs[IDX4(RHO_STAR_HLL_FLUXGF, i0, i1, i2)] =
              FDPart3tmp51 * (FDPart3tmp59 * FDPart3tmp63 * rescaledvlU0 + FDPart3tmp65 * FDPart3tmp66 * rescaledvrU0 + FDPart3tmp68 * FDPart3tmp69);
          auxevol_gfs[IDX4(TAU_TILDE_HLL_FLUXGF, i0, i1, i2)] =
              FDPart3tmp51 * (FDPart3tmp59 * (FDPart3tmp13 * FDPart3tmp60 * FDPart3tmp72 - FDPart3tmp63 * rescaledvlU0) +
                              FDPart3tmp65 * (FDPart3tmp13 * FDPart3tmp60 * FDPart3tmp75 - FDPart3tmp66 * rescaledvrU0) -
                              FDPart3tmp69 * (FDPart3tmp61 * FDPart3tmp76 - FDPart3tmp61 * FDPart3tmp77 + FDPart3tmp68));
          auxevol_gfs[IDX4(RESCALEDSTILDE_FLUX_HLLD0GF, i0, i1, i2)] =
              FDPart3tmp51 *
              (FDPart3tmp103 * (FDPart3tmp100 * FDPart3tmp89 + FDPart3tmp101 * FDPart3tmp93 + FDPart3tmp102 * FDPart3tmp87 * FDPart3tmp96 +
                                FDPart3tmp75 * FDPart3tmp80) -
               FDPart3tmp69 *
                   (FDPart3tmp13 * FDPart3tmp60 * alpha_face *
                        (FDPart3tmp106 * FDPart3tmp88 + FDPart3tmp107 * FDPart3tmp92 + FDPart3tmp75 * FDPart3tmp96 + FDPart3tmp76 * FDPart3tmp80) -
                    FDPart3tmp62 *
                        (FDPart3tmp104 * FDPart3tmp88 + FDPart3tmp105 * FDPart3tmp92 + FDPart3tmp72 * FDPart3tmp96 + FDPart3tmp77 * FDPart3tmp80)) +
               FDPart3tmp98 * (FDPart3tmp72 * FDPart3tmp80 + FDPart3tmp86 * FDPart3tmp89 + FDPart3tmp87 * FDPart3tmp95 * FDPart3tmp96 +
                               FDPart3tmp91 * FDPart3tmp93));
          auxevol_gfs[IDX4(RESCALEDSTILDE_FLUX_HLLD1GF, i0, i1, i2)] =
              FDPart3tmp51 * (FDPart3tmp103 * (FDPart3tmp100 * FDPart3tmp111 * FDPart3tmp87 + FDPart3tmp101 * FDPart3tmp110 +
                                               FDPart3tmp102 * FDPart3tmp89 + FDPart3tmp108 * FDPart3tmp75) -
                              FDPart3tmp69 * (FDPart3tmp13 * FDPart3tmp60 * alpha_face *
                                                  (FDPart3tmp106 * FDPart3tmp111 + FDPart3tmp107 * FDPart3tmp109 + FDPart3tmp108 * FDPart3tmp76 +
                                                   FDPart3tmp75 * FDPart3tmp88) -
                                              FDPart3tmp62 * (FDPart3tmp104 * FDPart3tmp111 + FDPart3tmp105 * FDPart3tmp109 +
                                                              FDPart3tmp108 * FDPart3tmp77 + FDPart3tmp72 * FDPart3tmp88)) +
                              FDPart3tmp98 * (FDPart3tmp108 * FDPart3tmp72 + FDPart3tmp110 * FDPart3tmp91 +
                                              FDPart3tmp111 * FDPart3tmp86 * FDPart3tmp87 + FDPart3tmp89 * FDPart3tmp95));
          auxevol_gfs[IDX4(RESCALEDSTILDE_FLUX_HLLD2GF, i0, i1, i2)] =
              FDPart3tmp51 * (FDPart3tmp103 * (FDPart3tmp100 * FDPart3tmp110 + FDPart3tmp101 * FDPart3tmp114 * FDPart3tmp87 +
                                               FDPart3tmp102 * FDPart3tmp93 + FDPart3tmp113 * FDPart3tmp75) -
                              FDPart3tmp69 * (FDPart3tmp13 * FDPart3tmp60 * alpha_face *
                                                  (FDPart3tmp106 * FDPart3tmp109 + FDPart3tmp107 * FDPart3tmp114 + FDPart3tmp113 * FDPart3tmp76 +
                                                   FDPart3tmp75 * FDPart3tmp92) -
                                              FDPart3tmp62 * (FDPart3tmp104 * FDPart3tmp109 + FDPart3tmp105 * FDPart3tmp114 +
                                                              FDPart3tmp113 * FDPart3tmp77 + FDPart3tmp72 * FDPart3tmp92)) +
                              FDPart3tmp98 * (FDPart3tmp110 * FDPart3tmp86 + FDPart3tmp113 * FDPart3tmp72 +
                                              FDPart3tmp114 * FDPart3tmp87 * FDPart3tmp91 + FDPart3tmp93 * FDPart3tmp95));
        }

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS + 1; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS + 1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS + 1; i2++)
} // END FUNCTION calculate_HLL_fluxes_direction_0
