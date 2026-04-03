#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Add source and connection terms to rho_star, tau_tilde, rescaledStildeD, and optional composition variables
 */
void calculate_all_source_terms__rfm__Spherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                REAL *restrict xx[3], const ghl_eos_parameters *restrict eos, REAL *restrict auxevol_gfs,
                                                const REAL *restrict in_gfs, REAL *restrict rhs_gfs) {
#include "set_CodeParameters.h"
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    MAYBE_UNUSED const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      MAYBE_UNUSED const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        MAYBE_UNUSED const REAL xx0 = xx[0][i0];

        ghl_primitive_quantities prims;
        ghl_initialize_primitives(auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)], auxevol_gfs[IDX4(PGF, i0, i1, i2)], NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                                  NAN, NAN, &prims);

        double h, cs2;
        ghl_compute_h_and_cs2(eos, &prims, &h, &cs2);
        /*
         * NRPy-Generated GF Access/FD Code, Step 1 of 2:
         * Read gridfunction(s) from main memory and compute FD stencils as needed.
         */
        const REAL aDD00 = in_gfs[IDX4(ADD00GF, i0, i1, i2)];
        const REAL aDD01 = in_gfs[IDX4(ADD01GF, i0, i1, i2)];
        const REAL aDD02 = in_gfs[IDX4(ADD02GF, i0, i1, i2)];
        const REAL aDD11 = in_gfs[IDX4(ADD11GF, i0, i1, i2)];
        const REAL aDD12 = in_gfs[IDX4(ADD12GF, i0, i1, i2)];
        const REAL aDD22 = in_gfs[IDX4(ADD22GF, i0, i1, i2)];
        const REAL alpha_i2m2 = in_gfs[IDX4(ALPHAGF, i0, i1, i2 - 2)];
        const REAL alpha_i2m1 = in_gfs[IDX4(ALPHAGF, i0, i1, i2 - 1)];
        const REAL alpha_i1m2 = in_gfs[IDX4(ALPHAGF, i0, i1 - 2, i2)];
        const REAL alpha_i1m1 = in_gfs[IDX4(ALPHAGF, i0, i1 - 1, i2)];
        const REAL alpha_i0m2 = in_gfs[IDX4(ALPHAGF, i0 - 2, i1, i2)];
        const REAL alpha_i0m1 = in_gfs[IDX4(ALPHAGF, i0 - 1, i1, i2)];
        const REAL alpha = in_gfs[IDX4(ALPHAGF, i0, i1, i2)];
        const REAL alpha_i0p1 = in_gfs[IDX4(ALPHAGF, i0 + 1, i1, i2)];
        const REAL alpha_i0p2 = in_gfs[IDX4(ALPHAGF, i0 + 2, i1, i2)];
        const REAL alpha_i1p1 = in_gfs[IDX4(ALPHAGF, i0, i1 + 1, i2)];
        const REAL alpha_i1p2 = in_gfs[IDX4(ALPHAGF, i0, i1 + 2, i2)];
        const REAL alpha_i2p1 = in_gfs[IDX4(ALPHAGF, i0, i1, i2 + 1)];
        const REAL alpha_i2p2 = in_gfs[IDX4(ALPHAGF, i0, i1, i2 + 2)];
        const REAL cf_i2m2 = in_gfs[IDX4(CFGF, i0, i1, i2 - 2)];
        const REAL cf_i2m1 = in_gfs[IDX4(CFGF, i0, i1, i2 - 1)];
        const REAL cf_i1m2 = in_gfs[IDX4(CFGF, i0, i1 - 2, i2)];
        const REAL cf_i1m1 = in_gfs[IDX4(CFGF, i0, i1 - 1, i2)];
        const REAL cf_i0m2 = in_gfs[IDX4(CFGF, i0 - 2, i1, i2)];
        const REAL cf_i0m1 = in_gfs[IDX4(CFGF, i0 - 1, i1, i2)];
        const REAL cf = in_gfs[IDX4(CFGF, i0, i1, i2)];
        const REAL cf_i0p1 = in_gfs[IDX4(CFGF, i0 + 1, i1, i2)];
        const REAL cf_i0p2 = in_gfs[IDX4(CFGF, i0 + 2, i1, i2)];
        const REAL cf_i1p1 = in_gfs[IDX4(CFGF, i0, i1 + 1, i2)];
        const REAL cf_i1p2 = in_gfs[IDX4(CFGF, i0, i1 + 2, i2)];
        const REAL cf_i2p1 = in_gfs[IDX4(CFGF, i0, i1, i2 + 1)];
        const REAL cf_i2p2 = in_gfs[IDX4(CFGF, i0, i1, i2 + 2)];
        const REAL hDD00_i2m2 = in_gfs[IDX4(HDD00GF, i0, i1, i2 - 2)];
        const REAL hDD00_i2m1 = in_gfs[IDX4(HDD00GF, i0, i1, i2 - 1)];
        const REAL hDD00_i1m2 = in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2)];
        const REAL hDD00_i1m1 = in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2)];
        const REAL hDD00_i0m2 = in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2)];
        const REAL hDD00_i0m1 = in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2)];
        const REAL hDD00 = in_gfs[IDX4(HDD00GF, i0, i1, i2)];
        const REAL hDD00_i0p1 = in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2)];
        const REAL hDD00_i0p2 = in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2)];
        const REAL hDD00_i1p1 = in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2)];
        const REAL hDD00_i1p2 = in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2)];
        const REAL hDD00_i2p1 = in_gfs[IDX4(HDD00GF, i0, i1, i2 + 1)];
        const REAL hDD00_i2p2 = in_gfs[IDX4(HDD00GF, i0, i1, i2 + 2)];
        const REAL hDD01_i2m2 = in_gfs[IDX4(HDD01GF, i0, i1, i2 - 2)];
        const REAL hDD01_i2m1 = in_gfs[IDX4(HDD01GF, i0, i1, i2 - 1)];
        const REAL hDD01_i1m2 = in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2)];
        const REAL hDD01_i1m1 = in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2)];
        const REAL hDD01_i0m2 = in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2)];
        const REAL hDD01_i0m1 = in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2)];
        const REAL hDD01 = in_gfs[IDX4(HDD01GF, i0, i1, i2)];
        const REAL hDD01_i0p1 = in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2)];
        const REAL hDD01_i0p2 = in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2)];
        const REAL hDD01_i1p1 = in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2)];
        const REAL hDD01_i1p2 = in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2)];
        const REAL hDD01_i2p1 = in_gfs[IDX4(HDD01GF, i0, i1, i2 + 1)];
        const REAL hDD01_i2p2 = in_gfs[IDX4(HDD01GF, i0, i1, i2 + 2)];
        const REAL hDD02_i2m2 = in_gfs[IDX4(HDD02GF, i0, i1, i2 - 2)];
        const REAL hDD02_i2m1 = in_gfs[IDX4(HDD02GF, i0, i1, i2 - 1)];
        const REAL hDD02_i1m2 = in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2)];
        const REAL hDD02_i1m1 = in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2)];
        const REAL hDD02_i0m2 = in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2)];
        const REAL hDD02_i0m1 = in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2)];
        const REAL hDD02 = in_gfs[IDX4(HDD02GF, i0, i1, i2)];
        const REAL hDD02_i0p1 = in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2)];
        const REAL hDD02_i0p2 = in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2)];
        const REAL hDD02_i1p1 = in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2)];
        const REAL hDD02_i1p2 = in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2)];
        const REAL hDD02_i2p1 = in_gfs[IDX4(HDD02GF, i0, i1, i2 + 1)];
        const REAL hDD02_i2p2 = in_gfs[IDX4(HDD02GF, i0, i1, i2 + 2)];
        const REAL hDD11_i2m2 = in_gfs[IDX4(HDD11GF, i0, i1, i2 - 2)];
        const REAL hDD11_i2m1 = in_gfs[IDX4(HDD11GF, i0, i1, i2 - 1)];
        const REAL hDD11_i1m2 = in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2)];
        const REAL hDD11_i1m1 = in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2)];
        const REAL hDD11_i0m2 = in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2)];
        const REAL hDD11_i0m1 = in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2)];
        const REAL hDD11 = in_gfs[IDX4(HDD11GF, i0, i1, i2)];
        const REAL hDD11_i0p1 = in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2)];
        const REAL hDD11_i0p2 = in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2)];
        const REAL hDD11_i1p1 = in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2)];
        const REAL hDD11_i1p2 = in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2)];
        const REAL hDD11_i2p1 = in_gfs[IDX4(HDD11GF, i0, i1, i2 + 1)];
        const REAL hDD11_i2p2 = in_gfs[IDX4(HDD11GF, i0, i1, i2 + 2)];
        const REAL hDD12_i2m2 = in_gfs[IDX4(HDD12GF, i0, i1, i2 - 2)];
        const REAL hDD12_i2m1 = in_gfs[IDX4(HDD12GF, i0, i1, i2 - 1)];
        const REAL hDD12_i1m2 = in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2)];
        const REAL hDD12_i1m1 = in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2)];
        const REAL hDD12_i0m2 = in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2)];
        const REAL hDD12_i0m1 = in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2)];
        const REAL hDD12 = in_gfs[IDX4(HDD12GF, i0, i1, i2)];
        const REAL hDD12_i0p1 = in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2)];
        const REAL hDD12_i0p2 = in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2)];
        const REAL hDD12_i1p1 = in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2)];
        const REAL hDD12_i1p2 = in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2)];
        const REAL hDD12_i2p1 = in_gfs[IDX4(HDD12GF, i0, i1, i2 + 1)];
        const REAL hDD12_i2p2 = in_gfs[IDX4(HDD12GF, i0, i1, i2 + 2)];
        const REAL hDD22_i2m2 = in_gfs[IDX4(HDD22GF, i0, i1, i2 - 2)];
        const REAL hDD22_i2m1 = in_gfs[IDX4(HDD22GF, i0, i1, i2 - 1)];
        const REAL hDD22_i1m2 = in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2)];
        const REAL hDD22_i1m1 = in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2)];
        const REAL hDD22_i0m2 = in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2)];
        const REAL hDD22_i0m1 = in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2)];
        const REAL hDD22 = in_gfs[IDX4(HDD22GF, i0, i1, i2)];
        const REAL hDD22_i0p1 = in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2)];
        const REAL hDD22_i0p2 = in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2)];
        const REAL hDD22_i1p1 = in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2)];
        const REAL hDD22_i1p2 = in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2)];
        const REAL hDD22_i2p1 = in_gfs[IDX4(HDD22GF, i0, i1, i2 + 1)];
        const REAL hDD22_i2p2 = in_gfs[IDX4(HDD22GF, i0, i1, i2 + 2)];
        const REAL trK = in_gfs[IDX4(TRKGF, i0, i1, i2)];
        const REAL vetU0_i2m2 = in_gfs[IDX4(VETU0GF, i0, i1, i2 - 2)];
        const REAL vetU0_i2m1 = in_gfs[IDX4(VETU0GF, i0, i1, i2 - 1)];
        const REAL vetU0_i1m2 = in_gfs[IDX4(VETU0GF, i0, i1 - 2, i2)];
        const REAL vetU0_i1m1 = in_gfs[IDX4(VETU0GF, i0, i1 - 1, i2)];
        const REAL vetU0_i0m2 = in_gfs[IDX4(VETU0GF, i0 - 2, i1, i2)];
        const REAL vetU0_i0m1 = in_gfs[IDX4(VETU0GF, i0 - 1, i1, i2)];
        const REAL vetU0 = in_gfs[IDX4(VETU0GF, i0, i1, i2)];
        const REAL vetU0_i0p1 = in_gfs[IDX4(VETU0GF, i0 + 1, i1, i2)];
        const REAL vetU0_i0p2 = in_gfs[IDX4(VETU0GF, i0 + 2, i1, i2)];
        const REAL vetU0_i1p1 = in_gfs[IDX4(VETU0GF, i0, i1 + 1, i2)];
        const REAL vetU0_i1p2 = in_gfs[IDX4(VETU0GF, i0, i1 + 2, i2)];
        const REAL vetU0_i2p1 = in_gfs[IDX4(VETU0GF, i0, i1, i2 + 1)];
        const REAL vetU0_i2p2 = in_gfs[IDX4(VETU0GF, i0, i1, i2 + 2)];
        const REAL vetU1_i2m2 = in_gfs[IDX4(VETU1GF, i0, i1, i2 - 2)];
        const REAL vetU1_i2m1 = in_gfs[IDX4(VETU1GF, i0, i1, i2 - 1)];
        const REAL vetU1_i1m2 = in_gfs[IDX4(VETU1GF, i0, i1 - 2, i2)];
        const REAL vetU1_i1m1 = in_gfs[IDX4(VETU1GF, i0, i1 - 1, i2)];
        const REAL vetU1_i0m2 = in_gfs[IDX4(VETU1GF, i0 - 2, i1, i2)];
        const REAL vetU1_i0m1 = in_gfs[IDX4(VETU1GF, i0 - 1, i1, i2)];
        const REAL vetU1 = in_gfs[IDX4(VETU1GF, i0, i1, i2)];
        const REAL vetU1_i0p1 = in_gfs[IDX4(VETU1GF, i0 + 1, i1, i2)];
        const REAL vetU1_i0p2 = in_gfs[IDX4(VETU1GF, i0 + 2, i1, i2)];
        const REAL vetU1_i1p1 = in_gfs[IDX4(VETU1GF, i0, i1 + 1, i2)];
        const REAL vetU1_i1p2 = in_gfs[IDX4(VETU1GF, i0, i1 + 2, i2)];
        const REAL vetU1_i2p1 = in_gfs[IDX4(VETU1GF, i0, i1, i2 + 1)];
        const REAL vetU1_i2p2 = in_gfs[IDX4(VETU1GF, i0, i1, i2 + 2)];
        const REAL vetU2_i2m2 = in_gfs[IDX4(VETU2GF, i0, i1, i2 - 2)];
        const REAL vetU2_i2m1 = in_gfs[IDX4(VETU2GF, i0, i1, i2 - 1)];
        const REAL vetU2_i1m2 = in_gfs[IDX4(VETU2GF, i0, i1 - 2, i2)];
        const REAL vetU2_i1m1 = in_gfs[IDX4(VETU2GF, i0, i1 - 1, i2)];
        const REAL vetU2_i0m2 = in_gfs[IDX4(VETU2GF, i0 - 2, i1, i2)];
        const REAL vetU2_i0m1 = in_gfs[IDX4(VETU2GF, i0 - 1, i1, i2)];
        const REAL vetU2 = in_gfs[IDX4(VETU2GF, i0, i1, i2)];
        const REAL vetU2_i0p1 = in_gfs[IDX4(VETU2GF, i0 + 1, i1, i2)];
        const REAL vetU2_i0p2 = in_gfs[IDX4(VETU2GF, i0 + 2, i1, i2)];
        const REAL vetU2_i1p1 = in_gfs[IDX4(VETU2GF, i0, i1 + 1, i2)];
        const REAL vetU2_i1p2 = in_gfs[IDX4(VETU2GF, i0, i1 + 2, i2)];
        const REAL vetU2_i2p1 = in_gfs[IDX4(VETU2GF, i0, i1, i2 + 1)];
        const REAL vetU2_i2p2 = in_gfs[IDX4(VETU2GF, i0, i1, i2 + 2)];
        static const REAL FDPart1_Rational_2_3 = 2.0 / 3.0;
        static const REAL FDPart1_Rational_1_12 = 1.0 / 12.0;
        const REAL alpha_dD0 = invdxx0 * (FDPart1_Rational_1_12 * (alpha_i0m2 - alpha_i0p2) + FDPart1_Rational_2_3 * (-alpha_i0m1 + alpha_i0p1));
        const REAL alpha_dD1 = invdxx1 * (FDPart1_Rational_1_12 * (alpha_i1m2 - alpha_i1p2) + FDPart1_Rational_2_3 * (-alpha_i1m1 + alpha_i1p1));
        const REAL alpha_dD2 = invdxx2 * (FDPart1_Rational_1_12 * (alpha_i2m2 - alpha_i2p2) + FDPart1_Rational_2_3 * (-alpha_i2m1 + alpha_i2p1));
        const REAL cf_dD0 = invdxx0 * (FDPart1_Rational_1_12 * (cf_i0m2 - cf_i0p2) + FDPart1_Rational_2_3 * (-cf_i0m1 + cf_i0p1));
        const REAL cf_dD1 = invdxx1 * (FDPart1_Rational_1_12 * (cf_i1m2 - cf_i1p2) + FDPart1_Rational_2_3 * (-cf_i1m1 + cf_i1p1));
        const REAL cf_dD2 = invdxx2 * (FDPart1_Rational_1_12 * (cf_i2m2 - cf_i2p2) + FDPart1_Rational_2_3 * (-cf_i2m1 + cf_i2p1));
        const REAL hDD_dD000 = invdxx0 * (FDPart1_Rational_1_12 * (hDD00_i0m2 - hDD00_i0p2) + FDPart1_Rational_2_3 * (-hDD00_i0m1 + hDD00_i0p1));
        const REAL hDD_dD001 = invdxx1 * (FDPart1_Rational_1_12 * (hDD00_i1m2 - hDD00_i1p2) + FDPart1_Rational_2_3 * (-hDD00_i1m1 + hDD00_i1p1));
        const REAL hDD_dD002 = invdxx2 * (FDPart1_Rational_1_12 * (hDD00_i2m2 - hDD00_i2p2) + FDPart1_Rational_2_3 * (-hDD00_i2m1 + hDD00_i2p1));
        const REAL hDD_dD010 = invdxx0 * (FDPart1_Rational_1_12 * (hDD01_i0m2 - hDD01_i0p2) + FDPart1_Rational_2_3 * (-hDD01_i0m1 + hDD01_i0p1));
        const REAL hDD_dD011 = invdxx1 * (FDPart1_Rational_1_12 * (hDD01_i1m2 - hDD01_i1p2) + FDPart1_Rational_2_3 * (-hDD01_i1m1 + hDD01_i1p1));
        const REAL hDD_dD012 = invdxx2 * (FDPart1_Rational_1_12 * (hDD01_i2m2 - hDD01_i2p2) + FDPart1_Rational_2_3 * (-hDD01_i2m1 + hDD01_i2p1));
        const REAL hDD_dD020 = invdxx0 * (FDPart1_Rational_1_12 * (hDD02_i0m2 - hDD02_i0p2) + FDPart1_Rational_2_3 * (-hDD02_i0m1 + hDD02_i0p1));
        const REAL hDD_dD021 = invdxx1 * (FDPart1_Rational_1_12 * (hDD02_i1m2 - hDD02_i1p2) + FDPart1_Rational_2_3 * (-hDD02_i1m1 + hDD02_i1p1));
        const REAL hDD_dD022 = invdxx2 * (FDPart1_Rational_1_12 * (hDD02_i2m2 - hDD02_i2p2) + FDPart1_Rational_2_3 * (-hDD02_i2m1 + hDD02_i2p1));
        const REAL hDD_dD110 = invdxx0 * (FDPart1_Rational_1_12 * (hDD11_i0m2 - hDD11_i0p2) + FDPart1_Rational_2_3 * (-hDD11_i0m1 + hDD11_i0p1));
        const REAL hDD_dD111 = invdxx1 * (FDPart1_Rational_1_12 * (hDD11_i1m2 - hDD11_i1p2) + FDPart1_Rational_2_3 * (-hDD11_i1m1 + hDD11_i1p1));
        const REAL hDD_dD112 = invdxx2 * (FDPart1_Rational_1_12 * (hDD11_i2m2 - hDD11_i2p2) + FDPart1_Rational_2_3 * (-hDD11_i2m1 + hDD11_i2p1));
        const REAL hDD_dD120 = invdxx0 * (FDPart1_Rational_1_12 * (hDD12_i0m2 - hDD12_i0p2) + FDPart1_Rational_2_3 * (-hDD12_i0m1 + hDD12_i0p1));
        const REAL hDD_dD121 = invdxx1 * (FDPart1_Rational_1_12 * (hDD12_i1m2 - hDD12_i1p2) + FDPart1_Rational_2_3 * (-hDD12_i1m1 + hDD12_i1p1));
        const REAL hDD_dD122 = invdxx2 * (FDPart1_Rational_1_12 * (hDD12_i2m2 - hDD12_i2p2) + FDPart1_Rational_2_3 * (-hDD12_i2m1 + hDD12_i2p1));
        const REAL hDD_dD220 = invdxx0 * (FDPart1_Rational_1_12 * (hDD22_i0m2 - hDD22_i0p2) + FDPart1_Rational_2_3 * (-hDD22_i0m1 + hDD22_i0p1));
        const REAL hDD_dD221 = invdxx1 * (FDPart1_Rational_1_12 * (hDD22_i1m2 - hDD22_i1p2) + FDPart1_Rational_2_3 * (-hDD22_i1m1 + hDD22_i1p1));
        const REAL hDD_dD222 = invdxx2 * (FDPart1_Rational_1_12 * (hDD22_i2m2 - hDD22_i2p2) + FDPart1_Rational_2_3 * (-hDD22_i2m1 + hDD22_i2p1));
        const REAL vetU_dD00 = invdxx0 * (FDPart1_Rational_1_12 * (vetU0_i0m2 - vetU0_i0p2) + FDPart1_Rational_2_3 * (-vetU0_i0m1 + vetU0_i0p1));
        const REAL vetU_dD01 = invdxx1 * (FDPart1_Rational_1_12 * (vetU0_i1m2 - vetU0_i1p2) + FDPart1_Rational_2_3 * (-vetU0_i1m1 + vetU0_i1p1));
        const REAL vetU_dD02 = invdxx2 * (FDPart1_Rational_1_12 * (vetU0_i2m2 - vetU0_i2p2) + FDPart1_Rational_2_3 * (-vetU0_i2m1 + vetU0_i2p1));
        const REAL vetU_dD10 = invdxx0 * (FDPart1_Rational_1_12 * (vetU1_i0m2 - vetU1_i0p2) + FDPart1_Rational_2_3 * (-vetU1_i0m1 + vetU1_i0p1));
        const REAL vetU_dD11 = invdxx1 * (FDPart1_Rational_1_12 * (vetU1_i1m2 - vetU1_i1p2) + FDPart1_Rational_2_3 * (-vetU1_i1m1 + vetU1_i1p1));
        const REAL vetU_dD12 = invdxx2 * (FDPart1_Rational_1_12 * (vetU1_i2m2 - vetU1_i2p2) + FDPart1_Rational_2_3 * (-vetU1_i2m1 + vetU1_i2p1));
        const REAL vetU_dD20 = invdxx0 * (FDPart1_Rational_1_12 * (vetU2_i0m2 - vetU2_i0p2) + FDPart1_Rational_2_3 * (-vetU2_i0m1 + vetU2_i0p1));
        const REAL vetU_dD21 = invdxx1 * (FDPart1_Rational_1_12 * (vetU2_i1m2 - vetU2_i1p2) + FDPart1_Rational_2_3 * (-vetU2_i1m1 + vetU2_i1p1));
        const REAL vetU_dD22 = invdxx2 * (FDPart1_Rational_1_12 * (vetU2_i2m2 - vetU2_i2p2) + FDPart1_Rational_2_3 * (-vetU2_i2m1 + vetU2_i2p1));

        /*
         * NRPy-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        const REAL FDPart3tmp0 = (1.0 / (xx0));
        const REAL FDPart3tmp2 = pow(fabs(cf), -3.0);
        const REAL FDPart3tmp5 = cos(xx1);
        const REAL FDPart3tmp6 = sin(xx1);
        const REAL FDPart3tmp10 = ((alpha) * (alpha));
        const REAL FDPart3tmp29 = hDD00 + 1;
        const REAL FDPart3tmp30 = (1.0 / 3.0) * trK;
        const REAL FDPart3tmp34 = hDD22 + 1;
        const REAL FDPart3tmp36 = hDD11 + 1;
        const REAL FDPart3tmp55 = ((xx0) * (xx0));
        const REAL FDPart3tmp77 = 2 * vetU2;
        const REAL FDPart3tmp98 = ((cf) * (cf) * (cf));
        const REAL FDPart3tmp105 = cf * xx0;
        const REAL FDPart3tmp110 = 1.0 * cf_dD0 * xx0;
        const REAL FDPart3tmp128 = 1.0 * cf;
        const REAL FDPart3tmp142 = 1.0 * cf_dD1;
        const REAL FDPart3tmp143 = 0.5 * cf;
        const REAL FDPart3tmp150 = 1.0 * cf_dD2;
        const REAL FDPart3tmp1 = 2 * FDPart3tmp0;
        const REAL FDPart3tmp3 = FDPart3tmp2 * alpha * rhob * u4Ut;
        const REAL FDPart3tmp7 = (1.0 / (FDPart3tmp6));
        const REAL FDPart3tmp12 = FDPart3tmp10 * h * rhob * ((u4Ut) * (u4Ut));
        const REAL FDPart3tmp17 = (1.0 / (FDPart3tmp10));
        const REAL FDPart3tmp19 = FDPart3tmp10 * h * rhob * ((u4Ut) * (u4Ut)) - P;
        const REAL FDPart3tmp28 = (1.0 / ((cf) * (cf)));
        const REAL FDPart3tmp38 = FDPart3tmp29 * FDPart3tmp34 * FDPart3tmp36 - FDPart3tmp29 * ((hDD12) * (hDD12)) -
                                  FDPart3tmp34 * ((hDD01) * (hDD01)) - FDPart3tmp36 * ((hDD02) * (hDD02)) + 2 * hDD01 * hDD02 * hDD12;
        const REAL FDPart3tmp40 = FDPart3tmp10 * ((cf) * (cf));
        const REAL FDPart3tmp57 = (1.0 / (FDPart3tmp55));
        const REAL FDPart3tmp70 = FDPart3tmp6 * hDD12;
        const REAL FDPart3tmp84 = FDPart3tmp36 * vetU1 + hDD01 * vetU0 + hDD12 * vetU2;
        const REAL FDPart3tmp89 = FDPart3tmp6 * cf;
        const REAL FDPart3tmp90 = FDPart3tmp34 * vetU2 + hDD02 * vetU0 + hDD12 * vetU1;
        const REAL FDPart3tmp94 = FDPart3tmp29 * vetU0 + hDD01 * vetU1 + hDD02 * vetU2;
        const REAL FDPart3tmp96 = FDPart3tmp5 * cf;
        const REAL FDPart3tmp109 = 0.5 * FDPart3tmp105;
        const REAL FDPart3tmp120 = 2 * cf_dD0 * xx0;
        const REAL FDPart3tmp135 = FDPart3tmp5 * vetU2;
        const REAL FDPart3tmp9 = FDPart3tmp0 * FDPart3tmp3 * rescaledvU1;
        const REAL FDPart3tmp20 = FDPart3tmp17 * FDPart3tmp19;
        const REAL FDPart3tmp23 = FDPart3tmp0 * FDPart3tmp7;
        const REAL FDPart3tmp24 = FDPart3tmp12 * rescaledvU2 + P * vetU2;
        const REAL FDPart3tmp26 = FDPart3tmp7 * vetU2;
        const REAL FDPart3tmp31 = FDPart3tmp28 * FDPart3tmp30;
        const REAL FDPart3tmp39 = FDPart3tmp12 * FDPart3tmp38;
        const REAL FDPart3tmp41 = FDPart3tmp38 * ((vetU0) * (vetU0));
        const REAL FDPart3tmp43 = (1.0 / (FDPart3tmp38));
        const REAL FDPart3tmp45 = FDPart3tmp28 * xx0;
        const REAL FDPart3tmp50 = FDPart3tmp38 * vetU0;
        const REAL FDPart3tmp56 = FDPart3tmp28 * FDPart3tmp55;
        const REAL FDPart3tmp62 = FDPart3tmp38 * ((vetU1) * (vetU1));
        const REAL FDPart3tmp73 = FDPart3tmp38 * vetU1;
        const REAL FDPart3tmp79 = FDPart3tmp55 * ((FDPart3tmp6) * (FDPart3tmp6));
        const REAL FDPart3tmp80 = (1.0 / ((FDPart3tmp6) * (FDPart3tmp6)));
        const REAL FDPart3tmp82 = FDPart3tmp38 * ((vetU2) * (vetU2));
        const REAL FDPart3tmp99 = FDPart3tmp19 * FDPart3tmp38 * FDPart3tmp98 * alpha;
        const REAL FDPart3tmp111 = FDPart3tmp109 * hDD_dD120 - FDPart3tmp110 * hDD12;
        const REAL FDPart3tmp121 = 0.5 * FDPart3tmp105 * hDD_dD010 - 0.5 * FDPart3tmp120 * hDD01;
        const REAL FDPart3tmp125 = 0.5 * FDPart3tmp105 * hDD_dD020 - 0.5 * FDPart3tmp120 * hDD02;
        const REAL FDPart3tmp138 = -1.0 * FDPart3tmp70 * cf_dD1 + 0.5 * FDPart3tmp89 * hDD02 + 0.5 * FDPart3tmp89 * hDD_dD121;
        const REAL FDPart3tmp140 = 0.5 * FDPart3tmp6 * cf * hDD_dD021 - 1.0 * FDPart3tmp6 * cf_dD1 * hDD02 - 0.5 * FDPart3tmp70 * cf;
        const REAL FDPart3tmp147 = 0.5 * FDPart3tmp6 * (FDPart3tmp29 * cf - FDPart3tmp36 * cf + cf * hDD_dD011 - 2 * cf_dD1 * hDD01);
        const REAL FDPart3tmp149 = 0.5 * cf * hDD_dD012 - 0.5 * cf * (FDPart3tmp5 * hDD02 + FDPart3tmp70) - 1.0 * cf_dD2 * hDD01;
        const REAL FDPart3tmp151 = FDPart3tmp128 * FDPart3tmp6 * hDD02;
        const REAL FDPart3tmp152 = FDPart3tmp128 * FDPart3tmp5 * hDD12;
        const REAL FDPart3tmp154 = -0.5 * FDPart3tmp34 * FDPart3tmp96 + 0.5 * FDPart3tmp36 * FDPart3tmp5 * cf + 0.5 * FDPart3tmp89 * hDD01 +
                                   0.5 * cf * hDD_dD122 - 1.0 * cf_dD2 * hDD12;
        const REAL FDPart3tmp155 = 0.5 * FDPart3tmp29 * FDPart3tmp6 * cf - 0.5 * FDPart3tmp34 * FDPart3tmp89 + 0.5 * FDPart3tmp96 * hDD01 +
                                   0.5 * cf * hDD_dD022 - 1.0 * cf_dD2 * hDD02;
        const REAL FDPart3tmp14 = FDPart3tmp12 * rescaledvU0 + P * vetU0;
        const REAL FDPart3tmp16 = FDPart3tmp12 * rescaledvU1 + P * vetU1;
        const REAL FDPart3tmp21 = FDPart3tmp20 * vetU0;
        const REAL FDPart3tmp25 = FDPart3tmp17 * FDPart3tmp24;
        const REAL FDPart3tmp42 =
            FDPart3tmp39 * ((rescaledvU0) * (rescaledvU0)) + P * (FDPart3tmp40 * (FDPart3tmp34 * FDPart3tmp36 - ((hDD12) * (hDD12))) - FDPart3tmp41);
        const REAL FDPart3tmp44 = FDPart3tmp17 * FDPart3tmp43;
        const REAL FDPart3tmp46 = FDPart3tmp30 * FDPart3tmp45 * hDD01 + FDPart3tmp45 * aDD01;
        const REAL FDPart3tmp49 = FDPart3tmp12 * FDPart3tmp38 * rescaledvU0;
        const REAL FDPart3tmp61 = FDPart3tmp20 * FDPart3tmp57;
        const REAL FDPart3tmp63 =
            FDPart3tmp39 * ((rescaledvU1) * (rescaledvU1)) + P * (FDPart3tmp40 * (FDPart3tmp29 * FDPart3tmp34 - ((hDD02) * (hDD02))) - FDPart3tmp62);
        const REAL FDPart3tmp71 = FDPart3tmp30 * FDPart3tmp56 * FDPart3tmp70 + FDPart3tmp56 * FDPart3tmp6 * aDD12;
        const REAL FDPart3tmp83 =
            FDPart3tmp39 * ((rescaledvU2) * (rescaledvU2)) + P * (FDPart3tmp40 * (FDPart3tmp29 * FDPart3tmp36 - ((hDD01) * (hDD01))) - FDPart3tmp82);
        const REAL FDPart3tmp91 = FDPart3tmp24 * FDPart3tmp38;
        const REAL FDPart3tmp112 = 2 * FDPart3tmp73;
        const REAL FDPart3tmp116 = 2 * FDPart3tmp50;
        const REAL FDPart3tmp133 = FDPart3tmp2 * FDPart3tmp23 * FDPart3tmp43 / (FDPart3tmp98 * alpha);
        const REAL FDPart3tmp18 = FDPart3tmp14 * FDPart3tmp17;
        const REAL FDPart3tmp52 = FDPart3tmp49 * rescaledvU1 + P * (FDPart3tmp40 * (-FDPart3tmp34 * hDD01 + hDD02 * hDD12) - FDPart3tmp50 * vetU1);
        const REAL FDPart3tmp64 = FDPart3tmp44 * FDPart3tmp57;
        const REAL FDPart3tmp66 = FDPart3tmp30 * FDPart3tmp45 * FDPart3tmp6 * hDD02 + FDPart3tmp45 * FDPart3tmp6 * aDD02;
        const REAL FDPart3tmp68 = FDPart3tmp49 * rescaledvU2 + P * (FDPart3tmp40 * (-FDPart3tmp36 * hDD02 + hDD01 * hDD12) - FDPart3tmp50 * vetU2);
        const REAL FDPart3tmp75 = FDPart3tmp12 * FDPart3tmp38 * rescaledvU1 * rescaledvU2 +
                                  P * (FDPart3tmp40 * (-FDPart3tmp29 * hDD12 + hDD01 * hDD02) - FDPart3tmp73 * vetU2);
        const REAL FDPart3tmp85 = FDPart3tmp16 * FDPart3tmp38;
        const REAL FDPart3tmp97 = FDPart3tmp14 * FDPart3tmp38;
        const REAL FDPart3tmp100 = FDPart3tmp14 * hDD01 + FDPart3tmp16 * FDPart3tmp36 + FDPart3tmp19 * FDPart3tmp84 + FDPart3tmp24 * hDD12;
        const REAL FDPart3tmp102 = FDPart3tmp14 * hDD02 + FDPart3tmp16 * hDD12 + FDPart3tmp19 * FDPart3tmp90 + FDPart3tmp24 * FDPart3tmp34;
        const REAL FDPart3tmp106 = FDPart3tmp14 * FDPart3tmp29 + FDPart3tmp16 * hDD01 + FDPart3tmp19 * FDPart3tmp94 + FDPart3tmp24 * hDD02;
        const REAL FDPart3tmp127 = FDPart3tmp116 * FDPart3tmp14 + FDPart3tmp19 * FDPart3tmp41 + FDPart3tmp42;
        const REAL FDPart3tmp131 = FDPart3tmp19 * FDPart3tmp82 + FDPart3tmp77 * FDPart3tmp91 + FDPart3tmp83;
        const REAL FDPart3tmp132 = FDPart3tmp112 * FDPart3tmp16 + FDPart3tmp19 * FDPart3tmp62 + FDPart3tmp63;
        const REAL FDPart3tmp53 = FDPart3tmp0 * FDPart3tmp21 * vetU1 + FDPart3tmp0 * FDPart3tmp44 * FDPart3tmp52;
        const REAL FDPart3tmp58 = FDPart3tmp16 * FDPart3tmp17 * FDPart3tmp57;
        const REAL FDPart3tmp69 = FDPart3tmp21 * FDPart3tmp23 * vetU2 + FDPart3tmp23 * FDPart3tmp44 * FDPart3tmp68;
        const REAL FDPart3tmp76 = FDPart3tmp26 * FDPart3tmp61 * vetU1 + FDPart3tmp64 * FDPart3tmp7 * FDPart3tmp75;
        const REAL FDPart3tmp95 = FDPart3tmp29 * FDPart3tmp52 + FDPart3tmp63 * hDD01 + FDPart3tmp75 * hDD02 + FDPart3tmp85 * FDPart3tmp94;
        const REAL FDPart3tmp104 = FDPart3tmp102 * FDPart3tmp38 * cf;
        const REAL FDPart3tmp113 = FDPart3tmp19 * FDPart3tmp73 * vetU2 + FDPart3tmp75;
        const REAL FDPart3tmp117 = FDPart3tmp19 * FDPart3tmp50 * vetU1 + FDPart3tmp52;
        const REAL FDPart3tmp123 = FDPart3tmp19 * FDPart3tmp50 * vetU2 + FDPart3tmp68;
        const REAL FDPart3tmp88 = FDPart3tmp36 * FDPart3tmp63 + FDPart3tmp52 * hDD01 + FDPart3tmp75 * hDD12 + FDPart3tmp84 * FDPart3tmp85;
        const REAL FDPart3tmp93 = FDPart3tmp34 * FDPart3tmp83 + FDPart3tmp68 * hDD02 + FDPart3tmp75 * hDD12 + FDPart3tmp90 * FDPart3tmp91;
        const REAL FDPart3tmp108 = FDPart3tmp100 * FDPart3tmp38 * cf;
        const REAL FDPart3tmp114 = FDPart3tmp112 * FDPart3tmp24 + FDPart3tmp113;
        const REAL FDPart3tmp115 = FDPart3tmp113 + FDPart3tmp77 * FDPart3tmp85;
        const REAL FDPart3tmp118 = FDPart3tmp116 * FDPart3tmp16 + FDPart3tmp117;
        const REAL FDPart3tmp122 = FDPart3tmp112 * FDPart3tmp14 + FDPart3tmp117;
        const REAL FDPart3tmp124 = FDPart3tmp116 * FDPart3tmp24 + FDPart3tmp123;
        const REAL FDPart3tmp126 = FDPart3tmp123 + FDPart3tmp77 * FDPart3tmp97;
        const REAL FDPart3tmp136 = FDPart3tmp106 * FDPart3tmp38 * FDPart3tmp89;
        rhs_gfs[IDX4(RHO_STARGF, i0, i1, i2)] = -FDPart3tmp1 * FDPart3tmp3 * rescaledvU0 - FDPart3tmp5 * FDPart3tmp7 * FDPart3tmp9;
        rhs_gfs[IDX4(TAU_TILDEGF, i0, i1, i2)] =
            -FDPart3tmp1 * (FDPart3tmp14 * FDPart3tmp2 - FDPart3tmp3 * rescaledvU0) +
            FDPart3tmp2 * alpha *
                (FDPart3tmp46 * (FDPart3tmp1 * FDPart3tmp18 * vetU1 + FDPart3tmp53) +
                 FDPart3tmp46 * (FDPart3tmp1 * FDPart3tmp16 * FDPart3tmp17 * vetU0 + FDPart3tmp53) +
                 FDPart3tmp66 * (FDPart3tmp1 * FDPart3tmp18 * FDPart3tmp26 + FDPart3tmp69) +
                 FDPart3tmp66 * (FDPart3tmp1 * FDPart3tmp25 * FDPart3tmp7 * vetU0 + FDPart3tmp69) +
                 FDPart3tmp71 * (FDPart3tmp58 * FDPart3tmp7 * FDPart3tmp77 + FDPart3tmp76) +
                 FDPart3tmp71 * (2 * FDPart3tmp25 * FDPart3tmp57 * FDPart3tmp7 * vetU1 + FDPart3tmp76) + alpha_dD0 * (-FDPart3tmp18 - FDPart3tmp21) +
                 alpha_dD1 * (-FDPart3tmp0 * FDPart3tmp16 * FDPart3tmp17 - FDPart3tmp0 * FDPart3tmp20 * vetU1) +
                 alpha_dD2 * (-FDPart3tmp0 * FDPart3tmp20 * FDPart3tmp26 - FDPart3tmp23 * FDPart3tmp25) +
                 (FDPart3tmp28 * aDD00 + FDPart3tmp29 * FDPart3tmp31) *
                     (2 * FDPart3tmp18 * vetU0 + FDPart3tmp20 * ((vetU0) * (vetU0)) + FDPart3tmp42 * FDPart3tmp44) +
                 (FDPart3tmp31 * (FDPart3tmp55 * hDD11 + FDPart3tmp55) + FDPart3tmp56 * aDD11) *
                     (2 * FDPart3tmp58 * vetU1 + FDPart3tmp61 * ((vetU1) * (vetU1)) + FDPart3tmp63 * FDPart3tmp64) +
                 (FDPart3tmp28 * FDPart3tmp79 * aDD22 + FDPart3tmp31 * (FDPart3tmp79 * hDD22 + FDPart3tmp79)) *
                     (FDPart3tmp25 * FDPart3tmp57 * FDPart3tmp77 * FDPart3tmp80 + FDPart3tmp61 * FDPart3tmp80 * ((vetU2) * (vetU2)) +
                      FDPart3tmp64 * FDPart3tmp80 * FDPart3tmp83)) -
            FDPart3tmp5 * FDPart3tmp7 * (FDPart3tmp0 * FDPart3tmp16 * FDPart3tmp2 - FDPart3tmp9);
        rhs_gfs[IDX4(RESCALEDSTILDED0GF, i0, i1, i2)] =
            FDPart3tmp133 *
            (FDPart3tmp6 * (FDPart3tmp100 * FDPart3tmp73 * cf + FDPart3tmp104 * vetU2 + FDPart3tmp104 * (-vetU2 + vetU_dD20 * xx0) +
                            FDPart3tmp105 * FDPart3tmp106 * FDPart3tmp38 * vetU_dD00 + FDPart3tmp108 * (-vetU1 + vetU_dD10 * xx0) +
                            FDPart3tmp111 * FDPart3tmp114 + FDPart3tmp111 * FDPart3tmp115 + FDPart3tmp118 * FDPart3tmp121 +
                            FDPart3tmp121 * FDPart3tmp122 + FDPart3tmp124 * FDPart3tmp125 + FDPart3tmp125 * FDPart3tmp126 +
                            0.5 * FDPart3tmp127 * xx0 * (-2 * FDPart3tmp29 * cf_dD0 + cf * hDD_dD000) +
                            FDPart3tmp131 * (FDPart3tmp109 * hDD_dD220 - FDPart3tmp110 * FDPart3tmp34 - FDPart3tmp128 * FDPart3tmp34 +
                                             FDPart3tmp128 * hDD22 + FDPart3tmp128) +
                            FDPart3tmp132 * (FDPart3tmp109 * hDD_dD110 - FDPart3tmp110 * FDPart3tmp36 - FDPart3tmp128 * FDPart3tmp36 +
                                             FDPart3tmp128 * hDD11 + FDPart3tmp128) -
                            FDPart3tmp99 * alpha_dD0 * xx0) +
             FDPart3tmp88 * FDPart3tmp89 + FDPart3tmp89 * FDPart3tmp93 -
             2 * FDPart3tmp89 * (FDPart3tmp29 * FDPart3tmp42 + FDPart3tmp52 * hDD01 + FDPart3tmp68 * hDD02 + FDPart3tmp94 * FDPart3tmp97) -
             FDPart3tmp95 * FDPart3tmp96);
        rhs_gfs[IDX4(RESCALEDSTILDED1GF, i0, i1, i2)] =
            FDPart3tmp133 *
            (FDPart3tmp100 * FDPart3tmp38 * FDPart3tmp89 * vetU_dD11 + FDPart3tmp100 * FDPart3tmp50 * FDPart3tmp89 + FDPart3tmp104 * FDPart3tmp135 +
             FDPart3tmp104 * (-FDPart3tmp135 + FDPart3tmp6 * vetU_dD21) - FDPart3tmp106 * FDPart3tmp73 * FDPart3tmp89 +
             FDPart3tmp114 * FDPart3tmp138 + FDPart3tmp115 * FDPart3tmp138 + FDPart3tmp118 * FDPart3tmp147 + FDPart3tmp122 * FDPart3tmp147 +
             FDPart3tmp124 * FDPart3tmp140 + FDPart3tmp126 * FDPart3tmp140 +
             FDPart3tmp127 * FDPart3tmp6 * (-FDPart3tmp128 * hDD01 - FDPart3tmp142 * FDPart3tmp29 + 0.5 * cf * hDD_dD001) +
             FDPart3tmp131 * (-FDPart3tmp128 * FDPart3tmp34 * FDPart3tmp5 + FDPart3tmp128 * FDPart3tmp5 * hDD22 + FDPart3tmp128 * FDPart3tmp5 -
                              FDPart3tmp142 * FDPart3tmp34 * FDPart3tmp6 + 0.5 * FDPart3tmp89 * hDD_dD221) +
             FDPart3tmp132 * FDPart3tmp6 * (FDPart3tmp128 * hDD01 - FDPart3tmp142 * FDPart3tmp36 + FDPart3tmp143 * hDD_dD111) +
             FDPart3tmp136 * vetU_dD01 - FDPart3tmp6 * FDPart3tmp99 * alpha_dD1 - FDPart3tmp88 * FDPart3tmp96 - FDPart3tmp89 * FDPart3tmp95 -
             FDPart3tmp89 * (FDPart3tmp36 * FDPart3tmp52 + FDPart3tmp42 * hDD01 + FDPart3tmp68 * hDD12 + FDPart3tmp84 * FDPart3tmp97) +
             FDPart3tmp93 * FDPart3tmp96);
        rhs_gfs[IDX4(RESCALEDSTILDED2GF, i0, i1, i2)] =
            FDPart3tmp133 *
            (FDPart3tmp102 * FDPart3tmp50 * FDPart3tmp89 + FDPart3tmp102 * FDPart3tmp73 * FDPart3tmp96 + FDPart3tmp104 * vetU_dD22 +
             FDPart3tmp106 * FDPart3tmp38 * cf * vetU_dD02 - FDPart3tmp108 * FDPart3tmp135 + FDPart3tmp108 * vetU_dD12 +
             FDPart3tmp114 * FDPart3tmp154 + FDPart3tmp115 * FDPart3tmp154 + FDPart3tmp118 * FDPart3tmp149 + FDPart3tmp122 * FDPart3tmp149 +
             FDPart3tmp124 * FDPart3tmp155 + FDPart3tmp126 * FDPart3tmp155 +
             FDPart3tmp127 * (-FDPart3tmp150 * FDPart3tmp29 - FDPart3tmp151 + 0.5 * cf * hDD_dD002) +
             FDPart3tmp131 * (FDPart3tmp143 * hDD_dD222 - FDPart3tmp150 * FDPart3tmp34 + FDPart3tmp151 + FDPart3tmp152) +
             FDPart3tmp132 * (-FDPart3tmp150 * FDPart3tmp36 - FDPart3tmp152 + 0.5 * cf * hDD_dD112) - FDPart3tmp136 * vetU2 -
             FDPart3tmp89 * (FDPart3tmp29 * FDPart3tmp68 + FDPart3tmp75 * hDD01 + FDPart3tmp83 * hDD02 + FDPart3tmp91 * FDPart3tmp94) -
             FDPart3tmp89 * (FDPart3tmp34 * FDPart3tmp68 + FDPart3tmp42 * hDD02 + FDPart3tmp52 * hDD12 + FDPart3tmp90 * FDPart3tmp97) -
             FDPart3tmp96 * (FDPart3tmp36 * FDPart3tmp75 + FDPart3tmp68 * hDD01 + FDPart3tmp83 * hDD12 + FDPart3tmp84 * FDPart3tmp91) -
             FDPart3tmp99 * alpha_dD2);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION calculate_all_source_terms__rfm__Spherical
