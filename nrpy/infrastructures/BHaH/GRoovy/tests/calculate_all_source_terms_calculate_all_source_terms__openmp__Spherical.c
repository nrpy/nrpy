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
        ghl_initialize_primitives(auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)], auxevol_gfs[IDX4(PGF, i0, i1, i2)], NAN, NAN, NAN, NAN, NAN, NAN, NAN,
                                  auxevol_gfs[IDX4(SGF, i0, i1, i2)], auxevol_gfs[IDX4(YEGF, i0, i1, i2)],
                                  auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)], &prims);

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
        const REAL FDPart3tmp7 = cos(xx1);
        const REAL FDPart3tmp8 = sin(xx1);
        const REAL FDPart3tmp14 = ((alpha) * (alpha));
        const REAL FDPart3tmp33 = hDD00 + 1;
        const REAL FDPart3tmp34 = (1.0 / 3.0) * trK;
        const REAL FDPart3tmp38 = hDD22 + 1;
        const REAL FDPart3tmp40 = hDD11 + 1;
        const REAL FDPart3tmp59 = ((xx0) * (xx0));
        const REAL FDPart3tmp81 = 2 * vetU2;
        const REAL FDPart3tmp102 = ((cf) * (cf) * (cf));
        const REAL FDPart3tmp109 = cf * xx0;
        const REAL FDPart3tmp114 = 1.0 * cf_dD0 * xx0;
        const REAL FDPart3tmp132 = 1.0 * cf;
        const REAL FDPart3tmp146 = 1.0 * cf_dD1;
        const REAL FDPart3tmp147 = 0.5 * cf;
        const REAL FDPart3tmp154 = 1.0 * cf_dD2;
        const REAL FDPart3tmp1 = 2 * FDPart3tmp0;
        const REAL FDPart3tmp3 = FDPart3tmp2 * alpha * u4Ut;
        const REAL FDPart3tmp9 = (1.0 / (FDPart3tmp8));
        const REAL FDPart3tmp16 = FDPart3tmp14 * h * rhob * ((u4Ut) * (u4Ut));
        const REAL FDPart3tmp21 = (1.0 / (FDPart3tmp14));
        const REAL FDPart3tmp23 = FDPart3tmp14 * h * rhob * ((u4Ut) * (u4Ut)) - P;
        const REAL FDPart3tmp32 = (1.0 / ((cf) * (cf)));
        const REAL FDPart3tmp42 = FDPart3tmp33 * FDPart3tmp38 * FDPart3tmp40 - FDPart3tmp33 * ((hDD12) * (hDD12)) -
                                  FDPart3tmp38 * ((hDD01) * (hDD01)) - FDPart3tmp40 * ((hDD02) * (hDD02)) + 2 * hDD01 * hDD02 * hDD12;
        const REAL FDPart3tmp44 = FDPart3tmp14 * ((cf) * (cf));
        const REAL FDPart3tmp61 = (1.0 / (FDPart3tmp59));
        const REAL FDPart3tmp74 = FDPart3tmp8 * hDD12;
        const REAL FDPart3tmp88 = FDPart3tmp40 * vetU1 + hDD01 * vetU0 + hDD12 * vetU2;
        const REAL FDPart3tmp93 = FDPart3tmp8 * cf;
        const REAL FDPart3tmp94 = FDPart3tmp38 * vetU2 + hDD02 * vetU0 + hDD12 * vetU1;
        const REAL FDPart3tmp98 = FDPart3tmp33 * vetU0 + hDD01 * vetU1 + hDD02 * vetU2;
        const REAL FDPart3tmp100 = FDPart3tmp7 * cf;
        const REAL FDPart3tmp113 = 0.5 * FDPart3tmp109;
        const REAL FDPart3tmp124 = 2 * cf_dD0 * xx0;
        const REAL FDPart3tmp139 = FDPart3tmp7 * vetU2;
        const REAL FDPart3tmp10 = FDPart3tmp7 * FDPart3tmp9;
        const REAL FDPart3tmp24 = FDPart3tmp21 * FDPart3tmp23;
        const REAL FDPart3tmp27 = FDPart3tmp0 * FDPart3tmp9;
        const REAL FDPart3tmp28 = FDPart3tmp16 * rescaledvU2 + P * vetU2;
        const REAL FDPart3tmp30 = FDPart3tmp9 * vetU2;
        const REAL FDPart3tmp35 = FDPart3tmp32 * FDPart3tmp34;
        const REAL FDPart3tmp43 = FDPart3tmp16 * FDPart3tmp42;
        const REAL FDPart3tmp45 = FDPart3tmp42 * ((vetU0) * (vetU0));
        const REAL FDPart3tmp47 = (1.0 / (FDPart3tmp42));
        const REAL FDPart3tmp49 = FDPart3tmp32 * xx0;
        const REAL FDPart3tmp54 = FDPart3tmp42 * vetU0;
        const REAL FDPart3tmp60 = FDPart3tmp32 * FDPart3tmp59;
        const REAL FDPart3tmp66 = FDPart3tmp42 * ((vetU1) * (vetU1));
        const REAL FDPart3tmp77 = FDPart3tmp42 * vetU1;
        const REAL FDPart3tmp83 = FDPart3tmp59 * ((FDPart3tmp8) * (FDPart3tmp8));
        const REAL FDPart3tmp84 = (1.0 / ((FDPart3tmp8) * (FDPart3tmp8)));
        const REAL FDPart3tmp86 = FDPart3tmp42 * ((vetU2) * (vetU2));
        const REAL FDPart3tmp103 = FDPart3tmp102 * FDPart3tmp23 * FDPart3tmp42 * alpha;
        const REAL FDPart3tmp115 = FDPart3tmp113 * hDD_dD120 - FDPart3tmp114 * hDD12;
        const REAL FDPart3tmp125 = 0.5 * FDPart3tmp109 * hDD_dD010 - 0.5 * FDPart3tmp124 * hDD01;
        const REAL FDPart3tmp129 = 0.5 * FDPart3tmp109 * hDD_dD020 - 0.5 * FDPart3tmp124 * hDD02;
        const REAL FDPart3tmp142 = -1.0 * FDPart3tmp74 * cf_dD1 + 0.5 * FDPart3tmp93 * hDD02 + 0.5 * FDPart3tmp93 * hDD_dD121;
        const REAL FDPart3tmp144 = -0.5 * FDPart3tmp74 * cf + 0.5 * FDPart3tmp8 * cf * hDD_dD021 - 1.0 * FDPart3tmp8 * cf_dD1 * hDD02;
        const REAL FDPart3tmp151 = 0.5 * FDPart3tmp8 * (FDPart3tmp33 * cf - FDPart3tmp40 * cf + cf * hDD_dD011 - 2 * cf_dD1 * hDD01);
        const REAL FDPart3tmp153 = 0.5 * cf * hDD_dD012 - 0.5 * cf * (FDPart3tmp7 * hDD02 + FDPart3tmp74) - 1.0 * cf_dD2 * hDD01;
        const REAL FDPart3tmp155 = FDPart3tmp132 * FDPart3tmp8 * hDD02;
        const REAL FDPart3tmp156 = FDPart3tmp132 * FDPart3tmp7 * hDD12;
        const REAL FDPart3tmp158 = -0.5 * FDPart3tmp100 * FDPart3tmp38 + 0.5 * FDPart3tmp40 * FDPart3tmp7 * cf + 0.5 * FDPart3tmp93 * hDD01 +
                                   0.5 * cf * hDD_dD122 - 1.0 * cf_dD2 * hDD12;
        const REAL FDPart3tmp159 = 0.5 * FDPart3tmp100 * hDD01 + 0.5 * FDPart3tmp33 * FDPart3tmp8 * cf - 0.5 * FDPart3tmp38 * FDPart3tmp93 +
                                   0.5 * cf * hDD_dD022 - 1.0 * cf_dD2 * hDD02;
        const REAL FDPart3tmp5 = FDPart3tmp3 * rescaledvU0 * rhob;
        const REAL FDPart3tmp12 = FDPart3tmp0 * FDPart3tmp3 * rescaledvU1 * rhob;
        const REAL FDPart3tmp18 = FDPart3tmp16 * rescaledvU0 + P * vetU0;
        const REAL FDPart3tmp20 = FDPart3tmp16 * rescaledvU1 + P * vetU1;
        const REAL FDPart3tmp25 = FDPart3tmp24 * vetU0;
        const REAL FDPart3tmp29 = FDPart3tmp21 * FDPart3tmp28;
        const REAL FDPart3tmp46 =
            FDPart3tmp43 * ((rescaledvU0) * (rescaledvU0)) + P * (FDPart3tmp44 * (FDPart3tmp38 * FDPart3tmp40 - ((hDD12) * (hDD12))) - FDPart3tmp45);
        const REAL FDPart3tmp48 = FDPart3tmp21 * FDPart3tmp47;
        const REAL FDPart3tmp50 = FDPart3tmp34 * FDPart3tmp49 * hDD01 + FDPart3tmp49 * aDD01;
        const REAL FDPart3tmp53 = FDPart3tmp16 * FDPart3tmp42 * rescaledvU0;
        const REAL FDPart3tmp65 = FDPart3tmp24 * FDPart3tmp61;
        const REAL FDPart3tmp67 =
            FDPart3tmp43 * ((rescaledvU1) * (rescaledvU1)) + P * (FDPart3tmp44 * (FDPart3tmp33 * FDPart3tmp38 - ((hDD02) * (hDD02))) - FDPart3tmp66);
        const REAL FDPart3tmp75 = FDPart3tmp34 * FDPart3tmp60 * FDPart3tmp74 + FDPart3tmp60 * FDPart3tmp8 * aDD12;
        const REAL FDPart3tmp87 =
            FDPart3tmp43 * ((rescaledvU2) * (rescaledvU2)) + P * (FDPart3tmp44 * (FDPart3tmp33 * FDPart3tmp40 - ((hDD01) * (hDD01))) - FDPart3tmp86);
        const REAL FDPart3tmp95 = FDPart3tmp28 * FDPart3tmp42;
        const REAL FDPart3tmp116 = 2 * FDPart3tmp77;
        const REAL FDPart3tmp120 = 2 * FDPart3tmp54;
        const REAL FDPart3tmp137 = FDPart3tmp2 * FDPart3tmp27 * FDPart3tmp47 / (FDPart3tmp102 * alpha);
        const REAL FDPart3tmp22 = FDPart3tmp18 * FDPart3tmp21;
        const REAL FDPart3tmp56 = FDPart3tmp53 * rescaledvU1 + P * (FDPart3tmp44 * (-FDPart3tmp38 * hDD01 + hDD02 * hDD12) - FDPart3tmp54 * vetU1);
        const REAL FDPart3tmp68 = FDPart3tmp48 * FDPart3tmp61;
        const REAL FDPart3tmp70 = FDPart3tmp34 * FDPart3tmp49 * FDPart3tmp8 * hDD02 + FDPart3tmp49 * FDPart3tmp8 * aDD02;
        const REAL FDPart3tmp72 = FDPart3tmp53 * rescaledvU2 + P * (FDPart3tmp44 * (-FDPart3tmp40 * hDD02 + hDD01 * hDD12) - FDPart3tmp54 * vetU2);
        const REAL FDPart3tmp79 = FDPart3tmp16 * FDPart3tmp42 * rescaledvU1 * rescaledvU2 +
                                  P * (FDPart3tmp44 * (-FDPart3tmp33 * hDD12 + hDD01 * hDD02) - FDPart3tmp77 * vetU2);
        const REAL FDPart3tmp89 = FDPart3tmp20 * FDPart3tmp42;
        const REAL FDPart3tmp101 = FDPart3tmp18 * FDPart3tmp42;
        const REAL FDPart3tmp104 = FDPart3tmp18 * hDD01 + FDPart3tmp20 * FDPart3tmp40 + FDPart3tmp23 * FDPart3tmp88 + FDPart3tmp28 * hDD12;
        const REAL FDPart3tmp106 = FDPart3tmp18 * hDD02 + FDPart3tmp20 * hDD12 + FDPart3tmp23 * FDPart3tmp94 + FDPart3tmp28 * FDPart3tmp38;
        const REAL FDPart3tmp110 = FDPart3tmp18 * FDPart3tmp33 + FDPart3tmp20 * hDD01 + FDPart3tmp23 * FDPart3tmp98 + FDPart3tmp28 * hDD02;
        const REAL FDPart3tmp131 = FDPart3tmp120 * FDPart3tmp18 + FDPart3tmp23 * FDPart3tmp45 + FDPart3tmp46;
        const REAL FDPart3tmp135 = FDPart3tmp23 * FDPart3tmp86 + FDPart3tmp81 * FDPart3tmp95 + FDPart3tmp87;
        const REAL FDPart3tmp136 = FDPart3tmp116 * FDPart3tmp20 + FDPart3tmp23 * FDPart3tmp66 + FDPart3tmp67;
        const REAL FDPart3tmp57 = FDPart3tmp0 * FDPart3tmp25 * vetU1 + FDPart3tmp0 * FDPart3tmp48 * FDPart3tmp56;
        const REAL FDPart3tmp62 = FDPart3tmp20 * FDPart3tmp21 * FDPart3tmp61;
        const REAL FDPart3tmp73 = FDPart3tmp25 * FDPart3tmp27 * vetU2 + FDPart3tmp27 * FDPart3tmp48 * FDPart3tmp72;
        const REAL FDPart3tmp80 = FDPart3tmp30 * FDPart3tmp65 * vetU1 + FDPart3tmp68 * FDPart3tmp79 * FDPart3tmp9;
        const REAL FDPart3tmp99 = FDPart3tmp33 * FDPart3tmp56 + FDPart3tmp67 * hDD01 + FDPart3tmp79 * hDD02 + FDPart3tmp89 * FDPart3tmp98;
        const REAL FDPart3tmp108 = FDPart3tmp106 * FDPart3tmp42 * cf;
        const REAL FDPart3tmp117 = FDPart3tmp23 * FDPart3tmp77 * vetU2 + FDPart3tmp79;
        const REAL FDPart3tmp121 = FDPart3tmp23 * FDPart3tmp54 * vetU1 + FDPart3tmp56;
        const REAL FDPart3tmp127 = FDPart3tmp23 * FDPart3tmp54 * vetU2 + FDPart3tmp72;
        const REAL FDPart3tmp92 = FDPart3tmp40 * FDPart3tmp67 + FDPart3tmp56 * hDD01 + FDPart3tmp79 * hDD12 + FDPart3tmp88 * FDPart3tmp89;
        const REAL FDPart3tmp97 = FDPart3tmp38 * FDPart3tmp87 + FDPart3tmp72 * hDD02 + FDPart3tmp79 * hDD12 + FDPart3tmp94 * FDPart3tmp95;
        const REAL FDPart3tmp112 = FDPart3tmp104 * FDPart3tmp42 * cf;
        const REAL FDPart3tmp118 = FDPart3tmp116 * FDPart3tmp28 + FDPart3tmp117;
        const REAL FDPart3tmp119 = FDPart3tmp117 + FDPart3tmp81 * FDPart3tmp89;
        const REAL FDPart3tmp122 = FDPart3tmp120 * FDPart3tmp20 + FDPart3tmp121;
        const REAL FDPart3tmp126 = FDPart3tmp116 * FDPart3tmp18 + FDPart3tmp121;
        const REAL FDPart3tmp128 = FDPart3tmp120 * FDPart3tmp28 + FDPart3tmp127;
        const REAL FDPart3tmp130 = FDPart3tmp101 * FDPart3tmp81 + FDPart3tmp127;
        const REAL FDPart3tmp140 = FDPart3tmp110 * FDPart3tmp42 * FDPart3tmp93;
        rhs_gfs[IDX4(RHO_STARGF, i0, i1, i2)] = -FDPart3tmp1 * FDPart3tmp5 - FDPart3tmp10 * FDPart3tmp12;
        rhs_gfs[IDX4(TAU_TILDEGF, i0, i1, i2)] =
            -FDPart3tmp1 * (FDPart3tmp18 * FDPart3tmp2 - FDPart3tmp5) - FDPart3tmp10 * (FDPart3tmp0 * FDPart3tmp2 * FDPart3tmp20 - FDPart3tmp12) +
            FDPart3tmp2 * alpha *
                (FDPart3tmp50 * (FDPart3tmp1 * FDPart3tmp22 * vetU1 + FDPart3tmp57) +
                 FDPart3tmp50 * (FDPart3tmp1 * FDPart3tmp20 * FDPart3tmp21 * vetU0 + FDPart3tmp57) +
                 FDPart3tmp70 * (FDPart3tmp1 * FDPart3tmp22 * FDPart3tmp30 + FDPart3tmp73) +
                 FDPart3tmp70 * (FDPart3tmp1 * FDPart3tmp29 * FDPart3tmp9 * vetU0 + FDPart3tmp73) +
                 FDPart3tmp75 * (FDPart3tmp62 * FDPart3tmp81 * FDPart3tmp9 + FDPart3tmp80) +
                 FDPart3tmp75 * (2 * FDPart3tmp29 * FDPart3tmp61 * FDPart3tmp9 * vetU1 + FDPart3tmp80) + alpha_dD0 * (-FDPart3tmp22 - FDPart3tmp25) +
                 alpha_dD1 * (-FDPart3tmp0 * FDPart3tmp20 * FDPart3tmp21 - FDPart3tmp0 * FDPart3tmp24 * vetU1) +
                 alpha_dD2 * (-FDPart3tmp0 * FDPart3tmp24 * FDPart3tmp30 - FDPart3tmp27 * FDPart3tmp29) +
                 (FDPart3tmp32 * aDD00 + FDPart3tmp33 * FDPart3tmp35) *
                     (2 * FDPart3tmp22 * vetU0 + FDPart3tmp24 * ((vetU0) * (vetU0)) + FDPart3tmp46 * FDPart3tmp48) +
                 (FDPart3tmp35 * (FDPart3tmp59 * hDD11 + FDPart3tmp59) + FDPart3tmp60 * aDD11) *
                     (2 * FDPart3tmp62 * vetU1 + FDPart3tmp65 * ((vetU1) * (vetU1)) + FDPart3tmp67 * FDPart3tmp68) +
                 (FDPart3tmp32 * FDPart3tmp83 * aDD22 + FDPart3tmp35 * (FDPart3tmp83 * hDD22 + FDPart3tmp83)) *
                     (FDPart3tmp29 * FDPart3tmp61 * FDPart3tmp81 * FDPart3tmp84 + FDPart3tmp65 * FDPart3tmp84 * ((vetU2) * (vetU2)) +
                      FDPart3tmp68 * FDPart3tmp84 * FDPart3tmp87));
        rhs_gfs[IDX4(RESCALEDSTILDED0GF, i0, i1, i2)] =
            FDPart3tmp137 *
            (-FDPart3tmp100 * FDPart3tmp99 +
             FDPart3tmp8 * (-FDPart3tmp103 * alpha_dD0 * xx0 + FDPart3tmp104 * FDPart3tmp77 * cf + FDPart3tmp108 * vetU2 +
                            FDPart3tmp108 * (-vetU2 + vetU_dD20 * xx0) + FDPart3tmp109 * FDPart3tmp110 * FDPart3tmp42 * vetU_dD00 +
                            FDPart3tmp112 * (-vetU1 + vetU_dD10 * xx0) + FDPart3tmp115 * FDPart3tmp118 + FDPart3tmp115 * FDPart3tmp119 +
                            FDPart3tmp122 * FDPart3tmp125 + FDPart3tmp125 * FDPart3tmp126 + FDPart3tmp128 * FDPart3tmp129 +
                            FDPart3tmp129 * FDPart3tmp130 + 0.5 * FDPart3tmp131 * xx0 * (-2 * FDPart3tmp33 * cf_dD0 + cf * hDD_dD000) +
                            FDPart3tmp135 * (FDPart3tmp113 * hDD_dD220 - FDPart3tmp114 * FDPart3tmp38 - FDPart3tmp132 * FDPart3tmp38 +
                                             FDPart3tmp132 * hDD22 + FDPart3tmp132) +
                            FDPart3tmp136 * (FDPart3tmp113 * hDD_dD110 - FDPart3tmp114 * FDPart3tmp40 - FDPart3tmp132 * FDPart3tmp40 +
                                             FDPart3tmp132 * hDD11 + FDPart3tmp132)) +
             FDPart3tmp92 * FDPart3tmp93 + FDPart3tmp93 * FDPart3tmp97 -
             2 * FDPart3tmp93 * (FDPart3tmp101 * FDPart3tmp98 + FDPart3tmp33 * FDPart3tmp46 + FDPart3tmp56 * hDD01 + FDPart3tmp72 * hDD02));
        rhs_gfs[IDX4(RESCALEDSTILDED1GF, i0, i1, i2)] =
            FDPart3tmp137 *
            (-FDPart3tmp100 * FDPart3tmp92 + FDPart3tmp100 * FDPart3tmp97 - FDPart3tmp103 * FDPart3tmp8 * alpha_dD1 +
             FDPart3tmp104 * FDPart3tmp42 * FDPart3tmp93 * vetU_dD11 + FDPart3tmp104 * FDPart3tmp54 * FDPart3tmp93 + FDPart3tmp108 * FDPart3tmp139 +
             FDPart3tmp108 * (-FDPart3tmp139 + FDPart3tmp8 * vetU_dD21) - FDPart3tmp110 * FDPart3tmp77 * FDPart3tmp93 +
             FDPart3tmp118 * FDPart3tmp142 + FDPart3tmp119 * FDPart3tmp142 + FDPart3tmp122 * FDPart3tmp151 + FDPart3tmp126 * FDPart3tmp151 +
             FDPart3tmp128 * FDPart3tmp144 + FDPart3tmp130 * FDPart3tmp144 +
             FDPart3tmp131 * FDPart3tmp8 * (-FDPart3tmp132 * hDD01 - FDPart3tmp146 * FDPart3tmp33 + 0.5 * cf * hDD_dD001) +
             FDPart3tmp135 * (-FDPart3tmp132 * FDPart3tmp38 * FDPart3tmp7 + FDPart3tmp132 * FDPart3tmp7 * hDD22 + FDPart3tmp132 * FDPart3tmp7 -
                              FDPart3tmp146 * FDPart3tmp38 * FDPart3tmp8 + 0.5 * FDPart3tmp93 * hDD_dD221) +
             FDPart3tmp136 * FDPart3tmp8 * (FDPart3tmp132 * hDD01 - FDPart3tmp146 * FDPart3tmp40 + FDPart3tmp147 * hDD_dD111) +
             FDPart3tmp140 * vetU_dD01 - FDPart3tmp93 * FDPart3tmp99 -
             FDPart3tmp93 * (FDPart3tmp101 * FDPart3tmp88 + FDPart3tmp40 * FDPart3tmp56 + FDPart3tmp46 * hDD01 + FDPart3tmp72 * hDD12));
        rhs_gfs[IDX4(RESCALEDSTILDED2GF, i0, i1, i2)] =
            FDPart3tmp137 *
            (FDPart3tmp100 * FDPart3tmp106 * FDPart3tmp77 -
             FDPart3tmp100 * (FDPart3tmp40 * FDPart3tmp79 + FDPart3tmp72 * hDD01 + FDPart3tmp87 * hDD12 + FDPart3tmp88 * FDPart3tmp95) -
             FDPart3tmp103 * alpha_dD2 + FDPart3tmp106 * FDPart3tmp54 * FDPart3tmp93 + FDPart3tmp108 * vetU_dD22 +
             FDPart3tmp110 * FDPart3tmp42 * cf * vetU_dD02 - FDPart3tmp112 * FDPart3tmp139 + FDPart3tmp112 * vetU_dD12 +
             FDPart3tmp118 * FDPart3tmp158 + FDPart3tmp119 * FDPart3tmp158 + FDPart3tmp122 * FDPart3tmp153 + FDPart3tmp126 * FDPart3tmp153 +
             FDPart3tmp128 * FDPart3tmp159 + FDPart3tmp130 * FDPart3tmp159 +
             FDPart3tmp131 * (-FDPart3tmp154 * FDPart3tmp33 - FDPart3tmp155 + 0.5 * cf * hDD_dD002) +
             FDPart3tmp135 * (FDPart3tmp147 * hDD_dD222 - FDPart3tmp154 * FDPart3tmp38 + FDPart3tmp155 + FDPart3tmp156) +
             FDPart3tmp136 * (-FDPart3tmp154 * FDPart3tmp40 - FDPart3tmp156 + 0.5 * cf * hDD_dD112) - FDPart3tmp140 * vetU2 -
             FDPart3tmp93 * (FDPart3tmp101 * FDPart3tmp94 + FDPart3tmp38 * FDPart3tmp72 + FDPart3tmp46 * hDD02 + FDPart3tmp56 * hDD12) -
             FDPart3tmp93 * (FDPart3tmp33 * FDPart3tmp72 + FDPart3tmp79 * hDD01 + FDPart3tmp87 * hDD02 + FDPart3tmp95 * FDPart3tmp98));
        rhs_gfs[IDX4(YE_STARGF, i0, i1, i2)] = -FDPart3tmp1 * FDPart3tmp5 * Ye - FDPart3tmp10 * FDPart3tmp12 * Ye;
        rhs_gfs[IDX4(S_STARGF, i0, i1, i2)] =
            -FDPart3tmp0 * FDPart3tmp10 * FDPart3tmp3 * S * rescaledvU1 - FDPart3tmp1 * FDPart3tmp3 * S * rescaledvU0;

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION calculate_all_source_terms__rfm__Spherical
