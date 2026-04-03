#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Compute flux divergences for rho_star, tau_tilde, rescaledStildeD, and optional composition variables.
 */
void calculate_flux_divergences__rfm__Spherical(const int flux_dirn, const commondata_struct *restrict commondata,
                                                const params_struct *restrict params, REAL *restrict xx[3], const ghl_eos_parameters *restrict eos,
                                                REAL *restrict auxevol_gfs, const REAL *restrict in_gfs, REAL *restrict rhs_gfs) {
#include "set_CodeParameters.h"
  const REAL invdxi[3] = {invdxx0, invdxx1, invdxx2};
  const REAL invdx = invdxi[flux_dirn];

  const int xdir = (flux_dirn == 0);
  const int ydir = (flux_dirn == 1);
  const int zdir = (flux_dirn == 2);
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    MAYBE_UNUSED const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      MAYBE_UNUSED const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        MAYBE_UNUSED const REAL xx0 = xx[0][i0];

        const int index = IDX3(i0, i1, i2);
        const int indexp1 = IDX3(i0 + xdir, i1 + ydir, i2 + zdir);

        REAL ReU[3];
        REAL ReUdD[3][3];
        REAL ReUD[3][3];
        REAL ReUDdD[3][3][3];

        REAL rescaledStildefluxUD[3][3];
        REAL rescaledtau_tildefluxU[3];
        REAL rescaledrho_starfluxU[3];

        REAL rescaledYe_starfluxU[3];

        REAL rescaledS_starfluxU[3];

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
        const REAL alpha = in_gfs[IDX4(ALPHAGF, i0, i1, i2)];
        const REAL cf = in_gfs[IDX4(CFGF, i0, i1, i2)];
        const REAL hDD00 = in_gfs[IDX4(HDD00GF, i0, i1, i2)];
        const REAL hDD01 = in_gfs[IDX4(HDD01GF, i0, i1, i2)];
        const REAL hDD02 = in_gfs[IDX4(HDD02GF, i0, i1, i2)];
        const REAL hDD11 = in_gfs[IDX4(HDD11GF, i0, i1, i2)];
        const REAL hDD12 = in_gfs[IDX4(HDD12GF, i0, i1, i2)];
        const REAL hDD22 = in_gfs[IDX4(HDD22GF, i0, i1, i2)];
        const REAL vetU0 = in_gfs[IDX4(VETU0GF, i0, i1, i2)];
        const REAL vetU1 = in_gfs[IDX4(VETU1GF, i0, i1, i2)];
        const REAL vetU2 = in_gfs[IDX4(VETU2GF, i0, i1, i2)];

        /*
         * NRPy-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        const REAL FDPart3tmp0 = sin(xx1);
        const REAL FDPart3tmp1 = cos(xx1);
        const REAL FDPart3tmp2 = (1.0 / (xx0));
        const REAL FDPart3tmp3 = (1.0 / ((xx0) * (xx0)));
        const REAL FDPart3tmp10 = pow(fabs(cf), -3.0);
        const REAL FDPart3tmp12 = alpha * u4Ut;
        const REAL FDPart3tmp19 = hDD00 + 1;
        const REAL FDPart3tmp22 = hDD22 + 1;
        const REAL FDPart3tmp24 = hDD11 + 1;
        const REAL FDPart3tmp5 = (1.0 / (FDPart3tmp0));
        const REAL FDPart3tmp8 = FDPart3tmp1 / ((FDPart3tmp0) * (FDPart3tmp0));
        const REAL FDPart3tmp13 = FDPart3tmp12 * rescaledvU0 * rhob;
        const REAL FDPart3tmp15 = ((alpha) * (alpha)) * h * ((u4Ut) * (u4Ut));
        const REAL FDPart3tmp20 = FDPart3tmp19 * vetU0 + hDD01 * vetU1 + hDD02 * vetU2;
        const REAL FDPart3tmp26 = FDPart3tmp19 * FDPart3tmp22 * FDPart3tmp24 - FDPart3tmp19 * ((hDD12) * (hDD12)) -
                                  FDPart3tmp22 * ((hDD01) * (hDD01)) - FDPart3tmp24 * ((hDD02) * (hDD02)) + 2 * hDD01 * hDD02 * hDD12;
        const REAL FDPart3tmp30 = ((alpha) * (alpha)) * ((cf) * (cf));
        const REAL FDPart3tmp40 = FDPart3tmp24 * vetU1 + hDD01 * vetU0 + hDD12 * vetU2;
        const REAL FDPart3tmp41 = FDPart3tmp22 * vetU2 + hDD02 * vetU0 + hDD12 * vetU1;
        const REAL FDPart3tmp56 = FDPart3tmp10 * FDPart3tmp12 * S;
        const REAL FDPart3tmp7 = -FDPart3tmp3 * FDPart3tmp5;
        const REAL FDPart3tmp9 = -FDPart3tmp2 * FDPart3tmp8;
        const REAL FDPart3tmp16 = FDPart3tmp15 * rescaledvU0 * rhob;
        const REAL FDPart3tmp36 = FDPart3tmp15 * rhob;
        const REAL FDPart3tmp39 = FDPart3tmp10 / (FDPart3tmp26 * alpha * ((cf) * (cf)));
        const REAL FDPart3tmp43 = FDPart3tmp12 * rescaledvU1 * rhob;
        const REAL FDPart3tmp51 = FDPart3tmp12 * rescaledvU2 * rhob;
        const REAL FDPart3tmp17 = FDPart3tmp16 + P * vetU0;
        const REAL FDPart3tmp37 = FDPart3tmp26 * FDPart3tmp36;
        const REAL FDPart3tmp52 = FDPart3tmp36 * rescaledvU2 + P * vetU2;
        const REAL FDPart3tmp27 = FDPart3tmp17 * FDPart3tmp26;
        const REAL FDPart3tmp32 =
            FDPart3tmp16 * FDPart3tmp26 * rescaledvU1 + P * (-FDPart3tmp26 * vetU0 * vetU1 + FDPart3tmp30 * (-FDPart3tmp22 * hDD01 + hDD02 * hDD12));
        const REAL FDPart3tmp34 =
            FDPart3tmp16 * FDPart3tmp26 * rescaledvU2 + P * (-FDPart3tmp26 * vetU0 * vetU2 + FDPart3tmp30 * (-FDPart3tmp24 * hDD02 + hDD01 * hDD12));
        const REAL FDPart3tmp38 = FDPart3tmp37 * ((rescaledvU0) * (rescaledvU0)) +
                                  P * (-FDPart3tmp26 * ((vetU0) * (vetU0)) + FDPart3tmp30 * (FDPart3tmp22 * FDPart3tmp24 - ((hDD12) * (hDD12))));
        const REAL FDPart3tmp45 = FDPart3tmp36 * rescaledvU1 + P * vetU1;
        const REAL FDPart3tmp48 = FDPart3tmp26 * FDPart3tmp36 * rescaledvU1 * rescaledvU2 +
                                  P * (-FDPart3tmp26 * vetU1 * vetU2 + FDPart3tmp30 * (-FDPart3tmp19 * hDD12 + hDD01 * hDD02));
        const REAL FDPart3tmp49 = FDPart3tmp37 * ((rescaledvU1) * (rescaledvU1)) +
                                  P * (-FDPart3tmp26 * ((vetU1) * (vetU1)) + FDPart3tmp30 * (FDPart3tmp19 * FDPart3tmp22 - ((hDD02) * (hDD02))));
        const REAL FDPart3tmp54 = FDPart3tmp26 * FDPart3tmp52;
        const REAL FDPart3tmp55 = FDPart3tmp37 * ((rescaledvU2) * (rescaledvU2)) +
                                  P * (-FDPart3tmp26 * ((vetU2) * (vetU2)) + FDPart3tmp30 * (FDPart3tmp19 * FDPart3tmp24 - ((hDD01) * (hDD01))));
        const REAL FDPart3tmp47 = FDPart3tmp26 * FDPart3tmp45;
        ReU[0] = 1;
        ReUdD[0][0] = 0;
        ReUD[0][0] = 1;
        ReUDdD[0][0][0] = 0;
        ReUDdD[0][0][1] = 0;
        ReUDdD[0][0][2] = 0;
        ReUdD[0][1] = 0;
        ReUD[0][1] = xx0;
        ReUDdD[0][1][0] = 1;
        ReUDdD[0][1][1] = 0;
        ReUDdD[0][1][2] = 0;
        ReUdD[0][2] = 0;
        ReUD[0][2] = FDPart3tmp0 * xx0;
        ReUDdD[0][2][0] = FDPart3tmp0;
        ReUDdD[0][2][1] = FDPart3tmp1 * xx0;
        ReUDdD[0][2][2] = 0;
        ReU[1] = FDPart3tmp2;
        ReUdD[1][0] = -FDPart3tmp3;
        ReUD[1][0] = FDPart3tmp2;
        ReUDdD[1][0][0] = -FDPart3tmp3;
        ReUDdD[1][0][1] = 0;
        ReUDdD[1][0][2] = 0;
        ReUdD[1][1] = 0;
        ReUD[1][1] = 1;
        ReUDdD[1][1][0] = 0;
        ReUDdD[1][1][1] = 0;
        ReUDdD[1][1][2] = 0;
        ReUdD[1][2] = 0;
        ReUD[1][2] = FDPart3tmp0;
        ReUDdD[1][2][0] = 0;
        ReUDdD[1][2][1] = FDPart3tmp1;
        ReUDdD[1][2][2] = 0;
        ReU[2] = FDPart3tmp2 * FDPart3tmp5;
        ReUdD[2][0] = FDPart3tmp7;
        ReUD[2][0] = FDPart3tmp2 * FDPart3tmp5;
        ReUDdD[2][0][0] = FDPart3tmp7;
        ReUDdD[2][0][1] = FDPart3tmp9;
        ReUDdD[2][0][2] = 0;
        ReUdD[2][1] = FDPart3tmp9;
        ReUD[2][1] = FDPart3tmp5;
        ReUDdD[2][1][0] = 0;
        ReUDdD[2][1][1] = -FDPart3tmp8;
        ReUDdD[2][1][2] = 0;
        ReUdD[2][2] = 0;
        ReUD[2][2] = 1;
        ReUDdD[2][2][0] = 0;
        ReUDdD[2][2][1] = 0;
        ReUDdD[2][2][2] = 0;
        rescaledtau_tildefluxU[0] = FDPart3tmp10 * (-FDPart3tmp13 + FDPart3tmp17);
        rescaledrho_starfluxU[0] = FDPart3tmp10 * FDPart3tmp13;
        rescaledStildefluxUD[0][0] =
            FDPart3tmp39 * (FDPart3tmp19 * FDPart3tmp38 + FDPart3tmp20 * FDPart3tmp27 + FDPart3tmp32 * hDD01 + FDPart3tmp34 * hDD02);
        rescaledStildefluxUD[0][1] =
            FDPart3tmp39 * (FDPart3tmp24 * FDPart3tmp32 + FDPart3tmp27 * FDPart3tmp40 + FDPart3tmp34 * hDD12 + FDPart3tmp38 * hDD01);
        rescaledStildefluxUD[0][2] =
            FDPart3tmp39 * (FDPart3tmp22 * FDPart3tmp34 + FDPart3tmp27 * FDPart3tmp41 + FDPart3tmp32 * hDD12 + FDPart3tmp38 * hDD02);
        rescaledtau_tildefluxU[1] = FDPart3tmp10 * (-FDPart3tmp43 + FDPart3tmp45);
        rescaledrho_starfluxU[1] = FDPart3tmp10 * FDPart3tmp43;
        rescaledStildefluxUD[1][0] =
            FDPart3tmp39 * (FDPart3tmp19 * FDPart3tmp32 + FDPart3tmp20 * FDPart3tmp47 + FDPart3tmp48 * hDD02 + FDPart3tmp49 * hDD01);
        rescaledStildefluxUD[1][1] =
            FDPart3tmp39 * (FDPart3tmp24 * FDPart3tmp49 + FDPart3tmp32 * hDD01 + FDPart3tmp40 * FDPart3tmp47 + FDPart3tmp48 * hDD12);
        rescaledStildefluxUD[1][2] =
            FDPart3tmp39 * (FDPart3tmp22 * FDPart3tmp48 + FDPart3tmp32 * hDD02 + FDPart3tmp41 * FDPart3tmp47 + FDPart3tmp49 * hDD12);
        rescaledtau_tildefluxU[2] = FDPart3tmp10 * (-FDPart3tmp51 + FDPart3tmp52);
        rescaledrho_starfluxU[2] = FDPart3tmp10 * FDPart3tmp51;
        rescaledStildefluxUD[2][0] =
            FDPart3tmp39 * (FDPart3tmp19 * FDPart3tmp34 + FDPart3tmp20 * FDPart3tmp54 + FDPart3tmp48 * hDD01 + FDPart3tmp55 * hDD02);
        rescaledStildefluxUD[2][1] =
            FDPart3tmp39 * (FDPart3tmp24 * FDPart3tmp48 + FDPart3tmp34 * hDD01 + FDPart3tmp40 * FDPart3tmp54 + FDPart3tmp55 * hDD12);
        rescaledStildefluxUD[2][2] =
            FDPart3tmp39 * (FDPart3tmp22 * FDPart3tmp55 + FDPart3tmp34 * hDD02 + FDPart3tmp41 * FDPart3tmp54 + FDPart3tmp48 * hDD12);
        rescaledYe_starfluxU[0] = FDPart3tmp10 * FDPart3tmp13 * Ye;
        rescaledYe_starfluxU[1] = FDPart3tmp10 * FDPart3tmp43 * Ye;
        rescaledYe_starfluxU[2] = FDPart3tmp10 * FDPart3tmp51 * Ye;
        rescaledS_starfluxU[0] = FDPart3tmp56 * rescaledvU0;
        rescaledS_starfluxU[1] = FDPart3tmp56 * rescaledvU1;
        rescaledS_starfluxU[2] = FDPart3tmp56 * rescaledvU2;

        rhs_gfs[IDX4pt(RHO_STARGF, index)] +=
            -rescaledrho_starfluxU[flux_dirn] * ReUdD[flux_dirn][flux_dirn] +
            ReU[flux_dirn] * (auxevol_gfs[IDX4pt(RHO_STAR_HLL_FLUXGF, index)] - auxevol_gfs[IDX4pt(RHO_STAR_HLL_FLUXGF, indexp1)]) * invdx;

        rhs_gfs[IDX4pt(TAU_TILDEGF, index)] +=
            -rescaledtau_tildefluxU[flux_dirn] * ReUdD[flux_dirn][flux_dirn] +
            ReU[flux_dirn] * (auxevol_gfs[IDX4pt(TAU_TILDE_HLL_FLUXGF, index)] - auxevol_gfs[IDX4pt(TAU_TILDE_HLL_FLUXGF, indexp1)]) * invdx;

        REAL StildeD0_rhs;
        REAL StildeD1_rhs;
        REAL StildeD2_rhs;

        StildeD0_rhs = -rescaledStildefluxUD[flux_dirn][0] * ReUDdD[flux_dirn][0][flux_dirn] +
                       ReUD[flux_dirn][0] *
                           (auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD0GF, index)] - auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD0GF, indexp1)]) *
                           invdx;
        StildeD1_rhs = -rescaledStildefluxUD[flux_dirn][1] * ReUDdD[flux_dirn][1][flux_dirn] +
                       ReUD[flux_dirn][1] *
                           (auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD1GF, index)] - auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD1GF, indexp1)]) *
                           invdx;
        StildeD2_rhs = -rescaledStildefluxUD[flux_dirn][2] * ReUDdD[flux_dirn][2][flux_dirn] +
                       ReUD[flux_dirn][2] *
                           (auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD2GF, index)] - auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD2GF, indexp1)]) *
                           invdx;

        // Rescale the momentum RHSs back to the evolved variables.
        rhs_gfs[IDX4pt(RESCALEDSTILDED0GF, index)] += StildeD0_rhs * ReU[0];
        rhs_gfs[IDX4pt(RESCALEDSTILDED1GF, index)] += StildeD1_rhs * ReU[1];
        rhs_gfs[IDX4pt(RESCALEDSTILDED2GF, index)] += StildeD2_rhs * ReU[2];

        rhs_gfs[IDX4pt(YE_STARGF, index)] +=
            -rescaledYe_starfluxU[flux_dirn] * ReUdD[flux_dirn][flux_dirn] +
            ReU[flux_dirn] * (auxevol_gfs[IDX4pt(YE_STAR_HLL_FLUXGF, index)] - auxevol_gfs[IDX4pt(YE_STAR_HLL_FLUXGF, indexp1)]) * invdx;

        rhs_gfs[IDX4pt(S_STARGF, index)] +=
            -rescaledS_starfluxU[flux_dirn] * ReUdD[flux_dirn][flux_dirn] +
            ReU[flux_dirn] * (auxevol_gfs[IDX4pt(S_STAR_HLL_FLUXGF, index)] - auxevol_gfs[IDX4pt(S_STAR_HLL_FLUXGF, indexp1)]) * invdx;

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION calculate_flux_divergences__rfm__Spherical
