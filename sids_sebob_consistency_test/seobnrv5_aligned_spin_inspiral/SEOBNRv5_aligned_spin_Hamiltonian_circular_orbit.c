#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

/**
 * Evaluate SEOBNRv5 Hamiltonian's circular orbit conditions to compute conservative initial conditions.
 */
int SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit(const gsl_vector *restrict x, void *restrict params, gsl_vector *restrict f) {

  const REAL r = gsl_vector_get(x, 0);
  const REAL pphi = gsl_vector_get(x, 1);
  const REAL m1 = ((commondata_struct *restrict)params)->m1;
  const REAL m2 = ((commondata_struct *restrict)params)->m2;
  const REAL chi1 = ((commondata_struct *restrict)params)->chi1;
  const REAL chi2 = ((commondata_struct *restrict)params)->chi2;
  const REAL a6 = ((commondata_struct *restrict)params)->a6;
  const REAL dSO = ((commondata_struct *restrict)params)->dSO;
  const REAL initial_omega = ((commondata_struct *restrict)params)->initial_omega;
  const REAL tmp0 = ((r) * (r) * (r));
  const REAL tmp1 = m1 + m2;
  const REAL tmp16 = ((r) * (r));
  const REAL tmp23 = (1.0 / (r));
  const REAL tmp27 = ((pphi) * (pphi));
  const REAL tmp29 = ((m1) * (m1));
  const REAL tmp30 = ((m2) * (m2));
  const REAL tmp36 = ((r) * (r) * (r) * (r));
  const REAL tmp38 = ((pphi) * (pphi) * (pphi) * (pphi));
  const REAL tmp53 = ((r) * (r) * (r) * (r) * (r));
  const REAL tmp56 = log(r);
  const REAL tmp2 = (1.0 / ((tmp1) * (tmp1)));
  const REAL tmp5 = chi1 * m1 + chi2 * m2;
  const REAL tmp10 = (1.0 / (tmp0));
  const REAL tmp13 = (1.0 / ((tmp1) * (tmp1) * (tmp1)));
  const REAL tmp17 = (1.0 / (tmp16));
  const REAL tmp18 = (1.0 / ((tmp1) * (tmp1) * (tmp1) * (tmp1)));
  const REAL tmp19 = chi1 * m1 - chi2 * m2;
  const REAL tmp37 = (1.0 / (tmp36));
  const REAL tmp54 = (1.0 / (tmp53));
  const REAL tmp57 = 1120 * ((m1) * (m1) * (m1)) * ((m2) * (m2) * (m2)) / pow(tmp1, 6);
  const REAL tmp60 = 8 * r + 2 * tmp0 + 4 * tmp16 + 16;
  const REAL tmp78 = 7680 * tmp36;
  const REAL tmp80 = 2 * tmp23 + 1;
  const REAL tmp103 = 8 * r + 6 * tmp16 + 8;
  const REAL tmp6 = ((tmp5) * (tmp5));
  const REAL tmp20 = tmp19 * (m1 - m2);
  const REAL tmp24 = m1 * m2 * tmp2;
  const REAL tmp28 = tmp17 * tmp27;
  const REAL tmp32 = tmp18 * tmp29 * tmp30;
  const REAL tmp34 = tmp10 * tmp27;
  const REAL tmp40 = tmp5 / tmp1;
  const REAL tmp50 = ((tmp19) * (tmp19)) * tmp2;
  const REAL tmp94 = 2 * tmp10;
  const REAL tmp106 = 2 * pphi * tmp17;
  const REAL tmp108 = 4 * ((pphi) * (pphi) * (pphi)) * tmp37;
  const REAL tmp7 = tmp2 * tmp6;
  const REAL tmp15 = dSO * m1 * m2 * tmp13 * tmp5;
  const REAL tmp21 = -tmp13 * ((tmp5) * (tmp5) * (tmp5)) + tmp18 * tmp20 * tmp6;
  const REAL tmp25 = (23.0 / 32.0) * tmp24 - 3.0 / 32.0;
  const REAL tmp26 = -45.0 / 32.0 * tmp24 - 15.0 / 32.0;
  const REAL tmp31 = (109.0 / 192.0) * tmp18 * tmp29 * tmp30 - 177.0 / 32.0 * tmp24 - 5.0 / 64.0;
  const REAL tmp33 = -267.0 / 128.0 * tmp24 - 1591.0 / 768.0 * tmp32 + 59.0 / 256.0;
  const REAL tmp35 = (75.0 / 128.0) * tmp24 + (345.0 / 256.0) * tmp32 + 105.0 / 256.0;
  const REAL tmp41 = (11.0 / 32.0) * tmp24 + 3.0 / 32.0;
  const REAL tmp42 = 15.0 / 32.0 - 9.0 / 32.0 * tmp24;
  const REAL tmp43 = -1.0 / 32.0 * tmp24 + (103.0 / 192.0) * tmp32 + 5.0 / 64.0;
  const REAL tmp44 = -35.0 / 128.0 * tmp24 - 613.0 / 768.0 * tmp32 - 59.0 / 256.0;
  const REAL tmp45 = (75.0 / 256.0) * tmp18 * tmp29 * tmp30 - 45.0 / 128.0 * tmp24 - 105.0 / 256.0;
  const REAL tmp46 = tmp2 * tmp20;
  const REAL tmp51 = tmp13 * tmp20 * tmp5;
  const REAL tmp67 = tmp24 * tmp56;
  const REAL tmp68 = 588 * tmp24 + 1079;
  const REAL tmp93 = 2 * tmp34;
  const REAL tmp96 = 3 * tmp27 * tmp37;
  const REAL tmp98 = 4 * tmp38 * tmp54;
  const REAL tmp8 = tmp0 + tmp7 * (r + 2);
  const REAL tmp22 = (1.0 / 4.0) * tmp17 * tmp21;
  const REAL tmp47 = tmp40 * (tmp17 * tmp31 + tmp23 * tmp25 + tmp26 * tmp28 + tmp33 * tmp34 + tmp35 * tmp37 * tmp38 + 7.0 / 4.0) +
                     tmp46 * (tmp17 * tmp43 + tmp23 * tmp41 + tmp28 * tmp42 + tmp34 * tmp44 + tmp37 * tmp38 * tmp45 + 1.0 / 4.0);
  const REAL tmp52 = tmp50 * ((1.0 / 2.0) * tmp24 + 1.0 / 8.0) - 5.0 / 4.0 * tmp51 + (9.0 / 8.0) * tmp7;
  const REAL tmp55 = tmp50 * ((21.0 / 16.0) * tmp18 * tmp29 * tmp30 - 81.0 / 64.0 * tmp24 - 9.0 / 64.0) +
                     tmp51 * (117.0 / 32.0 - 39.0 / 16.0 * tmp24) + tmp7 * (-175.0 / 64.0 * tmp24 - 225.0 / 64.0);
  const REAL tmp70 = 756 * tmp24 + 1079;
  const REAL tmp75 = 336 * r + 756 * tmp24 + 407;
  const REAL tmp81 = tmp7 * tmp80;
  const REAL tmp9 = (1.0 / (tmp8));
  const REAL tmp48 = pphi * tmp10 * tmp15 + pphi * tmp22 + pphi * tmp47;
  const REAL tmp76 = 2048 * m1 * m2 * tmp2 * tmp56 * tmp75 + 28 * m1 * m2 * tmp2 * (1920 * a6 + 733955307463037.0 / 1000000000.0) -
                     7 * r * ((938918400156317.0 / 1000000000.0) * m1 * m2 * tmp2 - 185763092693281.0 / 1000000000.0 * tmp32 - 245760) -
                     270820329770593.0 / 50000000.0 * tmp32 - 3440640;
  const REAL tmp82 = tmp17 * tmp81 + 1;
  const REAL tmp72 = 7864320 * r * (-373313530533103.0 / 50000000000.0 * tmp24 - 3024 * tmp32 + 17264) + 31457280 * tmp0 * tmp70 +
                     62914560 * tmp16 * tmp68 - 99436138180993792.0 / 78125.0 * tmp24 - 24095439828988416.0 / 390625.0 * tmp32 +
                     15728640 * tmp36 * tmp70 + 5284823040 * tmp53 + 440653578240;
  const REAL tmp83 = (1.0 / (tmp82));
  const REAL tmp85 = (1.0 / (tmp16 + tmp81));
  const REAL tmp73 =
      (241555486248807.0 / 1000.0) * ((m1) * (m1) * (m1) * (m1)) * ((m2) * (m2) * (m2) * (m2)) / pow(tmp1, 8) +
      53760 * tmp24 *
          (7680 * a6 * (tmp36 + tmp60) + (113485217444961.0 / 1000000000.0) * r * (-tmp36 + tmp60) +
           (7402203300817.0 / 50000000000.0) * r * (7704 * r + 1926 * tmp0 + 3852 * tmp16 + 349 * tmp36 + 36400) +
           128 * r *
               ((33046962773603.0 / 2500000000.0) * r + (42646962773603.0 / 10000000000.0) * tmp0 + (852939255472061.0 / 100000000000.0) * tmp16 -
                137046962773603.0 / 20000000000.0 * tmp36 - 42153037226397.0 / 1250000000.0)) +
      67645734912 * tmp32 * ((tmp56) * (tmp56)) +
      7 * tmp32 *
          (-39321600 * a6 * (3 * r + 59) + (186464462028901.0 / 250000.0) * a6 + (122635399361987.0 / 1000.0) * r -
           308925070376879.0 / 50000.0 * tmp0 + (206478381132587.0 / 100000.0) * tmp16 - 308925070376879.0 / 100000.0 * tmp36 +
           (3566651379711.0 / 2500.0) * tmp53 + 276057889687011.0 / 1000.0) +
      13212057600 * tmp53 + tmp57 * (-163683964822551.0 / 1000000.0 * r - 3566651379711.0 / 200000.0 * tmp16 - 59449372951581.0 / 50000.0) +
      tmp67 * tmp72;
  const REAL tmp87 = -tmp80 * tmp85;
  const REAL tmp74 = (1.0 / (tmp73));
  const REAL tmp89 = tmp28 * tmp7 * tmp87 + tmp28 + 1;
  const REAL tmp79 = tmp17 * tmp7 + tmp37 * tmp52 + tmp54 * tmp55 + tmp74 * tmp76 * tmp78;
  const REAL tmp91 = sqrt(tmp79 * tmp83 * tmp89);
  const REAL tmp92 = (1.0 / sqrt(2 * tmp24 * (tmp48 * tmp9 + tmp91 - 1) + 1));
  const REAL tmp104 = tmp91 / tmp89;
  const REAL dHreal_dr_circ =
      tmp92 *
      (tmp104 * tmp82 *
           ((1.0 / 2.0) * tmp79 * tmp83 *
                (2 * tmp27 * tmp37 * tmp7 * tmp85 - tmp28 * tmp7 * tmp80 * (-2 * r + 2 * tmp17 * tmp2 * tmp6) / ((tmp16 + tmp81) * (tmp16 + tmp81)) -
                 tmp7 * tmp87 * tmp93 - tmp93) +
            (1.0 / 2.0) * tmp79 * tmp89 * (2 * tmp37 * tmp7 + tmp81 * tmp94) / ((tmp82) * (tmp82)) +
            (1.0 / 2.0) * tmp83 * tmp89 *
                (30720 * tmp0 * tmp74 * tmp76 - 4 * tmp52 * tmp54 - tmp7 * tmp94 +
                 tmp74 * tmp78 *
                     (2048 * tmp23 * tmp24 * tmp75 - 6572428801094219.0 / 1000000000.0 * tmp24 + (1300341648852967.0 / 1000000000.0) * tmp32 +
                      688128 * tmp67 + 1720320) +
                 tmp76 * tmp78 *
                     (-tmp23 * tmp24 * tmp72 - 135291469824 * tmp23 * tmp32 * tmp56 -
                      53760 * tmp24 *
                          (7680 * a6 * (4 * tmp0 + tmp103) + (113485217444961.0 / 1000000000.0) * r * (-4 * tmp0 + tmp103) +
                           (7402203300817.0 / 50000000000.0) * r * (7704 * r + 1396 * tmp0 + 5778 * tmp16 + 7704) +
                           128 * r *
                               ((852939255472061.0 / 50000000000.0) * r - 137046962773603.0 / 5000000000.0 * tmp0 +
                                (127940888320809.0 / 10000000000.0) * tmp16 + 33046962773603.0 / 2500000000.0) +
                           (23377610738487781.0 / 6250000000.0) * r + (26449610738487781.0 / 25000000000.0) * tmp0 +
                           (26449610738487797.0 / 12500000000.0) * tmp16 - 46945920007815877.0 / 50000000000.0 * tmp36 +
                           1805060293831937.0 / 625000000.0) -
                      7 * tmp32 *
                          (-117964800 * a6 + (206478381132587.0 / 50000.0) * r - 308925070376879.0 / 25000.0 * tmp0 -
                           926775211130637.0 / 50000.0 * tmp16 + (3566651379711.0 / 500.0) * tmp36 + 122635399361987.0 / 1000.0) -
                      66060288000 * tmp36 - tmp57 * (-3566651379711.0 / 100000.0 * r - 163683964822551.0 / 1000000.0) -
                      32768 * tmp67 *
                          (3840 * r * tmp68 + 1920 * tmp0 * tmp70 + 2880 * tmp16 * tmp70 - 1119940591599309.0 / 625000000.0 * tmp24 - 725760 * tmp32 +
                           806400 * tmp36 + 4143360)) /
                     ((tmp73) * (tmp73)) -
                 5 * tmp55 / pow(r, 6))) /
           tmp79 +
       tmp48 * (-3 * tmp16 - tmp7) / ((tmp8) * (tmp8)) +
       tmp9 * (-1.0 / 2.0 * pphi * tmp10 * tmp21 - 3 * pphi * tmp15 * tmp37 +
               pphi * (tmp40 * (-tmp17 * tmp25 - tmp26 * tmp93 - tmp31 * tmp94 - tmp33 * tmp96 - tmp35 * tmp98) +
                       tmp46 * (-tmp17 * tmp41 - tmp42 * tmp93 - tmp43 * tmp94 - tmp44 * tmp96 - tmp45 * tmp98))));
  const REAL dHreal_dpphi_circ = tmp92 * ((1.0 / 2.0) * tmp104 * (2 * pphi * tmp17 * tmp7 * tmp87 + tmp106) +
                                          tmp9 * (pphi * (tmp40 * (pphi * tmp33 * tmp94 + tmp106 * tmp26 + tmp108 * tmp35) +
                                                          tmp46 * (pphi * tmp44 * tmp94 + tmp106 * tmp42 + tmp108 * tmp45)) +
                                                  tmp10 * tmp15 + tmp22 + tmp47));

  gsl_vector_set(f, 0, dHreal_dr_circ);
  gsl_vector_set(f, 1, dHreal_dpphi_circ - initial_omega);
  return GSL_SUCCESS;
} // END FUNCTION SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit
