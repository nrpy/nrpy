#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Evaluate SEOBNRv5 Hamiltonian and needed derivatives to compute binary dynamics.
 */
int SEOBNRv5_aligned_spin_right_hand_sides(REAL t, const REAL *restrict y, REAL *restrict f, void *restrict params) {

  const REAL m1 = ((commondata_struct *restrict)params)->m1;
  const REAL m2 = ((commondata_struct *restrict)params)->m2;
  const REAL chi1 = ((commondata_struct *restrict)params)->chi1;
  const REAL chi2 = ((commondata_struct *restrict)params)->chi2;
  const REAL a6 = ((commondata_struct *restrict)params)->a6;
  const REAL dSO = ((commondata_struct *restrict)params)->dSO;
  const REAL r = y[0];
  const REAL prstar = y[2];
  const REAL pphi = y[3];
  const REAL tmp0 = ((m1) * (m1) * (m1));
  const REAL tmp1 = ((m2) * (m2) * (m2));
  const REAL tmp2 = m1 + m2;
  const REAL tmp5 = log(r);
  const REAL tmp9 = ((m1) * (m1));
  const REAL tmp10 = ((m2) * (m2));
  const REAL tmp14 = ((r) * (r) * (r) * (r));
  const REAL tmp16 = ((prstar) * (prstar) * (prstar) * (prstar));
  const REAL tmp19 = ((r) * (r) * (r));
  const REAL tmp22 = ((r) * (r));
  const REAL tmp26 = pow(prstar, 6);
  const REAL tmp32 = pow(prstar, 8);
  const REAL tmp35 = ((prstar) * (prstar));
  const REAL tmp53 = (1.0 / (r));
  const REAL tmp62 = ((r) * (r) * (r) * (r) * (r));
  const REAL tmp101 = pow(r, -5.0 / 2.0);
  const REAL tmp103 = pow(r, -9.0 / 2.0);
  const REAL tmp113 = ((pphi) * (pphi));
  const REAL tmp144 = ((pphi) * (pphi) * (pphi) * (pphi));
  const REAL tmp168 = 5 / pow(r, 6);
  const REAL tmp192 = ((prstar) * (prstar) * (prstar));
  const REAL tmp193 = ((prstar) * (prstar) * (prstar) * (prstar) * (prstar));
  const REAL tmp194 = pow(prstar, 7);
  const REAL tmp3 = pow(tmp2, -6);
  const REAL tmp6 = (1.0 / ((tmp2) * (tmp2)));
  const REAL tmp11 = (1.0 / ((tmp2) * (tmp2) * (tmp2) * (tmp2)));
  const REAL tmp15 = (1.0 / (tmp14));
  const REAL tmp20 = (1.0 / (tmp19));
  const REAL tmp23 = (1.0 / (tmp22));
  const REAL tmp27 = ((m1) * (m1) * (m1) * (m1)) * ((m2) * (m2) * (m2) * (m2)) / pow(tmp2, 8);
  const REAL tmp38 = chi1 * m1 + chi2 * m2;
  const REAL tmp44 = ((tmp5) * (tmp5));
  const REAL tmp66 = 8 * r + 2 * tmp19 + 4 * tmp22 + 16;
  const REAL tmp84 = 7680 * tmp14;
  const REAL tmp87 = chi1 * m1 - chi2 * m2;
  const REAL tmp90 = (1.0 / ((tmp2) * (tmp2) * (tmp2)));
  const REAL tmp98 = 7680 * tmp62;
  const REAL tmp106 = (1.0 / (tmp62));
  const REAL tmp115 = 2 * tmp53;
  const REAL tmp173 = 8 * r + 6 * tmp22 + 8;
  const REAL tmp195 = 4 * tmp192;
  const REAL tmp4 = tmp0 * tmp1 * tmp3;
  const REAL tmp8 = m1 * m2 * tmp6;
  const REAL tmp12 = tmp10 * tmp11 * tmp9;
  const REAL tmp39 = ((tmp38) * (tmp38));
  const REAL tmp45 = 5787938193408 * tmp44;
  const REAL tmp88 = tmp6 * ((tmp87) * (tmp87));
  const REAL tmp89 = tmp87 * (m1 - m2);
  const REAL tmp114 = tmp113 * tmp23;
  const REAL tmp116 = tmp115 + 1;
  const REAL tmp142 = tmp113 * tmp20;
  const REAL tmp146 = tmp38 / tmp2;
  const REAL tmp159 = 3 * tmp15;
  const REAL tmp161 = 2 * tmp20;
  const REAL tmp163 = 4 * tmp106;
  const REAL tmp198 = 2 * pphi * tmp23;
  const REAL tmp200 = 4 * ((pphi) * (pphi) * (pphi)) * tmp15;
  const REAL tmp13 = tmp12 * (118.40000000000001 * tmp5 - 179613660498019.0 / 100000000000.0) + (150579635104141.0 / 250000000000.0) * tmp4 +
                     tmp8 * (452542166996693.0 / 1000000000000.0 - 516952380952381.0 / 10000000000000.0 * tmp5);
  const REAL tmp18 = -131 * tmp12 + 10 * tmp4 + (115888805356193.0 / 1250000000000.0) * tmp8;
  const REAL tmp24 = 8 * m1 * m2 * tmp6 - 6 * tmp12;
  const REAL tmp28 = 188 * tmp0 * tmp1 * tmp3 - 447649163680617.0 / 5000000000000.0 * tmp12 - 14 * tmp27 - 84945530542609.0 / 2500000000000.0 * tmp8;
  const REAL tmp30 = 6 * tmp0 * tmp1 * tmp3 - 27.0 / 5.0 * tmp12 - 139150381847503.0 / 50000000000000.0 * tmp8;
  const REAL tmp33 = (166921011824161.0 / 50000000000000.0) * tmp12 - 6 * tmp27 + (342857142857143.0 / 100000000000000.0) * tmp4 +
                     (4343054718629.0 / 3125000000000.0) * tmp8;
  const REAL tmp40 = tmp39 * tmp6;
  const REAL tmp46 = -98886748396767.0 / 500000.0 * tmp12 - 7876452485916241.0 / 12500000.0 * tmp8 + 580530436787913.0 / 100000.0;
  const REAL tmp47 = -138240 * tmp12 - 53511513369581.0 / 20000000.0 * tmp8 - 52783413229329.0 / 10000000.0;
  const REAL tmp48 = ((1 - 0.49694878161693501 * tmp8) * (1 - 0.49694878161693501 * tmp8));
  const REAL tmp50 = 14700 * tmp8 + 42911;
  const REAL tmp54 = r * tmp8;
  const REAL tmp73 = tmp5 * tmp8;
  const REAL tmp74 = 588 * tmp8 + 1079;
  const REAL tmp92 = tmp38 * tmp89 * tmp90;
  const REAL tmp100 = tmp8 / pow(r, 7.0 / 2.0);
  const REAL tmp105 = tmp53 * tmp8;
  const REAL tmp135 = dSO * m1 * m2 * tmp38 * tmp90;
  const REAL tmp136 = tmp11 * tmp39 * tmp89 - ((tmp38) * (tmp38) * (tmp38)) * tmp90;
  const REAL tmp138 = (23.0 / 32.0) * tmp8 - 3.0 / 32.0;
  const REAL tmp139 = -45.0 / 32.0 * tmp8 - 15.0 / 32.0;
  const REAL tmp140 = (109.0 / 192.0) * tmp10 * tmp11 * tmp9 - 177.0 / 32.0 * tmp8 - 5.0 / 64.0;
  const REAL tmp141 = -1591.0 / 768.0 * tmp12 - 267.0 / 128.0 * tmp8 + 59.0 / 256.0;
  const REAL tmp143 = (345.0 / 256.0) * tmp12 + (75.0 / 128.0) * tmp8 + 105.0 / 256.0;
  const REAL tmp147 = (11.0 / 32.0) * tmp8 + 3.0 / 32.0;
  const REAL tmp148 = 15.0 / 32.0 - 9.0 / 32.0 * tmp8;
  const REAL tmp149 = (103.0 / 192.0) * tmp12 - 1.0 / 32.0 * tmp8 + 5.0 / 64.0;
  const REAL tmp150 = -613.0 / 768.0 * tmp12 - 35.0 / 128.0 * tmp8 - 59.0 / 256.0;
  const REAL tmp151 = (75.0 / 256.0) * tmp10 * tmp11 * tmp9 - 45.0 / 128.0 * tmp8 - 105.0 / 256.0;
  const REAL tmp152 = tmp6 * tmp89;
  const REAL tmp160 = 2 * tmp142;
  const REAL tmp178 = -1759846072320 * r - 4161798144000 * tmp8 - 12148770078720;
  const REAL tmp41 = tmp23 * tmp40;
  const REAL tmp49 = (3369809522332764779.0 / 78125.0) * tmp12 + (296393260946151.0 / 50.0) * tmp4 + (866182644304933.0 / 10.0) * tmp48 +
                     (3390101660860523312.0 / 78125.0) * tmp8 + 188440788778196;
  const REAL tmp51 = 49152 * r * (102574080 * tmp12 + (16368307925443.0 / 40000.0) * tmp8 - 105983591868019.0 / 50000.0) + 879923036160 * tmp19 +
                     283115520 * tmp22 * tmp50 - 1698693120 * tmp8 * (11592 * tmp8 + 69847);
  const REAL tmp55 = r * tmp50;
  const REAL tmp76 = 756 * tmp8 + 1079;
  const REAL tmp81 = 336 * r + 756 * tmp8 + 407;
  const REAL tmp93 = tmp40 * (3 * tmp8 + 45.0 / 16.0) + tmp88 * ((3.0 / 4.0) * tmp8 - 3.0 / 16.0) - 21.0 / 8.0 * tmp92;
  const REAL tmp94 = tmp40 * (-1171.0 / 64.0 * tmp8 - 861.0 / 64.0) + tmp88 * ((1.0 / 16.0) * tmp12 + (115.0 / 64.0) * tmp8 - 37.0 / 64.0) +
                     tmp92 * ((13.0 / 16.0) * tmp8 + 449.0 / 32.0);
  const REAL tmp109 = tmp40 * (-5 * tmp12 + (165.0 / 32.0) * tmp8 + 25.0 / 32.0) +
                      tmp88 * ((75.0 / 32.0) * m1 * m2 * tmp6 - 15.0 / 8.0 * tmp12 - 15.0 / 32.0) + tmp92 * ((45.0 / 8.0) * tmp8 - 5.0 / 16.0);
  const REAL tmp127 = tmp40 * (-175.0 / 64.0 * tmp8 - 225.0 / 64.0) +
                      tmp88 * ((21.0 / 16.0) * tmp10 * tmp11 * tmp9 - 81.0 / 64.0 * tmp8 - 9.0 / 64.0) + tmp92 * (117.0 / 32.0 - 39.0 / 16.0 * tmp8);
  const REAL tmp128 = (9.0 / 8.0) * tmp40 + tmp88 * ((1.0 / 2.0) * tmp8 + 1.0 / 8.0) - 5.0 / 4.0 * tmp92;
  const REAL tmp132 = tmp19 + tmp40 * (r + 2);
  const REAL tmp137 = (1.0 / 4.0) * tmp136 * tmp23;
  const REAL tmp153 = tmp146 * (tmp114 * tmp139 + tmp138 * tmp53 + tmp140 * tmp23 + tmp141 * tmp142 + tmp143 * tmp144 * tmp15 + 7.0 / 4.0) +
                      tmp152 * (tmp114 * tmp148 + tmp142 * tmp150 + tmp144 * tmp15 * tmp151 + tmp147 * tmp53 + tmp149 * tmp23 + 1.0 / 4.0);
  const REAL tmp167 = -tmp161 * tmp40;
  const REAL tmp42 = tmp41 + 1;
  const REAL tmp52 =
      r * tmp45 + r * tmp49 - 967680 * tmp19 * tmp47 - 9216 * tmp22 * tmp46 + tmp5 * tmp51 +
      55296 * tmp8 *
          (-426364516032331.0 / 10000000.0 * tmp12 + 14515200 * tmp4 - 31383302728516087.0 / 12500000.0 * tmp8 + 100201376401019.0 / 100000.0);
  const REAL tmp58 = -5041721180160 * tmp12 - 879923036160 * tmp22 - 283115520 * tmp55 - 253929149957443.0 / 10.0 * tmp8 + 104186110149937;
  const REAL tmp82 = 2048 * m1 * m2 * tmp5 * tmp6 * tmp81 + 28 * m1 * m2 * tmp6 * (1920 * a6 + 733955307463037.0 / 1000000000.0) -
                     7 * r * ((938918400156317.0 / 1000000000.0) * m1 * m2 * tmp6 - 185763092693281.0 / 1000000000.0 * tmp12 - 245760) -
                     270820329770593.0 / 50000000.0 * tmp12 - 3440640;
  const REAL tmp118 = tmp116 * tmp40 + tmp22;
  const REAL tmp125 = tmp116 * tmp41 + 1;
  const REAL tmp133 = (1.0 / (tmp132));
  const REAL tmp154 = pphi * tmp135 * tmp20 + pphi * tmp137 + pphi * tmp153;
  const REAL tmp181 = tmp5 * (5041721180160 * tmp12 + 2639769108480 * tmp22 + 566231040 * tmp55 + (12570860486740224.0 / 625.0) * tmp8 -
                              325581594218554368.0 / 3125.0);
  const REAL tmp43 = ((tmp42) * (tmp42));
  const REAL tmp59 = (182268054644921.0 / 100.0) * r * tmp12 - 197383821284628.0 / 5.0 * r + 133772083200 * tmp12 * tmp22 +
                     (400296247701391.0 / 5.0) * tmp12 + (258910106287381.0 / 100.0) * tmp22 * tmp8 + (510774533137571.0 / 100.0) * tmp22 +
                     (336524885906151.0 / 50.0) * tmp4 + tmp45 - tmp5 * tmp58 + (53058305272831.0 / 5.0) * tmp54 -
                     163418713120743.0 / 500000.0 * tmp55 - 60249543508726.0 / 5.0 * tmp8 + 275059053208689;
  const REAL tmp78 = 7864320 * r * (-3024 * tmp12 - 373313530533103.0 / 50000000000.0 * tmp8 + 17264) - 24095439828988416.0 / 390625.0 * tmp12 +
                     15728640 * tmp14 * tmp76 + 31457280 * tmp19 * tmp76 + 62914560 * tmp22 * tmp74 + 5284823040 * tmp62 -
                     99436138180993792.0 / 78125.0 * tmp8 + 440653578240;
  const REAL tmp95 = (1.0 / (tmp52));
  const REAL tmp107 = ((tmp42) * (tmp42) * (tmp42) * (tmp42));
  const REAL tmp111 = ((tmp52) * (tmp52));
  const REAL tmp119 = (1.0 / (tmp118));
  const REAL tmp126 = (1.0 / (tmp125));
  const REAL tmp155 = tmp133 * tmp154 - 1;
  const REAL tmp186 = -18432 * r * tmp46 + tmp181 - 2903040 * tmp22 * tmp47 + tmp45 + tmp49 + 11575876386816 * tmp5 + tmp51 * tmp53;
  const REAL tmp187 = 267544166400 * r * tmp12 + (510774533137571.0 / 50.0) * r + (182268054644921.0 / 100.0) * tmp12 - tmp178 * tmp5 +
                      11575876386816 * tmp5 * tmp53 - tmp53 * tmp58 + (258910106287381.0 / 50.0) * tmp54 + (29035754444081779.0 / 5000.0) * tmp8 -
                      26750842527187002873.0 / 500000.0;
  const REAL tmp201 = tmp133 * (pphi * (tmp146 * (pphi * tmp141 * tmp161 + tmp139 * tmp198 + tmp143 * tmp200) +
                                        tmp152 * (pphi * tmp150 * tmp161 + tmp148 * tmp198 + tmp151 * tmp200)) +
                                tmp135 * tmp20 + tmp137 + tmp153);
  const REAL tmp60 = (1.0 / (tmp59));
  const REAL tmp79 =
      67645734912 * tmp12 * tmp44 +
      7 * tmp12 *
          (-39321600 * a6 * (3 * r + 59) + (186464462028901.0 / 250000.0) * a6 + (122635399361987.0 / 1000.0) * r -
           308925070376879.0 / 100000.0 * tmp14 - 308925070376879.0 / 50000.0 * tmp19 + (206478381132587.0 / 100000.0) * tmp22 +
           (3566651379711.0 / 2500.0) * tmp62 + 276057889687011.0 / 1000.0) +
      (241555486248807.0 / 1000.0) * tmp27 +
      1120 * tmp4 * (-163683964822551.0 / 1000000.0 * r - 3566651379711.0 / 200000.0 * tmp22 - 59449372951581.0 / 50000.0) + 13212057600 * tmp62 +
      tmp73 * tmp78 +
      53760 * tmp8 *
          (7680 * a6 * (tmp14 + tmp66) + (113485217444961.0 / 1000000000.0) * r * (-tmp14 + tmp66) +
           (7402203300817.0 / 50000000000.0) * r * (7704 * r + 349 * tmp14 + 1926 * tmp19 + 3852 * tmp22 + 36400) +
           128 * r *
               ((33046962773603.0 / 2500000000.0) * r - 137046962773603.0 / 20000000000.0 * tmp14 + (42646962773603.0 / 10000000000.0) * tmp19 +
                (852939255472061.0 / 100000000000.0) * tmp22 - 42153037226397.0 / 1250000000.0));
  const REAL tmp96 = tmp59 * tmp95;
  const REAL tmp110 = (1.0 / ((tmp59) * (tmp59)));
  const REAL tmp121 = -tmp116 * tmp119;
  const REAL tmp80 = (1.0 / (tmp79));
  const REAL tmp112 = tmp110 * tmp111;
  const REAL tmp174 =
      tmp82 *
      (-tmp105 * tmp78 - 135291469824 * tmp12 * tmp5 * tmp53 -
       7 * tmp12 *
           (-117964800 * a6 + (206478381132587.0 / 50000.0) * r + (3566651379711.0 / 500.0) * tmp14 - 308925070376879.0 / 25000.0 * tmp19 -
            926775211130637.0 / 50000.0 * tmp22 + 122635399361987.0 / 1000.0) -
       66060288000 * tmp14 - 1120 * tmp4 * (-3566651379711.0 / 100000.0 * r - 163683964822551.0 / 1000000.0) -
       32768 * tmp73 *
           (3840 * r * tmp74 - 725760 * tmp12 + 806400 * tmp14 + 1920 * tmp19 * tmp76 + 2880 * tmp22 * tmp76 -
            1119940591599309.0 / 625000000.0 * tmp8 + 4143360) -
       53760 * tmp8 *
           (7680 * a6 * (tmp173 + 4 * tmp19) + (113485217444961.0 / 1000000000.0) * r * (tmp173 - 4 * tmp19) +
            (7402203300817.0 / 50000000000.0) * r * (7704 * r + 1396 * tmp19 + 5778 * tmp22 + 7704) +
            128 * r *
                ((852939255472061.0 / 50000000000.0) * r - 137046962773603.0 / 5000000000.0 * tmp19 + (127940888320809.0 / 10000000000.0) * tmp22 +
                 33046962773603.0 / 2500000000.0) +
            (23377610738487781.0 / 6250000000.0) * r - 46945920007815877.0 / 50000000000.0 * tmp14 + (26449610738487781.0 / 25000000000.0) * tmp19 +
            (26449610738487797.0 / 12500000000.0) * tmp22 + 1805060293831937.0 / 625000000.0)) /
      ((tmp79) * (tmp79));
  const REAL tmp202 = pphi * tmp121 * tmp41 + (1.0 / 2.0) * tmp198;
  const REAL tmp83 = tmp80 * tmp82;
  const REAL tmp123 = tmp114 * tmp121 * tmp40 + tmp114 + 1;
  const REAL tmp169 = tmp80 * (2048 * tmp105 * tmp81 + (1300341648852967.0 / 1000000000.0) * tmp12 + 688128 * tmp73 -
                               6572428801094219.0 / 1000000000.0 * tmp8 + 1720320);
  const REAL tmp85 = tmp41 + tmp83 * tmp84;
  const REAL tmp170 = tmp19 * tmp83;
  const REAL tmp86 = (1.0 / ((tmp85) * (tmp85)));
  const REAL tmp99 = tmp15 * tmp94 + tmp20 * tmp93 + tmp41 + tmp83 * tmp96 * tmp98;
  const REAL tmp108 = (1.0 / ((tmp85) * (tmp85) * (tmp85) * (tmp85)));
  const REAL tmp129 = tmp106 * tmp127 + tmp128 * tmp15 + tmp85;
  const REAL tmp124 = -113175085791863.0 / 10000000000000.0 * tmp100 * tmp26 + (29655068404873.0 / 20000000000000.0) * tmp101 * tmp32 * tmp8 +
                      (73721876495073.0 / 500000000000.0) * tmp103 * tmp16 * tmp8 + (121954868780449.0 / 1000000000000000.0) * tmp105 * tmp32 +
                      tmp106 * tmp107 * tmp108 * tmp109 * tmp112 * tmp16 + tmp123 + tmp13 * tmp15 * tmp16 + tmp16 * tmp18 * tmp20 +
                      tmp16 * tmp23 * tmp24 + tmp20 * tmp26 * tmp28 + tmp23 * tmp26 * tmp30 + tmp23 * tmp32 * tmp33 +
                      tmp35 * tmp43 * tmp52 * tmp53 * tmp60 * tmp86 * tmp99;
  const REAL tmp130 = tmp126 * tmp129;
  const REAL tmp175 = tmp108 * tmp109 * tmp112;
  const REAL tmp183 = tmp52 * tmp60 * tmp86;
  const REAL tmp131 = sqrt(tmp124 * tmp130);
  const REAL tmp203 = sqrt(tmp123 * tmp130);
  const REAL tmp157 = sqrt(2 * tmp8 * (tmp131 + tmp155) + 1);
  const REAL tmp191 = tmp131 / tmp124;
  const REAL tmp158 = (1.0 / (tmp157));
  const REAL Hreal = tmp157;
  const REAL xi = tmp85 * sqrt(r * tmp96) / tmp42;
  const REAL dHreal_dr =
      tmp158 *
      (tmp125 * tmp191 *
           ((1.0 / 2.0) * tmp124 * tmp126 * (-tmp127 * tmp168 - tmp128 * tmp163 + tmp167 + tmp169 * tmp84 + 30720 * tmp170 + tmp174 * tmp84) +
            (1.0 / 2.0) * tmp124 * tmp129 * (tmp116 * tmp161 * tmp40 + 2 * tmp15 * tmp40) / ((tmp125) * (tmp125)) +
            (1.0 / 2.0) * tmp130 *
                ((792225600543041.0 / 20000000000000.0) * m1 * m2 * tmp103 * tmp26 * tmp6 - 29655068404873.0 / 8000000000000.0 * tmp100 * tmp32 +
                 tmp106 * tmp107 * tmp108 * tmp109 * tmp110 * tmp16 * tmp52 *
                     (-36864 * r * tmp46 + tmp115 * tmp51 + (6739619044665529558.0 / 78125.0) * tmp12 + 2 * tmp181 - 5806080 * tmp22 * tmp47 +
                      (296393260946151.0 / 25.0) * tmp4 + 11575876386816 * tmp44 + (866182644304933.0 / 5.0) * tmp48 + 23151752773632 * tmp5 +
                      (6780203321721046624.0 / 78125.0) * tmp8 + 376881577556392) +
                 tmp106 * tmp107 * tmp108 * tmp109 * tmp111 * tmp16 *
                     (-535088332800 * r * tmp12 - 510774533137571.0 / 25.0 * r - 182268054644921.0 / 50.0 * tmp12 + 2 * tmp178 * tmp5 -
                      23151752773632 * tmp5 * tmp53 + 2 * tmp53 * tmp58 - 258910106287381.0 / 25.0 * tmp54 - 29035754444081779.0 / 2500.0 * tmp8 +
                      26750842527187002873.0 / 250000.0) /
                     ((tmp59) * (tmp59) * (tmp59)) +
                 tmp106 * tmp107 * tmp109 * tmp110 * tmp111 * tmp16 *
                     (-30720 * tmp14 * tmp169 - 30720 * tmp14 * tmp174 - 122880 * tmp170 + 8 * tmp20 * tmp39 * tmp6) /
                     ((tmp85) * (tmp85) * (tmp85) * (tmp85) * (tmp85)) -
                 tmp107 * tmp108 * tmp109 * tmp112 * tmp16 * tmp168 -
                 tmp113 * tmp116 * tmp23 * tmp39 * tmp6 * (-2 * r + 2 * tmp23 * tmp39 * tmp6) / ((tmp118) * (tmp118)) +
                 2 * tmp113 * tmp119 * tmp15 * tmp39 * tmp6 - tmp121 * tmp160 * tmp40 - tmp13 * tmp16 * tmp163 - 3 * tmp15 * tmp16 * tmp18 +
                 tmp15 * tmp16 * (-516952380952381.0 / 10000000000000.0 * tmp105 + 118.40000000000001 * tmp12 * tmp53) -
                 4 * tmp15 * tmp183 * tmp35 * tmp40 * tmp42 * tmp99 - tmp159 * tmp26 * tmp28 - 2 * tmp16 * tmp20 * tmp24 - tmp160 -
                 tmp161 * tmp26 * tmp30 - tmp161 * tmp32 * tmp33 - 121954868780449.0 / 1000000000000000.0 * tmp23 * tmp32 * tmp8 +
                 tmp35 * tmp43 * tmp52 * tmp53 * tmp60 * tmp86 *
                     (38400 * tmp14 * tmp83 * tmp96 - tmp159 * tmp93 - tmp163 * tmp94 + tmp167 + tmp169 * tmp96 * tmp98 + tmp174 * tmp96 * tmp98 +
                      tmp187 * tmp83 * tmp95 * tmp98 - tmp186 * tmp59 * tmp83 * tmp98 / tmp111) +
                 tmp35 * tmp43 * tmp52 * tmp53 * tmp60 * tmp99 *
                     (-15360 * tmp14 * tmp169 - 15360 * tmp14 * tmp174 - 61440 * tmp170 + 4 * tmp20 * tmp39 * tmp6) / ((tmp85) * (tmp85) * (tmp85)) +
                 tmp35 * tmp43 * tmp86 * tmp99 * (-tmp110 * tmp187 * tmp52 * tmp53 + tmp186 * tmp53 * tmp60 - tmp23 * tmp52 * tmp60) -
                 8 * tmp16 * tmp175 * tmp40 * ((tmp42) * (tmp42) * (tmp42)) / pow(r, 8) -
                 663496888455657.0 / 1000000000000.0 * tmp16 * tmp8 / pow(r, 11.0 / 2.0))) /
           tmp129 +
       tmp133 * (-pphi * tmp135 * tmp159 - 1.0 / 2.0 * pphi * tmp136 * tmp20 +
                 pphi * (tmp146 * (-tmp113 * tmp141 * tmp159 - tmp138 * tmp23 - tmp139 * tmp160 - tmp140 * tmp161 - tmp143 * tmp144 * tmp163) +
                         tmp152 * (-tmp113 * tmp150 * tmp159 - tmp144 * tmp151 * tmp163 - tmp147 * tmp23 - tmp148 * tmp160 - tmp149 * tmp161))) +
       tmp154 * (-3 * tmp22 - tmp40) / ((tmp132) * (tmp132)));
  const REAL dHreal_dprstar =
      (1.0 / 2.0) * tmp158 * tmp191 *
      (prstar * tmp115 * tmp183 * tmp43 * tmp99 - 339525257375589.0 / 5000000000000.0 * tmp100 * tmp193 +
       (29655068404873.0 / 2500000000000.0) * tmp101 * tmp194 * tmp8 + (73721876495073.0 / 125000000000.0) * tmp103 * tmp192 * tmp8 +
       (121954868780449.0 / 125000000000000.0) * tmp105 * tmp194 + tmp107 * tmp163 * tmp175 * tmp192 + tmp13 * tmp15 * tmp195 +
       tmp18 * tmp195 * tmp20 + 6 * tmp193 * tmp20 * tmp28 + 6 * tmp193 * tmp23 * tmp30 + 8 * tmp194 * tmp23 * tmp33 + tmp195 * tmp23 * tmp24);
  const REAL dHreal_dpphi = tmp158 * (tmp191 * tmp202 + tmp201);
  const REAL dHreal_dpphi_circ = (tmp201 + tmp202 * tmp203 / tmp123) / sqrt(2 * tmp8 * (tmp155 + tmp203) + 1);

  REAL flux[2];
  SEOBNRv5_aligned_spin_flux(y, Hreal, dHreal_dpphi, dHreal_dpphi_circ, flux, params);
  f[0] = xi * dHreal_dprstar;
  f[1] = dHreal_dpphi;
  f[2] = -xi * dHreal_dr + flux[0];
  f[3] = flux[1];
  return GSL_SUCCESS;
} // END FUNCTION SEOBNRv5_aligned_spin_right_hand_sides
