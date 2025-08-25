#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Evaluate SEOBNRv5 radial momentum condition.
 */
REAL SEOBNRv5_aligned_spin_radial_momentum_conditions(REAL x, void *restrict params) {

  const REAL m1 = ((commondata_struct *restrict)params)->m1;
  const REAL m2 = ((commondata_struct *restrict)params)->m2;
  const REAL chi1 = ((commondata_struct *restrict)params)->chi1;
  const REAL chi2 = ((commondata_struct *restrict)params)->chi2;
  const REAL a6 = ((commondata_struct *restrict)params)->a6;
  const REAL dSO = ((commondata_struct *restrict)params)->dSO;
  const REAL r = ((commondata_struct *restrict)params)->r;
  const REAL pphi = ((commondata_struct *restrict)params)->pphi;
  REAL prstar = x;
  const REAL tmp0 = ((r) * (r) * (r));
  const REAL tmp1 = m1 + m2;
  const REAL tmp11 = ((r) * (r) * (r) * (r));
  const REAL tmp25 = ((r) * (r));
  const REAL tmp30 = ((pphi) * (pphi));
  const REAL tmp33 = ((m1) * (m1));
  const REAL tmp34 = ((m2) * (m2));
  const REAL tmp41 = ((pphi) * (pphi) * (pphi) * (pphi));
  const REAL tmp42 = ((r) * (r) * (r) * (r) * (r));
  const REAL tmp57 = (1.0 / (r));
  const REAL tmp65 = ((m1) * (m1) * (m1));
  const REAL tmp66 = ((m2) * (m2) * (m2));
  const REAL tmp69 = log(r);
  const REAL tmp71 = ((prstar) * (prstar) * (prstar) * (prstar));
  const REAL tmp77 = pow(prstar, 6);
  const REAL tmp83 = pow(prstar, 8);
  const REAL tmp86 = ((prstar) * (prstar));
  const REAL tmp112 = 7680 * a6;
  const REAL tmp113 = 128 * r;
  const REAL tmp114 = (113485217444961.0 / 1000000000.0) * r;
  const REAL tmp116 = (7402203300817.0 / 50000000000.0) * r;
  const REAL tmp143 = pow(r, -5.0 / 2.0);
  const REAL tmp145 = pow(r, -9.0 / 2.0);
  const REAL tmp180 = pow(r, -6);
  const REAL tmp187 = -3566651379711.0 / 100000.0 * r - 163683964822551.0 / 1000000.0;
  const REAL tmp203 = pow(r, -11.0 / 2.0);
  const REAL tmp212 = pow(r, -8);
  const REAL tmp296 = pow(r, -7);
  const REAL tmp304 = 12 * r + 8;
  const REAL tmp366 = ((pphi) * (pphi) * (pphi));
  const REAL tmp378 = ((prstar) * (prstar) * (prstar));
  const REAL tmp379 = ((prstar) * (prstar) * (prstar) * (prstar) * (prstar));
  const REAL tmp380 = pow(prstar, 7);
  const REAL tmp6 = chi1 * m1 + chi2 * m2;
  const REAL tmp12 = (1.0 / (tmp11));
  const REAL tmp15 = (1.0 / ((tmp1) * (tmp1) * (tmp1)));
  const REAL tmp19 = (1.0 / (tmp0));
  const REAL tmp20 = (1.0 / ((tmp1) * (tmp1) * (tmp1) * (tmp1)));
  const REAL tmp21 = chi1 * m1 - chi2 * m2;
  const REAL tmp26 = (1.0 / (tmp25));
  const REAL tmp43 = (1.0 / (tmp42));
  const REAL tmp67 = pow(tmp1, -6);
  const REAL tmp78 = ((m1) * (m1) * (m1) * (m1)) * ((m2) * (m2) * (m2) * (m2)) / pow(tmp1, 8);
  const REAL tmp90 = ((tmp69) * (tmp69));
  const REAL tmp111 = 8 * r + 2 * tmp0 + 4 * tmp25 + 16;
  const REAL tmp130 = 7680 * tmp11;
  const REAL tmp140 = 7680 * tmp42;
  const REAL tmp157 = 2 * tmp57;
  const REAL tmp188 = tmp57 * tmp69;
  const REAL tmp190 = -117964800 * a6 + (206478381132587.0 / 50000.0) * r - 308925070376879.0 / 25000.0 * tmp0 + (3566651379711.0 / 500.0) * tmp11 -
                      926775211130637.0 / 50000.0 * tmp25 + 122635399361987.0 / 1000.0;
  const REAL tmp237 = 30720 * tmp11;
  const REAL tmp256 = 15360 * tmp11;
  const REAL tmp288 = 20 * tmp180;
  const REAL tmp357 = 15360 * tmp42;
  const REAL tmp3 = (1.0 / ((tmp1) * (tmp1)));
  const REAL tmp7 = ((tmp6) * (tmp6));
  const REAL tmp13 = 3 * tmp12;
  const REAL tmp22 = tmp21 * (m1 - m2);
  const REAL tmp31 = tmp19 * tmp30;
  const REAL tmp36 = 2 * tmp19;
  const REAL tmp37 = tmp20 * tmp33 * tmp34;
  const REAL tmp44 = 4 * tmp43;
  const REAL tmp46 = tmp6 / tmp1;
  const REAL tmp58 = tmp26 * tmp30;
  const REAL tmp68 = tmp65 * tmp66 * tmp67;
  const REAL tmp72 = tmp12 * tmp71;
  const REAL tmp91 = 5787938193408 * tmp90;
  const REAL tmp155 = tmp43 * tmp71;
  const REAL tmp158 = tmp157 + 1;
  const REAL tmp196 = 8 * r + 6 * tmp25 + 8;
  const REAL tmp218 = 23151752773632 * tmp57 * tmp69;
  const REAL tmp244 = 11575876386816 * tmp57 * tmp69;
  const REAL tmp254 = 4 * tmp19;
  const REAL tmp283 = 12 * tmp43;
  const REAL tmp285 = 6 * tmp12;
  const REAL tmp314 = tmp20 * ((tmp6) * (tmp6) * (tmp6) * (tmp6));
  const REAL tmp333 = 16 * tmp43;
  const REAL tmp363 = ((tmp1) * (tmp1)) / (m1 * m2);
  const REAL tmp368 = 4 * tmp12;
  const REAL tmp8 = tmp3 * tmp7;
  const REAL tmp17 = dSO * m1 * m2 * tmp15 * tmp6;
  const REAL tmp23 = -tmp15 * ((tmp6) * (tmp6) * (tmp6)) + tmp20 * tmp22 * tmp7;
  const REAL tmp27 = m1 * m2 * tmp3;
  const REAL tmp32 = 2 * tmp31;
  const REAL tmp52 = tmp22 * tmp3;
  const REAL tmp75 = 8 * m1 * m2 * tmp3 - 6 * tmp37;
  const REAL tmp109 = 7 * tmp37;
  const REAL tmp133 = ((tmp21) * (tmp21)) * tmp3;
  const REAL tmp134 = tmp15 * tmp22 * tmp6;
  const REAL tmp197 = (23377610738487781.0 / 6250000000.0) * r + (26449610738487781.0 / 25000000000.0) * tmp0 -
                      46945920007815877.0 / 50000000000.0 * tmp11 + tmp112 * (4 * tmp0 + tmp196) +
                      tmp113 * ((852939255472061.0 / 50000000000.0) * r - 137046962773603.0 / 5000000000.0 * tmp0 +
                                (127940888320809.0 / 10000000000.0) * tmp25 + 33046962773603.0 / 2500000000.0) +
                      tmp114 * (-4 * tmp0 + tmp196) + tmp116 * (7704 * r + 1396 * tmp0 + 5778 * tmp25 + 7704) +
                      (26449610738487797.0 / 12500000000.0) * tmp25 + 1805060293831937.0 / 625000000.0;
  const REAL tmp204 = 118.40000000000001 * tmp37;
  const REAL tmp208 = -2 * r + 2 * tmp26 * tmp3 * tmp7;
  const REAL tmp286 = tmp285 * tmp30;
  const REAL tmp287 = tmp283 * tmp30;
  const REAL tmp364 = 2 * pphi * tmp26;
  const REAL tmp376 = pphi * tmp254;
  const REAL tmp9 = tmp0 + tmp8 * (r + 2);
  const REAL tmp24 = (1.0 / 2.0) * tmp19 * tmp23;
  const REAL tmp28 = (23.0 / 32.0) * tmp27 - 3.0 / 32.0;
  const REAL tmp29 = -45.0 / 32.0 * tmp27 - 15.0 / 32.0;
  const REAL tmp35 = (109.0 / 192.0) * tmp20 * tmp33 * tmp34 - 177.0 / 32.0 * tmp27 - 5.0 / 64.0;
  const REAL tmp38 = -267.0 / 128.0 * tmp27 - 1591.0 / 768.0 * tmp37 + 59.0 / 256.0;
  const REAL tmp40 = (75.0 / 128.0) * tmp27 + (345.0 / 256.0) * tmp37 + 105.0 / 256.0;
  const REAL tmp47 = (11.0 / 32.0) * tmp27 + 3.0 / 32.0;
  const REAL tmp48 = 15.0 / 32.0 - 9.0 / 32.0 * tmp27;
  const REAL tmp49 = -1.0 / 32.0 * tmp27 + (103.0 / 192.0) * tmp37 + 5.0 / 64.0;
  const REAL tmp50 = -35.0 / 128.0 * tmp27 - 613.0 / 768.0 * tmp37 - 59.0 / 256.0;
  const REAL tmp51 = (75.0 / 256.0) * tmp20 * tmp33 * tmp34 - 45.0 / 128.0 * tmp27 - 105.0 / 256.0;
  const REAL tmp56 = (1.0 / 4.0) * tmp23 * tmp26;
  const REAL tmp63 = -3 * tmp25 - tmp8;
  const REAL tmp70 = tmp27 * (452542166996693.0 / 1000000000000.0 - 516952380952381.0 / 10000000000000.0 * tmp69) +
                     tmp37 * (118.40000000000001 * tmp69 - 179613660498019.0 / 100000000000.0) + (150579635104141.0 / 250000000000.0) * tmp68;
  const REAL tmp73 = (115888805356193.0 / 1250000000000.0) * tmp27 - 131 * tmp37 + 10 * tmp68;
  const REAL tmp79 =
      -84945530542609.0 / 2500000000000.0 * tmp27 - 447649163680617.0 / 5000000000000.0 * tmp37 + 188 * tmp65 * tmp66 * tmp67 - 14 * tmp78;
  const REAL tmp81 = -139150381847503.0 / 50000000000000.0 * tmp27 - 27.0 / 5.0 * tmp37 + 6 * tmp65 * tmp66 * tmp67;
  const REAL tmp84 = (4343054718629.0 / 3125000000000.0) * tmp27 + (166921011824161.0 / 50000000000000.0) * tmp37 +
                     (342857142857143.0 / 100000000000000.0) * tmp68 - 6 * tmp78;
  const REAL tmp87 = tmp26 * tmp8;
  const REAL tmp92 = -7876452485916241.0 / 12500000.0 * tmp27 - 98886748396767.0 / 500000.0 * tmp37 + 580530436787913.0 / 100000.0;
  const REAL tmp93 = -53511513369581.0 / 20000000.0 * tmp27 - 138240 * tmp37 - 52783413229329.0 / 10000000.0;
  const REAL tmp94 = ((1 - 0.49694878161693501 * tmp27) * (1 - 0.49694878161693501 * tmp27));
  const REAL tmp96 = 14700 * tmp27 + 42911;
  const REAL tmp117 = 53760 * tmp27;
  const REAL tmp118 = tmp27 * tmp69;
  const REAL tmp119 = 588 * tmp27 + 1079;
  const REAL tmp135 = tmp133 * ((3.0 / 4.0) * tmp27 - 3.0 / 16.0) - 21.0 / 8.0 * tmp134 + tmp8 * (3 * tmp27 + 45.0 / 16.0);
  const REAL tmp136 = tmp133 * ((115.0 / 64.0) * tmp27 + (1.0 / 16.0) * tmp37 - 37.0 / 64.0) + tmp134 * ((13.0 / 16.0) * tmp27 + 449.0 / 32.0) +
                      tmp8 * (-1171.0 / 64.0 * tmp27 - 861.0 / 64.0);
  const REAL tmp142 = tmp27 / pow(r, 7.0 / 2.0);
  const REAL tmp144 = tmp27 * tmp83;
  const REAL tmp146 = tmp27 * tmp71;
  const REAL tmp147 = tmp27 * tmp57;
  const REAL tmp151 = tmp133 * ((75.0 / 32.0) * m1 * m2 * tmp3 - 15.0 / 8.0 * tmp37 - 15.0 / 32.0) + tmp134 * ((45.0 / 8.0) * tmp27 - 5.0 / 16.0) +
                      tmp8 * ((165.0 / 32.0) * tmp27 - 5 * tmp37 + 25.0 / 32.0);
  const REAL tmp159 = tmp158 * tmp8;
  const REAL tmp164 = -tmp158 * tmp58;
  const REAL tmp168 = tmp133 * ((21.0 / 16.0) * tmp20 * tmp33 * tmp34 - 81.0 / 64.0 * tmp27 - 9.0 / 64.0) +
                      tmp134 * (117.0 / 32.0 - 39.0 / 16.0 * tmp27) + tmp8 * (-175.0 / 64.0 * tmp27 - 225.0 / 64.0);
  const REAL tmp169 = tmp133 * ((1.0 / 2.0) * tmp27 + 1.0 / 8.0) - 5.0 / 4.0 * tmp134 + (9.0 / 8.0) * tmp8;
  const REAL tmp174 = 2 * tmp8;
  const REAL tmp179 = -tmp36 * tmp8;
  const REAL tmp205 = tmp204 * tmp57 - 516952380952381.0 / 10000000000000.0 * tmp27 * tmp57;
  const REAL tmp221 = -1759846072320 * r - 4161798144000 * tmp27 - 12148770078720;
  const REAL tmp232 = tmp12 * tmp8;
  const REAL tmp255 = tmp254 * tmp8;
  const REAL tmp326 = tmp69 * (5279538216960 * r + 8323596288000 * tmp27 + 24297540157440);
  const REAL tmp10 = (1.0 / (tmp9));
  const REAL tmp53 = tmp46 * (-tmp13 * tmp30 * tmp38 - tmp26 * tmp28 - tmp29 * tmp32 - tmp35 * tmp36 - tmp40 * tmp41 * tmp44) +
                     tmp52 * (-tmp13 * tmp30 * tmp50 - tmp26 * tmp47 - tmp32 * tmp48 - tmp36 * tmp49 - tmp41 * tmp44 * tmp51);
  const REAL tmp60 = tmp46 * (tmp12 * tmp40 * tmp41 + tmp26 * tmp35 + tmp28 * tmp57 + tmp29 * tmp58 + tmp31 * tmp38 + 7.0 / 4.0) +
                     tmp52 * (tmp12 * tmp41 * tmp51 + tmp26 * tmp49 + tmp31 * tmp50 + tmp47 * tmp57 + tmp48 * tmp58 + 1.0 / 4.0);
  const REAL tmp62 = (1.0 / ((tmp9) * (tmp9)));
  const REAL tmp80 = tmp77 * tmp79;
  const REAL tmp82 = tmp77 * tmp81;
  const REAL tmp85 = tmp83 * tmp84;
  const REAL tmp88 = tmp87 + 1;
  const REAL tmp95 = (3390101660860523312.0 / 78125.0) * tmp27 + (3369809522332764779.0 / 78125.0) * tmp37 + (296393260946151.0 / 50.0) * tmp68 +
                     (866182644304933.0 / 10.0) * tmp94 + 188440788778196;
  const REAL tmp97 = 49152 * r * ((16368307925443.0 / 40000.0) * tmp27 + 102574080 * tmp37 - 105983591868019.0 / 50000.0) + 879923036160 * tmp0 +
                     283115520 * tmp25 * tmp96 - 1698693120 * tmp27 * (11592 * tmp27 + 69847);
  const REAL tmp100 = r * tmp96;
  const REAL tmp121 = 756 * tmp27 + 1079;
  const REAL tmp127 = 336 * r + 756 * tmp27 + 407;
  const REAL tmp160 = tmp159 + tmp25;
  const REAL tmp171 = tmp158 * tmp87 + 1;
  const REAL tmp175 = tmp12 * tmp174 + tmp159 * tmp36;
  const REAL tmp233 = 4 * tmp232;
  const REAL tmp372 = 2 * pphi * tmp87;
  const REAL tmp54 = -pphi * tmp13 * tmp17 - pphi * tmp24 + pphi * tmp53;
  const REAL tmp61 = pphi * tmp17 * tmp19 + pphi * tmp56 + pphi * tmp60;
  const REAL tmp64 = tmp62 * tmp63;
  const REAL tmp89 = ((tmp88) * (tmp88));
  const REAL tmp98 =
      r * tmp91 + r * tmp95 - 967680 * tmp0 * tmp93 - 9216 * tmp25 * tmp92 +
      55296 * tmp27 *
          (-31383302728516087.0 / 12500000.0 * tmp27 - 426364516032331.0 / 10000000.0 * tmp37 + 14515200 * tmp68 + 100201376401019.0 / 100000.0) +
      tmp69 * tmp97;
  const REAL tmp103 = -283115520 * tmp100 - 879923036160 * tmp25 - 253929149957443.0 / 10.0 * tmp27 - 5041721180160 * tmp37 + 104186110149937;
  const REAL tmp128 = 2048 * m1 * m2 * tmp127 * tmp3 * tmp69 + 28 * m1 * m2 * tmp3 * (1920 * a6 + 733955307463037.0 / 1000000000.0) -
                      7 * r * ((938918400156317.0 / 1000000000.0) * m1 * m2 * tmp3 - 185763092693281.0 / 1000000000.0 * tmp37 - 245760) -
                      270820329770593.0 / 50000000.0 * tmp37 - 3440640;
  const REAL tmp148 = ((tmp88) * (tmp88) * (tmp88) * (tmp88));
  const REAL tmp161 = (1.0 / (tmp160));
  const REAL tmp178 = (1.0 / (tmp171));
  const REAL tmp209 = (1.0 / ((tmp160) * (tmp160)));
  const REAL tmp213 = ((tmp88) * (tmp88) * (tmp88));
  const REAL tmp226 =
      566231040 * tmp100 + 2639769108480 * tmp25 + (12570860486740224.0 / 625.0) * tmp27 + 5041721180160 * tmp37 - 325581594218554368.0 / 3125.0;
  const REAL tmp234 = tmp86 * tmp88;
  const REAL tmp370 = pphi * (tmp46 * (pphi * tmp36 * tmp38 + tmp29 * tmp364 + tmp366 * tmp368 * tmp40) +
                              tmp52 * (pphi * tmp36 * tmp50 + tmp364 * tmp48 + tmp366 * tmp368 * tmp51)) +
                      tmp17 * tmp19 + tmp56 + tmp60;
  const REAL tmp104 = (53058305272831.0 / 5.0) * r * tmp27 + (182268054644921.0 / 100.0) * r * tmp37 - 197383821284628.0 / 5.0 * r -
                      163418713120743.0 / 500000.0 * tmp100 - tmp103 * tmp69 + (258910106287381.0 / 100.0) * tmp25 * tmp27 +
                      133772083200 * tmp25 * tmp37 + (510774533137571.0 / 100.0) * tmp25 - 60249543508726.0 / 5.0 * tmp27 +
                      (400296247701391.0 / 5.0) * tmp37 + (336524885906151.0 / 50.0) * tmp68 + tmp91 + 275059053208689;
  const REAL tmp123 = 240 * r * (-373313530533103.0 / 50000000000.0 * tmp27 - 3024 * tmp37 + 17264) + 960 * tmp0 * tmp121 + 480 * tmp11 * tmp121 +
                      1920 * tmp119 * tmp25 - 388422414769507.0 / 10000000.0 * tmp27 - 47061405915993.0 / 25000000.0 * tmp37 + 161280 * tmp42 +
                      13447680;
  const REAL tmp137 = (1.0 / (tmp98));
  const REAL tmp150 = ((tmp98) * (tmp98));
  const REAL tmp162 = tmp161 * tmp8;
  const REAL tmp173 = (1.0 / ((tmp171) * (tmp171)));
  const REAL tmp183 =
      688128 * tmp118 + 2048 * tmp127 * tmp147 - 6572428801094219.0 / 1000000000.0 * tmp27 + (1300341648852967.0 / 1000000000.0) * tmp37 + 1720320;
  const REAL tmp192 = 3840 * r * tmp119 + 1920 * tmp0 * tmp121 + 806400 * tmp11 + 2880 * tmp121 * tmp25 - 1119940591599309.0 / 625000000.0 * tmp27 -
                      725760 * tmp37 + 4143360;
  const REAL tmp207 = -tmp158 * tmp161;
  const REAL tmp222 = -258910106287381.0 / 25.0 * r * tmp27 - 535088332800 * r * tmp37 - 510774533137571.0 / 25.0 * r + 2 * tmp103 * tmp57 - tmp218 +
                      2 * tmp221 * tmp69 - 29035754444081779.0 / 2500.0 * tmp27 - 182268054644921.0 / 50.0 * tmp37 +
                      26750842527187002873.0 / 250000.0;
  const REAL tmp247 = (258910106287381.0 / 50.0) * r * tmp27 + 267544166400 * r * tmp37 + (510774533137571.0 / 50.0) * r - tmp103 * tmp57 -
                      tmp221 * tmp69 + tmp244 + (29035754444081779.0 / 5000.0) * tmp27 + (182268054644921.0 / 100.0) * tmp37 -
                      26750842527187002873.0 / 500000.0;
  const REAL tmp250 = tmp57 * tmp98;
  const REAL tmp277 = tmp10 * tmp61 - 1;
  const REAL tmp316 = tmp148 * tmp151;
  const REAL tmp324 = 16 * tmp212 * tmp213 * tmp8;
  const REAL tmp340 = tmp86 * tmp89;
  const REAL tmp349 = -5806080 * r * tmp93 + tmp157 * tmp226 + tmp244 - tmp26 * tmp97 + (4536836631887754816.0 / 390625.0) * tmp27 + tmp326 +
                      (56958767076537792.0 / 15625.0) * tmp37 + 11575876386816 * tmp57 - 334385531589837888.0 / 3125.0;
  const REAL tmp105 = (1.0 / (tmp104));
  const REAL tmp138 = tmp104 * tmp137;
  const REAL tmp149 = (1.0 / ((tmp104) * (tmp104)));
  const REAL tmp165 = tmp162 * tmp164 + tmp58 + 1;
  const REAL tmp223 = (1.0 / ((tmp104) * (tmp104) * (tmp104)));
  const REAL tmp228 = -36864 * r * tmp92 + tmp157 * tmp97 + 2 * tmp226 * tmp69 - 5806080 * tmp25 * tmp93 + (6780203321721046624.0 / 78125.0) * tmp27 +
                      (6739619044665529558.0 / 78125.0) * tmp37 + (296393260946151.0 / 25.0) * tmp68 + 23151752773632 * tmp69 +
                      11575876386816 * tmp90 + (866182644304933.0 / 5.0) * tmp94 + 376881577556392;
  const REAL tmp240 = -18432 * r * tmp92 + tmp226 * tmp69 - 2903040 * tmp25 * tmp93 + tmp57 * tmp97 + 11575876386816 * tmp69 + tmp91 + tmp95;
  const REAL tmp258 = tmp137 * tmp247;
  const REAL tmp261 = (1.0 / (tmp150));
  const REAL tmp310 = -tmp158 * tmp162;
  const REAL tmp313 = -tmp158 * tmp208 * tmp209;
  const REAL tmp350 = tmp103 * tmp26 - 2 * tmp221 * tmp57 - 11575876386816 * tmp26 * tmp69 + 11575876386816 * tmp26 +
                      (258910106287381.0 / 50.0) * tmp27 + 267544166400 * tmp37 + 1759846072320 * tmp69 + 510774533137571.0 / 50.0;
  const REAL tmp373 = tmp207 * tmp372 + tmp364;
  const REAL tmp106 = tmp105 * tmp57;
  const REAL tmp125 = tmp109 * (-39321600 * a6 * (3 * r + 59) + (186464462028901.0 / 250000.0) * a6 + (122635399361987.0 / 1000.0) * r -
                                308925070376879.0 / 50000.0 * tmp0 - 308925070376879.0 / 100000.0 * tmp11 + (206478381132587.0 / 100000.0) * tmp25 +
                                (3566651379711.0 / 2500.0) * tmp42 + 276057889687011.0 / 1000.0) +
                      tmp117 * (tmp112 * (tmp11 + tmp111) +
                                tmp113 * ((33046962773603.0 / 2500000000.0) * r + (42646962773603.0 / 10000000000.0) * tmp0 -
                                          137046962773603.0 / 20000000000.0 * tmp11 + (852939255472061.0 / 100000000000.0) * tmp25 -
                                          42153037226397.0 / 1250000000.0) +
                                tmp114 * (-tmp11 + tmp111) + tmp116 * (7704 * r + 1926 * tmp0 + 349 * tmp11 + 3852 * tmp25 + 36400)) +
                      32768 * tmp118 * tmp123 + 67645734912 * tmp37 * tmp90 + 13212057600 * tmp42 +
                      1120 * tmp68 * (-163683964822551.0 / 1000000.0 * r - 3566651379711.0 / 200000.0 * tmp25 - 59449372951581.0 / 50000.0) +
                      (241555486248807.0 / 1000.0) * tmp78;
  const REAL tmp152 = tmp149 * tmp150 * tmp151;
  const REAL tmp198 = -tmp109 * tmp190 - 66060288000 * tmp11 - tmp117 * tmp197 - 32768 * tmp118 * tmp192 - 32768 * tmp123 * tmp147 -
                      1120 * tmp187 * tmp68 - 135291469824 * tmp188 * tmp37;
  const REAL tmp229 = tmp105 * tmp98;
  const REAL tmp249 = -tmp149 * tmp247;
  const REAL tmp260 = tmp138 * tmp140;
  const REAL tmp319 = tmp222 * tmp223;
  const REAL tmp374 = (1.0 / 2.0) * tmp373;
  const REAL tmp107 = tmp106 * tmp98;
  const REAL tmp126 = (1.0 / (tmp125));
  const REAL tmp153 = tmp148 * tmp152;
  const REAL tmp186 = (1.0 / ((tmp125) * (tmp125)));
  const REAL tmp241 = tmp106 * tmp240;
  const REAL tmp243 = tmp105 * tmp26 * tmp98;
  const REAL tmp251 = tmp249 * tmp250;
  const REAL tmp263 = -tmp240 * tmp261;
  const REAL tmp308 = tmp198 *
                      (-132120576000 * tmp11 - 65536 * tmp118 * tmp192 - 65536 * tmp123 * tmp147 - 2240 * tmp187 * tmp68 -
                       270582939648 * tmp188 * tmp37 - 14 * tmp190 * tmp37 - 107520 * tmp197 * tmp27) /
                      ((tmp125) * (tmp125) * (tmp125));
  const REAL tmp330 = tmp149 * tmp228 * tmp98;
  const REAL tmp344 = 8 * tmp229 * tmp232;
  const REAL tmp129 = tmp126 * tmp128;
  const REAL tmp184 = tmp126 * tmp183;
  const REAL tmp252 = tmp241 - tmp243 + tmp251;
  const REAL tmp298 = tmp126 * (-2048 * tmp127 * tmp26 * tmp27 + 1376256 * tmp147);
  const REAL tmp306 =
      tmp186 *
      (32768 * m1 * m2 * tmp123 * tmp26 * tmp3 - 264241152000 * tmp0 -
       tmp109 *
           (-926775211130637.0 / 25000.0 * r + (3566651379711.0 / 125.0) * tmp0 - 926775211130637.0 / 25000.0 * tmp25 + 206478381132587.0 / 50000.0) -
       tmp117 * ((26449610738487797.0 / 3125000000.0) * r - 46945920007815877.0 / 6250000000.0 * tmp0 + tmp112 * (12 * tmp25 + tmp304) +
                 tmp113 * ((127940888320809.0 / 5000000000.0) * r - 411140888320809.0 / 5000000000.0 * tmp25 + 852939255472061.0 / 50000000000.0) +
                 tmp114 * (-12 * tmp25 + tmp304) + tmp116 * (11556 * r + 4188 * tmp25 + 7704) + (79348832215463343.0 / 12500000000.0) * tmp25 +
                 23377610738487781.0 / 3125000000.0) -
       32768 * tmp118 * (5760 * r * tmp121 + 3225600 * tmp0 + 5760 * tmp121 * tmp25 + 2257920 * tmp27 + 4143360) - 65536 * tmp147 * tmp192 +
       135291469824 * tmp20 * tmp26 * tmp33 * tmp34 * tmp69 - 135291469824 * tmp26 * tmp37 + (24966559657977.0 / 625.0) * tmp65 * tmp66 * tmp67);
  const REAL tmp309 = tmp128 * tmp308;
  const REAL tmp321 = tmp150 * tmp316 * tmp319;
  const REAL tmp131 = tmp129 * tmp130 + tmp87;
  const REAL tmp185 = tmp0 * tmp129;
  const REAL tmp200 = tmp128 * tmp186 * tmp198;
  const REAL tmp259 = tmp129 * tmp140;
  const REAL tmp299 = tmp0 * tmp184;
  const REAL tmp300 = tmp129 * tmp25;
  const REAL tmp301 = tmp183 * tmp186 * tmp198;
  const REAL tmp307 = tmp128 * tmp306;
  const REAL tmp345 = tmp11 * tmp184;
  const REAL tmp359 = 76800 * tmp11 * tmp129;
  const REAL tmp362 = tmp104 * tmp263 * tmp357;
  const REAL tmp132 = (1.0 / ((tmp131) * (tmp131)));
  const REAL tmp141 = tmp12 * tmp136 + tmp129 * tmp138 * tmp140 + tmp135 * tmp19 + tmp87;
  const REAL tmp154 = (1.0 / ((tmp131) * (tmp131) * (tmp131) * (tmp131)));
  const REAL tmp170 = tmp12 * tmp169 + tmp131 + tmp168 * tmp43;
  const REAL tmp201 = tmp130 * tmp184 + tmp130 * tmp200 - 5 * tmp168 * tmp180 - tmp169 * tmp44 + tmp179 + 30720 * tmp185;
  const REAL tmp236 = (1.0 / ((tmp131) * (tmp131) * (tmp131) * (tmp131) * (tmp131)));
  const REAL tmp238 = -tmp184 * tmp237 - 122880 * tmp185 + 8 * tmp19 * tmp3 * tmp7 - tmp200 * tmp237;
  const REAL tmp253 = (1.0 / ((tmp131) * (tmp131) * (tmp131)));
  const REAL tmp257 = -tmp184 * tmp256 - 61440 * tmp185 - tmp200 * tmp256 + tmp255;
  const REAL tmp264 = tmp104 * tmp259;
  const REAL tmp302 = tmp0 * tmp200;
  const REAL tmp347 = tmp11 * tmp200;
  const REAL tmp177 = tmp170 * tmp173 * tmp175;
  const REAL tmp202 = tmp178 * tmp201;
  const REAL tmp214 = tmp154 * tmp71;
  const REAL tmp230 = tmp132 * tmp141;
  const REAL tmp265 = 38400 * tmp11 * tmp129 * tmp138 - tmp13 * tmp135 - tmp136 * tmp44 + tmp179 + tmp184 * tmp260 + tmp200 * tmp260 +
                      tmp258 * tmp259 + tmp263 * tmp264;
  const REAL tmp268 = tmp170 * tmp178;
  const REAL tmp271 = (1.0 / (tmp170));
  const REAL tmp290 = (1.0 / ((tmp170) * (tmp170)));
  const REAL tmp352 = tmp132 * tmp340;
  const REAL tmp166 = tmp107 * tmp132 * tmp141 * tmp86 * tmp89 - 113175085791863.0 / 10000000000000.0 * tmp142 * tmp77 +
                      (29655068404873.0 / 20000000000000.0) * tmp143 * tmp144 + (73721876495073.0 / 500000000000.0) * tmp145 * tmp146 +
                      (121954868780449.0 / 1000000000000000.0) * tmp147 * tmp83 + tmp153 * tmp154 * tmp155 + tmp165 + tmp19 * tmp71 * tmp73 +
                      tmp19 * tmp80 + tmp26 * tmp71 * tmp75 + tmp26 * tmp82 + tmp26 * tmp85 + tmp70 * tmp72;
  const REAL tmp211 = tmp153 * tmp154 * tmp71;
  const REAL tmp231 = tmp229 * tmp230;
  const REAL tmp269 = (1.0 / 2.0) * tmp268;
  const REAL tmp317 = tmp154 * tmp155 * tmp316;
  const REAL tmp323 = 10 * tmp180 * tmp214;
  const REAL tmp325 = tmp151 * tmp214 * tmp324;
  const REAL tmp332 = tmp230 * tmp233 * tmp234;
  const REAL tmp336 = tmp236 * tmp238 * tmp71;
  const REAL tmp337 = 2 * tmp155 * tmp236 * tmp238;
  const REAL tmp342 = tmp253 * tmp257 * tmp340;
  const REAL tmp354 = tmp265 * tmp352;
  const REAL tmp382 = sqrt(tmp165 * tmp268);
  const REAL tmp167 = (1.0 / 2.0) * tmp166;
  const REAL tmp216 = tmp152 * tmp213 * tmp214 * tmp8;
  const REAL tmp272 = sqrt(tmp166 * tmp268);
  const REAL tmp343 = tmp141 * tmp342;
  const REAL tmp266 = -792225600543041.0 / 20000000000000.0 * m1 * m2 * tmp145 * tmp3 * tmp77 -
                      tmp105 * tmp132 * tmp265 * tmp57 * tmp86 * tmp89 * tmp98 - tmp105 * tmp141 * tmp253 * tmp257 * tmp57 * tmp86 * tmp89 * tmp98 -
                      2 * tmp12 * tmp161 * tmp3 * tmp30 * tmp7 - tmp12 * tmp205 * tmp71 + tmp13 * tmp80 - tmp132 * tmp141 * tmp252 * tmp86 * tmp89 +
                      (29655068404873.0 / 8000000000000.0) * tmp142 * tmp83 + (121954868780449.0 / 1000000000000000.0) * tmp144 * tmp26 +
                      (663496888455657.0 / 1000000000000.0) * tmp146 * tmp203 - tmp148 * tmp149 * tmp150 * tmp151 * tmp236 * tmp238 * tmp43 * tmp71 -
                      tmp148 * tmp149 * tmp151 * tmp154 * tmp228 * tmp43 * tmp71 * tmp98 -
                      tmp148 * tmp150 * tmp151 * tmp154 * tmp222 * tmp223 * tmp43 * tmp71 + tmp158 * tmp208 * tmp209 * tmp26 * tmp3 * tmp30 * tmp7 +
                      tmp174 * tmp207 * tmp31 + 5 * tmp180 * tmp211 + 2 * tmp19 * tmp71 * tmp75 + 8 * tmp212 * tmp216 + tmp231 * tmp233 * tmp234 +
                      tmp32 + tmp36 * tmp82 + tmp36 * tmp85 + tmp44 * tmp70 * tmp71 + 3 * tmp72 * tmp73;
  const REAL tmp273 = tmp272 / tmp166;
  const REAL tmp279 = 2 * tmp27 * (tmp272 + tmp277) + 1;
  const REAL tmp291 = tmp272 / ((tmp166) * (tmp166));
  const REAL tmp280 = tmp27 / pow(tmp279, 3.0 / 2.0);
  const REAL tmp281 = sqrt(tmp279);
  const REAL tmp375 = tmp10 * tmp370 + tmp273 * tmp374;
  const REAL tmp270 = tmp167 * tmp177 + tmp167 * tmp202 - tmp266 * tmp269;
  const REAL tmp275 = tmp171 * tmp271 * tmp273;
  const REAL tmp282 = (1.0 / (tmp281));
  const REAL tmp276 = tmp10 * tmp54 + tmp270 * tmp275 + tmp61 * tmp64;
  const REAL tmp293 = tmp171 * tmp270 * tmp271 * tmp291;
  const REAL dHreal_dr_dr =
      tmp363 *
      (-((tmp276) * (tmp276)) * tmp280 +
       tmp282 *
           (-6 * r * tmp61 * tmp62 +
            tmp10 * ((3.0 / 2.0) * pphi * tmp12 * tmp23 + pphi * tmp17 * tmp283 +
                     pphi * (tmp46 * (tmp28 * tmp36 + tmp285 * tmp35 + tmp286 * tmp29 + tmp287 * tmp38 + tmp288 * tmp40 * tmp41) +
                             tmp52 * (tmp285 * tmp49 + tmp286 * tmp48 + tmp287 * tmp50 + tmp288 * tmp41 * tmp51 + tmp36 * tmp47))) +
            ((tmp171) * (tmp171)) * ((tmp270) * (tmp270)) * tmp290 * tmp291 - tmp171 * tmp201 * tmp270 * tmp273 * tmp290 -
            tmp175 * tmp270 * tmp271 * tmp273 + tmp266 * tmp293 +
            tmp275 * (tmp166 * tmp173 * tmp175 * tmp201 + tmp167 * tmp170 * tmp173 * (-tmp159 * tmp285 - tmp283 * tmp8) +
                      tmp167 * tmp170 * tmp175 * (tmp159 * tmp254 + tmp233) / ((tmp171) * (tmp171) * (tmp171)) +
                      tmp167 * tmp178 *
                          (tmp130 * tmp298 + tmp130 * tmp307 + tmp130 * tmp309 + 30 * tmp168 * tmp296 + tmp169 * tmp288 + tmp256 * tmp301 +
                           tmp285 * tmp8 + 61440 * tmp299 + 92160 * tmp300 + 61440 * tmp302) -
                      tmp177 * tmp266 - tmp202 * tmp266 +
                      tmp269 * (-tmp105 * tmp240 * tmp332 +
                                tmp107 * tmp141 * tmp154 * tmp257 * tmp340 *
                                    (-92160 * tmp185 + 6 * tmp19 * tmp3 * tmp7 - 23040 * tmp345 - 23040 * tmp347) +
                                tmp107 * tmp141 * tmp253 * tmp340 *
                                    (-tmp128 * tmp256 * tmp306 - tmp128 * tmp256 * tmp308 - 12 * tmp232 - tmp237 * tmp301 - tmp256 * tmp298 -
                                     122880 * tmp299 - 184320 * tmp300 - 122880 * tmp302) +
                                tmp107 * tmp352 *
                                    (tmp104 * tmp263 * tmp359 + tmp129 * tmp247 * tmp263 * tmp357 + tmp135 * tmp283 + tmp136 * tmp288 +
                                     tmp137 * tmp259 * tmp350 + 153600 * tmp138 * tmp185 + tmp138 * tmp301 * tmp357 + 76800 * tmp138 * tmp345 +
                                     76800 * tmp138 * tmp347 + tmp184 * tmp258 * tmp357 + tmp184 * tmp362 + tmp200 * tmp258 * tmp357 +
                                     tmp200 * tmp362 + tmp228 * tmp240 * tmp264 / ((tmp98) * (tmp98) * (tmp98)) + tmp258 * tmp359 + tmp260 * tmp298 +
                                     tmp260 * tmp307 + tmp260 * tmp309 - tmp261 * tmp264 * tmp349 + tmp285 * tmp8) -
                                tmp132 * tmp141 * tmp234 * tmp252 * tmp255 - tmp132 * tmp234 * tmp265 * tmp344 -
                                tmp141 * tmp234 * tmp253 * tmp257 * tmp344 + tmp141 * tmp252 * tmp342 +
                                tmp141 * tmp352 *
                                    (-2 * tmp105 * tmp240 * tmp26 + tmp106 * tmp349 - tmp149 * tmp250 * tmp350 + tmp157 * tmp240 * tmp249 +
                                     tmp229 * tmp36 - tmp247 * tmp250 * tmp319 - 2 * tmp249 * tmp26 * tmp98) +
                                (207585478834111.0 / 16000000000000.0) * tmp144 * tmp145 + (121954868780449.0 / 500000000000000.0) * tmp144 * tmp19 +
                                tmp149 * tmp228 * tmp240 * tmp317 +
                                tmp149 * tmp317 * tmp98 *
                                    (-11612160 * r * tmp93 + tmp218 + 4 * tmp226 * tmp57 - 2 * tmp26 * tmp97 +
                                     (9073673263775509632.0 / 390625.0) * tmp27 + 2 * tmp326 + (113917534153075584.0 / 15625.0) * tmp37 +
                                     23151752773632 * tmp57 - 668771063179675776.0 / 3125.0) +
                                tmp150 * tmp223 * tmp317 *
                                    (-2 * tmp103 * tmp26 + 4 * tmp221 * tmp57 + 23151752773632 * tmp26 * tmp69 - 23151752773632 * tmp26 -
                                     258910106287381.0 / 25.0 * tmp27 - 535088332800 * tmp37 - 3519692144640 * tmp69 - 510774533137571.0 / 25.0) -
                                tmp150 * tmp319 * tmp325 - tmp152 * tmp324 * tmp336 +
                                tmp153 * tmp155 * tmp236 *
                                    (-61440 * tmp11 * tmp301 - 24 * tmp232 - tmp237 * tmp298 - tmp237 * tmp307 - tmp237 * tmp309 - 245760 * tmp299 -
                                     368640 * tmp300 - 245760 * tmp302) -
                                10 * tmp153 * tmp180 * tmp336 - 8 * tmp155 * tmp205 + tmp157 * tmp229 * tmp265 * tmp342 - tmp162 * tmp287 +
                                tmp164 * tmp209 * tmp8 * (-tmp255 - 2) - 7130030404887369.0 / 40000000000000.0 * tmp203 * tmp27 * tmp77 +
                                tmp208 * tmp209 * tmp233 * tmp30 + 30 * tmp211 * tmp296 + 2 * tmp228 * tmp317 * tmp319 * tmp98 +
                                tmp231 * tmp234 * tmp333 * tmp8 + 8 * tmp231 * tmp296 * tmp314 * tmp86 + tmp241 * tmp343 + tmp241 * tmp354 -
                                tmp243 * tmp343 - tmp243 * tmp354 - tmp249 * tmp332 * tmp98 + tmp251 * tmp343 + tmp251 * tmp354 +
                                tmp252 * tmp265 * tmp352 + tmp283 * tmp71 * tmp73 + tmp283 * tmp80 + tmp285 * tmp82 + tmp285 * tmp85 +
                                tmp286 * tmp310 + tmp286 + tmp288 * tmp70 * tmp71 - 4 * tmp31 * tmp313 * tmp8 - tmp316 * tmp323 * tmp330 +
                                tmp316 * tmp330 * tmp337 - tmp321 * tmp323 + tmp321 * tmp337 - tmp325 * tmp330 + 6 * tmp72 * tmp75 +
                                tmp72 * ((516952380952381.0 / 10000000000000.0) * m1 * m2 * tmp26 * tmp3 - tmp204 * tmp26) +
                                tmp164 * tmp208 * tmp8 * (-4 * r + 4 * tmp26 * tmp3 * tmp7) / ((tmp160) * (tmp160) * (tmp160)) +
                                tmp153 * tmp155 * tmp238 * (-153600 * tmp185 + 10 * tmp19 * tmp3 * tmp7 - 38400 * tmp345 - 38400 * tmp347) /
                                    pow(tmp131, 6) +
                                tmp150 * tmp222 * tmp317 *
                                    (-776730318862143.0 / 50.0 * r * tmp27 - 802632499200 * r * tmp37 - 1532323599412713.0 / 50.0 * r +
                                     3 * tmp103 * tmp57 - 34727629160448 * tmp188 + 3 * tmp221 * tmp69 - 87107263332245337.0 / 5000.0 * tmp27 -
                                     546804163934763.0 / 100.0 * tmp37 + 80252527581561008619.0 / 500000.0) /
                                    ((tmp104) * (tmp104) * (tmp104) * (tmp104)) +
                                104 * tmp216 / pow(r, 9) + 48 * tmp152 * tmp214 * tmp314 * tmp89 / pow(r, 11) +
                                (7298465773012227.0 / 2000000000000.0) * tmp146 / pow(r, 13.0 / 2.0))) +
            2 * tmp54 * tmp64 + tmp61 * tmp63 * (-tmp174 - 6 * tmp25) / ((tmp9) * (tmp9) * (tmp9))));
  const REAL dHreal_dr_dpphi =
      tmp363 * (-tmp276 * tmp280 * tmp375 + tmp282 * (tmp10 * (pphi * (tmp46 * (-pphi * tmp285 * tmp38 - tmp29 * tmp376 - tmp333 * tmp366 * tmp40) +
                                                                       tmp52 * (-pphi * tmp285 * tmp50 - tmp333 * tmp366 * tmp51 - tmp376 * tmp48)) -
                                                               tmp13 * tmp17 - tmp24 + tmp53) +
                                                      tmp275 * (tmp177 * tmp374 + tmp202 * tmp374 +
                                                                tmp269 * (4 * pphi * tmp12 * tmp162 - tmp310 * tmp376 + tmp313 * tmp372 - tmp376)) -
                                                      tmp293 * tmp373 + tmp293 * tmp374 + tmp370 * tmp64));
  const REAL xi = tmp131 * sqrt(r * tmp138) / tmp88;
  const REAL dHreal_dprstar =
      (1.0 / 2.0) * tmp273 * tmp282 *
      (prstar * tmp157 * tmp229 * tmp230 * tmp89 - 339525257375589.0 / 5000000000000.0 * tmp142 * tmp379 +
       (29655068404873.0 / 2500000000000.0) * tmp143 * tmp27 * tmp380 + (73721876495073.0 / 125000000000.0) * tmp145 * tmp27 * tmp378 +
       (121954868780449.0 / 125000000000000.0) * tmp147 * tmp380 + tmp153 * tmp154 * tmp378 * tmp44 + 6 * tmp19 * tmp379 * tmp79 +
       tmp254 * tmp378 * tmp73 + 4 * tmp26 * tmp378 * tmp75 + 6 * tmp26 * tmp379 * tmp81 + 8 * tmp26 * tmp380 * tmp84 + tmp368 * tmp378 * tmp70);
  const REAL Hreal = tmp281;
  const REAL Omega = tmp282 * tmp375;
  const REAL Omega_circ = (tmp10 * tmp370 + tmp374 * tmp382 / tmp165) / sqrt(2 * tmp27 * (tmp277 + tmp382) + 1);

  const REAL dLdr = -dHreal_dr_dr / dHreal_dr_dpphi;
  const REAL rdot_dyn = xi * dHreal_dprstar;
  REAL flux[2];
  const REAL y[4] = {r, 0., prstar, pphi};
  SEOBNRv5_aligned_spin_flux(y, Hreal, Omega, Omega_circ, flux, params);
  const REAL dLdt = flux[1];
  const REAL rdot_rad = dLdt / dLdr;
  const REAL prstar_condition = rdot_dyn - rdot_rad;
  return prstar_condition;
} // END FUNCTION SEOBNRv5_aligned_spin_radial_momentum_conditions
