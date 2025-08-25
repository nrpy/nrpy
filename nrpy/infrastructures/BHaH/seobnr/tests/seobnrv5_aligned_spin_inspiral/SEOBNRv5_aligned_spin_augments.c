#include "BHaH_defines.h"
/**
 * Evaluate SEOBNRv5 Hamiltonian and circular derivatives.
 */
void SEOBNRv5_aligned_spin_augments(commondata_struct *restrict commondata) {
#include "set_CodeParameters.h"
  const REAL tmp0 = ((prstar) * (prstar) * (prstar) * (prstar));
  const REAL tmp1 = ((r) * (r) * (r) * (r));
  const REAL tmp3 = ((m1) * (m1) * (m1));
  const REAL tmp4 = ((m2) * (m2) * (m2));
  const REAL tmp5 = m1 + m2;
  const REAL tmp8 = log(r);
  const REAL tmp12 = ((m1) * (m1));
  const REAL tmp13 = ((m2) * (m2));
  const REAL tmp16 = ((r) * (r) * (r));
  const REAL tmp18 = ((r) * (r));
  const REAL tmp20 = pow(prstar, 6);
  const REAL tmp22 = pow(prstar, 8);
  const REAL tmp30 = (1.0 / (r));
  const REAL tmp37 = ((r) * (r) * (r) * (r) * (r));
  const REAL tmp76 = 2 * pphi;
  const REAL tmp2 = (1.0 / (tmp1));
  const REAL tmp6 = pow(tmp5, -6);
  const REAL tmp9 = (1.0 / ((tmp5) * (tmp5)));
  const REAL tmp14 = (1.0 / ((tmp5) * (tmp5) * (tmp5) * (tmp5)));
  const REAL tmp17 = (1.0 / (tmp16));
  const REAL tmp19 = (1.0 / (tmp18));
  const REAL tmp21 = ((m1) * (m1) * (m1) * (m1)) * ((m2) * (m2) * (m2) * (m2)) / pow(tmp5, 8);
  const REAL tmp25 = chi1 * m1 + chi2 * m2;
  const REAL tmp38 = 8 * r + 2 * tmp16 + 4 * tmp18 + 16;
  const REAL tmp43 = chi1 * m1 - chi2 * m2;
  const REAL tmp46 = (1.0 / ((tmp5) * (tmp5) * (tmp5)));
  const REAL tmp50 = (1.0 / (tmp37));
  const REAL tmp53 = 2 * tmp30 + 1;
  const REAL tmp7 = tmp3 * tmp4 * tmp6;
  const REAL tmp11 = m1 * m2 * tmp9;
  const REAL tmp15 = tmp12 * tmp13 * tmp14;
  const REAL tmp32 = 5787938193408 * ((tmp8) * (tmp8));
  const REAL tmp44 = ((tmp43) * (tmp43)) * tmp9;
  const REAL tmp45 = tmp43 * (m1 - m2);
  const REAL tmp52 = ((pphi) * (pphi)) * tmp19;
  const REAL tmp64 = ((pphi) * (pphi)) * tmp17;
  const REAL tmp66 = ((pphi) * (pphi) * (pphi) * (pphi)) * tmp2;
  const REAL tmp67 = tmp25 / tmp5;
  const REAL tmp77 = tmp19 * tmp76;
  const REAL tmp79 = 4 * ((pphi) * (pphi) * (pphi)) * tmp2;
  const REAL tmp27 = ((tmp25) * (tmp25)) * tmp9;
  const REAL tmp33 = 14700 * tmp11 + 42911;
  const REAL tmp48 = tmp25 * tmp45 * tmp46;
  const REAL tmp60 = dSO * m1 * m2 * tmp17 * tmp25 * tmp46;
  const REAL tmp61 = (1.0 / 4.0) * tmp19 * (tmp14 * ((tmp25) * (tmp25)) * tmp45 - ((tmp25) * (tmp25) * (tmp25)) * tmp46);
  const REAL tmp62 = -45.0 / 32.0 * tmp11 - 15.0 / 32.0;
  const REAL tmp63 = -267.0 / 128.0 * tmp11 - 1591.0 / 768.0 * tmp15 + 59.0 / 256.0;
  const REAL tmp65 = (75.0 / 128.0) * tmp11 + (345.0 / 256.0) * tmp15 + 105.0 / 256.0;
  const REAL tmp68 = 15.0 / 32.0 - 9.0 / 32.0 * tmp11;
  const REAL tmp69 = -35.0 / 128.0 * tmp11 - 613.0 / 768.0 * tmp15 - 59.0 / 256.0;
  const REAL tmp70 = -45.0 / 128.0 * tmp11 + (75.0 / 256.0) * tmp12 * tmp13 * tmp14 - 105.0 / 256.0;
  const REAL tmp28 = tmp19 * tmp27;
  const REAL tmp36 =
      r * tmp32 +
      r * ((3390101660860523312.0 / 78125.0) * tmp11 + (3369809522332764779.0 / 78125.0) * tmp15 + (296393260946151.0 / 50.0) * tmp7 +
           (866182644304933.0 / 10.0) * ((1 - 0.49694878161693501 * tmp11) * (1 - 0.49694878161693501 * tmp11)) + 188440788778196) +
      55296 * tmp11 *
          (-31383302728516087.0 / 12500000.0 * tmp11 - 426364516032331.0 / 10000000.0 * tmp15 + 14515200 * tmp7 + 100201376401019.0 / 100000.0) -
      967680 * tmp16 * (-53511513369581.0 / 20000000.0 * tmp11 - 138240 * tmp15 - 52783413229329.0 / 10000000.0) -
      9216 * tmp18 * (-7876452485916241.0 / 12500000.0 * tmp11 - 98886748396767.0 / 500000.0 * tmp15 + 580530436787913.0 / 100000.0) +
      tmp8 * (49152 * r * ((16368307925443.0 / 40000.0) * tmp11 + 102574080 * tmp15 - 105983591868019.0 / 50000.0) -
              1698693120 * tmp11 * (11592 * tmp11 + 69847) + 879923036160 * tmp16 + 283115520 * tmp18 * tmp33);
  const REAL tmp40 = 756 * tmp11 + 1079;
  const REAL tmp54 = -tmp53 / (tmp18 + tmp27 * tmp53);
  const REAL tmp59 = (1.0 / (tmp16 + tmp27 * (r + 2)));
  const REAL tmp72 = tmp45 * tmp9 *
                         (tmp19 * (-1.0 / 32.0 * tmp11 + (103.0 / 192.0) * tmp15 + 5.0 / 64.0) + tmp30 * ((11.0 / 32.0) * tmp11 + 3.0 / 32.0) +
                          tmp52 * tmp68 + tmp64 * tmp69 + tmp66 * tmp70 + 1.0 / 4.0) +
                     tmp67 * (tmp19 * (-177.0 / 32.0 * tmp11 + (109.0 / 192.0) * tmp12 * tmp13 * tmp14 - 5.0 / 64.0) +
                              tmp30 * ((23.0 / 32.0) * tmp11 - 3.0 / 32.0) + tmp52 * tmp62 + tmp63 * tmp64 + tmp65 * tmp66 + 7.0 / 4.0);
  const REAL tmp35 =
      (53058305272831.0 / 5.0) * r * tmp11 + (182268054644921.0 / 100.0) * r * tmp15 - 163418713120743.0 / 500000.0 * r * tmp33 -
      197383821284628.0 / 5.0 * r + (258910106287381.0 / 100.0) * tmp11 * tmp18 - 60249543508726.0 / 5.0 * tmp11 + 133772083200 * tmp15 * tmp18 +
      (400296247701391.0 / 5.0) * tmp15 + (510774533137571.0 / 100.0) * tmp18 + tmp32 + (336524885906151.0 / 50.0) * tmp7 -
      tmp8 * (-283115520 * r * tmp33 - 253929149957443.0 / 10.0 * tmp11 - 5041721180160 * tmp15 - 879923036160 * tmp18 + 104186110149937) +
      275059053208689;
  const REAL tmp41 =
      7680 *
      (2048 * m1 * m2 * tmp8 * tmp9 * (336 * r + 756 * tmp11 + 407) + 28 * m1 * m2 * tmp9 * (1920 * a6 + 733955307463037.0 / 1000000000.0) -
       7 * r * ((938918400156317.0 / 1000000000.0) * m1 * m2 * tmp9 - 185763092693281.0 / 1000000000.0 * tmp15 - 245760) -
       270820329770593.0 / 50000000.0 * tmp15 - 3440640) /
      (32768 * tmp11 * tmp8 *
           (240 * r * (-373313530533103.0 / 50000000000.0 * tmp11 - 3024 * tmp15 + 17264) + 480 * tmp1 * tmp40 -
            388422414769507.0 / 10000000.0 * tmp11 - 47061405915993.0 / 25000000.0 * tmp15 + 960 * tmp16 * tmp40 +
            1920 * tmp18 * (588 * tmp11 + 1079) + 161280 * tmp37 + 13447680) +
       53760 * tmp11 *
           (7680 * a6 * (tmp1 + tmp38) + (113485217444961.0 / 1000000000.0) * r * (-tmp1 + tmp38) +
            (7402203300817.0 / 50000000000.0) * r * (7704 * r + 349 * tmp1 + 1926 * tmp16 + 3852 * tmp18 + 36400) +
            128 * r *
                ((33046962773603.0 / 2500000000.0) * r - 137046962773603.0 / 20000000000.0 * tmp1 + (42646962773603.0 / 10000000000.0) * tmp16 +
                 (852939255472061.0 / 100000000000.0) * tmp18 - 42153037226397.0 / 1250000000.0)) +
       67645734912 * tmp15 * ((tmp8) * (tmp8)) +
       7 * tmp15 *
           (-39321600 * a6 * (3 * r + 59) + (186464462028901.0 / 250000.0) * a6 + (122635399361987.0 / 1000.0) * r -
            308925070376879.0 / 100000.0 * tmp1 - 308925070376879.0 / 50000.0 * tmp16 + (206478381132587.0 / 100000.0) * tmp18 +
            (3566651379711.0 / 2500.0) * tmp37 + 276057889687011.0 / 1000.0) +
       (241555486248807.0 / 1000.0) * tmp21 + 13212057600 * tmp37 +
       1120 * tmp7 * (-163683964822551.0 / 1000000.0 * r - 3566651379711.0 / 200000.0 * tmp18 - 59449372951581.0 / 50000.0));
  const REAL tmp55 = tmp27 * tmp52 * tmp54 + tmp52 + 1;
  const REAL tmp73 = tmp59 * (pphi * tmp60 + pphi * tmp61 + pphi * tmp72) - 1;
  const REAL tmp80 = tmp59 * (pphi * (tmp45 * tmp9 * (tmp17 * tmp69 * tmp76 + tmp68 * tmp77 + tmp70 * tmp79) +
                                      tmp67 * (tmp17 * tmp63 * tmp76 + tmp62 * tmp77 + tmp65 * tmp79)) +
                              tmp60 + tmp61 + tmp72);
  const REAL tmp81 = (1.0 / 2.0) * tmp28 * tmp54 * tmp76 + (1.0 / 2.0) * tmp77;
  const REAL tmp42 = tmp1 * tmp41 + tmp28;
  const REAL tmp56 =
      ((prstar) * (prstar)) * tmp30 * tmp36 * ((tmp28 + 1) * (tmp28 + 1)) *
          (tmp17 * (tmp27 * (3 * tmp11 + 45.0 / 16.0) + tmp44 * ((3.0 / 4.0) * tmp11 - 3.0 / 16.0) - 21.0 / 8.0 * tmp48) +
           tmp2 * (tmp27 * (-1171.0 / 64.0 * tmp11 - 861.0 / 64.0) + tmp44 * ((115.0 / 64.0) * tmp11 + (1.0 / 16.0) * tmp15 - 37.0 / 64.0) +
                   tmp48 * ((13.0 / 16.0) * tmp11 + 449.0 / 32.0)) +
           tmp28 + tmp35 * tmp37 * tmp41 / tmp36) /
          (tmp35 * ((tmp42) * (tmp42))) +
      tmp0 * tmp17 * ((115888805356193.0 / 1250000000000.0) * tmp11 - 131 * tmp15 + 10 * tmp7) + tmp0 * tmp19 * (8 * m1 * m2 * tmp9 - 6 * tmp15) +
      tmp0 * tmp2 *
          (tmp11 * (452542166996693.0 / 1000000000000.0 - 516952380952381.0 / 10000000000000.0 * tmp8) +
           tmp15 * (118.40000000000001 * tmp8 - 179613660498019.0 / 100000000000.0) + (150579635104141.0 / 250000000000.0) * tmp7) +
      tmp0 * ((tmp36) * (tmp36)) * tmp50 * ((tmp28 + 1) * (tmp28 + 1) * (tmp28 + 1) * (tmp28 + 1)) *
          (tmp27 * ((165.0 / 32.0) * tmp11 - 5 * tmp15 + 25.0 / 32.0) + tmp44 * ((75.0 / 32.0) * m1 * m2 * tmp9 - 15.0 / 8.0 * tmp15 - 15.0 / 32.0) +
           tmp48 * ((45.0 / 8.0) * tmp11 - 5.0 / 16.0)) /
          (((tmp35) * (tmp35)) * ((tmp42) * (tmp42) * (tmp42) * (tmp42))) +
      (121954868780449.0 / 1000000000000000.0) * tmp11 * tmp22 * tmp30 +
      tmp17 * tmp20 *
          (-84945530542609.0 / 2500000000000.0 * tmp11 - 447649163680617.0 / 5000000000000.0 * tmp15 - 14 * tmp21 + 188 * tmp3 * tmp4 * tmp6) +
      tmp19 * tmp20 * (-139150381847503.0 / 50000000000000.0 * tmp11 - 27.0 / 5.0 * tmp15 + 6 * tmp3 * tmp4 * tmp6) +
      tmp19 * tmp22 *
          ((4343054718629.0 / 3125000000000.0) * tmp11 + (166921011824161.0 / 50000000000000.0) * tmp15 - 6 * tmp21 +
           (342857142857143.0 / 100000000000000.0) * tmp7) +
      tmp55 + (29655068404873.0 / 20000000000000.0) * tmp11 * tmp22 / pow(r, 5.0 / 2.0) -
      113175085791863.0 / 10000000000000.0 * tmp11 * tmp20 / pow(r, 7.0 / 2.0) +
      (73721876495073.0 / 500000000000.0) * tmp0 * tmp11 / pow(r, 9.0 / 2.0);
  const REAL tmp57 =
      (tmp2 * ((9.0 / 8.0) * tmp27 + tmp44 * ((1.0 / 2.0) * tmp11 + 1.0 / 8.0) - 5.0 / 4.0 * tmp48) + tmp42 +
       tmp50 * (tmp27 * (-175.0 / 64.0 * tmp11 - 225.0 / 64.0) + tmp44 * (-81.0 / 64.0 * tmp11 + (21.0 / 16.0) * tmp12 * tmp13 * tmp14 - 9.0 / 64.0) +
                tmp48 * (117.0 / 32.0 - 39.0 / 16.0 * tmp11))) /
      (tmp28 * tmp53 + 1);
  const REAL tmp58 = sqrt(tmp56 * tmp57);
  const REAL tmp82 = sqrt(tmp55 * tmp57);
  const REAL tmp75 = sqrt(2 * tmp11 * (tmp58 + tmp73) + 1);
  commondata->Hreal = tmp75;
  commondata->dHreal_dpphi = (tmp80 + tmp58 * tmp81 / tmp56) / tmp75;
  commondata->Omega_circ = (tmp80 + tmp81 * tmp82 / tmp55) / sqrt(2 * tmp11 * (tmp73 + tmp82) + 1);
} // END FUNCTION SEOBNRv5_aligned_spin_augments
