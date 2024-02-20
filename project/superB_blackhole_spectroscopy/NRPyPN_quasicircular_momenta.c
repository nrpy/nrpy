#include "BHaH_defines.h"
/*
 * Compute quasicircular momenta using validated expressions from NRPyPN.
 */
void NRPyPN_quasicircular_momenta(commondata_struct *restrict commondata) {
  // compute quasicircular parameters if commondata.p_t and commondata.p_r not
  // already set in the parfile (i.e., reset from their default values of -1.0).
  if (commondata->initial_p_t == -1.0 && commondata->initial_p_r == -1.0) {
    // In NRPyPN, q = m2/m1 = mass_M/mass_m. So mass_m = object 1, and mass_M is object 2.
    const REAL q = commondata->mass_ratio;
    const REAL r = commondata->initial_sep;
    const REAL chi1U0 = commondata->bbhxy_BH_M_chix;
    const REAL chi1U1 = commondata->bbhxy_BH_M_chiy;
    const REAL chi1U2 = commondata->bbhxy_BH_M_chiz;
    const REAL chi2U0 = commondata->bbhxy_BH_m_chix;
    const REAL chi2U1 = commondata->bbhxy_BH_m_chiy;
    const REAL chi2U2 = commondata->bbhxy_BH_m_chiz;

    const REAL mass_M = q / (1.0 + q);
    const REAL mass_m = 1.0 / (1.0 + q);
    // In NRPyPN, q = m2/m1 = mass_M/mass_m. So mass_m = object 1, and mass_M is object 2.
    const REAL m1 = mass_m;
    const REAL S1U0 = chi1U0 * mass_m * mass_m;
    const REAL S1U1 = chi1U1 * mass_m * mass_m;
    const REAL S1U2 = chi1U2 * mass_m * mass_m;
    const REAL m2 = mass_M;
    const REAL S2U0 = chi2U0 * mass_M * mass_M;
    const REAL S2U1 = chi2U1 * mass_M * mass_M;
    const REAL S2U2 = chi2U2 * mass_M * mass_M;

    REAL Pt, Pr;
    {
      const REAL tmp0 = m2 / m1;
      const REAL tmp5 = ((m2) * (m2)) / ((m1) * (m1));
      const REAL tmp6 = ((m2) * (m2) * (m2)) / ((m1) * (m1) * (m1));
      const REAL tmp10 = ((m2) * (m2) * (m2) * (m2)) / ((m1) * (m1) * (m1) * (m1));
      const REAL tmp12 = ((chi1U0) * (chi1U0));
      const REAL tmp14 = ((chi1U1) * (chi1U1));
      const REAL tmp16 = ((chi1U2) * (chi1U2));
      const REAL tmp22 = ((chi2U0) * (chi2U0));
      const REAL tmp24 = ((chi2U1) * (chi2U1));
      const REAL tmp26 = ((chi2U2) * (chi2U2));
      const REAL tmp31 = ((M_PI) * (M_PI));
      const REAL tmp45 = ((m2) * (m2) * (m2) * (m2) * (m2)) / ((m1) * (m1) * (m1) * (m1) * (m1));
      const REAL tmp1 = tmp0 + 1;
      const REAL tmp38 = tmp6 * (20 * tmp0 + 5);
      const REAL tmp40 = (15.0 / 4.0) * tmp5;
      const REAL tmp2 = (1.0 / ((tmp1) * (tmp1)));
      const REAL tmp8 = (1.0 / ((tmp1) * (tmp1) * (tmp1) * (tmp1)));
      const REAL tmp33 = (1.0 / ((tmp1) * (tmp1) * (tmp1)));
      const REAL tmp36 = 5 * tmp0 + 20;
      const REAL tmp46 = (1.0 / 32.0) / pow(tmp1, 6);
      const REAL tmp9 = chi1U2 * tmp8;
      const REAL tmp11 = chi2U2 * tmp8;
      const REAL tmp13 = (3.0 / 2.0) * tmp2;
      const REAL tmp15 = (3.0 / 4.0) * tmp2;
      const REAL tmp28 = (1.0 / 128.0) * tmp8;
      const REAL tmp39 = (1.0 / 8.0) * tmp38 * tmp8;
      const REAL tmp42 = (15.0 / 4.0) * tmp0 * tmp11;
      Pt = tmp0 * tmp2 *
           (1 + 2 / r +
            (-3 * chi1U0 * chi2U0 * tmp0 * tmp2 + chi1U1 * chi2U1 * tmp0 * tmp13 + chi1U2 * chi2U2 * tmp0 * tmp13 - tmp12 * tmp13 -
             tmp13 * tmp22 * tmp5 + tmp14 * tmp15 + tmp15 * tmp16 + tmp15 * tmp24 * tmp5 + tmp15 * tmp26 * tmp5 +
             (1.0 / 16.0) * tmp2 * (41 * tmp0 + 42 * tmp5 + 42)) /
                ((r) * (r)) +
            (chi1U0 * chi2U0 * tmp28 * (192 * tmp0 + 560 * tmp5 + 192 * tmp6) + chi1U1 * chi2U1 * tmp28 * (-864 * tmp0 - 1856 * tmp5 - 864 * tmp6) +
             chi1U2 * chi2U2 * tmp28 * (480 * tmp0 + 1064 * tmp5 + 480 * tmp6) + tmp12 * tmp28 * (472 * tmp5 - 640) +
             tmp14 * tmp28 * (-640 * tmp0 - 512 * tmp5 - 64) + tmp16 * tmp28 * (224 * tmp0 - 108 * tmp5 + 512) +
             tmp22 * tmp28 * (-640 * tmp10 + 472 * tmp5) + tmp24 * tmp28 * (-64 * tmp10 - 512 * tmp5 - 640 * tmp6) +
             tmp26 * tmp28 * (512 * tmp10 - 108 * tmp5 + 224 * tmp6) +
             tmp28 * (163 * tmp0 * tmp31 - 2636 * tmp0 + 480 * tmp10 + 326 * tmp31 * tmp5 + 163 * tmp31 * tmp6 - 6128 * tmp5 - 2636 * tmp6 + 480)) /
                ((r) * (r) * (r)) +
            (-1.0 / 4.0 * chi1U2 * tmp2 * (9 * tmp0 + 12) + (1.0 / 4.0) * chi2U2 * tmp2 * (-9 * tmp0 - 12 * tmp5)) / pow(r, 3.0 / 2.0) +
            ((1.0 / 16.0) * tmp11 * (-13 * tmp0 - 72 * tmp10 - 60 * tmp5 - 116 * tmp6) +
             (1.0 / 16.0) * tmp9 * (-116 * tmp0 - 60 * tmp5 - 13 * tmp6 - 72)) /
                pow(r, 5.0 / 2.0) +
            (chi1U0 * ((1.0 / 4.0) * chi2U0 * tmp0 * tmp9 * (15 * tmp0 + 30) + (1.0 / 4.0) * chi2U0 * tmp11 * tmp5 * (30 * tmp0 + 15)) +
             chi1U1 * (chi2U1 * tmp11 * tmp40 + chi2U1 * tmp40 * tmp9) - 1.0 / 8.0 * ((chi1U2) * (chi1U2) * (chi1U2)) * tmp36 * tmp8 +
             chi1U2 * ((15.0 / 4.0) * tmp22 * tmp5 * tmp8 * (2 * tmp0 + 3) - tmp24 * tmp40 * tmp8 * (tmp0 + 2) - tmp26 * tmp33 * tmp40 -
                       tmp46 * (670 * tmp0 + 145 * tmp10 + 103 * tmp45 + 252 * tmp5 - 27 * tmp6 + 348)) -
             ((chi2U2) * (chi2U2) * (chi2U2)) * tmp39 - 15.0 / 4.0 * chi2U2 * tmp0 * tmp16 * tmp33 -
             chi2U2 * tmp0 * tmp46 * (145 * tmp0 + 670 * tmp10 + 348 * tmp45 - 27 * tmp5 + 252 * tmp6 + 103) - chi2U2 * tmp24 * tmp39 +
             (1.0 / 2.0) * tmp11 * tmp22 * tmp38 + tmp12 * ((1.0 / 2.0) * tmp36 * tmp9 + tmp42 * (3 * tmp0 + 2)) +
             tmp14 * (-tmp42 * (2 * tmp0 + 1) + (1.0 / 8.0) * tmp9 * (-5 * tmp0 - 20))) /
                pow(r, 7.0 / 2.0)) /
           sqrt(r);
    }
    {
      const REAL tmp0 = (1.0 / (m1));
      const REAL tmp4 = (1.0 / (m2));
      const REAL tmp7 = ((m1) * (m1));
      const REAL tmp9 = ((m2) * (m2));
      const REAL tmp11 = (1.0 / (r));
      const REAL tmp14 = (1.0 / ((m1) * (m1) * (m1) * (m1)));
      const REAL tmp15 = ((m2) * (m2) * (m2) * (m2));
      const REAL tmp17 = ((m1) * (m1) * (m1));
      const REAL tmp19 = ((m2) * (m2) * (m2));
      const REAL tmp22 = (1.0 / ((r) * (r)));
      const REAL tmp24 = pow(r, -5.0 / 2.0);
      const REAL tmp31 = pow(r, -7.0 / 2.0);
      const REAL tmp32 = (1.0 / ((m1) * (m1) * (m1) * (m1) * (m1)));
      const REAL tmp33 = ((m2) * (m2) * (m2) * (m2) * (m2));
      const REAL tmp40 = (1.0 / ((r) * (r) * (r)));
      const REAL tmp41 = ((M_PI) * (M_PI));
      const REAL tmp42 = ((chi1U1) * (chi1U1));
      const REAL tmp44 = ((chi1U2) * (chi1U2));
      const REAL tmp45 = ((chi2U1) * (chi2U1));
      const REAL tmp48 = ((chi2U2) * (chi2U2));
      const REAL tmp52 = chi1U0 * chi2U0;
      const REAL tmp56 = chi1U1 * chi2U1;
      const REAL tmp57 = chi1U2 * chi2U2;
      const REAL tmp58 = ((chi1U0) * (chi1U0));
      const REAL tmp60 = ((chi2U0) * (chi2U0));
      const REAL tmp61 = pow(m1, -6);
      const REAL tmp62 = pow(m2, 6);
      const REAL tmp71 = (1.0 / ((r) * (r) * (r) * (r)));
      const REAL tmp78 = pow(r, -3.0 / 2.0);
      const REAL tmp122 = m1 + m2;
      const REAL tmp124 = m1 * m2;
      const REAL tmp127 = ((Pt) * (Pt));
      const REAL tmp128 = S1U2 * S2U2;
      const REAL tmp133 = 2 * Pt;
      const REAL tmp134 = (3.0 / 2.0) * Pt;
      const REAL tmp136 = S1U0 * S2U0;
      const REAL tmp138 = S1U1 * S2U1;
      const REAL tmp139 = 4 / ((r) * (r) * (r) * (r) * (r));
      const REAL tmp146 = ((S2U0) * (S2U0));
      const REAL tmp147 = ((S2U1) * (S2U1));
      const REAL tmp151 = ((S1U0) * (S1U0));
      const REAL tmp152 = ((S1U1) * (S1U1));
      const REAL tmp164 = ((Pt) * (Pt) * (Pt));
      const REAL tmp191 = ((Pt) * (Pt) * (Pt) * (Pt));
      const REAL tmp197 = S1U1 * S2U2 - S1U2 * S2U1;
      const REAL tmp199 = S1U0 * S2U2;
      const REAL tmp201 = S1U2 * S2U0;
      const REAL tmp266 = (9.0 / 2.0) * Pt;
      const REAL tmp303 = m1 - m2;
      const REAL tmp1 = m2 * tmp0;
      const REAL tmp5 = m1 * tmp4;
      const REAL tmp8 = (1.0 / (tmp7));
      const REAL tmp16 = tmp14 * tmp15;
      const REAL tmp18 = (1.0 / (tmp17));
      const REAL tmp34 = tmp32 * tmp33;
      const REAL tmp123 = (1.0 / (tmp122));
      const REAL tmp131 = (1.0 / (tmp9));
      const REAL tmp135 = 2 * tmp40;
      const REAL tmp149 = ((S2U2) * (S2U2)) + tmp147;
      const REAL tmp154 = ((S1U2) * (S1U2)) + tmp152;
      const REAL tmp158 = tmp0 * tmp4;
      const REAL tmp163 = 3 * tmp71;
      const REAL tmp166 = (3.0 / 4.0) * tmp127;
      const REAL tmp176 = (1.0 / (tmp19));
      const REAL tmp188 = ((tmp122) * (tmp122) * (tmp122) * (tmp122));
      const REAL tmp192 = (1.0 / (tmp15));
      const REAL tmp216 = pow(tmp122, 6);
      const REAL tmp229 = tmp133 * tmp40;
      const REAL tmp234 = (75.0 / 8.0) * tmp0 * tmp19;
      const REAL tmp238 = (3.0 / 16.0) * tmp191;
      const REAL tmp241 = (159.0 / 16.0) * tmp0 + (37.0 / 8.0) * tmp4;
      const REAL tmp249 = (75.0 / 8.0) * tmp17 * tmp4;
      const REAL tmp250 = m1 / tmp33;
      const REAL tmp252 = (37.0 / 8.0) * tmp0 + (159.0 / 16.0) * tmp4;
      const REAL tmp257 = Pt * tmp40;
      const REAL tmp270 = (3.0 / 2.0) * tmp199 - 3.0 / 2.0 * tmp201;
      const REAL tmp272 = -2 * tmp199 + 2 * tmp201;
      const REAL tmp273 = Pt * tmp22;
      const REAL tmp304 = tmp303 * (-S1U2 * tmp0 * tmp122 + S2U2 * tmp122 * tmp4) / ((tmp122) * (tmp122) * (tmp122));
      const REAL tmp2 = tmp1 + 1;
      const REAL tmp10 = tmp8 * tmp9;
      const REAL tmp20 = tmp18 * tmp19;
      const REAL tmp30 = (1.0 / 4.0) * tmp5;
      const REAL tmp49 = 5 * tmp1;
      const REAL tmp59 = 9 * tmp1;
      const REAL tmp69 = (1.0 / 4.0) * tmp1;
      const REAL tmp72 = (3.0 / 2.0) * tmp1;
      const REAL tmp75 = 2 * tmp1;
      const REAL tmp79 = 3 * tmp1;
      const REAL tmp125 = tmp123 * tmp124;
      const REAL tmp130 = (9.0 / 2.0) * tmp127 * tmp128;
      const REAL tmp141 = (1.0 / ((tmp122) * (tmp122)));
      const REAL tmp143 = tmp131 * tmp8;
      const REAL tmp150 = tmp146 + tmp149;
      const REAL tmp155 = tmp151 + tmp154;
      const REAL tmp157 = (1.0 / 2.0) * tmp128 + (1.0 / 2.0) * tmp136 + (1.0 / 2.0) * tmp138;
      const REAL tmp162 = -5 * tmp0 * tmp9;
      const REAL tmp167 = tmp158 * tmp166;
      const REAL tmp170 = m2 * tmp18;
      const REAL tmp174 = -5 * tmp4 * tmp7;
      const REAL tmp177 = m1 * tmp176;
      const REAL tmp179 = tmp5 + 1;
      const REAL tmp189 = tmp7 * tmp9 / tmp188;
      const REAL tmp193 = tmp14 * tmp188 * tmp192;
      const REAL tmp196 = (15.0 / 2.0) * S1U2 * tmp146;
      const REAL tmp203 = (3.0 / 2.0) * Pt * tmp199 - 3.0 / 2.0 * Pt * tmp201;
      const REAL tmp205 = S2U2 * (-4 * tmp146 + tmp149);
      const REAL tmp206 = 5 * S1U2 * tmp146;
      const REAL tmp207 = -2 * Pt * tmp199 + 2 * Pt * tmp201;
      const REAL tmp209 = (15.0 / 2.0) * S2U2 * tmp151;
      const REAL tmp211 = S1U2 * (-4 * tmp151 + tmp154);
      const REAL tmp212 = 5 * S2U2 * tmp151;
      const REAL tmp213 = (3.0 / 2.0) * tmp5;
      const REAL tmp217 = tmp17 * tmp19 / tmp216;
      const REAL tmp219 = tmp216 * tmp61 / tmp62;
      const REAL tmp222 = (47.0 / 8.0) * tmp0 + 5 * tmp4;
      const REAL tmp224 = Pt * tmp163;
      const REAL tmp225 = tmp18 * tmp4;
      const REAL tmp227 = tmp0 * tmp176;
      const REAL tmp246 = 5 * tmp0 + (47.0 / 8.0) * tmp4;
      const REAL tmp258 = 3 * tmp128 * tmp257;
      const REAL tmp263 = tmp134 * tmp158;
      const REAL tmp306 = (1.0 / 2.0) * tmp8;
      const REAL tmp308 = (1.0 / 2.0) * tmp131;
      const REAL tmp314 = tmp123 * tmp303;
      const REAL tmp26 = (1.0 / (4 * tmp1 + 4));
      const REAL tmp29 = (1.0 / (tmp2));
      const REAL tmp35 = ((tmp2) * (tmp2) * (tmp2) * (tmp2));
      const REAL tmp50 = 3 * tmp10 + 3;
      const REAL tmp87 = tmp10 * tmp60;
      const REAL tmp95 = pow(tmp2, -6);
      const REAL tmp97 = (1.0 / ((tmp2) * (tmp2) * (tmp2)));
      const REAL tmp101 = tmp20 * (12 * tmp1 + 3);
      const REAL tmp106 = (9.0 / 4.0) * tmp10;
      const REAL tmp142 = tmp124 * tmp141;
      const REAL tmp144 = ((tmp122) * (tmp122)) * tmp143;
      const REAL tmp171 = -5.0 / 8.0 * tmp127 * tmp170 + tmp166 * tmp8;
      const REAL tmp178 = -5.0 / 8.0 * tmp127 * tmp177 + tmp131 * tmp166;
      const REAL tmp184 = (9.0 / 4.0) * tmp127 * tmp158;
      const REAL tmp204 = S1U2 * tmp150;
      const REAL tmp210 = S2U2 * tmp155;
      const REAL tmp223 = (13.0 / 2.0) * tmp0 * tmp127 + tmp127 * tmp222 + (53.0 / 8.0) * tmp127 * tmp4;
      const REAL tmp228 = (1.0 / 4.0) * tmp143 * tmp191 + (1.0 / 2.0) * tmp191 * tmp225 + (1.0 / 2.0) * tmp191 * tmp227;
      const REAL tmp239 = tmp143 * tmp238 - tmp225 * tmp238 - tmp227 * tmp238;
      const REAL tmp243 = (9.0 / 2.0) * m2 * tmp8 + (23.0 / 4.0) * tmp0;
      const REAL tmp244 = -3.0 / 2.0 * m2 * tmp8 + (27.0 / 8.0) * tmp18 * tmp9;
      const REAL tmp247 = (53.0 / 8.0) * tmp0 * tmp127 + tmp127 * tmp246 + (13.0 / 2.0) * tmp127 * tmp4;
      const REAL tmp254 = (9.0 / 2.0) * m1 * tmp131 + (23.0 / 4.0) * tmp4;
      const REAL tmp255 = -3.0 / 2.0 * m1 * tmp131 + (27.0 / 8.0) * tmp176 * tmp7;
      const REAL tmp264 = (15.0 / 4.0) * tmp127 * tmp158;
      const REAL tmp274 = ((Pt) * (Pt) * (Pt) * (Pt) * (Pt)) * tmp219;
      const REAL tmp296 = tmp20 * (20 * tmp1 + 5);
      const REAL tmp297 = (15.0 / 4.0) * tmp10;
      const REAL tmp302 = tmp141 * (S1U2 + S2U2);
      const REAL tmp13 = (1.0 / ((tmp2) * (tmp2)));
      const REAL tmp36 = (1.0 / (tmp35));
      const REAL tmp43 = (9.0 / 4.0) * tmp29 * tmp5;
      const REAL tmp47 = 9 * tmp10 * tmp26;
      const REAL tmp98 = chi2U2 * tmp44 * tmp97;
      const REAL tmp115 = (1.0 / 32.0) * tmp95;
      const REAL tmp145 = tmp144 * (tmp142 + 3);
      const REAL tmp187 = tmp144 * (8 * tmp142 + 5);
      const REAL tmp190 = -20 * tmp142 - 3 * tmp189 + 5;
      const REAL tmp214 = (17.0 / 2.0) * tmp142 + (109.0 / 16.0) * tmp189 - 27.0 / 16.0;
      const REAL tmp215 = tmp144 * (tmp142 * ((1.0 / 64.0) * tmp41 - 335.0 / 48.0) - 23.0 / 8.0 * tmp189 - 25.0 / 8.0);
      const REAL tmp218 = 42 * tmp142 - 53 * tmp189 - 5 * tmp217 - 7;
      const REAL tmp240 = (7.0 / 16.0) * m2 * tmp191 * tmp32 - 9.0 / 16.0 * tmp14 * tmp191 + tmp239;
      const REAL tmp245 = -tmp127 * tmp241 - tmp127 * tmp243 + tmp127 * tmp244;
      const REAL tmp251 = -9.0 / 16.0 * tmp191 * tmp192 + (7.0 / 16.0) * tmp191 * tmp250 + tmp239;
      const REAL tmp256 = -tmp127 * tmp252 - tmp127 * tmp254 + tmp127 * tmp255;
      const REAL tmp260 = 4 * tmp164 * tmp193;
      const REAL tmp280 = (3.0 / 4.0) * tmp143 * tmp164 - 3.0 / 4.0 * tmp164 * tmp225 - 3.0 / 4.0 * tmp164 * tmp227;
      const REAL tmp284 = -tmp22 * tmp228 - tmp273 * (tmp143 * tmp164 + 2 * tmp164 * tmp225 + 2 * tmp164 * tmp227);
      const REAL tmp305 = M_PI * tmp142;
      const REAL tmp310 = S1U2 * tmp306 + S2U2 * tmp308;
      const REAL tmp311 = S1U2 * tmp306 - S2U2 * tmp308;
      const REAL tmp318 = S1U0 * tmp306 - S2U0 * tmp308;
      const REAL tmp321 = S1U1 * tmp306 - S2U1 * tmp308;
      const REAL tmp322 = S1U0 * tmp306 + S2U0 * tmp308;
      const REAL tmp323 = S1U1 * tmp306 + S2U1 * tmp308;
      const REAL tmp37 = (1.0 / 16.0) * tmp36;
      const REAL tmp54 = (1.0 / 2.0) * tmp13;
      const REAL tmp73 = tmp36 * (tmp49 + 2);
      const REAL tmp83 = chi2U2 * tmp36;
      const REAL tmp85 = (3.0 / 2.0) * tmp13;
      const REAL tmp86 = (3.0 / 4.0) * tmp13;
      const REAL tmp93 = (1.0 / 8.0) * tmp36;
      const REAL tmp105 = chi1U2 * tmp36;
      const REAL tmp116 = chi2U2 * tmp1 * tmp115;
      const REAL tmp118 = tmp36 * tmp45 * (tmp1 + 2);
      const REAL tmp119 = tmp36 * tmp87 * (tmp75 + 3);
      const REAL tmp290 = (1.0 / 128.0) * tmp36;
      const REAL tmp312 = ((tmp311) * (tmp311));
      const REAL tmp313 = ((tmp310) * (tmp310));
      const REAL tmp38 = chi2U2 * tmp37;
      const REAL tmp39 = chi1U2 * tmp37;
      const REAL tmp55 = tmp54 * (8 * tmp1 + tmp50);
      const REAL tmp68 = chi2U0 * tmp10 * tmp36;
      const REAL tmp74 = chi1U0 * tmp72 * tmp73;
      const REAL tmp82 = (1.0 / 4.0) * chi1U2 * tmp13;
      const REAL tmp100 = ((chi1U2) * (chi1U2) * (chi1U2)) * tmp93;
      const REAL tmp102 = ((chi2U2) * (chi2U2) * (chi2U2)) * tmp93;
      const REAL tmp103 = (1.0 / 2.0) * tmp60 * tmp83;
      const REAL tmp110 = (9.0 / 4.0) * tmp1 * tmp83;
      const REAL tmp113 = chi2U0 * tmp105 * tmp69;
      const REAL tmp294 =
          tmp290 * tmp42 * (-640 * tmp1 - 512 * tmp10 - 64) + tmp290 * tmp44 * (224 * tmp1 - 108 * tmp10 + 512) +
          tmp290 * tmp45 * (-512 * tmp10 - 64 * tmp16 - 640 * tmp20) + tmp290 * tmp48 * (-108 * tmp10 + 512 * tmp16 + 224 * tmp20) +
          tmp290 * tmp52 * (192 * tmp1 + 560 * tmp10 + 192 * tmp20) + tmp290 * tmp56 * (-864 * tmp1 - 1856 * tmp10 - 864 * tmp20) +
          tmp290 * tmp57 * (480 * tmp1 + 1064 * tmp10 + 480 * tmp20) + tmp290 * tmp58 * (472 * tmp10 - 640) +
          tmp290 * tmp60 * (472 * tmp10 - 640 * tmp16) +
          tmp290 * (163 * tmp1 * tmp41 - 2636 * tmp1 + 326 * tmp10 * tmp41 - 6128 * tmp10 + 480 * tmp16 + 163 * tmp20 * tmp41 - 2636 * tmp20 + 480);
      const REAL tmp77 = (3.0 / 2.0) * tmp68 * (tmp75 + 5);
      const REAL tmp90 = tmp1 * tmp56 * tmp85 + tmp1 * tmp57 * tmp85 + tmp10 * tmp45 * tmp86 + tmp10 * tmp48 * tmp86 - tmp13 * tmp52 * tmp79 +
                         tmp42 * tmp86 + tmp44 * tmp86 - tmp58 * tmp85 - tmp85 * tmp87;
      const REAL tmp104 = chi2U2 * tmp45 * tmp93;
      const REAL tmp114 = (1.0 / 4.0) * chi2U2 * tmp68;
      const REAL tmp285 = (1.0 / 4.0) * chi2U2 * tmp13 * (-12 * tmp10 - tmp59) - tmp82 * (tmp59 + 12);
      const REAL tmp287 = tmp38 * (-13 * tmp1 - 60 * tmp10 - 72 * tmp16 - 116 * tmp20) + tmp39 * (-116 * tmp1 - 60 * tmp10 - 13 * tmp20 - 72);
      const REAL tmp299 = (15.0 / 4.0) * tmp1 * tmp83;
      const REAL tmp120 =
          -tmp11 * tmp54 * (tmp49 + tmp50) + tmp22 * (tmp37 * (103 * tmp1 + 164 * tmp10 + 24 * tmp16 + 103 * tmp20 + 24) + tmp90) +
          tmp24 * ((3.0 / 16.0) * tmp1 * tmp83 * (34 * tmp1 + 30 * tmp10 + 16 * tmp20 + 13) + tmp39 * (90 * tmp1 + 102 * tmp10 + 39 * tmp20 + 48)) +
          tmp31 * (chi1U0 * (tmp113 * (tmp59 + 18) + tmp114 * (18 * tmp1 + 9)) + chi1U1 * (chi2U1 * tmp105 * tmp106 + chi2U1 * tmp106 * tmp83) +
                   chi1U2 * (-tmp106 * tmp118 - tmp106 * tmp48 * tmp97 -
                             tmp115 * (387 * tmp1 + 377 * tmp10 + 385 * tmp16 + 363 * tmp20 + 135 * tmp34 + 168) + (9.0 / 4.0) * tmp119) -
                   9.0 / 4.0 * tmp1 * tmp98 - tmp100 * (tmp79 + 12) - tmp101 * tmp102 + tmp101 * tmp103 - tmp101 * tmp104 -
                   tmp116 * (385 * tmp1 + 363 * tmp10 + 387 * tmp16 + 377 * tmp20 + 168 * tmp34 + 135) +
                   tmp42 * (chi1U2 * tmp93 * (-tmp79 - 12) - tmp110 * (tmp75 + 1)) +
                   tmp58 * ((1.0 / 2.0) * tmp105 * (tmp79 + 12) + tmp110 * (tmp79 + 2))) +
          tmp40 * (-chi1U2 * tmp1 * tmp38 * (127 * tmp1 + 96 * tmp10 + 96) + tmp1 * tmp52 * tmp93 * (187 * tmp1 + 120 * tmp10 + 120) -
                   1.0 / 32.0 * tmp10 * tmp36 * tmp48 * (14 * tmp1 + 27) * (tmp75 + 5) - tmp10 * tmp45 * tmp93 * (85 * tmp1 + 43 * tmp10 + 55) -
                   tmp36 * tmp56 * tmp69 * (95 * tmp1 + 54 * tmp10 + 54) + tmp37 * tmp58 * (180 * tmp1 + 155 * tmp10 + 76) +
                   tmp37 * tmp87 * (180 * tmp1 + 76 * tmp10 + 155) - tmp42 * tmp93 * (85 * tmp1 + 55 * tmp10 + 43) -
                   1.0 / 32.0 * tmp44 * tmp73 * (27 * tmp1 + 14) +
                   (1.0 / 384.0) * tmp95 *
                       (501 * tmp1 * tmp35 * tmp41 - 10976 * tmp1 - 40196 * tmp10 - 40196 * tmp16 - 59280 * tmp20 - 10976 * tmp34 -
                        480 * tmp61 * tmp62 - 480)) +
          tmp78 * (-chi2U2 * tmp13 * tmp69 * (4 * tmp1 + 3) + tmp82 * (-tmp79 - 4)) + 1;
      const REAL tmp288 = (1.0 / 16.0) * tmp13 * (41 * tmp1 + 42 * tmp10 + 42) + tmp90;
      const REAL tmp300 =
          chi1U0 * (tmp113 * (15 * tmp1 + 30) + tmp114 * (30 * tmp1 + 15)) + chi1U1 * (chi2U1 * tmp105 * tmp297 + chi2U1 * tmp297 * tmp83) +
          chi1U2 * (-tmp115 * (670 * tmp1 + 252 * tmp10 + 145 * tmp16 - 27 * tmp20 + 103 * tmp34 + 348) - tmp118 * tmp297 + (15.0 / 4.0) * tmp119 -
                    tmp297 * tmp48 * tmp97) -
          15.0 / 4.0 * tmp1 * tmp98 - tmp100 * (tmp49 + 20) - tmp102 * tmp296 + tmp103 * tmp296 - tmp104 * tmp296 -
          tmp116 * (145 * tmp1 - 27 * tmp10 + 670 * tmp16 + 252 * tmp20 + 348 * tmp34 + 103) +
          tmp42 * (chi1U2 * tmp93 * (-tmp49 - 20) - tmp299 * (tmp75 + 1)) + tmp58 * ((1.0 / 2.0) * tmp105 * (tmp49 + 20) + tmp299 * (tmp79 + 2));
      const REAL tmp121 = tmp120 * tmp78;
      const REAL tmp301 = pow(tmp121, 2.0 / 3.0);
      const double EulerGamma = 0.57721566490153286;
      Pr =
          (-32.0 / 5.0 * pow(tmp121, 10.0 / 3.0) * tmp189 *
               (((tmp120) * (tmp120)) * tmp40 *
                    (tmp142 * ((41.0 / 48.0) * tmp41 - 134543.0 / 7776.0) - 94403.0 / 3024.0 * tmp189 - 775.0 / 324.0 * tmp217 - 16 * M_PI * tmp302 -
                     31.0 / 6.0 * M_PI * tmp304 + (16.0 / 3.0) * tmp41 - 856.0 / 105.0 * log(16 * tmp301) - 1712.0 / 105.0 * EulerGamma +
                     6643739519.0 / 69854400.0) +
                pow(tmp121, 8.0 / 3.0) *
                    (tmp302 * ((13879.0 / 72.0) * tmp305 - 3485.0 / 96.0 * M_PI) + tmp304 * ((130583.0 / 2016.0) * tmp305 - 7163.0 / 672.0 * M_PI)) +
                pow(tmp121, 7.0 / 3.0) * (tmp302 * ((6172.0 / 189.0) * tmp142 - 2810.0 / 27.0 * tmp189 + 476645.0 / 6804.0) +
                                          tmp304 * ((1849.0 / 126.0) * tmp142 - 1501.0 / 36.0 * tmp189 + 9535.0 / 336.0) +
                                          M_PI * ((214745.0 / 1728.0) * tmp142 + (193385.0 / 3024.0) * tmp189 - 16285.0 / 504.0)) +
                pow(tmp121, 5.0 / 3.0) * (tmp142 * ((272.0 / 9.0) * tmp302 + (43.0 / 4.0) * tmp304 - 583.0 / 24.0 * M_PI) - 9.0 / 2.0 * tmp302 -
                                          13.0 / 16.0 * tmp304 - 1.0 / 4.0 * tmp310 * (1 - 3 * tmp142) * (9 * tmp312 + 3 * tmp313 + 1) -
                                          1.0 / 4.0 * tmp311 * tmp314 * (1 - tmp142) * (3 * tmp312 + 9 * tmp313 + 1) - 8191.0 / 672.0 * M_PI) +
                pow(tmp121, 4.0 / 3.0) * ((9271.0 / 504.0) * tmp142 + (65.0 / 18.0) * tmp189 + (287.0 / 48.0) * tmp310 * tmp311 * tmp314 +
                                          tmp312 * (287.0 / 96.0 - 12 * tmp142) + tmp313 * ((1.0 / 24.0) * tmp142 + 287.0 / 96.0) -
                                          89.0 / 48.0 * tmp314 * (tmp310 * tmp311 + tmp318 * tmp322 + tmp321 * tmp323) -
                                          ((7.0 / 24.0) * tmp142 + 89.0 / 96.0) * (tmp313 + ((tmp322) * (tmp322)) + ((tmp323) * (tmp323))) +
                                          (4 * tmp142 - 89.0 / 96.0) * (tmp312 + ((tmp318) * (tmp318)) + ((tmp321) * (tmp321))) - 44711.0 / 9072.0) +
                tmp121 * (-4 * tmp302 - 5.0 / 4.0 * tmp304 + 4 * M_PI) + tmp301 * (-35.0 / 12.0 * tmp142 - 1247.0 / 336.0) + 1) /
               (-S1U2 * tmp135 * (tmp1 * tmp134 + tmp133) +
                S1U2 * (-tmp135 * (Pt * (tmp167 + tmp171) + tmp158 * tmp164) -
                        tmp163 * (-Pt * (6 * m1 + (15.0 / 2.0) * m2) + Pt * (-11.0 / 2.0 * m2 + tmp162))) +
                S1U2 * (-tmp139 * (Pt * ((181.0 / 16.0) * tmp124 + tmp234 + (95.0 / 4.0) * tmp9) +
                                   Pt * ((473.0 / 16.0) * tmp124 + (21.0 / 2.0) * tmp7 + (63.0 / 4.0) * tmp9)) +
                        tmp223 * tmp224 - tmp224 * tmp245 + tmp228 * tmp229 - tmp229 * tmp240) -
                S2U2 * tmp135 * (tmp133 + tmp134 * tmp5) +
                S2U2 * (-tmp135 * (Pt * (tmp167 + tmp178) + tmp158 * tmp164) -
                        tmp163 * (Pt * (-11.0 / 2.0 * m1 + tmp174) - Pt * ((15.0 / 2.0) * m1 + 6 * m2))) +
                S2U2 * (-tmp139 * (Pt * ((181.0 / 16.0) * tmp124 + tmp249 + (95.0 / 4.0) * tmp7) +
                                   Pt * ((473.0 / 16.0) * tmp124 + (63.0 / 4.0) * tmp7 + (21.0 / 2.0) * tmp9)) +
                        tmp224 * tmp247 - tmp224 * tmp256 + tmp228 * tmp229 - tmp229 * tmp251) -
                m1 * tmp139 * ((5.0 / 2.0) * ((S2U2) * (S2U2)) - 7 * tmp146 * tmp5 - 2 * tmp146 + (5.0 / 2.0) * tmp147 + 3 * tmp150 * tmp5) -
                m2 * tmp139 * ((5.0 / 2.0) * ((S1U2) * (S1U2)) - 7 * tmp1 * tmp151 - 2 * tmp151 + (5.0 / 2.0) * tmp152 + tmp155 * tmp79) -
                tmp122 * tmp139 * (6 * tmp128 - 6 * tmp136 + 6 * tmp138) + tmp125 * tmp22 + tmp125 * ((1.0 / 2.0) * tmp127 * tmp145 * tmp22 - tmp40) +
                tmp125 * (-tmp127 * tmp187 * tmp40 + tmp163 * ((3.0 / 4.0) * tmp142 + 1.0 / 4.0) - 1.0 / 8.0 * tmp190 * tmp191 * tmp193 * tmp22) +
                tmp125 * (-1.0 / 16.0 * pow(Pt, 6) * tmp218 * tmp219 * tmp22 - tmp127 * tmp163 * tmp215 - tmp135 * tmp191 * tmp193 * tmp214 -
                          tmp139 * (tmp142 * (109.0 / 12.0 - 21.0 / 32.0 * tmp41) + 1.0 / 8.0)) +
                tmp130 * tmp131 * tmp71 + tmp130 * tmp71 * tmp8 -
                tmp131 * tmp139 *
                    (-Pt * tmp196 - S2U0 * tmp203 - S2U1 * tmp134 * tmp197 + tmp134 * tmp204 + tmp205 * (Pt * tmp30 + Pt) -
                     tmp72 * (-Pt * tmp204 + Pt * tmp206 - S2U0 * tmp207)) -
                tmp139 * tmp8 *
                    (-Pt * tmp209 + S1U0 * tmp203 + S1U1 * tmp134 * tmp197 + tmp134 * tmp210 - tmp211 * (-Pt * tmp69 - Pt) -
                     tmp213 * (-Pt * tmp210 + Pt * tmp212 + S1U0 * tmp207)) -
                3.0 / 2.0 * tmp142 * tmp71 *
                    (2 * ((S1U0 * tmp2 + S2U0 * tmp179) * (S1U0 * tmp2 + S2U0 * tmp179)) -
                     ((S1U1 * tmp2 + S2U1 * tmp179) * (S1U1 * tmp2 + S2U1 * tmp179)) -
                     ((S1U2 * tmp2 + S2U2 * tmp179) * (S1U2 * tmp2 + S2U2 * tmp179))) -
                3.0 / 2.0 * tmp158 * tmp71 * (-15.0 / 2.0 * tmp127 * tmp128 + 3 * tmp127 * tmp136 - 1.0 / 2.0 * tmp127 * tmp138 - tmp127 * tmp157) -
                tmp163 * ((9.0 / 4.0) * tmp127 * tmp131 * tmp146 - 3.0 / 8.0 * tmp127 * tmp146 * tmp177 + (1.0 / 4.0) * tmp127 * tmp147 * tmp177 -
                          tmp131 * tmp150 * tmp166 + tmp146 * tmp184 - tmp150 * tmp167) -
                tmp163 * (-3.0 / 8.0 * tmp127 * tmp151 * tmp170 + (9.0 / 4.0) * tmp127 * tmp151 * tmp8 + (1.0 / 4.0) * tmp127 * tmp152 * tmp170 +
                          tmp151 * tmp184 - tmp155 * tmp166 * tmp8 - tmp155 * tmp167) +
                (-tmp1 * tmp54 * tmp78 * (2 * tmp11 + tmp22 * tmp288 + tmp24 * tmp287 + tmp285 * tmp78 + tmp294 * tmp40 + tmp300 * tmp31 + 1) +
                 tmp1 * tmp13 *
                     (-tmp135 * tmp288 - tmp163 * tmp294 - 2 * tmp22 - 3.0 / 2.0 * tmp24 * tmp285 - 5.0 / 2.0 * tmp287 * tmp31 -
                      7.0 / 2.0 * tmp300 / pow(r, 9.0 / 2.0)) /
                     sqrt(r)) *
                    (Pt * tmp122 * tmp158 + S1U2 * tmp22 * (tmp72 + 2) +
                     S1U2 * (tmp22 * (Pt * (-5.0 / 4.0 * Pt * tmp170 + tmp134 * tmp8 + tmp263) + tmp171 + tmp264) +
                             tmp40 * (-6 * m1 - 13 * m2 + tmp162)) +
                     S1U2 *
                         (tmp22 * tmp240 - tmp223 * tmp40 + tmp245 * tmp40 - tmp257 * (13 * Pt * tmp0 + (53.0 / 4.0) * Pt * tmp4 + tmp133 * tmp222) +
                          tmp257 * (-tmp133 * tmp241 - tmp133 * tmp243 + tmp133 * tmp244) +
                          tmp273 * ((7.0 / 4.0) * m2 * tmp164 * tmp32 - 9.0 / 4.0 * tmp14 * tmp164 + tmp280) + tmp284 +
                          tmp71 * ((327.0 / 8.0) * tmp124 + tmp234 + (21.0 / 2.0) * tmp7 + (79.0 / 2.0) * tmp9)) +
                     S2U2 * tmp22 * (tmp213 + 2) +
                     S2U2 * (tmp22 * (Pt * (-5.0 / 4.0 * Pt * tmp177 + tmp131 * tmp134 + tmp263) + tmp178 + tmp264) +
                             tmp40 * (-13 * m1 - 6 * m2 + tmp174)) +
                     S2U2 *
                         (tmp22 * tmp251 - tmp247 * tmp40 + tmp256 * tmp40 - tmp257 * ((53.0 / 4.0) * Pt * tmp0 + 13 * Pt * tmp4 + tmp133 * tmp246) +
                          tmp257 * (-tmp133 * tmp252 - tmp133 * tmp254 + tmp133 * tmp255) +
                          tmp273 * (-9.0 / 4.0 * tmp164 * tmp192 + (7.0 / 4.0) * tmp164 * tmp250 + tmp280) + tmp284 +
                          tmp71 * ((327.0 / 8.0) * tmp124 + tmp249 + (79.0 / 2.0) * tmp7 + (21.0 / 2.0) * tmp9)) +
                     tmp125 * (-Pt * tmp11 * tmp145 + tmp260 * ((3.0 / 8.0) * tmp142 - 1.0 / 8.0)) +
                     tmp125 * ((1.0 / 2.0) * tmp11 * tmp164 * tmp190 * tmp193 + tmp187 * tmp273 +
                               6 * tmp274 * (-5.0 / 16.0 * tmp142 + (5.0 / 16.0) * tmp189 + 1.0 / 16.0)) +
                     tmp125 *
                         (8 * pow(Pt, 7) * pow(tmp122, 8) * ((35.0 / 128.0) * tmp142 - 35.0 / 64.0 * tmp189 + (35.0 / 128.0) * tmp217 - 5.0 / 128.0) /
                              (pow(m1, 8) * pow(m2, 8)) +
                          (3.0 / 8.0) * tmp11 * tmp218 * tmp274 + tmp214 * tmp22 * tmp260 + tmp215 * tmp229) -
                     tmp131 * tmp258 +
                     tmp131 * tmp71 *
                         (-S2U0 * tmp270 - 3.0 / 2.0 * S2U1 * tmp197 - tmp196 + (3.0 / 2.0) * tmp204 + tmp205 * (tmp30 + 1) -
                          tmp72 * (-S2U0 * tmp272 - tmp204 + tmp206)) +
                     (1.0 / 2.0) * tmp158 * tmp40 * (-15 * Pt * tmp128 + 6 * Pt * tmp136 - Pt * tmp138 - tmp133 * tmp157) - tmp258 * tmp8 +
                     tmp40 * (-3.0 / 4.0 * Pt * tmp146 * tmp177 + (1.0 / 2.0) * Pt * tmp147 * tmp177 - tmp131 * tmp134 * tmp150 +
                              tmp131 * tmp146 * tmp266 + tmp146 * tmp158 * tmp266 - tmp150 * tmp263) +
                     tmp40 * (-3.0 / 4.0 * Pt * tmp151 * tmp170 + (1.0 / 2.0) * Pt * tmp152 * tmp170 - tmp134 * tmp155 * tmp8 +
                              tmp151 * tmp158 * tmp266 + tmp151 * tmp266 * tmp8 - tmp155 * tmp263) +
                     tmp71 * tmp8 *
                         (S1U0 * tmp270 + (3.0 / 2.0) * S1U1 * tmp197 - tmp209 + (3.0 / 2.0) * tmp210 - tmp211 * (-tmp69 - 1) -
                          tmp213 * (S1U0 * tmp272 - tmp210 + tmp212)))) +
           tmp31 * ((1.0 / 4.0) * chi1U0 * chi2U1 * tmp10 * tmp36 * (-6 * tmp1 - 13) +
                    chi1U1 * (-chi1U0 * tmp36 * tmp69 * (tmp1 + 6) - chi2U0 * tmp36 * tmp69 * (13 * tmp1 + 6)) -
                    1.0 / 4.0 * chi2U1 * tmp68 * (6 * tmp1 + 1)) +
           tmp71 * (chi1U1 * chi2U2 * (-tmp74 + tmp77) + chi1U2 * (chi2U1 * tmp74 - chi2U1 * tmp77))) /
          (-1.0 / 2.0 * tmp11 * tmp5 * (-15 * tmp1 - 7 * tmp10 - 7) -
           1.0 / 8.0 * tmp13 * tmp22 * tmp5 * (229 * tmp1 + 363 * tmp10 + 47 * tmp16 + 229 * tmp20 + 47) - ((tmp2) * (tmp2)) * tmp5 -
           tmp24 * (chi1U2 * tmp29 * tmp30 * (11 * tmp1 + 4 * tmp10 + 12) + chi2U2 * tmp26 * (11 * tmp1 + 12 * tmp10 + 4)) -
           tmp31 * (tmp38 * (-357 * tmp1 - 1097 * tmp10 - 842 * tmp16 - 1486 * tmp20 - 144 * tmp34 - 53) +
                    tmp39 * tmp5 * (-842 * tmp1 - 1486 * tmp10 - 357 * tmp16 - 1097 * tmp20 - 53 * tmp34 - 144)) -
           tmp40 * (tmp13 * tmp52 * (tmp49 + tmp50) +
                    (1.0 / 48.0) * tmp36 * tmp5 *
                        (-2608 * tmp1 - 7324 * tmp10 - 7324 * tmp16 - 10161 * tmp20 - 2608 * tmp34 - 363 * tmp61 * tmp62 - 363) -
                    1.0 / 16.0 * tmp41 - tmp42 * tmp43 - tmp43 * tmp44 - tmp45 * tmp47 - tmp47 * tmp48 + tmp5 * tmp54 * tmp58 * (tmp10 + tmp59 + 9) +
                    tmp54 * tmp60 * (tmp1 + 9 * tmp10 + 9 * tmp20) + tmp55 * tmp56 + tmp55 * tmp57));
    }

    commondata->initial_p_t = Pt;
    commondata->initial_p_r = Pr;
    printf("p_t, p_r = %.15e %.15e\n", Pt, Pr);
  }
}
