#include "../BHaH_defines.h"
/**
 * Compute AUXEVOL grid functions at a single point.
 */
__device__ void auxevol_gfs_single_point__rfm__SinhSymTP(const size_t streamid, const REAL xx0, const REAL xx1, const REAL xx2,
                                                         REAL *restrict psi_background, REAL *restrict ADD_times_AUU) {
  // Load necessary parameters from params_struct
  const REAL AMAX = d_params[streamid].AMAX;
  const REAL SINHWAA = d_params[streamid].SINHWAA;
  const REAL bScale = d_params[streamid].bScale;

  // Load necessary parameters from commondata_struct
  const REAL P0_x = d_commondata.P0_x;
  const REAL P0_y = d_commondata.P0_y;
  const REAL P0_z = d_commondata.P0_z;
  const REAL P1_x = d_commondata.P1_x;
  const REAL P1_y = d_commondata.P1_y;
  const REAL P1_z = d_commondata.P1_z;
  const REAL S0_x = d_commondata.S0_x;
  const REAL S0_y = d_commondata.S0_y;
  const REAL S0_z = d_commondata.S0_z;
  const REAL S1_x = d_commondata.S1_x;
  const REAL S1_y = d_commondata.S1_y;
  const REAL S1_z = d_commondata.S1_z;
  const REAL bare_mass_0 = d_commondata.bare_mass_0;
  const REAL bare_mass_1 = d_commondata.bare_mass_1;
  const REAL zPunc = d_commondata.zPunc;

  const REAL tmp0 = (1.0 / (SINHWAA));
  const REAL tmp7 = cos(xx2);
  const REAL tmp8 = sin(xx1);
  const REAL tmp11 = sin(xx2);
  const REAL tmp1 = exp(tmp0) - exp(-tmp0);
  const REAL tmp3 = exp(tmp0 * xx0) - exp(-tmp0 * xx0);
  const REAL tmp18 = (1.0 / (tmp1));
  const REAL tmp4 = ((AMAX) * (AMAX)) * ((tmp3) * (tmp3)) / ((tmp1) * (tmp1));
  const REAL tmp21 = AMAX * tmp18 * tmp3 * tmp8;
  const REAL tmp5 = sqrt(((bScale) * (bScale)) + tmp4) * cos(xx1);
  const REAL tmp9 = tmp4 * ((tmp8) * (tmp8));
  const REAL tmp22 = tmp11 * tmp21;
  const REAL tmp25 = tmp21 * tmp7;
  const REAL tmp6 = tmp5 + zPunc;
  const REAL tmp10 = ((tmp7) * (tmp7)) * tmp9;
  const REAL tmp12 = ((tmp11) * (tmp11)) * tmp9;
  const REAL tmp15 = -tmp5 + zPunc;
  const REAL tmp60 = AMAX * S1_y * tmp18 * tmp3 * tmp7 * tmp8 - S1_x * tmp22;
  const REAL tmp65 = AMAX * S0_y * tmp18 * tmp3 * tmp7 * tmp8 - S0_x * tmp22;
  const REAL tmp70 = 3 * tmp25;
  const REAL tmp71 = 3 * tmp22;
  const REAL tmp73 = (3.0 / 2.0) * tmp22;
  const REAL tmp75 = (3.0 / 2.0) * tmp25;
  const REAL tmp78 = 6 * tmp11 * tmp7 * tmp9;
  const REAL tmp26 = S1_x * tmp6 - S1_z * tmp25;
  const REAL tmp51 = S0_y * tmp15 + S0_z * tmp22;
  const REAL tmp14 = tmp10 + tmp12 + ((tmp6) * (tmp6));
  const REAL tmp17 = tmp10 + tmp12 + ((tmp15) * (tmp15));
  const REAL tmp29 = -S0_x * tmp15 - S0_z * tmp25;
  const REAL tmp35 = -1.0 / 4.0 * P1_x * tmp25 - 1.0 / 4.0 * P1_y * tmp22 - 1.0 / 4.0 * P1_z * tmp6;
  const REAL tmp39 = -P0_z * tmp15;
  const REAL tmp49 = -S1_y * tmp6 + S1_z * tmp22;
  const REAL tmp19 = pow(tmp17, -3.0 / 2.0);
  const REAL tmp20 = pow(tmp14, -3.0 / 2.0);
  const REAL tmp24 = pow(tmp14, -5.0 / 2.0);
  const REAL tmp28 = pow(tmp17, -5.0 / 2.0);
  const REAL tmp40 = -1.0 / 4.0 * P0_x * tmp25 - 1.0 / 4.0 * P0_y * tmp22 - 1.0 / 4.0 * tmp39;
  const REAL tmp42 = 2 * tmp35;
  const REAL tmp27 = tmp24 * tmp26;
  const REAL tmp30 = tmp28 * tmp29;
  const REAL tmp36 = tmp24 * tmp35;
  const REAL tmp41 = tmp28 * tmp40;
  const REAL tmp44 = 2 * tmp40;
  const REAL tmp50 = tmp24 * tmp49;
  const REAL tmp52 = tmp28 * tmp51;
  const REAL tmp54 = P1_z * tmp20 * tmp6;
  const REAL tmp57 = -tmp20 * tmp6;
  const REAL tmp63 = tmp24 * (-3 * tmp5 - 3 * zPunc);
  const REAL tmp66 = tmp28 * (-3 * tmp5 + 3 * zPunc);
  const REAL tmp86 = -tmp15 * tmp19;
  const REAL tmp56 = P0_z * tmp15 * tmp19;
  const REAL tmp69 = -tmp15 * tmp66;
  const REAL tmp84 = (7.0 / 4.0) * tmp15 * tmp19;
  const REAL tmp95 = -tmp15 * tmp41;
  const REAL tmp46 = -tmp19 * tmp44 - tmp20 * tmp42;
  const REAL tmp83 = -P0_x * tmp19 * tmp25 - P0_y * tmp19 * tmp22 - P1_x * tmp20 * tmp25 - P1_y * tmp20 * tmp22 + (7.0 / 6.0) * P1_z * tmp57 +
                     tmp10 * tmp24 * tmp42 + tmp10 * tmp28 * tmp44 + tmp12 * tmp24 * tmp42 + tmp12 * tmp28 * tmp44 + (1.0 / 6.0) * tmp19 * tmp39 +
                     2 * tmp22 * tmp27 + 2 * tmp22 * tmp30 + 2 * tmp25 * tmp50 + 2 * tmp25 * tmp52 - 2.0 / 3.0 * tmp35 * tmp6 * tmp63 -
                     2.0 / 3.0 * tmp40 * tmp69 + tmp46 + (1.0 / 6.0) * tmp54 + (7.0 / 6.0) * tmp56 - 2.0 / 3.0 * tmp60 * tmp63 -
                     2.0 / 3.0 * tmp65 * tmp66;
  *psi_background = (1.0 / 2.0) * bare_mass_0 / sqrt(tmp17) + (1.0 / 2.0) * bare_mass_1 / sqrt(tmp14) + 1;
  *ADD_times_AUU = ((3 * AMAX * P0_x * tmp18 * tmp19 * tmp3 * tmp7 * tmp8 + 3 * AMAX * P1_x * tmp18 * tmp20 * tmp3 * tmp7 * tmp8 - 6 * tmp10 * tmp36 -
                     6 * tmp10 * tmp41 - 6 * tmp25 * tmp50 - 6 * tmp25 * tmp52 - tmp46) *
                    (3 * AMAX * P0_x * tmp18 * tmp19 * tmp3 * tmp7 * tmp8 + 3 * AMAX * P1_x * tmp18 * tmp20 * tmp3 * tmp7 * tmp8 - 6 * tmp10 * tmp36 -
                     6 * tmp10 * tmp41 - 6 * tmp25 * tmp50 - 6 * tmp25 * tmp52 - tmp46)) +
                   ((3 * AMAX * P0_y * tmp11 * tmp18 * tmp19 * tmp3 * tmp8 + 3 * AMAX * P1_y * tmp11 * tmp18 * tmp20 * tmp3 * tmp8 -
                     6 * tmp12 * tmp36 - 6 * tmp12 * tmp41 - 6 * tmp22 * tmp27 - 6 * tmp22 * tmp30 - tmp46) *
                    (3 * AMAX * P0_y * tmp11 * tmp18 * tmp19 * tmp3 * tmp8 + 3 * AMAX * P1_y * tmp11 * tmp18 * tmp20 * tmp3 * tmp8 -
                     6 * tmp12 * tmp36 - 6 * tmp12 * tmp41 - 6 * tmp22 * tmp27 - 6 * tmp22 * tmp30 - tmp46)) +
                   ((-7.0 / 2.0 * P1_z * tmp57 - 1.0 / 2.0 * tmp19 * tmp39 + tmp19 * tmp44 + tmp20 * tmp42 + tmp42 * tmp6 * tmp63 + tmp44 * tmp69 -
                     1.0 / 2.0 * tmp54 - 7.0 / 2.0 * tmp56 + 2 * tmp60 * tmp63 + 2 * tmp65 * tmp66) *
                    (-7.0 / 2.0 * P1_z * tmp57 - 1.0 / 2.0 * tmp19 * tmp39 + tmp19 * tmp44 + tmp20 * tmp42 + tmp42 * tmp6 * tmp63 + tmp44 * tmp69 -
                     1.0 / 2.0 * tmp54 - 7.0 / 2.0 * tmp56 + 2 * tmp60 * tmp63 + 2 * tmp65 * tmp66)) +
                   2 * ((P0_x * tmp19 * tmp73 + P0_y * tmp19 * tmp75 + P1_x * tmp20 * tmp73 + P1_y * tmp20 * tmp75 - tmp27 * tmp70 - tmp30 * tmp70 -
                         tmp36 * tmp78 - tmp41 * tmp78 - tmp50 * tmp71 - tmp52 * tmp71 + tmp83) *
                        (P0_x * tmp19 * tmp73 + P0_y * tmp19 * tmp75 + P1_x * tmp20 * tmp73 + P1_y * tmp20 * tmp75 - tmp27 * tmp70 - tmp30 * tmp70 -
                         tmp36 * tmp78 - tmp41 * tmp78 - tmp50 * tmp71 - tmp52 * tmp71 + tmp83)) +
                   2 * ((-P0_x * tmp84 - 1.0 / 4.0 * P0_x * tmp86 + P0_z * tmp19 * tmp75 - 1.0 / 4.0 * P1_x * tmp20 * tmp6 -
                         7.0 / 4.0 * P1_x * tmp57 + P1_z * tmp20 * tmp75 - tmp24 * tmp60 * tmp70 + tmp25 * tmp35 * tmp63 + tmp25 * tmp40 * tmp66 -
                         tmp28 * tmp65 * tmp70 - tmp36 * tmp6 * tmp70 + tmp49 * tmp63 + tmp51 * tmp66 - tmp70 * tmp95 + tmp83) *
                        (-P0_x * tmp84 - 1.0 / 4.0 * P0_x * tmp86 + P0_z * tmp19 * tmp75 - 1.0 / 4.0 * P1_x * tmp20 * tmp6 -
                         7.0 / 4.0 * P1_x * tmp57 + P1_z * tmp20 * tmp75 - tmp24 * tmp60 * tmp70 + tmp25 * tmp35 * tmp63 + tmp25 * tmp40 * tmp66 -
                         tmp28 * tmp65 * tmp70 - tmp36 * tmp6 * tmp70 + tmp49 * tmp63 + tmp51 * tmp66 - tmp70 * tmp95 + tmp83)) +
                   2 * ((-P0_y * tmp84 - 1.0 / 4.0 * P0_y * tmp86 + P0_z * tmp19 * tmp73 - 1.0 / 4.0 * P1_y * tmp20 * tmp6 -
                         7.0 / 4.0 * P1_y * tmp57 + P1_z * tmp20 * tmp73 + tmp22 * tmp35 * tmp63 + tmp22 * tmp40 * tmp66 - tmp24 * tmp60 * tmp71 +
                         tmp26 * tmp63 - tmp28 * tmp65 * tmp71 + tmp29 * tmp66 - tmp36 * tmp6 * tmp71 - tmp71 * tmp95 + tmp83) *
                        (-P0_y * tmp84 - 1.0 / 4.0 * P0_y * tmp86 + P0_z * tmp19 * tmp73 - 1.0 / 4.0 * P1_y * tmp20 * tmp6 -
                         7.0 / 4.0 * P1_y * tmp57 + P1_z * tmp20 * tmp73 + tmp22 * tmp35 * tmp63 + tmp22 * tmp40 * tmp66 - tmp24 * tmp60 * tmp71 +
                         tmp26 * tmp63 - tmp28 * tmp65 * tmp71 + tmp29 * tmp66 - tmp36 * tmp6 * tmp71 - tmp71 * tmp95 + tmp83));
} // END FUNCTION auxevol_gfs_single_point__rfm__SinhSymTP
