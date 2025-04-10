#include "../BHaH_defines.h"
/**
 * Compute AUXEVOL grid functions at a single point.
 */
void auxevol_gfs_single_point__rfm__HoleySinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                       const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict psi_background,
                                                       REAL *restrict ADD_times_AUU) {
  // Load necessary parameters from params_struct
  const REAL AMPL = params->AMPL;
  const REAL SINHW = params->SINHW;

  // Load necessary parameters from commondata_struct
  const REAL P0_x = commondata->P0_x;
  const REAL P0_y = commondata->P0_y;
  const REAL P0_z = commondata->P0_z;
  const REAL P1_x = commondata->P1_x;
  const REAL P1_y = commondata->P1_y;
  const REAL P1_z = commondata->P1_z;
  const REAL S0_x = commondata->S0_x;
  const REAL S0_y = commondata->S0_y;
  const REAL S0_z = commondata->S0_z;
  const REAL S1_x = commondata->S1_x;
  const REAL S1_y = commondata->S1_y;
  const REAL S1_z = commondata->S1_z;
  const REAL bare_mass_0 = commondata->bare_mass_0;
  const REAL bare_mass_1 = commondata->bare_mass_1;
  const REAL zPunc = commondata->zPunc;

  const REAL tmp0 = (1.0 / (SINHW));
  const REAL tmp8 = cos(xx2);
  const REAL tmp9 = sin(xx1);
  const REAL tmp12 = sin(xx2);
  const REAL tmp1 = exp(tmp0) - exp(-tmp0);
  const REAL tmp2 = (1.0 / (tmp1));
  const REAL tmp4 = exp(tmp0 * xx0) - exp(-tmp0 * xx0);
  const REAL tmp5 = AMPL * tmp2 * tmp4;
  const REAL tmp10 = ((AMPL) * (AMPL)) * ((tmp4) * (tmp4)) * ((tmp9) * (tmp9)) / ((tmp1) * (tmp1));
  const REAL tmp6 = tmp5 * cos(xx1);
  const REAL tmp11 = tmp10 * ((tmp8) * (tmp8));
  const REAL tmp13 = tmp10 * ((tmp12) * (tmp12));
  const REAL tmp78 = 6 * tmp10 * tmp12 * tmp8;
  const REAL tmp7 = tmp6 + zPunc;
  const REAL tmp16 = tmp6 - zPunc;
  const REAL tmp21 = tmp12 * tmp5 * tmp9;
  const REAL tmp24 = tmp5 * tmp8 * tmp9;
  const REAL tmp15 = tmp11 + tmp13 + ((tmp7) * (tmp7));
  const REAL tmp17 = tmp11 + tmp13 + ((tmp16) * (tmp16));
  const REAL tmp25 = S1_x * tmp7 - S1_z * tmp24;
  const REAL tmp28 = S0_x * tmp16 - S0_z * tmp24;
  const REAL tmp60 = AMPL * S1_y * tmp2 * tmp4 * tmp8 * tmp9 - S1_x * tmp21;
  const REAL tmp65 = AMPL * S0_y * tmp2 * tmp4 * tmp8 * tmp9 - S0_x * tmp21;
  const REAL tmp70 = 3 * tmp24;
  const REAL tmp71 = 3 * tmp21;
  const REAL tmp73 = (3.0 / 2.0) * tmp21;
  const REAL tmp75 = (3.0 / 2.0) * tmp24;
  const REAL tmp18 = pow(tmp17, -3.0 / 2.0);
  const REAL tmp19 = pow(tmp15, -3.0 / 2.0);
  const REAL tmp23 = pow(tmp15, -5.0 / 2.0);
  const REAL tmp27 = pow(tmp17, -5.0 / 2.0);
  const REAL tmp34 = -1.0 / 4.0 * P1_x * tmp24 - 1.0 / 4.0 * P1_y * tmp21 - 1.0 / 4.0 * P1_z * tmp7;
  const REAL tmp39 = -1.0 / 4.0 * P0_x * tmp24 - 1.0 / 4.0 * P0_y * tmp21 - 1.0 / 4.0 * P0_z * tmp16;
  const REAL tmp48 = -S1_y * tmp7 + S1_z * tmp21;
  const REAL tmp51 = -S0_y * tmp16 + S0_z * tmp21;
  const REAL tmp26 = tmp23 * tmp25;
  const REAL tmp29 = tmp27 * tmp28;
  const REAL tmp35 = tmp23 * tmp34;
  const REAL tmp40 = tmp27 * tmp39;
  const REAL tmp41 = 2 * tmp34;
  const REAL tmp43 = 2 * tmp39;
  const REAL tmp49 = tmp23 * tmp48;
  const REAL tmp52 = tmp27 * tmp51;
  const REAL tmp54 = P1_z * tmp19 * tmp7;
  const REAL tmp55 = -tmp16 * tmp18;
  const REAL tmp57 = -tmp19 * tmp7;
  const REAL tmp59 = P0_z * tmp16 * tmp18;
  const REAL tmp63 = tmp23 * (-3 * tmp6 - 3 * zPunc);
  const REAL tmp66 = tmp27 * (-3 * tmp6 + 3 * zPunc);
  const REAL tmp45 = -tmp18 * tmp43 - tmp19 * tmp41;
  const REAL tmp83 = -P0_x * tmp18 * tmp24 - P0_y * tmp18 * tmp21 + (7.0 / 6.0) * P0_z * tmp55 - P1_x * tmp19 * tmp24 - P1_y * tmp19 * tmp21 +
                     (7.0 / 6.0) * P1_z * tmp57 + tmp11 * tmp23 * tmp41 + tmp11 * tmp27 * tmp43 + tmp13 * tmp23 * tmp41 + tmp13 * tmp27 * tmp43 -
                     2.0 / 3.0 * tmp16 * tmp39 * tmp66 + 2 * tmp21 * tmp26 + 2 * tmp21 * tmp29 + 2 * tmp24 * tmp49 + 2 * tmp24 * tmp52 -
                     2.0 / 3.0 * tmp34 * tmp63 * tmp7 + tmp45 + (1.0 / 6.0) * tmp54 + (1.0 / 6.0) * tmp59 - 2.0 / 3.0 * tmp60 * tmp63 -
                     2.0 / 3.0 * tmp65 * tmp66;
  *psi_background = (1.0 / 2.0) * bare_mass_0 / sqrt(tmp17) + (1.0 / 2.0) * bare_mass_1 / sqrt(tmp15) + 1;
  *ADD_times_AUU = ((3 * AMPL * P0_x * tmp18 * tmp2 * tmp4 * tmp8 * tmp9 + 3 * AMPL * P1_x * tmp19 * tmp2 * tmp4 * tmp8 * tmp9 - 6 * tmp11 * tmp35 -
                     6 * tmp11 * tmp40 - 6 * tmp24 * tmp49 - 6 * tmp24 * tmp52 - tmp45) *
                    (3 * AMPL * P0_x * tmp18 * tmp2 * tmp4 * tmp8 * tmp9 + 3 * AMPL * P1_x * tmp19 * tmp2 * tmp4 * tmp8 * tmp9 - 6 * tmp11 * tmp35 -
                     6 * tmp11 * tmp40 - 6 * tmp24 * tmp49 - 6 * tmp24 * tmp52 - tmp45)) +
                   ((3 * AMPL * P0_y * tmp12 * tmp18 * tmp2 * tmp4 * tmp9 + 3 * AMPL * P1_y * tmp12 * tmp19 * tmp2 * tmp4 * tmp9 - 6 * tmp13 * tmp35 -
                     6 * tmp13 * tmp40 - 6 * tmp21 * tmp26 - 6 * tmp21 * tmp29 - tmp45) *
                    (3 * AMPL * P0_y * tmp12 * tmp18 * tmp2 * tmp4 * tmp9 + 3 * AMPL * P1_y * tmp12 * tmp19 * tmp2 * tmp4 * tmp9 - 6 * tmp13 * tmp35 -
                     6 * tmp13 * tmp40 - 6 * tmp21 * tmp26 - 6 * tmp21 * tmp29 - tmp45)) +
                   ((-7.0 / 2.0 * P0_z * tmp55 - 7.0 / 2.0 * P1_z * tmp57 + tmp16 * tmp43 * tmp66 + tmp18 * tmp43 + tmp19 * tmp41 +
                     tmp41 * tmp63 * tmp7 - 1.0 / 2.0 * tmp54 - 1.0 / 2.0 * tmp59 + 2 * tmp60 * tmp63 + 2 * tmp65 * tmp66) *
                    (-7.0 / 2.0 * P0_z * tmp55 - 7.0 / 2.0 * P1_z * tmp57 + tmp16 * tmp43 * tmp66 + tmp18 * tmp43 + tmp19 * tmp41 +
                     tmp41 * tmp63 * tmp7 - 1.0 / 2.0 * tmp54 - 1.0 / 2.0 * tmp59 + 2 * tmp60 * tmp63 + 2 * tmp65 * tmp66)) +
                   2 * ((P0_x * tmp18 * tmp73 + P0_y * tmp18 * tmp75 + P1_x * tmp19 * tmp73 + P1_y * tmp19 * tmp75 - tmp26 * tmp70 - tmp29 * tmp70 -
                         tmp35 * tmp78 - tmp40 * tmp78 - tmp49 * tmp71 - tmp52 * tmp71 + tmp83) *
                        (P0_x * tmp18 * tmp73 + P0_y * tmp18 * tmp75 + P1_x * tmp19 * tmp73 + P1_y * tmp19 * tmp75 - tmp26 * tmp70 - tmp29 * tmp70 -
                         tmp35 * tmp78 - tmp40 * tmp78 - tmp49 * tmp71 - tmp52 * tmp71 + tmp83)) +
                   2 * ((-1.0 / 4.0 * P0_x * tmp16 * tmp18 - 7.0 / 4.0 * P0_x * tmp55 + P0_z * tmp18 * tmp75 - 1.0 / 4.0 * P1_x * tmp19 * tmp7 -
                         7.0 / 4.0 * P1_x * tmp57 + P1_z * tmp19 * tmp75 - tmp16 * tmp40 * tmp70 - tmp23 * tmp60 * tmp70 + tmp24 * tmp34 * tmp63 +
                         tmp24 * tmp39 * tmp66 - tmp27 * tmp65 * tmp70 - tmp35 * tmp7 * tmp70 + tmp48 * tmp63 + tmp51 * tmp66 + tmp83) *
                        (-1.0 / 4.0 * P0_x * tmp16 * tmp18 - 7.0 / 4.0 * P0_x * tmp55 + P0_z * tmp18 * tmp75 - 1.0 / 4.0 * P1_x * tmp19 * tmp7 -
                         7.0 / 4.0 * P1_x * tmp57 + P1_z * tmp19 * tmp75 - tmp16 * tmp40 * tmp70 - tmp23 * tmp60 * tmp70 + tmp24 * tmp34 * tmp63 +
                         tmp24 * tmp39 * tmp66 - tmp27 * tmp65 * tmp70 - tmp35 * tmp7 * tmp70 + tmp48 * tmp63 + tmp51 * tmp66 + tmp83)) +
                   2 * ((-1.0 / 4.0 * P0_y * tmp16 * tmp18 - 7.0 / 4.0 * P0_y * tmp55 + P0_z * tmp18 * tmp73 - 1.0 / 4.0 * P1_y * tmp19 * tmp7 -
                         7.0 / 4.0 * P1_y * tmp57 + P1_z * tmp19 * tmp73 - tmp16 * tmp40 * tmp71 + tmp21 * tmp34 * tmp63 + tmp21 * tmp39 * tmp66 -
                         tmp23 * tmp60 * tmp71 + tmp25 * tmp63 - tmp27 * tmp65 * tmp71 + tmp28 * tmp66 - tmp35 * tmp7 * tmp71 + tmp83) *
                        (-1.0 / 4.0 * P0_y * tmp16 * tmp18 - 7.0 / 4.0 * P0_y * tmp55 + P0_z * tmp18 * tmp73 - 1.0 / 4.0 * P1_y * tmp19 * tmp7 -
                         7.0 / 4.0 * P1_y * tmp57 + P1_z * tmp19 * tmp73 - tmp16 * tmp40 * tmp71 + tmp21 * tmp34 * tmp63 + tmp21 * tmp39 * tmp66 -
                         tmp23 * tmp60 * tmp71 + tmp25 * tmp63 - tmp27 * tmp65 * tmp71 + tmp28 * tmp66 - tmp35 * tmp7 * tmp71 + tmp83));
} // END FUNCTION auxevol_gfs_single_point__rfm__HoleySinhSpherical
