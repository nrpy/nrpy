#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include <math.h>
#include <stdint.h>

#include "BHaH_defines.h"

/**
 *  Computes the 10 unique components of the KerrSchild_Cartesian metric $g_{mu nu}$ for a photon particle.
 *     @param commondata Struct containing global spacetime parameters.
 *     @param f_local Thread-local array containing the 1D flattened state vector.
 *     @param metric_local Thread-local array where the symmetric metric components are stored.
 */
BHAH_HD_INLINE void g4DD_metric_KerrSchild_Cartesian(const commondata_struct *restrict commondata, const double *restrict f_local,
                                                     double *restrict metric_local) {
#include "set_CodeParameters.h"
  //==========================================
  // METRIC EVALUATION & THREAD-LOCAL UNPACKING
  //==========================================
  // Extract spatial coordinates $x^i$ and compute $g_{mu nu}$.
  // Unpack position coordinates $x^i$ from the thread-local state vector.
  // Evaluated at compile time for state vector size: 9
  const double x = f_local[1];
  const double y = f_local[2];
  const double z = f_local[3];

  const REAL tmp0 = ((a_spin) * (a_spin));
  const REAL tmp4 = ((z) * (z));
  const REAL tmp6 = (1.0 / 2.0) * tmp4 + (1.0 / 2.0) * ((x) * (x)) + (1.0 / 2.0) * ((y) * (y)) +
                    (1.0 / 2.0) * sqrt(4 * tmp0 * tmp4 + ((-tmp0 + tmp4 + ((x) * (x)) + ((y) * (y))) * (-tmp0 + tmp4 + ((x) * (x)) + ((y) * (y)))));
  const REAL tmp7 = -1.0 / 2.0 * tmp0 + tmp6;
  const REAL tmp12 = (1.0 / 2.0) * tmp0 + tmp6;
  const REAL tmp8 = 2 * M_scale / (tmp0 * tmp4 + ((tmp7) * (tmp7)));
  const REAL tmp10 = sqrt(tmp7);
  const REAL tmp13 = (1.0 / (tmp12));
  const REAL tmp9 = pow(tmp7, 3.0 / 2.0) * tmp8;
  const REAL tmp11 = a_spin * y + tmp10 * x;
  const REAL tmp15 = -a_spin * x + tmp10 * y;
  const REAL tmp16 = tmp7 * tmp8 * z;
  const REAL tmp17 = tmp9 / ((tmp12) * (tmp12));
  metric_local[0] = tmp9 - 1;
  metric_local[1] = tmp11 * tmp13 * tmp9;
  metric_local[2] = tmp13 * tmp15 * tmp9;
  metric_local[3] = tmp16;
  metric_local[4] = ((tmp11) * (tmp11)) * tmp17 + 1;
  metric_local[5] = tmp11 * tmp15 * tmp17;
  metric_local[6] = tmp11 * tmp13 * tmp16;
  metric_local[7] = ((tmp15) * (tmp15)) * tmp17 + 1;
  metric_local[8] = tmp13 * tmp15 * tmp16;
  metric_local[9] = tmp10 * tmp4 * tmp8 + 1;

} // END FUNCTION: g4DD_metric_KerrSchild_Cartesian

#include "BHaH_defines.h"

/**
 *  Computes the 40 unique Christoffel symbols $Gamma^alpha_{mu nu}$ for the KerrSchild_Cartesian metric.
 *     @param commondata Struct containing global spacetime parameters.
 *     @param f_local Thread-local array containing the 1D flattened state vector.
 *     @param Gamma_local Thread-local array where the connection components are stored.
 */
BHAH_HD_INLINE void connections_KerrSchild_Cartesian(const commondata_struct *restrict commondata, const double *restrict f_local,
                                                     double *restrict Gamma_local) {
#include "set_CodeParameters.h"
  //==========================================
  // CONNECTION EVALUATION & THREAD-LOCAL UNPACKING
  //==========================================
  // Extract spatial coordinates $x^i$ and compute $Gamma^alpha_{mu nu}$.
  // Unpack position coordinates $x^i$ from the thread-local state vector.
  // Evaluated at compile time for state vector size: 9
  const double x = f_local[1];
  const double y = f_local[2];
  const double z = f_local[3];

  const REAL tmp1 = ((a_spin) * (a_spin));
  const REAL tmp2 = ((z) * (z));
  const REAL tmp22 = 2 * x;
  const REAL tmp37 = ((M_scale) * (M_scale) * (M_scale));
  const REAL tmp42 = ((M_scale) * (M_scale));
  const REAL tmp74 = 4 * z;
  const REAL tmp86 = 2 * y;
  const REAL tmp93 = 2 * z;
  const REAL tmp6 = -tmp1 + tmp2 + ((x) * (x)) + ((y) * (y));
  const REAL tmp7 = sqrt(4 * tmp1 * tmp2 + ((tmp6) * (tmp6)));
  const REAL tmp8 = (1.0 / (tmp7));
  const REAL tmp12 = (1.0 / 2.0) * tmp2 + (1.0 / 2.0) * tmp7 + (1.0 / 2.0) * ((x) * (x)) + (1.0 / 2.0) * ((y) * (y));
  const REAL tmp9 = tmp6 * tmp8;
  const REAL tmp13 = -1.0 / 2.0 * tmp1 + tmp12;
  const REAL tmp32 = (1.0 / 2.0) * tmp1 + tmp12;
  const REAL tmp94 = tmp8 * (tmp1 * tmp74 + tmp6 * tmp93);
  const REAL tmp10 = (3.0 / 2.0) * tmp9 * x + (3.0 / 2.0) * x;
  const REAL tmp14 = sqrt(tmp13);
  const REAL tmp16 = ((tmp13) * (tmp13));
  const REAL tmp23 = tmp22 * tmp9 + tmp22;
  const REAL tmp24 = pow(tmp13, 5.0 / 2.0);
  const REAL tmp33 = (1.0 / ((tmp32) * (tmp32) * (tmp32)));
  const REAL tmp36 = pow(tmp13, 7.0 / 2.0);
  const REAL tmp43 = ((tmp13) * (tmp13) * (tmp13));
  const REAL tmp45 = (1.0 / (tmp32));
  const REAL tmp46 = (1.0 / ((tmp32) * (tmp32)));
  const REAL tmp47 = pow(tmp13, 3.0 / 2.0);
  const REAL tmp54 = (1.0 / ((tmp32) * (tmp32) * (tmp32) * (tmp32)));
  const REAL tmp73 = ((tmp13) * (tmp13) * (tmp13) * (tmp13)) * z;
  const REAL tmp84 = (3.0 / 2.0) * tmp9 * y + (3.0 / 2.0) * y;
  const REAL tmp87 = tmp86 * tmp9 + tmp86;
  const REAL tmp108 = tmp9 * x + x;
  const REAL tmp112 = (1.0 / 4.0) * tmp94 + (1.0 / 2.0) * z;
  const REAL tmp116 = (1.0 / 2.0) * tmp94 + z;
  const REAL tmp123 = (1.0 / 2.0) * tmp9 * y + (1.0 / 2.0) * y;
  const REAL tmp129 = tmp9 * y + y;
  const REAL tmp137 = (1.0 / 2.0) * tmp9 * x + (1.0 / 2.0) * x;
  const REAL tmp15 = 2 * tmp14;
  const REAL tmp17 = tmp1 * tmp2 + tmp16;
  const REAL tmp28 = -a_spin * x + tmp14 * y;
  const REAL tmp48 = 2 * tmp47;
  const REAL tmp61 = tmp16 * tmp2;
  const REAL tmp124 = (1.0 / (tmp14));
  const REAL tmp170 = -tmp93 - tmp94;
  const REAL tmp18 = (1.0 / (tmp17));
  const REAL tmp25 = (1.0 / ((tmp17) * (tmp17)));
  const REAL tmp29 = ((tmp28) * (tmp28));
  const REAL tmp31 = a_spin * y + tmp14 * x;
  const REAL tmp38 = (1.0 / ((tmp17) * (tmp17) * (tmp17)));
  const REAL tmp138 = -a_spin + tmp124 * tmp137 * y;
  const REAL tmp140 = tmp28 * tmp45;
  const REAL tmp162 = 4 * tmp24 * tmp45;
  const REAL tmp181 = tmp28 * tmp46;
  const REAL tmp225 = tmp123 * tmp124 * y + tmp14;
  const REAL tmp19 = M_scale * tmp18;
  const REAL tmp55 = ((tmp31) * (tmp31));
  const REAL tmp59 = tmp25 * tmp42;
  const REAL tmp65 = tmp29 * tmp46;
  const REAL tmp78 = tmp29 * tmp43 * tmp46;
  const REAL tmp97 = M_scale * tmp25;
  const REAL tmp115 = tmp31 * tmp45;
  const REAL tmp126 = a_spin + tmp123 * tmp124 * x;
  const REAL tmp155 = tmp31 * tmp46;
  const REAL tmp164 = tmp124 * tmp137 * x + tmp14;
  const REAL tmp20 = tmp15 * tmp19;
  const REAL tmp35 = tmp29 * tmp31 * tmp33;
  const REAL tmp40 = 8 * tmp37 * tmp38;
  const REAL tmp49 = tmp19 * tmp48;
  const REAL tmp60 = 4 * tmp59;
  const REAL tmp69 = tmp46 * tmp55;
  const REAL tmp71 = tmp29 * tmp54 * tmp55;
  const REAL tmp79 = tmp43 * tmp46 * tmp55;
  const REAL tmp98 = tmp97 * (-tmp1 * tmp93 - tmp13 * (tmp93 + tmp94));
  const REAL tmp103 = tmp13 * tmp19;
  const REAL tmp109 = tmp23 * tmp97;
  const REAL tmp133 = tmp87 * tmp97;
  const REAL tmp159 = 4 * tmp19 * tmp45;
  const REAL tmp215 = 2 * tmp124 * tmp19 * tmp2;
  const REAL tmp235 = tmp19 * tmp74;
  const REAL tmp21 = tmp10 * tmp20;
  const REAL tmp41 = tmp2 * tmp36 * tmp40;
  const REAL tmp44 = tmp2 * tmp20 + 1;
  const REAL tmp50 = tmp46 * tmp49;
  const REAL tmp52 = tmp45 * tmp49;
  const REAL tmp56 = tmp49 - 1;
  const REAL tmp72 = 16 * tmp37 * tmp38 * tmp71;
  const REAL tmp85 = tmp20 * tmp84;
  const REAL tmp95 = tmp20 * ((3.0 / 4.0) * tmp94 + (3.0 / 2.0) * z);
  const REAL tmp99 = tmp48 * tmp98;
  const REAL tmp102 = tmp40 * tmp73;
  const REAL tmp104 = tmp103 * tmp93;
  const REAL tmp111 = tmp109 * tmp16 * tmp93;
  const REAL tmp113 = tmp103 * tmp112;
  const REAL tmp134 = 2 * tmp133 * tmp24;
  const REAL tmp139 = 2 * tmp109 * tmp24;
  const REAL tmp145 = tmp133 * tmp16 * tmp93;
  const REAL tmp151 = tmp33 * tmp49;
  const REAL tmp157 = 4 * tmp19 * tmp47;
  const REAL tmp168 = 4 * tmp103;
  const REAL tmp182 = 4 * tmp14 * tmp19;
  const REAL tmp197 = tmp129 * tmp19 * tmp93;
  const REAL tmp208 = tmp108 * tmp19 * tmp93;
  const REAL tmp218 = tmp109 * tmp2 * tmp48;
  const REAL tmp231 = tmp133 * tmp2 * tmp48;
  const REAL tmp236 = tmp112 * tmp215 + tmp14 * tmp235 + tmp15 * tmp2 * tmp98;
  const REAL tmp237 = tmp116 * tmp235;
  const REAL tmp238 = tmp13 * tmp74 * tmp98;
  const REAL tmp241 = tmp112 * tmp235 * tmp45;
  const REAL tmp26 = -2 * M_scale * tmp23 * tmp24 * tmp25 + tmp21;
  const REAL tmp51 = tmp29 * tmp50 + 1;
  const REAL tmp57 = tmp50 * tmp55 + 1;
  const REAL tmp88 = -2 * M_scale * tmp24 * tmp25 * tmp87 + tmp85;
  const REAL tmp100 = tmp95 + tmp99;
  const REAL tmp106 = tmp44 * tmp60;
  const REAL tmp118 = tmp31 * tmp50;
  const REAL tmp142 = tmp28 * tmp50;
  const REAL tmp152 = tmp151 * tmp55;
  const REAL tmp166 = tmp10 * tmp14 * tmp159 * tmp31 - tmp108 * tmp155 * tmp157 - tmp109 * tmp162 * tmp31 + tmp159 * tmp164 * tmp47;
  const REAL tmp180 = tmp157 * tmp28 * tmp31 * tmp33;
  const REAL tmp184 = tmp157 * tmp181;
  const REAL tmp187 = tmp151 * tmp29;
  const REAL tmp199 = -tmp104 * tmp129 * tmp155;
  const REAL tmp200 = tmp113 * tmp181 * tmp22;
  const REAL tmp201 = tmp113 * tmp155 * tmp86;
  const REAL tmp203 = tmp151 * tmp170 * tmp28 * tmp31;
  const REAL tmp205 = tmp155 * tmp28 * tmp95;
  const REAL tmp206 = tmp155 * tmp28 * tmp99;
  const REAL tmp211 = -tmp104 * tmp108 * tmp181;
  const REAL tmp219 = tmp137 * tmp215 - tmp218;
  const REAL tmp232 = tmp123 * tmp215 - tmp231;
  const REAL tmp239 = tmp168 + tmp237 + tmp238;
  const REAL tmp246 = tmp56 * tmp60 * tmp61;
  const REAL tmp58 = tmp51 * tmp57;
  const REAL tmp119 = tmp113 * tmp22 * tmp45 + tmp115 * tmp95 + tmp115 * tmp99 - tmp116 * tmp118;
  const REAL tmp131 = -tmp118 * tmp129;
  const REAL tmp143 = -tmp108 * tmp142 + tmp138 * tmp52 - tmp139 * tmp140 + tmp140 * tmp21;
  const REAL tmp146 = tmp113 * tmp45 * tmp86 - tmp116 * tmp142 + tmp140 * tmp95 + tmp140 * tmp99;
  const REAL tmp154 = tmp118 * (tmp124 * tmp137 * tmp22 + tmp15) - tmp139 * tmp69 - tmp152 * tmp23 + tmp21 * tmp69;
  const REAL tmp171 = tmp112 * tmp155 * tmp168 * x + tmp152 * tmp170 + tmp69 * tmp95 + tmp69 * tmp99;
  const REAL tmp174 = -tmp152 * tmp87;
  const REAL tmp179 = tmp118 * (2 * a_spin + tmp123 * tmp124 * tmp22);
  const REAL tmp188 = -tmp187 * tmp23;
  const REAL tmp191 = tmp142 * (-2 * a_spin + 2 * tmp124 * tmp137 * y);
  const REAL tmp195 = tmp104 * tmp126 * tmp45;
  const REAL tmp210 = tmp104 * tmp138 * tmp45;
  const REAL tmp222 = -tmp134 * tmp65 + tmp142 * (tmp123 * tmp124 * tmp86 + tmp15) - tmp187 * tmp87 + tmp65 * tmp85;
  const REAL tmp226 = -tmp129 * tmp184 - tmp133 * tmp162 * tmp28 + tmp14 * tmp159 * tmp28 * tmp84 + tmp159 * tmp225 * tmp47;
  const REAL tmp227 = tmp112 * tmp168 * tmp181 * y + tmp170 * tmp187 + tmp65 * tmp95 + tmp65 * tmp99;
  const REAL tmp243 =
      -tmp103 * tmp116 * tmp155 * tmp74 + tmp115 * tmp237 + tmp115 * tmp238 - tmp137 * tmp215 + tmp14 * tmp241 * x + tmp168 * tmp31 * tmp45 + tmp218;
  const REAL tmp244 =
      -tmp103 * tmp116 * tmp181 * tmp74 - tmp123 * tmp215 + tmp14 * tmp241 * y + tmp140 * tmp237 + tmp140 * tmp238 + tmp168 * tmp28 * tmp45 + tmp231;
  const REAL tmp64 = tmp57 * tmp60 * tmp61;
  const REAL tmp68 = tmp51 * tmp60 * tmp61;
  const REAL tmp80 = -pow(tmp13, 9.0 / 2.0) * tmp72 + tmp43 * tmp56 * tmp60 * tmp71 + tmp51 * tmp60 * tmp79 - tmp56 * tmp58 + tmp57 * tmp60 * tmp78;
  const REAL tmp120 = -2 * M_scale * tmp108 * tmp18 * z + tmp111 + tmp119;
  const REAL tmp144 = tmp115 * tmp134 - tmp115 * tmp85 - tmp126 * tmp52 - tmp131 + tmp143;
  const REAL tmp147 = -2 * M_scale * tmp129 * tmp18 * z + tmp145 + tmp146;
  const REAL tmp172 = -4 * M_scale * tmp108 * tmp13 * tmp18 * tmp31 * tmp46 * z + 4 * M_scale * tmp108 * tmp18 * tmp31 * tmp45 * z +
                      4 * M_scale * tmp13 * tmp164 * tmp18 * tmp45 * z - tmp109 * tmp115 * tmp16 * tmp74 - tmp171;
  const REAL tmp185 = tmp10 * tmp181 * tmp182 * tmp31 - 4 * tmp109 * tmp181 * tmp24 * tmp31 + tmp134 * tmp69 + tmp138 * tmp155 * tmp157 +
                      tmp164 * tmp184 - tmp174 - tmp179 - tmp180 * tmp23 - tmp69 * tmp85;
  const REAL tmp186 = -tmp134 * tmp69 + tmp174 + tmp179 + tmp69 * tmp85;
  const REAL tmp192 = -tmp139 * tmp65 + tmp188 + tmp191 + tmp21 * tmp65;
  const REAL tmp193 = -tmp115 * tmp134 + tmp115 * tmp85 + tmp126 * tmp52 + tmp131 + tmp143;
  const REAL tmp207 = tmp115 * tmp145 - tmp115 * tmp197 - tmp195 - tmp199 + tmp200 + tmp201 + tmp203 + tmp205 + tmp206;
  const REAL tmp213 = tmp111 * tmp140 - tmp140 * tmp208 - tmp210 - tmp211;
  const REAL tmp220 = -tmp111 + tmp119 + tmp208;
  const REAL tmp228 = -4 * M_scale * tmp129 * tmp13 * tmp18 * tmp28 * tmp46 * z + 4 * M_scale * tmp129 * tmp18 * tmp28 * tmp45 * z +
                      4 * M_scale * tmp13 * tmp18 * tmp225 * tmp45 * z - tmp133 * tmp140 * tmp16 * tmp74 - tmp227;
  const REAL tmp229 = tmp126 * tmp184 - 4 * tmp133 * tmp155 * tmp24 * tmp28 + tmp139 * tmp65 + tmp155 * tmp157 * tmp225 +
                      tmp155 * tmp182 * tmp28 * tmp84 - tmp180 * tmp87 - tmp188 - tmp191 - tmp21 * tmp65;
  const REAL tmp233 = -tmp145 + tmp146 + tmp197;
  const REAL tmp214 = -tmp207 - tmp213;
  const REAL tmp221 = -tmp111 * tmp140 + tmp140 * tmp208 + tmp207 + tmp210 + tmp211;
  const REAL tmp234 = -tmp115 * tmp145 + tmp115 * tmp197 + tmp195 + tmp199 + tmp200 + tmp201 + tmp203 + tmp205 + tmp206 + tmp213;
  const REAL tmp81 = (1.0 / 2.0) / (16 * ((M_scale) * (M_scale) * (M_scale) * (M_scale)) * ((tmp13) * (tmp13) * (tmp13) * (tmp13) * (tmp13)) * tmp2 *
                                        tmp29 * tmp54 * tmp55 / ((tmp17) * (tmp17) * (tmp17) * (tmp17)) +
                                    4 * M_scale * tmp13 * tmp18 * z *
                                        (tmp24 * tmp51 * tmp59 * tmp69 * tmp74 + tmp24 * tmp57 * tmp59 * tmp65 * tmp74 - tmp72 * tmp73) +
                                    16 * tmp2 * tmp29 * tmp36 * tmp37 * tmp38 * tmp54 * tmp55 * tmp56 - tmp44 * tmp80 - tmp56 * tmp64 * tmp65 -
                                    tmp56 * tmp68 * tmp69 - tmp58 * tmp60 * tmp61);
  const REAL tmp82 = tmp81 * (4 * tmp16 * tmp2 * tmp25 * tmp31 * tmp42 * tmp45 * tmp51 + 4 * tmp25 * tmp29 * tmp31 * tmp33 * tmp42 * tmp43 * tmp44 -
                              tmp31 * tmp44 * tmp51 * tmp52 - tmp35 * tmp41);
  const REAL tmp92 = tmp81 * (4 * tmp16 * tmp2 * tmp25 * tmp28 * tmp42 * tmp45 * tmp57 + 4 * tmp25 * tmp28 * tmp33 * tmp42 * tmp43 * tmp44 * tmp55 -
                              tmp28 * tmp33 * tmp41 * tmp55 - tmp28 * tmp44 * tmp52 * tmp57);
  const REAL tmp105 = tmp81 * (-tmp102 * tmp71 - tmp104 * tmp58 + 4 * tmp24 * tmp25 * tmp29 * tmp42 * tmp46 * tmp57 * z +
                               4 * tmp24 * tmp25 * tmp42 * tmp46 * tmp51 * tmp55 * z);
  const REAL tmp107 = tmp81 * (-tmp106 * tmp43 * tmp71 + 16 * tmp2 * tmp29 * tmp36 * tmp37 * tmp38 * tmp54 * tmp55 + tmp44 * tmp51 * tmp57 -
                               tmp64 * tmp65 - tmp68 * tmp69);
  const REAL tmp245 =
      tmp81 * (-tmp118 * tmp28 * tmp44 * tmp56 - tmp155 * tmp28 * tmp41 + 4 * tmp16 * tmp2 * tmp25 * tmp28 * tmp31 * tmp42 * tmp46 * tmp56 +
               4 * tmp25 * tmp28 * tmp31 * tmp42 * tmp43 * tmp44 * tmp46);
  const REAL tmp247 = tmp81 * (-tmp106 * tmp78 + 16 * tmp2 * tmp29 * tmp36 * tmp37 * tmp38 * tmp46 - tmp246 * tmp65 + tmp44 * tmp51 * tmp56 - tmp68);
  const REAL tmp249 = tmp81 * (-tmp102 * tmp35 - tmp104 * tmp115 * tmp51 * tmp56 + 4 * tmp24 * tmp25 * tmp29 * tmp31 * tmp33 * tmp42 * tmp56 * z +
                               4 * tmp24 * tmp25 * tmp31 * tmp42 * tmp45 * tmp51 * z);
  const REAL tmp250 = tmp81 * (-tmp106 * tmp79 + 16 * tmp2 * tmp36 * tmp37 * tmp38 * tmp46 * tmp55 - tmp246 * tmp69 + tmp44 * tmp56 * tmp57 - tmp64);
  const REAL tmp251 = tmp81 * (-tmp102 * tmp28 * tmp33 * tmp55 - tmp104 * tmp140 * tmp56 * tmp57 +
                               4 * tmp24 * tmp25 * tmp28 * tmp33 * tmp42 * tmp55 * tmp56 * z + 4 * tmp24 * tmp25 * tmp28 * tmp42 * tmp45 * tmp57 * z);
  const REAL tmp252 = -tmp80 * tmp81;
  Gamma_local[0] = -tmp100 * tmp105 - tmp26 * tmp82 - tmp88 * tmp92;
  Gamma_local[1] = -tmp105 * tmp120 + tmp107 * tmp26 + tmp144 * tmp92;
  Gamma_local[2] = -tmp105 * tmp147 + tmp107 * tmp88 - tmp144 * tmp82;
  Gamma_local[3] = tmp100 * tmp107 + tmp120 * tmp82 + tmp147 * tmp92;
  Gamma_local[4] = tmp105 * tmp172 + tmp107 * tmp166 + tmp154 * tmp82 + tmp185 * tmp92;
  Gamma_local[5] = tmp105 * tmp214 + tmp107 * tmp193 + tmp186 * tmp82 + tmp192 * tmp92;
  Gamma_local[6] = tmp105 * tmp219 + tmp107 * tmp220 + tmp171 * tmp82 + tmp221 * tmp92;
  Gamma_local[7] = tmp105 * tmp228 + tmp107 * tmp226 + tmp222 * tmp92 + tmp229 * tmp82;
  Gamma_local[8] = tmp105 * tmp232 + tmp107 * tmp233 + tmp227 * tmp92 + tmp234 * tmp82;
  Gamma_local[9] = tmp105 * tmp236 + tmp107 * tmp239 + tmp243 * tmp82 + tmp244 * tmp92;
  Gamma_local[10] = -tmp100 * tmp249 - tmp245 * tmp88 - tmp247 * tmp26;
  Gamma_local[11] = -tmp120 * tmp249 + tmp144 * tmp245 + tmp26 * tmp82;
  Gamma_local[12] = -tmp144 * tmp247 - tmp147 * tmp249 + tmp82 * tmp88;
  Gamma_local[13] = tmp100 * tmp82 + tmp120 * tmp247 + tmp147 * tmp245;
  Gamma_local[14] = tmp154 * tmp247 + tmp166 * tmp82 + tmp172 * tmp249 + tmp185 * tmp245;
  Gamma_local[15] = tmp186 * tmp247 + tmp192 * tmp245 + tmp193 * tmp82 + tmp214 * tmp249;
  Gamma_local[16] = tmp171 * tmp247 + tmp219 * tmp249 + tmp220 * tmp82 + tmp221 * tmp245;
  Gamma_local[17] = tmp222 * tmp245 + tmp226 * tmp82 + tmp228 * tmp249 + tmp229 * tmp247;
  Gamma_local[18] = tmp227 * tmp245 + tmp232 * tmp249 + tmp233 * tmp82 + tmp234 * tmp247;
  Gamma_local[19] = tmp236 * tmp249 + tmp239 * tmp82 + tmp243 * tmp247 + tmp244 * tmp245;
  Gamma_local[20] = -tmp100 * tmp251 - tmp245 * tmp26 - tmp250 * tmp88;
  Gamma_local[21] = -tmp120 * tmp251 + tmp144 * tmp250 + tmp26 * tmp92;
  Gamma_local[22] = -tmp144 * tmp245 - tmp147 * tmp251 + tmp88 * tmp92;
  Gamma_local[23] = tmp100 * tmp92 + tmp120 * tmp245 + tmp147 * tmp250;
  Gamma_local[24] = tmp154 * tmp245 + tmp166 * tmp92 + tmp172 * tmp251 + tmp185 * tmp250;
  Gamma_local[25] = tmp186 * tmp245 + tmp192 * tmp250 + tmp193 * tmp92 + tmp214 * tmp251;
  Gamma_local[26] = tmp171 * tmp245 + tmp219 * tmp251 + tmp220 * tmp92 + tmp221 * tmp250;
  Gamma_local[27] = tmp222 * tmp250 + tmp226 * tmp92 + tmp228 * tmp251 + tmp229 * tmp245;
  Gamma_local[28] = tmp227 * tmp250 + tmp232 * tmp251 + tmp233 * tmp92 + tmp234 * tmp245;
  Gamma_local[29] = tmp236 * tmp251 + tmp239 * tmp92 + tmp243 * tmp245 + tmp244 * tmp250;
  Gamma_local[30] = -tmp100 * tmp252 - tmp249 * tmp26 - tmp251 * tmp88;
  Gamma_local[31] = tmp105 * tmp26 - tmp120 * tmp252 + tmp144 * tmp251;
  Gamma_local[32] = tmp105 * tmp88 - tmp144 * tmp249 - tmp147 * tmp252;
  Gamma_local[33] = tmp100 * tmp105 + tmp120 * tmp249 + tmp147 * tmp251;
  Gamma_local[34] = tmp105 * tmp166 + tmp154 * tmp249 + tmp172 * tmp252 + tmp185 * tmp251;
  Gamma_local[35] = tmp105 * tmp193 + tmp186 * tmp249 + tmp192 * tmp251 + tmp214 * tmp252;
  Gamma_local[36] = tmp105 * tmp220 + tmp171 * tmp249 + tmp219 * tmp252 + tmp221 * tmp251;
  Gamma_local[37] = tmp105 * tmp226 + tmp222 * tmp251 + tmp228 * tmp252 + tmp229 * tmp249;
  Gamma_local[38] = tmp105 * tmp233 + tmp227 * tmp251 + tmp232 * tmp252 + tmp234 * tmp249;
  Gamma_local[39] = tmp105 * tmp239 + tmp236 * tmp252 + tmp243 * tmp249 + tmp244 * tmp251;

} // END FUNCTION: connections_KerrSchild_Cartesian

/**
 * Interpolate piecewise analytic/numerical spacetime tensors for one photon chunk.
 *
 * The caller supplies an active numerical time window, a spatial interpolation
 * context, and one chunk of photon states in the same Structure-of-Arrays bundle
 * layout used by the analytic geodesic interpolation kernel. This CPU wrapper
 * parallelizes over rays, evaluates `KerrSchild_Cartesian` directly for times
 * below `t_metric_0`, evaluates `KerrSchild_Cartesian` directly for times above
 * `t_metric_1`, and otherwise reconstructs one mixed temporal stencil by
 * combining spatial interpolation on the mapped numerical slices with one static
 * analytic fill on the missing stencil edge before performing temporal
 * interpolation at the photon coordinate time.
 *
 * The design goal is to let all photons in the chunk reuse the same mapped
 * numerical-spacetime payload window rather than loading numerical grids
 * independently ray-by-ray.
 *
 * @param[in] commondata Common runtime parameters.
 * @param[in] params Generated BHaH grid parameters for the mapped numerical data.
 * @param[in] spatial_context Trusted azimuthal-symmetry spatial interpolation context.
 * @param[in] numerical_window Active mapped numerical time-window manager.
 * @param[in] d_f_bundle Photon state bundle.
 * @param[out] d_metric_bundle Destination metric bundle.
 * @param[out] d_connection_bundle Destination Christoffel bundle, or NULL.
 * @param chunk_size Number of active rays in the chunk.
 * @param stream_idx Analytic-kernel compatibility argument; ignored on CPU.
 *
 * @pre `numerical_window` must already map all slices needed by the active chunk.
 */
void numerical_interpolation(const commondata_struct *restrict commondata, const params_struct *restrict params,
                             const azimuthal_symmetry_spatial_lagrange_context_struct *restrict spatial_context,
                             const NumericalTimeWindowManager *restrict numerical_window, const double *restrict d_f_bundle,
                             double *restrict d_metric_bundle, double *restrict d_connection_bundle, const long int chunk_size,
                             const int stream_idx) {
  (void)stream_idx;

#define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
#define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
#define IDX_CONN(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
#define G4_SLICE(base_ptr, slice_idx, comp_idx) ((base_ptr)[(slice_idx) * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + (comp_idx)])
#define GAMMA_SLICE(base_ptr, slice_idx, comp_idx) ((base_ptr)[(slice_idx) * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT + (comp_idx)])

  const int temporal_half_width = commondata->numerical_spacetime_temporal_interp_order;
  const REAL t_metric_0 = (REAL)commondata->t_metric_0;
  const REAL t_metric_1 = (REAL)commondata->t_metric_1;
  // The mapped numerical window and the temporal helper must agree on the
  // centered temporal stencil width before any variable-length arrays are
  // sized from runtime data.
  if (temporal_half_width < 0 || temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH ||
      temporal_half_width != numerical_window->temporal_interp_half_width || !isfinite((double)t_metric_0) || !isfinite((double)t_metric_1) ||
      t_metric_0 >= t_metric_1) {
#pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_connection_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++)
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
    } // END LOOP: for i over rays after invalid temporal order
    return;
  } // END IF: runtime temporal interpolation half-width was invalid
  const int temporal_num_points = 2 * temporal_half_width + 1;
  if (temporal_num_points != numerical_window->temporal_interp_num_points) {
#pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_connection_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++)
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
    } // END LOOP: for i over rays after inconsistent temporal stencil size
    return;
  } // END IF: runtime temporal stencil size did not match the mapped numerical window

#pragma omp parallel for
  for (long int i = 0; i < chunk_size; i++) {
    double f_local[9];
    for (int comp = 0; comp < 9; comp++)
      f_local[comp] = d_f_bundle[IDX_F(comp, i)];
    const REAL t = (REAL)f_local[0];
    const REAL x = (REAL)f_local[1];
    const REAL y = (REAL)f_local[2];
    const REAL z = (REAL)f_local[3];
    int ray_failed = 0;
    uint64_t available_slice_indices[temporal_num_points];
    REAL available_slice_times[temporal_num_points];
    const double *available_slice_payloads[temporal_num_points];
    int num_available_slices = 0;
    REAL missing_slice_times[temporal_num_points];
    int num_missing_slices = 0;
    REAL g4dd_available[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_available[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];
    REAL full_slice_times[temporal_num_points];
    REAL g4dd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];
    REAL g4dd_local[TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_local[TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];
    REAL g4dd_missing_local[TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT] = {0};
    REAL gamma4udd_missing_local[TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT] = {0};

    // Step 1: Dispatch directly to the static analytic metrics when this ray
    // is fully outside the numerical time interval.
    if (t < t_metric_0) {
      g4DD_metric_KerrSchild_Cartesian(commondata, f_local, g4dd_local);
      if (d_connection_bundle != NULL)
        connections_KerrSchild_Cartesian(commondata, f_local, gamma4udd_local);
    } else if (t > t_metric_1) {
      g4DD_metric_KerrSchild_Cartesian(commondata, f_local, g4dd_local);
      if (d_connection_bundle != NULL)
        connections_KerrSchild_Cartesian(commondata, f_local, gamma4udd_local);
    } else {
      // Step 2: Recover one adaptive numerical stencil from the slot-level
      // time window shared by the whole chunk.
      const int window_status = time_window_manager_numerical_stencil_for_time(
          numerical_window, (double)t, temporal_half_width, available_slice_indices, available_slice_times, available_slice_payloads,
          &num_available_slices, missing_slice_times, &num_missing_slices);
      if (window_status != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS || num_available_slices <= 0 ||
          num_available_slices + num_missing_slices != temporal_num_points) {
        ray_failed = 1;
      } else {
        // Step 3: Interpolate only the available mapped numerical slices in
        // space at the photon position.
        const int spatial_status = azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical(
            spatial_context, commondata, params, x, y, z, num_available_slices, available_slice_payloads, g4dd_available, gamma4udd_available);
        if (spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
          ray_failed = 1;
        else {
          int missing_is_upper = 0;
          int missing_is_lower = 0;

          // Step 4: If the centered stencil is only partially mapped, fill
          // the missing edge with one static analytic metric evaluation.
          if (num_missing_slices > 0) {
            missing_is_upper = missing_slice_times[0] > available_slice_times[num_available_slices - 1];
            missing_is_lower = missing_slice_times[num_missing_slices - 1] < available_slice_times[0];

            if (missing_is_upper == missing_is_lower) {
              ray_failed = 1;
            } else {
              f_local[0] = (double)missing_slice_times[0];
              if (missing_is_upper) {
                g4DD_metric_KerrSchild_Cartesian(commondata, f_local, g4dd_missing_local);
                connections_KerrSchild_Cartesian(commondata, f_local, gamma4udd_missing_local);
              } else {
                g4DD_metric_KerrSchild_Cartesian(commondata, f_local, g4dd_missing_local);
                connections_KerrSchild_Cartesian(commondata, f_local, gamma4udd_missing_local);
              } // END ELSE: lower missing stencil edge uses analytical_metric_0
              f_local[0] = (double)t;
            } // END ELSE: missing stencil edge classification was usable
          } // END IF: one analytic stencil edge must be synthesized

          if (!ray_failed) {
            // Step 5: Reconstruct the full ordered temporal stencil expected
            // by the temporal interpolation helper.
            if (num_missing_slices > 0 && missing_is_lower) {
              for (int s = 0; s < num_missing_slices; s++) {
                full_slice_times[s] = missing_slice_times[s];
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, s, comp) = g4dd_missing_local[comp];
                } // END LOOP: for comp over missing analytic metric components
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, s, comp) = gamma4udd_missing_local[comp];
                } // END LOOP: for comp over missing analytic Christoffel components
              } // END LOOP: for s over lower missing stencil nodes
              for (int s = 0; s < num_available_slices; s++) {
                const int full_slot = num_missing_slices + s;
                full_slice_times[full_slot] = available_slice_times[s];
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, full_slot, comp) = G4_SLICE(g4dd_available, s, comp);
                } // END LOOP: for comp over mapped numerical metric components after lower analytic padding
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, full_slot, comp) = GAMMA_SLICE(gamma4udd_available, s, comp);
                } // END LOOP: for comp over mapped numerical Christoffel components after lower analytic padding
              } // END LOOP: for s over available numerical stencil nodes after lower analytic padding
            } else {
              for (int s = 0; s < num_available_slices; s++) {
                full_slice_times[s] = available_slice_times[s];
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, s, comp) = G4_SLICE(g4dd_available, s, comp);
                } // END LOOP: for comp over mapped numerical metric components
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, s, comp) = GAMMA_SLICE(gamma4udd_available, s, comp);
                } // END LOOP: for comp over mapped numerical Christoffel components
              } // END LOOP: for s over available numerical stencil nodes
              for (int s = 0; s < num_missing_slices; s++) {
                const int full_slot = num_available_slices + s;
                full_slice_times[full_slot] = missing_slice_times[s];
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, full_slot, comp) = g4dd_missing_local[comp];
                } // END LOOP: for comp over upper missing analytic metric components
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
                  GAMMA_SLICE(gamma4udd_slices, full_slot, comp) = gamma4udd_missing_local[comp];
                } // END LOOP: for comp over upper missing analytic Christoffel components
              } // END LOOP: for s over upper missing stencil nodes
            } // END ELSE: upper missing edge or fully numerical centered stencil
          } // END IF: full ordered stencil reconstruction remained valid

          if (!ray_failed) {
            // Step 6: Interpolate the reconstructed per-slice tensor bundles
            // in physical time to the photon coordinate time.
            const int temporal_status =
                temporal_lagrange_interpolation(commondata, full_slice_times, g4dd_slices, gamma4udd_slices, t, g4dd_local, gamma4udd_local);
            if (temporal_status != TEMPORAL_LAGRANGE_INTERP_SUCCESS)
              ray_failed = 1;
          } // END IF: reconstructed stencil was ready for temporal interpolation
        } // END ELSE: spatial interpolation succeeded for the mapped numerical stencil subset
      } // END ELSE: adaptive stencil query succeeded for this mixed ray
    } // END ELSE: photon required numerical or mixed interpolation

    if (ray_failed) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_connection_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++)
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
      continue;
    } // END IF: at least one interpolation stage failed for this ray

    for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
      d_metric_bundle[IDX_METRIC(comp, i)] = (double)g4dd_local[comp];
    if (d_connection_bundle != NULL)
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++)
        d_connection_bundle[IDX_CONN(comp, i)] = (double)gamma4udd_local[comp];
  } // END LOOP: for i over rays in chunk

#undef IDX_F
#undef IDX_METRIC
#undef IDX_CONN
#undef G4_SLICE
#undef GAMMA_SLICE
} // END FUNCTION: numerical_interpolation
