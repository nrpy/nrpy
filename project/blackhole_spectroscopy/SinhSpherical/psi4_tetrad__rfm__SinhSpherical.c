#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/*
 * Compute quasiKinnersley tetrad for psi4, with use_metric_to_construct_unit_normal=False
 */
void psi4_tetrad__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL cf,
                                     const REAL hDD00, const REAL hDD01, const REAL hDD02, const REAL hDD11, const REAL hDD12, const REAL hDD22,
                                     REAL *mre4U0, REAL *mre4U1, REAL *mre4U2, REAL *mre4U3, REAL *mim4U0, REAL *mim4U1, REAL *mim4U2, REAL *mim4U3,
                                     REAL *n4U0, REAL *n4U1, REAL *n4U2, REAL *n4U3, REAL *restrict xx[3], const int i0, const int i1, const int i2) {
#include "../set_CodeParameters.h"
  const REAL xx0 = xx[0][i0];
  const REAL xx1 = xx[1][i1];
  const REAL xx2 = xx[2][i2];
  // Compute tetrads:
  const REAL tmp0 = sin(xx2);
  const REAL tmp2 = (1.0 / (SINHW));
  const REAL tmp17 = cos(xx2);
  const REAL tmp26 = sin(xx1);
  const REAL tmp33 = cos(xx1);
  const REAL tmp45 = 2 * hDD01;
  const REAL tmp62 = hDD00 + 1;
  const REAL tmp67 = hDD22 + 1;
  const REAL tmp70 = hDD11 + 1;
  const REAL tmp1 = ((tmp0) * (tmp0));
  const REAL tmp3 = 2 * tmp2;
  const REAL tmp8 = exp(tmp2) - exp(-tmp2);
  const REAL tmp25 = (1.0 / ((cf) * (cf)));
  const REAL tmp27 = ((tmp26) * (tmp26));
  const REAL tmp38 = tmp26 * tmp33;
  const REAL tmp61 = (1.0 / (tmp26));
  const REAL tmp66 = hDD02 * hDD12 * tmp45;
  const REAL tmp68 = ((hDD01) * (hDD01)) * tmp67;
  const REAL tmp71 = ((hDD02) * (hDD02)) * tmp70;
  const REAL tmp73 = ((hDD12) * (hDD12)) * tmp62;
  const REAL tmp4 = exp(tmp3 * xx0);
  const REAL tmp6 = exp(tmp3);
  const REAL tmp10 = exp(tmp2 * xx0);
  const REAL tmp11 = exp(-tmp2 * xx0);
  const REAL tmp28 = ((AMPL) * (AMPL)) / ((tmp8) * (tmp8));
  const REAL tmp75 = tmp62 * tmp67 * tmp70;
  const REAL tmp5 = (1.0 / (tmp4 - 1));
  const REAL tmp7 = tmp6 - 1;
  const REAL tmp12 = tmp10 - tmp11;
  const REAL tmp19 = 1 - tmp4;
  const REAL tmp20 = 1 - tmp6;
  const REAL tmp42 = tmp10 * tmp2 + tmp11 * tmp2;
  const REAL tmp14 = tmp12 * exp(tmp2 * (xx0 - 1)) / tmp8;
  const REAL tmp21 = tmp20 / tmp19;
  const REAL tmp29 = ((tmp12) * (tmp12)) * tmp28;
  const REAL tmp35 = SINHW / (tmp4 + 1);
  const REAL tmp43 = tmp28 * ((tmp42) * (tmp42));
  const REAL tmp47 = tmp12 * tmp28 * tmp42;
  const REAL tmp64 = ((tmp4 + 1) * (tmp4 + 1));
  const REAL tmp18 = tmp14 * ((tmp17) * (tmp17));
  const REAL tmp41 = tmp25 * (hDD11 * tmp29 + tmp29);
  const REAL tmp44 = tmp25 * (hDD00 * tmp43 + tmp43);
  const REAL tmp53 = hDD12 * tmp29;
  const REAL tmp16 = tmp14 * tmp5 * tmp7;
  const REAL tmp31 = tmp25 * (hDD22 * tmp27 * tmp29 + tmp27 * tmp29);
  const REAL tmp36 = tmp14 * tmp35 * tmp7;
  const REAL tmp23 = tmp1 * tmp16 + tmp18 * tmp21;
  const REAL tmp37 = tmp1 * tmp27 * tmp36 + tmp18 * tmp27 * tmp35 * tmp7 + ((tmp33) * (tmp33)) * tmp36;
  const REAL tmp40 = tmp1 * tmp14 * tmp21 * tmp38 - tmp16 * tmp38 + tmp18 * tmp21 * tmp38;
  const REAL tmp50 = tmp0 * tmp14 * tmp17 * tmp21 - tmp0 * tmp16 * tmp17;
  const REAL tmp32 = tmp23 / sqrt(((tmp23) * (tmp23)) * tmp31);
  const REAL tmp46 = tmp25 * tmp40;
  const REAL tmp48 = tmp37 * tmp47;
  const REAL tmp76 = ((cf) * (cf)) * tmp23 *
                     sqrt(((tmp19) * (tmp19) * (tmp19) * (tmp19)) * tmp27 * tmp64 * (tmp66 - tmp68 - tmp71 - tmp73 + tmp75) *
                          exp(6 * tmp2 * (1 - xx0)) / (((SINHW) * (SINHW)) * pow(cf, 6) * pow(tmp20, 6))) *
                     exp(tmp3 * (xx0 - 1)) * fabs(AMPL) / (-tmp66 + tmp68 + tmp71 + tmp73 - tmp75);
  const REAL tmp52 = tmp26 * tmp32;
  const REAL tmp54 = tmp46 * tmp53;
  const REAL tmp56 = hDD02 * tmp25 * tmp48;
  const REAL tmp77 = ((tmp20) * (tmp20)) * tmp76;
  const REAL tmp83 = tmp35 * tmp5 * ((tmp7) * (tmp7)) * tmp76 * (hDD01 * tmp67 - hDD02 * hDD12);
  const REAL tmp57 = -tmp32 * (tmp31 * tmp32 * tmp50 + tmp52 * tmp54 + tmp52 * tmp56) + tmp50;
  const REAL tmp78 = tmp37 * tmp77 / ((tmp19) * (tmp19));
  const REAL tmp87 = -((SINHW) * (SINHW)) * tmp40 * tmp77 * (((hDD12) * (hDD12)) - tmp67 * tmp70) / tmp64 + tmp37 * tmp83;
  const REAL tmp58 = tmp26 * tmp57;
  const REAL tmp80 =
      tmp35 * tmp40 * tmp5 * tmp61 * ((tmp7) * (tmp7)) * tmp76 * (hDD01 * hDD12 - hDD02 * tmp70) + tmp61 * tmp78 * (-hDD01 * hDD02 + hDD12 * tmp62);
  const REAL tmp84 = -tmp40 * tmp83 + tmp78 * (((hDD02) * (hDD02)) - tmp62 * tmp67);
  const REAL tmp60 = (1.0 / sqrt(tmp31 * ((tmp57) * (tmp57)) + ((tmp37) * (tmp37)) * tmp44 + ((tmp40) * (tmp40)) * tmp41 + tmp45 * tmp46 * tmp48 +
                                 2 * tmp54 * tmp58 + 2 * tmp56 * tmp58));
  const REAL tmp99 = M_SQRT1_2 * tmp60;
  const REAL tmp82 = tmp26 * tmp60 * tmp80;
  const REAL tmp86 = tmp25 * tmp60 * tmp84;
  const REAL tmp89 = tmp47 * tmp60 * tmp87;
  const REAL tmp90 = tmp60 * (hDD01 * tmp46 * tmp89 + hDD01 * tmp48 * tmp86 + hDD02 * tmp25 * tmp58 * tmp89 + tmp31 * tmp57 * tmp60 * tmp80 +
                              tmp37 * tmp44 * tmp60 * tmp87 + tmp40 * tmp41 * tmp60 * tmp84 + tmp53 * tmp58 * tmp86 + tmp54 * tmp82 + tmp56 * tmp82);
  const REAL tmp91 = -tmp37 * tmp90 + tmp87;
  const REAL tmp92 = -tmp40 * tmp90 + tmp84;
  const REAL tmp95 = -tmp32 * (hDD02 * tmp25 * tmp47 * tmp52 * tmp87 + tmp25 * tmp52 * tmp53 * tmp84 + tmp31 * tmp32 * tmp80) - tmp57 * tmp90 + tmp80;
  const REAL tmp96 = tmp25 * tmp47 * tmp91;
  const REAL tmp97 = 2 * tmp26 * tmp95;
  const REAL tmp98 = M_SQRT1_2 / sqrt(hDD02 * tmp96 * tmp97 + tmp25 * tmp53 * tmp92 * tmp97 + tmp31 * ((tmp95) * (tmp95)) +
                                      tmp41 * ((tmp92) * (tmp92)) + tmp44 * ((tmp91) * (tmp91)) + tmp45 * tmp92 * tmp96);
  *mim4U0 = 0;
  *mim4U1 = 0;
  *mim4U2 = 0;
  *mim4U3 = M_SQRT1_2 * tmp32;
  *mre4U0 = 0;
  *mre4U1 = tmp91 * tmp98;
  *mre4U2 = tmp92 * tmp98;
  *mre4U3 = tmp95 * tmp98;
  *n4U0 = M_SQRT1_2;
  *n4U1 = -tmp37 * tmp99;
  *n4U2 = -tmp40 * tmp99;
  *n4U3 = -tmp57 * tmp99;
}
