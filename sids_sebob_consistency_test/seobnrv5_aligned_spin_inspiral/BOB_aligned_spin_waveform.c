#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Calculate the BOB 22 mode.
 */
void BOB_aligned_spin_waveform(const REAL t, commondata_struct *restrict commondata, REAL *restrict waveform) {

  const REAL m1 = commondata->m1;
  const REAL m2 = commondata->m2;
  const REAL chi1 = commondata->chi1;
  const REAL chi2 = commondata->chi2;
  const REAL omega_qnm = commondata->omega_qnm;
  const REAL tau_qnm = commondata->tau_qnm;
  const REAL t_0 = commondata->t_attach;
  // compute
  const REAL tmp0 = m1 + m2;
  const REAL tmp2 = ((m1) * (m1));
  const REAL tmp3 = ((m2) * (m2));
  const REAL tmp21 = (1.0 / (omega_qnm));
  const REAL tmp28 = (1.0 / 2.0) * omega_qnm;
  const REAL tmp1 = m1 * m2 / ((tmp0) * (tmp0));
  const REAL tmp4 = tmp2 * tmp3 / ((tmp0) * (tmp0) * (tmp0) * (tmp0));
  const REAL tmp5 = ((m1) * (m1) * (m1)) * ((m2) * (m2) * (m2)) / pow(tmp0, 6);
  const REAL tmp6 = ((m1) * (m1) * (m1) * (m1)) * ((m2) * (m2) * (m2) * (m2)) / pow(tmp0, 8);
  const REAL tmp7 = tmp2 + tmp3;
  const REAL tmp8 = chi1 * tmp2 + chi2 * tmp3;
  const REAL tmp22 = 2 * tmp21;
  const REAL tmp9 = tmp8 / tmp7;
  const REAL tmp10 = ((tmp8) * (tmp8) * (tmp8) * (tmp8)) / ((tmp7) * (tmp7) * (tmp7) * (tmp7));
  const REAL tmp11 = ((tmp8) * (tmp8) * (tmp8)) / ((tmp7) * (tmp7) * (tmp7));
  const REAL tmp12 = ((tmp8) * (tmp8)) / ((tmp7) * (tmp7));
  const REAL tmp13 = tmp1 * tmp9;
  const REAL tmp14 = tmp1 * tmp11;
  const REAL tmp16 = tmp4 * tmp9;
  const REAL tmp17 = tmp5 * tmp9;
  const REAL tmp18 = tmp12 * tmp4;
  const REAL tmp19 = -647517837725651.0 / 2500000000000000.0 * tmp1 * tmp12 + (31709602351033533.0 / 200000000000000000.0) * tmp1 +
                     (381474417039561.0 / 50000000000000000.0) * tmp10 + (1673164620878479.0 / 50000000000000000.0) * tmp11 +
                     (871517604568457.0 / 20000000000000000.0) * tmp12 - 392132032264821.0 / 3125000000000000.0 * tmp13 -
                     24194837236629313.0 / 200000000000000000.0 * tmp14 + (16973430239436997.0 / 20000000000000000.0) * tmp16 -
                     7502911609839309.0 / 4000000000000000.0 * tmp17 + (9714093262519423.0 / 20000000000000000.0) * tmp18 -
                     7134801214011141.0 / 50000000000000000.0 * tmp4 + (6698610724189451.0 / 4000000000000000.0) * tmp5 -
                     5893523296177077.0 / 2000000000000000.0 * tmp6 + (4566965972049467.0 / 100000000000000000.0) * tmp9 +
                     1342707196092513.0 / 10000000000000000.0;
  const REAL tmp20 = ((tmp19) * (tmp19) * (tmp19) * (tmp19));
  const REAL tmp23 = tmp19 * tmp22;
  const REAL tmp31 = ((tmp19) * (tmp19));
  const REAL tmp35 = tmp28 / tmp19;
  const REAL tmp24 = 2 * log(tmp23);
  const REAL tmp33 = tau_qnm * tmp21 * tmp31;
  const REAL tmp25 = tanh(tmp24);
  const REAL tmp26 = (t - t_0 + tau_qnm * tmp24) / tau_qnm;
  const REAL tmp27 = tmp20 + ((1.0 / 16.0) * ((omega_qnm) * (omega_qnm) * (omega_qnm) * (omega_qnm)) - tmp20) * (-tmp25 + tanh(tmp26)) / (1 - tmp25);
  const REAL tmp29 = pow(tmp27, 0.25);
  const REAL tmp36 = tmp28 * tmp29 / tmp31;
  const REAL h = (1.0 / 4.0) * tmp1 * (1.0 / sqrt(tmp27)) *
                 ((-647517837725651.0 / 1250000000000000.0 * tmp1 * tmp12 + (31709602351033533.0 / 100000000000000000.0) * tmp1 +
                   (381474417039561.0 / 25000000000000000.0) * tmp10 + (1673164620878479.0 / 25000000000000000.0) * tmp11 +
                   (871517604568457.0 / 10000000000000000.0) * tmp12 - 392132032264821.0 / 1562500000000000.0 * tmp13 -
                   24194837236629313.0 / 100000000000000000.0 * tmp14 + (16973430239436997.0 / 10000000000000000.0) * tmp16 -
                   7502911609839309.0 / 2000000000000000.0 * tmp17 + (9714093262519423.0 / 10000000000000000.0) * tmp18 -
                   7134801214011141.0 / 25000000000000000.0 * tmp4 + (6698610724189451.0 / 2000000000000000.0) * tmp5 -
                   5893523296177077.0 / 1000000000000000.0 * tmp6 + (4566965972049467.0 / 50000000000000000.0) * tmp9 +
                   1342707196092513.0 / 5000000000000000.0) *
                  (-647517837725651.0 / 1250000000000000.0 * tmp1 * tmp12 + (31709602351033533.0 / 100000000000000000.0) * tmp1 +
                   (381474417039561.0 / 25000000000000000.0) * tmp10 + (1673164620878479.0 / 25000000000000000.0) * tmp11 +
                   (871517604568457.0 / 10000000000000000.0) * tmp12 - 392132032264821.0 / 1562500000000000.0 * tmp13 -
                   24194837236629313.0 / 100000000000000000.0 * tmp14 + (16973430239436997.0 / 10000000000000000.0) * tmp16 -
                   7502911609839309.0 / 2000000000000000.0 * tmp17 + (9714093262519423.0 / 10000000000000000.0) * tmp18 -
                   7134801214011141.0 / 25000000000000000.0 * tmp4 + (6698610724189451.0 / 2000000000000000.0) * tmp5 -
                   5893523296177077.0 / 1000000000000000.0 * tmp6 + (4566965972049467.0 / 50000000000000000.0) * tmp9 +
                   1342707196092513.0 / 5000000000000000.0)) *
                 cosh(tmp24) *
                 fabs(-217072339408107.0 / 250000000000000.0 * tmp1 - 8493901280736431.0 / 100000000000000000.0 * tmp11 -
                      10413107147836477.0 / 500000000000000000.0 * tmp12 - 217891933479267.0 / 125000000000000.0 * tmp13 +
                      (860293461381563.0 / 2000000000000000.0) * tmp14 + (7194264161892297.0 / 1000000000000000.0) * tmp16 -
                      6678807011156761.0 / 500000000000000.0 * tmp17 + (3099447225891283.0 / 5000000000000000.0) * tmp18 +
                      (12440404909323101.0 / 1000000000000000.0) * tmp4 - 4687585958426211.0 / 100000000000000.0 * tmp5 +
                      (3598984888018441.0 / 50000000000000.0) * tmp6 + (1869399114678491.0 / 10000000000000000.0) * tmp9 +
                      14670966347991181.0 / 10000000000000000.0) /
                 cosh(tmp26);
  const REAL phi = omega_qnm * tau_qnm * (-atan2(tmp19, tmp28) + atan2(tmp29, tmp28)) +
                   0.5 * omega_qnm * tau_qnm * log((1 - tmp23) * (tmp22 * tmp29 + 1) / ((tmp23 + 1) * (-tmp22 * tmp29 + 1))) -
                   4 * tmp33 * (-atan2(tmp19, tmp22 * tmp31) + atan2(tmp29, tmp22 * tmp31)) -
                   2.0 * tmp33 * log((1 - tmp35) * (tmp36 + 1) / ((1 - tmp36) * (tmp35 + 1)));

  waveform[0] = h;
  waveform[1] = phi;
} // END FUNCTION BOB_aligned_spin_waveform
