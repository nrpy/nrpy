#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Calculate the BOB informed NQC amplitudes and phases.
 */
void BOB_aligned_spin_NQC_rhs(commondata_struct *restrict commondata, REAL *restrict amps, REAL *restrict omegas) {

  const REAL m1 = commondata->m1;
  const REAL m2 = commondata->m2;
  const REAL chi1 = commondata->chi1;
  const REAL chi2 = commondata->chi2;
  const REAL omega_qnm = commondata->omega_qnm;
  const REAL tau_qnm = commondata->tau_qnm;
  // compute
  const REAL tmp0 = m1 + m2;
  const REAL tmp2 = ((m1) * (m1));
  const REAL tmp3 = ((m2) * (m2));
  const REAL tmp1 = m1 * m2 / ((tmp0) * (tmp0));
  const REAL tmp4 = tmp2 * tmp3 / ((tmp0) * (tmp0) * (tmp0) * (tmp0));
  const REAL tmp5 = ((m1) * (m1) * (m1)) * ((m2) * (m2) * (m2)) / pow(tmp0, 6);
  const REAL tmp6 = ((m1) * (m1) * (m1) * (m1)) * ((m2) * (m2) * (m2) * (m2)) / pow(tmp0, 8);
  const REAL tmp7 = tmp2 + tmp3;
  const REAL tmp8 = chi1 * tmp2 + chi2 * tmp3;
  const REAL tmp9 = tmp8 / tmp7;
  const REAL tmp10 = ((tmp8) * (tmp8) * (tmp8)) / ((tmp7) * (tmp7) * (tmp7));
  const REAL tmp11 = ((tmp8) * (tmp8)) / ((tmp7) * (tmp7));
  const REAL tmp18 = ((tmp8) * (tmp8) * (tmp8) * (tmp8)) / ((tmp7) * (tmp7) * (tmp7) * (tmp7));
  const REAL tmp12 = tmp1 * tmp9;
  const REAL tmp13 = tmp1 * tmp10;
  const REAL tmp14 = tmp4 * tmp9;
  const REAL tmp15 = tmp5 * tmp9;
  const REAL tmp16 = tmp11 * tmp4;
  const REAL tmp17 = tmp1 * fabs(-217072339408107.0 / 250000000000000.0 * tmp1 - 8493901280736431.0 / 100000000000000000.0 * tmp10 -
                                 10413107147836477.0 / 500000000000000000.0 * tmp11 - 217891933479267.0 / 125000000000000.0 * tmp12 +
                                 (860293461381563.0 / 2000000000000000.0) * tmp13 + (7194264161892297.0 / 1000000000000000.0) * tmp14 -
                                 6678807011156761.0 / 500000000000000.0 * tmp15 + (3099447225891283.0 / 5000000000000000.0) * tmp16 +
                                 (12440404909323101.0 / 1000000000000000.0) * tmp4 - 4687585958426211.0 / 100000000000000.0 * tmp5 +
                                 (3598984888018441.0 / 50000000000000.0) * tmp6 + (1869399114678491.0 / 10000000000000000.0) * tmp9 +
                                 14670966347991181.0 / 10000000000000000.0);
  const REAL tmp20 = -647517837725651.0 / 2500000000000000.0 * tmp1 * tmp11 + (31709602351033533.0 / 200000000000000000.0) * tmp1 +
                     (1673164620878479.0 / 50000000000000000.0) * tmp10 + (871517604568457.0 / 20000000000000000.0) * tmp11 -
                     392132032264821.0 / 3125000000000000.0 * tmp12 - 24194837236629313.0 / 200000000000000000.0 * tmp13 +
                     (16973430239436997.0 / 20000000000000000.0) * tmp14 - 7502911609839309.0 / 4000000000000000.0 * tmp15 +
                     (9714093262519423.0 / 20000000000000000.0) * tmp16 + (381474417039561.0 / 50000000000000000.0) * tmp18 -
                     7134801214011141.0 / 50000000000000000.0 * tmp4 + (6698610724189451.0 / 4000000000000000.0) * tmp5 -
                     5893523296177077.0 / 2000000000000000.0 * tmp6 + (4566965972049467.0 / 100000000000000000.0) * tmp9 +
                     1342707196092513.0 / 10000000000000000.0;
  const REAL tmp22 = -647517837725651.0 / 1250000000000000.0 * tmp1 * tmp11 + (31709602351033533.0 / 100000000000000000.0) * tmp1 +
                     (1673164620878479.0 / 25000000000000000.0) * tmp10 + (871517604568457.0 / 10000000000000000.0) * tmp11 -
                     392132032264821.0 / 1562500000000000.0 * tmp12 - 24194837236629313.0 / 100000000000000000.0 * tmp13 +
                     (16973430239436997.0 / 10000000000000000.0) * tmp14 - 7502911609839309.0 / 2000000000000000.0 * tmp15 +
                     (9714093262519423.0 / 10000000000000000.0) * tmp16 + (381474417039561.0 / 25000000000000000.0) * tmp18 -
                     7134801214011141.0 / 25000000000000000.0 * tmp4 + (6698610724189451.0 / 2000000000000000.0) * tmp5 -
                     5893523296177077.0 / 1000000000000000.0 * tmp6 + (4566965972049467.0 / 50000000000000000.0) * tmp9 +
                     1342707196092513.0 / 5000000000000000.0;
  const REAL tmp21 = ((tmp20) * (tmp20) * (tmp20) * (tmp20));
  const REAL tmp23 = tmp17 * ((tmp22) * (tmp22)) / ((tau_qnm) * (tau_qnm));
  const REAL tmp25 = 2 * log(2 * tmp20 / omega_qnm);
  const REAL tmp24 = (1.0 / sqrt(tmp21)) * tmp23;
  const REAL tmp26 = cosh(tmp25);
  const REAL tmp27 = sinh(tmp25);
  const REAL tmp28 = tanh(tmp25);
  const REAL tmp31 = (1.0 / 16.0) * ((omega_qnm) * (omega_qnm) * (omega_qnm) * (omega_qnm)) - tmp21;
  const REAL tmp29 = 1 - ((tmp28) * (tmp28));
  const REAL tmp30 = 1 - tmp28;
  const REAL tmp32 = tmp29 * tmp31 / tmp30;
  const REAL tmp33 = 0.25 * pow(tmp21, -1.5) * tmp23 * tmp32;
  const REAL h_t_attach = tmp17;
  const REAL hdot_t_attach = 0;
  const REAL hddot_t_attach = 0.1875 * pow(tmp21, -2.5) * tmp23 * ((tmp29) * (tmp29)) * ((tmp31) * (tmp31)) / ((tmp30) * (tmp30)) -
                              1.0 / 4.0 * tmp24 + (1.0 / 2.0) * tmp24 * ((tmp27) * (tmp27)) / ((tmp26) * (tmp26)) + tmp28 * tmp33 +
                              tmp27 * tmp33 / tmp26;
  const REAL w_t_attach = tmp22;
  const REAL wdot_t_attach = 0.5 * pow(tmp21, -0.75) * tmp32 / tau_qnm;

  amps[0] = h_t_attach;
  amps[1] = hdot_t_attach;
  amps[2] = hddot_t_attach;
  omegas[0] = w_t_attach;
  omegas[1] = wdot_t_attach;
} // END FUNCTION BOB_aligned_spin_NQC_rhs
