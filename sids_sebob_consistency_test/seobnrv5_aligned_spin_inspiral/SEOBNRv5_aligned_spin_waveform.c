#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include <complex.h>
/**
 * Calculate the SEOBNRv5 22 mode.
 */
double complex SEOBNRv5_aligned_spin_waveform(REAL *restrict dynamics, commondata_struct *restrict commondata) {

  double complex gamma_22;
  const REAL m1 = commondata->m1;
  const REAL m2 = commondata->m2;
  const REAL chi1 = commondata->chi1;
  const REAL chi2 = commondata->chi2;
  const REAL phi = dynamics[PHI];
  const REAL Hreal = dynamics[H];
  const REAL Omega = dynamics[OMEGA];
  const REAL Omega_circ = dynamics[OMEGA_CIRC];
  // compute
  const REAL khat2 = 2 * Hreal * Omega;

  gamma_22 = SEOBNRv5_aligned_spin_gamma_wrapper(3., -2. * khat2);
  const double complex tmp0 = ((Omega) * (Omega));
  const double complex tmp1 = pow(Omega, 10.0 / 3.0);
  const double complex tmp2 = log(4 * cbrt(Omega));
  const double complex tmp3 = m1 + m2;
  const double complex tmp7 = pow(Omega, 8.0 / 3.0);
  const double complex tmp15 = pow(Omega, 5.0 / 3.0);
  const double complex tmp25 = Hreal * Omega;
  const double complex tmp10 = (1.0 / 2.0) * chi1 + (1.0 / 2.0) * chi2;
  const double complex tmp11 = (1.0 / 2.0) * chi1 - 1.0 / 2.0 * chi2;
  const double complex tmp12 = (1.0 / (tmp3));
  const double complex tmp16 = ((m1) * (m1)) * ((m2) * (m2)) / ((tmp3) * (tmp3) * (tmp3) * (tmp3));
  const double complex tmp23 = ((m1) * (m1) * (m1)) * ((m2) * (m2) * (m2)) / pow(tmp3, 6);
  const double complex tmp5 = (1.0 / ((tmp3) * (tmp3)));
  const double complex tmp13 = tmp12 * (m1 - m2);
  const double complex tmp17 = ((tmp10) * (tmp10));
  const double complex tmp18 = ((tmp11) * (tmp11));
  const double complex tmp6 = m1 * m2 * tmp5;
  const double complex tmp14 = tmp11 * tmp13;
  const double complex tmp19 = 2 * tmp6;
  const double complex tmp21 = tmp10 * (1 - tmp19) + tmp14;
  const double complex h22 =
      -1.0092530088080638 * I * M_PI * gamma_22 * tmp0 * tmp6 * (1 + (1.0 / 2.0) * ((tmp3) * (tmp3)) * (((Hreal) * (Hreal)) - 1) / (m1 * m2)) *
      (m1 * tmp12 + m2 * tmp12) *
      ((pow(Omega, 7.0 / 3.0) * (((tmp10) * (tmp10) * (tmp10)) * (tmp6 + 1.0 / 3.0) + tmp10 * tmp18 * (-4 * tmp16 - 3 * tmp6 + 1) +
                                 tmp10 * (-245717.0 / 63504.0 * tmp16 + (50803.0 / 63504.0) * tmp23 + (74749.0 / 5292.0) * tmp6 + 18733.0 / 15876.0) +
                                 ((tmp11) * (tmp11) * (tmp11)) * tmp13 * (1.0 / 3.0 - 4.0 / 3.0 * tmp6) + tmp14 * tmp17 * (tmp19 + 1) +
                                 tmp14 * ((97865.0 / 63504.0) * tmp16 + (50140.0 / 3969.0) * tmp6 + 18733.0 / 15876.0)) +
        pow(Omega, 4.0 / 3.0) * (tmp10 * tmp14 + (19583.0 / 42336.0) * tmp16 + (1.0 / 2.0) * tmp17 + tmp18 * (1.0 / 2.0 - tmp19) -
                                 33025.0 / 21168.0 * tmp6 - 20555.0 / 10584.0) +
        pow(Omega, 2.0 / 3.0) * ((55.0 / 84.0) * tmp6 - 43.0 / 42.0) +
        Omega * (-1.0 / 3.0 * chi1 - 1.0 / 3.0 * chi2 + (2.0 / 3.0) * m1 * m2 * tmp10 * tmp5 - 2.0 / 3.0 * tmp14) +
        tmp0 * (tmp10 * tmp14 * (89.0 / 126.0 - 781.0 / 252.0 * tmp6) - 6292061.0 / 3259872.0 * tmp16 +
                tmp17 * ((10.0 / 9.0) * tmp16 - 1817.0 / 504.0 * tmp6 + 89.0 / 252.0) +
                tmp18 * (-27.0 / 14.0 * tmp16 - 457.0 / 504.0 * tmp6 + 89.0 / 252.0) - 428.0 / 105.0 * tmp2 + (10620745.0 / 39118464.0) * tmp23 -
                48993925.0 / 9779616.0 * tmp6 + (41.0 / 192.0) * ((M_PI) * (M_PI)) * tmp6 + 39665487339946151125741.0 / 3820162500000000000000.0) +
        tmp1 * ((439877.0 / 55566.0) * tmp2 - 42670673890617980628777401.0 / 1668646980000000000000000.0) +
        tmp15 * (tmp10 * ((209.0 / 126.0) * tmp16 + (49.0 / 18.0) * tmp6 - 34.0 / 21.0) + tmp14 * (-19.0 / 42.0 * tmp6 - 34.0 / 21.0)) +
        tmp6 * (-411 * tmp1 + (106.0 / 5.0) * tmp7) +
        tmp7 * ((9202.0 / 2205.0) * tmp2 - 0.125 * ((tmp21) * (tmp21) * (tmp21) * (tmp21)) + 0.86701625094482238 * ((tmp21) * (tmp21)) -
                52376996519382937867657.0 / 6257426175000000000000000.0) +
        1) *
       (pow(Omega, 7.0 / 3.0) * (((tmp10) * (tmp10) * (tmp10)) * (tmp6 + 1.0 / 3.0) + tmp10 * tmp18 * (-4 * tmp16 - 3 * tmp6 + 1) +
                                 tmp10 * (-245717.0 / 63504.0 * tmp16 + (50803.0 / 63504.0) * tmp23 + (74749.0 / 5292.0) * tmp6 + 18733.0 / 15876.0) +
                                 ((tmp11) * (tmp11) * (tmp11)) * tmp13 * (1.0 / 3.0 - 4.0 / 3.0 * tmp6) + tmp14 * tmp17 * (tmp19 + 1) +
                                 tmp14 * ((97865.0 / 63504.0) * tmp16 + (50140.0 / 3969.0) * tmp6 + 18733.0 / 15876.0)) +
        pow(Omega, 4.0 / 3.0) * (tmp10 * tmp14 + (19583.0 / 42336.0) * tmp16 + (1.0 / 2.0) * tmp17 + tmp18 * (1.0 / 2.0 - tmp19) -
                                 33025.0 / 21168.0 * tmp6 - 20555.0 / 10584.0) +
        pow(Omega, 2.0 / 3.0) * ((55.0 / 84.0) * tmp6 - 43.0 / 42.0) +
        Omega * (-1.0 / 3.0 * chi1 - 1.0 / 3.0 * chi2 + (2.0 / 3.0) * m1 * m2 * tmp10 * tmp5 - 2.0 / 3.0 * tmp14) +
        tmp0 * (tmp10 * tmp14 * (89.0 / 126.0 - 781.0 / 252.0 * tmp6) - 6292061.0 / 3259872.0 * tmp16 +
                tmp17 * ((10.0 / 9.0) * tmp16 - 1817.0 / 504.0 * tmp6 + 89.0 / 252.0) +
                tmp18 * (-27.0 / 14.0 * tmp16 - 457.0 / 504.0 * tmp6 + 89.0 / 252.0) - 428.0 / 105.0 * tmp2 + (10620745.0 / 39118464.0) * tmp23 -
                48993925.0 / 9779616.0 * tmp6 + (41.0 / 192.0) * ((M_PI) * (M_PI)) * tmp6 + 39665487339946151125741.0 / 3820162500000000000000.0) +
        tmp1 * ((439877.0 / 55566.0) * tmp2 - 42670673890617980628777401.0 / 1668646980000000000000000.0) +
        tmp15 * (tmp10 * ((209.0 / 126.0) * tmp16 + (49.0 / 18.0) * tmp6 - 34.0 / 21.0) + tmp14 * (-19.0 / 42.0 * tmp6 - 34.0 / 21.0)) +
        tmp6 * (-411 * tmp1 + (106.0 / 5.0) * tmp7) +
        tmp7 * ((9202.0 / 2205.0) * tmp2 - 0.125 * ((tmp21) * (tmp21) * (tmp21) * (tmp21)) + 0.86701625094482238 * ((tmp21) * (tmp21)) -
                52376996519382937867657.0 / 6257426175000000000000000.0) +
        1)) *
      cexp(-2 * I * phi) *
      cexp(I * (((Hreal) * (Hreal) * (Hreal)) * ((Omega) * (Omega) * (Omega)) * (-2203.0 / 81.0 + (1712.0 / 315.0) * ((M_PI) * (M_PI))) +
                ((Hreal) * (Hreal)) * tmp0 * (tmp10 * ((8.0 / 3.0) * tmp6 - 4.0 / 3.0) - 4.0 / 3.0 * tmp14 + (428.0 / 105.0) * M_PI) -
                24 * tmp15 * tmp6 + (7.0 / 3.0) * tmp25)) *
      cexp(2 * M_PI * tmp25) * cexp(4 * I * tmp25 * log(8 * Omega * cexp(-1.0 / 2.0))) / pow(Omega_circ, 4.0 / 3.0);

  return h22;
} // END FUNCTION SEOBNRv5_aligned_spin_waveform
