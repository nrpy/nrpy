#include "BHaH_defines.h"
/*
 * // Compute at a single point (th,ph) the spin-weight -2 spherical harmonic Y_{s=-2, l,m}(th,ph)
 * // Manual "inline void" of this function results in compilation error with clang.
 *
 */
void spin_weight_minus2_sph_harmonics(const int l, const int m, const REAL th, const REAL ph, REAL *restrict reYlmswm2_l_m,
                                      REAL *restrict imYlmswm2_l_m) {
  switch (l) {
  case 0:
    switch (m) {
    case 0: {
      // l = 0, m = 0:
      *reYlmswm2_l_m = 0;
      *imYlmswm2_l_m = 0;
    }
      return;
    } // END switch(l == 0)
  case 1:
    switch (m) {
    case -1: {
      // l = 1, m = -1:
      *reYlmswm2_l_m = 0;
      *imYlmswm2_l_m = 0;
    }
      return;
    case 0: {
      // l = 1, m = 0:
      *reYlmswm2_l_m = 0;
      *imYlmswm2_l_m = 0;
    }
      return;
    case 1: {
      // l = 1, m = 1:
      *reYlmswm2_l_m = 0;
      *imYlmswm2_l_m = 0;
    }
      return;
    } // END switch(l == 1)
  case 2:
    switch (m) {
    case -2: {
      // l = 2, m = -2:
      const REAL tmp1 = (1.0 / 2.0) * sqrt(5) *
                        ((sin((1.0 / 2.0) * th)) * (sin((1.0 / 2.0) * th)) * (sin((1.0 / 2.0) * th)) * (sin((1.0 / 2.0) * th))) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp1 * cos(2 * ph);
      *imYlmswm2_l_m = -tmp1 * sin(2 * ph);
    }
      return;
    case -1: {
      // l = 2, m = -1:
      const REAL tmp1 = sqrt(5) * ((sin((1.0 / 2.0) * th)) * (sin((1.0 / 2.0) * th)) * (sin((1.0 / 2.0) * th))) * cos((1.0 / 2.0) * th) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp1 * cos(ph);
      *imYlmswm2_l_m = -tmp1 * sin(ph);
    }
      return;
    case 0: {
      // l = 2, m = 0:
      *reYlmswm2_l_m = (1.0 / 12.0) * sqrt(30) * (3.0 / 4.0 - 3.0 / 4.0 * cos(2 * th)) / sqrt(M_PI);
      *imYlmswm2_l_m = 0;
    }
      return;
    case 1: {
      // l = 2, m = 1:
      const REAL tmp1 = sqrt(5) * sin((1.0 / 2.0) * th) * ((cos((1.0 / 2.0) * th)) * (cos((1.0 / 2.0) * th)) * (cos((1.0 / 2.0) * th))) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp1 * cos(ph);
      *imYlmswm2_l_m = tmp1 * sin(ph);
    }
      return;
    case 2: {
      // l = 2, m = 2:
      const REAL tmp1 = (1.0 / 2.0) * sqrt(5) *
                        ((cos((1.0 / 2.0) * th)) * (cos((1.0 / 2.0) * th)) * (cos((1.0 / 2.0) * th)) * (cos((1.0 / 2.0) * th))) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp1 * cos(2 * ph);
      *imYlmswm2_l_m = tmp1 * sin(2 * ph);
    }
      return;
    } // END switch(l == 2)
  case 3:
    switch (m) {
    case -3: {
      // l = 3, m = -3:
      const REAL tmp2 =
          (1.0 / 2.0) * sqrt(42) *
          ((sin((1.0 / 2.0) * th)) * (sin((1.0 / 2.0) * th)) * (sin((1.0 / 2.0) * th)) * (sin((1.0 / 2.0) * th)) * (sin((1.0 / 2.0) * th))) *
          cos((1.0 / 2.0) * th) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp2 * cos(3 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(3 * ph);
    }
      return;
    case -2: {
      // l = 3, m = -2:
      const REAL tmp1 = sin((1.0 / 2.0) * th);
      const REAL tmp2 = (1.0 / 2.0) * sqrt(7) * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) * (5 - 6 * ((tmp1) * (tmp1))) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp2 * cos(2 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(2 * ph);
    }
      return;
    case -1: {
      // l = 3, m = -1:
      const REAL tmp0 = (1.0 / 64.0) * sqrt(70) * (sin(th) + 4 * sin(2 * th) - 3 * sin(3 * th)) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp0 * cos(ph);
      *imYlmswm2_l_m = -tmp0 * sin(ph);
    }
      return;
    case 0: {
      // l = 3, m = 0:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      *reYlmswm2_l_m =
          (1.0 / 2.0) * sqrt(210) * (1 - ((tmp1) * (tmp1))) * pow(sin((1.0 / 2.0) * th), 6) / (sqrt(M_PI) * ((tmp1) * (tmp1) * (tmp1) * (tmp1)));
      *imYlmswm2_l_m = 0;
    }
      return;
    case 1: {
      // l = 3, m = 1:
      const REAL tmp0 = (1.0 / 64.0) * sqrt(70) * (-sin(th) + 4 * sin(2 * th) + 3 * sin(3 * th)) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp0 * cos(ph);
      *imYlmswm2_l_m = tmp0 * sin(ph);
    }
      return;
    case 2: {
      // l = 3, m = 2:
      const REAL tmp1 = cos(th);
      const REAL tmp2 = (1.0 / 8.0) * sqrt(7) * ((tmp1 + 1) * (tmp1 + 1)) * (3 * tmp1 - 2) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp2 * cos(2 * ph);
      *imYlmswm2_l_m = tmp2 * sin(2 * ph);
    }
      return;
    case 3: {
      // l = 3, m = 3:
      const REAL tmp1 = (1.0 / 16.0) * sqrt(42) * ((cos(th) + 1) * (cos(th) + 1)) * sin(th) / sqrt(M_PI);
      *reYlmswm2_l_m = -tmp1 * cos(3 * ph);
      *imYlmswm2_l_m = -tmp1 * sin(3 * ph);
    }
      return;
    } // END switch(l == 3)
  case 4:
    switch (m) {
    case -4: {
      // l = 4, m = -4:
      const REAL tmp2 = 3 * sqrt(7) * pow(sin((1.0 / 2.0) * th), 8) / (sqrt(M_PI) * ((tan((1.0 / 2.0) * th)) * (tan((1.0 / 2.0) * th))));
      *reYlmswm2_l_m = tmp2 * cos(4 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(4 * ph);
    }
      return;
    case -3: {
      // l = 4, m = -3:
      const REAL tmp1 = (3.0 / 16.0) * sqrt(14) * ((1 - cos(th)) * (1 - cos(th))) * (sin(th) + sin(2 * th)) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp1 * cos(3 * ph);
      *imYlmswm2_l_m = -tmp1 * sin(3 * ph);
    }
      return;
    case -2: {
      // l = 4, m = -2:
      const REAL tmp1 = cos(th);
      const REAL tmp2 = ((1 - tmp1) * (1 - tmp1));
      const REAL tmp3 = (3.0 / 8.0) * tmp2 * (21 * tmp1 + 7 * tmp2 - 6) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp3 * cos(2 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(2 * ph);
    }
      return;
    case -1: {
      // l = 4, m = -1:
      const REAL tmp1 = sin((1.0 / 2.0) * th);
      const REAL tmp2 = (3.0 / 2.0) * M_SQRT2 * ((tmp1) * (tmp1) * (tmp1)) *
                        (28 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 35 * ((tmp1) * (tmp1)) + 10) * cos((1.0 / 2.0) * th) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp2 * cos(ph);
      *imYlmswm2_l_m = -tmp2 * sin(ph);
    }
      return;
    case 0: {
      // l = 4, m = 0:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      *reYlmswm2_l_m = (3.0 / 10.0) * sqrt(10) * (15 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 40 * ((tmp1) * (tmp1)) + 15) *
                       pow(sin((1.0 / 2.0) * th), 8) / (sqrt(M_PI) * pow(tmp1, 6));
      *imYlmswm2_l_m = 0;
    }
      return;
    case 1: {
      // l = 4, m = 1:
      const REAL tmp0 = tan((1.0 / 2.0) * th);
      const REAL tmp1 = (3.0 / 32.0) * M_SQRT2 * ((cos(th) + 1) * (cos(th) + 1) * (cos(th) + 1)) *
                        (10 * ((tmp0) * (tmp0) * (tmp0) * (tmp0)) - 15 * ((tmp0) * (tmp0)) + 3) * sin(th) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp1 * cos(ph);
      *imYlmswm2_l_m = tmp1 * sin(ph);
    }
      return;
    case 2: {
      // l = 4, m = 2:
      const REAL tmp1 = cos(th);
      const REAL tmp2 = (3.0 / 8.0) * ((tmp1 + 1) * (tmp1 + 1)) * (7 * tmp1 + 7 * ((1 - tmp1) * (1 - tmp1)) - 6) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp2 * cos(2 * ph);
      *imYlmswm2_l_m = tmp2 * sin(2 * ph);
    }
      return;
    case 3: {
      // l = 4, m = 3:
      const REAL tmp2 = cos((1.0 / 2.0) * th);
      const REAL tmp3 =
          (3.0 / 2.0) * sqrt(14) * ((tmp2) * (tmp2) * (tmp2) * (tmp2) * (tmp2)) * (3 - 4 * ((tmp2) * (tmp2))) * sin((1.0 / 2.0) * th) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp3 * cos(3 * ph);
      *imYlmswm2_l_m = tmp3 * sin(3 * ph);
    }
      return;
    case 4: {
      // l = 4, m = 4:
      const REAL tmp2 = 3 * sqrt(7) * pow(sin((1.0 / 2.0) * th), 8) / (sqrt(M_PI) * pow(tan((1.0 / 2.0) * th), 6));
      *reYlmswm2_l_m = tmp2 * cos(4 * ph);
      *imYlmswm2_l_m = tmp2 * sin(4 * ph);
    }
      return;
    } // END switch(l == 4)
  case 5:
    switch (m) {
    case -5: {
      // l = 5, m = -5:
      const REAL tmp2 =
          sqrt(330) * pow(sin((1.0 / 2.0) * th), 10) / (sqrt(M_PI) * ((tan((1.0 / 2.0) * th)) * (tan((1.0 / 2.0) * th)) * (tan((1.0 / 2.0) * th))));
      *reYlmswm2_l_m = tmp2 * cos(5 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(5 * ph);
    }
      return;
    case -4: {
      // l = 5, m = -4:
      const REAL tmp1 = cos(th);
      const REAL tmp2 = (1.0 / 32.0) * sqrt(33) * ((1 - tmp1) * (1 - tmp1) * (1 - tmp1)) * (14 * tmp1 + 5 * cos(2 * th) + 9) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp2 * cos(4 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(4 * ph);
    }
      return;
    case -3: {
      // l = 5, m = -3:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = sqrt(66) * (((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 7 * ((tmp2) * (tmp2)) + 7) * pow(sin((1.0 / 2.0) * th), 10) /
                        (sqrt(M_PI) * ((tmp2) * (tmp2) * (tmp2) * (tmp2) * (tmp2)));
      *reYlmswm2_l_m = tmp3 * cos(3 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(3 * ph);
    }
      return;
    case -2: {
      // l = 5, m = -2:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = pow(tmp2, 6);
      const REAL tmp4 = (1.0 / 2.0) * sqrt(11) * (21 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 63 * ((tmp2) * (tmp2)) - tmp3 + 35) *
                        pow(sin((1.0 / 2.0) * th), 10) / (sqrt(M_PI) * tmp3);
      *reYlmswm2_l_m = tmp4 * cos(2 * ph);
      *imYlmswm2_l_m = -tmp4 * sin(2 * ph);
    }
      return;
    case -1: {
      // l = 5, m = -1:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      const REAL tmp2 = sqrt(77) * (-pow(tmp1, 6) + 9 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 15 * ((tmp1) * (tmp1)) + 5) *
                        pow(sin((1.0 / 2.0) * th), 10) / (sqrt(M_PI) * pow(tmp1, 7));
      *reYlmswm2_l_m = tmp2 * cos(ph);
      *imYlmswm2_l_m = -tmp2 * sin(ph);
    }
      return;
    case 0: {
      // l = 5, m = 0:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      *reYlmswm2_l_m = (1.0 / 2.0) * sqrt(2310) * (-pow(tmp1, 6) + 5 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 5 * ((tmp1) * (tmp1)) + 1) *
                       pow(sin((1.0 / 2.0) * th), 10) / (sqrt(M_PI) * pow(tmp1, 8));
      *imYlmswm2_l_m = 0;
    }
      return;
    case 1: {
      // l = 5, m = 1:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      const REAL tmp2 = (1.0 / 7.0) * sqrt(77) * (-35 * pow(tmp1, 6) + 105 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 63 * ((tmp1) * (tmp1)) + 7) *
                        pow(sin((1.0 / 2.0) * th), 10) / (sqrt(M_PI) * pow(tmp1, 9));
      *reYlmswm2_l_m = tmp2 * cos(ph);
      *imYlmswm2_l_m = tmp2 * sin(ph);
    }
      return;
    case 2: {
      // l = 5, m = 2:
      const REAL tmp2 = cos((1.0 / 2.0) * th);
      const REAL tmp3 = sin((1.0 / 2.0) * th);
      const REAL tmp4 = (1.0 / 2.0) * sqrt(11) * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) *
                        (pow(tmp2, 6) - 119 * pow(tmp3, 6) + 105 * ((tmp3) * (tmp3) * (tmp3) * (tmp3)) - 21 * ((tmp3) * (tmp3))) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp4 * cos(2 * ph);
      *imYlmswm2_l_m = tmp4 * sin(2 * ph);
    }
      return;
    case 3: {
      // l = 5, m = 3:
      const REAL tmp2 = sin((1.0 / 2.0) * th);
      const REAL tmp3 =
          sqrt(66) * tmp2 * (-15 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) + 9 * ((tmp2) * (tmp2)) - 1) *
          ((cos((1.0 / 2.0) * th)) * (cos((1.0 / 2.0) * th)) * (cos((1.0 / 2.0) * th)) * (cos((1.0 / 2.0) * th)) * (cos((1.0 / 2.0) * th))) /
          sqrt(M_PI);
      *reYlmswm2_l_m = tmp3 * cos(3 * ph);
      *imYlmswm2_l_m = tmp3 * sin(3 * ph);
    }
      return;
    case 4: {
      // l = 5, m = 4:
      const REAL tmp1 = cos(th);
      const REAL tmp2 = (1.0 / 32.0) * sqrt(33) * ((tmp1 + 1) * (tmp1 + 1) * (tmp1 + 1)) * (14 * tmp1 - 5 * cos(2 * th) - 9) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp2 * cos(4 * ph);
      *imYlmswm2_l_m = tmp2 * sin(4 * ph);
    }
      return;
    case 5: {
      // l = 5, m = 5:
      const REAL tmp2 = sqrt(330) * pow(sin((1.0 / 2.0) * th), 10) / (sqrt(M_PI) * pow(tan((1.0 / 2.0) * th), 7));
      *reYlmswm2_l_m = -tmp2 * cos(5 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(5 * ph);
    }
      return;
    } // END switch(l == 5)
  case 6:
    switch (m) {
    case -6: {
      // l = 6, m = -6:
      const REAL tmp2 = (3.0 / 2.0) * sqrt(715) * pow(sin((1.0 / 2.0) * th), 12) /
                        (sqrt(M_PI) * ((tan((1.0 / 2.0) * th)) * (tan((1.0 / 2.0) * th)) * (tan((1.0 / 2.0) * th)) * (tan((1.0 / 2.0) * th))));
      *reYlmswm2_l_m = tmp2 * cos(6 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(6 * ph);
    }
      return;
    case -5: {
      // l = 6, m = -5:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 =
          sqrt(2145) * (2 - ((tmp2) * (tmp2))) * pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * ((tmp2) * (tmp2) * (tmp2) * (tmp2) * (tmp2)));
      *reYlmswm2_l_m = tmp3 * cos(5 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(5 * ph);
    }
      return;
    case -4: {
      // l = 6, m = -4:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 4.0) * sqrt(390) * (6 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 32 * ((tmp2) * (tmp2)) + 28) *
                        pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * pow(tmp2, 6));
      *reYlmswm2_l_m = tmp3 * cos(4 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(4 * ph);
    }
      return;
    case -3: {
      // l = 6, m = -3:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = 3 * sqrt(13) * (-pow(tmp2, 6) + 12 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 28 * ((tmp2) * (tmp2)) + 14) *
                        pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * pow(tmp2, 7));
      *reYlmswm2_l_m = tmp3 * cos(3 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(3 * ph);
    }
      return;
    case -2: {
      // l = 6, m = -2:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = pow(tmp2, 8);
      const REAL tmp4 = (1.0 / 2.0) * sqrt(13) *
                        (-32 * pow(tmp2, 6) + 168 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 224 * ((tmp2) * (tmp2)) + tmp3 + 70) *
                        pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * tmp3);
      *reYlmswm2_l_m = tmp4 * cos(2 * ph);
      *imYlmswm2_l_m = -tmp4 * sin(2 * ph);
    }
      return;
    case -1: {
      // l = 6, m = -1:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      const REAL tmp2 = sqrt(130) * (pow(tmp1, 8) - 14 * pow(tmp1, 6) + 42 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 35 * ((tmp1) * (tmp1)) + 7) *
                        pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * pow(tmp1, 9));
      *reYlmswm2_l_m = tmp2 * cos(ph);
      *imYlmswm2_l_m = -tmp2 * sin(ph);
    }
      return;
    case 0: {
      // l = 6, m = 0:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      *reYlmswm2_l_m = sqrt(1365) * (pow(tmp1, 8) - 8 * pow(tmp1, 6) + 15 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 8 * ((tmp1) * (tmp1)) + 1) *
                       pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * pow(tmp1, 10));
      *imYlmswm2_l_m = 0;
    }
      return;
    case 1: {
      // l = 6, m = 1:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      const REAL tmp2 = (1.0 / 8.0) * sqrt(130) *
                        (56 * pow(tmp1, 8) - 280 * pow(tmp1, 6) + 336 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 112 * ((tmp1) * (tmp1)) + 8) *
                        pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * pow(tmp1, 11));
      *reYlmswm2_l_m = tmp2 * cos(ph);
      *imYlmswm2_l_m = tmp2 * sin(ph);
    }
      return;
    case 2: {
      // l = 6, m = 2:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 2.0) * sqrt(13) *
                        (70 * pow(tmp2, 8) - 224 * pow(tmp2, 6) + 168 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 32 * ((tmp2) * (tmp2)) + 1) *
                        pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * pow(tmp2, 12));
      *reYlmswm2_l_m = tmp3 * cos(2 * ph);
      *imYlmswm2_l_m = tmp3 * sin(2 * ph);
    }
      return;
    case 3: {
      // l = 6, m = 3:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (3.0 / 4.0) * sqrt(13) * (56 * pow(tmp2, 6) - 112 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) + 48 * ((tmp2) * (tmp2)) - 4) *
                        pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * pow(tmp2, 11));
      *reYlmswm2_l_m = tmp3 * cos(3 * ph);
      *imYlmswm2_l_m = tmp3 * sin(3 * ph);
    }
      return;
    case 4: {
      // l = 6, m = 4:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 4.0) * sqrt(390) * (28 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 32 * ((tmp2) * (tmp2)) + 6) *
                        pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * pow(tmp2, 10));
      *reYlmswm2_l_m = tmp3 * cos(4 * ph);
      *imYlmswm2_l_m = tmp3 * sin(4 * ph);
    }
      return;
    case 5: {
      // l = 6, m = 5:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 4.0) * sqrt(2145) * (8 * ((tmp2) * (tmp2)) - 4) * pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * pow(tmp2, 9));
      *reYlmswm2_l_m = tmp3 * cos(5 * ph);
      *imYlmswm2_l_m = tmp3 * sin(5 * ph);
    }
      return;
    case 6: {
      // l = 6, m = 6:
      const REAL tmp2 = (3.0 / 2.0) * sqrt(715) * pow(sin((1.0 / 2.0) * th), 12) / (sqrt(M_PI) * pow(tan((1.0 / 2.0) * th), 8));
      *reYlmswm2_l_m = tmp2 * cos(6 * ph);
      *imYlmswm2_l_m = tmp2 * sin(6 * ph);
    }
      return;
    } // END switch(l == 6)
  case 7:
    switch (m) {
    case -7: {
      // l = 7, m = -7:
      const REAL tmp2 = (1.0 / 2.0) * sqrt(30030) * pow(sin((1.0 / 2.0) * th), 14) /
                        (sqrt(M_PI) * ((tan((1.0 / 2.0) * th)) * (tan((1.0 / 2.0) * th)) * (tan((1.0 / 2.0) * th)) * (tan((1.0 / 2.0) * th)) *
                                       (tan((1.0 / 2.0) * th))));
      *reYlmswm2_l_m = tmp2 * cos(7 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(7 * ph);
    }
      return;
    case -6: {
      // l = 7, m = -6:
      const REAL tmp1 = cos(th);
      const REAL tmp2 =
          (1.0 / 128.0) * sqrt(2145) * ((1 - tmp1) * (1 - tmp1) * (1 - tmp1) * (1 - tmp1)) * ((tmp1 + 1) * (tmp1 + 1)) * (7 * tmp1 + 2) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp2 * cos(6 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(6 * ph);
    }
      return;
    case -5: {
      // l = 7, m = -5:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 2.0) * sqrt(330) * (10 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 45 * ((tmp2) * (tmp2)) + 36) *
                        pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tmp2, 7));
      *reYlmswm2_l_m = tmp3 * cos(5 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(5 * ph);
    }
      return;
    case -4: {
      // l = 7, m = -4:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 4.0) * sqrt(330) * (-10 * pow(tmp2, 6) + 90 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 180 * ((tmp2) * (tmp2)) + 84) *
                        pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tmp2, 8));
      *reYlmswm2_l_m = tmp3 * cos(4 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(4 * ph);
    }
      return;
    case -3: {
      // l = 7, m = -3:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 2.0) * sqrt(30) *
                        (5 * pow(tmp2, 8) - 90 * pow(tmp2, 6) + 360 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 420 * ((tmp2) * (tmp2)) + 126) *
                        pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tmp2, 9));
      *reYlmswm2_l_m = tmp3 * cos(3 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(3 * ph);
    }
      return;
    case -2: {
      // l = 7, m = -2:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = pow(tmp2, 10);
      const REAL tmp4 = (1.0 / 2.0) * sqrt(15) *
                        (45 * pow(tmp2, 8) - 360 * pow(tmp2, 6) + 840 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 630 * ((tmp2) * (tmp2)) - tmp3 + 126) *
                        pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * tmp3);
      *reYlmswm2_l_m = tmp4 * cos(2 * ph);
      *imYlmswm2_l_m = -tmp4 * sin(2 * ph);
    }
      return;
    case -1: {
      // l = 7, m = -1:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      const REAL tmp2 =
          (1.0 / 2.0) * sqrt(10) *
          (-9 * pow(tmp1, 10) + 180 * pow(tmp1, 8) - 840 * pow(tmp1, 6) + 1260 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 630 * ((tmp1) * (tmp1)) + 84) *
          pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tmp1, 11));
      *reYlmswm2_l_m = tmp2 * cos(ph);
      *imYlmswm2_l_m = -tmp2 * sin(ph);
    }
      return;
    case 0: {
      // l = 7, m = 0:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      *reYlmswm2_l_m = (1.0 / 4.0) * sqrt(35) *
                       (-36 * pow(tmp1, 10) + 420 * pow(tmp1, 8) - 1260 * pow(tmp1, 6) + 1260 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) -
                        420 * ((tmp1) * (tmp1)) + 36) *
                       pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tmp1, 12));
      *imYlmswm2_l_m = 0;
    }
      return;
    case 1: {
      // l = 7, m = 1:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      const REAL tmp2 =
          (1.0 / 2.0) * sqrt(10) *
          (-84 * pow(tmp1, 10) + 630 * pow(tmp1, 8) - 1260 * pow(tmp1, 6) + 840 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 180 * ((tmp1) * (tmp1)) + 9) *
          pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tmp1, 13));
      *reYlmswm2_l_m = tmp2 * cos(ph);
      *imYlmswm2_l_m = tmp2 * sin(ph);
    }
      return;
    case 2: {
      // l = 7, m = 2:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 =
          (1.0 / 2.0) * sqrt(15) *
          (-126 * pow(tmp2, 10) + 630 * pow(tmp2, 8) - 840 * pow(tmp2, 6) + 360 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 45 * ((tmp2) * (tmp2)) + 1) *
          pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tmp2, 14));
      *reYlmswm2_l_m = tmp3 * cos(2 * ph);
      *imYlmswm2_l_m = tmp3 * sin(2 * ph);
    }
      return;
    case 3: {
      // l = 7, m = 3:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 2.0) * sqrt(30) *
                        (-126 * pow(tmp2, 8) + 420 * pow(tmp2, 6) - 360 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) + 90 * ((tmp2) * (tmp2)) - 5) *
                        pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tmp2, 13));
      *reYlmswm2_l_m = tmp3 * cos(3 * ph);
      *imYlmswm2_l_m = tmp3 * sin(3 * ph);
    }
      return;
    case 4: {
      // l = 7, m = 4:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 4.0) * sqrt(330) * (-84 * pow(tmp2, 6) + 180 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 90 * ((tmp2) * (tmp2)) + 10) *
                        pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tmp2, 12));
      *reYlmswm2_l_m = tmp3 * cos(4 * ph);
      *imYlmswm2_l_m = tmp3 * sin(4 * ph);
    }
      return;
    case 5: {
      // l = 7, m = 5:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 2.0) * sqrt(330) * (-36 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) + 45 * ((tmp2) * (tmp2)) - 10) *
                        pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tmp2, 11));
      *reYlmswm2_l_m = tmp3 * cos(5 * ph);
      *imYlmswm2_l_m = tmp3 * sin(5 * ph);
    }
      return;
    case 6: {
      // l = 7, m = 6:
      const REAL tmp1 = cos(th);
      const REAL tmp2 =
          (1.0 / 128.0) * sqrt(2145) * ((1 - tmp1) * (1 - tmp1)) * ((tmp1 + 1) * (tmp1 + 1) * (tmp1 + 1) * (tmp1 + 1)) * (7 * tmp1 - 2) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp2 * cos(6 * ph);
      *imYlmswm2_l_m = tmp2 * sin(6 * ph);
    }
      return;
    case 7: {
      // l = 7, m = 7:
      const REAL tmp2 = (1.0 / 2.0) * sqrt(30030) * pow(sin((1.0 / 2.0) * th), 14) / (sqrt(M_PI) * pow(tan((1.0 / 2.0) * th), 9));
      *reYlmswm2_l_m = -tmp2 * cos(7 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(7 * ph);
    }
      return;
    } // END switch(l == 7)
  case 8:
    switch (m) {
    case -8: {
      // l = 8, m = -8:
      const REAL tmp2 = sqrt(34034) * pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tan((1.0 / 2.0) * th), 6));
      *reYlmswm2_l_m = tmp2 * cos(8 * ph);
      *imYlmswm2_l_m = -tmp2 * sin(8 * ph);
    }
      return;
    case -7: {
      // l = 8, m = -7:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 4.0) * sqrt(34034) * (10 - 6 * ((tmp2) * (tmp2))) * pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp2, 7));
      *reYlmswm2_l_m = tmp3 * cos(7 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(7 * ph);
    }
      return;
    case -6: {
      // l = 8, m = -6:
      const REAL tmp1 = cos(th);
      const REAL tmp3 = (1.0 / 128.0) * sqrt(255255) * ((1 - tmp1) * (1 - tmp1) * (1 - tmp1) * (1 - tmp1)) * ((tmp1 + 1) * (tmp1 + 1)) *
                        (tmp1 + cos(2 * th) + 1) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp3 * cos(6 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(6 * ph);
    }
      return;
    case -5: {
      // l = 8, m = -5:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 20.0) * sqrt(24310) *
                        (-20 * pow(tmp2, 6) + 150 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 270 * ((tmp2) * (tmp2)) + 120) *
                        pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp2, 9));
      *reYlmswm2_l_m = tmp3 * cos(5 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(5 * ph);
    }
      return;
    case -4: {
      // l = 8, m = -4:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 10.0) * sqrt(1870) *
                        (15 * pow(tmp2, 8) - 200 * pow(tmp2, 6) + 675 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 720 * ((tmp2) * (tmp2)) + 210) *
                        pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp2, 10));
      *reYlmswm2_l_m = tmp3 * cos(4 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(4 * ph);
    }
      return;
    case -3: {
      // l = 8, m = -3:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 =
          (1.0 / 2.0) * sqrt(1122) *
          (-pow(tmp2, 10) + 25 * pow(tmp2, 8) - 150 * pow(tmp2, 6) + 300 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 210 * ((tmp2) * (tmp2)) + 42) *
          pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp2, 11));
      *reYlmswm2_l_m = tmp3 * cos(3 * ph);
      *imYlmswm2_l_m = -tmp3 * sin(3 * ph);
    }
      return;
    case -2: {
      // l = 8, m = -2:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = pow(tmp2, 12);
      const REAL tmp4 = (1.0 / 2.0) * sqrt(17) *
                        (-60 * pow(tmp2, 10) + 675 * pow(tmp2, 8) - 2400 * pow(tmp2, 6) + 3150 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) -
                         1512 * ((tmp2) * (tmp2)) + tmp3 + 210) *
                        pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * tmp3);
      *reYlmswm2_l_m = tmp4 * cos(2 * ph);
      *imYlmswm2_l_m = -tmp4 * sin(2 * ph);
    }
      return;
    case -1: {
      // l = 8, m = -1:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      const REAL tmp2 = (1.0 / 2.0) * sqrt(1190) *
                        (pow(tmp1, 12) - 27 * pow(tmp1, 10) + 180 * pow(tmp1, 8) - 420 * pow(tmp1, 6) + 378 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) -
                         126 * ((tmp1) * (tmp1)) + 12) *
                        pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp1, 13));
      *reYlmswm2_l_m = tmp2 * cos(ph);
      *imYlmswm2_l_m = -tmp2 * sin(ph);
    }
      return;
    case 0: {
      // l = 8, m = 0:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      *reYlmswm2_l_m = 3 * sqrt(595) *
                       (pow(tmp1, 12) - 16 * pow(tmp1, 10) + 70 * pow(tmp1, 8) - 112 * pow(tmp1, 6) + 70 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) -
                        16 * ((tmp1) * (tmp1)) + 1) *
                       pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp1, 14));
      *imYlmswm2_l_m = 0;
    }
      return;
    case 1: {
      // l = 8, m = 1:
      const REAL tmp1 = tan((1.0 / 2.0) * th);
      const REAL tmp2 = (1.0 / 20.0) * sqrt(1190) *
                        (120 * pow(tmp1, 12) - 1260 * pow(tmp1, 10) + 3780 * pow(tmp1, 8) - 4200 * pow(tmp1, 6) +
                         1800 * ((tmp1) * (tmp1) * (tmp1) * (tmp1)) - 270 * ((tmp1) * (tmp1)) + 10) *
                        pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp1, 15));
      *reYlmswm2_l_m = tmp2 * cos(ph);
      *imYlmswm2_l_m = tmp2 * sin(ph);
    }
      return;
    case 2: {
      // l = 8, m = 2:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 2.0) * sqrt(17) *
                        (210 * pow(tmp2, 12) - 1512 * pow(tmp2, 10) + 3150 * pow(tmp2, 8) - 2400 * pow(tmp2, 6) +
                         675 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 60 * ((tmp2) * (tmp2)) + 1) *
                        pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp2, 16));
      *reYlmswm2_l_m = tmp3 * cos(2 * ph);
      *imYlmswm2_l_m = tmp3 * sin(2 * ph);
    }
      return;
    case 3: {
      // l = 8, m = 3:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 12.0) * sqrt(1122) *
                        (252 * pow(tmp2, 10) - 1260 * pow(tmp2, 8) + 1800 * pow(tmp2, 6) - 900 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) +
                         150 * ((tmp2) * (tmp2)) - 6) *
                        pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp2, 15));
      *reYlmswm2_l_m = tmp3 * cos(3 * ph);
      *imYlmswm2_l_m = tmp3 * sin(3 * ph);
    }
      return;
    case 4: {
      // l = 8, m = 4:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 10.0) * sqrt(1870) *
                        (210 * pow(tmp2, 8) - 720 * pow(tmp2, 6) + 675 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) - 200 * ((tmp2) * (tmp2)) + 15) *
                        pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp2, 14));
      *reYlmswm2_l_m = tmp3 * cos(4 * ph);
      *imYlmswm2_l_m = tmp3 * sin(4 * ph);
    }
      return;
    case 5: {
      // l = 8, m = 5:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 20.0) * sqrt(24310) * (120 * pow(tmp2, 6) - 270 * ((tmp2) * (tmp2) * (tmp2) * (tmp2)) + 150 * ((tmp2) * (tmp2)) - 20) *
                        pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp2, 13));
      *reYlmswm2_l_m = tmp3 * cos(5 * ph);
      *imYlmswm2_l_m = tmp3 * sin(5 * ph);
    }
      return;
    case 6: {
      // l = 8, m = 6:
      const REAL tmp1 = cos(th);
      const REAL tmp2 = 1 - tmp1;
      const REAL tmp3 =
          (1.0 / 128.0) * sqrt(255255) * ((tmp2) * (tmp2)) * ((tmp1 + 1) * (tmp1 + 1) * (tmp1 + 1) * (tmp1 + 1)) * (tmp2 + cos(2 * th)) / sqrt(M_PI);
      *reYlmswm2_l_m = tmp3 * cos(6 * ph);
      *imYlmswm2_l_m = tmp3 * sin(6 * ph);
    }
      return;
    case 7: {
      // l = 8, m = 7:
      const REAL tmp2 = tan((1.0 / 2.0) * th);
      const REAL tmp3 = (1.0 / 4.0) * sqrt(34034) * (10 * ((tmp2) * (tmp2)) - 6) * pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tmp2, 11));
      *reYlmswm2_l_m = tmp3 * cos(7 * ph);
      *imYlmswm2_l_m = tmp3 * sin(7 * ph);
    }
      return;
    case 8: {
      // l = 8, m = 8:
      const REAL tmp2 = sqrt(34034) * pow(sin((1.0 / 2.0) * th), 16) / (sqrt(M_PI) * pow(tan((1.0 / 2.0) * th), 10));
      *reYlmswm2_l_m = tmp2 * cos(8 * ph);
      *imYlmswm2_l_m = tmp2 * sin(8 * ph);
    }
      return;
    } // END switch(l == 8)
  }   // END switch blocks

  fprintf(stderr, "ERROR: SpinWeight_minus2_SphHarmonics handles only l=[0,swm2sh_maximum_l_mode_generated=8] and only m=[-l,+l] is defined.\n");
  fprintf(stderr, "       You chose l=%d and m=%d, which is out of these bounds.\n", l, m);
  exit(1);
}
