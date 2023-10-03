/* TwoPunctures:  File  "TwoPunctures.h"*/
#ifndef __TWOPUNCTURES_H__
#define __TWOPUNCTURES_H__

#include "stdbool.h" // for the bools below
#ifndef REAL
#define REAL double // for the REALs below
#endif

#define StencilSize 19
#define N_PlaneRelax 1
#define NRELAX 200
#define Step_Relax 1

typedef struct __derivs__ {
  REAL *d0, *d1, *d2, *d3, *d11, *d12, *d13, *d22, *d23, *d33;
} derivs;

#endif // #ifndef __TWOPUNCTURES_H__
