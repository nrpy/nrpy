/*
TwoPunctures creates initial for two puncture black holes using a
single domain spectral method.  This method is described in
Marcus Ansorg, Bernd Brügmann, Wolfgang Tichy, "A single-domain
spectral method for black hole puncture data", PRD 70, 064011 (2004),
arXiv:gr-qc/0404056.

Code originally from Einstein Toolkit's EinsteinInitialData
arrangement:
https://bitbucket.org/einsteintoolkit/einsteininitialdata

License: Lesser GNU Public License, version 2.0+
*/

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
