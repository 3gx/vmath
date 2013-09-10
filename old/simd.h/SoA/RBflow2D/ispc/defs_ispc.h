#pragma once

#include "dims.h"

#ifdef SINGLE
typedef float real;
#else
typedef double real;
#endif

inline real SQR(const real x) { return x*x;}


#define __out

struct Data
{
  int n2m, n3m;
  real ren, pec, dt;
  real al, ga, ro, beta;

  /* local array */
  real (*q2  )[M2], (*q3  )[M2];                /* q2[M2][M3], q3[M2][M3] */
  real (*pr  )[M2], (*dens)[M2], (*hro )[M2];   /* pr[M2][M3], dens[M2][M3], etc.. */
  real (*rhs1)[M2], (*rhs2)[M2], (*rhs3)[M2];
  real (*ru1 )[M2], (*ru2 )[M2], (*ru3 )[M2], (*ruro)[M2];
  real (*qcap)[M2], (*dph )[M2+16];

  real *denbs, *denbn;

  /* grid parameters */

  real dx2,  dx3;
  real dx2q, dx3q;
  real *rc, *rm, *g2rc, *g2rm;
  real *zz, *zm, *g3rc, *g3rm;
      
  /******* QUANTITIES FOR DERIVATIVES ********/
      
  real *udx2rm, *udx2c, *udx2m;
  real *usg2rc;
  real *udx3c, *udx3m;
      
  /******* Grid indices**************************************/

  int *jmv,*jmc,*jpv,*jpc;
  int *jmhv;
  int *kmc,*kpc,*kmv,*kpv,*kup,*kum;

  /******* Metric coefficients *******************************/

  real *ap2j, *ac2j, *am2j;
  real *ap2je, *ac2je, *am2je;
  real *ap3j, *ac3j, *am3j;
  real *apscj, *acscj, *amscj;

  real *ap3ck, *ac3ck, *am3ck;
  real *ap3sk, *ac3sk, *am3sk;
  real *ap3ssk, *ac3ssk, *am3ssk;
      
  /******* Variables for FFT and Poisson solver****************/

  real *ak2;
  real *amphj, *apphj, *acphj;
  real *amphk, *apphk, *acphk;

};

extern uniform Data data;
