#pragma once

/************************************************************************
 *   this subroutine performs the inversion of the q3 momentum equation
 *   by a factored implicit scheme, only the derivatives 11,22,33 of q3
 *   are treated implicitly
 *       direction x2
 ************************************************************************/
void solq3j()
{
#if 0  /* scalar */
  Tridiag<real> mat[NX];
  real u[NX], d[NX];
  
  const real betadx = beta*al;
  for (int jc = 0; jc < n2m; jc++)
  {
    mat[jc].a =     - am3j(jc)*betadx;
    mat[jc].b = 1.0 - ac3j(jc)*betadx;
    mat[jc].c =     - ap3j(jc)*betadx;
  }

  Tridiag<real>::Reduce(n2m, mat, u);

  const double t0 = rtc();
  for (int kc = 0; kc < n3m; kc++)
  {
    for (int jc = 0; jc < n2m; jc++)
      d[jc] = rhs(jc,kc);

    Tridiag<real>::Solve(n2m, mat, u, d);

    for (int jc = 0; jc < n2m; jc++)
      rhs(jc,kc) = d[jc];
  }
  mytimer[60] += rtc() - t0;

#else /*vector */

  static Tridiag<vreal> mat[NX];

  const vreal betadx = -beta*al;
  assert(n2m % vreal::WIDTH == 0);

#pragma omp parallel for
  for (int jc = 0; jc < n2m; jc += vreal::WIDTH)
  {
    const vreal a =              VLDA(am3j(jc))*betadx; 
    const vreal b = vreal(1.0) + VLDA(ac3j(jc))*betadx; 
    const vreal c =              VLDA(ap3j(jc))*betadx; 

#define UNROLL(ch) {\
  mat [jc+ch].a = SIMD::broadcast<ch>(a); \
  mat [jc+ch].b = SIMD::broadcast<ch>(b); \
  mat [jc+ch].c = SIMD::broadcast<ch>(c); }

    assert(vreal::WIDTH >= 2);
    assert(vreal::WIDTH <= 16);
    if (vreal::WIDTH > 0) { UNROLL(0); UNROLL(1);}
    if (vreal::WIDTH > 2) { UNROLL(2); UNROLL(3);}
    if (vreal::WIDTH > 4) { UNROLL(4); UNROLL(5); UNROLL(6); UNROLL(7);}
    if (vreal::WIDTH > 8) { 
      UNROLL( 8); UNROLL( 9); UNROLL(10); UNROLL(11);
      UNROLL(12); UNROLL(13); UNROLL(14); UNROLL(15);
    }
#undef UNROLL
  }

  solveX(n2m, n3m, rhs2, rhs2, mat);

#endif

}

/************************************************************************
 *   this subroutine performs the inversion of the q3 momentum equation
 *   by a factored implicit scheme, only the derivatives 11,22,33 of q3
 *   are treated implicitly
 *       direction x3
 ************************************************************************/

void solq3k()
{
#if 0 /* scalar */
  Tridiag<real> mat[NY];
  real cnst[NY], d[NY];

  const real betadx = beta*al;

  for (int kc = 1; kc < n3-1; kc++)
  {
    const real ackl_b = 1.0 - ac3ck(kc)*betadx;
    cnst[kc]   = 1.0/ackl_b;
    mat [kc].a = -am3ck(kc)*betadx * cnst[kc];
    mat [kc].b = 1.0;
    mat [kc].c = -ap3ck(kc)*betadx * cnst[kc];
  }
  d  [0]   = d  [n3-1]   = 0.0;
  mat[0].a = mat[n3-1].a = 0.0;
  mat[0].b = mat[n3-1].b = 1.0;
  mat[0].c = mat[n3-1].c = 0.0;

  Tridiag<real>::Reduce(n3, mat);

  const double t0 = rtc();
  for (int jc = 0; jc < n2m; jc++)
  {
    for (int kc = 1; kc < n3-1; kc++)
      d[kc] = rhs(jc,kc) * cnst[kc];

    Tridiag<real>::Solve(n3, mat, d);

    for (int kc = 0; kc < n3; kc++)
      q3(jc,kc) += d[kc];
  }
  mytimer[61] += rtc() - t0;

#else /* vector */

  static Tridiag<vreal> mat[NY];
  static vreal cnst[NY];

  const vreal betadx = beta*al;
  assert(n3m % vreal::WIDTH == 0);

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc += vreal::WIDTH)
  {
    const vreal ackl_b = vreal(1.0) - VLDA(ac3ck(kc))*betadx;
    const vreal inv =  vreal(1.0)/ackl_b;
    const vreal a   = -VLDA(am3ck(kc)) * betadx * inv;
    const vreal c   = -VLDA(ap3ck(kc)) * betadx * inv;

    for (int ch = 0; ch < vreal::WIDTH; ch++)
      mat[kc+ch].b = 1.0;

#define UNROLL(ch) {\
  mat [kc+ch].a = SIMD::broadcast<ch>(a); \
  mat [kc+ch].c = SIMD::broadcast<ch>(c); \
  cnst[kc+ch]   = SIMD::broadcast<ch>(inv); }

    assert(vreal::WIDTH >= 2);
    assert(vreal::WIDTH <= 16);
    if (vreal::WIDTH > 0) { UNROLL(0); UNROLL(1);}
    if (vreal::WIDTH > 2) { UNROLL(2); UNROLL(3);}
    if (vreal::WIDTH > 4) { UNROLL(4); UNROLL(5); UNROLL(6); UNROLL(7);}
    if (vreal::WIDTH > 8) { 
      UNROLL( 8); UNROLL( 9); UNROLL(10); UNROLL(11);
      UNROLL(12); UNROLL(13); UNROLL(14); UNROLL(15);
    }
#undef UNROLL
  }

  mat[0].a = mat[n3-1].a = 0.0;
  mat[0].b = mat[n3-1].b = 1.0;
  mat[0].c = mat[n3-1].c = 0.0;

  solveY(n2m, n3, cnst, rhs2, q3, mat, true);

#endif
}
