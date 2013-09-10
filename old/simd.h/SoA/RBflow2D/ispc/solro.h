#pragma once

/************************************************************************
 *   this subroutine performs the inversion of the q3 momentum equation
 *   by a factored implicit scheme, only the derivatives 11,22,33 of q3
 *   are treated implicitly
 *       direction x2
 ************************************************************************/
void solroj()
{
#if 0
  Tridiag<real> mat[NX];
  real u[NX], d[NX];

  const real betadx=0.5*al*dt/pec;

  for (int jc = 0; jc < n2m; jc++)
  {
    mat[jc].a =     - amscj(jc)*betadx;
    mat[jc].b = 1.0 - acscj(jc)*betadx;
    mat[jc].c =     - apscj(jc)*betadx;
  }

  Tridiag<real>::Reduce(n2m, mat, u);

  for (int kc = 0; kc < n3m; kc++)
  {
    for (int jc = 0; jc < n2m; jc++)
      d[jc] = rhs(jc,kc);

    Tridiag<real>::Solve(n2m, mat, u, d);

    for (int jc = 0; jc < n2m; jc++)
      rhs(jc,kc) = d[jc];
  }
#else
  static Tridiag<SIMD::real> mat[NX];
  
  const SIMD::real betadx = -0.5*al*dt/pec;
  assert(n2m % SIMD::WIDTH == 0);

#pragma omp parallel for
  for (int jc = 0; jc < n2m; jc += SIMD::WIDTH)
  {
    const SIMD::real a =                   SIMD::real(&amscj(jc))*betadx; 
    const SIMD::real b = SIMD::real(1.0) + SIMD::real(&acscj(jc))*betadx; 
    const SIMD::real c =                   SIMD::real(&apscj(jc))*betadx; 
#define UNROLL(ch) {\
  mat [jc+ch].a = SIMD::broadcast<ch>(a); \
  mat [jc+ch].b = SIMD::broadcast<ch>(b); \
  mat [jc+ch].c = SIMD::broadcast<ch>(c); }

    assert(SIMD::WIDTH <= 8);
    if (SIMD::WIDTH > 0) { UNROLL(0); UNROLL(1);}
    if (SIMD::WIDTH > 2) { UNROLL(2); UNROLL(3);}
    if (SIMD::WIDTH > 4) { UNROLL(4); UNROLL(5); UNROLL(6); UNROLL(7);}
#undef UNROLL
  }

  solveX(n2m, n3m, rhs3, rhs3, mat);

#endif
}

/************************************************************************
 *   this subroutine performs the inversion of the q3 momentum equation
 *   by a factored implicit scheme, only the derivatives 11,22,33 of q3
 *   are treated implicitly
 *       direction x3
 ************************************************************************/
void solrok()
{
#if 0
  Tridiag<real> mat[NY];
  real cnst[NY], d[NY];

  const real betadx=0.5*al*dt/pec;

  for (int kc = 0; kc < n3m; kc++)
  {
    const real ackl_b=1.-ac3ssk(kc)*betadx;
    cnst[kc] = 1.0/ackl_b;
    mat [kc].a = -am3ssk(kc)*betadx * cnst[kc];
    mat [kc].b = 1.0;
    mat [kc].c = -ap3ssk(kc)*betadx * cnst[kc];
  }

  Tridiag<real>::Reduce(n3m, mat);

  for (int jc = 0; jc < n2m; jc++)
  {
    for (int kc = 0; kc < n3m; kc++)
      d[kc] = rhs(jc,kc) * cnst[kc];

    Tridiag<real>::Solve(n3m, mat, d);

    for (int kc = 0; kc < n3m; kc++)
      dens(jc,kc) += d[kc];
  }
#else
  static Tridiag<SIMD::real> mat[NY];
  static SIMD::real cnst[NY];

  const SIMD::real betadx=0.5*al*dt/pec;
  assert(n3m % SIMD::WIDTH == 0);

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc += SIMD::WIDTH)
  {
    const SIMD::real ackl_b = SIMD::real(1.0) - SIMD::real(&ac3ssk(kc))*betadx;
    const SIMD::real inv =  SIMD::real(1.0)/ackl_b;
    const SIMD::real a   = -SIMD::real(&am3ssk(kc)) * betadx * inv;
    const SIMD::real c   = -SIMD::real(&ap3ssk(kc)) * betadx * inv;

    for (int ch = 0; ch < SIMD::WIDTH; ch++)
      mat[kc+ch].b = 1.0;

#define UNROLL(ch) {\
  mat [kc+ch].a = SIMD::broadcast<ch>(a); \
  mat [kc+ch].c = SIMD::broadcast<ch>(c); \
  cnst[kc+ch]   = SIMD::broadcast<ch>(inv); }

    assert(SIMD::WIDTH <= 8);
    if (SIMD::WIDTH > 0) { UNROLL(0); UNROLL(1);}
    if (SIMD::WIDTH > 2) { UNROLL(2); UNROLL(3);}
    if (SIMD::WIDTH > 4) { UNROLL(4); UNROLL(5); UNROLL(6); UNROLL(7);}
#undef UNROLL
  }

  solveY(n2m, n3m, cnst, rhs3, dens, mat);

#endif
}
