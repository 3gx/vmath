#pragma once

void solq2j()  /* passed */
{

#if 0
  Tridiag<real> mat[NX];
  real u[NX], d[NX];

  const real betadx=beta*al;

  for (int jc = 0; jc < n2m; jc++)
  {
    mat[jc].a =     - am2j(jc)*betadx;
    mat[jc].b = 1.0 - ac2j(jc)*betadx;
    mat[jc].c =     - ap2j(jc)*betadx;
  }

  Tridiag<real>::Reduce(n2m, mat, u);

  static int iter = 0; 
  static unsigned long long cnt = 0;
  const unsigned long long t1 = TIMER();
  const float FLOPS = 20;

  for (int kc = 0; kc < n3m; kc++)
  {
    for (int jc = 0; jc < n2m; jc++)
      d[jc] = rhs(jc,kc);
          
    Tridiag<real>::Solve(n2m, mat, u, d);

    for (int jc = 0; jc < n2m; jc++)
      rhs(jc,kc) = d[jc];
  }
  
  cnt += TIMER() - t1;
  iter++;
  const double dt = (double)cnt/CCLOCK;
  if (iter % 10 == 0)
    fprintf(stderr, " solq2jI: GFLOP/s= %g \n", FLOPS*208*208*iter/dt/1e9);

#else
  static Tridiag<SIMD::real> mat[NX] __attribute__((aligned(64)));

  const SIMD::real betadx=-beta*al;
  assert(n2m % SIMD::WIDTH == 0);

#pragma omp parallel for
  for (int jc = 0; jc < n2m; jc += SIMD::WIDTH)
  {
    const SIMD::real a =                   SIMD::real(&am2j(jc))*betadx; 
    const SIMD::real b = SIMD::real(1.0) + SIMD::real(&ac2j(jc))*betadx; 
    const SIMD::real c =                   SIMD::real(&ap2j(jc))*betadx; 
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

  double t0 = rtc();
  solveX(n2m, n3m, rhs1, rhs1, mat);
  mytimer[60] += rtc() - t0;

#endif
}

void solq2k()  /* checked */
{
#if 0

  Tridiag<real> mat[NY];
  real cnst[NY], d[NY];

  const real betadx = beta*al;
  for (int kc = 0; kc < n3m; kc++)
  {
    const real ackl_b = 1.0 - ac3sk(kc)*betadx;
    cnst[kc] = 1.0/ackl_b;
    mat [kc].a = -am3sk(kc)*betadx * cnst[kc];
    mat [kc].b =  1.0;
    mat [kc].c = -ap3sk(kc)*betadx * cnst[kc];
  }

  Tridiag<real>::Reduce(n3m, mat);
  
  static int iter = 0; 
  static unsigned long long cnt = 0;
  const unsigned long long t1 = TIMER();
  const float FLOPS = 10;

  for (int jc = 0; jc < n2m; jc++)
  {
    for (int kc = 0; kc < n3m; kc++)
      d[kc] = rhs(jc,kc) * cnst[kc];

    Tridiag<real>::Solve(n3m, mat, d);
    
    for (int kc = 0; kc < n3m; kc++)
      q2(jc,kc) += d[kc];
  }

  for (int kc = 0; kc < n3; kc++)
    q2(n2-1, kc) = q2(0, kc);
  
  cnt += TIMER() - t1;
  iter++;
  const double dt = (double)cnt/CCLOCK;
  if (iter % 10 == 0)
    fprintf(stderr, " solq2kI: GFLOP/s= %g \n", FLOPS*208*208*iter/dt/1e9);

#else
  static Tridiag<SIMD::real> mat[NY] __attribute__((aligned(64))) ;
  static SIMD::real cnst[NY] __attribute__((aligned(64))) ;
  
  const SIMD::real betadx = beta*al;
  assert(n3m % SIMD::WIDTH == 0);

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc += SIMD::WIDTH)
  {
    const SIMD::real ackl_b = SIMD::real(1.0) - SIMD::real(&ac3sk(kc))*betadx;
    const SIMD::real inv =  SIMD::real(1.0)/ackl_b;
    const SIMD::real a   = -SIMD::real(&am3sk(kc)) * betadx * inv;
    const SIMD::real c   = -SIMD::real(&ap3sk(kc)) * betadx * inv;

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

  double t0 = rtc();
  solveY(n2m, n3m, cnst, rhs1, q2, mat);
  mytimer[61] += rtc() - t0;

#endif

}
