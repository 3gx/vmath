#pragma once

/************************************************************************
 *   this subroutine performs the inversion of the q3 momentum equation
 *   by a factored implicit scheme, only the derivatives 11,22,33 of q3
 *   are treated implicitly
 *       direction x2
 ************************************************************************/
void solroj()
{
  static Tridiag<vreal> mat[NX];
  
  const vreal betadx = -0.5*al*dt/pec;
  assert(n2m % vreal::WIDTH == 0);

#pragma omp parallel for
  for (int jc = 0; jc < n2m; jc += vreal::WIDTH)
  {
    const vreal a =                   VLDA(amscj(jc))*betadx; 
    const vreal b = vreal(1.0) + VLDA(acscj(jc))*betadx; 
    const vreal c =                   VLDA(apscj(jc))*betadx; 
#define UNROLL(ch) {\
  mat [jc+ch].a = SIMD::broadcast<ch>(a); \
  mat [jc+ch].b = SIMD::broadcast<ch>(b); \
  mat [jc+ch].c = SIMD::broadcast<ch>(c); }

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

  solveX(n2m, n3m, rhs3, rhs3, mat);

}

/************************************************************************
 *   this subroutine performs the inversion of the q3 momentum equation
 *   by a factored implicit scheme, only the derivatives 11,22,33 of q3
 *   are treated implicitly
 *       direction x3
 ************************************************************************/
void solrok()
{
  static Tridiag<vreal> mat[NY];
  static vreal cnst[NY];

  const vreal betadx=0.5*al*dt/pec;
  assert(n3m % vreal::WIDTH == 0);

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc += vreal::WIDTH)
  {
    const vreal ackl_b = vreal(1.0) - VLDA(ac3ssk(kc))*betadx;
    const vreal inv =  vreal(1.0)/ackl_b;
    const vreal a   = -VLDA(am3ssk(kc)) * betadx * inv;
    const vreal c   = -VLDA(ap3ssk(kc)) * betadx * inv;

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

  solveY(n2m, n3m, cnst, rhs3, dens, mat);

}
