#pragma once

void solq2j()  /* passed */
{
  static Tridiag<vreal> mat[NX] __attribute__((aligned(64)));

  const vreal betadx = -beta*al;
  assert(n2m % vreal::WIDTH == 0);

#pragma omp parallel for
  for (int jc = 0; jc < n2m; jc += vreal::WIDTH)
  {
    const vreal a =              VLDA(am2j(jc))*betadx;
    const vreal b = vreal(1.0) + VLDA(ac2j(jc))*betadx; 
    const vreal c =              VLDA(ap2j(jc))*betadx; 
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

  solveX(n2m, n3m, rhs1, rhs1, mat);
}

void solq2k()  /* checked */
{
  static Tridiag<vreal> mat[NY] __attribute__((aligned(64))) ;
  
  const vreal betadx = beta*al;
  assert(n3m % vreal::WIDTH == 0);

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc += vreal::WIDTH)
  {
    const vreal a =            - VLDA(am3sk(kc)) * betadx;
    const vreal b = vreal(1.0) - VLDA(ac3sk(kc)) * betadx;
    const vreal c =            - VLDA(ap3sk(kc)) * betadx;

#define UNROLL(ch) {\
  mat [kc+ch].a = SIMD::broadcast<ch>(a); \
  mat [kc+ch].b = SIMD::broadcast<ch>(b); \
  mat [kc+ch].c = SIMD::broadcast<ch>(c);}

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

  solveY(n2m, n3m, rhs1, q2, mat);
}
