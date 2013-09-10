#pragma once

void divg_scal(const int kbeg, const int kend, const int jbeg, const int jend)
{
  const real usdtal = 1.0/(dt*al);

#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {
    const int kp = kc + 1;
    const real c1 = udx3m(kc);
    for (int jc = jbeg; jc < jend; jc++)
    {
      const int jp = jpv(jc) - 1;
      const real dqcap = (q2(jp,kc) - q2(jc,kc))*udx2m(jc) + (q3(jc,kp)-q3(jc,kc))*c1;
      dph(jc,kc)=dqcap*usdtal;
    }
  }

}

void divg()
{
  assert(n2m % SIMD::WIDTH == 0);
  assert(n3m % SIMD::WIDTH == 0);
  divg_scal(0, n3m, n2m-SIMD::WIDTH, n2m);

  const SIMD::real usdtal = 1.0/(dt*al);

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc++)
  {
    const int kp = kc + 1;
    const SIMD::real c1 = udx3m(kc);
    for (int jc = 0; jc < n2m-SIMD::WIDTH; jc += SIMD::WIDTH)
    {
      const int jp = jc+1;
      const SIMD::real dqcap = 
        (LDU(q2(jp,kc)) - LDA(q2(jc,kc)))*LDA(udx2m(jc)) +
        (LDA(q3(jc,kp)) - LDA(q3(jc,kc)))*c1;

      STA(dph(jc,kc)) = dqcap * usdtal;
    }
  }

}
