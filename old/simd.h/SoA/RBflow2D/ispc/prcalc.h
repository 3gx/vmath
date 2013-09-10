#pragma once

/***********************************************************************
 *   this subroutine performs the calculation of the pressure.
 *   this depends on the fractional step
 ***********************************************************************/
void prcalc_scal(const int kbeg, const int kend, const int jbeg, const int jend)
{
  //    the pressure is evaluated at the center of the box.
  //
  //    p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}

  const real be = al*beta;

#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {
    const int kp = kpv(kc)-1;
    const int km = kmv(kc)-1;
    const real c1 = amphk(kc);
    for (int jc = jbeg; jc < jend; jc++)
    {
      const int jm = jmv(jc)-1;
      const int jp = jpv(jc)-1;
      pr(jc,kc) += dph(jc,kc) - be*(
          (dph(jc,kp)*apphk(kc) + dph(jc,kc)*acphk (kc) + dph(jc,km)*c1       ) +
          (dph(jp,kc)*apphj(jc) + dph(jc,kc)*acphjj(jc) + dph(jm,kc)*amphj(jc)) );

    }
  }
}
void prcalc()
{
  assert(n2m % SIMD::WIDTH == 0);
  assert(n3m % SIMD::WIDTH == 0);

  prcalc_scal(0,n3m, 0,  SIMD::WIDTH);
  prcalc_scal(0,n3m, n2m-SIMD::WIDTH, n2m);

  //    the pressure is evaluated at the center of the box.
  //
  //    p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}

  const SIMD::real be = al*beta;

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc++)
  {
    const int kp = kpv(kc)-1;
    const int km = kmv(kc)-1;
    const SIMD::real ap3 = apphk(kc);
    const SIMD::real ac3 = acphk(kc);
    const SIMD::real am3 = amphk(kc);
    for (int jc = SIMD::WIDTH; jc < n2m-SIMD::WIDTH; jc += SIMD::WIDTH)
    {
      const int jm = jc-1;
      const int jp = jc+1;

      STA(pr(jc,kc)) += LDA(dph(jc,kc)) - be*(
          (LDA(dph(jc,kp))*ap3 + 
           LDA(dph(jc,kc))*ac3 +
           LDA(dph(jc,km))*am3) +

          (LDU(dph(jp,kc))*LDA(apphj (jc)) +
           LDA(dph(jc,kc))*LDA(acphjj(jc)) +
           LDU(dph(jm,kc))*LDA(amphj (jc)) ) );

    }
  }
}
