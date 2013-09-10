#pragma once

/***********************************************************************
 *     this subroutine calculates the solenoidal vel field.
 *       q(n+1)=qhat-grad(dph)*dt ,  pr=dph
 *    third order runge-kutta is used.
 ***********************************************************************/

void updvp_scal(const int kbeg, const int kend, const int jbeg, const int jend)
{
  asm("#test1");
#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {
    asm("#test2");
    const int km = kmv(kc)-1;
    const real usukm = al*dt*udx3c(kc);

    //RO  Compute the q2 and q3 velocity components
    //RO  by updating with the gradient of dph

    for (int jc = jbeg; jc < jend; jc++)
    {
      asm("#test3");
      const int jm = jmv(jc)-1;
      const real  usurm = al*dt*udx2c(jc);
      q2(jc,kc) -= (dph(jc,kc)-dph(jm,kc))*usurm;
      q3(jc,kc) -= (dph(jc,kc)-dph(jc,km))*usukm;
      asm("#test4");
    }
  }

}   /* 8 flops */

void updvp()
{
  assert(n2m % SIMD::WIDTH == 0);
  assert(n3m % SIMD::WIDTH == 0);

  updvp_scal(0, n3m, 0, SIMD::WIDTH);

  const SIMD::real c0 = al*dt;

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc++)
  {
    const int km = kmv(kc)-1;
    const SIMD::real usukm = al*dt*udx3c(kc);

    //RO  Compute the q2 and q3 velocity components
    //RO  by updating with the gradient of dph

    for (int jc = SIMD::WIDTH; jc < n2m; jc += SIMD::WIDTH)
    {
      const int jm = jc-1;
      const SIMD::real usurm = c0*LDA(udx2c(jc));
      STA(q2(jc,kc)) -= (LDA(dph(jc,kc)) - LDU(dph(jm,kc)))*usurm;
      STA(q3(jc,kc)) -= (LDA(dph(jc,kc)) - LDA(dph(jc,km)))*usukm;
    }
  }

  for (int kc = 0; kc < NY; kc++)
    q2(n2-1,kc) = q2(0, kc);
  for (int jc = 0; jc < NX; jc++)
    q3(jc,0) = q3(jc,n3-1) = 0.0;
}   /* 8 flops */
