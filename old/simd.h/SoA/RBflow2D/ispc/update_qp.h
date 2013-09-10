#pragma once

void update_qp_scal(const int kbeg, const int kend, const int jbeg, const int jend)
{
#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {
    {
      /***********************************************************************
       *     this subroutine calculates the solenoidal vel field.
       *       q(n+1)=qhat-grad(dph)*dt ,  pr=dph
       *    third order runge-kutta is used.
       ***********************************************************************/
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
      }
    }

    /***********************************************************************
     *   this subroutine performs the calculation of the pressure.
     *   this depends on the fractional step
     ***********************************************************************/
    {
      const real be = al*beta;
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

}   /* 8 flops */

void update_qp()
{
  assert(n2m % SIMD::WIDTH == 0);
  assert(n3m % SIMD::WIDTH == 0);

  update_qp_scal(0,n3m, 0,  SIMD::WIDTH);
  update_qp_scal(0,n3m, n2m-SIMD::WIDTH, n2m);


#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc++)
  {

    {
      const SIMD::real c0 = al*dt;
      const int km = kmv(kc)-1;
      const SIMD::real usukm = al*dt*udx3c(kc);

      //RO  Compute the q2 and q3 velocity components
      //RO  by updating with the gradient of dph

      for (int jc = SIMD::WIDTH; jc < n2m-SIMD::WIDTH; jc += SIMD::WIDTH)
      {
        const int jm = jc-1;
        const SIMD::real usurm = c0*LDA(udx2c(jc));
        STA(q2(jc,kc)) -= (LDA(dph(jc,kc)) - LDU(dph(jm,kc)))*usurm;
        STA(q3(jc,kc)) -= (LDA(dph(jc,kc)) - LDA(dph(jc,km)))*usukm;
      }
    }

    {
      const SIMD::real be = al*beta;
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

  for (int kc = 0; kc < NY; kc++)
    q2(n2-1,kc) = q2(0, kc);
  for (int jc = 0; jc < NX; jc++)
    q3(jc,0) = q3(jc,n3-1) = 0.0;

}
