#pragma once

template<typename REAL>
inline void update_qp_ij(
    const int jc, const int kc,
    const int jm, const int km,
    const int jp, const int kp)
{
  /***********************************************************************
   *     this subroutine calculates the solenoidal vel field.
   *       q(n+1)=qhat-grad(dph)*dt ,  pr=dph
   *    third order runge-kutta is used.
   ***********************************************************************/
  const REAL usukm = al*dt*udx3c(kc);

  //RO  Compute the q2 and q3 velocity components
  //RO  by updating with the gradient of dph

  const REAL usurm = REAL(al*dt)*LDA(udx2c(jc));
  REF(q2(jc,kc)) -= (LDA(dph(jc,kc)) - LDU(dph(jm,kc)))*usurm;
  REF(q3(jc,kc)) -= (LDA(dph(jc,kc)) - LDA(dph(jc,km)))*usukm;

  /***********************************************************************
   *   this subroutine performs the calculation of the pressure.
   *   this depends on the fractional step
   ***********************************************************************/
  const REAL be = al*beta;
  REF(pr(jc,kc)) += LDA(dph(jc,kc)) - be*(
      (LDA(dph(jc,kp)) * REAL(apphk(kc)) +
       LDA(dph(jc,kc)) * REAL(acphk(kc)) +
       LDA(dph(jc,km)) * REAL(amphk(kc)) ) +
      (LDU(dph(jp,kc)) * LDA(apphj (jc)) + 
       LDA(dph(jc,kc)) * LDA(acphjj(jc)) + 
       LDU(dph(jm,kc)) * LDA(amphj (jc))) );
}

void update_qp()
{
  assert(n2m % vreal::WIDTH == 0);

#pragma omp parallel 
  {
    const int nt = omp_get_num_threads();
#pragma omp for schedule(dynamic, 2*nt) nowait
    for (int kc = 0; kc < n3m; kc++)
    {
      const int kp = kpv(kc) - 1;
      const int km = kmv(kc) - 1;
      foreach_warp(jc, 0, n2m)
        if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
          foreach_lane(ch)
          {
            const int jp = jpv(jc+ch) - 1;
            const int jm = jmv(jc+ch) - 1;
            update_qp_ij<sreal>(jc+ch, kc, jm, km, jp, kp);
          }
        else
          update_qp_ij<vreal>(jc, kc, jc-1, km, jc+1, kp);
    }
  }

  for (int kc = 0; kc < n3; kc++)
    q2(n2-1,kc) = q2(0, kc);
  for (int jc = 0; jc < n2; jc++)
    q3(jc,0) = q3(jc,n3-1) = 0.0;

}

