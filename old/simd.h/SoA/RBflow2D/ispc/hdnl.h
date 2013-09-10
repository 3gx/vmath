#pragma once

void hdnl2scal(const int kbeg, const int kend, const int jbeg, const int jend)
{
  /*    h term for the q2 momentum equation at i+1/2,j,k+1/2     */

#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {
    const int kmm=kmv(kc)-1;
    const int kpp=kpv(kc)-1;
    const int kp  = kc+1;
    const real c1 = (-0.25)*udx3m(kc);
    for (int jc = jbeg; jc < jend; jc++)
    {
      const int jm=jmv(jc)-1;
      const int jp=jpv(jc)-1;

      /*     q2 q2 term
       *                 d  q_r q_r 
       *                ------------
       *                 d   r      
       */
      const real h22 = (-0.25)*udx2c(jc)*(SQR(q2(jp,kc)+q2(jc,kc)) - SQR(q2(jc,kc)+q2(jm,kc)));

      /*
       *     q2 q3 term
       *                 d  q_x q_r 
       *                -----------
       *                 d   x      
       */
      const real h23 = c1*(
          (q3(jc,kp) + q3(jm,kp)) * (q2(jc,kpp) + q2(jc,kc )) -
          (q3(jc,kc) + q3(jm,kc)) * (q2(jc,kc ) + q2(jc,kmm)) );
      dph(jc, kc) = h22 + h23;
    }
  }
}  /* 16 flops */

void hdnl2()
{
  assert(n2m % SIMD::WIDTH == 0);
  assert(n3m % SIMD::WIDTH == 0);

  hdnl2scal(0, n3m, 0,  SIMD::WIDTH);
  hdnl2scal(0, n3m, n2m-SIMD::WIDTH, n2m);

  /*    h term for the q2 momentum equation at i+1/2,j,k+1/2     */

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc++)
  {
    const int kmm=kmv(kc)-1;
    const int kpp=kpv(kc)-1;
    const int kp  = kc+1;
    const SIMD::real c1 = (-0.25)*udx3m(kc);

    for (int jc = SIMD::WIDTH; jc < n2m-SIMD::WIDTH; jc += SIMD::WIDTH)
    {
      const int jm=jc-1;
      const int jp=jc+1;

      /*     q2 q2 term
       *                 d  q_r q_r 
       *                ------------
       *                 d   r      
       */
      const SIMD::real h22 = SIMD::real(-0.25)*LDA(udx2c(jc)) * (
        (SQR(LDU(q2(jp,kc))+LDA(q2(jc,kc)))) - SQR(LDA(q2(jc,kc))+LDU(q2(jm,kc))) );

      /*
       *     q2 q3 term
       *                 d  q_x q_r 
       *                -----------
       *                 d   x      
       */
      const SIMD::real h23 = c1 * (
          (LDA(q3(jc,kp)) + LDU(q3(jm,kp))) * (LDA(q2(jc,kpp)) + LDA(q2(jc, kc))) -
          (LDA(q3(jc,kc)) + LDU(q3(jm,kc))) * (LDA(q2(jc, kc)) + LDA(q2(jc,kmm))) );

      STA(dph(jc,kc)) = h22 + h23;
    }
  }
}  /* 16 flops */

/*************************************************************/
/*************************************************************/
/*************************************************************/

void hdnl3scal(const int kbeg, const int kend, const int jbeg, const int jend)
{

  /*     h term for the q3 momentum equation at i+1/2,j+1/2,k   */
#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {
    const int km = kmv(kc) - 1;
    const int kp = kc + 1;
    const real c1 = (-0.25) * udx3c(kc);
    for (int jc = jbeg; jc < jend; jc++)
    {
      const int jmm = jmv(jc)-1;
      const int jpp = jpv(jc)-1;
      const int  jp = jc + 1;
      /*   q3 q2 term
       *                d  q_x q_r 
       *             -----------
       *                d   r      
       */
      const real h32 = (-0.25) * udx2c(jc) * (
          (q2(jp,kc) + q2(jp,km)) * (q3(jpp,kc) + q3(jc,  kc)) -
          (q2(jc,kc) + q2(jc,km)) * (q3(jc, kc) + q3(jmm, kc)) );
      /*   q3 q3 term
       *                 d  q_x q_x 
       *                -----------
       *                 d   x      
       */
      const real h33 = c1 * (SQR(q3(jc,kp) + q3(jc,kc)) - SQR(q3(jc,kc)+q3(jc,km)));

      /* add the buoyancy term */

      const real densit = 0.5 * (dens(jc,kc) + dens(jc,km));

      qcap(jc,kc) = h32 + h33 + densit;
    }
  }
} /* 19 flops */

void hdnl3()
{
  assert(n2m % SIMD::WIDTH == 0);
  assert(n3m % SIMD::WIDTH == 0);

  hdnl3scal(1, n3m, 0,  SIMD::WIDTH);
  hdnl3scal(1, n3m, n2m-SIMD::WIDTH, n2m);

  /*     h term for the q3 momentum equation at i+1/2,j+1/2,k   */

#pragma omp parallel for
  for (int kc = 1; kc < n3m; kc++)
  {
    const int km = kmv(kc) - 1;
    const int kp = kc + 1;
    const SIMD::real c1 = (-0.25) * udx3c(kc);
    for (int jc = SIMD::WIDTH; jc < n2m-SIMD::WIDTH; jc += SIMD::WIDTH)
    {
      const int jmm = jc-1;
      const int jpp = jc+1;
      const int  jp = jc+1;
      
      /*   q3 q2 term
       *                d  q_x q_r 
       *             -----------
       *                d   r      
       */
      const SIMD::real h32 = SIMD::real(-0.25) * LDA(udx2c(jc)) * (
          (LDU(q2(jp,kc)) + LDU(q2(jp,km)))*(LDU(q3(jpp,kc)) + LDA(q3(jc ,kc)))  -
          (LDA(q2(jc,kc)) + LDA(q2(jc,km)))*(LDA(q3(jc ,kc)) + LDU(q3(jmm,kc))) );
      
      /*   q3 q3 term
       *                 d  q_x q_x 
       *                -----------
       *                 d   x      
       */
      const SIMD::real h33 = c1 * (
          SQR(LDA(q3(jc,kp))+LDA(q3(jc,kc))) - SQR(LDA(q3(jc,kc))+LDA(q3(jc,km))) );

      /* add the buoyancy term */

      const SIMD::real densit =  SIMD::real(0.5) * (LDA(dens(jc,kc)) + LDA(dens(jc,km)));

      STA(qcap(jc,kc)) = h32 + h33 + densit;
    }
  }
} /* 19 flops */

/*************************************************************/
/*************************************************************/
/*************************************************************/

void hdnlro_scal(const int kbeg, const int kend, const int jbeg, const int jend)
{
  /*     h term for the q3 momentum equation at i+1/2,j+1/2,k    */

#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {
    const int km  = kmv(kc)-1;
    const int kpp = kpv(kc)-1;
    const int kp  = kc + 1;
    const real c1 = (-0.5)*udx3m(kc);
    for (int jc = jbeg; jc < jend; jc++)
    {
      const int jmm = jmv(jc)-1;
      const int jpp = jpv(jc)-1;
      const int jp  = jc+1;
      /*    rho q2 term
       *                d  rho q_r 
       *             -----------
       *                d   r      
       */
      const real h32 = (-0.5)*udx2m(jc) * (
          q2(jp,kc) * (dens(jpp,kc) + dens(jc, kc)) -
          q2(jc,kc) * (dens(jc ,kc) + dens(jmm,kc)) );
      /*    rho q3 term
       *                 d  rho q_x 
       *                -----------
       *                 d   x      
       */
      const real h33 = c1 * (
          q3(jc,kp) * (dens(jc,kpp) + dens(jc,kc)) - 
          q3(jc,kc) * (dens(jc,kc ) + dens(jc,km)) );

      hro(jc,kc) = h32 + h33;
    }
  }
}  /* 14 flops */

void hdnlro()
{
  assert(n2m % SIMD::WIDTH == 0);
  assert(n3m % SIMD::WIDTH == 0);

  hdnlro_scal(0, n3m, 0,  SIMD::WIDTH);
  hdnlro_scal(0, n3m, n2m-SIMD::WIDTH, n2m);

  /*     h term for the q3 momentum equation at i+1/2,j+1/2,k    */

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc++)
  {
    const int km  = kmv(kc)-1;
    const int kpp = kpv(kc)-1;
    const int kp  = kc + 1;
    const SIMD::real c1 = (-0.5)*udx3m(kc);
    for (int jc = SIMD::WIDTH; jc < n2m-SIMD::WIDTH; jc += SIMD::WIDTH)
    {
      const int jmm = jc-1;
      const int jpp = jc+1;
      const int jp  = jc+1;
      /*    rho q2 term
       *                d  rho q_r 
       *             -----------
       *                d   r      
       */
      const SIMD::real h32 = SIMD::real(-0.5)*LDA(udx2m(jc)) * (
          LDU(q2(jp,kc)) * (LDU(dens(jpp,kc)) + LDU(dens(jc,kc))) -
          LDA(q2(jc,kc)) * (LDU(dens(jmm,kc)) + LDU(dens(jc,kc))) );

      /*    rho q3 term
       *                 d  rho q_x 
       *                -----------
       *                 d   x      
       */
      const SIMD::real h33 = c1 * (
          LDA(q3(jc,kp))*(LDA(dens(jc,kpp)) + LDA(dens(jc,kc))) -
          LDA(q3(jc,kc))*(LDA(dens(jc,km )) + LDA(dens(jc,kc))) );

      STA(hro(jc,kc)) = h32 + h33;
    }
  }
}  /* 14 flops */
