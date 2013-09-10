#pragma once

template<typename REAL>
inline void hdnl2_ij(
    const int jc,  const int kc,
    const int jm,  const int km,
    const int jp,  const int kp,
    const int jmm, const int kmm,
    const int jpp, const int kpp)
{
  /*     q2 q2 term
   *                 d  q_r q_r 
   *                ------------
   *                 d   r      
   */
  const REAL h22 = REAL(-0.25)*LDA(udx2c(jc))*
    ( SQR(LDU(q2(jpp,kc))+LDA(q2(jc,kc))) - SQR(LDA(q2(jc,kc))+LDU(q2(jmm,kc))) );

  /*
   *     q2 q3 term
   *                 d  q_x q_r 
   *                -----------
   *                 d   x      
   */
  const REAL h23 = REAL(-0.25*udx3m(kc))*(
      (LDA(q3(jc,kp)) + LDU(q3(jmm,kp))) * (LDA(q2(jc,kpp)) + LDA(q2(jc,kc))) -
      (LDA(q3(jc,kc)) + LDU(q3(jmm,kc))) * (LDA(q2(jc,kc )) + LDA(q2(jc,kmm))) );
  REF(dph(jc, kc)) = h22 + h23;
}

template<typename REAL>
inline void hdnl3_ij(
    const int jc,  const int kc,
    const int jm,  const int km,
    const int jp,  const int kp,
    const int jmm, const int kmm,
    const int jpp, const int kpp)
{
  /*   q3 q2 term
   *                d  q_x q_r 
   *             -----------
   *                d   r      
   */
  const REAL h32 = REAL(-0.25) * LDA(udx2c(jc)) * (
      (LDU(q2(jp,kc)) + LDU(q2(jp,kmm))) * (LDU(q3(jpp,kc)) + LDA(q3(jc,  kc))) -
      (LDA(q2(jc,kc)) + LDA(q2(jc,kmm))) * (LDA(q3(jc, kc)) + LDU(q3(jmm, kc))) );
  /*   q3 q3 term
   *                 d  q_x q_x 
   *                -----------
   *                 d   x      
   */
  const REAL h33 = REAL(-0.25*udx3c(kc)) * 
    (SQR(LDA(q3(jc,kp)) + LDA(q3(jc,kc))) - SQR(LDA(q3(jc,kc))+LDA(q3(jc,kmm))));

  /* add the buoyancy term */

  const REAL densit = REAL(0.5) * (LDA(dens(jc,kc)) + LDA(dens(jc,kmm)));

  REF(qcap(jc,kc)) = h32 + h33 + densit;
}

template<typename REAL>
inline void hdnlr_ij(
    const int jc,  const int kc,
    const int jm,  const int km,
    const int jp,  const int kp,
    const int jmm, const int kmm,
    const int jpp, const int kpp)
{
  const REAL h32 = REAL(-0.5)*LDA(udx2m(jc)) * (
      LDU(q2(jp,kc)) * (LDU(dens(jpp,kc)) + LDA(dens(jc, kc))) -
      LDA(q2(jc,kc)) * (LDA(dens(jc ,kc)) + LDU(dens(jmm,kc))) );
  /*    rho q3 term
   *                 d  rho q_x 
   *                -----------
   *                 d   x      
   */
  const REAL h33 = REAL(-0.5*udx3m(kc)) * (
      LDA(q3(jc,kp)) * (LDA(dens(jc,kpp)) + LDA(dens(jc,kc))) - 
      LDA(q3(jc,kc)) * (LDA(dens(jc,kc )) + LDA(dens(jc,km))) );

  REF(hro(jc,kc)) = h32 + h33;
}

void hdnl()
{
  assert(n2m % vreal::WIDTH == 0);

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc++)
  {
    const int kmm=kmv(kc)-1;
    const int kpp=kpv(kc)-1;
    const int kp  = kc+1;
    const int km  = kc-1;
    for (int jc = 0; jc < n2m; jc += vreal::WIDTH)
    {
      const int jp  = jc+1;
      const int jm  = jc-1;
      if (jc < vreal::WIDTH || jc >= n2m-vreal::WIDTH)
        for (int ch = 0; ch < vreal::WIDTH; ch++)
        {
          const int jmm = jmv(jc+ch)-1;
          const int jpp = jpv(jc+ch)-1;
          hdnl2_ij<sreal>(jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
          hdnl3_ij<sreal>(jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
          hdnlr_ij<sreal>(jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
        }
      else
      {
        hdnl2_ij<vreal>(jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
        hdnl3_ij<vreal>(jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
        hdnlr_ij<vreal>(jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
      }
    }
  }
}

