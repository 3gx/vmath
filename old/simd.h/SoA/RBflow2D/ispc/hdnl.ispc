#include "defs_ispc.h"

inline static void 
hdnl_all_ij(
    const int kc, const int jc,
    const int km, const int jm,
    const int kp, const int jp,
    const int kmm, const int jmm,
    const int kpp, const int jpp)
{
  const uniform Data &d = data;
  const real ca1 = (-0.25) * d.udx3m[kc];  /* y */
  const real cb1 = (-0.25) * d.udx3c[kc];  /* y */
  const real cc1 = (-0.5 ) * d.udx3m[kc];  /* y */

  /*     q2 q2 term
   *                 d  q_r q_r 
   *                ------------
   *                 d   r      
   */
  const real ha22 = (-0.25)*d.udx2c[jc]*
    (SQR(d.q2[kc][jpp] + d.q2[kc][jc]) - SQR(d.q2[kc][jc] + d.q2[kc][jmm]));  /* y */

  /*
   *     q2 q3 term
   *                 d  q_x q_r 
   *                -----------
   *                 d   x      
   */
  const real ha23 = ca1*(
      (d.q3[kp][jc] + d.q3[kp][jmm]) * (d.q2[kpp][jc] + d.q2[kc ][jc]) -
      (d.q3[kc][jc] + d.q3[kc][jmm]) * (d.q2[kc ][jc] + d.q2[kmm][jc]) );  /* y */
  d.dph[kc][jc] = ha22 + ha23; /* y */

  /*     h term for the q3 momentum equation at i+1/2,j+1/2,k   */

  /*   q3 q2 term
   *                d  q_x q_r 
   *             -----------
   *                d   r      
   */
  const real hb32 = (-0.25) * d.udx2c[jc] * (
      (d.q2[kc][jp] + d.q2[kmm][jp]) * (d.q3[kc][jpp] + d.q3[kc][jc ]) -
      (d.q2[kc][jc] + d.q2[kmm][jc]) * (d.q3[kc][jc ] + d.q3[kc][jmm]) ); /* y */
  /*   q3 q3 term
   *                 d  q_x q_x 
   *                -----------
   *                 d   x      
   */
  const real hb33 = cb1 * 
    (SQR(d.q3[kp][jc] + d.q3[kc][jc]) - SQR(d.q3[kc][jc] + d.q3[kmm][jc])); /* y */

  /* add the buoyancy term */

  const real densit = 0.5 * (d.dens[kc][jc] + d.dens[kmm][jc]); /* y */

  d.qcap[kc][jc] = hb32 + hb33 + densit;  /* y */

  /*     h term for the q3 momentum equation at i+1/2,j+1/2,k    */

  /*    rho q2 term
   *                d  rho q_r 
   *             -----------
   *                d   r      
   */
  const real hc32 = (-0.5)*d.udx2m[jc] * (
      d.q2[kc][jp] * (d.dens[kc][jpp] + d.dens[kc][jc ]) -
      d.q2[kc][jc] * (d.dens[kc][jc ] + d.dens[kc][jmm]) );  /* y */
  /*    rho q3 term
   *                 d  rho q_x 
   *                -----------
   *                 d   x      
   */
  const real hc33 = cc1 * (
      d.q3[kp][jc] * (d.dens[kpp][jc] + d.dens[kc ][jc]) -
      d.q3[kc][jc] * (d.dens[kc ][jc] + d.dens[kmm][jc]) );  /* y */

  d.hro[kc][jc] = hc32 + hc33;  /* y */
} 

task void
hdnl_all_step()
{
  const uniform int n2m = data.n2m;
  const uniform int n3m = data.n3m;
  const uniform int k_per_task = (n3m-1)/taskCount + 1;
  const uniform int kbeg =  k_per_task * taskIndex;
  const uniform int kend = min(kbeg + k_per_task, n3m);

  for (uniform int kc = kbeg; kc < kend; kc++)
  {
    /*    h term for the q2 momentum equation at i+1/2,j,k+1/2     */
    const int kp  = kc + 1;
    const int km  = kc - 1;
    const int kmm = kc == 0     ? kc : km;
    const int kpp = kc == n3m-1 ? kc : kp;

#if 1
    foreach (jc = 0 ... n2m)
    {
      const int jm = jc-1;
      const int jp = jc+1;

      cif (jc < programCount || jc >= n2m - programCount)
      {
        const int jmm = jc == 0       ? n2m-1 : jm;
        const int jpp = jc == n2m - 1 ?  0    : jp;
        hdnl_all_ij(kc, jc, km, jm, kp, jp, kmm, jmm, kpp, jpp);
      }
      else
        hdnl_all_ij(kc, jc, km, jm, kp, jp, kmm, jm, kpp, jp);
    }
#else
    foreach (jc = 0 ... programCount)
    {
      const int jm = jc-1;
      const int jp = jc+1;
      const int jmm = jc == 0       ? n2m-1 : jm;
      const int jpp = jc == n2m - 1 ?  0    : jp;

      hdnl_all_ij(kc, jc, km, jm, kp, jp, kmm, jmm, kpp, jpp);
    }

    /* this one permits vectorized accesses */ 
    foreach (jc = programCount ... n2m - programCount)
    {
      const int jm = jc-1;
      const int jp = jc+1;
      hdnl_all_ij(kc, jc, km, jm, kp, jp, kmm, jm, kpp, jp);
    }

    foreach (jc = n2m - programCount ... n2m)
    {
      const int jm = jc-1;
      const int jp = jc+1;

      const int jmm = jc == 0       ? n2m-1 : jm;
      const int jpp = jc == n2m - 1 ?  0    : jp;
      hdnl_all_ij(kc, jc, km, jm, kp, jp, kmm, jmm, kpp, jpp);
    }
#endif
  }
}


  export void
hdnl_all()
{   
  launch[128] hdnl_all_step();
}
