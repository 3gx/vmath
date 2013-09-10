#include "defs_ispc.h"

#define d data

inline static void
invtr2_ij(
    const real alre,
    const int kc, const int jc,
    const int km, const int jm,
    const int kp, const int jp,
    const int kmm, const int jmm,
    const int kpp, const int jpp)
{
  /*
   *  compute the rhs of the factored equation
   *  everything at i+1/2,j,k+1/2
   *
   *    points inside the flowfield
   */
  const real udx2 = d.al*d.udx2c[jc]; /* y */

  /*
   *   22 second derivative of q2
   */
  const real d22q2 = 
    d.q2[kc][jpp] * d.ap2je[jc] +  
    d.q2[kc][jc ] * d.ac2je[jc] +
    d.q2[kc][jmm] * d.am2je[jc];  /* y */

  /*
   *   33 second derivative of q2
   */
  const real d33q2 =
    d.q2[kpp][jc] * d.ap3sk[kc] +
    d.q2[kc ][jc] * d.ac3sk[kc] +
    d.q2[kmm][jc] * d.am3sk[kc]; /* y */

  /*
   *    viscid terms
   */
  const real dcq2 = d22q2 + d33q2;

  /*
   *   component of grad(pr) along 2 direction
   */
  const real dpx22 = (d.pr[kc][jc] - d.pr[kc][jmm]) * udx2;

  d.rhs1[kc][jc] = (d.ga*d.dph[kc][jc] + d.ro*d.ru2[kc][jc] + alre*dcq2 - dpx22) * d.dt; /* y */

  d.ru2[kc][jc] = d.dph[kc][jc]; /* y */
} /* 21 flops */

inline static void
invtr3_ij(
    const real alre,
    const int kc, const int jc,
    const int km, const int jm,
    const int kp, const int jp,
    const int kmm, const int jmm,
    const int kpp, const int jpp)
{
  /*  compute the rhs of the factored equation
   *  everything at i+1/2,j+1/2,k
   */

  const real udx3 = d.al*d.udx3c[kc];  /* y */
  const real aap = d.ap3j[jc];
  const real aam = d.am3j[jc];
  const real aac = d.ac3j[jc];

  /*   22 second derivatives of q3 */
  const real dq32=aam*d.q3[kc][jmm] + aac*d.q3[kc][jc] + aap*d.q3[kc][jpp]; /* y */

  /*   33 second derivatives of q3 */
  const real dq33=d.q3[kp][jc]*d.ap3ck[kc] + d.q3[kc][jc]*d.ac3ck[kc] + d.q3[kmm][jc]*d.am3ck[kc]; /* y */

  /*   viscous terms    */
  const real dcq3=dq32+dq33; /* y */

  /*  component of grad(pr) along x3 direction */
  const real dpx33=(d.pr[kc][jc] - d.pr[kmm][jc])*udx3; /* y */

  /*======================================================= */
  d.rhs2[kc][jc]=(d.ga*d.qcap[kc][jc] + d.ro*d.ru3[kc][jc] + alre*dcq3-dpx33)*d.dt; /* y */
  /*======================================================= */

  /*  updating of the non-linear terms */
  d.ru3[kc][jc] = d.qcap[kc][jc];  /* y */
}  /* 20 flops */

inline static void
invtrro_ij(
    const real alpec,
    const int kc, const int jc,
    const int km, const int jm,
    const int kp, const int jp,
    const int kmm, const int jmm,
    const int kpp, const int jpp)
{

  const real aap = d.apscj[jc];
  const real aam = d.amscj[jc];
  const real aac = d.acscj[jc];

  /*   22 second derivatives of dens  */
  const real dq32 = d.dens[kc][jpp]*aap+d.dens[kc][jc]*aac + d.dens[kc][jmm]*aam; /* y */

  /*   33 second derivatives of dens   */
  const real dq33=
    d.dens[kpp][jc] * d.ap3ssk[kc] + 
    d.dens[kc ][jc] * d.ac3ssk[kc] +
    d.dens[kmm][jc] * d.am3ssk[kc];  /* y */

  const real dcq3=dq32+dq33; /* y */

  /*    right hand side of the density equation */
  /*m===========================================================*/
  d.rhs3[kc][jc]=(d.ga*d.hro[kc][jc]+d.ro*d.ruro[kc][jc] + alpec*dcq3)*d.dt;  /* y */
  /*m===========================================================*/

  /*    updating of the non-linear terms */
  d.ruro[kc][jc]=d.hro[kc][jc];  /* y */
}

inline static void
invtrroL_ij(
    const real alpec,
    const real del1, 
    const real del2, 
    const real fcder,
    const int kc, const int jc,
    const int km, const int jm,
    const int kp, const int jp,
    const int kmm, const int jmm,
    const int kpp, const int jpp)
{
  /*
   *       LOWER HOT WALL      
   */

  const real aap = d.apscj[jc];
  const real aam = d.amscj[jc];
  const real aac = d.acscj[jc];

  /*   22 second derivatives of dens */
  const real dq32 = d.dens[kc][jpp]*aap + d.dens[kc][jc]*aac + d.dens[kc][jmm]*aam; /* y */

  /*   33 second derivatives of dens */
  const real dq33=(d.dens[kp][jc]*del1-d.dens[kc][jc]*(del1+del2) + d.denbs[jc]*del2) *fcder; /* y */
  const real dcq3=dq32+dq33;
  /*    right hand side of the density equation */
  //m=======================================================
  d.rhs3[kc][jc] = (d.ga*d.hro[kc][jc] + d.ro*d.ruro[kc][jc] + alpec*dcq3)*d.dt; /* y */
  //cm=======================================================
  //
  //   updating of the non-linear terms
  //
  d.ruro[kc][jc] = d.hro[kc][jc];  /* y */
}

inline static void
invtrroU_ij(
    const real alpec,
    const real del1,
    const real del2,
    const real fcder,
    const int kc, const int jc,
    const int km, const int jm,
    const int kp, const int jp,
    const int kmm, const int jmm,
    const int kpp, const int jpp)
{
  /*
   *       UPPER COLD WALL
   */     


  const real aap = d.apscj[jc];
  const real aam = d.amscj[jc];
  const real aac = d.acscj[jc];

  /*   22 second derivatives of dens  */
  const real dq32 = d.dens[kc][jpp]*aap + d.dens[kc][jc]*aac + d.dens[kc][jmm] * aam; /* y */

  /*   33 second derivatives of dens  */
  const real dq33 = (d.dens[km][jc]*del1 - d.dens[kc][jc]*(del1+del2) + d.denbn[jc]*del2) * fcder; /* y */
  const real dcq3 = dq32+dq33; /* y */

  //
  //    right hand side of the density equation
  //
  //m========================================================
  d.rhs3[kc][jc]=(d.ga*d.hro[kc][jc] + d.ro*d.ruro[kc][jc] + alpec*dcq3)*d.dt; /* y */
  //m========================================================
  //
  //    updating of the non-linear terms
  //
  d.ruro[kc][jc]=d.hro[kc][jc]; /* y */
}

  task void
invtr_all_step()
{
  const uniform int n3m = data.n3m;
  const uniform int n2m = data.n2m;
  const uniform int k_per_task = (n3m-1)/taskCount + 1;
  const uniform int kbeg =  k_per_task * taskIndex;
  const uniform int kend = min(kbeg + k_per_task, n3m);

  const real alre  = data.al/data.ren;
  const real alpec = data.al/data.pec;
  for (uniform int kc = kbeg; kc < kend; kc++)
  {
    const int kp  = kc + 1;
    const int km  = kc - 1;
    const int kmm = kc == 0     ? kc : km;
    const int kpp = kc == n3m-1 ? kc : kp;

    /*******************/ 
    /*     invtr2      */
    /*******************/ 

    foreach (jc = 0 ... n2m)
    {
      const int jm = jc-1;
      const int jp = jc+1;
      cif (jc < programCount || jc >= n2m - programCount)
      {
        const int jmm = jc == 0       ? n2m-1 : jm;
        const int jpp = jc == n2m - 1 ?  0    : jp;
        invtr2_ij(alre, kc, jc, km, jm, kp, jp, kmm, jmm, kpp, jpp);
      }
      else
        invtr2_ij(alre, kc, jc, km, jm, kp, jp, kmm, jm, kpp, jp);
    }

    /*******************/ 
    /*     invtr3      */
    /*******************/ 

    foreach (jc = 0 ... n2m)
    {
      const int jm = jc-1;
      const int jp = jc+1;
      cif (jc < programCount || jc >= n2m - programCount)
      {
        const int jmm = jc == 0       ? n2m-1 : jm;
        const int jpp = jc == n2m - 1 ?  0    : jp;
        invtr3_ij(alre, kc, jc, km, jm, kp, jp, kmm, jmm, kpp, jpp);
      }
      else
        invtr3_ij(alre, kc, jc, km, jm, kp, jp, kmm, jm, kpp, jp);
    }

    /*******************/ 
    /*      invtro     */
    /*******************/ 

    if (kc >= 1 && kc <= n3m-2)
    {
      foreach (jc = 0 ... n2m)
      {
        const int jm = jc-1;
        const int jp = jc+1;
        cif (jc < programCount || jc >= n2m - programCount)
        {
          const int jmm = jc == 0       ? n2m-1 : jm;
          const int jpp = jc == n2m - 1 ?  0    : jp;
          invtrro_ij(alpec, kc, jc, km, jm, kp, jp, kmm, jmm, kpp, jpp);
        }
        else
          invtrro_ij(alpec, kc, jc, km, jm, kp, jp, kmm, jm,  kpp, jp);
      }
    }
    else if (kc == 0)
    {
      const real del1 = d.zm[0] - d.zz[0];
      const real del2 = d.zm[1] - d.zm[0];
      const real fcder = 2.0/(del1*del2*(del1+del2));

      foreach (jc = 0 ... n2m)
      {
        const int jm = jc-1;
        const int jp = jc+1;

        cif (jc < programCount || jc >= n2m - programCount)
        {
          const int jmm = jc == 0       ? n2m-1 : jm;
          const int jpp = jc == n2m - 1 ?  0    : jp;
          invtrroL_ij(alpec, del1, del2, fcder, kc, jc, km, jm, kp, jp, kmm, jmm, kpp, jpp);
        }
        else
          invtrroL_ij(alpec, del1, del2, fcder, kc, jc, km, jm, kp, jp, kmm, jm,  kpp, jp);
      }
    }
    else 
    {
      const real  del1 = d.zz[d.n3m  ] - d.zm[d.n3m-1];
      const real  del2 = d.zm[d.n3m-1] - d.zm[d.n3m-2];
      const real fcder = 2.0/(del1*del2*(del1+del2));

      foreach (jc = 0 ... n2m)
      {
        const int jm = jc-1;
        const int jp = jc+1;

        cif (jc < programCount || jc >= n2m - programCount)
        {
          const int jmm = jc == 0       ? n2m-1 : jm;
          const int jpp = jc == n2m - 1 ?  0    : jp;
          invtrroU_ij(alpec, del1, del2, fcder, kc, jc, km, jm, kp, jp, kmm, jmm, kpp, jpp);
        }
        else
          invtrroU_ij(alpec, del1, del2, fcder, kc, jc, km, jm, kp, jp, kmm, jm,  kpp, jp);
      }
    }
  }  /* kc loop */
}

  export void
invtr_all()
{   
  launch[128] invtr_all_step();
}


