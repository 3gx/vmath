#pragma once

template<typename REAL>
inline void invtr2_ij(
    const REAL alre,
    const int jc, const int kc,
    const int jm, const int km,
    const int jp, const int kp,
    const int jmm, const int kmm,
    const int jpp, const int kpp)
{
  const REAL udx2 = REAL(al)*LDA(udx2c(jc));

  /*
   *   22 second derivative of q2
   */
  const REAL d22q2 = 
    LDU(q2(jpp,kc)) * LDA(ap2je(jc)) + 
    LDA(q2(jc ,kc)) * LDA(ac2je(jc)) +
    LDU(q2(jmm,kc)) * LDA(am2je(jc));

  /*
   *   33 second derivative of q2
   */
  const REAL d33q2 =
    LDA(q2(jc,kpp)) * REAL(ap3sk(kc)) +
    LDA(q2(jc,kc )) * REAL(ac3sk(kc)) +
    LDA(q2(jc,kmm)) * REAL(am3sk(kc));

  /*
   *    viscid terms
   */
  const REAL dcq2 = d22q2 + d33q2;

  /*
   *   component of grad(pr) along 2 direction
   */
  const REAL dpx22 = (LDA(pr(jc,kc)) - LDU(pr(jmm,kc))) * udx2;

  REF(rhs1(jc,kc)) = REAL(dt)*(
      REAL(ga)*LDA(dph(jc,kc)) + 
      REAL(ro)*LDA(ru2(jc,kc)) + 
      alre*dcq2 - dpx22);

  REF(ru2(jc,kc)) = LDA(dph(jc,kc));
}

template<typename REAL>
inline void invtr3_ij(
    const REAL alre,
    const int jc, const int kc,
    const int jm, const int km,
    const int jp, const int kp,
    const int jmm, const int kmm,
    const int jpp, const int kpp)
{
  const REAL udx3 = al*udx3c(kc);
  const REAL aap  = LDA(ap3j(jc));
  const REAL aam  = LDA(am3j(jc));
  const REAL aac  = LDA(ac3j(jc));

  /*   22 second derivatives of q3 */
  const REAL dq32=aam*LDU(q3(jmm,kc))+aac*LDA(q3(jc,kc))+aap*LDU(q3(jpp,kc));

  /*   33 second derivatives of q3 */
  const REAL dq33=
    LDA(q3(jc,kp )) * REAL(ap3ck(kc)) + 
    LDA(q3(jc,kc )) * REAL(ac3ck(kc)) + 
    LDA(q3(jc,kmm)) * REAL(am3ck(kc));

  /*   viscous terms    */
  const REAL dcq3=dq32+dq33;

  /*  component of grad(pr) along x3 direction */
  const REAL dpx33=(LDA(pr(jc,kc)) - LDA(pr(jc,kmm)))*udx3;

  /*======================================================= */
  REF(rhs2(jc,kc)) = REAL(dt)*(
      REAL(ga) * LDA(qcap(jc,kc)) + 
      REAL(ro) * LDA(ru3 (jc,kc))
      + alre*dcq3-dpx33);
  /*======================================================= */

  /*  updating of the non-linear terms */
  REF(ru3(jc,kc)) = LDA(qcap(jc,kc));
}  /* 20 flops */

template<typename REAL>
inline void invtrro_ij(
    const REAL alpec,
    const int jc, const int kc,
    const int jm, const int km,
    const int jp, const int kp,
    const int jmm, const int kmm,
    const int jpp, const int kpp)
{
  const REAL aap = LDA(apscj(jc));
  const REAL aam = LDA(amscj(jc));
  const REAL aac = LDA(acscj(jc));

  /*   22 second derivatives of dens  */
  const REAL dq32 = LDU(dens(jpp,kc))*aap+LDA(dens(jc,kc))*aac + LDU(dens(jmm,kc))*aam;

  /*   33 second derivatives of dens   */
  const REAL dq33=
    LDA(dens(jc,kpp)) * REAL(ap3ssk(kc)) + 
    LDA(dens(jc,kc )) * REAL(ac3ssk(kc)) + 
    LDA(dens(jc,kmm)) * REAL(am3ssk(kc));

  const REAL dcq3=dq32+dq33;

  /*    right hand side of the density equation */
  /*m===========================================================*/
  REF(rhs3(jc,kc)) = REAL(dt)*(
      REAL(ga)*LDA(hro (jc,kc)) + 
      REAL(ro)*LDA(ruro(jc,kc)) + alpec*dcq3);
  /*m===========================================================*/

  /*    updating of the non-linear terms */
  REF(ruro(jc,kc)) = LDA(hro(jc,kc));
}

template<typename REAL>
inline void invtrroL_ij(
    const REAL alpec,
    const REAL del1, 
    const REAL del2, 
    const REAL fcder,
    const int jc, const int kc,
    const int jm, const int km,
    const int jp, const int kp,
    const int jmm, const int kmm,
    const int jpp, const int kpp)
{
  const REAL aap = LDA(apscj(jc));
  const REAL aam = LDA(amscj(jc));
  const REAL aac = LDA(acscj(jc));

  /*   22 second derivatives of dens */
  const REAL dq32 = LDU(dens(jpp,kc))*aap + LDA(dens(jc,kc))*aac + LDU(dens(jmm,kc))*aam;

  /*   33 second derivatives of dens */
  const REAL dq33=(LDA(dens(jc,kp))*del1-LDA(dens(jc,kc))*(del1+del2) + LDA(denbs(jc))*del2) *fcder;
  const REAL dcq3=dq32+dq33;
  /*    right hand side of the density equation */
  //m=======================================================
  REF(rhs3(jc,kc)) = REAL(dt)*(
      REAL(ga)*LDA(hro (jc,kc)) + 
      REAL(ro)*LDA(ruro(jc,kc)) + alpec*dcq3);
  //cm=======================================================
  //
  //   updating of the non-linear terms
  //
  REF(ruro(jc,kc)) = LDA(hro(jc,kc));
}

template<typename REAL>
inline void invtrroU_ij(
    const REAL alpec,
    const REAL del1, 
    const REAL del2, 
    const REAL fcder,
    const int jc, const int kc,
    const int jm, const int km,
    const int jp, const int kp,
    const int jmm, const int kmm,
    const int jpp, const int kpp)
{
  const REAL aap = LDA(apscj(jc));
  const REAL aam = LDA(amscj(jc));
  const REAL aac = LDA(acscj(jc));

  /*   22 second derivatives of dens  */
  const REAL dq32 = LDU(dens(jpp,kc))*aap + LDA(dens(jc,kc))*aac + LDU(dens(jmm,kc))*aam;

  /*   33 second derivatives of dens  */
  const REAL dq33 = (LDA(dens(jc,km))*del1 - LDA(dens(jc,kc))*(del1+del2) + LDA(denbn(jc))*del2) * fcder;
  const REAL dcq3 = dq32+dq33;

  //
  //    right hand side of the density equation
  //
  //m========================================================
  REF(rhs3(jc,kc)) = REAL(dt)*(
      REAL(ga)*LDA(hro (jc,kc)) + 
      REAL(ro)*LDA(ruro(jc,kc)) + alpec*dcq3);
  //m========================================================
  //
  //    updating of the non-linear terms
  //
  REF(ruro(jc,kc)) = LDA(hro(jc,kc));
}

void invtr()
{
  assert(n2m % vreal::WIDTH == 0);

  const real alre  = al/ren;
  const real alpec = al/pec;

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc++)
  {
    const int kp  = kc + 1;
    const int km  = kc - 1;
    const int kmm = kc == 0     ? kc : km;
    const int kpp = kc == n3m-1 ? kc : kp;

    /* invtr2 */

    foreach_warp(jc, 0, n2m)
    {
      const int jp  = jc+1;
      const int jm  = jc-1;
      if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
        foreach_lane(ch)
        {
          const int jmm = jmv(jc+ch)-1;
          const int jpp = jpv(jc+ch)-1;
          invtr2_ij<sreal>(alre, jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
        }
      else
        invtr2_ij<vreal>(alre, jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
    }

    /* invtr3 */

    foreach_warp(jc, 0, n2m)
    {
      const int jp  = jc+1;
      const int jm  = jc-1;
      if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
        foreach_lane(ch)
        {
          const int jmm = jmv(jc+ch)-1;
          const int jpp = jpv(jc+ch)-1;
          invtr3_ij<sreal>(alre, jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
        }
      else
          invtr3_ij<vreal>(alre, jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
    }

    /* invtro */

    if (kc >= 1 && kc <= n3m-2)
    {
      foreach_warp(jc, 0, n2m)
      {
        const int jp  = jc+1;
        const int jm  = jc-1;
        if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
          foreach_lane(ch)
          {
            const int jmm = jmv(jc+ch)-1;
            const int jpp = jpv(jc+ch)-1;
            invtrro_ij<sreal>(alpec, jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
          }
        else
          invtrro_ij<vreal>(alpec, jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
      }
    }
    else if (kc == 0)
    {
      const real del1 = zm(0) - zz(0);
      const real del2 = zm(1) - zm(0);
      const real fcder = 2.0/(del1*del2*(del1+del2));
      foreach_warp(jc, 0, n2m)
      {
        const int jp  = jc+1;
        const int jm  = jc-1;
        if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
          foreach_lane(ch)
          {
            const int jmm = jmv(jc+ch)-1;
            const int jpp = jpv(jc+ch)-1;
            invtrroL_ij<sreal>(alpec, del1, del2, fcder, jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
          }
        else
          invtrroL_ij<vreal>(alpec, del1, del2, fcder, jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
      }
    }
    else
    {
      const real  del1 = zz(n3m  ) - zm(n3m-1);
      const real  del2 = zm(n3m-1) - zm(n3m-2);
      const real fcder = 2.0/(del1*del2*(del1+del2));
      foreach_warp(jc, 0, n2m)
      {
        const int jp  = jc+1;
        const int jm  = jc-1;
        if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
          foreach_lane(ch)
          {
            const int jmm = jmv(jc+ch)-1;
            const int jpp = jpv(jc+ch)-1;
            invtrroU_ij<sreal>(alpec, del1, del2, fcder, jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
          }
        else
          invtrroU_ij<vreal>(alpec, del1, del2, fcder, jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
      }
    }

  } /* k-loop */
}

void hdnl_and_invtr()
{
  assert(n2m % vreal::WIDTH == 0);

  const real alre  = al/ren;
  const real alpec = al/pec;

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc++)
  {
    const int kp  = kc + 1;
    const int km  = kc - 1;
    const int kmm = kc == 0     ? kc : km;
    const int kpp = kc == n3m-1 ? kc : kp;

    /* hdnl */
    
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


    /* invtr2 */

    foreach_warp(jc, 0, n2m)
    {
      const int jp  = jc+1;
      const int jm  = jc-1;
      if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
        foreach_lane(ch)
        {
          const int jmm = jmv(jc+ch)-1;
          const int jpp = jpv(jc+ch)-1;
          invtr2_ij<sreal>(alre, jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
        }
      else
        invtr2_ij<vreal>(alre, jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
    }

    /* invtr3 */

    foreach_warp(jc, 0, n2m)
    {
      const int jp  = jc+1;
      const int jm  = jc-1;
      if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
        foreach_lane(ch)
        {
          const int jmm = jmv(jc+ch)-1;
          const int jpp = jpv(jc+ch)-1;
          invtr3_ij<sreal>(alre, jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
        }
      else
          invtr3_ij<vreal>(alre, jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
    }

    /* invtro */

    if (kc >= 1 && kc <= n3m-2)
    {
      foreach_warp(jc, 0, n2m)
      {
        const int jp  = jc+1;
        const int jm  = jc-1;
        if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
          foreach_lane(ch)
          {
            const int jmm = jmv(jc+ch)-1;
            const int jpp = jpv(jc+ch)-1;
            invtrro_ij<sreal>(alpec, jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
          }
        else
          invtrro_ij<vreal>(alpec, jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
      }
    }
    else if (kc == 0)
    {
      const real del1 = zm(0) - zz(0);
      const real del2 = zm(1) - zm(0);
      const real fcder = 2.0/(del1*del2*(del1+del2));
      foreach_warp(jc, 0, n2m)
      {
        const int jp  = jc+1;
        const int jm  = jc-1;
        if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
          foreach_lane(ch)
          {
            const int jmm = jmv(jc+ch)-1;
            const int jpp = jpv(jc+ch)-1;
            invtrroL_ij<sreal>(alpec, del1, del2, fcder, jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
          }
        else
          invtrroL_ij<vreal>(alpec, del1, del2, fcder, jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
      }
    }
    else
    {
      const real  del1 = zz(n3m  ) - zm(n3m-1);
      const real  del2 = zm(n3m-1) - zm(n3m-2);
      const real fcder = 2.0/(del1*del2*(del1+del2));
      foreach_warp(jc, 0, n2m)
      {
        const int jp  = jc+1;
        const int jm  = jc-1;
        if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
          foreach_lane(ch)
          {
            const int jmm = jmv(jc+ch)-1;
            const int jpp = jpv(jc+ch)-1;
            invtrroU_ij<sreal>(alpec, del1, del2, fcder, jc+ch, kc, jm+ch, km, jp+ch, kp, jmm, kmm, jpp, kpp);
          }
        else
          invtrroU_ij<vreal>(alpec, del1, del2, fcder, jc, kc, jm, km, jp, kp, jm, kmm, jp, kpp);
      }
    }

  } /* k-loop */
}
