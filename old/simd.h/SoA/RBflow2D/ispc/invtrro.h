#pragma once

/************************************************************************
 *                       SUBROUTINE INVTRRO
 *   This subroutine performs the computation of he scalar field.
 *   For details see the introduction of INVTR1
 ************************************************************************/

void invtrro_scal(const int kbeg, const int kend, const int jbeg, const int jend)
{
  const real alpec = al/pec;

  /*  compute the rhs of the factored equation
   *  everything at i+1/2,j+1/2,k+1/2
   */
#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {
    if (kc >= 1 && kc <= n3m-2)
    {
      const int km = kmv(kc)-1;
      const int kp = kpv(kc)-1;
      for (int jc = jbeg; jc < jend; jc++)
      {
        const int jmm = jmv(jc)-1;
        const int jpp = jpv(jc)-1;

        const real aap = apscj(jc);
        const real aam = amscj(jc);
        const real aac = acscj(jc);

        /*   22 second derivatives of dens  */
        const real dq32 = dens(jpp,kc)*aap+dens(jc,kc)*aac + dens(jmm,kc)*aam;

        /*   33 second derivatives of dens   */
        const real dq33=dens(jc,kp)*ap3ssk(kc) + dens(jc,kc)*ac3ssk(kc) + dens(jc,km)*am3ssk(kc);

        const real dcq3=dq32+dq33;

        /*    right hand side of the density equation */
        /*m===========================================================*/
        rhs3(jc,kc)=(ga*hro(jc,kc)+ro*ruro(jc,kc) + alpec*dcq3)*dt;
        /*m===========================================================*/

        /*    updating of the non-linear terms */
        ruro(jc,kc)=hro(jc,kc);
      }
    }

    {
      /*
       *       LOWER HOT WALL      
       */

      const real del1 = zm(0) - zz(0);
      const real del2 = zm(1) - zm(0);
      const real fcder = 2.0/(del1*del2*(del1+del2));
      if (kc == 0)
      {
        const int kp = kc + 1;
        for (int jc = jbeg; jc < jend; jc++)
        {
          const int jmm = jmv(jc)-1;
          const int jpp = jpv(jc)-1;

          const real aap = apscj(jc);
          const real aam = amscj(jc);
          const real aac = acscj(jc);

          /*   22 second derivatives of dens */
          const real dq32 = dens(jpp,kc)*aap + dens(jc,kc)*aac + dens(jmm,kc)*aam;

          /*   33 second derivatives of dens */
          const real dq33=(dens(jc,kp)*del1-dens(jc,kc)*(del1+del2) + denbs(jc)*del2) *fcder;
          const real dcq3=dq32+dq33;
          /*    right hand side of the density equation */
          //m=======================================================
          rhs3(jc,kc)=(ga*hro(jc,kc)+ro*ruro(jc,kc) + alpec*dcq3)*dt;
          //cm=======================================================
          //
          //   updating of the non-linear terms
          //
          ruro(jc,kc)=hro(jc,kc);
        }
      }
    }
    
    {
      //
      //       UPPER COLD WALL
      //     
      //
      const real  del1 = zz(n3 -1) - zm(n3m-1);
      const real  del2 = zm(n3m-1) - zm(n3m-2);
      const real fcder = 2.0/(del1*del2*(del1+del2));

      if (kc == n3m-1)
      {
        const int km = kc - 1;
        for (int jc = jbeg; jc < jend; jc++)
        {
          const int jmm = jmv(jc)-1;
          const int jpp = jpv(jc)-1;
          const real aap = apscj(jc);
          const real aam = amscj(jc);
          const real aac = acscj(jc);

          /*   22 second derivatives of dens  */
          const real dq32 = dens(jpp,kc)*aap + dens(jc,kc)*aac + dens(jmm,kc)*aam;

          /*   33 second derivatives of dens  */
          const real dq33 = (dens(jc,km)*del1 - dens(jc,kc)*(del1+del2) + denbn(jc)*del2) * fcder;
          const real dcq3 = dq32+dq33;

          //
          //    right hand side of the density equation
          //
          //m========================================================
          rhs3(jc,kc)=(ga*hro(jc,kc) + ro*ruro(jc,kc) + alpec*dcq3)*dt;
          //m========================================================
          //
          //    updating of the non-linear terms
          //
          ruro(jc,kc)=hro(jc,kc);
        }
      }
    }
  }
}

void invtrro()
{
  assert(n2m % SIMD::WIDTH == 0);
  assert(n3m % SIMD::WIDTH == 0);

  const int offset = 1;
  
  invtrro_scal(         0,     offset, 0,  n2m);
  invtrro_scal(n3m-offset, n3m,        0,  n2m);
  invtrro_scal(    offset, n3m-offset, 0,  SIMD::WIDTH);
  invtrro_scal(    offset, n3m-offset, n2m-SIMD::WIDTH, n2m);

  const SIMD::real alpec = al/pec;
  const SIMD::real ga(this->ga);
  const SIMD::real ro(this->ro);
  const SIMD::real dt(this->dt);

  /*  compute the rhs of the factored equation
   *  everything at i+1/2,j+1/2,k+1/2
   */
#pragma omp parallel for
  for (int kc = offset; kc < n3m-offset; kc++)
  {
    const int km = kmv(kc)-1;
    const int kp = kpv(kc)-1;
    const SIMD::real ap3 = ap3ssk(kc);
    const SIMD::real ac3 = ac3ssk(kc);
    const SIMD::real am3 = am3ssk(kc);

    for (int jc = SIMD::WIDTH; jc < n2m - SIMD::WIDTH; jc += SIMD::WIDTH)
    {
      const int jmm = jc+1;
      const int jpp = jc-1;

      const SIMD::real aap(&apscj(jc));
      const SIMD::real aam(&amscj(jc));
      const SIMD::real aac(&acscj(jc));

      /*   22 second derivatives of dens  */

      const SIMD::real dq32 = 
        LDU(dens(jpp,kc))*aap + 
        LDA(dens(jc ,kc))*aac + 
        LDU(dens(jmm,kc))*aam;

      /*   33 second derivatives of dens   */
      const SIMD::real dq33 = 
        LDA(dens(jc,kp))*ap3 +
        LDA(dens(jc,kc))*ac3 + 
        LDA(dens(jc,km))*am3;

      const SIMD::real dcq3=dq32 + dq33;
      
      const SIMD::real &hro  = STA(this->hro (jc,kc));
      __out SIMD::real &ruro = STA(this->ruro(jc,kc));
      __out SIMD::real &rhs  = STA(this->rhs3(jc,kc));

      /*    right hand side of the density equation */
      /*m===========================================================*/
      rhs = (ga*hro + ro*ruro + alpec*dcq3)*dt;
      /*m===========================================================*/

      /*    updating of the non-linear terms */
      ruro = hro;
    }
  }

  solroj();
  solrok();
}
