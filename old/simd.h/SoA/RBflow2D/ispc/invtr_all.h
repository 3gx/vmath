#pragma once

void invtr_all_scal(const int kbeg, const int kend, const int jbeg, const int jend)
{

#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {

    /*
     *  compute the rhs of the factored equation
     *  everything at i+1/2,j,k+1/2
     *
     *    points inside the flowfield
     */
    {
      const real alre = al/ren;
#if 0
      const int km = kmv(kc)-1;
      const int kp = kpv(kc)-1;
#else
      const int km = kc == 0     ?   0   : kc - 1;
      const int kp = kc == n3m-1 ? n3m-1 : kc + 1;
#endif
      for (int jc = jbeg; jc < jend; jc++)
      {
#if 0
        const int jm = jmv(jc)-1;
        const int jp = jpv(jc)-1;
#else
        const int jm = jc == 0     ?  n2m-1 : jc - 1;
        const int jp = jc == n2m-1 ?  0     : jc + 1;
#endif

        const real udx2 = al*udx2c(jc);

        /*
         *   22 second derivative of q2
         */
        const real d22q2 = 
          q2(jp,kc) * ap2je(jc) + 
          q2(jc,kc) * ac2je(jc) +
          q2(jm,kc) * am2je(jc);

        /*
         *   33 second derivative of q2
         */
        const real d33q2 =
          q2(jc,kp) * ap3sk(kc) +
          q2(jc,kc) * ac3sk(kc) +
          q2(jc,km) * am3sk(kc);

        /*
         *    viscid terms
         */
        const real dcq2 = d22q2 + d33q2;

        /*
         *   component of grad(pr) along 2 direction
         */
        const real dpx22 = (pr(jc,kc) - pr(jm,kc)) * udx2;

        rhs1(jc,kc) = (ga*dph(jc,kc) + ro*ru2(jc,kc) + alre*dcq2 - dpx22) * dt;

        ru2(jc,kc) = dph(jc,kc);
      }
    } /* 21 flops */


#if 1
    /*  compute the rhs of the factored equation
     *  everything at i+1/2,j+1/2,k
     */

    {
      const real alre=al/ren;

      const int  km   = kmv(kc)-1;
      const int  kp   = kc+1;
      const real udx3 = al*udx3c(kc);
      for (int jc = jbeg; jc < jend; jc++)
      {
        const real jmm =  jmv(jc)-1;
        const real jpp =  jpv(jc)-1;
        const real aap = ap3j(jc);
        const real aam = am3j(jc);
        const real aac = ac3j(jc);

        /*   22 second derivatives of q3 */
        const real dq32=aam*q3(jmm,kc)+aac*q3(jc,kc)+aap*q3(jpp,kc);

        /*   33 second derivatives of q3 */
        const real dq33=q3(jc,kp)*ap3ck(kc) + q3(jc,kc)*ac3ck(kc) + q3(jc,km)*am3ck(kc);

        /*   viscous terms    */
        const real dcq3=dq32+dq33;

        /*  component of grad(pr) along x3 direction */
        const real dpx33=(pr(jc,kc)-pr(jc,km))*udx3;

        /*======================================================= */
        rhs2(jc,kc)=(ga*qcap(jc,kc)+ro*ru3(jc,kc) + alre*dcq3-dpx33)*dt;
        /*======================================================= */

        /*  updating of the non-linear terms */
        ru3(jc,kc)=qcap(jc,kc);
      }
    }  /* 20 flops */

    /************************************************************************
     *                       SUBROUTINE INVTRRO
     *   This subroutine performs the computation of he scalar field.
     *   For details see the introduction of INVTR1
     ************************************************************************/

    /*  compute the rhs of the factored equation
     *  everything at i+1/2,j+1/2,k+1/2
     */
    {
      const real alpec = al/pec;
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
#endif
  }
}

void invtr_all()
{


  assert(n2m % SIMD::WIDTH == 0);
  
  int nthreads = 1;
#pragma omp parallel
  nthreads = omp_get_num_threads();

  const int offset = nthreads;

  invtr_all_scal(         0,     offset, 0,  n2m);
  invtr_all_scal(n3m-offset, n3m,        0,  n2m);
  invtr_all_scal(    offset, n3m-offset, 0,  SIMD::WIDTH);
  invtr_all_scal(    offset, n3m-offset, n2m-SIMD::WIDTH, n2m);

  const SIMD::real al(this->al);
  const SIMD::real ga(this->ga);
  const SIMD::real ro(this->ro);
  const SIMD::real dt(this->dt);


#pragma omp parallel for
  for (int kc = offset; kc < n3m-offset; kc++)
  {

    const int km = kc == 0     ?   0   : kc - 1;
    const int kp = kc == n3m-1 ? n3m-1 : kc + 1;

    /*  compute the rhs of the factored equation
     *  everything at i+1/2,j,k+1/2
     *
     *    points inside the flowfield
     */
    {
      const SIMD::real alre = al/SIMD::real(ren);
#if 0
      const SIMD::real al(this->al);
      const SIMD::real ga(this->ga);
      const SIMD::real ro(this->ro);
      const SIMD::real dt(this->dt);

      const int km = kc == 0     ?   0   : kc - 1;
      const int kp = kc == n3m-1 ? n3m-1 : kc + 1;
#endif
      const SIMD::real ap3 = ap3sk(kc);
      const SIMD::real ac3 = ac3sk(kc);
      const SIMD::real am3 = am3sk(kc);
      for (int jc = SIMD::WIDTH; jc < n2m - SIMD::WIDTH; jc += SIMD::WIDTH)
      {
        const int jm = jc - 1;
        const int jp = jc + 1;

        const SIMD::real udx2 = LDA(udx2c(jc)) * al;

        const SIMD::real am2 = LDA(am2je(jc));
        const SIMD::real ac2 = LDA(ac2je(jc));
        const SIMD::real ap2 = LDA(ap2je(jc));

        const SIMD::real q2mc = LDU(q2(jm,kc));
        const SIMD::real q2pc = LDU(q2(jp,kc));
        const SIMD::real q2cp = LDA(q2(jc,kp));
        const SIMD::real q2cc = LDA(q2(jc,kc));
        const SIMD::real q2cm = LDA(q2(jc,km));

        /*
         *   22 second derivative of q2
         */

        const SIMD::real d22q2 = q2pc*ap2 + q2cc*ac2 + q2mc*am2;

        /*
         *   33 second derivative of q2
         */

        const SIMD::real d33q2 = q2cp*ap3 + q2cc*ac3 + q2cm*am3;

        /*
         *    viscid terms
         */
        const SIMD::real dcq2 = d22q2 + d33q2;

        /*
         *   component of grad(pr) along 2 direction
         */
        const SIMD::real pcc(&pr(jc,kc), true);
        const SIMD::real pmc(&pr(jm,kc), false);
        const SIMD::real dpx22 = (pcc - pmc) * udx2;

        const SIMD::real &dph = STA(this->dph(jc,kc));
        __out SIMD::real &ru2 = STA(this->ru2(jc,kc));
        __out SIMD::real &rhs = STA(this->rhs1(jc,kc));

        rhs = (ga*dph + ro*ru2 + alre*dcq2 - dpx22) * dt;
        ru2 = dph;
      }
    } /* 21 flops */

    /*  compute the rhs of the factored equation
     *  everything at i+1/2,j+1/2,k
     */
    {
      const SIMD::real alre = al/SIMD::real(ren);
#if 0
      const SIMD::real ga(this->ga);
      const SIMD::real ro(this->ro);
      const SIMD::real dt(this->dt);
#if 0
      const int  km   = kmv(kc)-1;
      const int  kp   = kc+1;
#else
      const int km = kc == 0     ?   0   : kc - 1;
      const int kp = kc == n3m-1 ? n3m-1 : kc + 1;
#endif
#endif
      const SIMD::real udx3 = udx3c(kc) * al;
      const SIMD::real  ap3 = ap3ck(kc);
      const SIMD::real  ac3 = ac3ck(kc);
      const SIMD::real  am3 = am3ck(kc);

      for (int jc = SIMD::WIDTH; jc < n2m-SIMD::WIDTH; jc += SIMD::WIDTH)
      {
        const int jmm =  jc-1;
        const int jpp =  jc+1;

        const SIMD::real aap(&ap3j(jc));
        const SIMD::real aam(&am3j(jc));
        const SIMD::real aac(&ac3j(jc));

        /*   22 second derivatives of q3 */
        const SIMD::real dq32 =
          aam*LDU(q3(jmm,kc)) + 
          aac*LDA(q3(jc, kc)) +
          aap*LDU(q3(jpp,kc));

        /*   33 second derivatives of q3 */
        const SIMD::real dq33 = 
          ap3 * LDA(q3(jc,kp)) + 
          ac3 * LDA(q3(jc,kc)) +
          am3 * LDA(q3(jc,km));

        /*   viscous terms    */
        const SIMD::real dcq3=dq32+dq33;

        /*  component of grad(pr) along x3 direction */
        const SIMD::real dpx33 = (LDA(pr(jc,kc)) - LDA(pr(jc,km)))*udx3;

        const SIMD::real &qcap = STA(this->qcap(jc,kc));
        __out SIMD::real &ru3  = STA(this->ru3 (jc,kc));
        __out SIMD::real &rhs  = STA(this->rhs2(jc,kc));

        /*======================================================= */
        rhs = (ga*qcap + ro*ru3 + alre*dcq3 - dpx33)*dt;
        /*======================================================= */

        /*  updating of the non-linear terms */
        ru3 = qcap;
      }
    }  /* 20 flops */

    /*  compute the rhs of the factored equation
     *  everything at i+1/2,j+1/2,k+1/2
     */
    { 
      const SIMD::real alpec = al/SIMD::real(pec);
#if 0
      const SIMD::real ga(this->ga);
      const SIMD::real ro(this->ro);
      const SIMD::real dt(this->dt);
#if 1
      const int km = kmv(kc)-1;
      const int kp = kpv(kc)-1;
#else
      const int km = kc == 0     ?   0   : kc - 1;
      const int kp = kc == n3m-1 ? n3m-1 : kc + 1;
#endif
#endif
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

  } 

#if 0
  double t0 = rtc();
  solq2j();
  mytimer[49] += rtc() - t0;

  t0 = rtc();
  solq2k(); 
  mytimer[50] += rtc() - t0;

#if 1
  solq3j();
  solq3k();

  solroj();
  solrok();
#endif
#endif

}  
