#pragma once

void invtr2scal(const int kbeg, const int kend, const int jbeg, const int jend)
{
  const real alre = al/ren;
#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {
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
}

void invtr2()  /* check */
{

  /*
   *  compute the rhs of the factored equation
   *  everything at i+1/2,j,k+1/2
   *
   *    points inside the flowfield
   */

  assert(n2m % SIMD::WIDTH == 0);
  
  invtr2scal(0, n3m,               0,     SIMD::WIDTH);
  invtr2scal(0, n3m, n2m-SIMD::WIDTH, n2m            );

  const SIMD::real alre = al/ren;
  const SIMD::real al(this->al);
  const SIMD::real ga(this->ga);
  const SIMD::real ro(this->ro);
  const SIMD::real dt(this->dt);

#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    double t0 = rtc();
#pragma omp for
    for (int kc = 0; kc < n3m; kc++)
    {
      const int km = kc == 0     ?   0   : kc - 1;
      const int kp = kc == n3m-1 ? n3m-1 : kc + 1;
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
    if (tid == 0)
      mytimer[39] += rtc() - t0;
  }

  double t0 = rtc();
  solq2j();
  mytimer[49] += rtc() - t0;

  t0 = rtc();
  solq2k(); 
  mytimer[50] += rtc() - t0;

}  

