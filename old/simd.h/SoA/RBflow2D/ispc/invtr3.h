#pragma once

void invtr3scal(const int kbeg, const int kend, const int jbeg, const int jend)
{
  const real alre=al/ren;

  /*  compute the rhs of the factored equation
   *  everything at i+1/2,j+1/2,k
   */

#pragma omp parallel for
  for (int kc = kbeg; kc < kend; kc++)
  {
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

}

void invtr3()  /* check */
{

  invtr3scal(1, n3m,               0,     SIMD::WIDTH);
  invtr3scal(1, n3m, n2m-SIMD::WIDTH, n2m            );

  const SIMD::real alre = al/ren;
  const SIMD::real ga(this->ga);
  const SIMD::real ro(this->ro);
  const SIMD::real dt(this->dt);
#pragma omp parallel for
  for (int kc = 1; kc < n3m; kc++)
  {
    const int  km   = kmv(kc)-1;
    const int  kp   = kc+1;
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


  solq3j();  /* 20 flops */
  solq3k();  /* 10 flops */
}  /* 50 flops */
